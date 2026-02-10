// Author: Martin C. Frith 2025
// SPDX-License-Identifier: BSD-3-Clause

// See [Frith2025]: "Simple and thorough detection of related
// sequences with position-varying probabilities of substitutions,
// insertions, and deletions", MC Frith 2025

#include "dummer-util.hh"
#include "tantan-wrapper.hh"
#include "can_i_haz_simd.hh"

#include <algorithm>
#include <array>
#include <iomanip>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>
#include <queue>

#define OPT_e 10
#define OPT_s 2
#define OPT_m 0 // mask-filter doesnt work well on DNA
#define OPT_t 1000
#define OPT_l 1000
#define OPT_b 100

//#define DUMMER_SCORING
#define EVALUE

#ifdef DOUBLE
typedef double Float;
typedef SimdDbl SimdFloat;
const int simdLen = simdDblLen;
#else
typedef float Float;
typedef SimdFlt SimdFloat;
const int simdLen = simdFltLen;
#endif

const Float STOP_CODON_PROB = 0.01;
const Float FRAMESHIFT1_MULTIPLIER = 0.01;
const Float FRAMESHIFT2_MULTIPLIER = 0.005;
#define XX ((FRAMESHIFT1_MULTIPLIER + FRAMESHIFT2_MULTIPLIER) * 1)

int simdRoundUp(int x) {  // lowest multiple of simdLen that is >= x
  return x - 1 - (x - 1) % simdLen + simdLen;
}

Float simdHorizontalMax(SimdFloat x) {  // assuming it doesn't need to be fast
  Float y[simdLen];
  simdStore(y, x);
  return *std::max_element(y, y + simdLen);
}

SimdFloat simdPowersFwd(Float x) {
  Float a[simdLen];
  a[0] = x;
  for (int i = 1; i < simdLen; ++i) a[i] = a[i-1] * x;
  return simdLoad(a);
}

SimdFloat simdPowersRev(Float x) {
  Float a[simdLen];
  a[simdLen-1] = x;
  for (int i = simdLen-1; i > 0; --i) a[i-1] = a[i] * x;
  return simdLoad(a);
}

// Only consider similarities that are local maxima.  If 2
// similarities have identical 1st anchor coordinates, and their 2nd
// anchor coordinates are closer than this, omit the lower-scoring one.
const int minSeparation = 32;  // xxx ???

// down-scale probabilities by this amount, to delay overflow:
const Float scale = 1.0 / (1<<30) / (1<<30) / (1<<3); // sqrt[min normal float]
const int shift = 63;  // add this to scores, to undo the scaling

int verbosity = 0;

const int nonLetterWidth = 8;  // number of non-letter values per position

struct Params { // TODO: maybe SIMD order
  Float alpha_prime[3];
  Float beta_prime[3];
  Float delta_prime[3];
  Float epsilon_prime[3];
  Float enter_match_probability;
};

struct Profile {  // position-specific (insert, delete, letter) probabilities
  Float *values;  // probabilities or probability ratios
  std::vector<Params> values_v2;
  std::vector<Float> bg_probs;
  int width;   // number of values per position
  int length;  // number of positions
  size_t nameIdx;
  size_t consensusSequenceIdx;
  double gumbelKendAnchored, gumbelKbegAnchored, gumbelKmidAnchored;
  void *debug;
};

struct Sequence {
  size_t nameIdx;
  int length;
};

struct Contig {
  int start;
  int length;
};

struct SegmentPair {
  int start1, start2, length;
};

struct InitialSimilarity {
  Float probRatio;
  int anchor2;  // 2nd anchor coordinate (don't need to store the 1st one)
};

struct AlignedSimilarity {
  double probRatio;
  int anchor1, anchor2;
  Float wEndAnchored;
  std::vector<SegmentPair> alignment;
};

struct FinalSimilarity {
  double probRatio;
  size_t profileNum;
  size_t strandNum;
  int anchor1, anchor2;
  int start1, start2;
  std::vector<char> alignedSequences;
};

int simBeg2(const AlignedSimilarity &x) {
  return x.alignment.empty() ? x.anchor2 : x.alignment[0].start2;
}

int simEnd2(const AlignedSimilarity &x) {
  return x.alignment.empty() ?
    x.anchor2 : x.alignment.back().start2 + x.alignment.back().length;
}

double mean(const double *x, int n) {
  double s = 0;
  for (int i = 0; i < n; ++i) s += x[i];
  return s / n;
}

int numOfDigits(int x) {
  int n = 0;
  do ++n; while (x /= 10);
  return n;
}

const char *getAlphabet(int alphabetSize) {
  assert(alphabetSize == 20 || alphabetSize == 4);
  return alphabetSize == 20 ? "ACDEFGHIKLMNPQRSTVWYUO"  // 20 + 2 amino acids
    :    alphabetSize ==  4 ? "ACGT" : 0;
}

char complement(char c) {
  return 3 ^ c;
}

void reverseComplement(char *beg, char *end) {
  while (beg < end) {
    char c = *--end;
    *end = complement(*beg);
    *beg++ = complement(c);
  }
}

std::istream &readContig(std::istream &in, Sequence &sequence, Contig &contig,
			 std::vector<char> &vec, const char *charToNumber) {
  if (contig.length == 0) {
    char x;
    if (!(in >> x)) return in;
    if (x != '>') return fail(in, "bad sequence data: no '>'");
    std::string line, word;
    getline(in, line);
    std::istringstream iss(line);
    if (!(iss >> word)) return fail(in, "bad sequence data: no name");
    sequence.nameIdx = vec.size();
    sequence.length = 0;
    const char *name = word.c_str();
    vec.insert(vec.end(), name, name + word.size() + 1);
    if (verbosity > 0) std::cerr << "Sequence: " << name << "\n";
  }

  size_t seqIdx = vec.size();
  std::streambuf *buf = in.rdbuf();
  int c = buf->sgetc();

  while (c != std::streambuf::traits_type::eof() && c != '>') {
    if (charToNumber[c] < 125) break;  // found a contig symbol
    if (c > ' ') ++sequence.length;  // skip over non-contig symbols
    c = buf->snextc();
  }

  while (c != std::streambuf::traits_type::eof() && charToNumber[c] < 126) {
    if (c > ' ') {
        vec.push_back(charToNumber[c]);
    }
    c = buf->snextc();
  }

  size_t seqLen = vec.size() - seqIdx;
  if (seqLen > INT_MAX - 2 * simdLen) return fail(in, "sequence is too long!");
  contig.start = sequence.length;
  contig.length = seqLen;
  sequence.length += seqLen;
  return in;
}

char profileLetter(const char *alphabet, char letterCode) {
  return alphabet[letterCode & 31] + (letterCode & 32);  // upper/lowercase
}

void addAlignedProfile(std::vector<char> &gappedSeq,
		       const std::vector<SegmentPair> &alignment,
		       const char *alphabet, const char *consensusSequence) {
  int pos1 = alignment[0].start1;
  int pos2 = alignment[0].start2;
  for (auto a : alignment) {
    for (; pos1 < a.start1; ++pos1) {
      gappedSeq.push_back(profileLetter(alphabet, consensusSequence[pos1]));
    }
    gappedSeq.insert(gappedSeq.end(), a.start2 - pos2, '-');
    for (; pos1 < a.start1 + a.length; ++pos1) {
      gappedSeq.push_back(profileLetter(alphabet, consensusSequence[pos1]));
    }
    pos2 = a.start2 + a.length;
  }
}

char seqLetter(const char *alphabet, const char *sequence,
	       const char *maskedSequence, int position) {
  char c = sequence[position];
  return alphabet[c] + (maskedSequence[position] > c) * 32;  // upper/lowercase
}

void addAlignedSequence(std::vector<char> &gappedSeq,
			const std::vector<SegmentPair> &alignment,
			const char *alphabet, const char *sequence,
			const char *maskedSequence) {
  int pos1 = alignment[0].start1;
  int pos2 = alignment[0].start2;
  for (auto a : alignment) {
    gappedSeq.insert(gappedSeq.end(), a.start1 - pos1, '-');
    for (; pos2 < a.start2; ++pos2) {
      gappedSeq.push_back(seqLetter(alphabet, sequence, maskedSequence, pos2));
    }
    for (; pos2 < a.start2 + a.length; ++pos2) {
      gappedSeq.push_back(seqLetter(alphabet, sequence, maskedSequence, pos2));
    }
    pos1 = a.start1 + a.length;
  }
}

int strandPosition(size_t strandNum, int seqLength, int position) {
  return (strandNum % 2) ? seqLength - position : position;
}

Float xx_null;
void printSimilarity(const char *names, Profile p, Sequence s,
		     const FinalSimilarity &sim, double evalue) {
  char strand = "+-"[sim.strandNum % 2];
  const char *seq = sim.alignedSequences.data();
  int length = sim.alignedSequences.size() / 2;
  int span1 = length - std::count(seq, seq + length, '-');
  int span2 = length - std::count(seq + length, seq + length * 2, '-');
  int start2 = strandPosition(sim.strandNum, s.length, sim.start2);
  int anchor2 = strandPosition(sim.strandNum, s.length, sim.anchor2);
  int w1 = std::max(strlen(names + p.nameIdx), strlen(names + s.nameIdx));
  int w2 = std::max(numOfDigits(sim.start1), numOfDigits(start2));
  int w3 = std::max(numOfDigits(span1), numOfDigits(span2));
  int w4 = std::max(numOfDigits(p.length), numOfDigits(s.length));
  std::cout << "a score=" << (log2(sim.probRatio)+shift) << " E=" << evalue
	    << " anchor=" << sim.anchor1 << "," << anchor2 << "\n";
  std::cout << "s " << std::left << std::setw(w1) << names + p.nameIdx << " "
	    << std::right << std::setw(w2) << sim.start1 << " "
	    << std::setw(w3) << span1 << " " << '+' << " "
	    << std::setw(w4) << p.length << " ";
  std::cout.write(seq, length);
  std::cout << "\n";
  std::cout << "s " << std::left << std::setw(w1) << names + s.nameIdx << " "
	    << std::right << std::setw(w2) << start2 << " "
	    << std::setw(w3) << span2 << " " << strand << " "
	    << std::setw(w4) << s.length << " ";
  std::cout.write(seq + length, length);
  std::cout << "\n\n";
}

void addForwardMatch(std::vector<SegmentPair> &alignment, int pos1, int pos2) {
  if (!alignment.empty()) {
    SegmentPair &x = alignment.back();
    if (x.start1 + x.length == pos1 && x.start2 + x.length == pos2) {
      ++x.length;
    }
    if (x.start1 + x.length > pos1 || x.start2 + x.length > pos2) return;
  }
  SegmentPair sp = {pos1, pos2, 1};
  alignment.push_back(sp);
}

void addReverseMatch(std::vector<SegmentPair> &alignment, int pos1, int pos2) {
  if (!alignment.empty()) {
    SegmentPair &x = alignment.back();
    if (x.start1 - 1 == pos1 && x.start2 - 1 == pos2) {
      --x.start1;
      --x.start2;
      ++x.length;
    }
    if (x.start1 <= pos1 || x.start2 <= pos2) return;
  }
  SegmentPair sp = {pos1, pos2, 1};
  alignment.push_back(sp);
}

using DP_2D = std::vector<std::vector<Float>>;
DP_2D make_dp_table(size_t rows, size_t cols) {
    return DP_2D(rows, std::vector<Float>(cols));
}

struct metadata {
  Float metric = 0;
  void *dest = nullptr;
  int dest_i, dest_j;

  metadata(Float init) : metric(init) {
    assert(!isnan(init));
  }

  metadata add_cost(Float cost) const {
    metadata ret = *this;
    ret.metric += cost;
    assert(!isnan(cost));

    return ret;
  };

  constexpr bool operator < (const metadata& r) const noexcept {
    return this->metric < r.metric;
  }
};

using DP_2Dv2 = std::vector<std::vector<metadata>>;
DP_2Dv2 make_dp_table_v2(size_t rows, size_t cols) {
    return DP_2Dv2(rows, std::vector<metadata>(cols, 0));
}

inline double dp_access_safe(DP_2D &dp, int i, int j, double err = 0) {
  if(0 <= i && i < dp.size() && 0 <= j && j < dp[0].size()) {
    return dp[i][j];
  }
  return err;
}

struct {
  size_t r = -1, c = -1;
  std::array<DP_2D, 2> W;
  DP_2D X;
  std::array<DP_2D, 3> Y, Z;

  void resize_if_need(size_t r, size_t c, bool reset_w = true) {
    if(reset_w) {
      W = {
        make_dp_table(r, c),
        make_dp_table(r, c)
      };
    }
    
    if(r != this->r || c != this->c) {
      X = make_dp_table(r, c);
      Y = {
        make_dp_table(r, c),
        make_dp_table(r, c),
        make_dp_table(r, c)
      };
      Z = {
        make_dp_table(r, c),
        make_dp_table(r, c),
        make_dp_table(r, c)
      };
    }
  }
} scratch_v2;
std::vector<std::pair<int, Float>> decoded;

template <typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U>& a,
                          const std::pair<T, U>& b)
{
    return { a.first + b.first, b.second};
}
void addForwardAlignment(std::vector<SegmentPair> &alignment,
			 Profile profile, const char *sequence,
			 int sequenceLength, const Float *scratch,
			 int iBeg, int jBeg, double half) {
  /*
  TODO memory optimization: per-diagonal indexing using row (profile) indexes
  0 0 0 0
  1 1 1
  2 2
  3

  base cases:
    idx = 0 X and Y are 0 (only insertions allowed)
    idx = row #: X and Z are 0 (only deletions allowed)
  */

  //scratch_v2.resize_if_need(std::min(profile.length + 1, sequenceLength + 1), std::min(profile.length + 1, sequenceLength + 1));
  
  std::array<DP_2Dv2, 2> W{
      make_dp_table_v2(profile.length + 1, sequenceLength + 1),
      make_dp_table_v2(profile.length + 1, sequenceLength + 1)
  };
  DP_2Dv2 X = make_dp_table_v2(profile.length + 1, sequenceLength);
  std::array<DP_2Dv2, 3> Y{
      make_dp_table_v2(profile.length + 1, sequenceLength),
      make_dp_table_v2(profile.length + 1, sequenceLength),
      make_dp_table_v2(profile.length + 1, sequenceLength)
  };
  std::array<DP_2Dv2, 3> Z{
      make_dp_table_v2(profile.length + 1, sequenceLength),
      make_dp_table_v2(profile.length + 1, sequenceLength),
      make_dp_table_v2(profile.length + 1, sequenceLength)
  };

  addForwardMatch(alignment, iBeg, jBeg);
  //std::cout << "anchor is " << iBeg << " and " << jBeg << std::endl;

  auto dp_access_safe_v2 = [&](DP_2Dv2 &dp, int i, int j) {
    if(0 <= i && i < dp.size() && 0 <= j && j < dp[0].size()) {
      metadata ret(dp[i][j].metric);
      ret.dest = &dp;
      ret.dest_i = i;
      ret.dest_j = j;
      return ret;
    }

    metadata ret(0 /* end here */);
    return ret;
  };

  for(int radius = std::min(profile.length - iBeg, sequenceLength - 1 - jBeg); radius >= 0; radius--) {
    //std::cout << "radius " << radius << std::endl;
    for(int i = iBeg, j = jBeg + radius; i <= profile.length && j >= jBeg; i++, j--) {
      if(j >= sequenceLength) continue;

      std::array<Float, 4> w;
      const Params &params_cur = profile.values_v2[i];
      const Params &params_later = profile.values_v2[i + 1];
      const Float *params = profile.values + (i) * profile.width;
      const Float *params_emission_probabilities = params + 4;
      
      Float emit_prob = -INFINITY;
      if(j + 1 < sequenceLength) {
        auto [emitNum, divisor] = decoded[j + 1]; // upto j emitted alr
        emit_prob = log(params_emission_probabilities[emitNum] * divisor);
        X[i][j] = dp_access_safe_v2(W[0], i + 1, j + 3).add_cost(emit_prob);
      } else {
        X[i][j] = metadata(-INFINITY);
      }

      Y[0][i][j] = std::max(dp_access_safe_v2(W[0], i + 1, j + 0), dp_access_safe_v2(Y[0], i + 1, j + 0).add_cost(log(params_later.epsilon_prime[0])));
      Y[1][i][j] = std::max(dp_access_safe_v2(W[0], i + 1, j + 2), dp_access_safe_v2(Y[0], i + 1, j + 2).add_cost(log(params_later.epsilon_prime[1]))).add_cost(2 * log(0.25));
      Y[2][i][j] = std::max(dp_access_safe_v2(W[0], i + 1, j + 1), dp_access_safe_v2(Y[0], i + 1, j + 1).add_cost(log(params_later.epsilon_prime[2]))).add_cost(1 * log(0.25));
      
      Z[0][i][j] = std::max(dp_access_safe_v2(W[0], i, j + 3), dp_access_safe_v2(Z[0], i, j + 3).add_cost(log(params_cur.beta_prime[0]))).add_cost(emit_prob);
      for(int k = 1; k <= 2; k++) {
        Z[k][i][j] = std::max(dp_access_safe_v2(W[0], i, j + k), dp_access_safe_v2(Z[0], i, j + k).add_cost(log(params_cur.beta_prime[k]))).add_cost(k * log(0.25));
      }

      W[0][i][j] = std::max({
        //metadata(0), // end here
        dp_access_safe_v2(X, i, j).add_cost(log(params_cur.enter_match_probability)),
        dp_access_safe_v2(Y[0], i, j).add_cost(log(params_cur.delta_prime[0])),
        dp_access_safe_v2(Y[1], i, j).add_cost(log(params_cur.delta_prime[1])),
        dp_access_safe_v2(Y[2], i, j).add_cost(log(params_cur.delta_prime[2])),
        dp_access_safe_v2(Z[0], i, j).add_cost(log(params_cur.alpha_prime[0])),
        dp_access_safe_v2(Z[1], i, j).add_cost(log(params_cur.alpha_prime[1])),
        dp_access_safe_v2(Z[2], i, j).add_cost(log(params_cur.alpha_prime[2]))});
    }
  }
  
  bool is_print = false;
  int i = iBeg, j = jBeg;
  auto cur = &X[iBeg][jBeg];
  while(true) {
    std::cout << i << " " << j << " " << (cur->metric) << std::endl;
    if(is_print) {
      addForwardMatch(alignment, i, j + 3);
    }

    if(cur->dest != nullptr) {
      is_print = cur->dest == &X;
      i = cur->dest_i;
      j = cur->dest_j;
      cur = &((*((DP_2Dv2*)cur->dest))[i][j]);
    } else {
      break;
    }
    
  } 

}

void addReverseAlignment(std::vector<SegmentPair> &alignment,
			 Profile profile, const char *sequence,
			 int sequenceLength, const Float *scratch,
			 int iEnd, int jEnd, double half) {
  long rowSize = simdRoundUp(sequenceLength + 1) + simdLen;
  const char *seq = sequence + jEnd;

  for (int size = 16; ; size *= 2) {
    int iBeg = std::max(iEnd - size, -1);
    int jBeg = std::max(jEnd - size, -1);
    int jLen = jEnd - jBeg;
    std::vector<double> scratch2(jLen * 2);  // use "double" to avoid underflow
    double *X = scratch2.data() + jLen;
    double *Y = X + jLen;
    double wSum = 0;

    for (int i = iEnd-1; i >= iBeg; --i) {
      const Float *Xforward = (i >= 0) ? scratch + i * rowSize + jEnd : 0;
      double a = profile.values[(i+1) * profile.width + 0];
      double b = profile.values[(i+1) * profile.width + 1];
      double d = profile.values[(i+1) * profile.width + 2];
      double e = profile.values[(i+1) * profile.width + 3];
      const Float *S = (i >= 0) ? profile.values + i * profile.width + 4 : 0;

      double x = (i == iEnd-1) ? scale : 0;
      double z = 0;
      for (int j = -1; j >= -jLen; --j) {
	double y = Y[j];
	double w = x + d * y + a * z;  // this is: W[i+1][jEnd+j+1]
	wSum += w;
	Y[j] = w + e * y;
	x = X[j];
	z = w + b * z;
	if (i < 0 || jEnd + j < 0) continue;
	X[j] = S[seq[j]] * w;
	if (Xforward[j] * w > half * scale) {
	  addReverseMatch(alignment, i, jEnd + j);
	}
      }

      if (wSum >= half) return;
    }

    // this line should be unnecessary, except for problems such as overflow:
    if (iBeg <= 0 && jBeg <= 0) break;
  }
}

bool maybeLocalMaximum(Profile profile, const char *sequence,
		       int sequenceLength, const Float *scratch,
		       int anchor1, int anchor2, Float wMidAnchored) {
  Float X[minSeparation * 2 - 1];
  Float Y[minSeparation * 2 - 1];
  long rowSize = simdRoundUp(sequenceLength + 1) + simdLen;

  int iBeg = std::max(anchor1 - minSeparation + 1, 0);
  int iEnd = std::min(anchor1 + minSeparation - 1, profile.length);
  int jBeg = std::max(anchor2 - minSeparation + 1, 0);
  int jEnd = std::min(anchor2 + minSeparation - 1, sequenceLength);

  const char *seq = sequence + jBeg;
  const Float *Xfrom = scratch + rowSize * anchor1 + jBeg;
  const Float *Yfrom = scratch + rowSize * (profile.length + 1) + jBeg;

  for (int i = anchor1 + 1; i <= iEnd; ++i) {
    const Float *Wbackward = scratch + i * rowSize + jBeg;
    Float a = profile.values[i * profile.width + 0];
    Float b = profile.values[i * profile.width + 1];
    Float d = profile.values[i * profile.width + 2];
    Float e = profile.values[i * profile.width + 3];
    const Float *S = profile.values + i * profile.width + 4;

    Float x = 0;
    Float z = 0;
    for (int j = 0; j <= jEnd - jBeg; ++j) {
      Float y = Yfrom[j];
      Float w = x + y + z + scale;
      if (w * Wbackward[j] > wMidAnchored) return false;  // found higher score
      x = Xfrom[j];
      X[j] = S[seq[j]] * w;
      Y[j] = d * w + e * y;
      z = a * w + b * z;
    }

    Xfrom = X;
    Yfrom = Y;
  }

  Float *W = X;
  const Float *Wfrom =
    (anchor1 < profile.length) ? scratch + rowSize * (anchor1 + 1) + jBeg : Y;
  std::fill_n(Y, minSeparation * 2 - 1, 0);

  for (int i = anchor1; i >= iBeg; --i) {
    const Float *Xforward = scratch + i * rowSize + jBeg;
    Float a = profile.values[i * profile.width + 0];
    Float b = profile.values[i * profile.width + 1];
    Float d = profile.values[i * profile.width + 2];
    Float e = profile.values[i * profile.width + 3];
    const Float *S = profile.values + i * profile.width + 4;

    Float wOld = 0;
    Float z = 0;
    for (int j = jEnd - jBeg; j >= 0; --j) {
      Float y = Y[j];
      Float t = S[seq[j]];
      Float w = t * wOld + d * y + a * z + scale;
      if (i < anchor1 && w * (Xforward[j] / t) > wMidAnchored) return false;
      wOld = Wfrom[j];
      W[j] = w;
      Y[j] = w + e * y;
      z = w + b * z;
    }

    Wfrom = W;
  }

  return true;  // maybe there is no higher score nearby
}

void addMidAnchored(std::vector<AlignedSimilarity> &similarities,
		    Profile profile, const char *sequence,
		    int sequenceLength, const Float *scratch,
		    int anchor1, int anchor2,
		    Float wBegAnchored, Float wEndAnchored) {
  Float wMidAnchored = wEndAnchored * wBegAnchored;
  // this local maximum check makes it faster when there are many similarities:
  // if (!maybeLocalMaximum(profile, sequence, sequenceLength, scratch,
	// 		 anchor1, anchor2, wMidAnchored)) return;
  AlignedSimilarity s = {wMidAnchored / scale, anchor1, anchor2, wEndAnchored};
  addForwardAlignment(s.alignment, profile, sequence, sequenceLength,
		      scratch, anchor1, anchor2, wBegAnchored / 2);
  similarities.push_back(s);
}

void finishMidAnchored(AlignedSimilarity &s,
		       Profile profile, const char *sequence,
		       int sequenceLength, const Float *scratch) {
  // reverse(s.alignment.begin(), s.alignment.end());
  // addReverseAlignment(s.alignment, profile, sequence, sequenceLength,
	// 	      scratch, s.anchor1, s.anchor2, s.wEndAnchored / 2);
  // reverse(s.alignment.begin(), s.alignment.end());
}

bool isLess(const AlignedSimilarity &a, const AlignedSimilarity &b) {
  return simBeg2(a) < simBeg2(b);
}

bool isOverlapping(const std::vector<SegmentPair> &alignment1,
		   const std::vector<SegmentPair> &alignment2) {
  for (const auto &i : alignment1) {
    for (const auto &j : alignment2) {
      if (i.start1 - i.start2 == j.start1 - j.start2 &&
	  i.start1 + i.length > j.start1 && i.start1 < j.start1 + j.length)
	return true;
    }
  }
  return false;
}

void nonredundantize(std::vector<AlignedSimilarity> &similarities) {
  sort(similarities.begin(), similarities.end(), isLess);

  size_t k = 0;
  for (size_t i = 0; i < similarities.size(); ++i) {
    AlignedSimilarity &x = similarities[i];
    int end = simEnd2(x);
    for (size_t j = i + 1; j < similarities.size(); ++j) {
      AlignedSimilarity &y = similarities[j];
      if (simBeg2(y) >= end) break;
      if (isOverlapping(x.alignment, y.alignment)) {
	if (x.probRatio < y.probRatio) {
	  x.probRatio = 0;
	} else {
	  y.probRatio = 0;
	}
      }
    }
    if (x.probRatio > 0) {
      std::swap(similarities[k], similarities[i]);
      ++k;
    }
  }

  similarities.resize(k);
}

int updateInitialSimilarities(InitialSimilarity *sims, int count,
			      int anchor2, Float probRatio) {
  int i = 0;
  int j = 0;
  while (i < count && sims[i].anchor2 <= anchor2 - minSeparation) ++i;
  while (i < count && sims[i].probRatio > probRatio) sims[j++] = sims[i++];
  sims[j].probRatio = probRatio;
  sims[j].anchor2 = anchor2;
  return j + 1;
}

void setCharToNumber(char *charToNumber, const char *alphabet) {
  for (int i = 0; alphabet[i]; ++i) {
    int c = alphabet[i];
    charToNumber[toupper(c)] = charToNumber[tolower(c)] = i;
  }
}

/*
    vibe coded!!
    might not be accurate

    todo: rewrite
*/
std::unordered_map<char, std::vector<std::string>> aa2codons;
std::unordered_map<char, std::vector<std::string>> &build_standard_genetic_code() {
    // Amino acids are single-letter codes.
    // DNA codons (T not U).
    if(aa2codons.size() > 0) {
      return aa2codons;
    }

    aa2codons['A'] = {"GCT","GCC","GCA","GCG"};                 // Ala
    aa2codons['R'] = {"CGT","CGC","CGA","CGG","AGA","AGG"};     // Arg
    aa2codons['N'] = {"AAT","AAC"};                             // Asn
    aa2codons['D'] = {"GAT","GAC"};                             // Asp
    aa2codons['C'] = aa2codons['U'] = {"TGT","TGC"};            // Cys
    aa2codons['Q'] = {"CAA","CAG"};                             // Gln
    aa2codons['E'] = {"GAA","GAG"};                             // Glu
    aa2codons['G'] = {"GGT","GGC","GGA","GGG"};                 // Gly
    aa2codons['H'] = {"CAT","CAC"};                             // His
    aa2codons['I'] = {"ATT","ATC","ATA"};                       // Ile
    aa2codons['L'] = {"TTA","TTG","CTT","CTC","CTA","CTG"};     // Leu
    aa2codons['K'] = aa2codons['O'] = {"AAA","AAG"};            // Lys
    aa2codons['M'] = {"ATG"};                                   // Met
    aa2codons['F'] = {"TTT","TTC"};                             // Phe
    aa2codons['P'] = {"CCT","CCC","CCA","CCG"};                 // Pro
    aa2codons['S'] = {"TCT","TCC","TCA","TCG","AGT","AGC"};     // Ser
    aa2codons['T'] = {"ACT","ACC","ACA","ACG"};                 // Thr
    aa2codons['W'] = {"TGG"};                                   // Trp
    aa2codons['Y'] = {"TAT","TAC"};                             // Tyr
    aa2codons['V'] = {"GTT","GTC","GTA","GTG"};                 // Val
    aa2codons['*'] = {"TAA","TAG","TGA"};

    return aa2codons;
}

static void normalize(std::unordered_map<char,double>& m) {
    double s = 0.0;
    for (auto &kv : m) s += kv.second;
    if (s <= 0) return;
    for (auto &kv : m) kv.second /= s;

    // hardcode probabilities
    m['A'] = 0.25;
    m['T'] = 0.25;
    m['G'] = 0.25;
    m['C'] = 0.25;
}
struct NucDist {
    std::unordered_map<char, double> overall{{'A',0},{'C',0},{'G',0},{'T',0}};
};

NucDist infer_nucleotide_distribution_equal_synonyms(
    const std::unordered_map<char, double>& aaFreq
) {
    auto aa2codons = build_standard_genetic_code();
    NucDist out;

    // Distribute each amino acid's probability equally across its codons.
    double sm = 0;
    for (auto &kv : aaFreq) {
        char aa = (char)toupper((unsigned char)kv.first);
        auto it = aa2codons.find(aa);
        sm += kv.second;
        //std::cout << kv.first << " probs " << kv.second << std::endl;

        const std::vector<std::string>& codons = it->second;
        double perCodon = kv.second / (double)codons.size();

        for (const std::string& codon : codons) {
            char b1 = codon[0], b2 = codon[1], b3 = codon[2];

            out.overall[b1] += perCodon;
            out.overall[b2] += perCodon;
            out.overall[b3] += perCodon;
        }
    }

    //std::cout << "assert " << sm << " == 1" << std::endl;
    // At this point:
    // - overall sums to 3 (because each codon contributes 3 bases) after AA normalization,
    // - pos1/pos2/pos3 each sum to 1.
    normalize(out.overall);

    return out;
}

char translate(const char* dna, int i) {
    // Codon table (DNA codons → single-letter amino acid)
    static const std::unordered_map<std::string, char> codonTable = {
        {"TTT",'F'}, {"TTC",'F'}, {"TTA",'L'}, {"TTG",'L'},
        {"CTT",'L'}, {"CTC",'L'}, {"CTA",'L'}, {"CTG",'L'},
        {"ATT",'I'}, {"ATC",'I'}, {"ATA",'I'}, {"ATG",'M'},
        {"GTT",'V'}, {"GTC",'V'}, {"GTA",'V'}, {"GTG",'V'},

        {"TCT",'S'}, {"TCC",'S'}, {"TCA",'S'}, {"TCG",'S'},
        {"CCT",'P'}, {"CCC",'P'}, {"CCA",'P'}, {"CCG",'P'},
        {"ACT",'T'}, {"ACC",'T'}, {"ACA",'T'}, {"ACG",'T'},
        {"GCT",'A'}, {"GCC",'A'}, {"GCA",'A'}, {"GCG",'A'},

        {"TAT",'Y'}, {"TAC",'Y'}, {"TAA",'*'}, {"TAG",'*'},
        {"CAT",'H'}, {"CAC",'H'}, {"CAA",'Q'}, {"CAG",'Q'},
        {"AAT",'N'}, {"AAC",'N'}, {"AAA",'K'}, {"AAG",'K'},
        {"GAT",'D'}, {"GAC",'D'}, {"GAA",'E'}, {"GAG",'E'},

        {"TGT",'C'}, {"TGC",'C'}, {"TGA",'*'}, {"TGG",'W'},
        {"CGT",'R'}, {"CGC",'R'}, {"CGA",'R'}, {"CGG",'R'},
        {"AGT",'S'}, {"AGC",'S'}, {"AGA",'R'}, {"AGG",'R'},
        {"GGT",'G'}, {"GGC",'G'}, {"GGA",'G'}, {"GGG",'G'}
    };

    // Make sure we can read 3 characters
    // if (!dna || dna[i] == '\0' || dna[i+1] == '\0' || dna[i+2] == '\0') {
    //   std::cout << "error reading " << std::endl;
    //   return '?';
    // }

    std::string codon;
    codon += dna[i];
    codon += dna[i+1];
    codon += dna[i+2];

    // Convert to uppercase (in case input isn't)
    // for (char& c : codon)
    //     c = toupper(c);

    auto it = codonTable.find(codon);
    return (it != codonTable.end()) ? it->second : '?';
}
Float log_sum_exp(Float a, Float b) {
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  Float m = std::max(a, b);
  return m + log2(exp2(a - m) + exp2(b - m));
}
// vibe-coded section end


void findSimilarities(std::vector<AlignedSimilarity> &similarities,
		      Profile profile, const char *sequence,
		      int sequenceLength, Float *scratch,
		      Float minProbRatio) {
    //sequenceLength = std::min(sequenceLength, 5000);
    //std::cout << "findSimilarities called on " << (sequenceLength) << " x " << profile.length << std::endl;
    // need original characters to do DNA -> protein
    // TODO: remove this hack
    const char *alphabet = getAlphabet(profile.width - nonLetterWidth);
    std::string sequence_decompressed;
    bool contains_prot = false;
    for(int i = 0; i < sequenceLength; i++) {
        assert(sequence[i] < 22);
        sequence_decompressed += alphabet[sequence[i]];
        contains_prot |= (alphabet[sequence[i]] != 'A' && alphabet[sequence[i]] != 'C' && alphabet[sequence[i]] != 'G' && alphabet[sequence[i]] != 'T');
    }
    if(!contains_prot) {
        //std::cout << sequence_decompressed << std::endl;
        //exit(0);
    } 

    char charToNumber[256];
    setCharToNumber(charToNumber, alphabet);

    
    auto &dump = build_standard_genetic_code();

    std::priority_queue<AlignedSimilarity, std::vector<AlignedSimilarity>, decltype([](auto a, auto b) {
      return a.probRatio > b.probRatio;
    })> pq;

    // MEGA HACK to avoid retranslating every pHMM (unsure if actually helps)
//    if(!similarities.size()) {
      decoded = std::vector<std::pair<int, Float>>(sequenceLength, {INT_MIN, NAN});
      
      auto translate_wrapper = [&](std::string &dna, int i) -> std::pair<int, Float> {
          NucDist dist = *reinterpret_cast<NucDist*>(profile.debug);
          char translated = translate(dna.c_str(), i);

  #ifdef DUMMER_SCORING
          Float divisor = 1;
  #else
          //Float divisor = aa2codons.at(translated).size() * dist.overall[dna[i]] * dist.overall[dna[i + 1]] * dist.overall[dna[i + 2]]; // remove background probability 0.25^3, and also dedupliate emission
          //Float divisor = aa2codons.at(translated).size();
          Float divisor = 1;
  #endif
          if(translated == '*') {
            return {-1, STOP_CODON_PROB / divisor}; // for now
          }

          int toNum = charToNumber[translated];
          assert(0 <= toNum && toNum < strlen(alphabet));
          return {toNum, 1 / divisor};
      };
      for(int i = 0; i < sequenceLength - 2; i++) {
        decoded[i] = translate_wrapper(sequence_decompressed, i);
      }
  //  }
    
    std::vector<Float> dp(sequenceLength + 1), denom_multiplier_l(sequenceLength + 1), denom_multiplier_r(sequenceLength + 3, 1);
    

    // add compensate and merge denominators to actual/real denom
    for (int i = sequenceLength - 1; i >= 0; i--) {
#ifdef DUMMER_SCORING
      if(i + 3 < sequenceLength) {
        auto [emitNum, divisor] = decoded[i + 1];
        assert(profile.bg_probs[emitNum + 4] * divisor != 0);
        dp[i] = log2(profile.bg_probs[emitNum + 4] * divisor) + dp[i + 3];
      } else {
        dp[i] = 0;
      }
#else
      Float t1 = -INFINITY;
      Float t2;
      if (i + 3 < sequenceLength) {
        auto [emitNum, divisor] = decoded[i + 1];
        t1 = log2(1 - XX)
          + log2(profile.bg_probs[emitNum + 4] * divisor)
          + dp[i + 3];
        
      } else {
        t1 = log2(1 - XX);
      }
      
      t2 = log2(XX * 0.25) + dp[i + 1];
      
      
      dp[i] = log_sum_exp(t1, t2);
#endif
      denom_multiplier_r[i] = exp2(-(dp[i] - (i + 1 < sequenceLength ? dp[i + 1] : 0)));
    }

    std::vector dp_r(dp);
    //std::cout << "dp is " << xx_null << " in place of " << (-2 * sequenceLength) << std::endl;
    for (int i = 0; i < sequenceLength; i++) {
#ifdef DUMMER_SCORING
      if(i >= 3) {
        auto [emitNum, divisor] = decoded[i - 2];
        dp[i] = log2(profile.bg_probs[emitNum + 4] * divisor) + dp[i - 3];
      } else {
        dp[i] = log2(1.0 / 3.0);
      }
#else
      Float t1 = -INFINITY;
      Float t2 = -INFINITY;
      if (i - 2 >= 0) {
        auto [emitNum, divisor] = decoded[i - 2];
        t1 = log2(1 - XX)
          + log2(profile.bg_probs[emitNum + 4] * divisor)
          + (i - 3 >= 0 ? dp[i - 3] : log2(1.0 / 3.0));
      } else {
        t1 = log2(1.0 / 3.0);
      }

      if(i > 0) {
        t2 = log2(XX * 0.25) + dp[i - 1];
      }

      auto cur = log_sum_exp(t1, t2);
      auto delta = dp[sequenceLength - 1 - i] - cur; 
      dp[i] = cur;

#endif

      
      denom_multiplier_l[i] = exp2(-(dp[i] - (i - 1 >= 0 ? dp[i - 1] : 0)));
      // double anti_overflow = delta;
      // if(true) {
      //   std::cout << "i: " << i << " delta: " << cur << " - " << og_dp[i] << std::endl;
      // }

      //std::cout << dp[i] << ' ' << dp_r[i] << " = " << (dp[i] + dp_r[i]) << std::endl;
    }
#ifdef DUMMER_SCORING
    xx_null = dp[sequenceLength - 1];
#else
    xx_null = dp[sequenceLength - 1];
#endif
    //xx_null = log_sum_exp(log_sum_exp(dp[0], dp[1]), dp[2]) - log2(3);
    auto distribute1 = pow(2, -(xx_null / sequenceLength));
    auto distribute2 = distribute1 * distribute1;
    auto distribute3 = distribute2 * distribute1;
    //auto half = pow(2, (-2 * sequenceLength - xx_null) / 2.0);
    
    
    //std::cout << "redist " << distribute1 << std::endl;

    scratch_v2.resize_if_need(profile.length + 1, sequenceLength);
    auto &W = scratch_v2.W;
    auto &X = scratch_v2.X;
    auto &Y = scratch_v2.Y;
    auto &Z = scratch_v2.Z;

    Float mx = 0;
    AlignedSimilarity sel;
     
    for(int i = profile.length; i >= 0; i--) {
      const Float *params = profile.values + (i) * profile.width;
      const Float *params_emission_probabilities = params + 4;

      const Params &params_cur = profile.values_v2[i];
      const Params &params_later = profile.values_v2[i + 1];

      Float codon_emit_probs = 0;
      for(int j = sequenceLength - 1; j >= 0; j--) {
        if(j + 3 < sequenceLength) {
          auto [emitNum, divisor] = decoded[j + 1]; // upto j emitted alr
          codon_emit_probs = params_emission_probabilities[emitNum] * divisor;
          X[i][j] = codon_emit_probs * distribute3 * dp_access_safe(W[1], i + 1, j + 3);
        } else {
          X[i][j] = 0;
        }
        Y[0][i][j] = dp_access_safe(W[1], i + 1, j + 0) + params_later.epsilon_prime[0] * dp_access_safe(Y[0], i + 1, j + 0);
        Y[1][i][j] = dp_access_safe(W[1], i + 1, j + 2) + params_later.epsilon_prime[1] * 0.25 * 0.25 * distribute2 * dp_access_safe(Y[0], i + 1, j + 2);
        Y[2][i][j] = dp_access_safe(W[1], i + 1, j + 1) + params_later.epsilon_prime[2] * 0.25 * distribute1 * dp_access_safe(Y[0], i + 1, j + 1);
      
        Z[0][i][j] = dp_access_safe(W[1], i, j + 3) + params_cur.beta_prime[0] * codon_emit_probs * distribute3 * dp_access_safe(Z[0], i, j + 3);
        Z[1][i][j] = dp_access_safe(W[1], i, j + 1) + params_cur.beta_prime[1] * 0.25 * distribute1 * dp_access_safe(Z[0], i, j + 1);
        Z[2][i][j] = dp_access_safe(W[1], i, j + 2) + params_cur.beta_prime[2] * 0.25 * 0.25 * distribute2 * dp_access_safe(Z[0], i, j + 2);

        int alr_emitted = sequenceLength - 1 - j;
        auto one = exp2(-(xx_null / sequenceLength) * alr_emitted + dp_r[j]);

        if(i == 300) {
          //std::cout << j << " one is " << one << std::endl;
        }
        W[1][i][j] = dp_access_safe(X, i, j) * params_cur.enter_match_probability +
          dp_access_safe(Y[0], i, j) * params_cur.delta_prime[0] +
          dp_access_safe(Y[1], i, j) * params_cur.delta_prime[1] * 0.25 * 0.25 * distribute2 +
          dp_access_safe(Y[2], i, j) * params_cur.delta_prime[2] * 0.25 * distribute1 +
          dp_access_safe(Z[0], i, j) * params_cur.alpha_prime[0] * codon_emit_probs * distribute3 +
          dp_access_safe(Z[1], i, j) * params_cur.alpha_prime[1] * 0.25 * distribute1 +
          dp_access_safe(Z[2], i, j) * params_cur.alpha_prime[2] * 0.25 * 0.25 * distribute2 +
          one * scale;
      }
    }

    for(int i = 0; i <= profile.length; i++ ) {
      const Float *params = profile.values + (i) * profile.width;

      const Params &params_cur = profile.values_v2[i];
      const Float *params_emission_probabilities = params + 4;
      Float codon_emit_probs = 0;
      for(int j = 0; j < sequenceLength; j++) {
        std::array<Float, 4> w;
        

        for(int w_i = 1; w_i <= 3; w_i++) {
          int alr_emitted = j - w_i;
          Float one = 1;
          if(alr_emitted >= 0) {
            one = exp2(-(xx_null / sequenceLength) * (alr_emitted + 1) + dp[alr_emitted]);
          }
          w[w_i] = dp_access_safe(W[0], i, j - w_i) + (j - w_i >= -1 ? (scale * one) : 0);
        }

        if(j - 2 >= 0) {
          auto [emitNum, divisor] = decoded[j - 2];
          codon_emit_probs = params_emission_probabilities[emitNum] * divisor;
          X[i][j] = params_cur.enter_match_probability * codon_emit_probs * distribute3 * w[3];
        } else {
          X[i][j] = 0;
        }       
        
        Z[0][i][j] = params_cur.alpha_prime[0] * codon_emit_probs * distribute3 * w[3] +
                      params_cur.beta_prime[0] * codon_emit_probs * distribute3 * dp_access_safe(Z[0], i, j - 3) + 
                      params_cur.beta_prime[1] * codon_emit_probs * distribute3 * dp_access_safe(Z[1], i, j - 3) + 
                      params_cur.beta_prime[2] * codon_emit_probs * distribute3 * dp_access_safe(Z[2], i, j - 3);
        Z[1][i][j] = params_cur.alpha_prime[1] * 0.25 * distribute1 * w[1];
        Z[2][i][j] = params_cur.alpha_prime[2] * 0.25 * 0.25 * distribute2 * w[2];

        W[0][i][j] += Z[0][i][j] + Z[1][i][j] + Z[2][i][j];

        // Z[i][j] already computed at this point
        auto one = exp2(-(xx_null / sequenceLength) * (j + 1) + dp[j]);
        if(sequenceLength >= 2000 && i == 300) {
          //std::cout << j << " one is den: " << (-(xx_null / sequenceLength) * (j + 1)) << " num: " << dp[j] << std::endl;
        }
        w[0] = W[0][i][j] + scale * one;
        
        Y[0][i][j] = params_cur.delta_prime[0] * w[0] +
                      params_cur.epsilon_prime[0] * dp_access_safe(Y[0], i - 1, j) +
                      params_cur.epsilon_prime[1] * dp_access_safe(Y[1], i - 1, j) + 
                      params_cur.epsilon_prime[2] * dp_access_safe(Y[2], i - 1, j);
        Y[1][i][j] = params_cur.delta_prime[1] * 0.25 * 0.25 * distribute2 * w[2];
        Y[2][i][j] = params_cur.delta_prime[2] * 0.25 * distribute1 * w[1];

        //store at next hmm state
        if(i + 1 <= profile.length) W[0][i + 1][j] += X[i][j] + Y[0][i][j] + Y[1][i][j] + Y[2][i][j];

        
        auto wEndAnchored = W[0][i][j];
        auto wBegAnchored = W[1][i][j];

        Float wMidAnchored = wEndAnchored * wBegAnchored / scale;
        if(i == 456 - 1 && j == 19820) {
            std::cout << "debug sum of probabilities of ending at phmm idx 456 " << wEndAnchored << std::endl;
        }
        if(i == 9 - 1 && j == 18491 - 1) {
            std::cout << "debug sum of probabilities of starting at phmm idx 9 " << wBegAnchored << std::endl;
        }

        // if(i == 67) {
        //   std::cout << "debug f " << dp[j] << " " << dp_r[j] << " to " << dp[sequenceLength - 1] << " -> " << updated_denominator << std::endl;
        // }

        // assert(!std::isnan(wMidAnchored));
        AlignedSimilarity s = {wMidAnchored, i, j, wEndAnchored};

        //AlignedSimilarity s = {wBegAnchored, i, j, wEndAnchored};
        if(s.probRatio > mx) {
          mx = s.probRatio;
          sel = s;
        }
      }
      //W[reverse][i][sequenceLength] = scale; // boundary
    }

    if(minProbRatio >= 0) {
      //std::cout << "adding " << log2(mx) << std::endl;
      addMidAnchored(similarities, profile, sequence, sequenceLength, scratch, sel.anchor1, sel.anchor2, sel.probRatio * scale / sel.wEndAnchored, sel.wEndAnchored);
    } else {
        int i = sel.anchor1, j = sel.anchor2;
        AlignedSimilarity b = sel;
        //std::cout << log(sel.probRatio) << std::endl;
        b.probRatio = W[0][i][j];                                 
        similarities.push_back(b);
        b.probRatio = W[1][i][j];
        similarities.push_back(b);
        similarities.push_back(sel);
    }
    //std::cout << "findSimilarities finished with " << (sel.probRatio * scale) << std::endl;
}

int contigToSequencePos(Contig contig, size_t strandNum, int posInContig) {
  return contig.start + strandPosition(strandNum, contig.length, posInContig);
}

void findFinalSimilarities(std::vector<FinalSimilarity> &similarities,
			   Profile profile, const char *charVec,
			   size_t seqIdx, size_t maskedSeqIdx,
			   Contig contig, Float *scratch,
			   size_t profileNum, size_t strandNum,
			   Float minProbRatio) {
  const char *alphabet = getAlphabet(profile.width - nonLetterWidth);
  const char *profileSeq = charVec + profile.consensusSequenceIdx;
  const char *sequence = charVec + seqIdx;
  const char *maskedSequence = charVec + maskedSeqIdx;

  std::vector<AlignedSimilarity> sims;
  findSimilarities(sims, profile, maskedSequence, contig.length, scratch,
		   minProbRatio);

  for (const auto &x : sims) {
    int anchor2 = contigToSequencePos(contig, strandNum, x.anchor2);
    FinalSimilarity s = {x.probRatio, profileNum, strandNum,
			 x.anchor1, anchor2, x.anchor1, anchor2};
    if (!x.alignment.empty()) {
      s.start1 = x.alignment[0].start1;
      s.start2 = contigToSequencePos(contig, strandNum, x.alignment[0].start2);
      addAlignedProfile(s.alignedSequences, x.alignment, alphabet, profileSeq);
      addAlignedSequence(s.alignedSequences, x.alignment, alphabet,
			 sequence, maskedSequence);
    }
    similarities.push_back(s);
  }
}

double methodOfMomentsLambda(const double *scores, int n, double meanScore) {
  double pi = 3.1415926535897932;
  double s = 0;
  for (int i = 0; i < n; ++i) {
    s += (scores[i] - meanScore) * (scores[i] - meanScore);
  }
  double variance = s / n;  // apparently, method of moments doesn't use n-1
  return pi / sqrt(6 * variance);
}

double methodOfMomentsK(double meanScore, double lambda, double seqLength) {
  double euler = 0.57721566490153286;
  return exp(lambda * meanScore - euler) / seqLength;
}

double methodOfLmomentsLambda(const double *sortedScores, int n,
			      double meanScore) {
  double s = 0;
  for (int i = 0; i < n; ++i) {
    s += i * sortedScores[i];
  }
  double d = 0.5 * n * (n-1);  // !!! avoids int overflow
  return log(2.0) / (s / d - meanScore);
}

double shouldBe0(const double *scores, int scoreCount, double lambda) {
  double x = 0;
  double y = 0;
  double z = 0;
  for (int i = 0; i < scoreCount; ++i) {
    x += scores[i];
    y += exp(-lambda * scores[i]);
    z += scores[i] * exp(-lambda * scores[i]);
  }
  return 1 / lambda - x / scoreCount + z / y;
}

double maximumLikelihoodLambda(const double *scores, int n) {
  double lo = 1;
  double hi = 1;
  double x, y;
  do {
    lo /= 2;
    hi *= 2;
    x = shouldBe0(scores, n, lo);
    y = shouldBe0(scores, n, hi);
  } while ((x < 0 && y < 0) || (x > 0 && y > 0));
  double gap = hi - lo;
  while (1) {  // bisection method to find lambda that makes shouldBe0 = 0
    gap /= 2;
    double mid = lo + gap;
    if (mid <= lo) return lo;
    double z = shouldBe0(scores, n, mid);
    if ((x < 0 && z <= 0) || (x > 0 && z >= 0)) lo = mid;
  }
}

double maximumLikelihoodK(const double *scores, int n, double lambda,
			  double seqLength) {
  double s = 0;
  for (int i = 0; i < n; ++i) {
    s += exp(-lambda * scores[i]);
  }
  return n / (s * seqLength);
}

void methodOfMomentsGumbel(double &lambda, double &k, double &kSimple,
			   const double *scores, int n, double seqLength) {
  double meanScore = mean(scores, n);
  lambda = methodOfMomentsLambda(scores, n, meanScore);
  k = methodOfMomentsK(meanScore, lambda, seqLength);
  kSimple = methodOfMomentsK(meanScore, 1, seqLength);
}

void methodOfLmomentsGumbel(double &lambda, double &k,
			    const double *scores, int n, double seqLength) {
  double meanScore = mean(scores, n);
  lambda = methodOfLmomentsLambda(scores, n, meanScore);
  k = methodOfMomentsK(meanScore, lambda, seqLength);
}

void maximumLikelihoodGumbel(double &lambda, double &k, double &kSimple,
			     const double *scores, int n, double seqLength) {
  lambda = maximumLikelihoodLambda(scores, n);
  k = maximumLikelihoodK(scores, n, lambda, seqLength);
  kSimple = maximumLikelihoodK(scores, n, 1, seqLength);
}

void estimateGumbel(double &mmLambda, double &mmK, double &mmKsimple,
		    double &mlLambda, double &mlK, double &mlKsimple,
		    double &lmLambda, double &lmK,
		    double *scores, int n, double seqLength) {
  std::sort(scores, scores + n);
  methodOfMomentsGumbel(mmLambda, mmK, mmKsimple, scores, n, seqLength);
  maximumLikelihoodGumbel(mlLambda, mlK, mlKsimple, scores, n, seqLength);
  methodOfLmomentsGumbel(lmLambda, lmK, scores, n, seqLength);
}

void estimateK(Profile &profile, const Float *letterFreqs,
	       char *sequence, int sequenceLength, int border,
	       int numOfSequences, Float *scratch, int printVerbosity) {
  std::mt19937_64 randGen;
  int alphabetSize = profile.width - nonLetterWidth;
  std::discrete_distribution<> dist(letterFreqs, letterFreqs + alphabetSize);

  std::vector<double> scores(numOfSequences * 3);
  double *endScores = scores.data();
  double *begScores = endScores + numOfSequences;
  double *midScores = begScores + numOfSequences;

  if (printVerbosity > 1) {
    std::cout << "#trial\tend-anchored\t\tstart-anchored\t\tmid-anchored\n\
#\tprofPos\tseqPos\tscore\tprofPos\tseqPos\tscore\tprofPos\tseqPos\tscore"
	      << std::endl;
  }

  auto alphabet = getAlphabet(20);
  char charToNumber[256];
  setCharToNumber(charToNumber, alphabet);
  for (int i = 0; i < numOfSequences; ++i) {
    // should be "< sequenceLength", but kept for pseudo-random reproducibility
#if 0
    for (int j = 0; j <= sequenceLength; ++j) sequence[j] = dist(randGen);
#else
    static std::bernoulli_distribution frameshiftDist(XX);
    static const char bases[] = {'A', 'C', 'G', 'T'};
    static std::uniform_int_distribution<int> distDNA(0, 3);
    static std::uniform_int_distribution<int> distOffset(0, 2);

    int offset = distOffset(randGen);
    for(int j = 0; j < offset; j++) {
      sequence[j] = charToNumber[bases[distDNA(randGen)]];
    }
    for (int j = offset; j <= sequenceLength; j += 3) {
      bool shouldFS = frameshiftDist(randGen);
      if(shouldFS) {
        sequence[j] = charToNumber[bases[distDNA(randGen)]];
        j -= 2;
        continue;
      }

      int x = dist(randGen);
      std::discrete_distribution<> dist2(0, (int)aa2codons[alphabet[x]].size());
      auto &xx = aa2codons[alphabet[x]][dist2(randGen)];
      for(int k = 0; k < 3; k++) {
        if(j + k <= sequenceLength) {
          sequence[j + k] = charToNumber[xx[k]];
        }
      }
    }
#endif

    for (int j = 0; j < border; ++j) sequence[sequenceLength+j] = sequence[j];
    std::vector<AlignedSimilarity> sims;
    findSimilarities(sims, profile, sequence, sequenceLength + border,
		     scratch, -2);
    endScores[i] = log(sims[0].probRatio);
    begScores[i] = log(sims[1].probRatio);
    midScores[i] = log(sims[2].probRatio);
    if (printVerbosity > 1) {
      std::cout << (i+1) << "\t"
		<< sims[0].anchor1 << "\t" << sims[0].anchor2 << "\t"
		<< log2(sims[0].probRatio)+shift << "\t"
		<< sims[1].anchor1 << "\t" << sims[1].anchor2 << "\t"
		<< log2(sims[1].probRatio)+shift << "\t"
		<< sims[2].anchor1 << "\t" << sims[2].anchor2 << "\t"
		<< log2(sims[2].probRatio)+shift << std::endl;
    }
  }

  double MMendL, MMendK, MMendKsimple, MLendL, MLendK, MLendKsimple;
  double LMendL, LMendK;
  estimateGumbel(MMendL, MMendK, MMendKsimple, MLendL, MLendK, MLendKsimple,
		 LMendL, LMendK, endScores, numOfSequences, sequenceLength);

  double MMbegL, MMbegK, MMbegKsimple, MLbegL, MLbegK, MLbegKsimple;
  double LMbegL, LMbegK;
  estimateGumbel(MMbegL, MMbegK, MMbegKsimple, MLbegL, MLbegK, MLbegKsimple,
		 LMbegL, LMbegK, begScores, numOfSequences, sequenceLength);

  double MMmidL, MMmidK, MMmidKsimple, MLmidL, MLmidK, MLmidKsimple;
  double LMmidL, LMmidK;
  estimateGumbel(MMmidL, MMmidK, MMmidKsimple, MLmidL, MLmidK, MLmidKsimple,
		 LMmidL, LMmidK, midScores, numOfSequences, sequenceLength);

  double s = scale;

  if (printVerbosity > 1) {
    std::cout << "#\tend-\tstart-\tmid-anchored\n";

    std::cout << "#lamMM\t" << MMendL << "\t" << MMbegL << "\t"
	      << MMmidL << "\n"

	      << "#kMM\t" << MMendK / pow(s, MMendL) << "\t"
	      << MMbegK / pow(s, MMbegL) << "\t"
	      << MMmidK / pow(s, MMmidL) << "\n"

	      << "#kMM1\t" << MMendKsimple/scale << "\t" << MMbegKsimple/scale
	      << "\t" << MMmidKsimple/scale << "\n";

    std::cout << "#lamML\t" << MLendL << "\t" << MLbegL << "\t"
	      << MLmidL << "\n"

	      << "#kML\t" << MLendK / pow(s, MLendL) << "\t"
	      << MLbegK / pow(s, MLbegL) << "\t"
	      << MLmidK / pow(s, MLmidL) << "\n"

	      << "#kML1\t" << MLendKsimple/scale << "\t" << MLbegKsimple/scale
	      << "\t" << MLmidKsimple/scale << "\n";

    std::cout << "#lamLM\t" << LMendL << "\t" << LMbegL << "\t"
	      << LMmidL << "\n"

	      << "#kLM\t" << LMendK / pow(s, LMendL) << "\t"
	      << LMbegK / pow(s, LMbegL) << "\t"
	      << LMmidK / pow(s, LMmidL) << "\n";
  } else if (printVerbosity > 0) {
    std::cout << "# K: " << MMendKsimple/scale << " "
	      << MMbegKsimple/scale << " " << MMmidKsimple/scale << "\n";
  } else {
    std::cout << "# K: " << MMmidKsimple/scale << "\n";
  }

  static Float lambdas = 0, n = 0;
  std::cout << "Lambda: " << MMmidL << "\n";
  n++, lambdas += MMmidL;
  std::cout << "Avg Lambda: " << (lambdas / n) << "\n";

  profile.gumbelKendAnchored = MMendKsimple;
  profile.gumbelKbegAnchored = MMbegKsimple;
  profile.gumbelKmidAnchored = MMmidKsimple;
}

int intFromText(const char *text) {
  long x = strtol(text, 0, 0);
  if (x > INT_MAX || x < INT_MIN) return -1;
  return x;
}

double probFromText(const char *text) {
  if (*text == '*') return 0;
  double d = strtod(text, 0);
  return exp(-d);
}

void normalize(Float *x, int n) {
  double sum = 0;
  for (int i = 0; i < n; ++i) sum += x[i];
  assert(sum > 0);
  for (int i = 0; i < n; ++i) x[i] /= sum;
}

double meanOfLogs(const Float *x, int n) {
  double m = 1;
  for (int i = 0; i < n; ++i) m *= x[i];
  return log(m) / n;
}

double myMean(const Float *values, int length, int step, int meanType,
	      Float *valuesForMedian, const float *tantanProbs) {
  double mean = 0;
  int n = 0;
  for (int i = 0; i < length; ++i) {
    if (tantanProbs[i] >= 0.5) continue;
    double v = values[i * step];
    // Geometric mean is bad for zero (or very low) probabilities
    // All letter probs in Dfam-curated_only 3.9 and Pfam-A 38.0 are > 1e-6
    if (meanType == 'G') mean += log(std::max(v, 1e-6));  // geometric mean
    if (meanType == 'A') mean += v;                       // arithmetic mean
    if (meanType == 'M') valuesForMedian[n] = v;          // median
    ++n;
  }
  assert(n > 0);
  if (meanType == 'G') return exp(mean / n);
  if (meanType == 'A') return mean / n;
  std::sort(valuesForMedian, valuesForMedian + n);
  return valuesForMedian[n / 2];
}

void filterLetterProbabilities(Float *letterProbs, int length, int step,
			       double stdDev, bool keepNonvaryingTerm) {
  const double sqrt2pi = 2.5066282746310005;
  const double inv2var = 0.5 / (stdDev * stdDev);
  const int gaussianLimit = ceil(stdDev * 8);  // truncate Gaussian tails
  const int alphabetSize = step - nonLetterWidth;
  std::vector<double> meanLogProbs(length);
  std::vector<double> values(length);

  for (int i = 0; i < length; ++i) {
    meanLogProbs[i] = meanOfLogs(letterProbs + i * step, alphabetSize);
  }  // at each position in the profile, calculate: mean(log(letter prob))

  for (int k = 0; k < alphabetSize; ++k) {
    for (int i = 0; i < length; ++i) {
      double prob = letterProbs[i * step + k];
      values[i] = log(prob) - meanLogProbs[i];  // apply filter to this
    }

    double addItBack = keepNonvaryingTerm ? mean(values.data(), length) : 0.0;
    for (int i = 0; i < length; ++i) {
      double sum = 0;
      for (int j = -gaussianLimit; j <= gaussianLimit; ++j) {
	// xxx this treats the profile as circular (wrapping around at
	// the edges), which is rarely appropriate, but ensures no
	// change in average value:
	int x = (i + j) % length;
	if (x < 0) x += length;
	sum += values[x] * exp(-1.0 * j * j * inv2var);
      }
      sum /= stdDev * sqrt2pi;
      letterProbs[i * step + k] = exp(values[i] - sum + addItBack);
    }
  }

  for (int i = 0; i < length; ++i) {
    normalize(letterProbs + i * step, alphabetSize);
  }
}

int finalizeProfile(Profile &p, char *consensusSequence,
		    int backgroundProbsType, bool isMask,
		    double filterStdDev, bool keepNonvaryingTerm) {
  int alphabetSize = p.width - nonLetterWidth;
  std::vector<float> tantanProbs(p.length);
  std::vector<Float> valuesForMedian(p.length);
  Float *end = p.values + p.width * p.length;

  if (end[3] <= 0) {
    // set the final epsilon to the geometric mean of the other epsilons
    end[3] = myMean(p.values + p.width + 3, p.length - 1, p.width, 'G',
		    valuesForMedian.data(), tantanProbs.data());
  }

  if (filterStdDev > 0) {
    filterLetterProbabilities(p.values + 4, p.length, p.width,
			      filterStdDev, keepNonvaryingTerm);
  } else if (isMask && (alphabetSize == 4 || alphabetSize == 20)) {
    calcTantanProbabilities((const unsigned char *)consensusSequence, p.length,
			    alphabetSize > 4, tantanProbs.data());
  }

  double sumOfMeans = 0;
  for (int k = 4; k < 4 + alphabetSize; ++k) {
    double mean = myMean(p.values + k, p.length, p.width, backgroundProbsType,
			 valuesForMedian.data(), tantanProbs.data());
    end[k] = mean;
    sumOfMeans += mean;
  }

  p.bg_probs.push_back(0);
  p.bg_probs.push_back(0);
  p.bg_probs.push_back(0);
  p.bg_probs.push_back(1); // hack
  for (int k = 4; k < 4 + alphabetSize; ++k) {
    end[k] /= sumOfMeans;
    p.bg_probs.push_back((1 - STOP_CODON_PROB) * end[k]);
  }

  std::unordered_map<char, double> dist;
  char charToNumber[256];
  const char *alphabet = getAlphabet(p.width - nonLetterWidth);
  setCharToNumber(charToNumber, alphabet);
  for(auto c : std::string(alphabet)) {
    if(c == 'O') c = 'K';
    if(c == 'U') c = 'C';
    dist[c] = end[4 + charToNumber[c]] * (1 - STOP_CODON_PROB);
  }
  dist['*'] = STOP_CODON_PROB;

  auto ret = new NucDist(infer_nucleotide_distribution_equal_synonyms(dist));
  p.debug = ret;
  //std::cout << ret->overall.at('A') << ' ' << ret->overall.at('C') << ' ' << ret->overall.at('G') << ' ' << ret->overall.at('T') << '\n';
  

  for (int i = 0; ; ++i) {
    p.values_v2.push_back({0});

    Float *probs = p.values + i * p.width;
    double alpha = probs[0];
    double beta = probs[1];

    double alphaFS1 = alpha * FRAMESHIFT1_MULTIPLIER;
    double alphaFS2 = alpha * FRAMESHIFT2_MULTIPLIER;
    p.values_v2.rbegin()->alpha_prime[0] = alpha * (1 - beta);
    p.values_v2.rbegin()->alpha_prime[1] = alphaFS1 * (1 - beta);
    p.values_v2.rbegin()->alpha_prime[2] = alphaFS2 * (1 - beta);

    // all beta, are equal for now
    for(int i = 0; i <= 2; i++) {
      p.values_v2.rbegin()->beta_prime[i] = beta;
    }

    probs[0] = alpha;

    if (i == p.length) break;

    double delta = probs[2];
    double epsilon = probs[3];
    double epsilon1 = probs[p.width + 3];

    double deltaFS1 = delta * FRAMESHIFT1_MULTIPLIER;
    double deltaFS2 = delta * FRAMESHIFT2_MULTIPLIER;
    p.values_v2.rbegin()->delta_prime[0] = delta * (1 - epsilon1);
    p.values_v2.rbegin()->delta_prime[1] = deltaFS1 * (1 - epsilon1);
    p.values_v2.rbegin()->delta_prime[2] = deltaFS2 * (1 - epsilon1);

    for(int i = 0; i <= 2; i++) {
      p.values_v2.rbegin()->epsilon_prime[i] = epsilon * (1 - epsilon1) / (1 - epsilon);
    }

    p.values_v2.rbegin()->enter_match_probability = (1 - alpha - alphaFS1 - alphaFS2 - delta - deltaFS1 - deltaFS2);

    
    double c = (1 - alpha - delta);
    if (epsilon >= 1) return 0;
    probs[2] = delta;
    probs[3] = 1; // workaround -1 emitted by stop codon
    Float minVal = c;
    for (int k = 4; k < 4 + alphabetSize; ++k) {
      if (tantanProbs[i] >= 0.5) probs[k] = end[k];
      double p = probs[k];
      probs[k] = ((1 -STOP_CODON_PROB) /* minus stop codon */ * p);

      minVal = std::min(minVal, probs[k]);
    }
    if (alphabetSize == 20) {
      probs[4 + 20] = probs[4 + 1];  // selenocysteine = cysteine
      probs[4 + 21] = probs[4 + 8];  // pyrrolysine = lysine
    }
    probs[4 + alphabetSize + 2] = minVal;  // for masked sequence letters
    if (tantanProbs[i] >= 0.5) consensusSequence[i] |= 32;
  }

  // extra padding
  p.values_v2.push_back({0});

  return 1;
}

int readProfiles(std::istream &in, std::vector<Profile> &profiles,
		 std::vector<Float> &values, std::vector<char> &charVec,
		 int backgroundProbsType, bool isMask,
		 double filterStdDev, bool keepNonvaryingTerm) {
  Profile profile = {0};
  int state = 0;
  std::string line, word;
  while (getline(in, line)) {
    std::istringstream iss(line);
    iss >> word;
    switch (state) {
    case 0:
      if (word == "NAME") {
	profile.nameIdx = charVec.size();
	iss >> word;
	const char *name = word.c_str();
	charVec.insert(charVec.end(), name, name + word.size() + 1);
	profile.consensusSequenceIdx = charVec.size();
      } else if (word == "HMM") {
	++state;
      }
      break;
    case 1:
      ++state;
      break;
    case 2:
      if (word != "COMPO") ++state;
      break;
    case 3:
      {
	iss >> word;
	double MtoI = probFromText(word.c_str());
	iss >> word;
	double MtoD = probFromText(word.c_str());
	iss >> word >> word;
	double ItoI = probFromText(word.c_str());
	iss >> word >> word;
	double DtoD = probFromText(word.c_str());
	if (!iss) return 0;
	if (MtoI > 1 || MtoD > 1 || ItoI > 1 || DtoD > 1) return 0;
	values.push_back(MtoI);
	values.push_back(ItoI);
	values.push_back(MtoD);
	values.push_back(DtoD);
      }
      ++state;
      break;
    case 4:
      if (word == "//") {
	if (profile.length < 2) return 0;
	values.insert(values.end(), profile.width - 4, 0.0);
	profiles.push_back(profile);
	profile.width = profile.length = 0;
	state = 0;
      } else {
	int k = 0;
	while (iss >> word && strchr(word.c_str(), '.')) {  // xxx "*"?
	  double prob = probFromText(word.c_str());
	  if (prob > 1) return 0;
	  values.push_back(prob);
	  ++k;
	}
	values.insert(values.end(), nonLetterWidth - 4, 0.0);  // extra letters
	if (k == 0) return 0;
	if (profile.width > 0 && k + nonLetterWidth != profile.width) return 0;
	profile.width = k + nonLetterWidth;
	profile.length += 1;
	if (profile.length + 1 > INT_MAX / profile.width) return 0;
	const Float *letterProbs = &values[values.size() - profile.width + 4];
	const Float *m = std::max_element(letterProbs, letterProbs + k);
	charVec.push_back(m - letterProbs);  // consensus sequence
	state = 2;
      }
    }
  }

  Float *v = &values[0];
  for (auto &p : profiles) {
    p.values = v;
    char *consensus = &charVec[p.consensusSequenceIdx];
    if (!finalizeProfile(p, consensus, backgroundProbsType, isMask,
			 filterStdDev, keepNonvaryingTerm)) return 0;
    v += p.width * (p.length + 1);
  }

  return state == 0;
}

Float *resizeMem(Float *v, size_t &size,
		 int profileLength, int sequenceLength) {
  long rowSize = simdRoundUp(sequenceLength + 1) + simdLen;
  if (rowSize > LONG_MAX / (profileLength+2)) {
    std::cerr << "too big combination of sequence and profile\n";
    return 0;
  }
  size_t s = rowSize * (profileLength+2);
  if (s > size) {
    size = s;
    free(v);
    v = (Float *)aligned_alloc(simdLen * sizeof(Float), s * sizeof(Float));
    // this memory allocation doesn't get "free"-d at the end: that is ok!
    if (!v) std::cerr << "failed to allocate memory for " << s << " numbers\n";
  }
  return v;
}

void makeMaskedSequence(char *sequence, int length, int alphabetSize) {
  std::vector<float> tantanProbs(length);
  calcTantanProbabilities((const unsigned char *)sequence, length,
			  alphabetSize > 4, tantanProbs.data());
  int mask = alphabetSize + 2;
  for (int i = 0; i < length; ++i) {
    sequence[length + i] = (tantanProbs[i] < 0.5) ? sequence[i] : mask;
  }
}

int main(int argc, char* argv[]) {
  double evalueOpt = OPT_e;
  int strandOpt = OPT_s;
  int maskOpt = OPT_m;
  double filterStdDev = 0;
  bool keepNonvaryingTerm = false;
  int randomSeqNum = OPT_t;
  int randomSeqLen = OPT_l;
  int border = OPT_b;
  int backgroundProbsType = 'G';

  const char help[] = "\
usage: dummer profiles.hmm [sequences.fa]\n\
\n\
Find similarities between sequences and profiles.   A profile is a set of\n\
position-specific letter, deletion, and insertion probabilities.\n\
\n\
Options:\n\
  -h, --help        show this help message and exit\n\
  -V, --version     show version and exit\n\
  -v, --verbose     show progress messages\n\
  -e E, --evalue E  find similarities with E-value <= this (default: "
    STR(OPT_e) ")\n\
  -s S, --strand S  DNA strand: 0=reverse, 1=forward, 2=both (default: "
    STR(OPT_s) ")\n\
  -m M, --mask M    mask simple regions of:\n\
                    0=neither, 1=profile, 2=sequence, 3=both (default: "
    STR(OPT_m) ")\n\
\n\
Options for low-cut/high-pass filter on position-specific letter probabilities:\n\
  -d D, --dev D     standard deviation for Gaussian filter\n\
  -D D, --Dev D     same as above, but keep the non-varying component\n\
\n\
Options for random sequences:\n\
  -t T, --trials T  generate this many random sequences (default: "
    STR(OPT_t) ")\n\
  -l L, --length L  length of each random sequence (default: "
    STR(OPT_l) ")\n\
  -b B, --border B  add this size border to each random sequence (default: "
    STR(OPT_b) ")\n\
\n\
Options for background letter probabilities:\n\
  --barithmetic     arithmetic mean of position-specific probabilities\n\
  --bgeometric      geometric mean of position-specific probabilities (default)\n\
  --bmedian         median of position-specific probabilities\n\
";

  const char sOpts[] = "hVve:s:m:d:D:t:l:b:";

  static struct option lOpts[] = {
    {"help",    no_argument,       0, 'h'},
    {"version", no_argument,       0, 'V'},
    {"verbose", no_argument,       0, 'v'},
    {"evalue",  required_argument, 0, 'e'},
    {"strand",  required_argument, 0, 's'},
    {"mask",    required_argument, 0, 'm'},
    {"dev",     required_argument, 0, 'd'},
    {"Dev",     required_argument, 0, 'D'},
    {"trials",  required_argument, 0, 't'},
    {"length",  required_argument, 0, 'l'},
    {"border",  required_argument, 0, 'b'},
    {"barithmetic", no_argument,   0, 'A'},
    {"bgeometric",  no_argument,   0, 'G'},
    {"bmedian",     no_argument,   0, 'M'},
    {0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      std::cout << help;
      return 0;
    case 'V':
      std::cout << "DUMMER "
#include "version.hh"
	"\n";
      return 0;
    case 'v':
      ++verbosity;
      break;
    case 'e':
      evalueOpt = strtod(optarg, 0);
      if (evalueOpt < 0) return badOpt();
      break;
    case 's':
      strandOpt = intFromText(optarg);
      if (strandOpt < 0 || strandOpt > 2) return badOpt();
      break;
    case 'm':
      maskOpt = intFromText(optarg);
      if (maskOpt < 0 || maskOpt > 3) return badOpt();
      break;
    case 'd':
      filterStdDev = strtod(optarg, 0);
      // too low: discretized Gaussian problems; too high: overflow or slow
      if (filterStdDev < 2 || filterStdDev > 1000) return badOpt();
      break;
    case 'D':
      filterStdDev = strtod(optarg, 0);
      if (filterStdDev < 2 || filterStdDev > 1000) return badOpt();
      keepNonvaryingTerm = true;
      break;
    case 't':
      randomSeqNum = intFromText(optarg);
      if (randomSeqNum < 1) return badOpt();
      break;
    case 'l':
      randomSeqLen = intFromText(optarg);
      if (randomSeqLen < 1 || randomSeqLen > INT_MAX - 2 * simdLen)
	return badOpt();
      break;
    case 'b':
      border = intFromText(optarg);
      if (border < 0) return badOpt();
      break;
    case 'A':
      backgroundProbsType = 'A';
      break;
    case 'G':
      backgroundProbsType = 'G';
      break;
    case 'M':
      backgroundProbsType = 'M';
      break;
    case '?':
      std::cerr << help;
      return 1;
    }
  }

  if (filterStdDev > 0) maskOpt &= 2;  // filtering turns off profile-masking

  if (argc - optind < 1 || argc - optind > 2) {
    std::cerr << help;
    return 1;
  }

  if (border > INT_MAX - 2 * simdLen - randomSeqLen) {
    return err("sequence + border is too big");
  }

  std::vector<char> charVec;
  std::vector<Float> profileValues;
  std::vector<Profile> profiles;

  {
    std::ifstream file;
    std::istream &in = openFile(file, argv[optind]);
    if (!file) return 1;
    if (!readProfiles(in, profiles, profileValues, charVec,
		      backgroundProbsType, maskOpt & 1,
		      filterStdDev, keepNonvaryingTerm)) {
      return err("can't read the profile data");
    }
  }

  size_t numOfProfiles = profiles.size();

  int maxProfileLength = 0;
  for (size_t i = 0; i < numOfProfiles; ++i) {
    maxProfileLength = std::max(maxProfileLength, profiles[i].length);
  }

  size_t seqIdx = charVec.size();
  charVec.resize(seqIdx + simdRoundUp(randomSeqLen + border + 1));
  Float *scratch = 0;
  size_t scratchSize = 0;
  scratch = resizeMem(scratch, scratchSize,
		      maxProfileLength, randomSeqLen + border);
  if (!scratch) return 1;

  std::cout << "# DUMMER "
#include "version.hh"
    "\n";
  std::cout << "# Bytes per floating-point number: " << sizeof(Float) << "\n";
  if (filterStdDev > 0)
    std::cout << "# Filtering position-specific letter probabilities: std dev "
	      << filterStdDev << "\n";
  std::cout << "# Background letter probabilities: "
	    << (backgroundProbsType == 'A' ? "arithmetic mean" :
		backgroundProbsType == 'G' ? "geometric mean" : "median")
	    << " of foreground probabilities\n";
  std::cout << "# Random sequences: trials=" << randomSeqNum
	    << " length=" << randomSeqLen << " border=" << border << "\n";
  if (maskOpt & 1) std::cout << "# Masking simple regions in profiles\n";
  if (argc - optind > 1) {
    if (maskOpt & 2) std::cout << "# Masking simple regions in sequences\n";
    if (evalueOpt > 0) std::cout << "# E-value <= " << evalueOpt << "\n";
    if (strandOpt < 2)
      std::cout << "# Strand: " << (strandOpt ? "forward" : "reverse") << "\n";
  }

  int printVerbosity = (argc - optind < 2) * 2 + (evalueOpt <= 0);

  for (auto &p : profiles) {
    std::cout << "\n";
    std::cout << "# Profile name: " << &charVec[p.nameIdx] << "\n";
    std::cout << "# Profile length: " << p.length << "\n";
    if (maskOpt & 1) {
      int maskCount = 0;
      const char *consensus = &charVec[p.consensusSequenceIdx];
      for (int i = 0; i < p.length; ++i) maskCount += (consensus[i] > 31);
      std::cout << "# Positions masked by tantan: " << maskCount << "\n";
    }
    const Float *bgProbs = p.values + p.width * p.length + 4;
    std::cout << "# Background letter probabilities:";
    for (int j = 0; j < p.width - nonLetterWidth; ++j)
      std::cout << " " << bgProbs[j];
    std::cout << std::endl;

      char charToNumber[256];
      setCharToNumber(charToNumber, getAlphabet(p.width - nonLetterWidth));

      NucDist dist = *reinterpret_cast<NucDist*>(p.debug);
      Float bgProbsDNA[256] = {0};
      bgProbsDNA[charToNumber['A']] = dist.overall['A'];
      bgProbsDNA[charToNumber['C']] = dist.overall['C'];
      bgProbsDNA[charToNumber['G']] = dist.overall['G'];
      bgProbsDNA[charToNumber['T']] = dist.overall['T'];
#ifdef EVALUE
      estimateK(p, bgProbs, &charVec[seqIdx], randomSeqLen,
      border, randomSeqNum, scratch, printVerbosity);
#endif
  }

  if (argc - optind < 2 || numOfProfiles < 1) return 0;
  std::cout << std::endl;

  int width = profiles[0].width;
  for (size_t i = 1; i < numOfProfiles; ++i) {
    if (profiles[i].width != width) width = 0;
  }
  int alphabetSize = width - nonLetterWidth;
  const char *alphabet = getAlphabet(alphabetSize);
  if (!alphabet) {
    return err("the profiles should be all protein, or all nucleotide");
  }
  char charToNumber[256];
  memset(charToNumber, 127, 256);
  memset(charToNumber, 125, ' '+1);  // map "space characters" (<= ' ') to 125
  charToNumber['>'] = 126;  // record separator for FASTA-format sequences
  setCharToNumber(charToNumber, alphabet);
  if (alphabetSize == 4) setCharToNumber(charToNumber, "ACGU");  // set U = T
  if (alphabetSize != 4) strandOpt = 1;

  charVec.resize(seqIdx);
  std::vector<Sequence> sequences;
  std::vector<FinalSimilarity> similarities;
  size_t totSequenceLength = 0;

  std::ifstream file;
  std::istream &in = openFile(file, argv[optind + 1]);
  if (!file) return 1;
  Sequence sequence;
  Contig contig = {0, 0};
  while (readContig(in, sequence, contig, charVec, charToNumber)) {
    if (contig.length == 0) {
      sequences.push_back(sequence);
      continue;
    }
    seqIdx = charVec.size() - contig.length;
    size_t maskedSeqIdx = (maskOpt & 2) ? charVec.size() : seqIdx;
    // The algorithms need one arbitrary letter past the end
    // Then round up to a multiple of the SIMD length
    charVec.resize(maskedSeqIdx + simdRoundUp(contig.length + 1));
    scratch = resizeMem(scratch, scratchSize, maxProfileLength, contig.length);
    if (!scratch) return 1;
    totSequenceLength += contig.length;
    if (strandOpt == 2) totSequenceLength += contig.length;
    char *seq = &charVec[seqIdx];
    for (int s = 0; s < 2; ++s) {
      if (s != strandOpt) {
	if (maskOpt & 2) makeMaskedSequence(seq, contig.length, alphabetSize);
	size_t strandNum = sequences.size() * 2 + s;
	for (size_t j = 0; j < numOfProfiles; ++j) {
	  Profile p = profiles[j];
	  Float minProbRatio = (evalueOpt > 0) ?
	    p.gumbelKmidAnchored * totSequenceLength / evalueOpt * scale : -1;
	  if (verbosity > 1)
	    std::cerr << "Profile: " << &charVec[p.nameIdx] << "\n";
	  findFinalSimilarities(similarities, p, charVec.data(),
				seqIdx, maskedSeqIdx,
				contig, scratch, j, strandNum, minProbRatio);
	}
      }
      reverseComplement(seq, seq + contig.length);
    }
    charVec.resize(seqIdx);
  }

  std::cout << "# Total sequence length: " << totSequenceLength << "\n";

  std::cout.precision(3);
  for (size_t i = 0; i < similarities.size(); ++i) {
    Profile p = profiles[similarities[i].profileNum];
    Sequence s = sequences[similarities[i].strandNum / 2];
    double k = (evalueOpt > 0) ? p.gumbelKmidAnchored :
      (i % 3 == 0) ? p.gumbelKendAnchored :
      (i % 3 == 1) ? p.gumbelKbegAnchored : p.gumbelKmidAnchored;
    double evalue = k * totSequenceLength / similarities[i].probRatio;
    if (evalueOpt <= 0 && i % 3 == 0) std::cout << "\n";
    if (evalueOpt > 0 && evalue > evalueOpt) continue;
    printSimilarity(charVec.data(), p, s, similarities[i], evalue);
  }

  return 0;
}
