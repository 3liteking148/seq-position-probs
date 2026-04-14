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
#include <memory>
#include <queue>

#define OPT_e 100
#define OPT_s 2
#define OPT_m 3
#define OPT_t 1000
#define OPT_l 1000
#define OPT_b 100
#define OPT_x 100 // 0 to enable greedy mode

#define EVALUE
#define ALIGN

// using through BATH heuristic pipeline
#define PIPELINE_MODE

// uncomment to use codon probabilities instead of base probabilities to generate random sequences
#define ESTIMATOR_USE_RANDOM_CODONS

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
#define BACKGROUND_FRAMESHIFT_RATE (FRAMESHIFT1_MULTIPLIER + FRAMESHIFT2_MULTIPLIER)

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

  Float log2_alpha_prime[3];
  Float log2_beta_prime[3];
  Float log2_delta_prime[3];
  Float log2_epsilon_prime[3];
  Float log2_enter_match_probability;
};

struct Profile {  // position-specific (insert, delete, letter) probabilities
  Float *values;  // probabilities or probability ratios
  std::vector<Params> values_v2;
  std::vector<Float> bg_probs;
  std::vector<Float> dp, dp_r;
  Float not_align_probs;
  int width;   // number of values per position
  int length;  // number of positions
  size_t nameIdx;
  size_t consensusSequenceIdx;
  double gumbelKendAnchored, gumbelKbegAnchored, gumbelKmidAnchored, lambda;
  void *debug;
};

struct Sequence {
  size_t nameIdx;
  int length;

#ifdef PIPELINE_MODE
  size_t w_start, w_end, true_length;
  std::string target_profile;
  bool is_plus;
#endif
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

  bool operator<(const AlignedSimilarity &other) const {
    return this->probRatio < other.probRatio;
  }
  bool operator>(const AlignedSimilarity &other) const {
    return this->probRatio > other.probRatio;
  }
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
  return alphabetSize == 20 ? "ACDEFGHIKLMNPQRSTVWYUO?"  // 20 + 2 amino acids
    :    alphabetSize ==  4 ? "ACGT" : 0;
}

char complement(char c) {
  // Map DNA bases correctly using the protein alphabet indices
  // "ACDEFGHIKLMNPQRSTVWYUO?" -> A=0, C=1, G=5, T=16
  switch (c) {
    case 0:  return 16; // A -> T
    case 16: return 0;  // T -> A
    case 1:  return 5;  // C -> G
    case 5:  return 1;  // G -> C
    default: return c;  // Fallback for masked '?' or unrecognized chars
  }
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
#ifdef PIPELINE_MODE
    size_t slash = word.find('/');
    std::string chr = word.substr(0, slash);

    std::string range = word.substr(slash + 1);
    size_t dash = range.find('-');

    sequence.w_start = std::stoi(range.substr(0, dash));
    sequence.w_end = std::stoi(range.substr(dash + 1));

    std::string length, profile, strand;
    if (!(iss >> length >> profile >> strand)) return fail(in, "bad filtered sequence data: no true length, profile, or strand");

    sequence.true_length = stoll(length.substr(length.find('=') + 1));
    sequence.target_profile = profile.substr(profile.find('=') + 1);
    sequence.is_plus = strand.find("plus_strand") != std::string::npos;
    // keep sequence name
    word = chr;
#endif
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

Float not_align_probs;
void printSimilarity(const char *names, Profile p, Sequence s,
		     const FinalSimilarity &sim, double evalue) {
  if (std::isnan(evalue)) {
    return;
  }
#ifdef PIPELINE_MODE
  char strand = "+-"[!s.is_plus];
#else
  char strand = "+-"[sim.strandNum % 2];
#endif
  const char *seq = sim.alignedSequences.data();
  int length = sim.alignedSequences.size() / 2;
  int span1 = length - std::count(seq, seq + length, '-');
  int span2 = length - std::count(seq + length, seq + length * 2, '-');
  int start2 = strandPosition(sim.strandNum, s.length, sim.start2);
#ifdef PIPELINE_MODE
  if (s.is_plus) {
    start2 = s.w_start - 1 + start2;
  } else {
    start2 = s.true_length - s.w_end + start2;
  }
  int reportSeqLength = s.true_length;
#else
  int reportSeqLength = s.length;
#endif
  int anchor2 = strandPosition(sim.strandNum, reportSeqLength, sim.anchor2);
  int w1 = std::max(strlen(names + p.nameIdx), strlen(names + s.nameIdx));
  int w2 = std::max(numOfDigits(sim.start1), numOfDigits(start2));
  int w3 = std::max(numOfDigits(span1), numOfDigits(span2));
  int w4 = std::max(numOfDigits(p.length), numOfDigits(reportSeqLength));
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
	    << std::setw(w4) << reportSeqLength << " ";
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

struct metadata;
using DP_2Dv2 = std::vector<std::vector<metadata>>;
struct metadata {
  Float metric = -INFINITY;
  void *dest = nullptr;
  int dest_i, dest_j;

  metadata(Float init) : metric(init) {
    assert(!isnan(init));
  }

  metadata(DP_2Dv2 &dp, int i, int j) {
    if(0 <= i && i < dp.size() && 0 <= j && j < dp[0].size()) {
      this->metric = dp[i][j].metric;
      this->dest = &dp;
      this->dest_i = i;
      this->dest_j = j;
      return;
    }
  }

  metadata& add_cost(Float cost) {
    this->metric += cost;
    return *this;
  };

  void push_to(DP_2Dv2 &dp, int i, int j) {
    if(0 <= i && i < dp.size() && 0 <= j && j < dp[0].size()) {
      dp[i][j] = std::max(dp[i][j], *this);
    }
  }

  constexpr bool operator < (const metadata& r) const noexcept {
    return this->metric < r.metric;
  }
};


DP_2Dv2 make_dp_table_v2(size_t rows, size_t cols) {
    return DP_2Dv2(rows, std::vector<metadata>(cols, -INFINITY));
}

std::vector<std::pair<int, Float>> decoded;

template <typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U>& a,
                          const std::pair<T, U>& b)
{
    return { a.first + b.first, b.second};
}

enum class TraceState : uint8_t {
    NONE = 0, W, X, Y0, Y1, Y2, Z0, Z1, Z2
};

struct DP_Cell {
    Float metric = -INFINITY;
    TraceState backpointer = TraceState::NONE;

    inline void update(Float new_metric, TraceState state) {
        if (new_metric > metric) {
            metric = new_metric;
            backpointer = state;
        }
    }
};

struct DP_Bundle {
    DP_Cell W, X, Y0, Y1, Y2, Z0, Z1, Z2;
};

void addForwardAlignment(std::vector<SegmentPair> &alignment,
     Profile profile, const char *sequence,
     int sequenceLength, const Float *scratch,
     int iBeg, int jBeg, double half) {

  size_t rows = profile.length + 2;
  size_t cols = sequenceLength + 4;
  std::vector<DP_Bundle> dp(rows * cols);

  auto get_dp = [&](int i, int j) -> DP_Bundle& {
      return dp[i * cols + j];
  };

  auto null_prob = profile.not_align_probs;
  auto distribute1 = -(null_prob / sequenceLength);
  auto distribute2 = distribute1 * 2;
  auto distribute3 = distribute1 * 3;

  get_dp(iBeg, jBeg).W.metric = 0;

  Float cur_max_metric = -INFINITY;
  int max_dest_i = iBeg, max_dest_j = jBeg;
  TraceState max_dest_state = TraceState::W;

  int max_radius = (profile.length + 2) + (sequenceLength + 4);
  std::vector<int> min_i_per_rad(max_radius, INT_MAX);
  std::vector<int> max_i_per_rad(max_radius, -1);

  min_i_per_rad[0] = iBeg;
  max_i_per_rad[0] = iBeg;

  auto update_bounds = [&](int rad, int i_val) {
      if (rad < max_radius) {
          if (i_val < min_i_per_rad[rad]) min_i_per_rad[rad] = i_val;
          if (i_val > max_i_per_rad[rad]) max_i_per_rad[rad] = i_val;
      }
  };

  for (int radius = 0; radius < max_radius; ++radius) {
    int r_min_i = min_i_per_rad[radius];
    int r_max_i = max_i_per_rad[radius];

    if (r_min_i > r_max_i) {
        bool more_work = false;
        for (int r = radius + 1; r <= radius + 4 && r < max_radius; ++r) {
            if (min_i_per_rad[r] <= max_i_per_rad[r]) { more_work = true; break; }
        }
        if (!more_work) break;
        continue;
    }

    for (int i = r_min_i; i <= r_max_i; ++i) {
      int j = jBeg + radius - (i - iBeg);

      if (i > profile.length || j >= sequenceLength || j < 0) continue;

      DP_Bundle& cell = get_dp(i, j);
      cell.W.update(cell.X.metric, TraceState::X);
      cell.W.update(cell.Y0.metric, TraceState::Y0);
      cell.W.update(cell.Y1.metric, TraceState::Y1);
      cell.W.update(cell.Y2.metric, TraceState::Y2);
      cell.W.update(cell.Z0.metric, TraceState::Z0);
      cell.W.update(cell.Z1.metric, TraceState::Z1);
      cell.W.update(cell.Z2.metric, TraceState::Z2);

      Float w0_score = cell.W.metric;
      if (w0_score < cur_max_metric - OPT_x) continue;

      int to_emit_by_null = sequenceLength - 1 - j;
      Float one = (-(null_prob / sequenceLength) * to_emit_by_null + (profile.dp_r[j + 1]));

      if (w0_score + one > cur_max_metric) {
          cur_max_metric = w0_score + one;
          max_dest_i = i;
          max_dest_j = j;
          max_dest_state = TraceState::W;
      }

      const Params &params_cur = profile.values_v2[i];
      const Params &params_later = (i + 1 <= profile.length) ? profile.values_v2[i + 1] : params_cur;
      const Float *params_emission_probabilities = profile.values + i * profile.width + 4;

      Float codon_emit_probs = -INFINITY;
      if(j + 3 < sequenceLength) {
        auto [emitNum, divisor] = decoded[j + 1]; // upto j emitted alr
        codon_emit_probs = log2(params_emission_probabilities[emitNum] * divisor);
      }

      if (i + 1 <= profile.length && j + 3 < sequenceLength) {
          get_dp(i + 1, j + 3).X.update(w0_score + params_cur.log2_enter_match_probability + codon_emit_probs + distribute3, TraceState::W);
          update_bounds(radius + 4, i + 1);
      }

      if (i + 1 <= profile.length) {
          get_dp(i + 1, j).Y0.update(w0_score + params_cur.log2_delta_prime[0], TraceState::W);
          update_bounds(radius + 1, i + 1);

          if (j + 2 < sequenceLength) {
              get_dp(i + 1, j + 2).Y1.update(w0_score + params_cur.log2_delta_prime[1] + log2(0.25) * 2 + distribute2, TraceState::W);
              update_bounds(radius + 3, i + 1);
          }
          if (j + 1 < sequenceLength) {
              get_dp(i + 1, j + 1).Y2.update(w0_score + params_cur.log2_delta_prime[2] + log2(0.25) * 1 + distribute1, TraceState::W);
              update_bounds(radius + 2, i + 1);
          }

          get_dp(i + 1, j).Y0.update(cell.Y0.metric + params_later.log2_epsilon_prime[0], TraceState::Y0);
          get_dp(i + 1, j).Y0.update(cell.Y1.metric + params_later.log2_epsilon_prime[1], TraceState::Y1);
          get_dp(i + 1, j).Y0.update(cell.Y2.metric + params_later.log2_epsilon_prime[2], TraceState::Y2);
      }

      if (j + 3 < sequenceLength) {
          get_dp(i, j + 3).Z0.update(w0_score + params_cur.log2_alpha_prime[0] + codon_emit_probs + distribute3, TraceState::W);
          update_bounds(radius + 3, i);
      }
      if (j + 1 < sequenceLength) {
          get_dp(i, j + 1).Z1.update(w0_score + params_cur.log2_alpha_prime[1] + log2(0.25) * 1 + distribute1, TraceState::W);
          update_bounds(radius + 1, i);
      }
      if (j + 2 < sequenceLength) {
          get_dp(i, j + 2).Z2.update(w0_score + params_cur.log2_alpha_prime[2] + log2(0.25) * 2 + distribute2, TraceState::W);
          update_bounds(radius + 2, i);
      }

      if (j + 3 < sequenceLength) {
          get_dp(i, j + 3).Z0.update(cell.Z0.metric + params_cur.log2_beta_prime[0] + codon_emit_probs + distribute3, TraceState::Z0);
          get_dp(i, j + 3).Z0.update(cell.Z1.metric + params_cur.log2_beta_prime[1] + codon_emit_probs + distribute3, TraceState::Z1);
          get_dp(i, j + 3).Z0.update(cell.Z2.metric + params_cur.log2_beta_prime[2] + codon_emit_probs + distribute3, TraceState::Z2);
      }
    }
  }

  std::vector<std::pair<int, int>> path;
  int curr_i = max_dest_i;
  int curr_j = max_dest_j;
  TraceState curr_state = max_dest_state;

  while (true) {
      DP_Bundle& cell = get_dp(curr_i, curr_j);
      TraceState bp = TraceState::NONE;

      switch (curr_state) {
          case TraceState::W: bp = cell.W.backpointer; break;
          case TraceState::X:  bp = cell.X.backpointer;  break;
          case TraceState::Y0: bp = cell.Y0.backpointer; break;
          case TraceState::Y1: bp = cell.Y1.backpointer; break;
          case TraceState::Y2: bp = cell.Y2.backpointer; break;
          case TraceState::Z0: bp = cell.Z0.backpointer; break;
          case TraceState::Z1: bp = cell.Z1.backpointer; break;
          case TraceState::Z2: bp = cell.Z2.backpointer; break;
          default: break;
      }

      if (bp == TraceState::NONE) break;

      bool is_print = (bp == TraceState::X && curr_state == TraceState::W);
      if (is_print) {
          path.emplace_back(curr_i - 1, curr_j - 2);
      }

      if (curr_state == TraceState::W) {
          curr_state = bp;
      } else if (curr_state == TraceState::X) {
          curr_i -= 1; curr_j -= 3; curr_state = bp;
      } else if (curr_state == TraceState::Y0) {
          curr_i -= 1; curr_j -= 0; curr_state = bp;
      } else if (curr_state == TraceState::Y1) {
          curr_i -= 1; curr_j -= 2; curr_state = bp;
      } else if (curr_state == TraceState::Y2) {
          curr_i -= 1; curr_j -= 1; curr_state = bp;
      } else if (curr_state == TraceState::Z0) {
          curr_i -= 0; curr_j -= 3; curr_state = bp;
      } else if (curr_state == TraceState::Z1) {
          curr_i -= 0; curr_j -= 1; curr_state = bp;
      } else if (curr_state == TraceState::Z2) {
          curr_i -= 0; curr_j -= 2; curr_state = bp;
      }
  }
  std::reverse(path.begin(), path.end());
  for (const auto& p : path) {
    addForwardMatch(alignment, p.first, p.second);
  }
}

void addReverseAlignment(std::vector<SegmentPair> &alignment,
     Profile profile, const char *sequence,
     int sequenceLength, const Float *scratch,
     int iEnd, int jEnd, double half) {

  size_t rows = profile.length + 2;
  size_t cols = sequenceLength + 4;
  std::vector<DP_Bundle> dp(rows * cols);

  auto get_dp = [&](int i, int j) -> DP_Bundle& {
      return dp[i * cols + j];
  };

  auto null_prob = profile.not_align_probs;
  auto distribute1 = -(null_prob / sequenceLength);
  auto distribute2 = distribute1 * 2;
  auto distribute3 = distribute1 * 3;

  get_dp(iEnd, jEnd).W.metric = 0;

  Float cur_max_metric = -INFINITY;
  int max_dest_i = iEnd, max_dest_j = jEnd;
  TraceState max_dest_state = TraceState::W;

  int max_radius = iEnd + jEnd + 4;
  std::vector<int> min_i_per_rad(max_radius, INT_MAX);
  std::vector<int> max_i_per_rad(max_radius, -1);

  min_i_per_rad[0] = iEnd;
  max_i_per_rad[0] = iEnd;

  auto update_bounds = [&](int rad, int i_val) {
      if (rad < max_radius) {
          if (i_val < min_i_per_rad[rad]) min_i_per_rad[rad] = i_val;
          if (i_val > max_i_per_rad[rad]) max_i_per_rad[rad] = i_val;
      }
  };

  for (int radius = 0; radius < max_radius; ++radius) {
    int r_min_i = min_i_per_rad[radius];
    int r_max_i = max_i_per_rad[radius];

    if (r_min_i > r_max_i) {
        bool more_work = false;
        for (int r = radius + 1; r <= radius + 4 && r < max_radius; ++r) {
            if (min_i_per_rad[r] <= max_i_per_rad[r]) { more_work = true; break; }
        }
        if (!more_work) break;
        continue;
    }

    for (int i = r_max_i; i >= r_min_i; --i) {
      int j = jEnd - radius + (iEnd - i);
      if (i < 0 || j < 0 || j >= sequenceLength) continue;

      DP_Bundle& cell = get_dp(i, j);
      const Params &params_cur = profile.values_v2[i];
      const Float *params_emission_probabilities = profile.values + i * profile.width + 4;

      // moved to before check
      cell.W.update(cell.Y0.metric + params_cur.log2_delta_prime[0], TraceState::Y0);

      Float w_score = cell.W.metric;
      if (w_score < cur_max_metric - OPT_x) continue;

      Float codon_emit_probs = -INFINITY;
      if(j - 2 >= 0) {
        auto [emitNum, divisor] = decoded[j - 2];
        codon_emit_probs = log2(params_emission_probabilities[emitNum] * divisor);

        if (j - 3 >= 0) {
            get_dp(i, j - 3).W.update(cell.X.metric + params_cur.log2_enter_match_probability + codon_emit_probs + distribute3, TraceState::X);
            update_bounds(radius + 3, i);
        }
      }

      if (i - 1 >= 0) {
          get_dp(i - 1, j).Y0.update(cell.Y0.metric + params_cur.log2_epsilon_prime[0], TraceState::Y0);
          get_dp(i - 1, j).Y1.update(cell.Y0.metric + params_cur.log2_epsilon_prime[1], TraceState::Y0);
          get_dp(i - 1, j).Y2.update(cell.Y0.metric + params_cur.log2_epsilon_prime[2], TraceState::Y0);
          update_bounds(radius + 1, i - 1);
      }
      if (j - 2 >= 0) {
          get_dp(i, j - 2).W.update(cell.Y1.metric + params_cur.log2_delta_prime[1] + log2(0.25) * 2 + distribute2, TraceState::Y1);
          update_bounds(radius + 2, i);
      }
      if (j - 1 >= 0) {
          get_dp(i, j - 1).W.update(cell.Y2.metric + params_cur.log2_delta_prime[2] + log2(0.25) * 1 + distribute1, TraceState::Y2);
          update_bounds(radius + 1, i);
      }

      int to_emit_by_null = j + 1;
      Float one = (-(null_prob / sequenceLength) * to_emit_by_null + (profile.dp[j]));
      if (w_score + one > cur_max_metric) {
          cur_max_metric = w_score + one;
          max_dest_i = i;
          max_dest_j = j;
          max_dest_state = TraceState::W;
      }

      // from W
      if (i - 1 >= 0) {
          get_dp(i - 1, j).X.update(w_score, TraceState::W);
          get_dp(i - 1, j).Y0.update(w_score, TraceState::W);
          get_dp(i - 1, j).Y1.update(w_score, TraceState::W);
          get_dp(i - 1, j).Y2.update(w_score, TraceState::W);
          update_bounds(radius + 1, i - 1);
      }
      cell.Z0.update(w_score, TraceState::W);
      cell.Z1.update(w_score, TraceState::W);
      cell.Z2.update(w_score, TraceState::W);

      if (j - 3 >= 0) {
          get_dp(i, j - 3).Z0.update(cell.Z0.metric + params_cur.log2_beta_prime[0] + codon_emit_probs + distribute3, TraceState::Z0);
          get_dp(i, j - 3).Z1.update(cell.Z0.metric + params_cur.log2_beta_prime[1] + codon_emit_probs + distribute3, TraceState::Z0);
          get_dp(i, j - 3).Z2.update(cell.Z0.metric + params_cur.log2_beta_prime[2] + codon_emit_probs + distribute3, TraceState::Z0);
          get_dp(i, j - 3).W.update(cell.Z0.metric + params_cur.log2_alpha_prime[0] + codon_emit_probs + distribute3, TraceState::Z0);
          update_bounds(radius + 3, i);
      }
      if (j - 1 >= 0) {
          get_dp(i, j - 1).W.update(cell.Z1.metric + params_cur.log2_alpha_prime[1] + log2(0.25) * 1 + distribute1, TraceState::Z1);
          update_bounds(radius + 1, i);
      }
      if (j - 2 >= 0) {
          get_dp(i, j - 2).W.update(cell.Z2.metric + params_cur.log2_alpha_prime[2] + log2(0.25) * 2 + distribute2, TraceState::Z2);
          update_bounds(radius + 2, i);
      }
    }
  }

  std::vector<std::pair<int, int>> path;
  int curr_i = max_dest_i, curr_j = max_dest_j;
  TraceState curr_state = max_dest_state;

  while (true) {
    if (curr_i < 0 || curr_i >= rows || curr_j < 0 || curr_j >= cols) {
      break;
    }

    DP_Bundle& cell = get_dp(curr_i, curr_j);
    TraceState bp = TraceState::NONE;

    switch (curr_state) {
      case TraceState::W: bp = cell.W.backpointer; break;
      case TraceState::X:  bp = cell.X.backpointer;  break;
      case TraceState::Y0: bp = cell.Y0.backpointer; break;
      case TraceState::Y1: bp = cell.Y1.backpointer; break;
      case TraceState::Y2: bp = cell.Y2.backpointer; break;
      case TraceState::Z0: bp = cell.Z0.backpointer; break;
      case TraceState::Z1: bp = cell.Z1.backpointer; break;
      case TraceState::Z2: bp = cell.Z2.backpointer; break;
      default: break;
    }

    if (bp == TraceState::NONE) break;

    bool is_print = curr_state == TraceState::X;
    if (is_print) {
      path.emplace_back(curr_i, curr_j - 2);
    }

    if (curr_state == TraceState::W) {
      if (bp == TraceState::X || bp == TraceState::Z0) {
        curr_j += 3;
      } else if (bp == TraceState::Y1 || bp == TraceState::Z2) {
        curr_j += 2;
      } else if (bp == TraceState::Y2 || bp == TraceState::Z1) {
        curr_j += 1;
      }
    } else if (curr_state == TraceState::X || curr_state == TraceState::Y0 ||
               curr_state == TraceState::Y1 || curr_state == TraceState::Y2) {
      curr_i += 1;
    } else if (curr_state == TraceState::Z0 || curr_state == TraceState::Z1 ||
              curr_state == TraceState::Z2) {
       if (bp == TraceState::Z0) {
         curr_j += 3;
       }
    }
    curr_state = bp;
  }
  std::reverse(path.begin(), path.end());
  for (const auto& p : path) {
    addReverseMatch(alignment, p.first, p.second);
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
#ifdef ALIGN
  addForwardAlignment(s.alignment, profile, sequence, sequenceLength,
		      scratch, anchor1, anchor2, wBegAnchored / 2);
#endif
  similarities.push_back(s);
}

void finishMidAnchored(AlignedSimilarity &s,
		       Profile profile, const char *sequence,
		       int sequenceLength, const Float *scratch) {
  reverse(s.alignment.begin(), s.alignment.end());
#ifdef ALIGN
  addReverseAlignment(s.alignment, profile, sequence, sequenceLength,
		      scratch, s.anchor1, s.anchor2, s.wEndAnchored / 2);
#endif
  reverse(s.alignment.begin(), s.alignment.end());
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
    aa2codons['?'] = {"???"};  // masked

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
Float log2_sum_exp(Float a, Float b) {
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  Float m = std::max(a, b);
  return m + log2(exp2(a - m) + exp2(b - m));
}

template <typename T, bool Rolling = false>
class FlatMatrix {
  std::vector<T> data;
  size_t cols;
  size_t logical_rows;

public:
  FlatMatrix() : cols(0), logical_rows(0) {}

  void resize(size_t r, size_t c, T init = T()) {
    logical_rows = r;
    cols = c;
    if constexpr (Rolling) {
      data.resize(2 * c);
      data.assign(2 * c, init);
    } else {
      data.resize(r * c);
      data.assign(r * c, init);
    }
  }

  inline T& operator()(size_t i, size_t j) {
    // assert(0 <= i && i < logical_rows);
    // assert(0 <= j && j < cols);
    if constexpr (Rolling) {
      return data[(i & 1) * cols + j];
    } else {
      return data[i * cols + j];
    }
  }

  inline const T& operator()(size_t i, size_t j) const {
    if constexpr (Rolling) {
      return data[(i & 1) * cols + j];
    } else {
      return data[i * cols + j];
    }
  }

  // Crucial for DP: clear the current row before calculating it
  // so data from the "previous-previous" row doesn't pollute your maximums
  inline void clear_row(size_t i, T init_val = T()) {
    size_t actual_row = Rolling ? (i & 1) : i;
    auto row_start = data.begin() + (actual_row * cols);
    std::fill(row_start, row_start + cols, init_val);
  }
};

// vibe-coded section end

FlatMatrix<Float, true> W0;
FlatMatrix<Float> W1;
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
        assert(sequence[i] <= 22);
        sequence_decompressed += alphabet[sequence[i]];
        contains_prot |= (alphabet[sequence[i]] != 'A' && alphabet[sequence[i]] != 'C' && alphabet[sequence[i]] != 'G' && alphabet[sequence[i]] != 'T');
    }

    char charToNumber[256];
    setCharToNumber(charToNumber, alphabet);

    // MEGA HACK to avoid retranslating every pHMM (unsure if actually helps)
//    if(!similarities.size()) {
      decoded = std::vector<std::pair<int, Float>>(sequenceLength, {INT_MIN, NAN});

      auto translate_wrapper = [&](std::string &dna, int i) -> std::pair<int, Float> {
          NucDist dist = *reinterpret_cast<NucDist*>(profile.debug);
          char translated = translate(dna.c_str(), i);

          Float divisor = aa2codons.at(translated).size();
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

    std::vector<Float> dp(sequenceLength + 1);
    for (int i = sequenceLength - 1; i >= 0; i--) {
      Float t1 = -INFINITY;
      Float t2 = -INFINITY;
      if (i + 2 < sequenceLength) {
        auto [emitNum, divisor] = decoded[i];
        t1 = log2(1 - BACKGROUND_FRAMESHIFT_RATE)
          + log2(profile.bg_probs[emitNum + 4] * divisor)
          + dp[i + 3];
      } else {
        t1 = log2(1 - BACKGROUND_FRAMESHIFT_RATE)
          + log2(0.25) * (sequenceLength - i);
      }

      t2 = log2(BACKGROUND_FRAMESHIFT_RATE * 0.25) + dp[i + 1];
      dp[i] = log2_sum_exp(t1, t2);
    }

    std::vector dp_r(dp);
    for (int i = 0; i < sequenceLength; i++) {
      Float t1 = log2(1 - BACKGROUND_FRAMESHIFT_RATE);
      Float t2 = -INFINITY;
      if (i >= 3) {
        auto [emitNum, divisor] = decoded[i - 2];
        t1 += log2(profile.bg_probs[emitNum + 4] * divisor)
          + dp[i - 3];
      } else if (i == 2) {
        auto [emitNum, divisor] = decoded[i - 2];
        t1 += log2(profile.bg_probs[emitNum + 4] * divisor)
          + log2(1.0 / 3.0);
      } else {
        t1 += log2(0.25) * (i + 1)
          + log2(1.0 / 3.0);
      }

      t2 = log2(BACKGROUND_FRAMESHIFT_RATE * 0.25) + (i > 0 ? dp[i - 1] : log2(1.0 / 3.0));

      auto cur = log2_sum_exp(t1, t2);
      dp[i] = cur;
    }

    profile.dp = dp;
    profile.dp_r = dp_r;


    not_align_probs = dp_r[0]; // all at lower right
    for (int i = 0; i < sequenceLength; i++) {
      not_align_probs = log2_sum_exp(not_align_probs, dp[i] + dp_r[i + 1]);
    }
    profile.not_align_probs = not_align_probs;
    auto distribute1 = pow(2, -(not_align_probs / sequenceLength));
    auto distribute2 = distribute1 * distribute1;
    auto distribute3 = distribute2 * distribute1;

    W0.resize(profile.length + 1, sequenceLength + 4);
    W1.resize(profile.length + 2, sequenceLength + 4);

    std::vector<Float> Y0_next(sequenceLength + 4, 0.0);
    std::vector<Float> Y1_next(sequenceLength + 4, 0.0);
    std::vector<Float> Y2_next(sequenceLength + 4, 0.0);

    std::vector<Float> Y0_curr(sequenceLength + 4, 0.0);
    std::vector<Float> Y1_curr(sequenceLength + 4, 0.0);
    std::vector<Float> Y2_curr(sequenceLength + 4, 0.0);
    std::vector<Float> one(sequenceLength + 4, 0.0);

    for(int i = profile.length; i >= 0; i--) {
      const Params &params_cur = profile.values_v2[i];
      const Params &params_later = profile.values_v2[i + 1];
      const Float *params_emission_probabilities = profile.values + (i) * profile.width + 4;

      Float Z0_ring[4] = {0, 0, 0, 0};
      Float Z1_ring[4] = {0, 0, 0, 0};
      Float Z2_ring[4] = {0, 0, 0, 0};

      Float codon_emit_probs = 0;

      for(int j = sequenceLength - 1; j >= 0; j--) {
        // Ring buffer indices:
        // r_0 is current j. r_3 is j+3 (computed 3 iterations ago)
        int r_0 = j & 3;
        int r_1 = (j + 1) & 3;
        int r_2 = (j + 2) & 3;
        int r_3 = (j + 3) & 3;

        if(j + 3 < sequenceLength) {
          auto [emitNum, divisor] = decoded[j + 1];
          codon_emit_probs = params_emission_probabilities[emitNum] * divisor;
        }

        if (i == profile.length) {
          one[j] = exp2(-(not_align_probs / sequenceLength) * (sequenceLength - 1 - j) + dp_r[j + 1]);
        }

        Float w_val =
          W1(i + 1, j + 3) /* X[i+1][j+3] */ * params_cur.enter_match_probability * codon_emit_probs * distribute3 +
          Y0_next[j + 0] * params_cur.delta_prime[0] +
          Y1_next[j + 2] * params_cur.delta_prime[1] * 0.25 * 0.25 * distribute2 +
          Y2_next[j + 1] * params_cur.delta_prime[2] * 0.25 * distribute1 +
          Z0_ring[r_3]  * params_cur.alpha_prime[0] * codon_emit_probs * distribute3 +
          Z1_ring[r_1]  * params_cur.alpha_prime[1] * 0.25 * distribute1 +
          Z2_ring[r_2]  * params_cur.alpha_prime[2] * 0.25 * 0.25 * distribute2 +
          one[j] * scale;

        W1(i, j) = w_val;


        Y0_curr[j] = w_val + params_later.epsilon_prime[0] * Y0_next[j];
        Y1_curr[j] = w_val + params_later.epsilon_prime[1] * Y0_next[j];
        Y2_curr[j] = w_val + params_later.epsilon_prime[2] * Y0_next[j];

        Float z0_future = Z0_ring[r_3];
        Z0_ring[r_0] = w_val + params_cur.beta_prime[0] * codon_emit_probs * distribute3 * z0_future;
        Z1_ring[r_0] = w_val + params_cur.beta_prime[1] * codon_emit_probs * distribute3 * z0_future;
        Z2_ring[r_0] = w_val + params_cur.beta_prime[2] * codon_emit_probs * distribute3 * z0_future;
      }

      std::swap(Y0_curr, Y0_next);
      std::swap(Y1_curr, Y1_next);
      std::swap(Y2_curr, Y2_next);
    }

    fill(Y0_next.begin(), Y0_next.end(), 0);
    fill(Y1_next.begin(), Y1_next.end(), 0);
    fill(Y2_next.begin(), Y2_next.end(), 0);
    std::vector opt_profile_position(sequenceLength, (AlignedSimilarity){-INFINITY});
    for(int i = 0; i <= profile.length; i++ ) {
      const Float *params = profile.values + (i) * profile.width;

      const Params &params_cur = profile.values_v2[i];
      const Float *params_emission_probabilities = params + 4;
      Float codon_emit_probs = 0;

      // Ring buffer indices:
      // r_0 is current j. r_3 is j+3 (computed 3 iterations ago)
      Float Z0_ring[4] = {0, 0, 0, 0};
      Float Z1_ring[4] = {0, 0, 0, 0};
      Float Z2_ring[4] = {0, 0, 0, 0};

      for(int j = 0; j < sequenceLength; j++) {
        int last_idx_emitted_by_null = j;
        int cnt_emitted_by_null = last_idx_emitted_by_null + 1;
        if(i == 0) {
          one[j] = exp2(-(not_align_probs / sequenceLength) * cnt_emitted_by_null + dp[last_idx_emitted_by_null]);
        }

        std::array<Float, 4> w = {};
        for(int w_i = 1; w_i <= 3; w_i++) {
          if (j - w_i == -1) {
            Float one_val = (1.0 / 3.0);  // probability of having frame only, no emissions have happened yet
            w[w_i] = scale * one_val;
          } else if (j - w_i >= 0) {
            w[w_i] = W0(i, j - w_i);
          }
        }

        int r_0 = j & 3;
        int r_1 = (j - 1) & 3;
        int r_2 = (j - 2) & 3;
        int r_3 = (j - 3) & 3;

        Float X_ij = 0;
        if(j - 2 >= 0) {
          auto [emitNum, divisor] = decoded[j - 2];
          codon_emit_probs = params_emission_probabilities[emitNum] * divisor;
          X_ij = params_cur.enter_match_probability * codon_emit_probs * distribute3 * w[3];
        }
        Z0_ring[r_0] = params_cur.alpha_prime[0] * codon_emit_probs * distribute3 * w[3] +
                      params_cur.beta_prime[0] * codon_emit_probs * distribute3 * Z0_ring[r_3] +
                      params_cur.beta_prime[1] * codon_emit_probs * distribute3 * Z1_ring[r_3] +
                      params_cur.beta_prime[2] * codon_emit_probs * distribute3 * Z2_ring[r_3];
        Z1_ring[r_0] = params_cur.alpha_prime[1] * 0.25 * distribute1 * w[1];
        Z2_ring[r_0] = params_cur.alpha_prime[2] * 0.25 * 0.25 * distribute2 * w[2];

        w[0] = W0(i, j);
        w[0] += Z0_ring[r_0] + Z1_ring[r_0] + Z2_ring[r_0] + one[j] * scale;
        W0(i, j) = w[0];

        // Z[i][j] already computed at this point
        Y0_curr[j] = params_cur.delta_prime[0] * w[0] +
                      params_cur.epsilon_prime[0] * Y0_next[j] +
                      params_cur.epsilon_prime[1] * Y1_next[j] +
                      params_cur.epsilon_prime[2] * Y2_next[j];
        Y1_curr[j] = params_cur.delta_prime[1] * 0.25 * 0.25 * distribute2 * w[2];
        Y2_curr[j] = params_cur.delta_prime[2] * 0.25 * distribute1 * w[1];

        //store at next hmm state
        if(i + 1 <= profile.length) W0(i + 1, j) += X_ij + Y0_curr[j] + Y1_curr[j] + Y2_curr[j];


        auto wEndAnchored = w[0];
        auto wBegAnchored = W1(i, j);

        Float wMidAnchored = wEndAnchored * wBegAnchored / scale;
        if(i == 456 - 1 && j == 19820) {
            std::cout << "debug sum of probabilities of ending at phmm idx 456 " << wEndAnchored << std::endl;
        }
        if(i == 9 - 1 && j == 18491 - 1) {
            std::cout << "debug sum of probabilities of starting at phmm idx 9 " << wBegAnchored << std::endl;
        }

        AlignedSimilarity s = {wMidAnchored, i, j, wEndAnchored};
        opt_profile_position[j] = std::max(opt_profile_position[j], s);
      }

      W0.clear_row(i);
      std::swap(Y0_curr, Y0_next);
      std::swap(Y1_curr, Y1_next);
      std::swap(Y2_curr, Y2_next);
    }

    if(minProbRatio >= 0) {
      std::ranges::sort(opt_profile_position, std::greater<>());
      std::vector<bool>  aligned(sequenceLength);
      for (auto &aligned_similarity : opt_profile_position) {
        if (aligned_similarity.probRatio >= minProbRatio && !aligned[aligned_similarity.anchor2]) {
          std::cout << "adding " << aligned_similarity.anchor2 << " of score " << aligned_similarity.probRatio << std::endl;
          double evalue = profile.gumbelKmidAnchored * sequenceLength / pow(aligned_similarity.probRatio, profile.lambda);
          std::cout << "E: " << evalue << std::endl;
          addMidAnchored(similarities, profile, sequence, sequenceLength, scratch, aligned_similarity.anchor1, aligned_similarity.anchor2, aligned_similarity.probRatio * scale / aligned_similarity.wEndAnchored, aligned_similarity.wEndAnchored);
          auto &x = similarities.back();
          finishMidAnchored(x, profile, sequence, sequenceLength, scratch);

          // dumb heuristic
          int startIdx = std::max(simBeg2(x) - 3 * profile.length, 0);
          int endIdx   = std::min(simEnd2(x) + 3 * profile.length, sequenceLength - 1);
          std::cout << startIdx << " " << endIdx << std::endl;
          std::fill(aligned.begin() + startIdx, aligned.begin() + endIdx, true);
        }
      }
    } else {
        auto sel = *std::max_element(opt_profile_position.begin(), opt_profile_position.end());
        int i = sel.anchor1, j = sel.anchor2;
        AlignedSimilarity b = sel;
        //std::cout << log(sel.probRatio) << std::endl;
        b.probRatio = 0;
        similarities.push_back(b);
        b.probRatio = 0;
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
#ifdef ESTIMATOR_USE_RANDOM_CODONS
    static std::bernoulli_distribution frameshiftDist(BACKGROUND_FRAMESHIFT_RATE);
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
#else
    for (int j = 0; j <= sequenceLength; ++j) sequence[j] = dist(randGen);
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

  profile.gumbelKendAnchored = MMendK;
  profile.gumbelKbegAnchored = MMbegK;
  profile.gumbelKmidAnchored = MMmidK;
  profile.lambda = MMmidL;
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
  Float minVal = 1;
  for (int k = 4; k < 4 + alphabetSize; ++k) {
    end[k] /= sumOfMeans;
    p.bg_probs.push_back((1 - STOP_CODON_PROB) * end[k]);
    minVal = std::min(minVal, (1 - STOP_CODON_PROB) * end[k]);
  }
  p.bg_probs.push_back(0);
  p.bg_probs.push_back(0);
  p.bg_probs.push_back(1.0 / 64.0);

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

    for (int i = 0; i <= 2; i++) {
      p.values_v2.rbegin()->log2_alpha_prime[i] = log2(p.values_v2.rbegin()->alpha_prime[i]);
      p.values_v2.rbegin()->log2_beta_prime[i] = log2(p.values_v2.rbegin()->beta_prime[i]);
    }

    probs[0] = alpha;

    if (i == p.length) break;

    double delta = probs[2];
    double epsilon = probs[3];
    double epsilon1 = probs[p.width + 3];

    double deltaFS1 = FRAMESHIFT2_MULTIPLIER;
    double deltaFS2 = FRAMESHIFT1_MULTIPLIER;
    p.values_v2.rbegin()->delta_prime[0] = delta * (1 - epsilon1);
    p.values_v2.rbegin()->delta_prime[1] = deltaFS1 * (1 - epsilon1);
    p.values_v2.rbegin()->delta_prime[2] = deltaFS2 * (1 - epsilon1);

    for(int i = 0; i <= 2; i++) {
      p.values_v2.rbegin()->epsilon_prime[i] = epsilon * (1 - epsilon1) / (1 - epsilon);
    }

    p.values_v2.rbegin()->enter_match_probability = (1 - alpha - alphaFS1 - alphaFS2 - delta - deltaFS1 - deltaFS2);

    for (int i = 0; i <= 2; i++) {
      p.values_v2.rbegin()->log2_delta_prime[i] = log2(p.values_v2.rbegin()->delta_prime[i]);
      p.values_v2.rbegin()->log2_epsilon_prime[i] = log2(p.values_v2.rbegin()->epsilon_prime[i]);
    }
    p.values_v2.rbegin()->log2_enter_match_probability = log2(p.values_v2.rbegin()->enter_match_probability);

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
    probs[4 + alphabetSize + 2] = 1.0 / 64.0;  // for masked sequence letters
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
  std::string seq2;
  for(int i = 0; i < length; i++) {
    char val = '\0';
    switch(sequence[i]) {
      case 0:
        val = 0;
        break;
      case 1:
        val = 1;
        break;
      case 5:
        val = 2;
        break;
      case 16:
        val = 3;
        break;
      default:
        assert(0);
    }
    seq2 += val;
  }

  calcTantanProbabilities((const unsigned char *)seq2.c_str(), length,
			  false, tantanProbs.data());
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

#ifdef EVALUE
#ifdef ESTIMATOR_USE_RANDOM_CODONS
      estimateK(p, bgProbs, &charVec[seqIdx], randomSeqLen,
      border, randomSeqNum, scratch, printVerbosity);
#else
    NucDist dist = *reinterpret_cast<NucDist*>(p.debug);
    Float bgProbsDNA[256] = {0};
    bgProbsDNA[charToNumber['A']] = dist.overall['A'];
    bgProbsDNA[charToNumber['C']] = dist.overall['C'];
    bgProbsDNA[charToNumber['G']] = dist.overall['G'];
    bgProbsDNA[charToNumber['T']] = dist.overall['T'];
    estimateK(p, bgProbsDNA, &charVec[seqIdx], randomSeqLen,
      border, randomSeqNum, scratch, printVerbosity);
#endif
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
#ifdef PIPELINE_MODE
  strandOpt = 1;
#endif

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
      (std::pow(p.gumbelKmidAnchored * totSequenceLength / evalueOpt, 1.0 / p.lambda)) : -1;
	  if (verbosity > 1)
	    std::cerr << "Profile: " << &charVec[p.nameIdx] << "\n";
#ifdef PIPELINE_MODE
	  if (!strcmp(&charVec[p.nameIdx], sequence.target_profile.c_str()))
#endif
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
    double evalue = k * totSequenceLength / pow(similarities[i].probRatio, p.lambda);
    if (evalueOpt <= 0 && i % 3 == 0) std::cout << "\n";
    if (evalueOpt > 0 && evalue > evalueOpt) continue;
    printSimilarity(charVec.data(), p, s, similarities[i], evalue);
  }

  return 0;
}
