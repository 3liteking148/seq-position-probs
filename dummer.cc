// Author: Martin C. Frith 2025
// SPDX-License-Identifier: BSD-3-Clause

// See [Frith2025]: "Simple and thorough detection of related
// sequences with position-varying probabilities of substitutions,
// insertions, and deletions", MC Frith 2025

#include "dummer-util.hh"
#include "can_i_haz_simd.hh"

#include <algorithm>
#include <iomanip>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>

#define OPT_e 10.0
#define OPT_s 2
#define OPT_t 1000
#define OPT_l 500
#define OPT_b 100

#ifdef DOUBLE
typedef double Float;
typedef SimdDbl SimdFloat;
const int simdLen = simdDblLen;
#else
typedef float Float;
typedef SimdFlt SimdFloat;
const int simdLen = simdFltLen;
#endif

int simdRoundUp(int x) {  // lowest multiple of simdLen that is >= x
  return x - 1 - (x - 1) % simdLen + simdLen;
}

Float simdHorizontalMax(SimdFloat x) {  // assuming it doesn't need to be fast
  Float y[simdLen];
  simdStore(y, x);
  return *std::max_element(y, y + simdLen);
}

SimdFloat simdLookupTable(const Float *table) {  // assume 4 items in the table
  Float s[simdLen];
  for (int i = 0; i < simdLen; ++i) s[i] = table[i % 4];
  return simdLoad(s);
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

struct Profile {  // position-specific (insert, delete, letter) probabilities
  Float *values;  // probabilities or probability ratios
  int width;  // values per position: 4 + alphabetSize + 2
  int length;  // number of positions
  size_t nameIdx;
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

struct Triple {
  double endAnchored, begAnchored, midAnchored;
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
    if (c > ' ') vec.push_back(charToNumber[c]);
    c = buf->snextc();
  }

  size_t seqLen = vec.size() - seqIdx;
  if (seqLen > INT_MAX - 2 * simdLen) return fail(in, "sequence is too long!");
  contig.start = sequence.length;
  contig.length = seqLen;
  sequence.length += seqLen;
  return in;
}

int consensusLetter(Profile profile, int position) {
  const Float *S = profile.values + position * profile.width + 4;
  return std::max_element(S, S + profile.width - 6) - S;
}

void addAlignedProfile(std::vector<char> &gappedSeq,
		       const std::vector<SegmentPair> &alignment,
		       const char *alphabet, Profile profile) {
  int pos1 = alignment[0].start1;
  int pos2 = alignment[0].start2;
  for (auto a : alignment) {
    for (; pos1 < a.start1; ++pos1) {
      gappedSeq.push_back(alphabet[consensusLetter(profile, pos1)]);
    }
    gappedSeq.insert(gappedSeq.end(), a.start2 - pos2, '-');
    for (; pos1 < a.start1 + a.length; ++pos1) {
      gappedSeq.push_back(alphabet[consensusLetter(profile, pos1)]);
    }
    pos2 = a.start2 + a.length;
  }
}

void addAlignedSequence(std::vector<char> &gappedSeq,
			const std::vector<SegmentPair> &alignment,
			const char *alphabet, const char *sequence) {
  int pos1 = alignment[0].start1;
  int pos2 = alignment[0].start2;
  for (auto a : alignment) {
    gappedSeq.insert(gappedSeq.end(), a.start1 - pos1, '-');
    for (; pos2 < a.start2; ++pos2) {
      gappedSeq.push_back(alphabet[sequence[pos2]]);
    }
    for (; pos2 < a.start2 + a.length; ++pos2) {
      gappedSeq.push_back(alphabet[sequence[pos2]]);
    }
    pos1 = a.start1 + a.length;
  }
}

int strandPosition(size_t strandNum, int seqLength, int position) {
  return (strandNum % 2) ? seqLength - position : position;
}

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
  std::cout << "a score=" << log2(sim.probRatio)+shift << " E=" << evalue
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

void addForwardAlignment(std::vector<SegmentPair> &alignment,
			 Profile profile, const char *sequence,
			 int sequenceLength, const Float *scratch,
			 int iBeg, int jBeg, double half) {
  long rowSize = simdRoundUp(sequenceLength + 1) + simdLen;
  const char *seq = sequence + jBeg;

  for (int size = 16; ; size *= 2) {
    int iEnd = std::min(iBeg + size, profile.length + 1);
    int jEnd = std::min(jBeg + size, sequenceLength + 1);
    int jLen = jEnd - jBeg;
    std::vector<double> scratch2(jLen * 2);  // use "double" to avoid underflow
    double *X = scratch2.data();
    double *Y = X + jLen;
    double wSum = 0;

    for (int i = iBeg; i < iEnd; ++i) {
      const Float *Wbackward = scratch + (i+1) * rowSize + jBeg + 1;
      double a = profile.values[i * profile.width + 0];
      double b = profile.values[i * profile.width + 1];
      double d = profile.values[i * profile.width + 2];
      double e = profile.values[i * profile.width + 3];
      const Float *S = profile.values + i * profile.width + 4;

      double x = (i == iBeg) ? scale : 0;
      double z = 0;
      for (int j = 0; j < jLen; ++j) {
	double y = Y[j];
	double w = x + y + z;
	wSum += w;
	Y[j] = d * w + e * y;
	x = X[j];
	z = a * w + b * z;
	if (i == profile.length || jBeg + j == sequenceLength) continue;
	X[j] = S[seq[j]] * w;
	if (X[j] * Wbackward[j] > half * scale) {
	  addForwardMatch(alignment, i, jBeg + j);
	}
      }

      if (wSum >= half) return;
    }

    // this line should be unnecessary, except for problems such as overflow:
    if (iEnd >= profile.length && jEnd >= sequenceLength) break;
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
  if (!maybeLocalMaximum(profile, sequence, sequenceLength, scratch,
			 anchor1, anchor2, wMidAnchored)) return;
  AlignedSimilarity s = {wMidAnchored / scale, anchor1, anchor2, wEndAnchored};
  addForwardAlignment(s.alignment, profile, sequence, sequenceLength,
		      scratch, anchor1, anchor2, wBegAnchored / 2);
  similarities.push_back(s);
}

void finishMidAnchored(AlignedSimilarity &s,
		       Profile profile, const char *sequence,
		       int sequenceLength, const Float *scratch) {
  reverse(s.alignment.begin(), s.alignment.end());
  addReverseAlignment(s.alignment, profile, sequence, sequenceLength,
		      scratch, s.anchor1, s.anchor2, s.wEndAnchored / 2);
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

void findSimilarities(std::vector<AlignedSimilarity> &similarities,
		      Profile profile, const char *sequence,
		      int sequenceLength, Float *scratch,
		      Float minProbRatio) {
  const int lookupType = (simdLen > 4 && profile.width <= 10) ? 1
    :                    (simdLen > 8)                        ? 2 : 0;
  const Float zero = 0;
  SimdFloat simdScale = simdFill(scale);
  const int seqEnd = simdRoundUp(sequenceLength + 1);
  const long rowSize = seqEnd + simdLen;
  SimdFloat S1, S2;  // lookup tables for each type of letter

  // Dynamic programming initialization for forward & backward algorithms:
  for (int i = 0; i < profile.length + 2; ++i) {
    std::fill_n(scratch + i * rowSize, simdLen, 0);
  }

  scratch += simdLen;  // it's convenient to set the origin after 1st zero pad

  Float *Y = scratch + rowSize * (profile.length + 1);

  if (verbosity > 1 && minProbRatio >= -1)
    std::cerr << "Backward algorithm...\n";

  // This is algorithm S5 from [Frith2025]:
  // x       <-  S[i](Q[j]) W[i+1,j+1]
  // W[i,j]  <-  x  +  d[i] Y[i+1,j]  +  a[i] Z[i,j+1]  +  1
  // Y[i,j]  <-  W[i,j]  +  e[i] Y[i+1,j]
  // Z[i,j]  <-  W[i,j]  +  b[i] Z[i,j+1]

  // Here, it is rearranged like this:
  // u       <-  x  +  d[i] Y[i+1,j]  +  scale
  // w       <-  u  +  a[i] Z[i,j+1]
  // Z[i,j]  <-  u  +  (a[i] + b[i]) Z[i,j+1]    (using SIMD "cumulation")

  Float wMax = 0;
  int iMax, jMax;

  for (int j = 0; j < seqEnd; ++j) Y[j] = 0;

  for (int i = profile.length; i >= 0; --i) {
    Float *W = scratch + i * rowSize;
    const Float *Wfrom = (i < profile.length) ? W + rowSize + 1 : Y;
    const Float *params = profile.values + i * profile.width;
    const Float *S = params + 4;
    if (lookupType == 1) {
      S1 = simdLookupTable(S);
    } else if (lookupType == 2) {
      S1 = simdLoad(S);
      S2 = simdLoad(S + simdLen);
    }
    SimdFloat a = simdFill(params[0]);  // insert open
    SimdFloat d = simdFill(params[2]);  // delete open
    SimdFloat e = simdFill(params[3]);  // delete extend
    SimdFloat g = simdFill(params[0] + params[1]);  // insert open + extend
    SimdFloat gPowers = simdPowersRev(params[0] + params[1]);
    SimdFloat wMaxHere = simdFill(zero);

    Float scaleNow[simdLen] = {0};
    std::fill_n(scaleNow, simdLen - (seqEnd - (sequenceLength+1)), scale);
    SimdFloat simdScaleNow = simdLoad(scaleNow);

    SimdFloat z = simdFill(zero);
    for (int j = seqEnd - simdLen; j >= 0; j -= simdLen) {
      const char *seq = sequence + j;
      SimdFloat t = (lookupType < 1) ? simdLookup(S, seq)
	: (lookupType < 2) ? simdLookup(S1, seq) : simdLookup(S1, S2, seq);
      SimdFloat x = simdMul(t, simdLoad(Wfrom+j));
      SimdFloat y = simdLoad(Y+j);
      SimdFloat u = simdAdd(simdAdd(x, simdMul(d, y)), simdScaleNow);
      SimdFloat s = simdCumulateRev(u, g);
      s = simdAdd(s, simdMul(gPowers, z));  // now s has the Z[i,j] values
      z = simdShiftRev(z, s);               // shift to get the Z[i,j+1] values
      SimdFloat w = simdAdd(u, simdMul(a, z));
      wMaxHere = simdMax(wMaxHere, w);
      simdStore(W+j, w);
      simdStore(Y+j, simdAdd(w, simdMul(e, y)));
      z = simdLowItem(s);
      simdScaleNow = simdScale;
    }

    Float theMaxHere = simdHorizontalMax(wMaxHere);
    if (theMaxHere > wMax) {
      wMax = theMaxHere;
      iMax = i;
    }
  }

  if (wMax > DBL_MAX) {
    std::cerr << "numbers overflowed to infinity: giving up this comparison\n";
    return;
  }

  const Float *maxRow = scratch + iMax * rowSize;
  jMax = std::find(maxRow, maxRow + seqEnd, wMax) - maxRow;

  AlignedSimilarity begAnchored = {wMax, iMax, jMax};
  addForwardAlignment(begAnchored.alignment, profile, sequence,
		      sequenceLength, scratch, iMax, jMax, wMax / 2);

  if (verbosity > 1 && minProbRatio >= -1)
    std::cerr << "Forward algorithm...\n";

  // This is algorithm 3 from [Frith2025]:
  // w       <-  X[i-1,j-1]  +  Y[i-1,j]  +  Z[i,j-1]  +  1
  // X[i,j]  <-  S[i](Q[j]) w
  // Y[i,j]  <-  d[i] w  + e[i] Y[i-1,j]
  // Z[i,j]  <-  a[i] w  + b[i] Z[i,j-1]

  // Here, it is rearranged like this:
  // u       <-  X[i-1,j-1]  +  Y[i-1,j]  +  scale
  // w       <-  u  +  a[i] Z[i,j-1]
  // Z[i,j]  <-  u  +  (a[i] + b[i]) Z[i,j-1]    (using SIMD "cumulation")

  wMax = 0;
  Float wMidMax = 0;
  Float wBegAnchored, wEndAnchored;
  int iMidMax, jMidMax;
  std::vector<SegmentPair> alignment;
  SimdFloat simdMinProbRatio = simdFill(minProbRatio);

  for (int j = 0; j < seqEnd; ++j) Y[j] = 0;

  for (int i = 0; i <= profile.length; ++i) {
    Float wMidMaxHere = wMidMax;
    Float *X = scratch + i * rowSize;
    const Float *Xfrom = (i > 0) ? X - rowSize - 1 : Y;
    const Float *Wbackward = X;  // the forward Xs overwrite the backward Ws
    const Float *params = profile.values + i * profile.width;
    const Float *S = params + 4;
    if (lookupType == 1) {
      S1 = simdLookupTable(S);
    } else if (lookupType == 2) {
      S1 = simdLoad(S);
      S2 = simdLoad(S + simdLen);
    }
    SimdFloat a = simdFill(params[0]);  // insert open
    SimdFloat d = simdFill(params[2]);  // delete open
    SimdFloat e = simdFill(params[3]);  // delete extend
    SimdFloat g = simdFill(params[0] + params[1]);  // insert open + extend
    SimdFloat gPowers = simdPowersFwd(params[0] + params[1]);

    InitialSimilarity hits[minSeparation];
    int hitCount = 0;
    int jOld = INT_MAX;

    SimdFloat z = simdFill(zero);
    for (int j = 0; j < seqEnd; j += simdLen) {
      SimdFloat y = simdLoad(Y + j);
      SimdFloat u = simdAdd(simdAdd(simdLoad(Xfrom + j), y), simdScale);
      SimdFloat s = simdCumulateFwd(u, g);
      s = simdAdd(s, simdMul(gPowers, z));  // now s has the Z[i,j] values
      z = simdShiftFwd(z, s);               // shift to get the Z[i,j-1] values
      SimdFloat w = simdAdd(u, simdMul(a, z));
      SimdFloat wBck = simdLoad(Wbackward + j);
      SimdFloat wMid = simdMul(w, wBck);
      const char *seq = sequence + j;
      SimdFloat t = (lookupType < 1) ? simdLookup(S, seq)
	: (lookupType < 2) ? simdLookup(S1, seq) : simdLookup(S1, S2, seq);
      simdStore(X+j, simdMul(t, w));
      simdStore(Y+j, simdAdd(simdMul(d, w), simdMul(e, y)));
      z = simdHighItem(s);

      Float ws[simdLen];
      Float wBcks[simdLen];
      Float wMids[simdLen];
      if (minProbRatio >= 0) {
	if (simdGe(wMid, simdMinProbRatio)) {
	  simdStore(ws, w);
	  simdStore(wBcks, wBck);
	  simdStore(wMids, wMid);
	  for (int k = 0; k < simdLen && k <= sequenceLength - j; ++k) {
	    if (wMids[k] >= minProbRatio) {
	      if (j+k - jOld >= minSeparation) {
		addMidAnchored(similarities, profile, sequence, sequenceLength,
			       scratch, i, jOld, wBegAnchored, wEndAnchored);
		jOld = INT_MAX;
	      }
	      hitCount = updateInitialSimilarities(hits, hitCount,
						   j+k, wMids[k]);
	      if (hitCount == 1) {
		jOld = j+k;
		wBegAnchored = wBcks[k];
		wEndAnchored = ws[k];
	      }
	    }
	  }
	}
      } else {
	simdStore(ws, w);
	simdStore(wBcks, wBck);
	simdStore(wMids, wMid);
	for (int k = 0; k < simdLen && k <= sequenceLength - j; ++k) {
	  if (ws[k] > wMax) {
	    wMax = ws[k];
	    iMax = i;
	    jMax = j+k;
	  }
	  if (wMids[k] > wMidMaxHere) {
	    wMidMaxHere = wMids[k];
	    wBegAnchored = wBcks[k];
	    wEndAnchored = ws[k];
	    jMidMax = j+k;
	  }
	}
      }
    }

    if (wMidMaxHere > wMidMax) {
      wMidMax = wMidMaxHere;
      iMidMax = i;
      if (minProbRatio >= -1) {
	alignment.clear();
	addForwardAlignment(alignment, profile, sequence, sequenceLength,
			    scratch, iMidMax, jMidMax, wBegAnchored / 2);
      }
    }

    if (jOld <= sequenceLength) {
      addMidAnchored(similarities, profile, sequence, sequenceLength,
		     scratch, i, jOld, wBegAnchored, wEndAnchored);
    }
  }

  if (minProbRatio >= 0) {
    if (verbosity > 1) std::cerr << "Initial similarities: "
				 << similarities.size() << "\n";
    nonredundantize(similarities);
    if (verbosity > 1) std::cerr << "Non-overlapping forward extensions: "
				 << similarities.size() << "\n";
    for (auto &x : similarities) {
      finishMidAnchored(x, profile, sequence, sequenceLength, scratch);
    }
    nonredundantize(similarities);
    if (verbosity > 1) std::cerr << "Non-overlapping similarities: "
				 << similarities.size() << "\n";
  } else {
    AlignedSimilarity midAnchored = {wMidMax / scale, iMidMax, jMidMax,
				     wEndAnchored, alignment};
    finishMidAnchored(midAnchored, profile, sequence, sequenceLength, scratch);

    AlignedSimilarity endAnchored = {wMax, iMax, jMax};
    addReverseAlignment(endAnchored.alignment, profile, sequence,
			sequenceLength, scratch, iMax, jMax, wMax / 2);
    reverse(endAnchored.alignment.begin(), endAnchored.alignment.end());

    similarities.push_back(endAnchored);
    similarities.push_back(begAnchored);
    similarities.push_back(midAnchored);
  }
}

int contigToSequencePos(Contig contig, size_t strandNum, int posInContig) {
  return contig.start + strandPosition(strandNum, contig.length, posInContig);
}

void findFinalSimilarities(std::vector<FinalSimilarity> &similarities,
			   Profile profile, const char *sequence,
			   Contig contig, Float *scratch,
			   size_t profileNum, size_t strandNum,
			   Float minProbRatio) {
  const char *alphabet =
    (profile.width == 10) ? "acgt" : "ACDEFGHIKLMNPQRSTVWYUO";

  std::vector<AlignedSimilarity> sims;
  findSimilarities(sims, profile, sequence, contig.length, scratch,
		   minProbRatio);

  for (const auto &x : sims) {
    int anchor2 = contigToSequencePos(contig, strandNum, x.anchor2);
    FinalSimilarity s = {x.probRatio, profileNum, strandNum,
			 x.anchor1, anchor2, x.anchor1, anchor2};
    if (!x.alignment.empty()) {
      s.start1 = x.alignment[0].start1;
      s.start2 = contigToSequencePos(contig, strandNum, x.alignment[0].start2);
      addAlignedProfile(s.alignedSequences, x.alignment, alphabet, profile);
      addAlignedSequence(s.alignedSequences, x.alignment, alphabet, sequence);
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

Triple estimateK(Profile profile, const Float *letterFreqs,
		 char *sequence, int sequenceLength, int border,
		 int numOfSequences, Float *scratch, int printVerbosity) {
  std::mt19937_64 randGen;
  int alphabetSize = profile.width - 6;
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

  for (int i = 0; i < numOfSequences; ++i) {
    // should be "< sequenceLength", but kept for pseudo-random reproducibility
    for (int j = 0; j <= sequenceLength; ++j) sequence[j] = dist(randGen);
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

  Triple h = {MMendKsimple, MMbegKsimple, MMmidKsimple};
  return h;
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

double geometricMean(const Float *values, int length, int step) {
  double mean = 0;
  for (int i = 0; i < length; ++i) {
    double v = values[i * step];
    if (v <= 0) std::cerr << "zero probability at position " << i << "\n";
    if (v <= 0) return 0;
    mean += log(v);
  }
  return exp(mean / length);
}

int finalizeProfile(Profile p) {
  Float *end = p.values + p.width * p.length;

  if (end[3] <= 0) {
    // set the final epsilon to the geometric mean of the other epsilons
    end[3] = geometricMean(p.values + p.width + 3, p.length - 1, p.width);
  }

  // set the background letter probabilities proportional to the
  // geometric mean of the foreground letter probabilities
  double sum = 0;
  for (int k = 4; k < p.width - 2; ++k) {
    double m = geometricMean(p.values + k, p.length, p.width);
    if (m <= 0) return 0;
    end[k] = m;
    sum += m;
  }
  for (int k = 4; k < p.width - 2; ++k) end[k] /= sum;

  for (int i = 0; ; ++i) {
    Float *probs = p.values + i * p.width;
    double alpha = probs[0];
    double beta = probs[1];
    probs[0] = alpha * (1 - beta);

    if (i == p.length) break;

    double delta = probs[2];
    double epsilon = probs[3];
    double epsilon1 = probs[p.width + 3];
    double c = (1 - alpha - delta);
    if (epsilon >= 1) return 0;
    probs[2] = delta * (1 - epsilon1);
    probs[3] = epsilon * (1 - epsilon1) / (1 - epsilon);
    for (int k = 4; k < p.width - 2; ++k) {
      probs[k] = c * probs[k] / end[k];
    }
    if (p.width == 26) {
      probs[4 + 20] = probs[4 + 1];  // selenocysteine = cysteine
      probs[4 + 21] = probs[4 + 8];  // pyrrolysine = lysine
    }
  }

  return 1;
}

int readProfiles(std::istream &in, std::vector<Profile> &profiles,
		 std::vector<Float> &values, std::vector<char> &names) {
  Profile profile = {0};
  int state = 0;
  std::string line, word;
  while (getline(in, line)) {
    std::istringstream iss(line);
    iss >> word;
    switch (state) {
    case 0:
      if (word == "NAME") {
	profile.nameIdx = names.size();
	iss >> word;
	const char *name = word.c_str();
	names.insert(names.end(), name, name + word.size() + 1);
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
	int k = 6;
	while (iss >> word && strchr(word.c_str(), '.')) {  // xxx "*"?
	  double prob = probFromText(word.c_str());
	  if (prob > 1) return 0;
	  values.push_back(prob);
	  ++k;
	}
	values.insert(values.end(), 2, 0.0);  // allow for 2 unusual letters
	if (k == 6) return 0;
	if (profile.width > 0 && k != profile.width) return 0;
	profile.width = k;
	profile.length += 1;
	if (profile.length + 1 > INT_MAX / profile.width) return 0;
	state = 2;
      }
    }
  }

  Float *v = &values[0];
  for (auto &p : profiles) {
    p.values = v;
    if (!finalizeProfile(p)) return 0;
    v += p.width * (p.length + 1);
  }

  return state == 0;
}

void setCharToNumber(char *charToNumber, const char *alphabet) {
  for (int i = 0; alphabet[i]; ++i) {
    int c = alphabet[i];
    charToNumber[toupper(c)] = charToNumber[tolower(c)] = i;
  }
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
  }
  return v;
}

int main(int argc, char* argv[]) {
  double evalueOpt = OPT_e;
  int strandOpt = OPT_s;
  int randomSeqNum = OPT_t;
  int randomSeqLen = OPT_l;
  int border = OPT_b;

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
\n\
Options for real sequences:\n\
  -e E, --evalue E  find similarities with E-value <= this (default: "
    STR(OPT_e) ")\n\
  -s S, --strand S  DNA strand: 0=reverse, 1=forward, 2=both (default: "
    STR(OPT_s) ")\n\
\n\
Options for random sequences:\n\
  -t T, --trials T  generate this many random sequences (default: "
    STR(OPT_t) ")\n\
  -l L, --length L  length of each random sequence (default: "
    STR(OPT_l) ")\n\
  -b B, --border B  add this size border to each random sequence (default: "
    STR(OPT_b) ")\n\
";

  const char sOpts[] = "hVve:s:t:l:b:";

  static struct option lOpts[] = {
    {"help",    no_argument,       0, 'h'},
    {"version", no_argument,       0, 'V'},
    {"verbose", no_argument,       0, 'v'},
    {"evalue",  required_argument, 0, 'e'},
    {"strand",  required_argument, 0, 's'},
    {"trials",  required_argument, 0, 't'},
    {"length",  required_argument, 0, 'l'},
    {"border",  required_argument, 0, 'b'},
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
    case '?':
      std::cerr << help;
      return 1;
    }
  }

  if (argc - optind < 1 || argc - optind > 2) {
    std::cerr << help;
    return 1;
  }

  if (border > INT_MAX - 2 * simdLen - randomSeqLen) {
    std::cerr << "sequence + border is too big\n";
    return 1;
  }

  std::vector<char> charVec;
  std::vector<Float> profileValues;
  std::vector<Profile> profiles;

  {
    std::ifstream file;
    std::istream &in = openFile(file, argv[optind]);
    if (!file) return 1;
    if (!readProfiles(in, profiles, profileValues, charVec)) {
      std::cerr << "can't read the profile data\n";
      return 1;
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
  std::cout << "# Random sequences: trials=" << randomSeqNum
	    << " length=" << randomSeqLen << " border=" << border << "\n";
  if (argc - optind > 1) {
    if (evalueOpt > 0) std::cout << "# E-value <= " << evalueOpt << "\n";
    if (strandOpt < 2)
      std::cout << "# Strand: " << (strandOpt ? "forward" : "reverse") << "\n";
  }

  int printVerbosity = (argc - optind < 2) * 2 + (evalueOpt <= 0);

  double totEndK = 0;
  double totBegK = 0;
  double totMidK = 0;

  for (auto p : profiles) {
    const Float *profileEnd = p.values + p.width * p.length;
    std::cout << "\n";
    std::cout << "# Profile name: " << &charVec[p.nameIdx] << "\n";
    std::cout << "# Profile length: " << p.length << "\n";
    std::cout << "# Background letter probabilities:";
    for (int j = 4; j < p.width - 2; ++j) std::cout << " " << profileEnd[j];
    std::cout << std::endl;
    Triple r = estimateK(p, profileEnd+4, &charVec[seqIdx], randomSeqLen,
			 border, randomSeqNum, scratch, printVerbosity);
    totEndK += r.endAnchored;
    totBegK += r.begAnchored;
    totMidK += r.midAnchored;
  }

  if (argc - optind < 2 || numOfProfiles < 1) return 0;
  std::cout << std::endl;

  int width = profiles[0].width;
  for (size_t i = 1; i < numOfProfiles; ++i) {
    if (profiles[i].width != width) width = 0;
  }
  char charToNumber[256];
  memset(charToNumber, 127, 256);
  if (width == 10) {
    setCharToNumber(charToNumber, "ACGT");
    setCharToNumber(charToNumber, "ACGU");
  } else if (width == 26) {
    setCharToNumber(charToNumber, "ACDEFGHIKLMNPQRSTVWYUO");
    strandOpt = 1;
  } else {
    std::cerr << "the profiles should be all protein, or all nucleotide\n";
    return 1;
  }
  memset(charToNumber, 125, ' '+1);  // map "space characters" (<= ' ') to 125
  charToNumber['>'] = 126;

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
    // The algorithms need one arbitrary letter past the end
    // Then round up to a multiple of the SIMD length
    charVec.resize(seqIdx + simdRoundUp(contig.length + 1));
    scratch = resizeMem(scratch, scratchSize, maxProfileLength, contig.length);
    if (!scratch) return 1;
    totSequenceLength += contig.length;
    if (strandOpt == 2) totSequenceLength += contig.length;
    Float minProbRatio =
      (evalueOpt > 0) ? totMidK * totSequenceLength / evalueOpt * scale : -1;
    for (int s = 0; s < 2; ++s) {
      if (s != strandOpt) {
	size_t strandNum = sequences.size() * 2 + s;
	for (size_t j = 0; j < numOfProfiles; ++j) {
	  if (verbosity > 1)
	    std::cerr << "Profile: " << &charVec[profiles[j].nameIdx] << "\n";
	  findFinalSimilarities(similarities, profiles[j], &charVec[seqIdx],
				contig, scratch, j, strandNum, minProbRatio);
	}
      }
      reverseComplement(&charVec[seqIdx], &charVec[seqIdx] + contig.length);
    }
    charVec.resize(seqIdx);
  }

  std::cout << "# Total sequence length: " << totSequenceLength << "\n";

  std::cout.precision(3);
  for (size_t i = 0; i < similarities.size(); ++i) {
    Profile p = profiles[similarities[i].profileNum];
    Sequence s = sequences[similarities[i].strandNum / 2];
    double k = (evalueOpt > 0) ? totMidK :
      (i % 3 == 0) ? totEndK : (i % 3 == 1) ? totBegK : totMidK;
    double evalue = k * totSequenceLength / similarities[i].probRatio;
    if (evalueOpt <= 0 && i % 3 == 0) std::cout << "\n";
    if (evalueOpt > 0 && evalue > evalueOpt) continue;
    printSimilarity(charVec.data(), p, s, similarities[i], evalue);
  }

  return 0;
}
