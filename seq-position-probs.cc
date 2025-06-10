// Author: Martin C. Frith 2025

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>

typedef double Float;

// down-scale probabilities by this amount, to delay overflow:
const Float scale = 1.0 / (1<<30) / (1<<30) / (1<<3); // sqrt[min normal float]
const int shift = 63;  // add this to scores, to undo the scaling

struct Profile {  // position-specific (insert, delete, letter) probabilities
  Float *values;  // probabilities or probability ratios
  int width;  // values per position: 4 + alphabetSize + 1
  int length;  // number of positions
  size_t nameIdx;
};

struct Sequence {
  size_t nameIdx;
  int length;
};

struct SegmentPair {
  int start1, start2, length;
};

struct Triple {
  double endAnchored;
  double begAnchored;
  double midAnchored;
};

struct RawSimilarity {
  double probRatio;
  int anchor1, anchor2;
  std::vector<SegmentPair> alignment;
};

struct Similarity {
  double probRatio;
  size_t profileNum;
  size_t strandNum;
  int anchor1, anchor2;
  int start1, start2;
  std::vector<char> alignedSequences;
};

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
  return (c > 3) ? c : 3 - c;
}

void reverseComplement(char *beg, char *end) {
  while (beg < end) {
    char c = *--end;
    *end = complement(*beg);
    *beg++ = complement(c);
  }
}

bool isDash(const char *text) {
  return text[0] == '-' && text[1] == 0;
}

std::istream &openFile(std::ifstream &file, const char *name) {
  if (isDash(name)) return std::cin;
  file.open(name);
  if (!file) std::cerr << "can't open file: " << name << "\n";
  return file;
}

std::istream &readSequence(std::istream &in, Sequence &sequence,
			   std::vector<char> &vec, const char *charToNumber) {
  sequence.nameIdx = vec.size();

  std::string line, word;
  while (getline(in, line)) {
    std::istringstream iss(line);
    char c;
    if (iss >> c) {
      if (c != '>' || !(iss >> word)) in.setstate(std::ios::failbit);
      const char *name = word.c_str();
      vec.insert(vec.end(), name, name + word.size() + 1);
      break;
    }
  }

  std::streambuf *buf = in.rdbuf();
  int c = buf->sgetc();
  while (c != std::streambuf::traits_type::eof() && c != '>') {
    if (c > ' ') vec.push_back(charToNumber[c]);
    c = buf->snextc();
  }

  sequence.length = vec.size() - sequence.nameIdx - word.size() - 1;
  vec.push_back(0);  // the algorithms need one arbitrary letter past the end
  return in;
}

int consensusLetter(Profile profile, int position) {
  const Float *S = profile.values + position * profile.width + 4;
  return std::max_element(S, S + profile.width - 5) - S;
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

void printSimilarity(const char *names, Profile p, Sequence s,
		     const Similarity &sim, double evalue) {
  char strand = "+-"[sim.strandNum % 2];
  const char *seq = sim.alignedSequences.data();
  int length = sim.alignedSequences.size() / 2;
  int span1 = length - std::count(seq, seq + length, '-');
  int span2 = length - std::count(seq + length, seq + length * 2, '-');
  int w1 = std::max(strlen(names + p.nameIdx), strlen(names + s.nameIdx));
  int w2 = std::max(numOfDigits(sim.start1), numOfDigits(sim.start2));
  int w3 = std::max(numOfDigits(span1), numOfDigits(span2));
  int w4 = std::max(numOfDigits(p.length), numOfDigits(s.length));
  std::cout << "a score=" << log2(sim.probRatio)+shift << " E=" << evalue
	    << " anchor=" << sim.anchor1 << "," << sim.anchor2 << "\n";
  std::cout << "s " << std::left << std::setw(w1) << names + p.nameIdx << " "
	    << std::right << std::setw(w2) << sim.start1 << " "
	    << std::setw(w3) << span1 << " " << '+' << " "
	    << std::setw(w4) << p.length << " ";
  std::cout.write(seq, length);
  std::cout << "\n";
  std::cout << "s " << std::left << std::setw(w1) << names + s.nameIdx << " "
	    << std::right << std::setw(w2) << sim.start2 << " "
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
			 int iBeg, int jBeg, Float half) {
  long rowSize = sequenceLength + 1;
  const char *seq = sequence + jBeg;

  for (int size = 16; ; size *= 2) {
    int iEnd = std::min(iBeg + size, profile.length + 1);
    int jEnd = std::min(jBeg + size, sequenceLength + 1);
    int jLen = jEnd - jBeg;
    std::vector<Float> scratch2(jLen * 2);
    Float *X = scratch2.data();
    Float *Y = X + jLen;
    Float wSum = 0;

    for (int i = iBeg; i < iEnd; ++i) {
      const Float *Wreverse = scratch + (i+1) * rowSize + jBeg + 1;
      Float a = profile.values[i * profile.width + 0];
      Float b = profile.values[i * profile.width + 1];
      Float d = profile.values[i * profile.width + 2];
      Float e = profile.values[i * profile.width + 3];
      const Float *S = profile.values + i * profile.width + 4;

      Float x = (i == iBeg) ? scale : 0;
      Float z = 0;
      for (int j = 0; j < jLen; ++j) {
	Float y = Y[j];
	Float w = x + y + z;
	wSum += w;
	Y[j] = d * w + e * y;
	x = X[j];
	z = a * w + b * z;
	if (i == profile.length || jBeg + j == sequenceLength) continue;
	X[j] = S[seq[j]] * w;
	if (X[j] * Wreverse[j] > half * scale) {
	  addForwardMatch(alignment, i, jBeg + j);
	}
      }

      if (wSum >= half) return;
    }
  }
}

void addReverseAlignment(std::vector<SegmentPair> &alignment,
			 Profile profile, const char *sequence,
			 int sequenceLength, const Float *scratch,
			 int iEnd, int jEnd, Float half) {
  long rowSize = sequenceLength + 1;
  const char *seq = sequence + jEnd;

  for (int size = 16; ; size *= 2) {
    int iBeg = std::max(iEnd - size, -1);
    int jBeg = std::max(jEnd - size, -1);
    int jLen = jEnd - jBeg;
    std::vector<Float> scratch2(jLen * 2);
    Float *X = scratch2.data() + jLen;
    Float *Y = X + jLen;
    Float wSum = 0;

    for (int i = iEnd-1; i >= iBeg; --i) {
      const Float *Xforward = (i >= 0) ? scratch + i * rowSize + jEnd : 0;
      Float a = profile.values[(i+1) * profile.width + 0];
      Float b = profile.values[(i+1) * profile.width + 1];
      Float d = profile.values[(i+1) * profile.width + 2];
      Float e = profile.values[(i+1) * profile.width + 3];
      const Float *S = (i >= 0) ? profile.values + i * profile.width + 4 : 0;

      Float x = (i == iEnd-1) ? scale : 0;
      Float z = 0;
      for (int j = -1; j >= -jLen; --j) {
	Float y = Y[j];
	Float w = x + d * y + a * z;  // this is: W[i+1][jEnd+j+1]
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
  }
}

void findRawSimilarities(std::vector<RawSimilarity> &similarities,
			 Profile profile, const char *sequence,
			 int sequenceLength, Float *scratch,
			 double minProbRatio) {
  long rowSize = sequenceLength + 1;
  // scratch has space for (sequenceLength + 1) * (profile.length + 2) values
  Float *Y = scratch + rowSize * (profile.length + 1);

  // Backward algorithm:

  Float maxBeg = 0;
  int iMaxBeg, jMaxBeg;

  for (int j = 0; j <= sequenceLength; ++j) Y[j] = 0;

  for (int i = profile.length; i >= 0; --i) {
    Float *W = scratch + i * rowSize;
    const Float *Wfrom = W + rowSize;
    Float a = profile.values[i * profile.width + 0];
    Float b = profile.values[i * profile.width + 1];
    Float d = profile.values[i * profile.width + 2];
    Float e = profile.values[i * profile.width + 3];
    const Float *S = profile.values + i * profile.width + 4;

    Float wOld = 0;
    Float z = 0;
    for (int j = sequenceLength; j >= 0; --j) {
      Float x = S[sequence[j]] * wOld;
      Float y = Y[j];
      Float w = x + d * y + a * z + scale;
      if (w > maxBeg) {
	maxBeg = w;
	iMaxBeg = i;
	jMaxBeg = j;
      }
      wOld = Wfrom[j];
      W[j] = w;
      Y[j] = w + e * y;
      z = w + b * z;
    }
  }

  RawSimilarity begAnchored = {maxBeg, iMaxBeg, jMaxBeg};
  addForwardAlignment(begAnchored.alignment, profile, sequence,
		      sequenceLength, scratch, iMaxBeg, jMaxBeg, maxBeg / 2);

  // Forward algorithm:

  Float maxEnd = 0;
  int iMaxEnd, jMaxEnd;

  Float maxMid = 0;
  Float wBegAnchored, wEndAnchored;
  int iMaxMid, jMaxMid;
  std::vector<SegmentPair> alignment;

  for (int j = 0; j <= sequenceLength; ++j) Y[j] = 0;

  for (int i = 0; i <= profile.length; ++i) {
    Float maxMidHere = maxMid;
    Float *W = scratch + i * rowSize;
    const Float *X = i ? W - rowSize : Y;
    Float a = profile.values[i * profile.width + 0];
    Float b = profile.values[i * profile.width + 1];
    Float d = profile.values[i * profile.width + 2];
    Float e = profile.values[i * profile.width + 3];
    const Float *S = profile.values + i * profile.width + 4;

    Float x = 0;
    Float z = 0;
    for (int j = 0; j <= sequenceLength; ++j) {
      Float y = Y[j];
      Float w = x + y + z + scale;
      Float wMid = w * W[j];
      if (w > maxEnd) {
	maxEnd = w;
	iMaxEnd = i;
	jMaxEnd = j;
      }
      if (wMid > maxMidHere) {
	maxMidHere = wMid;
	wBegAnchored = W[j];
	wEndAnchored = w;
	jMaxMid = j;
      }
      if (minProbRatio > 0 && wMid >= minProbRatio) {
	RawSimilarity raw = {wMid / scale, i, j};
	addReverseAlignment(raw.alignment, profile, sequence, sequenceLength,
			    scratch, i, j, w / 2);
	reverse(raw.alignment.begin(), raw.alignment.end());
	addForwardAlignment(raw.alignment, profile, sequence, sequenceLength,
			    scratch, i, j, W[j] / 2);
	similarities.push_back(raw);
      }
      x = X[j];
      W[j] = S[sequence[j]] * w;
      Y[j] = d * w + e * y;
      z = a * w + b * z;
    }

    if (maxMidHere > maxMid) {
      maxMid = maxMidHere;
      iMaxMid = i;
      if (minProbRatio < 0) {
	alignment.clear();
	addForwardAlignment(alignment, profile, sequence, sequenceLength,
			    scratch, iMaxMid, jMaxMid, wBegAnchored / 2);
      }
    }

  }

  if (minProbRatio <= 0) {
    reverse(alignment.begin(), alignment.end());
    addReverseAlignment(alignment, profile, sequence, sequenceLength,
			scratch, iMaxMid, jMaxMid, wEndAnchored / 2);
    reverse(alignment.begin(), alignment.end());
    RawSimilarity midAnchored = {maxMid / scale, iMaxMid, jMaxMid, alignment};

    RawSimilarity endAnchored = {maxEnd, iMaxEnd, jMaxEnd};
    addReverseAlignment(endAnchored.alignment, profile, sequence,
			sequenceLength, scratch, iMaxEnd, jMaxEnd, maxEnd / 2);
    reverse(endAnchored.alignment.begin(), endAnchored.alignment.end());

    similarities.push_back(endAnchored);
    similarities.push_back(begAnchored);
    similarities.push_back(midAnchored);
  }
}

void findSimilarities(std::vector<Similarity> &similarities,
		      Profile profile, const char *sequence,
		      int sequenceLength, Float *scratch,
		      size_t profileNum, size_t strandNum,
		      double minProbRatio) {
  const char *alphabet =
    (profile.width == 9) ? "acgtn" : "ACDEFGHIKLMNPQRSTVWYX";

  std::vector<RawSimilarity> raws;
  findRawSimilarities(raws, profile, sequence, sequenceLength, scratch,
		      minProbRatio);

  for (const auto &r : raws) {
    Similarity s = {r.probRatio, profileNum, strandNum, r.anchor1, r.anchor2,
		    r.anchor1, r.anchor2};
    if (!r.alignment.empty()) {
      s.start1 = r.alignment[0].start1;
      s.start2 = r.alignment[0].start2;
      addAlignedProfile(s.alignedSequences, r.alignment, alphabet, profile);
      addAlignedSequence(s.alignedSequences, r.alignment, alphabet, sequence);
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
		 int numOfSequences, Float *scratch) {
  std::mt19937_64 randGen;
  int alphabetSize = profile.width - 5;
  std::discrete_distribution<> dist(letterFreqs, letterFreqs + alphabetSize);

  std::vector<double> scores(numOfSequences * 3);
  double *endScores = scores.data();
  double *begScores = endScores + numOfSequences;
  double *midScores = begScores + numOfSequences;

  std::cout << "#trial\tend-anchored\t\tstart-anchored\t\tmid-anchored\n"
    "#\tprofPos\tseqPos\tscore\tprofPos\tseqPos\tscore\tprofPos\tseqPos\tscore"
	    << std::endl;

  for (int i = 0; i < numOfSequences; ++i) {
    for (int j = 0; j < sequenceLength; ++j) sequence[j] = dist(randGen);
    for (int j = 0; j < border; ++j) sequence[sequenceLength+j] = sequence[j];
    sequence[sequenceLength + border] = dist(randGen);  // arbitrary letter
    std::vector<RawSimilarity> sims;
    findRawSimilarities(sims, profile, sequence, sequenceLength + border,
			scratch, 0);
    endScores[i] = log(sims[0].probRatio);
    begScores[i] = log(sims[1].probRatio);
    midScores[i] = log(sims[2].probRatio);
    std::cout << (i+1) << "\t"
	      << sims[0].anchor1 << "\t" << sims[0].anchor2 << "\t"
	      << log2(sims[0].probRatio)+shift << "\t"
	      << sims[1].anchor1 << "\t" << sims[1].anchor2 << "\t"
	      << log2(sims[1].probRatio)+shift << "\t"
	      << sims[2].anchor1 << "\t" << sims[2].anchor2 << "\t"
	      << log2(sims[2].probRatio)+shift << std::endl;
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

  std::cout << "#\tend-\tstart-\tmid-anchored\n";

  std::cout << "#lamMM\t" << MMendL << "\t" << MMbegL << "\t" << MMmidL << "\n"

	    << "#kMM\t" << MMendK / pow(s, MMendL) << "\t"
	    << MMbegK / pow(s, MMbegL) << "\t"
	    << MMmidK / pow(s, MMmidL) << "\n"

	    << "#kMM1\t" << MMendKsimple/scale << "\t" << MMbegKsimple/scale
	    << "\t" << MMmidKsimple/scale << "\n";

  std::cout << "#lamML\t" << MLendL << "\t" << MLbegL << "\t" << MLmidL << "\n"

	    << "#kML\t" << MLendK / pow(s, MLendL) << "\t"
	    << MLbegK / pow(s, MLbegL) << "\t"
	    << MLmidK / pow(s, MLmidL) << "\n"

	    << "#kML1\t" << MLendKsimple/scale << "\t" << MLbegKsimple/scale
	    << "\t" << MLmidKsimple/scale << "\n";

  std::cout << "#lamLM\t" << LMendL << "\t" << LMbegL << "\t" << LMmidL << "\n"

	    << "#kLM\t" << LMendK / pow(s, LMendL) << "\t"
	    << LMbegK / pow(s, LMbegL) << "\t"
	    << LMmidK / pow(s, LMmidL) << std::endl;

  Triple h = {MLendKsimple, MLbegKsimple, MLmidKsimple};
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

  // set the final epsilon to the geometric mean of the other epsilons
  end[3] = geometricMean(p.values + p.width + 3, p.length - 1, p.width);

  // set the background letter probabilities proportional to the
  // geometric mean of the foreground letter probabilities
  double sum = 0;
  for (int k = 4; k < p.width - 1; ++k) {
    double m = geometricMean(p.values + k, p.length, p.width);
    if (m <= 0) return 0;
    end[k] = m;
    sum += m;
  }
  for (int k = 4; k < p.width - 1; ++k) end[k] /= sum;

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
    for (int k = 4; k < p.width - 1; ++k) {
      probs[k] = c * probs[k] / end[k];
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
	int k = 5;
	while (iss >> word && strchr(word.c_str(), '.')) {  // xxx "*"?
	  double prob = probFromText(word.c_str());
	  if (prob > 1) return 0;
	  values.push_back(prob);
	  ++k;
	}
	values.push_back(0);
	if (k == 5) return 0;
	if (profile.width > 0 && k != profile.width) return 0;
	profile.width = k;
	profile.length += 1;
	if (profile.length + 1 > INT_MAX / profile.width) return 0;
	state = 2;
      }
    }
  }

  Float *v = &values[0];
  for (size_t i = 0; i < profiles.size(); ++i) {
    profiles[i].values = v;
    if (!finalizeProfile(profiles[i])) return 0;
    v += profiles[i].width * (profiles[i].length + 1);
  }

  return state == 0;
}

void setCharToNumber(char *charToNumber, const char *alphabet) {
  for (int i = 0; alphabet[i]; ++i) {
    int c = alphabet[i];
    charToNumber[toupper(c)] = charToNumber[tolower(c)] = i;
  }
}

int resizeMem(std::vector<Float> &v, int profileLength, int sequenceLength) {
  long rowSize = sequenceLength + 1;
  if (rowSize > LONG_MAX / (profileLength+2)) {
    std::cerr << "too big combination of sequence and profile\n";
    return 0;
  }
  v.resize(rowSize * (profileLength+2));
  return 1;
}

int main(int argc, char* argv[]) {
  int border = 0;
  int strandOpt = 2;
  double evalueOpt = 0;

  const char help[] = "\
usage: seq-position-probs profile.hmm randomTrials randomLength [sequences.fa]\n\
\n\
Get similarity scores between profiles and random sequences, estimate Gumbel K\n\
parameters, and optionally get scores and E-values for real sequences.\n\
\n\
options:\n\
  -h, --help        show this help message and exit\n\
  -b B, --border B  add a border of this length to each random sequence\n\
  -s S, --strand S  strand: 0=reverse, 1=forward, 2=both, ignored for protein\n\
                    (default: 2)\n\
";

  const char sOpts[] = "hb:e:s:";

  static struct option lOpts[] = {
    {"help",   no_argument,       0, 'h'},
    {"border", required_argument, 0, 'b'},
    {"evalue", required_argument, 0, 'e'},
    {"strand", required_argument, 0, 's'},
    {0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      std::cout << help;
      return 0;
    case 'b':
      border = intFromText(optarg);
      if (border < 0) {
	std::cerr << help;
	return 1;
      }
      break;
    case 'e':
      evalueOpt = strtod(optarg, 0);
      if (evalueOpt <= 0) {
	std::cerr << help;
	return 1;
      }
      break;
    case 's':
      strandOpt = intFromText(optarg);
      if (strandOpt < 0 || strandOpt > 2) {
	std::cerr << help;
	return 1;
      }
      break;
    case '?':
      std::cerr << help;
      return 1;
    }
  }

  if (argc - optind < 3 || argc - optind > 4) {
    std::cerr << help;
    return 1;
  }

  int numOfSequences = intFromText(argv[optind + 1]);
  int sequenceLength = intFromText(argv[optind + 2]);

  if (numOfSequences < 1) {
    std::cerr << "bad numOfSequences\n";
    return 1;
  }

  if (sequenceLength < 1 || sequenceLength == INT_MAX) {
    std::cerr << "bad sequenceLength\n";
    return 1;
  }

  if (INT_MAX - border <= sequenceLength) {
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
  charVec.resize(seqIdx + sequenceLength + border + 1);
  std::vector<Float> scratch;
  if (!resizeMem(scratch, maxProfileLength, sequenceLength + border)) return 1;

  std::cout << "# Length of random sequence: " << sequenceLength << "\n";

  double totEndK = 0;
  double totBegK = 0;
  double totMidK = 0;

  for (size_t i = 0; i < numOfProfiles; ++i) {
    Profile p = profiles[i];
    const Float *profileEnd = p.values + p.width * p.length;
    std::cout << "# Profile name: " << &charVec[p.nameIdx] << "\n";
    std::cout << "# Profile length: " << p.length << "\n";
    std::cout << "# Background letter probabilities:";
    for (int j = 4; j < p.width - 1; ++j) std::cout << " " << profileEnd[j];
    std::cout << "\n";
    Triple r = estimateK(p, profileEnd+4, &charVec[seqIdx], sequenceLength,
			 border, numOfSequences, &scratch[0]);
    totEndK += r.endAnchored;
    totBegK += r.begAnchored;
    totMidK += r.midAnchored;
  }

  if (argc - optind < 4 || numOfProfiles < 1) return 0;

  int width = profiles[0].width;
  for (size_t i = 1; i < numOfProfiles; ++i) {
    if (profiles[i].width != width) width = 0;
  }
  char charToNumber[256];
  memset(charToNumber, width - 5, 256);
  if (width == 9) {
    setCharToNumber(charToNumber, "ACGT");
    setCharToNumber(charToNumber, "ACGU");
  } else if (width == 25) {
    setCharToNumber(charToNumber, "ACDEFGHIKLMNPQRSTVWY");
    strandOpt = 1;
  } else {
    std::cerr << "the profiles should be all protein, or all nucleotide\n";
    return 1;
  }

  charVec.resize(seqIdx);
  std::vector<Sequence> sequences;
  std::vector<Similarity> similarities;
  size_t totSequenceLength = 0;

  std::ifstream file;
  std::istream &in = openFile(file, argv[optind + 3]);
  if (!file) return 1;
  Sequence sequence;
  for (size_t i = 0; readSequence(in, sequence, charVec, charToNumber); ++i) {
    if (!resizeMem(scratch, maxProfileLength, sequence.length)) return 1;
    seqIdx = charVec.size() - sequence.length - 1;
    totSequenceLength += sequence.length * (strandOpt / 2 + 1);
    double minProbRatio =
      (evalueOpt > 0) ? totMidK * totSequenceLength / evalueOpt * scale : -1;
    for (int s = 0; s < 2; ++s) {
      if (s != strandOpt) {
	size_t strandNum = i * 2 + s;
	for (size_t j = 0; j < numOfProfiles; ++j) {
	  findSimilarities(similarities, profiles[j], &charVec[seqIdx],
			   sequence.length, &scratch[0], j, strandNum,
			   minProbRatio);
	}
      }
      reverseComplement(&charVec[seqIdx], &charVec[seqIdx] + sequence.length);
    }
    sequences.push_back(sequence);
    charVec.resize(seqIdx);
  }

  std::cout << "# Total sequence length: " << totSequenceLength << "\n";

  double endKMN = totEndK * totSequenceLength;
  double begKMN = totBegK * totSequenceLength;
  double midKMN = totMidK * totSequenceLength;

  std::cout.precision(3);
  for (size_t i = 0; i < similarities.size(); ++i) {
    Profile p = profiles[similarities[i].profileNum];
    Sequence s = sequences[similarities[i].strandNum / 2];
    double kmn = (evalueOpt > 0) ? midKMN :
      (i % 3 == 0) ? endKMN : (i % 3 == 1) ? begKMN : midKMN;
    double evalue = kmn / similarities[i].probRatio;
    if (evalueOpt <= 0 && i % 3 == 0) std::cout << "\n";
    if (evalueOpt > 0 && evalue > evalueOpt) continue;
    printSimilarity(charVec.data(), p, s, similarities[i], evalue);
  }

  return 0;
}
