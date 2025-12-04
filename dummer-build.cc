// Author: Martin C. Frith 2025
//         Patrick Styll   2025
// SPDX-License-Identifier: BSD-3-Clause

#include "dummer-util.hh"

#include "priors/blocks9.hh"
#include "priors/wheeler4.hh"
#include "priors/gap-priors.hh"

#include <algorithm>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <getopt.h>

// Copied from HMMER 3.4:
#define OPT_symfrac 0.50
#define OPT_ere_aa  0.59
#define OPT_ere_nt  0.62
#define OPT_esigma  45.0

// Baum Welch parameters:
#define OPT_bw_maxiter 1000
#define OPT_bw_maxDiff 1e-6

// down-scale probabilities by this amount, to delay overflow:
const double scale1 = 1.0 / (1<<30) / (1<<30) / (1<<30) / (1<<30);
const double scale = scale1 * scale1 * scale1 * scale1;

int verbosity = 0;

struct DirichletMixture {
  const double *params;
  int componentCount;
};

struct MultipleAlignment {
  int sequenceCount;
  int alignmentLength;
  std::vector<unsigned char> alignment;
  std::vector<unsigned char> columnStatus;
  std::string name;
};

void setCharToNumber(unsigned char *charToNumber, const char *alphabet) {
  for (int i = 0; alphabet[i]; ++i) {
    int c = alphabet[i];
    charToNumber[toupper(c)] = charToNumber[tolower(c)] = i;
  }
}

void normalize(double *x, int n) {
  double s = std::accumulate(x, x + n, 0.0);
  assert(s > 0);
  for (int i = 0; i < n; ++i) x[i] /= s;
}

// some MSAs (if --pnone is used) may have zero probabilities
// log(0)/log(NaN) would produce errors
double geometricMean(const double *values, int length, int step) {
  const double minProb = 1e-6;
  double s = 0.0;
  for (int i = 0; i < length; ++i) {
    double v = values[i * step];
    assert(isfinite(v));
    v = std::max(v, minProb);
    s += log(v);
  }
  return exp(s / length);
}

void setBackgroundProbs(double *bgProbs, const double *probs,
			int alphabetSize, int profileLength) {
  for (int i = 0; i < alphabetSize; ++i) {
    bgProbs[i] = geometricMean(probs + i, profileLength, 7 + alphabetSize);
  }
  normalize(bgProbs, alphabetSize);
}

// Probabilities are estimated from counts and a prior probability
// distribution (pseudocounts).

// If we want probabilities fitted to the input multiple alignment, we
// should use option --enone: that makes maxCountSum effectively
// infinite, so it has no effect.

// However, we often want probabilities fitted to more-distantly
// related sequences than those in the input.  To do that, we
// downscale the counts (relative to the prior) when the total count
// is high (> maxCountSum).  This is slightly different from
// (supposedly, better than) HMMER's "entropy weighting".

double getProb(double count1, double count2,
	       double pseudocount1, double pseudocount2, double maxCountSum) {
  double s = count1 + count2;
  // avoid zero division when pnone removes priors
  if (s + (pseudocount1 + pseudocount2) == 0.0) return 0.5;
  double r = (s > maxCountSum) ? maxCountSum / s : 1.0;
  return (count1 * r + pseudocount1) / (s * r + (pseudocount1 + pseudocount2));
}

void applyDirichletMixture(DirichletMixture dmix,
			   int alphabetSize, double maxCountSum,
			   const double *counts, double *probEstimate) {
  int componentSize = 1 + alphabetSize;  // 1 mixture coefficient + alphas
  double countSum = std::accumulate(counts, counts + alphabetSize, 0.0);
  double rescale = (countSum > maxCountSum) ? maxCountSum / countSum : 1.0;
  std::vector<double> logs(dmix.componentCount);
  std::fill_n(probEstimate, alphabetSize, 0.0);

  for (int i = 0; i < dmix.componentCount; ++i) {
    if (dmix.componentCount == 1) break;  // allows for alpha values that = 0
    const double *alphas = dmix.params + i * componentSize + 1;
    double alphaSum = std::accumulate(alphas, alphas + alphabetSize, 0.0);
    logs[i] = lgamma(alphaSum) - lgamma(countSum * rescale + alphaSum + 1);
    for (int j = 0; j < alphabetSize; ++j) {
      logs[i] += lgamma(counts[j] * rescale + alphas[j]) - lgamma(alphas[j]);
    }
  }

  double maxLog = *std::max_element(logs.begin(), logs.end());

  for (int i = 0; i < dmix.componentCount; ++i) {
    double componentProb = dmix.params[i * componentSize];
    const double *alphas = dmix.params + i * componentSize + 1;
    double mul = componentProb * exp(logs[i] - maxLog);
    for (int j = 0; j < alphabetSize; ++j) {
      probEstimate[j] += mul * (counts[j] * rescale + alphas[j]);
    }
  }

  normalize(probEstimate, alphabetSize);
}

void gapCountsToProbs(const GapPriors &gp, double maxCountSum,
		      const double *counts, double *probs) {
  double match  = counts[0]; //gamma
  double notIns = counts[1]; //gamma + delta + tau = etap + oldGamma + epsilonp
  double delBeg = counts[2]; //delta
  double insBeg = counts[3]; //alpha
  double insEnd = insBeg;    //betap
  double insExt = counts[4]; //beta
  double delEnd = counts[5]; //epsilonp
  double delExt = counts[6]; //epsilon

  double gpNotIns = gp.match + gp.delStart;
  assert(notIns > 0);
  double a = getProb(insBeg, notIns, gp.insStart, gpNotIns, maxCountSum);
  double d = getProb(delBeg, match, gp.delStart, gp.match, maxCountSum);
  probs[1] = a;  // insertion start probability
  probs[0] = (1 - a) * (1 - d);
  probs[2] = (1 - a) * d;  // deletion start probability

  double b = getProb(insExt, insEnd, gp.insExtend, gp.insEnd, maxCountSum);
  probs[3] = 1 - b;
  probs[4] = b;  // insertion extend probability

  double e = getProb(delExt, delEnd, gp.delExtend, gp.delEnd, maxCountSum);
  probs[5] = 1 - e;
  probs[6] = e;  // deletion extend probability
}

void countsToProbs(DirichletMixture dmix, const GapPriors &gp,
		   int alphabetSize, double maxCountSum,
		   int profileLength, const double *counts, double *probs) {
  int countsPerPosition = 7 + alphabetSize;

  for (int i = 0; ; ++i) {
    const double *c = counts + i * countsPerPosition;
    /* */ double *p = probs  + i * countsPerPosition;

    gapCountsToProbs(gp, maxCountSum, c, p);

    if (i == profileLength) break;

    applyDirichletMixture(dmix, alphabetSize, maxCountSum, c + 7, p + 7);
  }

  // These gap probabilities are never used - set them similarly to HMMER:

  probs[5] = 1;  // delEnd
  probs[6] = 0;  // delExtend

  double *p = probs + profileLength * countsPerPosition;
  p[0] = 1 - p[1];  // match
  p[2] = 0;         // delStart

  setBackgroundProbs(p + 7, probs + 7, alphabetSize, profileLength);
}

double relativeEntropy(const double *probs,
		       int alphabetSize, int profileLength) {
  int probsPerPosition = 7 + alphabetSize;
  const double *bgProbs = probs + profileLength * probsPerPosition;

  double r = 0;
  for (int i = 0; i < profileLength; ++i) {
    for (int j = 0; j < alphabetSize; ++j) {
      if (probs[j] > 0) r += probs[j] * log2(probs[j] / bgProbs[j]);
    }
    probs += probsPerPosition;
  }
  return r;
}

double entropyWeight(DirichletMixture dmix, const GapPriors &gp,
		     int alphabetSize, double rangeEnd, double targetRelEnt,
		     int profileLength, const double *counts, double *probs) {
  double rangeBeg = 0;
  double rangeLen = rangeEnd;
  double maxCountSum = 1e37;  // start with no limit on total weight

  while (1) {
    countsToProbs(dmix, gp, alphabetSize, maxCountSum, profileLength,
		  counts, probs);
    if (verbosity > 1) {
      const double *bgProbs = probs + profileLength * (7 + alphabetSize) + 7;
      std::cerr << "Background letter probabilities:" << std::setprecision(3);
      for (int i = 0; i < alphabetSize; ++i) std::cerr << " " << bgProbs[i];
      std::cerr << std::setprecision(6) << "\n";
    }
    double r = relativeEntropy(probs + 7, alphabetSize, profileLength);
    if (verbosity) std::cerr << "Max total weight: " << maxCountSum
			     << "  Relative entropy: " << r << "\n";
    if (maxCountSum >= rangeEnd && r <= targetRelEnt) return rangeEnd;
    if (rangeLen < 0.01) return maxCountSum;
    if (r < targetRelEnt) rangeBeg = maxCountSum;
    rangeLen /= 2;
    maxCountSum = rangeBeg + rangeLen;
  }
}

std::istream &readGapPriors(std::istream &in, GapPriors &gp) {
  in >> gp.match >> gp.insStart >> gp.delStart
     >> gp.insEnd >> gp.insExtend >> gp.delEnd >> gp.delExtend;
  if (!in) std::cerr << "can't read the gap pseudocounts\n";
  return in;
}

std::istream &readDirichletMixture(std::istream &in, DirichletMixture &dmix,
				   std::vector<double> &params) {
  int alphabetSize = 0;
  int componentCount = 0;
  double d = 0;
  in >> alphabetSize >> componentCount;
  for (int i = 0; i < componentCount; ++i) {
    in >> d;
    params.push_back(d);
    for (int j = 0; j < alphabetSize; ++j) {
      in >> d;
      params.push_back(d);
    }
  }
  dmix.params = params.data();
  dmix.componentCount = componentCount;
  if (componentCount < 1) in.setstate(std::ios::failbit);
  if (!in) std::cerr << "can't read the Dirichlet mixture file\n";
  return in;
}

std::istream &checkAlignmentBlock(std::istream &in, MultipleAlignment &ma,
				  const std::vector<std::string> &lines,
				  const std::vector<std::string> &refLines,
				  int lineCount) {
  if (lineCount) {
    if (ma.sequenceCount && lineCount != ma.sequenceCount) {
      return fail(in, "alignment blocks should have same number of sequences");
    }
    ma.sequenceCount = lineCount;
  }

  if (!refLines.empty()) {
    size_t blockCount = ma.sequenceCount ? lines.size() / ma.sequenceCount : 0;
    if (refLines.size() != blockCount) fail(in, "bad RF lines");
    if (refLines.back().size() != lines.back().size()) {
      return fail(in, "RF length != aligned sequence length");
    }
  }

  return in;
}

std::istream &readMultipleAlignment(std::istream &in, MultipleAlignment &ma) {
  std::vector<std::string> lines;
  std::vector<std::string> refLines;
  std::string line, word;
  int lineCount = 0;
  ma.sequenceCount = 0;
  ma.alignmentLength = 0;
  ma.name = "Anon";

  while (getline(in, line)) {
    std::istringstream iss(line);
    iss >> word;
    if (!iss) {
      if (!checkAlignmentBlock(in, ma, lines, refLines, lineCount)) return in;
      lineCount = 0;
    } else if (word == "#=GF") {
      iss >> word;
      if (word == "ID") {
	iss >> ma.name;
      }
    } else if (word == "//") {
      break;
    } else if (word == "#=GC") {
      iss >> word;
      if (word != "RF") continue;
      iss >> word;
      if (!iss) return fail(in, "bad RF line");
      if (ma.sequenceCount && refLines.empty()) fail(in, "bad RF lines");
      refLines.push_back(word);
    } else if (word[0] != '#') {
      iss >> word;
      if (!iss) return fail(in, "bad alignment file");
      if (lineCount) {
	if (word.size() != lines.back().size()) {
	  return fail(in, "aligned sequences should have the same length");
	}
      } else {
	ma.alignmentLength += word.size();
      }
      lines.push_back(word);
      ++lineCount;
    }
  }

  if (!checkAlignmentBlock(in, ma, lines, refLines, lineCount)) return in;

  ma.alignment.resize(ma.sequenceCount * ma.alignmentLength);
  unsigned char *a = ma.alignment.data();

  for (int i = 0; i < ma.sequenceCount; ++i) {
    for (size_t j = 0; j < lines.size(); j += ma.sequenceCount) {
      for (size_t k = 0; k < lines[i + j].size(); ++k) {
	unsigned char c = lines[i + j][k];
	*a++ = c;
      }
    }
  }

  ma.columnStatus.resize(refLines.empty() ? 0 : ma.alignmentLength);
  unsigned char *s = ma.columnStatus.data();

  for (size_t j = 0; j < refLines.size(); ++j) {
    for (size_t k = 0; k < refLines[j].size(); ++k) {
      char c = refLines[j][k];
      *s++ = (c != '.' && c != '-' && c != '_' && c != '~');
    }
  }

  return in;
}

bool isProteinAlignment(const MultipleAlignment &ma) {
  int counts[256] = {0};
  for (auto x : ma.alignment) ++counts[x];

  int tot = 0;
  for (int x = 'A'; x <= 'Z'; ++x) {
    int y = tolower(x);
    tot += counts[x] + counts[y];
  }

  int dna = 0;
  const char dnaChars[] = "ACGTNU";
  for (int i = 0; dnaChars[i]; ++i) {
    int x = dnaChars[i];
    int y = tolower(x);
    dna += counts[x] + counts[y];
  }

  if (verbosity) std::cerr << "Nucleotide-like letters: " << dna << "\n"
			   << "Total letters: " << tot << "\n";

  return 1.0 * dna / tot < 0.9;  // xxx ???
}

void markEndGaps(MultipleAlignment &ma, int alphabetSize) {
  const int midGap = alphabetSize + 1;
  const int endGap = alphabetSize + 2;

  for (int i = 0; i < ma.sequenceCount; ++i) {
    unsigned char *seq = ma.alignment.data() + i * ma.alignmentLength;
    for (int j = 0; j < ma.alignmentLength && seq[j] == midGap; ++j) {
      seq[j] = endGap;
    }
    for (int j = ma.alignmentLength; j > 0 && seq[j-1] == midGap; --j) {
      seq[j-1] = endGap;
    }
  }
}

void makeSequenceWeights(const MultipleAlignment &ma, int alphabetSize,
			 bool isHand, double symfrac, double *weights) {
  const int midGap = alphabetSize + 1;
  const int endGap = alphabetSize + 2;
  std::vector<int> positionCounts(ma.sequenceCount);

  for (int i = 0; i < ma.alignmentLength; ++i) {
    int counts[32] = {0};
    const unsigned char *seq = ma.alignment.data() + i;
    for (int j = 0; j < ma.sequenceCount; ++j) {
      ++counts[seq[j * ma.alignmentLength]];
    }
    int totalCount = ma.sequenceCount - counts[endGap];
    int nonGapCount = totalCount - counts[midGap];

    if (isHand ? ma.columnStatus[i] : (nonGapCount >= symfrac * totalCount)) {
      int types = 0;
      for (int k = 0; k < alphabetSize; ++k) {
	types += (counts[k] > 0);
      }
      for (int j = 0; j < ma.sequenceCount; ++j) {
	int x = seq[j * ma.alignmentLength];
	if (x < alphabetSize) {
	  weights[j] += 1.0 * nonGapCount / (types * counts[x]);
	  ++positionCounts[j];
	}
      }
    }
  }

  for (int j = 0; j < ma.sequenceCount; ++j) {
    if (positionCounts[j]) weights[j] /= positionCounts[j];
  }

  double maxWeight = *std::max_element(weights, weights + ma.sequenceCount);
  if (maxWeight) {
    for (int j = 0; j < ma.sequenceCount; ++j) weights[j] /= maxWeight;
  }
}

void countEvents(const MultipleAlignment &ma, int alphabetSize,
		 bool isHand, double symfrac,
		 const double *weights, double weightSum, const GapPriors &gp,
		 std::vector<double> &allCounts, std::vector<int> &columns) {
  const int midGap = alphabetSize + 1;
  const int endGap = alphabetSize + 2;
  std::vector<char> states(ma.sequenceCount);
  double notIns = 0;
  double insBeg = 0;
  double insExt = 0;
  double delEnd = 0;

  for (int i = 0; i < ma.alignmentLength; ++i) {
    double counts[32] = {0};
    const unsigned char *seq = ma.alignment.data() + i;
    for (int j = 0; j < ma.sequenceCount; ++j) {
      counts[seq[j * ma.alignmentLength]] += weights[j];
    }
    double totalCount = weightSum - counts[endGap];
    double nonGapCount = totalCount - counts[midGap];

    if (isHand ? ma.columnStatus[i]
	:        (nonGapCount > 0 && nonGapCount >= symfrac * totalCount)) {
      columns.push_back(i);
      double delBeg = 0;
      double delExt = 0;
      for (int j = 0; j < ma.sequenceCount; ++j) {
	int x = seq[j * ma.alignmentLength];
	if (x <= alphabetSize) {  // an aligned letter in a "match" column
	  if (states[j] <= 'd')	notIns += weights[j];
	  if (states[j] == 'd') delEnd += weights[j];
	  states[j] = 'm';
	} else if (x == midGap) {  // a gap in a "match" column
	  if (states[j] != 'd') delBeg += weights[j];
	  // xxx assumes extension (not restart) of a deletion:
	  if (states[j] == 'd') delExt += weights[j];
	  states[j] = 'd';
	}
      }

      // Count adjacent inserts, or adjacent deletes, as extensions or restarts
      double iFrac = 1;  // fraction of adjacent inserts that are extensions
      double dFrac = 1;  // fraction of adjacent deletes that are extensions
      double iExt, dExt, iBeg, iEnd, dBeg, dEnd, nIns, iFracOld, dFracOld;
      do {
	iExt = insExt * iFrac;
	dExt = delExt * dFrac;
	iBeg = insBeg + (insExt - iExt);
	iEnd = iBeg;  // insert end count = insert begin count
	dBeg = delBeg + (delExt - dExt);
	dEnd = delEnd + (delExt - dExt);
	nIns = notIns + (delExt - dExt);
	// Estimate the insertion and deletion probabilities from the counts:
	double gpNotIns = gp.match + gp.delStart;
	double a = getProb(iBeg, nIns, gp.insStart, gpNotIns, 1e37);
	double d = getProb(dBeg, nonGapCount, gp.delStart, gp.match, 1e37);
	double b = getProb(iExt, iEnd, gp.insExtend, gp.insEnd, 1e37);
	double e = getProb(dExt, dEnd, gp.delExtend, gp.delEnd, 1e37);
	// Estimate the extension fractions from the probabilities:
	iFracOld = iFrac;
	dFracOld = dFrac;
	if (insExt > 0) iFrac = b / (b + a * (1 - b));
	if (delExt > 0) dFrac = e / (e + (1 - a) * d * (1 - e));
      } while (fabs(iFrac - iFracOld) > 1e-6 || fabs(dFrac - dFracOld) > 1e-6);

      allCounts.push_back(nonGapCount);
      allCounts.push_back(nIns);
      allCounts.push_back(dBeg);
      allCounts.push_back(iEnd);  // insert end count = insert begin count
      allCounts.push_back(iExt);
      allCounts.push_back(dEnd);
      allCounts.push_back(dExt);
      allCounts.insert(allCounts.end(), counts, counts + alphabetSize);
      insBeg = insExt = delEnd = 0;
      notIns = nonGapCount;
    } else {  // this position is defined as "insertion"
      for (int j = 0; j < ma.sequenceCount; ++j) {
	int x = seq[j * ma.alignmentLength];
	if (x <= alphabetSize) {  // a non-gap symbol
	  if (states[j] <= 'd') notIns += weights[j];
	  if (states[j] == 'd') delEnd += weights[j];
	  if (states[j] != 'i') insBeg += weights[j];
	  // xxx assumes extension (not restart) of an insertion:
	  if (states[j] == 'i') insExt += weights[j];
	  states[j] = 'i';
	}
      }
    }
  }

  allCounts.push_back(0.0);
  allCounts.push_back(notIns);
  allCounts.insert(allCounts.end(), 5, 0.0);  // final gap counts must be 0
  allCounts.insert(allCounts.end(), alphabetSize, 0.0);

  // We'll estimate probabilities from these counts, but
  // the counts assume the alignment is exactly correct and certain.
  // We will mitigate that by using a Baum-Welch algorithm.
}

void getSequenceWithoutGaps(
     const MultipleAlignment &ma, int alphabetSize,
     std::vector<unsigned char> &seqNoGap, std::vector<int> &seqLengths) {
  for (int i = 0; i < ma.sequenceCount; i++) {
    const unsigned char *seq = ma.alignment.data() + i * ma.alignmentLength;
    for (int j = 0; j < ma.alignmentLength; j++) {
      if (seq[j] <= alphabetSize) {
        seqNoGap.push_back(seq[j]);
        seqLengths[i]++;
      }
    }
  }
  seqNoGap.push_back(0);  // so it's safe to read one letter past the end
}

void forward(
     unsigned char *seq, int seqLength, std::vector<double> &probs,
     int profileLength, int width, double *wSum,
     std::vector<double> &X, std::vector<double> &Y, std::vector<double> &Z) {

  int cols = seqLength + 2;
  double aPrime, bPrime, dPrime, ePrime;  // parameters for BW
  double alphaProb, betaProb, deltaProb, epsilonProb, epsilonProb1;

  std::fill_n(&X[0], seqLength+1, 0.0);
  std::fill_n(&Y[1], seqLength+1, 0.0);

  for (int i = 1; i <= profileLength+1; i++) {

    alphaProb    = probs[(i-1) * width + 1];  // insStart
    betaProb     = probs[(i-1) * width + 4];  // insExt
    deltaProb    = probs[(i-1) * width + 2];  // delStart
    epsilonProb  = probs[(i-1) * width + 6];  // delExt

    epsilonProb1 = (i == profileLength+1) ? 0.0 : probs[i * width + 6];

    aPrime = alphaProb * (1.0 - betaProb);
    bPrime = betaProb;

    // can use arbitrary values for S_m, d_m, e_m
    dPrime = deltaProb   * (1 - epsilonProb1);
    ePrime = epsilonProb * (1 - epsilonProb1) / (1 - epsilonProb);
    if (!isfinite(ePrime)) ePrime = 0.0; // avoid NaN issues

    X[(i-1) * cols] = 0.0;
    double z = Z[i * cols] = 0.0;

    for (int j = 1; j <= seqLength+1; j++) {

      int letter = seq[j-1];

      // letter probs / background probs
      double S = (1 - alphaProb - deltaProb)
               * probs[(i-1) * width + (7 + letter)]
               / probs[profileLength * width + (7 + letter)];
      if (!isfinite(S)) S = 0.0; // if letter has 0 background probability

      double y = Y[(i-1) * cols + j];
      double w = X[(i-1) * cols + (j-1)] + y + z + scale;
      z = aPrime * w + bPrime * z;
      *wSum += w;

      X[i * cols + j] = S * w;
      Y[i * cols + j] = dPrime * w + ePrime * y;
      Z[i * cols + j] = z;
    }
  }
}

void backward(unsigned char *seq, int seqLength,
     std::vector<double> &probs, int profileLength, int width,
     std::vector<double> &Wbar, std::vector<double> &Ybar, std::vector<double> &Zbar) {

  int cols = seqLength + 2;
  double aPrime, bPrime, dPrime, ePrime;  // parameters for BW
  double alphaProb, betaProb, deltaProb, epsilonProb, epsilonProb1;

  std::fill_n(&Wbar[(profileLength+1) * cols + 1], seqLength+1, 0.0);
  std::fill_n(&Ybar[(profileLength+1) * cols + 0], seqLength+1, 0.0);

  for (int i = profileLength; i >= 0; i--) {

    alphaProb    = probs[i * width + 1];  // insStart
    betaProb     = probs[i * width + 4];  // insExt
    deltaProb    = probs[i * width + 2];  // delStart
    epsilonProb  = probs[i * width + 6];  // delExt

    epsilonProb1 = (i == profileLength) ? 0.0 : probs[(i+1) * width + 6];

    aPrime = alphaProb * (1.0 - betaProb);
    bPrime = betaProb;

    // can use arbitrary values for S_m, d_m, e_m
    dPrime = deltaProb   * (1 - epsilonProb1);
    ePrime = epsilonProb * (1 - epsilonProb1) / (1 - epsilonProb);
    if (!isfinite(ePrime)) ePrime = 0.0; // avoid NaN issues

    Wbar[(i+2) * cols - 1] = 0.0;
    double z = Zbar[(i+1) * cols - 1] = 0.0;

    for (int j = seqLength; j >= 0; j--) {

      int letter = seq[j];

      // letter probs / background probs
      double S = (1 - alphaProb - deltaProb)
               * probs[i * width + (7 + letter)]
               / probs[profileLength * width + (7 + letter)];
      if (!isfinite(S)) S = 0.0; // if letter has 0 background probability

      double x = S * Wbar[(i+1) * cols + (j+1)];
      double y = Ybar[(i+1) * cols + j];
      double w = x + dPrime * y + aPrime * z + scale;
      z = w + bPrime * z;

      Wbar[i * cols + j] = w;
      Ybar[i * cols + j] = w + ePrime * y;
      Zbar[i * cols + j] = z;
    }
  }
}

void calculateTransitionCounts(
     std::vector<double> &counts, int profileLength, int seqLength,
     std::vector<double> &Wbar, std::vector<double> &Ybar,
     std::vector<double> &Zbar, std::vector<double> &X,
     std::vector<double> &Y,    std::vector<double> &Z,
     double wt, int width, int alphabetSize, unsigned char *seq) {

  int cols = seqLength + 2;
  double alpha, beta, delta, epsilon;
  double gamma, betap, epsilonp, etap;
  double oldGamma = 0;

  // estimate counts for every position
  for (int i = 1; i <= profileLength+1; i++) {

    double emis[32] = {0};
    betap    = 0.0;  // expected count of (1 - beta)
    beta     = 0.0;
    epsilonp = 0.0;  // expected count of (1 - epsilon)
    epsilon  = 0.0;
    etap     = 0.0;  // expected count of (1 - eta)

    for (int j = 1; j <= seqLength+1; j++) {
      betap    += Z[i * cols + (j-1)] * Wbar[(i-1) * cols + (j-1)];
      beta     += Z[i * cols + (j-1)] * Zbar[(i-1) * cols + (j-1)];
      epsilonp += Y[(i-1) * cols + j] * Wbar[(i-1) * cols + (j-1)];
      epsilon  += Y[(i-1) * cols + j] * Ybar[(i-1) * cols + (j-1)];
      etap     += Wbar[(i-1) * cols + (j-1)];
      int letter = seq[j-1];
      emis[letter] += X[i * cols + j] * Wbar[i * cols + j];
    }

    beta  -= betap;
    assert(beta >= 0);

    epsilon  -= epsilonp;
    assert(epsilon >= 0);

    // expected count of alpha
    alpha = betap;

    gamma = std::accumulate(emis, emis + alphabetSize + 1, 0.0);

    double notIns = etap * scale + oldGamma + epsilonp;
    oldGamma = gamma;

    // update the HMM parameters
    counts[(i-1) * width + 0] += gamma    * wt;
    counts[(i-1) * width + 1] += notIns   * wt;

    counts[(i-1) * width + 3] += alpha    * wt;
    counts[(i-1) * width + 4] += beta     * wt;
    counts[(i-1) * width + 5] += epsilonp * wt;
    counts[(i-1) * width + 6] += epsilon  * wt;

    // emission probabilities
    if (i == profileLength+1) continue; // no emissions from the end state
    for (int letter = 0; letter < alphabetSize; letter++) {
      counts[(i-1) * width + (7 + letter)] += emis[letter] * wt;
    }
  }

  for (int i = 1; i <= profileLength; i++) {
    delta = counts[i * width + 5] + counts[i * width + 6]
      - counts[(i-1) * width + 6];
    delta = std::max(delta, 0.0);
    counts[(i-1) * width + 2] = delta;
  }
}

int determineTerminationCondition(double bwMaxDiff,
       std::vector<double> &probsOld, std::vector<double> &probsNew) {
  double maxDiff = 0.0;
  for (size_t i = 0; i < probsOld.size(); i++) {
    double diff = fabs(probsOld[i] - probsNew[i]);
    if (diff > maxDiff) maxDiff = diff;
  }
  return (maxDiff <= bwMaxDiff) ? 1 : 0; // converged or not
}

void baumWelch(std::vector<double> &counts, const MultipleAlignment &ma,
     int alphabetSize, DirichletMixture dmix, const GapPriors &gp,
     int profileLength, const double *weights,
     double maxIter, double bwMaxDiff) {

  // transitions + alphabet
  int width = 7 + alphabetSize;
  std::vector<double> probsOld(counts.size());
  std::vector<double> probsNew(counts.size());

  int rows = profileLength      + 2;
  int cols = ma.alignmentLength + 2;

  std::vector<double> Wbar(rows * cols, 0.0);
  std::vector<double> Ybar(rows * cols, 0.0);
  std::vector<double> Zbar(rows * cols, 0.0);
  std::vector<double> X(rows * cols, 0.0);
  std::vector<double> Y(rows * cols, 0.0);
  std::vector<double> Z(rows * cols, 0.0);

  /* Get gapless sequences and their lengths. */
  std::vector<unsigned char> seqsNoGap;
  std::vector<int>  seqsLengths(ma.sequenceCount, 0);
  getSequenceWithoutGaps(ma, alphabetSize, seqsNoGap, seqsLengths);

  countsToProbs(dmix, gp, alphabetSize, 1e37, profileLength,
      counts.data(), probsOld.data());

  int termCond;
  int num_iterations = 0;

  do {

    /* We aggregate the expected counts over all sequences. */
    for (size_t i = 0; i < counts.size(); i++) counts[i] = 0.0;

    unsigned char* seqNoGap = seqsNoGap.data();

    for (int idx = 0; idx < ma.sequenceCount; idx++) {
      int seqLength = seqsLengths[idx];

      /* Forward pass, calculate X, Y, Z and
         the aggregated v (i.e. sum of w-values). */
      double v = 0.0; // accumulate the sum of weights
      forward(seqNoGap, seqLength, probsOld,
        profileLength, width, &v, X, Y, Z);

      if (!isfinite(v)) {
	std::cerr
	  << "numbers overflowed to infinity in Baum-Welch: quitting\n";
	exit(1);
      }

      /* Backward pass, calculate Wbar, Ybar, Zbar. */
      backward(seqNoGap, seqLength, probsOld,
        profileLength, width, Wbar, Ybar, Zbar);

      double wt = weights[idx] / (v * scale); // weight for this sequence

      /* Calculate and update parameters in new parameter counts. */
      calculateTransitionCounts(counts, profileLength, seqLength,
          Wbar, Ybar, Zbar, X, Y, Z, wt, width, alphabetSize, seqNoGap);

      seqNoGap += seqLength;
    }

    /* Turn the parameter counts into probabilities via priors. */
    countsToProbs(dmix, gp, alphabetSize, 1e37, profileLength,
        counts.data(), probsNew.data());

    /* Determine termination condition.
       Then overwrite the old probabilities with the new ones. */
    termCond = determineTerminationCondition(bwMaxDiff, probsOld, probsNew);

    std::copy(probsNew.begin(), probsNew.end(), probsOld.begin());

    std::cerr << "\rIteration " << ++num_iterations << " / " << maxIter << ": "
              << (termCond ? "    converged" : "not converged");

  } while (!termCond && num_iterations < maxIter);

  std::cerr << "\nError tolerance: " << bwMaxDiff
            << ", iterations: " << num_iterations << ", "
            << (termCond ? "Baum-Welch converged" : "Baum-Welch not converged")
            << "\n";

}

void printProb(bool isCount, double prob) {
  if (isCount) {
    std::cout << "  " << std::left << std::setw(8) << prob << std::right;
  } else if (prob > 0) {
    std::cout << "  " << abs(log(prob));
  } else {
    std::cout << "        *";
  }
}

void printProfile(const double *probs, const int *columns,
		  const char *alphabet, int profileLength, char *argv[],
		  const MultipleAlignment &ma, double neff, bool isCounts) {
  int alphabetSize = strlen(alphabet);
  int width = alphabetSize + 7;
  const double *bgProbs = probs + profileLength * width + 7;

  std::cout << "HMMER3/f [DUMMER "
#include "version.hh"
    "]\n";
  std::cout << "NAME  " << ma.name << "\n";
  std::cout << "LENG  " << profileLength << "\n";
  std::cout << "ALPH  " << (alphabetSize > 4 ? "amino" : "DNA") << "\n";
  std::cout << "MAP   yes\n";
  std::cout << "NSEQ  " << ma.sequenceCount << "\n";
  std::cout << "EFFN  " << neff << "\n";
  std::cout << "BM   ";
  for (char **p = argv; p[1]; ++p) std::cout << " " << *p;
  std::cout << "\n";

  std::cout << "HMM     ";
  for (int i = 0; alphabet[i]; ++i) {
    std::cout << "     " << alphabet[i] << (isCounts ? "    " : "   ");
  }
  std::cout << "\n";

  std::cout << "           match   " <<
    (isCounts ? " notIns    delStart  insEnd    insExt    delEnd    delExt\n"
     :          "insStart delStart insEnd   insExt   delEnd   delExt\n");

  if (isCounts) std::cout << std::defaultfloat;
  std::cout.precision(isCounts ? 3 : 5);
  for (int i = 0; ; ++i) {
    const double *p = probs + i * width;
    std::cout << "        ";
    for (int j = 0; j < alphabetSize; ++j) printProb(isCounts, bgProbs[j]);
    std::cout << "\n";
    std::cout << "        ";
    for (int j = 0; j < 7; ++j) printProb(isCounts, p[j]);
    std::cout << "\n";
    if (i == profileLength) break;
    std::cout << std::setw(7) << i+1 << " ";
    for (int j = 0; j < alphabetSize; ++j) printProb(isCounts, p[7 + j]);
    std::cout << std::setw(7) << columns[i]+1 << "\n";
  }
  std::cout.precision(6);

  std::cout << "//\n";
}

int main(int argc, char* argv[]) {
  bool isHand = false;
  bool isCounts = false;
  double symfrac = OPT_symfrac;
  double ere     = 0;
  double esigma  = OPT_esigma;

  double bwMaxiter  = OPT_bw_maxiter;
  double bwMaxDiff = OPT_bw_maxDiff;
  bool   countOnly   = false;

  const char *dirichletMixtureFileName = 0;
  const char *gapPriorsFileName = 0;
  bool pnone = false;

  const char help[] = "\
usage: dummer-build alignments.stk\n\
\n\
Read multiple sequence alignments: write position-specific (-log) probabilities\n\
of letters, insertions, and deletions.\n\
\n\
Options:\n\
  -h, --help     show this help message and exit\n\
  -V, --version  show version and exit\n\
  -v, --verbose  show progress messages\n\
  --counts       output weighted counts instead of -log probs (implies --enone)\n\
\n\
Options for defining non-insert positions:\n\
  --hand         use alignment RF column annotation lines (ignore symfrac)\n\
  --symfrac F    minimum (weighted) non-gap fraction to define a non-insert\n\
                     position (default: "
    STR(OPT_symfrac) ")\n\
\n\
Options for effective sequence number:\n\
  --enone        ignore relative entropy: set maximum sequence weight to 1\n\
  --ere E        aim for this relative entropy per position (default:\n\
                     " STR(OPT_ere_aa) " for protein, "
    STR(OPT_ere_nt) " for nucleotide)\n\
  --esigma E     aim for this relative entropy (default: "
    STR(OPT_esigma) ")\n\
\n\
Baum-Welch options:\n\
  --maxiter N    maximum number of Baum-Welch iterations (default: "
    STR(OPT_bw_maxiter) ")\n\
  --maxdiff X    convergence threshold for Baum-Welch (default: "
    STR(OPT_bw_maxDiff) ")\n\
  --countonly    only count events, do not run Baum-Welch\n\
\n\
Prior probability options:\n\
  --dmix         Dirichlet mixture file (esl-mixdchlet format)\n\
  --gapprior     file with 7 gap pseudocounts:\n\
                     match insStart delStart insEnd insExtend delEnd delExtend\n\
  --pnone        don't use any prior probabilities (implies --enone)\n\
";

  const char sOpts[] = "hVv";

  static struct option lOpts[] = {
    {"help",      no_argument,       0, 'h'},
    {"version",   no_argument,       0, 'V'},
    {"verbose",   no_argument,       0, 'v'},
    {"symfrac",   required_argument, 0, 's'},
    {"counts",    no_argument,       0, 'C'},
    {"hand",      no_argument,       0, 'H'},
    {"enone",     no_argument,       0, 'n'},
    {"ere",       required_argument, 0, 'p'},
    {"esigma",    required_argument, 0, 'e'},
    {"maxiter",   required_argument, 0, 'N'},
    {"maxdiff",   required_argument, 0, 'X'},
    {"countonly", no_argument,       0, 'c'},
    {"dmix",      required_argument, 0, 'D'},
    {"gapprior",  required_argument, 0, 'G'},
    {"pnone",     no_argument,       0, 'P'},
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
    case 's':
      symfrac = strtod(optarg, 0);
      if (symfrac < 0 || symfrac > 1) return badOpt();
      break;
    case 'C':
      isCounts = true;
      break;
    case 'H':
      isHand = true;
      break;
    case 'n':
      esigma = 1e37;
      break;
    case 'p':
      ere = strtod(optarg, 0);
      if (ere <= 0) return badOpt();
      break;
    case 'e':
      esigma = strtod(optarg, 0);
      if (esigma <= 0) return badOpt();
      break;
    case 'N':
      bwMaxiter = strtol(optarg, 0, 0);
      if (bwMaxiter <= 0) return badOpt();
      break;
    case 'X':
      bwMaxDiff = strtod(optarg, 0);
      if (bwMaxDiff < 0) return badOpt();
      break;
    case 'c':
      countOnly = true;
      break;
    case 'D':
      dirichletMixtureFileName = optarg;
      break;
    case 'G':
      gapPriorsFileName = optarg;
      break;
    case 'P':
      pnone = true;
      esigma = 1e37;
      break;
    case '?':
      std::cerr << help;
      return 1;
    }
  }

  if (argc - optind != 1) {
    std::cerr << help;
    return 1;
  }

  DirichletMixture dmix;
  std::vector<double> dmixParameters;
  if (dirichletMixtureFileName) {
    std::ifstream file;
    std::istream &in = openFile(file, dirichletMixtureFileName);
    if (!file || !readDirichletMixture(in, dmix, dmixParameters)) return 1;
  }

  GapPriors gapPriors = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  if (gapPriorsFileName) {
    std::ifstream file;
    std::istream &in = openFile(file, gapPriorsFileName);
    if (!file || !readGapPriors(in, gapPriors)) return 1;
  }

  std::ifstream file;
  std::istream &in = openFile(file, argv[optind]);
  if (!file) return 1;

  unsigned char charToNumber[256];
  MultipleAlignment ma;

  unsigned long msa_count = 0;
  while (readMultipleAlignment(in, ma)) {
    std::cout << std::fixed;
    std::cerr << "MSA #" << ++msa_count << ": " << ma.name << std::endl;
    if (isHand && ma.columnStatus.empty()) {
      return err("missing RF annotation line: needed for --hand");
    }
    bool isProtein = isProteinAlignment(ma);
    const char *alphabet = isProtein ? "ACDEFGHIKLMNPQRSTVWY" : "ACGT";
    int alphabetSize = strlen(alphabet);
    memset(charToNumber, alphabetSize, 256);
    setCharToNumber(charToNumber, alphabet);
    setCharToNumber(charToNumber, isProtein ? "AUDEFGHIO" : "ACGU");
    charToNumber['-'] = charToNumber['.']
      = charToNumber['_'] = charToNumber['~'] = alphabetSize + 1;

    double myEre = (ere > 0) ? ere : isProtein ? OPT_ere_aa : OPT_ere_nt;

    if (dirichletMixtureFileName) {
      int dmixAlphabetSize = dmixParameters.size() / dmix.componentCount - 1;
      if (dmixAlphabetSize != alphabetSize) {
	return err("the Dirichlet mixture file has wrong alphabet size");
      }
    } else if (pnone) {
      dmixParameters.resize(1 + alphabetSize);
      dmixParameters[0] = 1;
      dmix.params = dmixParameters.data();
      dmix.componentCount = 1;
    } else {
      dmix.params = isProtein ? blocks9 : wheeler4;
      dmix.componentCount = (isProtein ? sizeof blocks9 : sizeof wheeler4)
	/ sizeof(double) / (1 + alphabetSize);
    }

    const GapPriors &gp = (gapPriorsFileName || pnone) ? gapPriors :
      isProtein ? mitchisonAaGapPriors : wheelerNtGapPriors;

    if (verbosity) std::cerr << "Alignment length: "
			     << ma.alignmentLength << "\n";

    for (auto &i : ma.alignment) i = charToNumber[i];

    markEndGaps(ma, alphabetSize);

    std::vector<double> weights(ma.sequenceCount);
    makeSequenceWeights(ma, alphabetSize, isHand, symfrac, weights.data());
    double weightSum = std::accumulate(weights.begin(), weights.end(), 0.0);

    if (verbosity > 1) {
      std::cerr << "Sequence weights:\n";
      std::cerr.precision(3);
      for (auto i : weights) std::cerr << i << "\n";
      std::cerr.precision(6);
    }

    if (verbosity) std::cerr << "Total weight: " << weightSum << "\n";

    std::vector<double> counts;
    std::vector<int> columns;
    countEvents(ma, alphabetSize, isHand, symfrac, weights.data(), weightSum,
		gp, counts, columns);
    int profileLength = columns.size();

    if (!countOnly) {
      baumWelch(counts, ma, alphabetSize, dmix, gp,
          profileLength, weights.data(), bwMaxiter, bwMaxDiff);
    }

    if (isCounts) {
      printProfile(counts.data(), columns.data(), alphabet, profileLength,
		   argv, ma, weightSum, isCounts);
      continue;
    }

    std::vector<double> probs(counts.size());

    double targetRelEnt = std::max(esigma, myEre * profileLength);
    if (verbosity) std::cerr << "Target relative entropy: "
			     << targetRelEnt << "\n";
    double neff = entropyWeight(dmix, gp, alphabetSize,
				weightSum, targetRelEnt,
				profileLength, counts.data(), probs.data());

    printProfile(probs.data(), columns.data(), alphabet, profileLength,
		 argv, ma, neff, isCounts);
  }
}
