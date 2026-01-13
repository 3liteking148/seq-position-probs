// Author: Martin C. Frith 2026
// SPDX-License-Identifier: BSD-3-Clause

#include "tantan/tantan.cc"

#include <math.h>

// Write the probability that each sequence position is repetitive
// into probabilities.
// The sequence should consist of small numbers.  If isProtein,
// 0,1,2,...,19 mean ACD...Y in that order.  Otherwise, 0,1,2,3 mean 4
// bases (in any order).
// Slightly larger numbers are also tolerated (e.g. for selenocysteine).
void calcTantanProbabilities(const unsigned char *sequence, int length,
			     bool isProtein, float *probabilities) {
  const int blosum62[20][20] = {
    { 4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-2,-1,-1,-1, 1, 0, 0,-3,-2},
    { 0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2},
    {-2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0,-1,-3,-4,-3},
    {-1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0,-1,-2,-3,-2},
    {-2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3},
    { 0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3, 0,-2,-2,-2, 0,-2,-3,-2,-3},
    {-2,-3,-1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1,-2,-3,-2, 2},
    {-1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-1, 3,-3,-1},
    {-1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0,-1,-2,-3,-2},
    {-1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-1, 1,-2,-1},
    {-1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1, 1,-1,-1},
    {-2,-3, 1, 0,-3, 0, 1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2},
    {-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2, 7,-1,-2,-1,-1,-2,-4,-3},
    {-1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0,-1,-2,-2,-1},
    {-1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2},
    { 1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2},
    { 0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 1, 5, 0,-2,-2},
    { 0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2, 0, 4,-3,-1},
    {-3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11, 2},
    {-2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7},
  };  // in alphabetical order of one-letter amino acid codes

  const int blosum62minScore = -4;
  const double blosum62lambda = 0.32403216734804408;

  double probRatioMatrix[22][22];  // >20 for unusual letters: enough for now
  const double *matrixRows[22];

  for (int i = 0; i < 22; ++i) {
    matrixRows[i] = probRatioMatrix[i];

    for (int j = 0; j < 22; ++j) {
      if (isProtein) {
	int score = (i < 20 && j < 20) ? blosum62[i][j] : blosum62minScore;
	probRatioMatrix[i][j] = exp(blosum62lambda * score);
      } else {
	// match score = +1, mismatch score = -1
	// => match ratio = alphabetSize-1, mismatch ratio = 1/(alphabetSize-1)
	probRatioMatrix[i][j] = (i < 4 && i == j) ? 3.0 : 1.0 / 3.0;
      }
    }
  }

  int maxRepeatOffset = isProtein ? 50 : 100;

  tantan::getProbabilities(sequence, sequence + length, maxRepeatOffset,
			   matrixRows, 0.005, 0.05, 0.9, 0, 0, probabilities);
}
