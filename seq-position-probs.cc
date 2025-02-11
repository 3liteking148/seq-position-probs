// Author: Martin C. Frith 2025

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <limits.h>
#include <math.h>
#include <stdlib.h>

typedef float Float;

struct Profile {
  Float *values;
  int width;  // values per position: 4 + alphabetSize
  int length;  // number of positions
  const char *name;
};

bool isDash(const char *text) {
  return text[0] == '-' && text[1] == 0;
}

std::istream &openFile(std::ifstream &file, const char *name) {
  if (isDash(name)) return std::cin;
  file.open(name);
  if (!file) std::cerr << "can't open file: " << name << "\n";
  return file;
}

void maxProbabilityRatios(Profile profile,
			  const unsigned char *sequence, int sequenceLength,
			  Float *scratch,
			  double &maxEndRatio, double &maxBegRatio,
			  double &maxMidRatio) {
  // scratch has space for (sequenceLength + 1) * (profile.length + 2) values
  Float *Y = scratch + (sequenceLength + 1) * (profile.length + 1);

  // Backward algorithm:

  for (int j = 0; j <= sequenceLength; ++j) {
    Y[j] = 0;
  }

  for (int i = profile.length; i >= 0; --i) {
    Float *W = scratch + i * (sequenceLength + 1);
    const Float *Wfrom = W + (sequenceLength + 1);
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
      Float w = x + d * y + a * z + 1;
      wOld = Wfrom[j];
      W[j] = w;
      Y[j] = w + e * y;
      z = w + b * z;
    }
  }

  // Forward algorithm:

  Float maxEnd = 0;
  Float maxBeg = 0;
  Float maxMid = 0;

  for (int j = 0; j <= sequenceLength; ++j) {
    Y[j] = 0;
  }

  for (int i = 0; i <= profile.length; ++i) {
    Float *W = scratch + i * (sequenceLength + 1);
    const Float *X = i ? W - (sequenceLength + 1) : Y;
    Float a = profile.values[i * profile.width + 0];
    Float b = profile.values[i * profile.width + 1];
    Float d = profile.values[i * profile.width + 2];
    Float e = profile.values[i * profile.width + 3];
    const Float *S = profile.values + i * profile.width + 4;

    Float x = 0;
    Float z = 0;
    for (int j = 0; j <= sequenceLength; ++j) {
      Float y = Y[j];
      Float w = x + y + z + 1;

      Float s = W[j];
      Float m = w * b;
      if (w > maxEnd) maxEnd = w;
      if (s > maxBeg) maxBeg = s;
      if (m > maxMid) maxMid = m;

      x = X[j];
      W[j] = S[sequence[j]] * w;
      Y[j] = d * w + e * y;
      z = a * w + b * z;
    }
  }

  maxEndRatio = maxEnd;
  maxBegRatio = maxBeg;
  maxMidRatio = maxMid;
}

void scoresOfRandomSequences(Profile profile, const Float *letterFreqs,
			     unsigned char *sequence, int sequenceLength,
			     int numOfSequences, Float *scratch) {
  std::mt19937_64 randGen;
  int alphabetSize = profile.width - 4;
  std::discrete_distribution<> dist(letterFreqs, letterFreqs + alphabetSize);

  double endSum = 0;
  double begSum = 0;
  double midSum = 0;

  std::cout << "#trial\tend-\tstart-\tmid-anchored" << std::endl;

  for (int i = 0; i < numOfSequences; ++i) {
    for (int j = 0; j <= sequenceLength; ++j) sequence[j] = dist(randGen);
    double maxEndRatio, maxBegRatio, maxMidRatio;
    maxProbabilityRatios(profile, sequence, sequenceLength, scratch,
			 maxEndRatio, maxBegRatio, maxMidRatio);
    endSum += 1 / maxEndRatio;
    begSum += 1 / maxBegRatio;
    midSum += 1 / maxMidRatio;
    std::cout << (i+1) << "\t" << log2(maxEndRatio) << "\t"
	      << log2(maxBegRatio) << "\t" << log2(maxMidRatio) << std::endl;
  }

  double kEnd = numOfSequences / (profile.length * sequenceLength * endSum);
  double kBeg = numOfSequences / (profile.length * sequenceLength * begSum);
  double kMid = numOfSequences / (profile.length * sequenceLength * midSum);
  std::cout.precision(3);
  std::cout << "#K\t" << kEnd << "\t" << kBeg << "\t" << kMid << "\n";
  std::cout.precision(6);
}

int intFromText(const char *text) {
  long x = strtol(text, 0, 0);
  if (x > INT_MAX || x < INT_MIN) return 0;
  return x;
}

double probFromText(const char *s) {
  if (*s == '*') return 0;
  double d = strtod(s, 0);
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
  std::cerr << "final epsilon: " << end[3] << "\n";

  // set the background letter probabilities proportional to the
  // geometric mean of the foreground letter probabilities
  double sum = 0;
  for (int k = 4; k < p.width; ++k) {
    double m = geometricMean(p.values + k, p.length, p.width);
    if (m <= 0) return 0;
    end[k] = m;
    sum += m;
  }
  for (int k = 4; k < p.width; ++k) {
    end[k] /= sum;
  }

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
    probs[2] = delta * (1 - epsilon1);
    probs[3] = epsilon * (1 - epsilon1) / (1 - epsilon);
    for (int k = 4; k < p.width; ++k) {
      probs[k] = c * probs[k] / end[k];
    }
  }

  return 1;
}

int readProfiles(std::istream &in, std::vector<Profile> &profiles,
		 std::vector<Float> &values, std::vector<char> &names) {
  Profile p = {0, 4, 0, 0};
  int state = 0;
  std::string line, word;
  while (getline(in, line)) {
    std::istringstream iss(line);
    iss >> word;
    switch (state) {
    case 0:
      if (word == "NAME") {
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
	values.insert(values.end(), p.width - 4, 0.0);
	profiles.push_back(p);
	p.width = 4;
	p.length = 0;
	state = 0;
      } else {
	int k = 4;
	while (iss >> word && strchr(word.c_str(), '.')) {  // xxx "*"?
	  double prob = probFromText(word.c_str());
	  if (prob > 1) return 0;
	  values.push_back(prob);
	  ++k;
	}
	if (p.width > 4 && k != p.width) return 0;
	if (k == 4) return 0;
	p.width = k;
	p.length += 1;
	if (p.length + 1 > INT_MAX / p.width) return 0;
	state = 2;
      }
    }
  }

  Float *v = &values[0];
  const char *n = &names[0];
  for (size_t i = 0; i < profiles.size(); ++i) {
    profiles[i].values = v;
    profiles[i].name = n;
    if (!finalizeProfile(profiles[i])) return 0;
    v += profiles[i].width * (profiles[i].length + 1);
    n += strlen(n) + 1;
  }

  return state == 0;
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "please give me: a-file.hmm numOfSequences sequenceLength\n";
    return 1;
  }

  int numOfSequences = intFromText(argv[2]);
  int sequenceLength = intFromText(argv[3]);

  if (numOfSequences < 1) {
    std::cerr << "bad numOfSequences\n";
    return 1;
  }

  if (sequenceLength < 1 || sequenceLength == INT_MAX) {
    std::cerr << "bad sequenceLength\n";
    return 1;
  }

  std::vector<char> profileNames;
  std::vector<Float> profileValues;
  std::vector<Profile> profiles;

  {
    std::ifstream file;
    std::istream &in = openFile(file, argv[1]);
    if (!file) return 1;
    if (!readProfiles(in, profiles, profileValues, profileNames)) {
      std::cerr << "can't read the profile data\n";
      return 1;
    }
  }

  size_t numOfProfiles = profiles.size();

  int maxProfileLength = 0;
  for (size_t i = 0; i < numOfProfiles; ++i) {
    maxProfileLength = std::max(maxProfileLength, profiles[i].length);
  }

  if (sequenceLength+1 > INT_MAX / (maxProfileLength+2)) {
    std::cerr << "too big combination of sequence and profile\n";
    return 1;
  }

  std::vector<unsigned char> sequence(sequenceLength+1);
  std::vector<Float> scratch((sequenceLength+1) * (maxProfileLength+2));

  std::cout << "# Sequence length: " << sequenceLength << "\n";

  for (size_t i = 0; i < numOfProfiles; ++i) {
    Profile p = profiles[i];
    const Float *profileEnd = p.values + p.width * p.length;
    std::cout << "# Profile name: " << p.name << "\n";
    std::cout << "# Profile length: " << p.length << "\n";
    std::cout << "# Background letter probabilities:";
    for (int j = 4; j < p.width; ++j) {
      std::cout << " " << profileEnd[j];
    }
    std::cout << "\n";
    scoresOfRandomSequences(p, profileEnd+4, &sequence[0], sequenceLength,
			    numOfSequences, &scratch[0]);
  }
}
