// Author: Martin C. Frith 2025

#include <algorithm>
#include <fstream>
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

struct Profile {  // position-specific (insert, delete, letter) probabilities
  Float *values;  // probabilities or probability ratios
  int width;  // values per position: 4 + alphabetSize + 1
  int length;  // number of positions
  size_t nameIdx;
  double endK;  // Gumbel K for end-anchored sum of alignment probabilities
  double begK;  // Gumbel K for start-anchored sum of alignment probabilities
  double midK;  // Gumbel K for mid-anchored sum of alignment probabilities
};

struct Sequence {
  size_t nameIdx;
  int length;
};

struct Result {
  double endAnchored;
  double begAnchored;
  double midAnchored;
};

void reverseComplement(char *beg, char *end) {
  while (beg < end) {
    char c = *--end;
    *end = 3 - *beg;
    *beg++ = 3 - c;
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

Result maxProbabilityRatios(Profile profile, const char *sequence,
			    int sequenceLength, Float *scratch) {
  long rowSize = sequenceLength + 1;
  // scratch has space for (sequenceLength + 1) * (profile.length + 2) values
  Float *Y = scratch + rowSize * (profile.length + 1);

  // Backward algorithm:

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

  for (int j = 0; j <= sequenceLength; ++j) Y[j] = 0;

  for (int i = 0; i <= profile.length; ++i) {
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
      Float w = x + y + z + 1;

      Float wBeg = W[j];
      Float wMid = w * wBeg;
      if (w > maxEnd) maxEnd = w;
      if (wBeg > maxBeg) maxBeg = wBeg;
      if (wMid > maxMid) maxMid = wMid;

      x = X[j];
      W[j] = S[sequence[j]] * w;
      Y[j] = d * w + e * y;
      z = a * w + b * z;
    }
  }

  Result result = {maxEnd, maxBeg, maxMid};
  return result;
}

void estimateK(Profile &profile, const Float *letterFreqs,
	       char *sequence, int sequenceLength, int numOfSequences,
	       Float *scratch) {
  std::mt19937_64 randGen;
  int alphabetSize = profile.width - 5;
  std::discrete_distribution<> dist(letterFreqs, letterFreqs + alphabetSize);

  double harmEnd = 0;
  double harmBeg = 0;
  double harmMid = 0;

  std::cout << "#trial\tend-\tstart-\tmid-anchored score" << std::endl;

  for (int i = 0; i < numOfSequences; ++i) {
    for (int j = 0; j <= sequenceLength; ++j) sequence[j] = dist(randGen);
    Result r = maxProbabilityRatios(profile, sequence, sequenceLength,
				    scratch);
    harmEnd += 1 / r.endAnchored;
    harmBeg += 1 / r.begAnchored;
    harmMid += 1 / r.midAnchored;
    std::cout << (i+1) << "\t" << log2(r.endAnchored) << "\t"
	      << log2(r.begAnchored) << "\t"
	      << log2(r.midAnchored) << std::endl;
  }

  profile.endK = numOfSequences / (sequenceLength * harmEnd);
  profile.begK = numOfSequences / (sequenceLength * harmBeg);
  profile.midK = numOfSequences / (sequenceLength * harmMid);
  std::cout.precision(3);
  std::cout << "#K\t" << profile.endK << "\t" << profile.begK << "\t"
	    << profile.midK << "\n";
  std::cout.precision(6);
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
  std::cerr << "final epsilon: " << end[3] << "\n";

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
  int strandOpt = 2;

  const char help[] = "\
usage: seq-position-probs profile.hmm randomTrials randomLength [sequences.fa]\n\
\n\
Get similarity scores between profiles and random sequences, estimate Gumbel K\n\
parameters, and optionally get scores and E-values for real sequences.\n\
\n\
options:\n\
  -h, --help        show this help message and exit\n\
  -s S, --strand S  strand: 0=reverse, 1=forward, 2=both, ignored for protein\n\
                    (default: 2)\n\
";

  const char sOpts[] = "hs:";

  static struct option lOpts[] = {
    {"help",   no_argument,       0, 'h'},
    {"strand", required_argument, 0, 's'},
    {0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      std::cout << help;
      return 0;
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
  charVec.resize(seqIdx + sequenceLength + 1);
  std::vector<Float> scratch;
  if (!resizeMem(scratch, maxProfileLength, sequenceLength)) return 1;

  std::cout << "# Length of random sequence: " << sequenceLength << "\n";

  double totEndK = 0;
  double totBegK = 0;
  double totMidK = 0;

  for (size_t i = 0; i < numOfProfiles; ++i) {
    Profile &p = profiles[i];
    const Float *profileEnd = p.values + p.width * p.length;
    std::cout << "# Profile name: " << &charVec[p.nameIdx] << "\n";
    std::cout << "# Profile length: " << p.length << "\n";
    std::cout << "# Background letter probabilities:";
    for (int j = 4; j < p.width - 1; ++j) std::cout << " " << profileEnd[j];
    std::cout << "\n";
    estimateK(p, profileEnd+4, &charVec[seqIdx], sequenceLength,
	      numOfSequences, &scratch[0]);
    totEndK += p.endK;
    totBegK += p.begK;
    totMidK += p.midK;
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
  std::vector<Result> results;
  size_t totSequenceLength = 0;

  std::ifstream file;
  std::istream &in = openFile(file, argv[optind + 3]);
  if (!file) return 1;
  Sequence sequence;
  while (readSequence(in, sequence, charVec, charToNumber)) {
    if (!resizeMem(scratch, maxProfileLength, sequence.length)) return 1;
    seqIdx = charVec.size() - sequence.length - 1;
    for (int s = 0; s < 2; ++s) {
      if (s != strandOpt) {
	totSequenceLength += sequence.length;
	for (size_t j = 0; j < numOfProfiles; ++j) {
	  results.push_back(maxProbabilityRatios(profiles[j], &charVec[seqIdx],
						 sequence.length,
						 &scratch[0]));
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

  std::cout << "#seq\tseqLen\tprofile\tproLen\tstrand\t"
    "EAscore\tE-value\tSAscore\tE-value\tMAscore\tE-value\n";
  std::cout.precision(3);
  const Result *r = results.empty() ? 0 : &results[0];
  for (size_t i = 0; i < sequences.size(); ++i) {
    Sequence s = sequences[i];
    for (int k = 0; k < 2; ++k) {
      if (k == strandOpt) continue;
      for (size_t j = 0; j < numOfProfiles; ++j) {
	Profile p = profiles[j];
	std::cout << &charVec[s.nameIdx] << "\t" << s.length << "\t"
		  << &charVec[p.nameIdx] << "\t" << p.length << "\t"
		  << "+-"[k] << "\t"
		  << log2(r->endAnchored) << "\t"
		  << endKMN / r->endAnchored << "\t"
		  << log2(r->begAnchored) << "\t"
		  << begKMN / r->begAnchored << "\t"
		  << log2(r->midAnchored) << "\t"
		  << midKMN / r->midAnchored << "\n";
	++r;
      }
    }
  }

  return 0;
}
