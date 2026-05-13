// Author: Martin C. Frith 2025
// SPDX-License-Identifier: BSD-3-Clause

// See [Frith2025]: "Simple and thorough detection of related
// sequences with position-varying probabilities of substitutions,
// insertions, and deletions", MC Frith 2025

#include "dummer-util.hh"
#include "tantan-wrapper.hh"
// clang-format off
#include "can_i_haz_simd.hh"
// clang-format on

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
#include <map>
#include <memory>
#include <optional>
#include <queue>

#include <Kokkos_SIMD.hpp>

#define OPT_e 10
#define OPT_s 2
#define OPT_m 3
#define OPT_t 30
#define OPT_l 1000
#define OPT_b 100
#define OPT_x 100 // 0 to enable greedy mode

#define EVALUE
#define ALIGN

// uncomment to enable I_1, I_2 edge to I_0
//#define ENABLE_FS_INSERT_EXTENSION

// uncomment to enable D_1, D_2 states
// D1, D2 cannot transition to another delete state for now (only to junction state)
#define ENABLE_FS_DELETE_STATES

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
using simd_t = Kokkos::Experimental::simd<Float>;
constexpr auto simdWidth = simd_t::size();

// for now they have to be the same?
const Float STOP_CODON_PROB = 0.001;
const Float BG_STOP_CODON_PROB = 0.001;

// reverse engineered from transmark
const Float FRAMESHIFT1_MULTIPLIER = 0.005; // 0.005 each for delete and insert
const Float FRAMESHIFT2_MULTIPLIER = (0.005 / 2);

#define BACKGROUND_FRAMESHIFT_RATE (0.01)

int simdRoundUp(int x) { // lowest multiple of simdLen that is >= x
    return x - 1 - (x - 1) % simdLen + simdLen;
}

Float simdHorizontalMax(SimdFloat x) { // assuming it doesn't need to be fast
    Float y[simdLen];
    simdStore(y, x);
    return *std::max_element(y, y + simdLen);
}

SimdFloat simdPowersFwd(Float x) {
    Float a[simdLen];
    a[0] = x;
    for (int i = 1; i < simdLen; ++i)
        a[i] = a[i - 1] * x;
    return simdLoad(a);
}

SimdFloat simdPowersRev(Float x) {
    Float a[simdLen];
    a[simdLen - 1] = x;
    for (int i = simdLen - 1; i > 0; --i)
        a[i - 1] = a[i] * x;
    return simdLoad(a);
}

// Only consider similarities that are local maxima.  If 2
// similarities have identical 1st anchor coordinates, and their 2nd
// anchor coordinates are closer than this, omit the lower-scoring one.
const int minSeparation = 32; // xxx ???

// down-scale probabilities by this amount, to delay overflow:
const Float scale = 1.0 / (1 << 30) / (1 << 30) / (1 << 3); // sqrt[min normal float]
const Float invScale = 1.0 / scale;
const int shift = 63; // add this to scores, to undo the scaling

int verbosity = 0;

const int nonLetterWidth = 9; // number of non-letter values per position

struct Params { // TODO: maybe SIMD order
    Float alpha_prime[3];
    Float beta_prime[3];
#ifdef ENABLE_FS_DELETE_STATES
    Float delta_prime[3];
#else
    Float delta_prime;
#endif
    Float epsilon_prime;
    Float enter_match_probability;
};

struct Profile {   // position-specific (insert, delete, letter) probabilities
    Float *values; // probabilities or probability ratios
    std::vector<Params> values_v2;
    std::vector<Float> bg_probs, log2_bg_probs;
    int width;  // number of values per position
    int length; // number of positions
    size_t nameIdx;
    size_t consensusSequenceIdx;
    double gumbelKendAnchored, gumbelKbegAnchored, gumbelKmidAnchored, lambda;
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
    int anchor2; // 2nd anchor coordinate (don't need to store the 1st one)
};

struct SequenceData {
    std::vector<uint8_t> decoded;
    std::string sequence;
    std::string maskedSequence;
    Contig contig;
    size_t strandNum;
};

struct SequenceRequest {
    std::shared_ptr<SequenceData> seqData;
    Float minProbRatio;

    bool operator<(const SequenceRequest &other) const {
        return seqData->decoded.size() < other.seqData->decoded.size();
    }

    bool operator>(const SequenceRequest &other) const {
        return seqData->decoded.size() > other.seqData->decoded.size();
    }
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
    return x.alignment.empty() ? x.anchor2 : x.alignment.back().start2 + x.alignment.back().length;
}

double mean(const double *x, int n) {
    double s = 0;
    for (int i = 0; i < n; ++i)
        s += x[i];
    return s / n;
}

int numOfDigits(int x) {
    int n = 0;
    do
        ++n;
    while (x /= 10);
    return n;
}

const char *getAlphabet(int alphabetSize) {
    assert(alphabetSize == 20 || alphabetSize == 4);
    return alphabetSize == 20  ? "ACDEFGHIKLMNPQRSTVWYUO?*" // 20 + 2 amino acids
           : alphabetSize == 4 ? "ACGT"
                               : 0;
}

char complement(char c) {
    // Map DNA bases correctly using the protein alphabet indices
    // "ACDEFGHIKLMNPQRSTVWYUO?*" -> A=0, C=1, G=5, T=16
    switch (c) {
    case 0:
        return 16; // A -> T
    case 16:
        return 0; // T -> A
    case 1:
        return 5; // C -> G
    case 5:
        return 1; // G -> C
    default:
        return c; // Fallback for masked '?' or unrecognized chars
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
        if (!(in >> x))
            return in;
        if (x != '>')
            return fail(in, "bad sequence data: no '>'");
        std::string line, word;
        getline(in, line);
        std::istringstream iss(line);
        if (!(iss >> word))
            return fail(in, "bad sequence data: no name");
#ifdef PIPELINE_MODE
        size_t slash = word.find('/');
        std::string chr = word.substr(0, slash);

        std::string range = word.substr(slash + 1);
        size_t dash = range.find('-');

        sequence.w_start = std::stoi(range.substr(0, dash));
        sequence.w_end = std::stoi(range.substr(dash + 1));

        std::string length, profile, strand;
        if (!(iss >> length >> profile >> strand))
            return fail(in, "bad filtered sequence data: no true length, profile, or strand");

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
        if (verbosity > 0)
            std::cerr << "Sequence: " << name << "\n";
    }

    size_t seqIdx = vec.size();
    std::streambuf *buf = in.rdbuf();
    int c = buf->sgetc();

    while (c != std::streambuf::traits_type::eof() && c != '>') {
        if (charToNumber[c] < 125)
            break; // found a contig symbol
        if (c > ' ')
            ++sequence.length; // skip over non-contig symbols
        c = buf->snextc();
    }

    while (c != std::streambuf::traits_type::eof() && charToNumber[c] < 126) {
        if (c > ' ') {
            vec.push_back(charToNumber[c]);
        }
        c = buf->snextc();
    }

    size_t seqLen = vec.size() - seqIdx;
    if (seqLen > INT_MAX - 2 * simdLen)
        return fail(in, "sequence is too long!");
    contig.start = sequence.length;
    contig.length = seqLen;
    sequence.length += seqLen;
    return in;
}

char profileLetter(const char *alphabet, char letterCode) {
    return alphabet[letterCode & 31] + (letterCode & 32); // upper/lowercase
}

void addAlignedProfile(std::vector<char> &gappedSeq, const std::vector<SegmentPair> &alignment,
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

char seqLetter(const char *alphabet, const char *sequence, const char *maskedSequence,
               int position) {
    char c = sequence[position];
    return alphabet[c] + (maskedSequence[position] > c) * 32; // upper/lowercase
}

void addAlignedSequence(std::vector<char> &gappedSeq, const std::vector<SegmentPair> &alignment,
                        const char *alphabet, const char *sequence, const char *maskedSequence) {
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


void printSimilarity(const char *names, Profile &p, Sequence s, const FinalSimilarity &sim,
                     double evalue) {
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
    std::cout << "a score=" << (log2(sim.probRatio) + shift) << " E=" << evalue
              << " anchor=" << sim.anchor1 << "," << anchor2 << "\n";
    std::cout << "s " << std::left << std::setw(w1) << names + p.nameIdx << " " << std::right
              << std::setw(w2) << sim.start1 << " " << std::setw(w3) << span1 << " " << '+' << " "
              << std::setw(w4) << p.length << " ";
    std::cout.write(seq, length);
    std::cout << "\n";
    std::cout << "s " << std::left << std::setw(w1) << names + s.nameIdx << " " << std::right
              << std::setw(w2) << start2 << " " << std::setw(w3) << span2 << " " << strand << " "
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
        if (x.start1 + x.length > pos1 || x.start2 + x.length > pos2)
            return;
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
        if (x.start1 <= pos1 || x.start2 <= pos2)
            return;
    }
    SegmentPair sp = {pos1, pos2, 1};
    alignment.push_back(sp);
}

template <typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U> &a, const std::pair<T, U> &b) {
    return {a.first + b.first, b.second};
}

using DP_Cell = Float;

template <typename T, bool Rolling = false> class FlatMatrix {
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
        } else {
            data.resize(r * c);
        }
    }

    void assign(size_t r, size_t c, T init = T()) {
        logical_rows = r;
        cols = c;
        if constexpr (Rolling) {
            data.assign(2 * c, init);
        } else {
            data.assign(r * c, init);
        }
    }

    inline T &operator()(size_t i, size_t j) {
        // assert(0 <= i && i < logical_rows);
        // assert(0 <= j && j < cols);
        if constexpr (Rolling) {
            return data[(i & 1) * cols + j];
        } else {
            return data[i * cols + j];
        }
    }

    inline const T &operator()(size_t i, size_t j) const {
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

    // Direct row pointer for hot loops — avoids repeated i*cols multiply
    inline T *row_ptr(size_t i) {
        if constexpr (Rolling) {
            return data.data() + (i & 1) * cols;
        } else {
            return data.data() + i * cols;
        }
    }
    inline const T *row_ptr(size_t i) const {
        if constexpr (Rolling) {
            return data.data() + (i & 1) * cols;
        } else {
            return data.data() + i * cols;
        }
    }
};

struct DPScratch {
    FlatMatrix<simd_t> W0, W1, X;
    FlatMatrix<simd_t> X_pfx, X_sfx;

    // Reusable temporary buffers for findSimilarities
    std::vector<simd_t> dp, dp_r;
    std::vector<simd_t> Y0_next, Y0_curr;
    std::vector<simd_t> one, one_sfx;
    std::vector<simd_t> left_side, right_side;
    std::array<std::vector<AlignedSimilarity>, simdWidth> opt_profile_position;
    std::array<std::vector<bool>, simdWidth> aligned;
    std::vector<uint8_t> transposed_decoded;
    std::vector<simd_t> bg_codon_probs;

    // SIMD anchor tracking — avoids scalar per-lane extraction in forward DP
    std::vector<simd_t> best_wMid;  // best wMidAnchored per j (SIMD)
    std::vector<simd_t> best_wEnd;  // corresponding wEndAnchored per j
    std::vector<simd_t> best_i;     // best profile position i per j (SIMD)
};


struct DP_Cell_v2 {
    Float metric;
    int i, j;
    bool emit = false;

    constexpr bool operator< (const DP_Cell_v2 &other) const {
        return metric < other.metric;
    }
};

void addForwardAlignment(int idx, size_t profileLength, size_t sequenceLength, std::vector<SegmentPair> &alignment, int iBeg, int jBeg,
                         DPScratch &scratch) {

    int i = iBeg, j = jBeg;
    while (i <= profileLength && j < sequenceLength) {
        auto choice = std::max({
            DP_Cell_v2{.metric=scratch.X(i, j)[idx] + scratch.X_sfx(i + 1, j + 3)[idx], .i=i + 1, .j=j + 3, .emit=true},
            DP_Cell_v2{.metric=scratch.X_sfx(i + 1, j)[idx], .i=i + 1, .j=j, .emit=false},
            DP_Cell_v2{.metric=scratch.X_sfx(i, j + 1)[idx], .i=i, .j=j + 1, .emit=false},
            DP_Cell_v2{.metric=scratch.left_side[j][idx], .i=INT_MAX, .j=INT_MAX, .emit=false},
        });

        if (choice.emit && j >= 2) {
            addForwardMatch(alignment, i, j - 2);
        }
        i = choice.i, j = choice.j;
    }
}

void addReverseAlignment(int idx, std::vector<SegmentPair> &alignment, int iEnd, int jEnd,
                         DPScratch &scratch) {
    int i = iEnd, j = jEnd;
    while (i >= 0 && j >= 0) {
        DP_Cell opt_succ = 0;
        if (i - 1 >= 0 && j - 3 >= 0) {
            opt_succ = scratch.X_pfx(i - 1, j - 3)[idx];
        }
        auto choice = std::max({
            DP_Cell_v2{.metric=scratch.X(i, j)[idx] + opt_succ, .i=i - 1, .j=j - 3, .emit=true},
            DP_Cell_v2{.metric=(i ? scratch.X_pfx(i - 1, j)[idx] : 0), .i=i - 1, .j=j, .emit=false},
            DP_Cell_v2{.metric=(j ? scratch.X_pfx(i, j - 1)[idx] : 0), .i=i, .j=j - 1, .emit=false},
            DP_Cell_v2{.metric=scratch.right_side[j][idx], .i=INT_MIN, .j=INT_MIN, .emit=false},
        });

        if (choice.emit && j >= 2) {
            addReverseMatch(alignment, i, j - 2);
        }
        i = choice.i, j = choice.j;
    }
}



void addMidAnchored(int idx, size_t profileLength, size_t sequenceLength, std::vector<AlignedSimilarity> &similarities, int anchor1, int anchor2,
                    Float wBegAnchored, Float wEndAnchored, DPScratch &scratch) {
    Float wMidAnchored = wEndAnchored * wBegAnchored;
    AlignedSimilarity s = {wMidAnchored / scale, anchor1, anchor2, wEndAnchored};
#ifdef ALIGN
    addForwardAlignment(idx, profileLength, sequenceLength, s.alignment, anchor1, anchor2, scratch);
#endif
    similarities.push_back(s);
}

void finishMidAnchored(int idx, AlignedSimilarity &s, DPScratch &scratch) {
    reverse(s.alignment.begin(), s.alignment.end());
#ifdef ALIGN
    addReverseAlignment(idx, s.alignment, s.anchor1, s.anchor2, scratch);
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
            if (i.start1 - i.start2 == j.start1 - j.start2 && i.start1 + i.length > j.start1 &&
                i.start1 < j.start1 + j.length)
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
            if (simBeg2(y) >= end)
                break;
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

int updateInitialSimilarities(InitialSimilarity *sims, int count, int anchor2, Float probRatio) {
    int i = 0;
    int j = 0;
    while (i < count && sims[i].anchor2 <= anchor2 - minSeparation)
        ++i;
    while (i < count && sims[i].probRatio > probRatio)
        sims[j++] = sims[i++];
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
    if (aa2codons.size() > 0) {
        return aa2codons;
    }

    aa2codons['A'] = {"GCT", "GCC", "GCA", "GCG"};               // Ala
    aa2codons['R'] = {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}; // Arg
    aa2codons['N'] = {"AAT", "AAC"};                             // Asn
    aa2codons['D'] = {"GAT", "GAC"};                             // Asp
    aa2codons['C'] = aa2codons['U'] = {"TGT", "TGC"};            // Cys
    aa2codons['Q'] = {"CAA", "CAG"};                             // Gln
    aa2codons['E'] = {"GAA", "GAG"};                             // Glu
    aa2codons['G'] = {"GGT", "GGC", "GGA", "GGG"};               // Gly
    aa2codons['H'] = {"CAT", "CAC"};                             // His
    aa2codons['I'] = {"ATT", "ATC", "ATA"};                      // Ile
    aa2codons['L'] = {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}; // Leu
    aa2codons['K'] = aa2codons['O'] = {"AAA", "AAG"};            // Lys
    aa2codons['M'] = {"ATG"};                                    // Met
    aa2codons['F'] = {"TTT", "TTC"};                             // Phe
    aa2codons['P'] = {"CCT", "CCC", "CCA", "CCG"};               // Pro
    aa2codons['S'] = {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}; // Ser
    aa2codons['T'] = {"ACT", "ACC", "ACA", "ACG"};               // Thr
    aa2codons['W'] = {"TGG"};                                    // Trp
    aa2codons['Y'] = {"TAT", "TAC"};                             // Tyr
    aa2codons['V'] = {"GTT", "GTC", "GTA", "GTG"};               // Val
    aa2codons['*'] = {"TAA", "TAG", "TGA"};
    aa2codons['?'] = {"???"}; // masked

    return aa2codons;
}

static void normalize(std::unordered_map<char, double> &m) {
    double s = 0.0;
    for (auto &kv : m)
        s += kv.second;
    if (s <= 0)
        return;
    for (auto &kv : m)
        kv.second /= s;

    // hardcode probabilities
    m['A'] = 0.25;
    m['T'] = 0.25;
    m['G'] = 0.25;
    m['C'] = 0.25;
}
struct NucDist {
    std::unordered_map<char, double> overall{{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
};

NucDist
infer_nucleotide_distribution_equal_synonyms(const std::unordered_map<char, double> &aaFreq) {
    auto aa2codons = build_standard_genetic_code();
    NucDist out;

    // Distribute each amino acid's probability equally across its codons.
    double sm = 0;
    for (auto &kv : aaFreq) {
        char aa = (char)toupper((unsigned char)kv.first);
        auto it = aa2codons.find(aa);
        sm += kv.second;
        // std::cout << kv.first << " probs " << kv.second << std::endl;

        const std::vector<std::string> &codons = it->second;
        double perCodon = kv.second / (double)codons.size();

        for (const std::string &codon : codons) {
            char b1 = codon[0], b2 = codon[1], b3 = codon[2];

            out.overall[b1] += perCodon;
            out.overall[b2] += perCodon;
            out.overall[b3] += perCodon;
        }
    }

    // std::cout << "assert " << sm << " == 1" << std::endl;
    //  At this point:
    //  - overall sums to 3 (because each codon contributes 3 bases) after AA normalization,
    //  - pos1/pos2/pos3 each sum to 1.
    normalize(out.overall);

    return out;
}

// Fast codon translation via flat lookup table (no heap alloc, no hashing)
// Build a flat 32768-entry codon table indexed by packed 5-bit-per-base key
static char codonTableFlat[32768];
static bool codonTableBuilt = false;

static void buildCodonTable() {
    if (codonTableBuilt) return;
    memset(codonTableFlat, '?', sizeof(codonTableFlat));
    static const char *codons[] = {
        "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
        "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
        "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
        "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
        "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
        "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
        "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
        "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"
    };

    static const char aaMap[] = {
        'F','F','L','L','L','L','L','L','I','I','I','M','V','V','V','V',
        'S','S','S','S','P','P','P','P','T','T','T','T','A','A','A','A',
        'Y','Y','*','*','H','H','Q','Q','N','N','K','K','D','D','E','E',
        'C','C','*','W','R','R','R','R','S','S','R','R','G','G','G','G'
    };
    for (int i = 0; i < 64; i++) {
        unsigned key = ((unsigned)(unsigned char)codons[i][0] & 0x1f)
                     | (((unsigned)(unsigned char)codons[i][1] & 0x1f) << 5)
                     | (((unsigned)(unsigned char)codons[i][2] & 0x1f) << 10);
        codonTableFlat[key] = aaMap[i];
    }
    // Also handle lowercase
    for (int i = 0; i < 64; i++) {
        char lc[3] = { (char)(codons[i][0] | 0x20), (char)(codons[i][1] | 0x20), (char)(codons[i][2] | 0x20) };
        unsigned key = ((unsigned)(unsigned char)lc[0] & 0x1f)
                     | (((unsigned)(unsigned char)lc[1] & 0x1f) << 5)
                     | (((unsigned)(unsigned char)lc[2] & 0x1f) << 10);
        codonTableFlat[key] = aaMap[i];
    }
    codonTableBuilt = true;
}

inline char translateFast(const char *dna, int i) {
    unsigned key = ((unsigned)(unsigned char)dna[i] & 0x1f)
                 | (((unsigned)(unsigned char)dna[i+1] & 0x1f) << 5)
                 | (((unsigned)(unsigned char)dna[i+2] & 0x1f) << 10);
    return codonTableFlat[key];
}

// Keep original for compatibility but mark as legacy
char translate(const char *dna, int i) {
    return translateFast(dna, i);
}

Float log2_sum_exp(Float a, Float b) {
    if (a == -INFINITY)
        return b;
    if (b == -INFINITY)
        return a;
    Float m = std::max(a, b);
    return m + log2(exp2(a - m) + exp2(b - m));
}

inline simd_t log2_sum_exp(simd_t a, simd_t b) {
    simd_t m = Kokkos::max(a, b);
    simd_t x = Kokkos::abs(a - b);

    return m + Kokkos::log2(simd_t(1.0) + Kokkos::exp2(-x));
}


// vibe-coded section end

// Single-pass decode: decompress + translate + charToNumber in one loop
std::vector<uint8_t> decodeSequence(const char *sequence, int sequenceLength, const char *alphabet,
                                    const char *charToNumber) {
    buildCodonTable();
    int n = sequenceLength;

    // Decompress sequence in-place to a stack buffer for the codon window
    // We only need a sliding window of 3 decompressed chars at a time
    std::vector<uint8_t> decoded(n, (uint8_t)INT_MIN);

    // Pre-decompress into a flat buffer (avoids per-char string append)
    // Use a local buffer instead of std::string for cache efficiency
    std::vector<char> decompressed(n);
    for (int i = 0; i < n; i++) {
        assert(sequence[i] <= 22);
        decompressed[i] = alphabet[(unsigned char)sequence[i]];
    }

    // Single pass: translate codons and map to numbers
    const char *dec = decompressed.data();
    for (int i = 0; i < n - 2; i++) {
        decoded[i] = charToNumber[(unsigned char)translateFast(dec, i)];
    }
    return decoded;
}

void findSimilarities(std::array<std::vector<AlignedSimilarity>, simdWidth> &similarities, const Profile &profile,
                      const std::array<std::vector<uint8_t>*, simdWidth> &decoded, std::array<Float, simdWidth> minProbRatio,
                      DPScratch &scratch, int activeCount) {
    assert(0 < activeCount && activeCount <= simdWidth);

    int maxSequenceLength = 0;
    for (int idx = 0; idx < activeCount; idx++) {
        maxSequenceLength = std::max(maxSequenceLength, (int)decoded[idx]->size());
    }

    int alphabetSize = profile.width - nonLetterWidth;
    int zero_idx = alphabetSize + 4; // new padded index that maps to 0.0

    scratch.transposed_decoded.assign((maxSequenceLength + 8) * simdWidth, zero_idx);
    uint8_t* transposed_base = scratch.transposed_decoded.data() + 4 * simdWidth;
    for (int idx = 0; idx < activeCount; idx++) {
        for (size_t j = 0; j < decoded[idx]->size(); j++) {
            transposed_base[j * simdWidth + idx] = (*decoded[idx])[j];
        }
    }

    alignas(64) Float rsl_tmp[simdWidth] = {};
    for (int idx = 0; idx < activeCount; idx++) rsl_tmp[idx] = (Float)decoded[idx]->size();
    simd_t realSeqLen = Kokkos::Experimental::simd_unchecked_load<simd_t>(rsl_tmp);

    auto &dp = scratch.dp;
    auto &dp_r = scratch.dp_r;
    dp.resize(maxSequenceLength + 4); dp_r.resize(maxSequenceLength + 4);

    dp_r[maxSequenceLength] = 0;
    for (int i = maxSequenceLength - 1; i >= 0; i--) {
        const Float *bg_probs_ptr = profile.log2_bg_probs.data() + 4;
        const char* indices = (const char*)&transposed_base[i * simdWidth];
        SimdFloat bg_raw = simdLookup(bg_probs_ptr, indices);
        simd_t bg_codon_emit_probs(bg_raw);

        Kokkos::Experimental::simd_mask<Float> msk = i + 2 < realSeqLen;
        Kokkos::Experimental::simd_mask<Float> msk2 = i < realSeqLen;
        simd_t full_codon = (Float)log2(1 - BACKGROUND_FRAMESHIFT_RATE) + bg_codon_emit_probs + dp_r[i + 3];
        simd_t partial_codon = (Float)log2(1 - BACKGROUND_FRAMESHIFT_RATE) + (Float)log2(0.25) * (realSeqLen - (Float)i);
        simd_t t1 = Kokkos::Experimental::condition(msk, full_codon, partial_codon);

        simd_t fs = (Float)log2(BACKGROUND_FRAMESHIFT_RATE * 0.25) + dp_r[i + 1];
        simd_t neg_inf_vec((Float)-INFINITY);
        simd_t t2 = Kokkos::Experimental::condition(msk2, fs, neg_inf_vec);

        dp_r[i] = log2_sum_exp(t1, t2);
    }


    dp[maxSequenceLength] = 0;
    const Float *log2_bg_probs_ptr = profile.log2_bg_probs.data() + 4;
    const Float log2_1_bg_fs = log2(1 - BACKGROUND_FRAMESHIFT_RATE);
    const Float log2_bg_fs_025 = log2(BACKGROUND_FRAMESHIFT_RATE * 0.25);
    const Float log2_025 = log2(0.25);
    const simd_t simd_log2_1_bg_fs(log2_1_bg_fs);
    const simd_t simd_log2_bg_fs_025(log2_bg_fs_025);
    const simd_t simd_log2_025(log2_025);
    const simd_t simd_neg_inf(-INFINITY);

    for (int i = 0; i < maxSequenceLength; i++) {
        // Vectorized bg emission lookup via transposed_base
        const char* indices = (const char*)&transposed_base[(i - 2) * simdWidth];
        SimdFloat bg_raw = simdLookup(log2_bg_probs_ptr, indices);
        simd_t bg_codon_emit_probs(bg_raw);

        // Mask for variable-length sequences (per-lane)
        Kokkos::Experimental::simd_mask<Float> msk_valid = (Float)i < realSeqLen;

        // t1: codon branch — i is scalar, so use plain if/else
        simd_t t1;
        if (i >= 3) {
            t1 = simd_log2_1_bg_fs + bg_codon_emit_probs + dp[i - 3];
        } else if (i == 2) {
            t1 = simd_log2_1_bg_fs + bg_codon_emit_probs;
        } else {
            t1 = simd_log2_1_bg_fs + simd_log2_025 * (Float)(i + 1);
        }

        // t2: frameshift branch
        simd_t t2 = simd_log2_bg_fs_025 + (i > 0 ? dp[i - 1] : simd_t(0));

        // log2_sum_exp and mask out-of-bounds lanes
        simd_t result = log2_sum_exp(t1, t2);
        dp[i] = Kokkos::Experimental::condition(msk_valid, result, simd_neg_inf);
    }

    alignas(64) Float dist1_tmp[simdWidth] = {0}, not_align_probs[simdWidth] = {0};

    for (int idx = 0; idx < activeCount; idx++) {
        if (rsl_tmp[idx] >= 3) {
            not_align_probs[idx] =
            log2_sum_exp(log2_sum_exp(dp[rsl_tmp[idx] - 1][idx],dp[rsl_tmp[idx] - 2][idx]),
                dp[rsl_tmp[idx] - 3][idx]
            );
        } else {
            // For very short sequences, treat as non‑alignable (e.g., -INF or 0)
            not_align_probs[idx] = -INFINITY;
        }

        Float d1 = exp2(-(not_align_probs[idx] / rsl_tmp[idx]));
        dist1_tmp[idx] = d1;
    }

    auto distribute1 = Kokkos::Experimental::simd_unchecked_load<simd_t>(dist1_tmp);
    auto distribute2 = distribute1 * distribute1;
    auto distribute3 = distribute2 * distribute1;

    distribute2 *= (Float)(0.25 * 0.25);
    distribute1 *= (Float)0.25;

    auto not_align_probs_simd = Kokkos::Experimental::simd_unchecked_load<simd_t>(not_align_probs);
    simd_t invRealSeqLen = (Float)1.0 / realSeqLen;
    simd_t nap_div_rsl = not_align_probs_simd * invRealSeqLen;

    scratch.W0.assign(profile.length + 1, maxSequenceLength + 4);
    scratch.W1.assign(profile.length + 2, maxSequenceLength + 4);

    const size_t bufSize = maxSequenceLength + 4;
    auto &Y0_next = scratch.Y0_next; Y0_next.assign(bufSize, 0.0);
    auto &Y0_curr = scratch.Y0_curr; Y0_curr.assign(bufSize, 0.0);

    auto &one = scratch.one; one.assign(bufSize, 0.0);

    for (int j = 0; j < maxSequenceLength; j++) {
        simd_t exponent = -nap_div_rsl * (realSeqLen - (Float)1.0 - (Float)j) + dp_r[j + 1];
        // Prevent the exponent from going into the subnormal range
        //exponent = Kokkos::max(exponent, simd_t(-125.0f));
        one[j] = Kokkos::exp2(exponent);
    }
    auto &one_sfx = scratch.one_sfx; one_sfx = one;
    auto &left_side = scratch.left_side; left_side.assign(one.size(), 0.0);
    auto &right_side = scratch.right_side; right_side.assign(one.size(), 0.0);
#ifdef ALIGN
    scratch.X.resize(profile.length + 2, maxSequenceLength);
    scratch.X_pfx.resize(profile.length + 2, maxSequenceLength + 4);
    scratch.X_sfx.resize(profile.length + 2, maxSequenceLength + 4);
#endif
        const Float *bg_probs_ptr = profile.bg_probs.data() + 4;
        scratch.bg_codon_probs.assign(maxSequenceLength + 8, simd_t(0.0));
        simd_t* bg_codon_probs_base = scratch.bg_codon_probs.data() + 4;
        for (int j = -4; j < maxSequenceLength + 4; j++) {
            const char* indices = (const char*)&transposed_base[j * simdWidth];
            bg_codon_probs_base[j] = simd_t(simdLookup(bg_probs_ptr, indices));
        }

        for (int i = profile.length; i >= 0; i--) {
            const Params &params_cur = profile.values_v2[i];
            const Float *params_emission_probabilities = profile.values + (i)*profile.width + 4;

            // Pre-calculate constants for this i
            const simd_t C_enter = params_cur.enter_match_probability * distribute3;
#ifdef ENABLE_FS_DELETE_STATES
            const simd_t C_delta0 = params_cur.delta_prime[0];
            const simd_t C_delta1 = params_cur.delta_prime[1] * distribute2;
            const simd_t C_delta2 = params_cur.delta_prime[2] * distribute1;
#else
            const simd_t C_delta0 = params_cur.delta_prime;
#endif
            const simd_t C_alpha0 = params_cur.alpha_prime[0] * distribute3;
            const simd_t C_alpha1 = params_cur.alpha_prime[1] * distribute1;
            const simd_t C_alpha2 = params_cur.alpha_prime[2] * distribute2;
            const simd_t C_beta0 = params_cur.beta_prime[0] * distribute3;
#ifdef ENABLE_FS_INSERT_EXTENSION
            const simd_t C_beta1 = params_cur.beta_prime[1] * distribute3;
            const simd_t C_beta2 = params_cur.beta_prime[2] * distribute3;
#endif
            simd_t Z0_ring[4] = {0, 0, 0, 0};
#ifdef ENABLE_FS_INSERT_EXTENSION
            simd_t Z1_ring[4] = {0, 0, 0, 0};
            simd_t Z2_ring[4] = {0, 0, 0, 0};
#endif

            // Raw row pointers — avoid repeated i*cols in inner loop
            simd_t *__restrict__ w1_row_i = scratch.W1.row_ptr(i);
            const simd_t *__restrict__ w1_row_ip1 = scratch.W1.row_ptr(i + 1);
            const simd_t C_eps0 = params_cur.epsilon_prime;
            const simd_t C_scale = scale;

            for (int j = maxSequenceLength - 1; j >= 0; j--) {
                int r_0 = j & 3;
                int r_1 = (j + 1) & 3;
                int r_2 = (j + 2) & 3;
                int r_3 = (j + 3) & 3;

                const char* indices = (const char*)&transposed_base[(j + 1) * simdWidth];
                SimdFloat codon_raw = simdLookup(params_emission_probabilities, indices);
                simd_t codon_emit_probs(codon_raw);
                simd_t bg_codon_emit_probs = bg_codon_probs_base[j + 1];

                simd_t w_val =
                    w1_row_ip1[j + 3] * codon_emit_probs * C_enter +
                    Y0_next[j + 0] * C_delta0 +
#ifdef ENABLE_FS_DELETE_STATES
                    w1_row_ip1[j + 2] * C_delta1 +
                    w1_row_ip1[j + 1] * C_delta2 +
#endif
                    Z0_ring[r_3] * bg_codon_emit_probs * C_alpha0
#ifdef ENABLE_FS_INSERT_EXTENSION
                    +
                    Z1_ring[r_1] * C_alpha1 +
                    Z2_ring[r_2] * C_alpha2
#endif
                 + one[j] * C_scale;
                w1_row_i[j] = w_val;
#ifdef ALIGN
                right_side[j] += w_val;
#endif

                Y0_curr[j] = Kokkos::fma(C_eps0, Y0_next[j], w_val);
                simd_t z0_future = Z0_ring[r_3] * bg_codon_emit_probs;
                Z0_ring[r_0] = Kokkos::fma(C_beta0, z0_future, w_val);
#ifdef ENABLE_FS_INSERT_EXTENSION
                Z1_ring[r_0] = Kokkos::fma(C_beta1, z0_future, w_val);
                Z2_ring[r_0] = Kokkos::fma(C_beta2, z0_future, w_val);
#endif
            }

            std::swap(Y0_curr, Y0_next);
        }

    std::fill(Y0_next.begin(), Y0_next.begin() + bufSize, simd_t(0.0));

    scratch.best_wMid.assign(maxSequenceLength + 4, simd_t(-INFINITY));
    scratch.best_wEnd.assign(maxSequenceLength + 4, simd_t(0.0));
    scratch.best_i.assign(maxSequenceLength + 4, simd_t(-1.0));
    for (int j = 0; j < maxSequenceLength; j++) {
        simd_t exponent = -nap_div_rsl * (Float)(j + 1) + dp[j];
        // Prevent the exponent from going into the subnormal range
        //exponent = Kokkos::max(exponent, simd_t(-125.0f));

        one[j] = Kokkos::exp2(exponent);
    }

#ifdef ALIGN
    // seems to be a clean way of getting expected value of null-sided junctions
    // TODO: verify logic
    for (int j = maxSequenceLength - 1; j >= 0; j--) {
        const char* indices = (const char*)&transposed_base[(j - 2) * simdWidth];
        SimdFloat bg_raw = simdLookup(bg_probs_ptr, indices);
        simd_t bg_codon_emit_probs(bg_raw);

        if (j - 3 >= 0) {
            right_side[j - 3] += (Float)(1 - BACKGROUND_FRAMESHIFT_RATE) * bg_codon_emit_probs *
                                 distribute3 * right_side[j];
        }

        if (j - 1 >= 0) {
            right_side[j - 1] += (Float)(BACKGROUND_FRAMESHIFT_RATE * 0.25) * distribute1 * right_side[j];
        }

        right_side[j] *= one[j]; // TODO: this one specifically (might be off by 1 idk)
        right_side[j] = Kokkos::max(right_side[j], Float(0.0));
    }
    for (int j = 1; j < maxSequenceLength; j++) {
        right_side[j] += right_side[j - 1];
    }
#endif

    for (int i = 0; i <= profile.length; i++) {
        const Params &params_cur = profile.values_v2[i];
        const Float *params_emission_probabilities = profile.values + (i)*profile.width + 4;

        // Pre-calculate constants
        const simd_t C_enter = params_cur.enter_match_probability * distribute3;
        const simd_t C_alpha0 = params_cur.alpha_prime[0];
        const simd_t C_alpha1 = params_cur.alpha_prime[1] * distribute1;
        const simd_t C_alpha2 = params_cur.alpha_prime[2] * distribute2;
        const simd_t C_beta0 = params_cur.beta_prime[0];
#ifdef ENABLE_FS_INSERT_EXTENSION
        const simd_t C_beta1 = params_cur.beta_prime[1];
        const simd_t C_beta2 = params_cur.beta_prime[2];
#endif
#ifdef ENABLE_FS_DELETE_STATES
        const simd_t C_delta0 = params_cur.delta_prime[0];
        const simd_t C_delta1 = params_cur.delta_prime[1] * distribute2;
        const simd_t C_delta2 = params_cur.delta_prime[2] * distribute1;
#else
        const simd_t C_delta0 = params_cur.delta_prime;
#endif
        const simd_t C_eps0 = params_cur.epsilon_prime;
        const simd_t C_scale = scale;

        simd_t Z0_ring[4] = {0, 0, 0, 0};
        simd_t Z1_ring[4] = {0, 0, 0, 0};
        simd_t Z2_ring[4] = {0, 0, 0, 0};

        simd_t one_val_scaled = (Float)(scale);

        // Raw row pointers — avoid repeated i*cols in inner loop
        simd_t *__restrict__ w0_row_i = scratch.W0.row_ptr(i);
        simd_t *__restrict__ w0_row_ip1 = (i + 1 <= profile.length) ? scratch.W0.row_ptr(i + 1) : nullptr;
        const simd_t *__restrict__ w1_row_ip1 = (i + 1 <= profile.length) ? scratch.W1.row_ptr(i + 1) : nullptr;
        const simd_t *__restrict__ w1_row_i = scratch.W1.row_ptr(i);
#ifdef ALIGN
        simd_t *__restrict__ x_row_i = scratch.X.row_ptr(i);
        simd_t *__restrict__ xpfx_row_i = scratch.X_pfx.row_ptr(i);
        const simd_t *__restrict__ xpfx_row_im1 = (i - 1 >= 0) ? scratch.X_pfx.row_ptr(i - 1) : nullptr;
#endif

        // Shift register for w[1..3] — avoids 3 matrix reads per iteration
        simd_t w_shift[3] = {0, 0, 0}; // w_shift[0]=w0(i,j-1), [1]=w0(i,j-2), [2]=w0(i,j-3)
        simd_t pfx_prev = simd_t(0.0); // X_pfx(i, j-1) rolling value

        const simd_t simd_invScale(invScale);

        for (int j = 0; j < maxSequenceLength; j++) {
            // w[1] = W0(i, j-1), w[2] = W0(i, j-2), w[3] = W0(i, j-3)
            simd_t w1, w2, w3;
            if (j == 0) {
                w1 = one_val_scaled; w2 = simd_t(0); w3 = simd_t(0);
            } else if (j == 1) {
                w1 = w_shift[0]; w2 = one_val_scaled; w3 = simd_t(0);
            } else if (j == 2) {
                w1 = w_shift[0]; w2 = w_shift[1]; w3 = one_val_scaled;
            } else {
                w1 = w_shift[0]; w2 = w_shift[1]; w3 = w_shift[2];
            }

            int r_0 = j & 3;
            int r_3 = (j - 3) & 3;

            const char* indices = (const char*)&transposed_base[(j - 2) * simdWidth];
            SimdFloat codon_raw = simdLookup(params_emission_probabilities, indices);

            simd_t codon_emit_probs(codon_raw);
            simd_t bg_codon_emit_probs = bg_codon_probs_base[j - 2];

            simd_t X_ij = C_enter * codon_emit_probs * w3;
            simd_t X_ij_EV = 0;
            if (w1_row_ip1)
                X_ij_EV = X_ij * w1_row_ip1[j] * simd_invScale;

#ifdef ALIGN
            x_row_i[j] = X_ij_EV;

            //
            simd_t opt_succ = 0;
            if (xpfx_row_im1 && j - 3 >= 0) {
                opt_succ = xpfx_row_im1[j - 3];
            }

            simd_t pfx_mx = 0;
            if (xpfx_row_im1)
                pfx_mx = xpfx_row_im1[j];
            if (j - 1 >= 0)
                pfx_mx = Kokkos::max(pfx_mx, pfx_prev);

            pfx_mx = Kokkos::max(pfx_mx, X_ij_EV + opt_succ);
            pfx_mx = Kokkos::max(pfx_mx, right_side[j]);
            xpfx_row_i[j] = pfx_mx;
            pfx_prev = pfx_mx;
            //
#endif

#ifdef ENABLE_FS_INSERT_EXTENSION
            Z0_ring[r_0] =
                bg_codon_emit_probs * distribute3 *
                (C_alpha0 * w3 + C_beta0 * Z0_ring[r_3] +
                    C_beta1 * Z1_ring[r_3] + C_beta2 * Z2_ring[r_3]);
#else
            Z0_ring[r_0] =
                bg_codon_emit_probs * distribute3 *
                (C_alpha0 * w3 + C_beta0 * Z0_ring[r_3]);
#endif
            Z1_ring[r_0] = C_alpha1 * w1;
            Z2_ring[r_0] = C_alpha2 * w2;

            simd_t w0 = w0_row_i[j];
            w0 += Z0_ring[r_0] + Z1_ring[r_0] + Z2_ring[r_0] + one[j] * C_scale;
            w0_row_i[j] = w0;
#ifdef ALIGN
            left_side[j] += w0;
#endif

            // Update shift register
            w_shift[2] = w_shift[1];
            w_shift[1] = w_shift[0];
            w_shift[0] = w0;


#ifdef ENABLE_FS_DELETE_STATES
            Y0_curr[j] = C_delta0 * w0 + C_eps0 * Y0_next[j];
            if (w0_row_ip1)
                w0_row_ip1[j] += X_ij + Y0_curr[j] + C_delta1 * w2 + C_delta2 * w1;
#else
            auto Y0_curr_j = C_delta0 * w0 + C_eps0 * Y0_next[j];
            if (w0_row_ip1)
                w0_row_ip1[j] += X_ij + Y0_curr_j;
#endif

            // SIMD anchor tracking — compute wMid in SIMD, update best values in parallel
            simd_t wBegAnchored = w1_row_i[j];
            simd_t wMidAnchored = w0 * wBegAnchored * simd_invScale;

            Kokkos::Experimental::simd_mask<Float> mask = (wMidAnchored > scratch.best_wMid[j]) && ((Float)j < realSeqLen);
            scratch.best_wMid[j] = Kokkos::Experimental::condition(mask, wMidAnchored, scratch.best_wMid[j]);
            scratch.best_wEnd[j] = Kokkos::Experimental::condition(mask, w0, scratch.best_wEnd[j]);
            scratch.best_i[j] = Kokkos::Experimental::condition(mask, simd_t((Float)i), scratch.best_i[j]);
        }

        std::swap(Y0_curr, Y0_next);
    }

#ifdef ALIGN
    for (int j = 0; j < maxSequenceLength; j++) {
        const char* indices = (const char*)&transposed_base[(j + 1) * simdWidth];
        SimdFloat bg_raw = simdLookup(bg_probs_ptr, indices);

        simd_t bg_codon_emit_probs(bg_raw);

        if (j + 3 < maxSequenceLength) {
            left_side[j + 3] += (Float)(1 - BACKGROUND_FRAMESHIFT_RATE) * bg_codon_emit_probs * distribute3 * left_side[j];
        }

        if (j + 1 < maxSequenceLength) {
            left_side[j + 1] += (Float)(BACKGROUND_FRAMESHIFT_RATE * 0.25) * distribute1 * left_side[j];
        }

        left_side[j] *= one_sfx[j];
        left_side[j] = Kokkos::max(left_side[j], Float(0.0));
    }
    for (int j = maxSequenceLength - 2; j >= 0; j--) {
        left_side[j] += left_side[j + 1];
    }

    for (int i = profile.length; i >= 0; i--) {
        simd_t *__restrict__ xsfx_row_i = scratch.X_sfx.row_ptr(i);
        const simd_t *__restrict__ xsfx_row_ip1 = (i + 1 <= profile.length) ? scratch.X_sfx.row_ptr(i + 1) : nullptr;
        const simd_t *__restrict__ x_row_i = scratch.X.row_ptr(i);

        simd_t opt_right_rolling = simd_t(0.0); // X_sfx(i, j+1) from previous iteration

        for (int j = maxSequenceLength - 1; j >= 0; j--) {
            // Guarded reads for X_sfx
            simd_t opt_succ = (xsfx_row_ip1 && j + 3 < maxSequenceLength)
                              ? xsfx_row_ip1[j + 3] : simd_t(0.0);

            simd_t opt_down = xsfx_row_ip1 ? xsfx_row_ip1[j] : simd_t(0.0);

            auto opt = Kokkos::max(opt_down, opt_right_rolling);
            opt = Kokkos::max(opt, left_side[j]);
            opt = Kokkos::max(opt, x_row_i[j] + opt_succ);
            xsfx_row_i[j] = opt;
            opt_right_rolling = opt;
        }
    }
#endif


    for (int idx = 0; idx < activeCount; idx++) {
        int realSequenceLength = rsl_tmp[idx];
        scratch.opt_profile_position[idx].assign(realSequenceLength, AlignedSimilarity(-INFINITY));
        for (int j = 0; j < realSequenceLength; j++) {
            Float best_prob = scratch.best_wMid[j][idx];
            if (best_prob > -INFINITY) {
                scratch.opt_profile_position[idx][j] = {
                    best_prob,
                    (int)scratch.best_i[j][idx],
                    j,
                    (Float)scratch.best_wEnd[j][idx]
                };
            }
        }

        if (minProbRatio[idx] >= 0) {
            std::ranges::sort(scratch.opt_profile_position[idx], std::greater<>());
            auto &aligned = scratch.aligned;
            aligned[idx].assign(realSequenceLength, false);
            for (auto &aligned_similarity : scratch.opt_profile_position[idx]) {
                if (aligned_similarity.probRatio >= minProbRatio[idx] &&
                    !aligned[idx][aligned_similarity.anchor2]) {
                    addMidAnchored(idx, profile.length, realSequenceLength, similarities[idx], aligned_similarity.anchor1,
                                   aligned_similarity.anchor2,
                                   aligned_similarity.probRatio * scale /
                                       aligned_similarity.wEndAnchored,
                                   aligned_similarity.wEndAnchored, scratch);
                    auto &x = similarities[idx].back();
                    finishMidAnchored(idx, x, scratch);
                    // dumb heuristic (4x length accounting for FS)
                    // todo: silence nuclear fallout
                    int startIdx = std::max(aligned_similarity.anchor2 - 12 * profile.length, 0);
                    int endIdx =
                        std::min(aligned_similarity.anchor2 + 12 * profile.length, realSequenceLength);
                    // std::cout << startIdx << " " << endIdx << std::endl;
                    std::fill(aligned[idx].begin() + startIdx, aligned[idx].begin() + endIdx, true);
                    }
            }
        } else {
            auto sel = *std::max_element(scratch.opt_profile_position[idx].begin(), scratch.opt_profile_position[idx].end());
            AlignedSimilarity b = sel;
            // std::cout << log(sel.probRatio) << std::endl;
            b.probRatio = 0;
            similarities[idx].push_back(b);
            b.probRatio = 0;
            similarities[idx].push_back(b);
            similarities[idx].push_back(sel);
        }
    }

}

int contigToSequencePos(Contig contig, size_t strandNum, int posInContig) {
    return contig.start + strandPosition(strandNum, contig.length, posInContig);
}

void findFinalSimilarities(std::vector<FinalSimilarity> &similarities, std::array<SequenceRequest, simdWidth> &req,
                           const Profile &profile, size_t profileNum, const char *charVec,
                           DPScratch &scratch, int activeCount) {
    const char *alphabet = getAlphabet(profile.width - nonLetterWidth);
    const char *profileSeq = charVec + profile.consensusSequenceIdx;

    std::array<std::vector<AlignedSimilarity>, simdWidth> sims;
    std::array<std::vector<uint8_t>*, simdWidth> decoded;
    std::array<Float, simdWidth> minProbRatio;
    for (int i = 0; i < activeCount; i++) {
        decoded[i] = &req[i].seqData->decoded;
        minProbRatio[i] = req[i].minProbRatio;
    }
    findSimilarities(sims, profile, decoded, minProbRatio, scratch, activeCount);

    for (int idx = 0; idx < activeCount; idx++) {
        const char *sequence = req[idx].seqData->sequence.c_str();
        const char *maskedSequence = req[idx].seqData->maskedSequence.c_str();
        for (const auto &x : sims[idx]) {
            int anchor2 = contigToSequencePos(req[idx].seqData->contig, req[idx].seqData->strandNum, x.anchor2);
            FinalSimilarity s = {x.probRatio, profileNum, req[idx].seqData->strandNum, x.anchor1,
                                 anchor2,     x.anchor1,  anchor2};
            if (!x.alignment.empty()) {
                s.start1 = x.alignment[0].start1;
                s.start2 =
                    contigToSequencePos(req[idx].seqData->contig, req[idx].seqData->strandNum, x.alignment[0].start2);
                addAlignedProfile(s.alignedSequences, x.alignment, alphabet, profileSeq);
                addAlignedSequence(s.alignedSequences, x.alignment, alphabet, sequence, maskedSequence);
            }
            similarities.push_back(s);
        }
    }

}

void findFinalSimilaritiesBatched(std::vector<FinalSimilarity> &similarities,
                                  std::vector<std::vector<SequenceRequest>> &allRequests,
                                  const std::vector<Profile> &profiles, const char *charVec,
                                  DPScratch &scratch) {
    for (size_t i = 0; i < profiles.size(); ++i) {
        auto &requests = allRequests[i];
        std::sort(requests.begin(), requests.end(), std::greater<>());

        int k = 0;
        std::array<SequenceRequest, simdWidth> curBatch;
        for (const auto &req : requests) {
            curBatch[k++] = req;
            if (k == simdWidth) {
                findFinalSimilarities(similarities, curBatch, profiles[i], i, charVec, scratch, simdWidth);
                k = 0;
            }
        }
        if (k > 0) {
            findFinalSimilarities(similarities, curBatch, profiles[i], i, charVec, scratch, k);
        }
    }
}

double methodOfMomentsLambda(const double *scores, int n, double meanScore) {
    double pi = 3.1415926535897932;
    double s = 0;
    for (int i = 0; i < n; ++i) {
        s += (scores[i] - meanScore) * (scores[i] - meanScore);
    }
    double variance = s / n; // apparently, method of moments doesn't use n-1
    return pi / sqrt(6 * variance);
}

double methodOfMomentsK(double meanScore, double lambda, double seqLength) {
    double euler = 0.57721566490153286;
    return exp(lambda * meanScore - euler) / seqLength;
}

double methodOfLmomentsLambda(const double *sortedScores, int n, double meanScore) {
    double s = 0;
    for (int i = 0; i < n; ++i) {
        s += i * sortedScores[i];
    }
    double d = 0.5 * n * (n - 1); // !!! avoids int overflow
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
    while (1) { // bisection method to find lambda that makes shouldBe0 = 0
        gap /= 2;
        double mid = lo + gap;
        if (mid <= lo)
            return lo;
        double z = shouldBe0(scores, n, mid);
        if ((x < 0 && z <= 0) || (x > 0 && z >= 0))
            lo = mid;
    }
}

double maximumLikelihoodK(const double *scores, int n, double lambda, double seqLength) {
    double s = 0;
    for (int i = 0; i < n; ++i) {
        s += exp(-lambda * scores[i]);
    }
    return n / (s * seqLength);
}

void methodOfMomentsGumbel(double &lambda, double &k, double &kSimple, const double *scores, int n,
                           double seqLength) {
    double meanScore = mean(scores, n);
    lambda = methodOfMomentsLambda(scores, n, meanScore);
    k = methodOfMomentsK(meanScore, lambda, seqLength);
    kSimple = methodOfMomentsK(meanScore, 1, seqLength);
}

void methodOfLmomentsGumbel(double &lambda, double &k, const double *scores, int n,
                            double seqLength) {
    double meanScore = mean(scores, n);
    lambda = methodOfLmomentsLambda(scores, n, meanScore);
    k = methodOfMomentsK(meanScore, lambda, seqLength);
}

void maximumLikelihoodGumbel(double &lambda, double &k, double &kSimple, const double *scores,
                             int n, double seqLength) {
    lambda = maximumLikelihoodLambda(scores, n);
    k = maximumLikelihoodK(scores, n, lambda, seqLength);
    kSimple = maximumLikelihoodK(scores, n, 1, seqLength);
}

void estimateGumbel(double &mmLambda, double &mmK, double &mmKsimple, double &mlLambda, double &mlK,
                    double &mlKsimple, double &lmLambda, double &lmK, double *scores, int n,
                    double seqLength) {
    std::sort(scores, scores + n);
    methodOfMomentsGumbel(mmLambda, mmK, mmKsimple, scores, n, seqLength);
    maximumLikelihoodGumbel(mlLambda, mlK, mlKsimple, scores, n, seqLength);
    methodOfLmomentsGumbel(lmLambda, lmK, scores, n, seqLength);
}

void estimateK(Profile &profile, const Float *letterFreqs, char *sequence, int sequenceLength,
               int border, int numOfSequences, int printVerbosity, DPScratch &scratch) {
    std::mt19937_64 randGen;
    int alphabetSize = profile.width - nonLetterWidth;
#ifdef ESTIMATOR_USE_RANDOM_CODONS
    std::discrete_distribution<> dist(letterFreqs,
                                      letterFreqs + alphabetSize); // TODO: no stop codons
#else
    std::discrete_distribution<> dist(letterFreqs, letterFreqs + alphabetSize);
#endif
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
        for (int j = 0; j < offset; j++) {
            sequence[j] = charToNumber[bases[distDNA(randGen)]];
        }
        for (int j = offset; j <= sequenceLength; j += 3) {
            bool shouldFS = frameshiftDist(randGen);
            if (shouldFS) {
                sequence[j] = charToNumber[bases[distDNA(randGen)]];
                j -= 2;
                continue;
            }
            int x = dist(randGen);
            auto &codons = aa2codons[alphabet[x]];
            std::discrete_distribution<> dist2(0, (int)codons.size());
            auto &xx = codons[dist2(randGen)];
            for (int k = 0; k < 3; k++) {
                if (j + k <= sequenceLength) {
                    sequence[j + k] = charToNumber[xx[k]];
                }
            }
        }
#else
        for (int j = 0; j <= sequenceLength; ++j)
            sequence[j] = dist(randGen);
#endif

        for (int j = 0; j < border; ++j)
            sequence[sequenceLength + j] = sequence[j];
        std::array<std::vector<AlignedSimilarity>, simdWidth> simsSIMD;
        std::array<std::vector<uint8_t>*, simdWidth> decoded;
        auto d = decodeSequence(sequence, sequenceLength + border, alphabet, charToNumber);
        decoded[0] = &d;
        findSimilarities(simsSIMD, profile, decoded, {-2}, scratch, 1);

        auto &sims = simsSIMD[0];
        endScores[i] = log(sims[0].probRatio);
        begScores[i] = log(sims[1].probRatio);
        midScores[i] = log(sims[2].probRatio);
        if (printVerbosity > 1) {
            std::cout << (i + 1) << "\t" << sims[0].anchor1 << "\t" << sims[0].anchor2 << "\t"
                      << log2(sims[0].probRatio) + shift << "\t" << sims[1].anchor1 << "\t"
                      << sims[1].anchor2 << "\t" << log2(sims[1].probRatio) + shift << "\t"
                      << sims[2].anchor1 << "\t" << sims[2].anchor2 << "\t"
                      << log2(sims[2].probRatio) + shift << std::endl;
        }
    }

    double MMendL, MMendK, MMendKsimple, MLendL, MLendK, MLendKsimple;
    double LMendL, LMendK;
    estimateGumbel(MMendL, MMendK, MMendKsimple, MLendL, MLendK, MLendKsimple, LMendL, LMendK,
                   endScores, numOfSequences, sequenceLength);

    double MMbegL, MMbegK, MMbegKsimple, MLbegL, MLbegK, MLbegKsimple;
    double LMbegL, LMbegK;
    estimateGumbel(MMbegL, MMbegK, MMbegKsimple, MLbegL, MLbegK, MLbegKsimple, LMbegL, LMbegK,
                   begScores, numOfSequences, sequenceLength);

    double MMmidL, MMmidK, MMmidKsimple, MLmidL, MLmidK, MLmidKsimple;
    double LMmidL, LMmidK;
    estimateGumbel(MMmidL, MMmidK, MMmidKsimple, MLmidL, MLmidK, MLmidKsimple, LMmidL, LMmidK,
                   midScores, numOfSequences, sequenceLength);

    double s = scale;

    if (printVerbosity > 1) {
        std::cout << "#\tend-\tstart-\tmid-anchored\n";

        std::cout << "#lamMM\t" << MMendL << "\t" << MMbegL << "\t" << MMmidL << "\n"

                  << "#kMM\t" << MMendK / pow(s, MMendL) << "\t" << MMbegK / pow(s, MMbegL) << "\t"
                  << MMmidK / pow(s, MMmidL) << "\n"

                  << "#kMM1\t" << MMendKsimple / scale << "\t" << MMbegKsimple / scale << "\t"
                  << MMmidKsimple / scale << "\n";

        std::cout << "#lamML\t" << MLendL << "\t" << MLbegL << "\t" << MLmidL << "\n"

                  << "#kML\t" << MLendK / pow(s, MLendL) << "\t" << MLbegK / pow(s, MLbegL) << "\t"
                  << MLmidK / pow(s, MLmidL) << "\n"

                  << "#kML1\t" << MLendKsimple / scale << "\t" << MLbegKsimple / scale << "\t"
                  << MLmidKsimple / scale << "\n";

        std::cout << "#lamLM\t" << LMendL << "\t" << LMbegL << "\t" << LMmidL << "\n"

                  << "#kLM\t" << LMendK / pow(s, LMendL) << "\t" << LMbegK / pow(s, LMbegL) << "\t"
                  << LMmidK / pow(s, LMmidL) << "\n";
    } else if (printVerbosity > 0) {
        std::cout << "# K: " << MMendKsimple / scale << " " << MMbegKsimple / scale << " "
                  << MMmidKsimple / scale << "\n";
    } else {
        std::cout << "# K: " << MMmidKsimple / scale << "\n";
    }

    static Float lambdas = 0, n = 0;
    std::cout << "# Lambda: " << MMmidL << "\n";
    n++, lambdas += MMmidL;
    std::cout << "# Avg Lambda: " << (lambdas / n) << "\n";

    profile.gumbelKendAnchored = MMendK;
    profile.gumbelKbegAnchored = MMbegK;
    profile.gumbelKmidAnchored = MMmidK;
    profile.lambda = MMmidL;
}

int intFromText(const char *text) {
    long x = strtol(text, 0, 0);
    if (x > INT_MAX || x < INT_MIN)
        return -1;
    return x;
}

double probFromText(const char *text) {
    if (*text == '*')
        return 0;
    double d = strtod(text, 0);
    return exp(-d);
}

void normalize(Float *x, int n) {
    double sum = 0;
    for (int i = 0; i < n; ++i)
        sum += x[i];
    assert(sum > 0);
    for (int i = 0; i < n; ++i)
        x[i] /= sum;
}

double meanOfLogs(const Float *x, int n) {
    double m = 1;
    for (int i = 0; i < n; ++i)
        m *= x[i];
    return log(m) / n;
}

double myMean(const Float *values, int length, int step, int meanType, Float *valuesForMedian,
              const float *tantanProbs) {
    double mean = 0;
    int n = 0;
    for (int i = 0; i < length; ++i) {
        if (tantanProbs[i] >= 0.5)
            continue;
        double v = values[i * step];
        // Geometric mean is bad for zero (or very low) probabilities
        // All letter probs in Dfam-curated_only 3.9 and Pfam-A 38.0 are > 1e-6
        if (meanType == 'G')
            mean += log(std::max(v, 1e-6)); // geometric mean
        if (meanType == 'A')
            mean += v; // arithmetic mean
        if (meanType == 'M')
            valuesForMedian[n] = v; // median
        ++n;
    }
    assert(n > 0);
    if (meanType == 'G')
        return exp(mean / n);
    if (meanType == 'A')
        return mean / n;
    std::sort(valuesForMedian, valuesForMedian + n);
    return valuesForMedian[n / 2];
}

void filterLetterProbabilities(Float *letterProbs, int length, int step, double stdDev,
                               bool keepNonvaryingTerm) {
    const double sqrt2pi = 2.5066282746310005;
    const double inv2var = 0.5 / (stdDev * stdDev);
    const int gaussianLimit = ceil(stdDev * 8); // truncate Gaussian tails
    const int alphabetSize = step - nonLetterWidth;
    std::vector<double> meanLogProbs(length);
    std::vector<double> values(length);

    for (int i = 0; i < length; ++i) {
        meanLogProbs[i] = meanOfLogs(letterProbs + i * step, alphabetSize);
    } // at each position in the profile, calculate: mean(log(letter prob))

    for (int k = 0; k < alphabetSize; ++k) {
        for (int i = 0; i < length; ++i) {
            double prob = letterProbs[i * step + k];
            values[i] = log(prob) - meanLogProbs[i]; // apply filter to this
        }

        double addItBack = keepNonvaryingTerm ? mean(values.data(), length) : 0.0;
        for (int i = 0; i < length; ++i) {
            double sum = 0;
            for (int j = -gaussianLimit; j <= gaussianLimit; ++j) {
                // xxx this treats the profile as circular (wrapping around at
                // the edges), which is rarely appropriate, but ensures no
                // change in average value:
                int x = (i + j) % length;
                if (x < 0)
                    x += length;
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

int finalizeProfile(Profile &p, char *consensusSequence, int backgroundProbsType, bool isMask,
                    double filterStdDev, bool keepNonvaryingTerm) {
    int alphabetSize = p.width - nonLetterWidth;
    std::vector<float> tantanProbs(p.length);
    std::vector<Float> valuesForMedian(p.length);
    Float *end = p.values + p.width * p.length;

    const char *alphabet = getAlphabet(alphabetSize);
    if (end[3] <= 0) {
        // set the final epsilon to the geometric mean of the other epsilons
        end[3] = myMean(p.values + p.width + 3, p.length - 1, p.width, 'G', valuesForMedian.data(),
                        tantanProbs.data());
    }

    if (filterStdDev > 0) {
        filterLetterProbabilities(p.values + 4, p.length, p.width, filterStdDev,
                                  keepNonvaryingTerm);
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
    p.bg_probs.push_back(0);
    for (int k = 4; k < 4 + alphabetSize; ++k) {
        end[k] /= sumOfMeans;
        end[k] *= (1 - BG_STOP_CODON_PROB) / aa2codons.at(alphabet[k - 4]).size();
        p.bg_probs.push_back(end[k]);
    }
    p.bg_probs.push_back(0);
    p.bg_probs.push_back(0);
    p.bg_probs.push_back(1.0 / 64.0);
    p.bg_probs.push_back(BG_STOP_CODON_PROB / 3.0);
    p.bg_probs.push_back(0); // For zero_idx padding

    for (int k = 0; k < alphabetSize + 9; ++k) {
        p.log2_bg_probs.push_back(p.bg_probs[k] > 0 ? log2(p.bg_probs[k]) : -INFINITY);
    }

    std::unordered_map<char, double> dist;
    char charToNumber[256];
    setCharToNumber(charToNumber, alphabet);
    for (auto c : std::string(alphabet)) {
        if (c == 'O')
            c = 'K';
        if (c == 'U')
            c = 'C';
        dist[c] = end[4 + charToNumber[c]] * (1 - STOP_CODON_PROB);
    }
    dist['*'] = STOP_CODON_PROB;

    p.values_v2.reserve(p.length + 5);
    for (int i = 0;; ++i) {
        p.values_v2.push_back({0});

        Float *probs = p.values + i * p.width;
        double alpha = probs[0];
        double beta = probs[1];

        double alphaFS1 = FRAMESHIFT1_MULTIPLIER;
        double alphaFS2 = FRAMESHIFT2_MULTIPLIER;
        p.values_v2.rbegin()->alpha_prime[0] = alpha * (1 - beta);
        p.values_v2.rbegin()->alpha_prime[1] = alphaFS1 * (1 - beta);
        p.values_v2.rbegin()->alpha_prime[2] = alphaFS2 * (1 - beta);

        p.values_v2.rbegin()->beta_prime[0] = beta;
        // all beta, are equal for now
        for (int j = 1; j <= 2; j++) { // perform actual insertion after frameshift
            p.values_v2.rbegin()->beta_prime[j] = 0;
        }

        double delta = probs[2];
        double epsilon = probs[3];
        if (i == p.length)
            break;

        double delta1 = probs[p.width + 2];
        double epsilon1 = probs[p.width + 3];

#ifdef ENABLE_FS_DELETE_STATES
        double deltaFS1 = FRAMESHIFT1_MULTIPLIER; // simulate delete
        double deltaFS2 = FRAMESHIFT2_MULTIPLIER;
        p.values_v2.rbegin()->delta_prime[0] = delta * (1 - epsilon1);
        p.values_v2.rbegin()->delta_prime[1] = deltaFS1 * (1 - epsilon1);
        p.values_v2.rbegin()->delta_prime[2] = deltaFS2 * (1 - epsilon1);

        p.values_v2.rbegin()->epsilon_prime = epsilon * (1 - epsilon1) / (1 - epsilon);
        p.values_v2.rbegin()->enter_match_probability = (1 - alpha - alphaFS1 - alphaFS2 - delta - deltaFS1 - deltaFS2);
#else
        p.values_v2.rbegin()->delta_prime = delta * (1 - epsilon1);
        p.values_v2.rbegin()->epsilon_prime = epsilon * (1 - epsilon1) / (1 - epsilon);
        p.values_v2.rbegin()->enter_match_probability = (1 - alpha - alphaFS1 - alphaFS2 - delta);
#endif


        double c = (1 - alpha - delta);
        if (epsilon >= 1)
            return 0;
        probs[2] = delta;
        probs[3] = epsilon;
        for (int k = 4; k < 4 + alphabetSize; ++k) {
            assert(alphabet[k - 4] != '*');
            if (tantanProbs[i] >= 0.5)
                probs[k] = (1 - BG_STOP_CODON_PROB) * end[k] / aa2codons.at(alphabet[k - 4]).size();
            double p = probs[k];
            probs[k] = ((1 - STOP_CODON_PROB) /* minus stop codon */ * p /
                        aa2codons.at(alphabet[k - 4]).size());
        }
        if (alphabetSize == 20) {
            probs[4 + 20] = probs[4 + 1]; // selenocysteine = cysteine
            probs[4 + 21] = probs[4 + 8]; // pyrrolysine = lysine
        }
        probs[4 + alphabetSize + 2] = 1.0 / 64.0; // for masked sequence letters
        probs[4 + alphabetSize + 3] = STOP_CODON_PROB / 3.0;
        probs[4 + alphabetSize + 4] = 0.0; // zero_idx padding
        if (tantanProbs[i] >= 0.5)
            consensusSequence[i] |= 32;
    }

    // extra padding
    p.values_v2.push_back({0});

    return 1;
}

int readProfiles(std::istream &in, std::vector<Profile> &profiles, std::vector<Float> &values,
                 std::vector<char> &charVec, int backgroundProbsType, bool isMask,
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
            if (word != "COMPO")
                ++state;
            break;
        case 3: {
            iss >> word;
            double MtoI = probFromText(word.c_str());
            iss >> word;
            double MtoD = probFromText(word.c_str());
            iss >> word >> word;
            double ItoI = probFromText(word.c_str());
            iss >> word >> word;
            double DtoD = probFromText(word.c_str());
            if (!iss)
                return 0;
            if (MtoI > 1 || MtoD > 1 || ItoI > 1 || DtoD > 1)
                return 0;
            values.push_back(MtoI);
            values.push_back(ItoI);
            values.push_back(MtoD);
            values.push_back(DtoD);
        }
            ++state;
            break;
        case 4:
            if (word == "//") {
                if (profile.length < 2)
                    return 0;
                values.insert(values.end(), profile.width - 4, 0.0);
                profiles.push_back(profile);
                profile.width = profile.length = 0;
                state = 0;
            } else {
                int k = 0;
                while (iss >> word && strchr(word.c_str(), '.')) { // xxx "*"?
                    double prob = probFromText(word.c_str());
                    if (prob > 1)
                        return 0;
                    values.push_back(prob);
                    ++k;
                }
                values.insert(values.end(), nonLetterWidth - 4, 0.0); // extra letters
                if (k == 0)
                    return 0;
                if (profile.width > 0 && k + nonLetterWidth != profile.width)
                    return 0;
                profile.width = k + nonLetterWidth;
                profile.length += 1;
                if (profile.length + 1 > INT_MAX / profile.width)
                    return 0;
                const Float *letterProbs = &values[values.size() - profile.width + 4];
                const Float *m = std::max_element(letterProbs, letterProbs + k);
                charVec.push_back(m - letterProbs); // consensus sequence
                state = 2;
            }
        }
    }

    Float *v = &values[0];
    for (auto &p : profiles) {
        p.values = v;
        char *consensus = &charVec[p.consensusSequenceIdx];
        if (!finalizeProfile(p, consensus, backgroundProbsType, isMask, filterStdDev,
                             keepNonvaryingTerm))
            return 0;
        v += p.width * (p.length + 1);
    }

    return state == 0;
}

Float *resizeMem(Float *v, size_t &size, int profileLength, int sequenceLength) {
    long rowSize = simdRoundUp(sequenceLength + 1) + simdLen;
    if (rowSize > LONG_MAX / (profileLength + 2)) {
        std::cerr << "too big combination of sequence and profile\n";
        return 0;
    }
    size_t s = rowSize * (profileLength + 2);
    if (s > size) {
        size = s;
        free(v);
        v = (Float *)aligned_alloc(simdLen * sizeof(Float), s * sizeof(Float));
        // this memory allocation doesn't get "free"-d at the end: that is ok!
        if (!v)
            std::cerr << "failed to allocate memory for " << s << " numbers\n";
    }
    return v;
}

void makeMaskedSequence(char *sequence, int length, int alphabetSize) {
    std::vector<float> tantanProbs(length);
    std::string seq2;
    for (int i = 0; i < length; i++) {
        char val = '\0';
        switch (sequence[i]) {
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

    calcTantanProbabilities((const unsigned char *)seq2.c_str(), length, false, tantanProbs.data());
    int mask = alphabetSize + 2;
    for (int i = 0; i < length; ++i) {
        sequence[length + i] = (tantanProbs[i] < 0.5) ? sequence[i] : mask;
    }
}

int main(int argc, char *argv[]) {
#if defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || defined(_M_X64)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
    build_standard_genetic_code(); // hack

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
  -e E, --evalue E  find similarities with E-value <= this (default: " STR(OPT_e) ")\n\
  -s S, --strand S  DNA strand: 0=reverse, 1=forward, 2=both (default: " STR(OPT_s) ")\n\
  -m M, --mask M    mask simple regions of:\n\
                    0=neither, 1=profile, 2=sequence, 3=both (default: " STR(OPT_m) ")\n\
\n\
Options for low-cut/high-pass filter on position-specific letter probabilities:\n\
  -d D, --dev D     standard deviation for Gaussian filter\n\
  -D D, --Dev D     same as above, but keep the non-varying component\n\
\n\
Options for random sequences:\n\
  -t T, --trials T  generate this many random sequences (default: " STR(OPT_t) ")\n\
  -l L, --length L  length of each random sequence (default: " STR(OPT_l) ")\n\
  -b B, --border B  add this size border to each random sequence (default: " STR(OPT_b) ")\n\
\n\
Options for background letter probabilities:\n\
  --barithmetic     arithmetic mean of position-specific probabilities\n\
  --bgeometric      geometric mean of position-specific probabilities (default)\n\
  --bmedian         median of position-specific probabilities\n\
";

    const char sOpts[] = "hVve:s:m:d:D:t:l:b:";

    static struct option lOpts[] = {{"help", no_argument, 0, 'h'},
                                    {"version", no_argument, 0, 'V'},
                                    {"verbose", no_argument, 0, 'v'},
                                    {"evalue", required_argument, 0, 'e'},
                                    {"strand", required_argument, 0, 's'},
                                    {"mask", required_argument, 0, 'm'},
                                    {"dev", required_argument, 0, 'd'},
                                    {"Dev", required_argument, 0, 'D'},
                                    {"trials", required_argument, 0, 't'},
                                    {"length", required_argument, 0, 'l'},
                                    {"border", required_argument, 0, 'b'},
                                    {"barithmetic", no_argument, 0, 'A'},
                                    {"bgeometric", no_argument, 0, 'G'},
                                    {"bmedian", no_argument, 0, 'M'},
                                    {0, 0, 0, 0}};

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
            if (evalueOpt < 0)
                return badOpt();
            break;
        case 's':
            strandOpt = intFromText(optarg);
            if (strandOpt < 0 || strandOpt > 2)
                return badOpt();
            break;
        case 'm':
            maskOpt = intFromText(optarg);
            if (maskOpt < 0 || maskOpt > 3)
                return badOpt();
            break;
        case 'd':
            filterStdDev = strtod(optarg, 0);
            // too low: discretized Gaussian problems; too high: overflow or slow
            if (filterStdDev < 2 || filterStdDev > 1000)
                return badOpt();
            break;
        case 'D':
            filterStdDev = strtod(optarg, 0);
            if (filterStdDev < 2 || filterStdDev > 1000)
                return badOpt();
            keepNonvaryingTerm = true;
            break;
        case 't':
            randomSeqNum = intFromText(optarg);
            if (randomSeqNum < 1)
                return badOpt();
            break;
        case 'l':
            randomSeqLen = intFromText(optarg);
            if (randomSeqLen < 1 || randomSeqLen > INT_MAX - 2 * simdLen)
                return badOpt();
            break;
        case 'b':
            border = intFromText(optarg);
            if (border < 0)
                return badOpt();
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

    if (filterStdDev > 0)
        maskOpt &= 2; // filtering turns off profile-masking

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
        if (!file)
            return 1;
        if (!readProfiles(in, profiles, profileValues, charVec, backgroundProbsType, maskOpt & 1,
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

    std::cout << "# DUMMER "
#include "version.hh"
                 "\n";
    std::cout << "# Bytes per floating-point number: " << sizeof(Float) << "\n";
    if (filterStdDev > 0)
        std::cout << "# Filtering position-specific letter probabilities: std dev " << filterStdDev
                  << "\n";
    std::cout << "# Background letter probabilities: "
              << (backgroundProbsType == 'A'   ? "arithmetic mean"
                  : backgroundProbsType == 'G' ? "geometric mean"
                                               : "median")
              << " of foreground probabilities\n";
    std::cout << "# Random sequences: trials=" << randomSeqNum << " length=" << randomSeqLen
              << " border=" << border << "\n";
    if (maskOpt & 1)
        std::cout << "# Masking simple regions in profiles\n";
    if (argc - optind > 1) {
        if (maskOpt & 2)
            std::cout << "# Masking simple regions in sequences\n";
        if (evalueOpt > 0)
            std::cout << "# E-value <= " << evalueOpt << "\n";
        if (strandOpt < 2)
            std::cout << "# Strand: " << (strandOpt ? "forward" : "reverse") << "\n";
    }

    int printVerbosity = (argc - optind < 2) * 2 + (evalueOpt <= 0);

    DPScratch scratch;

    for (auto &p : profiles) {
        std::cout << "\n";
        std::cout << "# Profile name: " << &charVec[p.nameIdx] << "\n";
        std::cout << "# Profile length: " << p.length << "\n";
        if (maskOpt & 1) {
            int maskCount = 0;
            const char *consensus = &charVec[p.consensusSequenceIdx];
            for (int i = 0; i < p.length; ++i)
                maskCount += (consensus[i] > 31);
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
        estimateK(p, bgProbs, &charVec[seqIdx], randomSeqLen, border, randomSeqNum, printVerbosity, scratch);
#else
        NucDist dist = *reinterpret_cast<NucDist *>(p.debug);
        Float bgProbsDNA[256] = {0};
        bgProbsDNA[charToNumber['A']] = dist.overall['A'];
        bgProbsDNA[charToNumber['C']] = dist.overall['C'];
        bgProbsDNA[charToNumber['G']] = dist.overall['G'];
        bgProbsDNA[charToNumber['T']] = dist.overall['T'];
        estimateK(p, bgProbsDNA, &charVec[seqIdx], randomSeqLen, border, randomSeqNum,
                  printVerbosity, scratch);
#endif
#endif
    }

    if (argc - optind < 2 || numOfProfiles < 1)
        return 0;
    std::cout << std::endl;

    int width = profiles[0].width;
    for (size_t i = 1; i < numOfProfiles; ++i) {
        if (profiles[i].width != width)
            width = 0;
    }
    int alphabetSize = width - nonLetterWidth;
    const char *alphabet = getAlphabet(alphabetSize);
    if (!alphabet) {
        return err("the profiles should be all protein, or all nucleotide");
    }
    char charToNumber[256];
    memset(charToNumber, 127, 256);
    memset(charToNumber, 125, ' ' + 1); // map "space characters" (<= ' ') to 125
    charToNumber['>'] = 126;            // record separator for FASTA-format sequences
    setCharToNumber(charToNumber, alphabet);
    if (alphabetSize == 4)
        setCharToNumber(charToNumber, "ACGU"); // set U = T
#ifdef PIPELINE_MODE
    strandOpt = 1;
#endif

    charVec.resize(seqIdx);
    std::vector<Sequence> sequences;
    std::vector<FinalSimilarity> similarities;
    size_t totSequenceLength = 0;

    std::ifstream file;
    std::istream &in = openFile(file, argv[optind + 1]);
    if (!file)
        return 1;
    Sequence sequence;
    Contig contig = {0, 0};
    std::vector<std::vector<SequenceRequest>> allRequests(numOfProfiles);
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
        totSequenceLength += contig.length;
        if (strandOpt == 2)
            totSequenceLength += contig.length;
        char *seq = &charVec[seqIdx];
        for (int s = 0; s < 2; ++s) {
            if (s != strandOpt) {
                if (maskOpt & 2)
                    makeMaskedSequence(seq, contig.length, alphabetSize);
                size_t strandNum = sequences.size() * 2 + s;
                std::vector<uint8_t> decoded =
                    decodeSequence(&charVec[maskedSeqIdx], contig.length, alphabet, charToNumber);
                std::shared_ptr<SequenceData> sd = std::make_shared<SequenceData>(
                    std::move(decoded), std::string(&charVec[seqIdx], contig.length),
                    std::string(&charVec[maskedSeqIdx], contig.length), contig, strandNum);

                for (size_t j = 0; j < numOfProfiles; ++j) {
                    const Profile &p = profiles[j];
#ifdef PIPELINE_MODE
                    if (!strcmp(&charVec[p.nameIdx], sequence.target_profile.c_str())) {
#endif
                        Float minProbRatio =
                            (evalueOpt > 0)
                                ? (std::pow(p.gumbelKmidAnchored * totSequenceLength / evalueOpt,
                                            1.0 / p.lambda))
                                : -1;
                        if (verbosity > 1)
                            std::cerr << "Profile: " << &charVec[p.nameIdx] << "\n";

                        allRequests[j].push_back({sd, minProbRatio});
#ifdef PIPELINE_MODE
                    }
#endif
                }
            }
#ifndef PIPELINE_MODE
            reverseComplement(seq, seq + contig.length);
#endif
        }
        charVec.resize(seqIdx);
    }

    findFinalSimilaritiesBatched(similarities, allRequests, profiles, charVec.data(), scratch);

    std::cout << "# Total sequence length: " << totSequenceLength << "\n";

    std::cout.precision(3);
    for (size_t i = 0; i < similarities.size(); ++i) {
        Profile p = profiles[similarities[i].profileNum];
        Sequence s = sequences[similarities[i].strandNum / 2];
        double k = (evalueOpt > 0) ? p.gumbelKmidAnchored
                   : (i % 3 == 0)  ? p.gumbelKendAnchored
                   : (i % 3 == 1)  ? p.gumbelKbegAnchored
                                   : p.gumbelKmidAnchored;
        double evalue = k * totSequenceLength / pow(similarities[i].probRatio, p.lambda);
        if (evalueOpt <= 0 && i % 3 == 0)
            std::cout << "\n";
        if (evalueOpt > 0 && evalue > evalueOpt)
            continue;
        printSimilarity(charVec.data(), p, s, similarities[i], evalue);
    }

    return 0;
}
