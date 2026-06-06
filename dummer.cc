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
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <functional>

#include <assert.h>
#include <ctype.h>
#include <filesystem>
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

#include <cereal/archives/binary.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/string.hpp>

#define XXH_INLINE_ALL
#include <xxhash.h>

#define OPT_e 10
#define OPT_s 2
#define OPT_m 3
#define OPT_t 1000
#define OPT_l 5000
#define OPT_b 100
#define OPT_x 0 // 0 to enable full DP mode

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

template <typename F>
inline simd_t gather_simd(F func) {
    alignas(64) Float arr[simdWidth];
    for (int k = 0; k < simdWidth; k++) arr[k] = func(k);
    return Kokkos::Experimental::simd_unchecked_load<simd_t>(arr);
}

// Memory-safe masked SIMD store to halt ghost zone zero-clobbering
inline void blend_store_simd(Float* ptr, simd_t vec, Kokkos::Experimental::simd_mask<Float> mask) {
    simd_t prev = Kokkos::Experimental::simd_unchecked_load<simd_t>(ptr);
    simd_t blended = Kokkos::Experimental::condition(mask, vec, prev);
    for (int k = 0; k < simdWidth; k++) ptr[k] = blended[k];
}

const Float STOP_CODON_PROB = 0.0005;
const Float BG_STOP_CODON_PROB = 0.046875; // 3/64

const Float INSERT1 = 0.0171;
const Float INSERT2 = 0.0018;
const Float DELETE1 = 0.0328;
const Float DELETE2 = 0.0083;

const Float BACKGROUND_FRAMESHIFT_RATE = INSERT1 + DELETE2;
const Float BACKGROUND_FRAMESHIFT_RATE_2 = INSERT2 + DELETE1;

#define TANTAN_MASK_THRESHOLD 0.5

int simdRoundUp(int x) { // lowest multiple of simdLen that is >= x
    return x - 1 - (x - 1) % simdLen + simdLen;
}

class ThreadPool {
public:
    ThreadPool(size_t numThreads) {
        for (size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back([this, i] {
                while (true) {
                    std::function<void(int)> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task(i);
                }
            });
        }
    }

    void enqueue(std::function<void(int)> task) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.push(std::move(task));
        }
        condition.notify_one();
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers)
            worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void(int)>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop = false;
};

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
    std::string name;
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

const size_t PADDING_SIZE = simdWidth + 32 /* block size */ + 3;

template <typename T, int K = 0>
class DiagonalMatrix {
public:
    std::vector<T> data;
    std::size_t padded_rows;
    std::size_t d_pad;
    std::size_t i_pad_front;

public:
    explicit DiagonalMatrix() : padded_rows(0), d_pad(PADDING_SIZE), i_pad_front(PADDING_SIZE) {}

    void resize(std::size_t r, std::size_t c) {
        // Flat 1-Row-per-Diagonal Layout with Padding
        // Align padded_rows to simdWidth for safe SIMD access
        std::size_t row_size = r + i_pad_front + PADDING_SIZE;
        padded_rows = ((row_size + simdWidth - 1) / simdWidth) * simdWidth;

        std::size_t num_diags;
        if constexpr (K > 0) {
            num_diags = K;
        } else {
            num_diags = r + c + d_pad + PADDING_SIZE;
        }

        if (r == 0 || c == 0) {
            data.clear();
            return;
        }

        std::size_t total_size = num_diags * padded_rows;
        if (data.size() < total_size) {
            data.resize(total_size);
        }
    }

    inline std::size_t get_d_index(std::ptrdiff_t d) const {
        std::size_t shifted_d = d + d_pad;
        if constexpr (K > 0) return shifted_d % K;
        return shifted_d;
    }

    inline T& operator()(std::ptrdiff_t i, std::ptrdiff_t j) {
        return data[get_d_index(i + j) * padded_rows + (i + i_pad_front)];
    }

    inline const T& operator()(std::ptrdiff_t i, std::ptrdiff_t j) const {
        return data[get_d_index(i + j) * padded_rows + (i + i_pad_front)];
    }

    inline T* data_ptr_diag_i(std::ptrdiff_t d, std::ptrdiff_t i) {
        return &data[get_d_index(d) * padded_rows + (i + i_pad_front)];
    }

    inline const T* data_ptr_diag_i(std::ptrdiff_t d, std::ptrdiff_t i) const {
        return &data[get_d_index(d) * padded_rows + (i + i_pad_front)];
    }

    inline simd_t load_simd(std::ptrdiff_t d, std::ptrdiff_t i) const {
        return Kokkos::Experimental::simd_unchecked_load<simd_t>(data_ptr_diag_i(d, i));
    }

    inline void store_simd(std::ptrdiff_t d, std::ptrdiff_t i, simd_t vec, Kokkos::Experimental::simd_mask<Float> mask) {
        blend_store_simd(data_ptr_diag_i(d, i), vec, mask);
    }

    inline T* get_diag_base_ptr(std::ptrdiff_t d) {
        return &data[get_d_index(d) * padded_rows + i_pad_front];
    }

    inline void safe_copy_to(std::ptrdiff_t d, std::ptrdiff_t i_start, std::size_t requested_count, T* dest) const {
        std::copy_n(data_ptr_diag_i(d, i_start), requested_count, dest);
    }
};


struct BlockDPScratch {
    DiagonalMatrix<Float, 0> W1, X, X_pfx;
    DiagonalMatrix<Float, 8> Y0, Z0_mat, W0;
#ifdef ENABLE_FS_INSERT_EXTENSION
    DiagonalMatrix<Float, 8> Z1_mat, Z2_mat;
#endif

    std::vector<Float> dp;
    std::vector<Float> dp_r;
    std::vector<Float> left_side;
    std::vector<Float> right_side;
    std::vector<Float> one;
    std::vector<Float> one_sfx;
    std::vector<Float> bg_codon_probs;
    std::vector<uint8_t> padded_decoded;

    std::vector<Float> best_wMid;
    std::vector<Float> best_wEnd;
    std::vector<int> best_i;

    std::vector<AlignedSimilarity> opt_profile_position;
    std::vector<bool> aligned;

    void resize(int maxProfileLength, int maxSequenceLength) {
        int r = maxProfileLength + 2;
        int c = maxSequenceLength + 8;

        W1.resize(r, c);
        X.resize(r, c);
        X_pfx.resize(r, c);

        Y0.resize(r, c);
        Z0_mat.resize(r, c);
#ifdef ENABLE_FS_INSERT_EXTENSION
        Z1_mat.resize(r, c);
        Z2_mat.resize(r, c);
#endif
        W0.resize(r, c);

        constexpr int NEG_PAD = 40; // B + 8
        constexpr int POS_PAD = 64; // 2 * B
        int pad_size = maxSequenceLength + NEG_PAD + POS_PAD;

        dp.assign(pad_size, 0);
        dp_r.assign(pad_size, 0);

        left_side.assign(pad_size, 0);
        right_side.assign(pad_size, 0);
        one.assign(pad_size, 0);
        one_sfx.assign(pad_size, 0);
        bg_codon_probs.assign(pad_size, 0);
        padded_decoded.assign(pad_size, 255);

        best_wMid.assign(c, -INFINITY);
        best_wEnd.assign(c, 0.0);
        best_i.assign(c, -1);
    }
};

struct DPScratch {
    BlockDPScratch blockScratch;
};


struct DP_Cell_v2 {
    Float metric;
    int i, j;
    bool emit = false;

    constexpr bool operator< (const DP_Cell_v2 &other) const {
        return metric < other.metric;
    }
};

void addForwardAlignment(size_t profileLength, size_t sequenceLength, std::vector<SegmentPair> &alignment, int iBeg, int jBeg,
                         BlockDPScratch &scratch) {
    constexpr int NEG_PAD = 40;
    Float* left_side_ptr = scratch.left_side.data() + NEG_PAD;

    int i = iBeg, j = jBeg;
    while (i <= profileLength && j < sequenceLength) {
        auto choice = std::max({
            DP_Cell_v2{.metric=scratch.X(i, j) + scratch.W1(i + 1, j + 3), .i=i + 1, .j=j + 3, .emit=true},
            DP_Cell_v2{.metric=scratch.W1(i + 1, j), .i=i + 1, .j=j, .emit=false},
            DP_Cell_v2{.metric=scratch.W1(i, j + 1), .i=i, .j=j + 1, .emit=false},
            DP_Cell_v2{.metric=left_side_ptr[j], .i=INT_MAX, .j=INT_MAX, .emit=false},
        });

        if (choice.emit && j >= 2) {
            addForwardMatch(alignment, i, j - 2);
        }
        i = choice.i, j = choice.j;
    }
}

void addReverseAlignment(std::vector<SegmentPair> &alignment, int iEnd, int jEnd,
                         BlockDPScratch &scratch) {
    constexpr int NEG_PAD = 40;
    Float* right_side_ptr = scratch.right_side.data() + NEG_PAD;

    int i = iEnd, j = jEnd;
    while (i >= 0 && j >= 0) {
        DP_Cell opt_succ = 0;
        if (i - 1 >= 0 && j - 3 >= 0) {
            opt_succ = scratch.X_pfx(i - 1, j - 3);
        }
        auto choice = std::max({
            DP_Cell_v2{.metric=scratch.X(i, j) + opt_succ, .i=i - 1, .j=j - 3, .emit=true},
            DP_Cell_v2{.metric=(i ? scratch.X_pfx(i - 1, j) : 0), .i=i - 1, .j=j, .emit=false},
            DP_Cell_v2{.metric=(j ? scratch.X_pfx(i, j - 1) : 0), .i=i, .j=j - 1, .emit=false},
            DP_Cell_v2{.metric=right_side_ptr[j], .i=INT_MIN, .j=INT_MIN, .emit=false},
        });

        if (choice.emit && j >= 2) {
            addReverseMatch(alignment, i, j - 2);
        }
        i = choice.i, j = choice.j;
    }
}



void addMidAnchored(size_t profileLength, size_t sequenceLength, std::vector<AlignedSimilarity> &similarities, int anchor1, int anchor2,
                    Float wBegAnchored, Float wEndAnchored, BlockDPScratch &scratch) {
    Float wMidAnchored = wEndAnchored * wBegAnchored;
    AlignedSimilarity s = {wMidAnchored / scale, anchor1, anchor2, wEndAnchored};
#ifdef ALIGN
    addForwardAlignment( profileLength, sequenceLength, s.alignment, anchor1, anchor2, scratch);
#endif
    similarities.push_back(s);
}

void finishMidAnchored(AlignedSimilarity &s, BlockDPScratch &scratch) {
    reverse(s.alignment.begin(), s.alignment.end());
#ifdef ALIGN
    addReverseAlignment(s.alignment, s.anchor1, s.anchor2, scratch);
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
    if (a == -INFINITY) return b;
    if (b == -INFINITY) return a;
    Float mx = std::max(a, b);
    return mx + std::log2(1.0 + std::exp2(std::min(a, b) - mx));
}

inline simd_t log2_sum_exp(simd_t a, simd_t b) {
    simd_t m = Kokkos::max(a, b);
    // When both a and b are -inf, a - b = NaN. Clamp to avoid NaN propagation:
    // min(abs(NaN), huge) would still be NaN, so use the fact that
    // m is -inf in that case and -inf + anything finite = -inf.
    simd_t diff = a - b;
    // Replace NaN lanes (from -inf - -inf) with 0: exp2(-0) = 1, harmless.
    // NaN comparison: NaN == NaN is false, so (diff == diff) is false for NaN lanes.
    Kokkos::Experimental::simd_mask<Float> valid = (diff == diff); // false for NaN
    simd_t x = Kokkos::Experimental::condition(valid, Kokkos::abs(diff), simd_t(0));

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

int contigToSequencePos(Contig contig, size_t strandNum, int posInContig) {
    return contig.start + strandPosition(strandNum, contig.length, posInContig);
}






#ifndef NDEBUG
#define OFFSET_ARRAY_ASSERT(cond) assert(cond)
#else
#define OFFSET_ARRAY_ASSERT(cond) ((void)0)
#endif

// A reusable wrapper for arrays allowing [-4...D_MAX-1] and [-1...I_MAX-1]
template<int D_MAX, int I_MAX>
struct OffsetArray {
    alignas(64) Float data[D_MAX + 4][I_MAX + 1];

    inline Float* operator[](int d) {
        OFFSET_ARRAY_ASSERT(-4 <= d && d < D_MAX);
        return &data[d + 4][1];
    }
    inline const Float* operator[](int d) const {
        OFFSET_ARRAY_ASSERT(-4 <= d && d < D_MAX);
        return &data[d + 4][1];
    }

    // SIMD operations
    inline simd_t load_simd(int d, int i) const {
        return Kokkos::Experimental::simd_unchecked_load<simd_t>(&(*this)[d][i]);
    }

    inline void store_simd(int d, int i, simd_t vec, Kokkos::Experimental::simd_mask<Float> mask) {
        blend_store_simd(&(*this)[d][i], vec, mask);
    }

};

void findSimilaritiesBlockDP(std::vector<AlignedSimilarity> &similarities, const Profile &profile,
                             const std::vector<uint8_t> &decoded, Float minProbRatio,
                             BlockDPScratch &scratch) {
    int maxSequenceLength = decoded.size();
    int alphabetSize = profile.width - nonLetterWidth;

    constexpr int NEG_PAD = 40;
    const int B = 32;
    scratch.resize(profile.length, maxSequenceLength);

    // Memory Alignment Fix: Inject the `zero_idx` into the unwritten padding lanes
    int zero_idx = alphabetSize + 4;
    scratch.padded_decoded.assign(scratch.padded_decoded.size(), zero_idx);

    std::copy(decoded.begin(), decoded.end(), scratch.padded_decoded.begin() + NEG_PAD);
    Float* one_ptr = scratch.one.data() + NEG_PAD;
    Float* one_sfx_ptr = scratch.one_sfx.data() + NEG_PAD;
    Float* right_side_ptr = scratch.right_side.data() + NEG_PAD;
    Float* bg_codon_probs_ptr = scratch.bg_codon_probs.data() + NEG_PAD;
    uint8_t* decoded_ptr = scratch.padded_decoded.data() + NEG_PAD;

    auto &dp = scratch.dp;
    auto &dp_r = scratch.dp_r;
    dp.assign(maxSequenceLength + 4, 0);
    dp_r.assign(maxSequenceLength + 4, 0);


    dp_r[maxSequenceLength] = 0;
    for (int i = maxSequenceLength - 1; i >= 0; i--) {
        const Float *bg_probs_ptr = profile.log2_bg_probs.data() + 4;
        Float bg_codon_emit_probs = bg_probs_ptr[decoded[i]];

        Float full_codon = log2(1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2) + bg_codon_emit_probs + dp_r[i + 3];
        Float partial_codon = log2(1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2) + log2(0.25) * (maxSequenceLength - i);
        Float t1 = (i + 2 < maxSequenceLength) ? full_codon : partial_codon;

        Float fs1 = log2(BACKGROUND_FRAMESHIFT_RATE * 0.25) + dp_r[i + 1];
        Float fs2 = log2(BACKGROUND_FRAMESHIFT_RATE_2 * 0.0625) + dp_r[i + 2];

        Float t_fs1 = (i < maxSequenceLength) ? fs1 : -INFINITY;
        Float t_fs2 = (i + 1 < maxSequenceLength) ? fs2 : -INFINITY;

        dp_r[i] = log2_sum_exp(t1, log2_sum_exp(t_fs1, t_fs2));
    }

    dp[maxSequenceLength] = 0;
    const Float *log2_bg_probs_ptr = profile.log2_bg_probs.data() + 4;
    const Float log2_1_bg_fs = log2(1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2);
    const Float log2_bg_fs_025 = log2(BACKGROUND_FRAMESHIFT_RATE * 0.25);
    const Float log2_bg_fs2_00625 = log2(BACKGROUND_FRAMESHIFT_RATE_2 * 0.0625);
    const Float log2_025 = log2(0.25);

    for (int i = 0; i < maxSequenceLength; i++) {
        Float bg_codon_emit_probs = (i >= 2) ? log2_bg_probs_ptr[decoded[i - 2]] : 0;

        Float t1;
        if (i >= 3) t1 = log2_1_bg_fs + bg_codon_emit_probs + dp[i - 3];
        else if (i == 2) t1 = log2_1_bg_fs + bg_codon_emit_probs;
        else t1 = log2_1_bg_fs + log2_025 * (i + 1);

        Float t2 = log2_bg_fs_025 + (i > 0 ? dp[i - 1] : 0);

        Float t3;
        if (i >= 2) t3 = log2_bg_fs2_00625 + dp[i - 2];
        else if (i == 1) t3 = log2_bg_fs2_00625 + 0;
        else t3 = -INFINITY;

        dp[i] = log2_sum_exp(t1, log2_sum_exp(t2, t3));
    }

    Float not_align_probs = -INFINITY;
    if (maxSequenceLength >= 3) {
        not_align_probs = log2_sum_exp(log2_sum_exp(dp[maxSequenceLength - 1], dp[maxSequenceLength - 2]), dp[maxSequenceLength - 3]);
    }

    Float distribute1 = exp2(-(not_align_probs / maxSequenceLength));
    Float distribute2 = distribute1 * distribute1;
    Float distribute3 = distribute2 * distribute1;
    distribute2 *= (0.25 * 0.25);
    distribute1 *= 0.25;

    Float invRealSeqLen = 1.0 / maxSequenceLength;
    Float nap_div_rsl = not_align_probs * invRealSeqLen;

    for (int j = 0; j < maxSequenceLength; j++) {
        Float exponent = -nap_div_rsl * (maxSequenceLength - 1.0 - j) + dp_r[j + 1];
        one_ptr[j] = exp2(exponent);
    }

    Float* left_side_ptr = scratch.left_side.data() + NEG_PAD;
    std::fill(scratch.left_side.begin(), scratch.left_side.end(), 0.0);
    std::fill(scratch.right_side.begin(), scratch.right_side.end(), 0.0);

    const Float *bg_probs_ptr = profile.bg_probs.data() + 4;
    for (int j = -4; j < maxSequenceLength + 4; j++) {
        if (j >= 0 && j < maxSequenceLength) {
            bg_codon_probs_ptr[j] = bg_probs_ptr[decoded[j]];
        }
    }

    std::fill(scratch.W1.data.begin(), scratch.W1.data.end(), 0.0);
    std::fill(scratch.Y0.data.begin(), scratch.Y0.data.end(), 0.0);
    std::fill(scratch.Z0_mat.data.begin(), scratch.Z0_mat.data.end(), 0.0);
#ifdef ENABLE_FS_INSERT_EXTENSION
    std::fill(scratch.Z1_mat.data.begin(), scratch.Z1_mat.data.end(), 0.0);
    std::fill(scratch.Z2_mat.data.begin(), scratch.Z2_mat.data.end(), 0.0);
#endif

    // BACKWARD DP BLOCK (No Tiling, Antidiagonal Sweep)
    int D_max_bwd = profile.length + maxSequenceLength + 3;
    for (int D = D_max_bwd; D >= 0; D--) {
        int min_i = std::max(0, D - maxSequenceLength + 1);
        int max_i = std::min(D, (int)profile.length);
        if (min_i > max_i) continue;

        // Fetch base pointers for required diagonals
        const Float* w1_j3_base = scratch.W1.get_diag_base_ptr(D + 4);
        const Float* y0_j_base  = scratch.Y0.get_diag_base_ptr(D + 1);
        const Float* z0_j3_base = scratch.Z0_mat.get_diag_base_ptr(D + 3);
#ifdef ENABLE_FS_DELETE_STATES
        const Float* w1_j2_base = scratch.W1.get_diag_base_ptr(D + 3);
        const Float* w1_j1_base = scratch.W1.get_diag_base_ptr(D + 2);
#endif
#ifdef ENABLE_FS_INSERT_EXTENSION
        const Float* z1_j1_base = scratch.Z1_mat.get_diag_base_ptr(D + 1);
        const Float* z2_j2_base = scratch.Z2_mat.get_diag_base_ptr(D + 2);
#endif

        // Store base pointers
        Float* w1_dst_base = scratch.W1.get_diag_base_ptr(D);
        Float* y0_dst_base = scratch.Y0.get_diag_base_ptr(D);
        Float* z0_dst_base = scratch.Z0_mat.get_diag_base_ptr(D);
#ifdef ENABLE_FS_INSERT_EXTENSION
        Float* z1_dst_base = scratch.Z1_mat.get_diag_base_ptr(D);
        Float* z2_dst_base = scratch.Z2_mat.get_diag_base_ptr(D);
#endif

        for (int i = min_i; i <= max_i; i += simdWidth) {
            simd_t is_valid = gather_simd([&](int k) -> Float {
                return (i + k > max_i) ? 0.0 : 1.0;
            });
            Kokkos::Experimental::simd_mask<Float> valid_mask = (is_valid > 0.5);

            simd_t C_enter = gather_simd([&](int k) -> Float {
                int r = i + k;
                return (r <= (int)profile.length) ? profile.values_v2[r].enter_match_probability * distribute3 : 0.0;
            });
#ifdef ENABLE_FS_DELETE_STATES
            simd_t C_delta0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].delta_prime[0] : 0.0; });
            simd_t C_delta1 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].delta_prime[1] * distribute2 : 0.0; });
            simd_t C_delta2 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].delta_prime[2] * distribute1 : 0.0; });
#else
            simd_t C_delta0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].delta_prime : 0.0; });
#endif
            simd_t C_alpha0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].alpha_prime[0] * distribute3 : 0.0; });
            simd_t C_alpha1 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].alpha_prime[1] * distribute1 : 0.0; });
            simd_t C_alpha2 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].alpha_prime[2] * distribute2 : 0.0; });
            simd_t C_beta0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].beta_prime[0] * distribute3 : 0.0; });
#ifdef ENABLE_FS_INSERT_EXTENSION
            simd_t C_beta1 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].beta_prime[1] * distribute3 : 0.0; });
            simd_t C_beta2 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].beta_prime[2] * distribute3 : 0.0; });
#endif
            simd_t C_eps0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].epsilon_prime : 0.0; });

            simd_t codon_emit_probs = gather_simd([&](int k) -> Float {
                int r = i + k;
                int c_j = D - r;
                if (r <= (int)profile.length && c_j + 1 < maxSequenceLength) {
                    return profile.values[r * profile.width + 4 + decoded_ptr[c_j + 1]];
                }
                return 0.0;
            });

            simd_t bg_codon_emit_probs = gather_simd([&](int k) -> Float {
                int r = i + k;
                int c_j = D - r;
                if (c_j + 1 < maxSequenceLength) {
                    return bg_codon_probs_ptr[c_j + 1];
                }
                return 0.0;
            });

            simd_t one_arr = gather_simd([&](int k) -> Float {
                int r = i + k;
                int c_j = D - r;
                if (c_j >= 0 && c_j < maxSequenceLength) {
                    return one_ptr[c_j];
                }
                return 0.0;
            });

            // LOAD MEMORY
            simd_t w1_i1_j3 = Kokkos::Experimental::simd_unchecked_load<simd_t>(w1_j3_base + i + 1);
            simd_t y0_i1_j  = Kokkos::Experimental::simd_unchecked_load<simd_t>(y0_j_base + i + 1);
#ifdef ENABLE_FS_DELETE_STATES
            simd_t w1_i1_j2 = Kokkos::Experimental::simd_unchecked_load<simd_t>(w1_j2_base + i + 1);
            simd_t w1_i1_j1 = Kokkos::Experimental::simd_unchecked_load<simd_t>(w1_j1_base + i + 1);
#endif
            simd_t z0_i_j3  = Kokkos::Experimental::simd_unchecked_load<simd_t>(z0_j3_base + i);
#ifdef ENABLE_FS_INSERT_EXTENSION
            simd_t z1_i_j1  = Kokkos::Experimental::simd_unchecked_load<simd_t>(z1_j1_base + i);
            simd_t z2_i_j2  = Kokkos::Experimental::simd_unchecked_load<simd_t>(z2_j2_base + i);
#endif

            // MATH
            simd_t w_val = w1_i1_j3 * codon_emit_probs * C_enter +
                          y0_i1_j * C_delta0 +
#ifdef ENABLE_FS_DELETE_STATES
                          w1_i1_j2 * C_delta1 +
                          w1_i1_j1 * C_delta2 +
#endif
                          z0_i_j3 * bg_codon_emit_probs * C_alpha0
#ifdef ENABLE_FS_INSERT_EXTENSION
                          + z1_i_j1 * C_alpha1 +
                          z2_i_j2 * C_alpha2
#endif
                          + one_arr * simd_t(scale);

            simd_t y0_val = C_eps0 * y0_i1_j + w_val;
            simd_t z0_future = z0_i_j3 * bg_codon_emit_probs;
            simd_t z0_val = C_beta0 * z0_future + w_val;
#ifdef ENABLE_FS_INSERT_EXTENSION
            simd_t z1_val = C_beta1 * z0_future + w_val;
            simd_t z2_val = C_beta2 * z0_future + w_val;
#endif

            // Store w_val to right_side scalar
            for (int k = 0; k < simdWidth; k++) {
                int r = i + k;
                if (r > max_i) break;
                int c_j = D - r;
                if (r > (int)profile.length || c_j >= maxSequenceLength) continue;
                right_side_ptr[c_j] += w_val[k];
            }

            // STORE MEMORY
            blend_store_simd(w1_dst_base + i, w_val, valid_mask);
            blend_store_simd(y0_dst_base + i, y0_val, valid_mask);
            blend_store_simd(z0_dst_base + i, z0_val, valid_mask);
#ifdef ENABLE_FS_INSERT_EXTENSION
            blend_store_simd(z1_dst_base + i, z1_val, valid_mask);
            blend_store_simd(z2_dst_base + i, z2_val, valid_mask);
#endif
        }
    }

    std::fill(scratch.Y0.data.begin(), scratch.Y0.data.end(), 0.0);
    std::fill(scratch.Z0_mat.data.begin(), scratch.Z0_mat.data.end(), 0.0);
#ifdef ENABLE_FS_INSERT_EXTENSION
    std::fill(scratch.Z1_mat.data.begin(), scratch.Z1_mat.data.end(), 0.0);
    std::fill(scratch.Z2_mat.data.begin(), scratch.Z2_mat.data.end(), 0.0);
#endif
    std::fill(scratch.W0.data.begin(), scratch.W0.data.end(), 0.0);
    std::fill(scratch.X.data.begin(), scratch.X.data.end(), 0.0);
    std::fill(scratch.X_pfx.data.begin(), scratch.X_pfx.data.end(), 0.0);
    std::fill(scratch.best_wMid.begin(), scratch.best_wMid.end(), -INFINITY);
    std::fill(scratch.best_wEnd.begin(), scratch.best_wEnd.end(), 0.0);
    std::fill(scratch.best_i.begin(), scratch.best_i.end(), -1);

    for (int j = 0; j < maxSequenceLength; j++) {
        Float exponent = -nap_div_rsl * (j + 1) + dp[j];
        one_sfx_ptr[j] = exp2(exponent);
    }

    for (int j = maxSequenceLength - 1; j >= 0; j--) {
        Float bg_codon_emit_probs = (j >= 2) ? bg_probs_ptr[decoded[j - 2]] : 0;

        if (j - 3 >= 0) right_side_ptr[j - 3] += (1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2) * bg_codon_emit_probs * distribute3 * right_side_ptr[j];
        if (j - 1 >= 0) right_side_ptr[j - 1] += (BACKGROUND_FRAMESHIFT_RATE * 0.25) * distribute1 * right_side_ptr[j];
        if (j - 2 >= 0) right_side_ptr[j - 2] += (BACKGROUND_FRAMESHIFT_RATE_2 * 0.0625) * distribute2 * right_side_ptr[j];
        right_side_ptr[j] *= one_sfx_ptr[j];
        right_side_ptr[j] = std::max(right_side_ptr[j], (Float)0.0);
    }
    for (int j = 1; j < maxSequenceLength; j++) {
        right_side_ptr[j] += right_side_ptr[j - 1];
    }

    // FORWARD DP BLOCK (No Tiling, Antidiagonal Sweep)
    int D_max_fwd = (int)profile.length + maxSequenceLength;
    for (int D = 0; D <= D_max_fwd; D++) {
        int min_i = std::max(0, D - maxSequenceLength + 1);
        int max_i = std::min(D, (int)profile.length);
        if (min_i > max_i) continue;

        // Fetch base pointers for required diagonals
        const Float* w1_base = scratch.W0.get_diag_base_ptr(D - 1);
        const Float* w2_base = scratch.W0.get_diag_base_ptr(D - 2);
        const Float* w3_base = scratch.W0.get_diag_base_ptr(D - 3);
        const Float* w1_bkwd_base = scratch.W1.get_diag_base_ptr(D + 1);
        const Float* z0_prev_base = scratch.Z0_mat.get_diag_base_ptr(D - 3);
#ifdef ENABLE_FS_INSERT_EXTENSION
        const Float* z1_prev_base = scratch.Z1_mat.get_diag_base_ptr(D - 3);
        const Float* z2_prev_base = scratch.Z2_mat.get_diag_base_ptr(D - 3);
#endif
        const Float* w0_prev_base = scratch.W0.get_diag_base_ptr(D);
        const Float* y0_prev_base = scratch.Y0.get_diag_base_ptr(D - 1);
        const Float* wBegAnchored_base = scratch.W1.get_diag_base_ptr(D);

        // Store base pointers
        Float* x_dst_base = scratch.X.get_diag_base_ptr(D);
        Float* z0_dst_base = scratch.Z0_mat.get_diag_base_ptr(D);
#ifdef ENABLE_FS_INSERT_EXTENSION
        Float* z1_dst_base = scratch.Z1_mat.get_diag_base_ptr(D);
        Float* z2_dst_base = scratch.Z2_mat.get_diag_base_ptr(D);
#endif
        Float* w0_dst_base = scratch.W0.get_diag_base_ptr(D);
        Float* y0_dst_base = scratch.Y0.get_diag_base_ptr(D);
        Float* w0_next_dst_base = scratch.W0.get_diag_base_ptr(D + 1);

        for (int i = min_i; i <= max_i; i += simdWidth) {
            simd_t is_valid = gather_simd([&](int k) -> Float {
                return (i + k > max_i) ? 0.0 : 1.0;
            });
            Kokkos::Experimental::simd_mask<Float> valid_mask = (is_valid > 0.5);

            simd_t C_enter = gather_simd([&](int k) -> Float {
                int r = i + k;
                return (r <= (int)profile.length) ? profile.values_v2[r].enter_match_probability * distribute3 : 0.0;
            });
            simd_t C_alpha0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].alpha_prime[0] : 0.0; });
            simd_t C_alpha1 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].alpha_prime[1] * distribute1 : 0.0; });
            simd_t C_alpha2 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].alpha_prime[2] * distribute2 : 0.0; });
            simd_t C_beta0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].beta_prime[0] : 0.0; });
#ifdef ENABLE_FS_INSERT_EXTENSION
            simd_t C_beta1 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].beta_prime[1] : 0.0; });
            simd_t C_beta2 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].beta_prime[2] : 0.0; });
#endif
#ifdef ENABLE_FS_DELETE_STATES
            simd_t C_delta0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].delta_prime[0] : 0.0; });
            simd_t C_delta1 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].delta_prime[1] * distribute2 : 0.0; });
            simd_t C_delta2 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].delta_prime[2] * distribute1 : 0.0; });
#else
            simd_t C_delta0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].delta_prime : 0.0; });
#endif
            simd_t C_eps0 = gather_simd([&](int k) -> Float { int r = i + k; return (r <= (int)profile.length) ? profile.values_v2[r].epsilon_prime : 0.0; });

            simd_t codon_emit_probs = gather_simd([&](int k) -> Float {
                int r = i + k;
                int c_j = D - r;
                if (r <= (int)profile.length && c_j >= 2 && c_j < maxSequenceLength) {
                    return profile.values[r * profile.width + 4 + decoded_ptr[c_j - 2]];
                }
                return 0.0;
            });

            simd_t bg_codon_emit_probs = gather_simd([&](int k) -> Float {
                int r = i + k;
                int c_j = D - r;
                if (c_j >= 2 && c_j < maxSequenceLength) {
                    return bg_codon_probs_ptr[c_j - 2];
                }
                return 0.0;
            });

            simd_t one_arr = gather_simd([&](int k) -> Float {
                int r = i + k;
                int c_j = D - r;
                if (c_j >= 0 && c_j < maxSequenceLength) {
                    return one_sfx_ptr[c_j];
                }
                return 0.0;
            });

            // Optimized LOAD MEMORY FROM DiagonalMatrix
            simd_t w1 = Kokkos::Experimental::simd_unchecked_load<simd_t>(w1_base + i);
            simd_t w2 = Kokkos::Experimental::simd_unchecked_load<simd_t>(w2_base + i);
            simd_t w3 = Kokkos::Experimental::simd_unchecked_load<simd_t>(w3_base + i);

            // Initial conditions for borders
            simd_t j_eq_0 = gather_simd([&](int k) -> Float { return (D - (i + k) == 0) ? 1.0 : 0.0; });
            simd_t j_eq_1 = gather_simd([&](int k) -> Float { return (D - (i + k) == 1) ? 1.0 : 0.0; });
            simd_t j_eq_2 = gather_simd([&](int k) -> Float { return (D - (i + k) == 2) ? 1.0 : 0.0; });

            w1 = Kokkos::Experimental::condition(j_eq_0 > 0.5, simd_t(scale), w1);
            w2 = Kokkos::Experimental::condition(j_eq_1 > 0.5, simd_t(scale), w2);
            w3 = Kokkos::Experimental::condition(j_eq_2 > 0.5, simd_t(scale), w3);

            simd_t X_ij = C_enter * codon_emit_probs * w3;
            simd_t w1_bkwd = Kokkos::Experimental::simd_unchecked_load<simd_t>(w1_bkwd_base + i + 1);
            simd_t X_ij_EV = X_ij * w1_bkwd * simd_t(invScale);

            simd_t z0_prev = Kokkos::Experimental::simd_unchecked_load<simd_t>(z0_prev_base + i);
#ifdef ENABLE_FS_INSERT_EXTENSION
            simd_t z1_prev = Kokkos::Experimental::simd_unchecked_load<simd_t>(z1_prev_base + i);
            simd_t z2_prev = Kokkos::Experimental::simd_unchecked_load<simd_t>(z2_prev_base + i);
#endif

#ifdef ENABLE_FS_INSERT_EXTENSION
            simd_t z0_val = bg_codon_emit_probs * simd_t(distribute3) * (C_alpha0 * w3 + C_beta0 * z0_prev + C_beta1 * z1_prev + C_beta2 * z2_prev);
#else
            simd_t z0_val = bg_codon_emit_probs * simd_t(distribute3) * (C_alpha0 * w3 + C_beta0 * z0_prev);
#endif
            simd_t z1_val = C_alpha1 * w1;
            simd_t z2_val = C_alpha2 * w2;

            simd_t w0_prev = Kokkos::Experimental::simd_unchecked_load<simd_t>(w0_prev_base + i);
            simd_t w0_val = w0_prev + z0_val + z1_val + z2_val + one_arr * simd_t(scale);

            simd_t y0_prev = Kokkos::Experimental::simd_unchecked_load<simd_t>(y0_prev_base + i - 1);
#ifdef ENABLE_FS_DELETE_STATES
            simd_t y0_val = C_delta0 * w0_val + C_eps0 * y0_prev;
            simd_t w0_next = X_ij + y0_val + C_delta1 * w2 + C_delta2 * w1;
#else
            simd_t y0_val = C_delta0 * w0_val + C_eps0 * y0_prev;
            simd_t w0_next = X_ij + y0_val;
#endif

            simd_t wBegAnchored = Kokkos::Experimental::simd_unchecked_load<simd_t>(wBegAnchored_base + i);
            simd_t wMidAnchored = w0_val * wBegAnchored * simd_t(invScale);

            // Scatter to global
            for (int k = 0; k < simdWidth; k++) {
                int r = i + k;
                if (r > max_i) break;
                int c_j = D - r;
                if (r > (int)profile.length || c_j >= maxSequenceLength) continue;

                Float w_mid = wMidAnchored[k];
                if (w_mid > scratch.best_wMid[c_j]) {
                    scratch.best_wMid[c_j] = w_mid;
                    scratch.best_wEnd[c_j] = w0_val[k];
                    scratch.best_i[c_j] = r;
                }
                left_side_ptr[c_j] += w0_val[k];
            }

            // Optimized STORE MEMORY
            blend_store_simd(x_dst_base + i, X_ij_EV, valid_mask);
            blend_store_simd(z0_dst_base + i, z0_val, valid_mask);
#ifdef ENABLE_FS_INSERT_EXTENSION
            blend_store_simd(z1_dst_base + i, z1_val, valid_mask);
            blend_store_simd(z2_dst_base + i, z2_val, valid_mask);
#endif
            blend_store_simd(w0_dst_base + i, w0_val, valid_mask);
            blend_store_simd(y0_dst_base + i, y0_val, valid_mask);

            blend_store_simd(w0_next_dst_base + i + 1, w0_next, valid_mask);
        }
    }

    for (int i = 0; i <= (int)profile.length; i++) {
        Float pfx_prev = 0.0;
        for (int j = 0; j < maxSequenceLength; j++) {
            Float opt_succ = 0.0;
            if (i - 1 >= 0 && j - 3 >= 0) {
                opt_succ = scratch.X_pfx(i - 1, j - 3);
            }

            Float pfx_mx = 0.0;
            if (i - 1 >= 0) {
                pfx_mx = scratch.X_pfx(i - 1, j);
            }
            if (j - 1 >= 0) {
                pfx_mx = std::max(pfx_mx, pfx_prev);
            }

            pfx_mx = std::max(pfx_mx, scratch.X(i, j) + opt_succ);
            pfx_mx = std::max(pfx_mx, right_side_ptr[j]);

            scratch.X_pfx(i, j) = pfx_mx;
            pfx_prev = pfx_mx;
        }
    }
    // =====================================================================

    for (int j = 0; j < maxSequenceLength; j++) {
        Float bg_codon_emit_probs = bg_codon_probs_ptr[j + 1];
        if (j + 3 < maxSequenceLength) left_side_ptr[j + 3] += (1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2) * bg_codon_emit_probs * distribute3 * left_side_ptr[j];
        if (j + 1 < maxSequenceLength) left_side_ptr[j + 1] += (BACKGROUND_FRAMESHIFT_RATE * 0.25) * distribute1 * left_side_ptr[j];
        if (j + 2 < maxSequenceLength) left_side_ptr[j + 2] += (BACKGROUND_FRAMESHIFT_RATE_2 * 0.0625) * distribute2 * left_side_ptr[j];
        left_side_ptr[j] *= one_ptr[j];
        left_side_ptr[j] = std::max(left_side_ptr[j], (Float)0.0);
    }
    for (int j = maxSequenceLength - 2; j >= 0; j--) {
        left_side_ptr[j] += left_side_ptr[j + 1];
    }

    for (int i = (int)profile.length; i >= 0; i--) {
        for (int j = maxSequenceLength - 1; j >= 0; j--) {
            Float opt_succ = (i + 1 <= profile.length && j + 3 < maxSequenceLength) ? scratch.W1(i + 1, j + 3) : 0;
            Float opt_down = (i + 1 <= profile.length) ? scratch.W1(i + 1, j) : 0;
            Float opt_right = (j + 1 < maxSequenceLength) ? scratch.W1(i, j + 1) : 0;

            Float opt = std::max(opt_down, opt_right);
            opt = std::max(opt, left_side_ptr[j]);
            opt = std::max(opt, scratch.X(i, j) + opt_succ);
            scratch.W1(i, j) = opt;
        }
    }

    scratch.opt_profile_position.assign(maxSequenceLength, AlignedSimilarity(-INFINITY));
    for (int j = 0; j < maxSequenceLength; j++) {
        Float best_prob = scratch.best_wMid[j];
        if (best_prob > -INFINITY) {
            scratch.opt_profile_position[j] = { best_prob, scratch.best_i[j], j, scratch.best_wEnd[j] };
        }
    }

    if (minProbRatio >= 0) {
        std::sort(scratch.opt_profile_position.begin(), scratch.opt_profile_position.end(), std::greater<>());
        std::cout << "# x-drop disabled (computed all cells)" << std::endl;
        scratch.aligned.assign(maxSequenceLength, false);
        for (auto &aligned_similarity : scratch.opt_profile_position) {
            if (aligned_similarity.probRatio >= minProbRatio &&
                !scratch.aligned[aligned_similarity.anchor2]) {
                addMidAnchored(profile.length, maxSequenceLength, similarities, aligned_similarity.anchor1,
                                    aligned_similarity.anchor2,
                                    aligned_similarity.probRatio * scale /
                                        aligned_similarity.wEndAnchored,
                                    aligned_similarity.wEndAnchored, scratch);
                auto &x = similarities.back();
                finishMidAnchored(x, scratch);

                int startIdx = std::max(aligned_similarity.anchor2 - 12 * profile.length, 0);
                int endIdx = std::min(aligned_similarity.anchor2 + 12 * profile.length, maxSequenceLength);
                std::fill(scratch.aligned.begin() + startIdx, scratch.aligned.begin() + endIdx, true);
            }
        }
    } else {
        auto sel = *std::max_element(scratch.opt_profile_position.begin(), scratch.opt_profile_position.end());
        AlignedSimilarity b = sel;
        b.probRatio = 0;
        similarities.push_back(b);
        similarities.push_back(b);
        similarities.push_back(sel);
    }
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
    for (int i = 0; i < activeCount; i++) {
        findSimilaritiesBlockDP(sims[i], profile, *decoded[i], minProbRatio[i], scratch.blockScratch);
    }

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
                                  ThreadPool &threadPool, std::vector<DPScratch> &threadScratches) {
    for (size_t i = 0; i < profiles.size(); ++i) {
        std::sort(allRequests[i].begin(), allRequests[i].end(), std::greater<>());
    }

    struct BatchJob {
        size_t profileIdx;
        size_t startRequestIdx;
        int activeCount;
    };

    std::vector<BatchJob> jobs;
    for (size_t i = 0; i < profiles.size(); ++i) {
        const auto &requests = allRequests[i];
        size_t n = requests.size();
        for (size_t start = 0; start < n; start += simdWidth) {
            int count = std::min(static_cast<size_t>(simdWidth), n - start);
            jobs.push_back({i, start, count});
        }
    }

    if (jobs.empty()) return;

    std::vector<std::vector<FinalSimilarity>> jobSimilarities(jobs.size());
    std::atomic<size_t> completedJobs(0);
    std::mutex mtx;
    std::condition_variable cv;

    for (size_t jobIdx = 0; jobIdx < jobs.size(); ++jobIdx) {
        threadPool.enqueue([&, jobIdx](int threadId) {
            DPScratch &threadScratch = threadScratches[threadId];
            const auto &job = jobs[jobIdx];
            std::array<SequenceRequest, simdWidth> curBatch;
            const auto &requests = allRequests[job.profileIdx];
            for (int k = 0; k < job.activeCount; ++k) {
                curBatch[k] = requests[job.startRequestIdx + k];
            }

            findFinalSimilarities(jobSimilarities[jobIdx], curBatch,
                                  profiles[job.profileIdx], job.profileIdx,
                                  charVec, threadScratch, job.activeCount);

            if (++completedJobs == jobs.size()) {
                std::lock_guard<std::mutex> lock(mtx);
                cv.notify_one();
            }
        });
    }

    if (!jobs.empty()) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]{ return completedJobs == jobs.size(); });
    }

    size_t totalSimilarities = 0;
    for (const auto &jobSims : jobSimilarities) {
        totalSimilarities += jobSims.size();
    }
    similarities.reserve(totalSimilarities);
    for (const auto &jobSims : jobSimilarities) {
        similarities.insert(similarities.end(), jobSims.begin(), jobSims.end());
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

static std::mutex g_cout_mutex;

class Hash128 {
public:
    Hash128() {
        XXH3_128bits_reset(state);
    }
    ~Hash128() {
        XXH3_freeState(state);
    }

    void add(const void* data, size_t size) {
        XXH3_128bits_update(state, data, size);
    }

    void add(const std::string& str) {
        add(str.data(), str.size());
    }

    template <typename T>
    requires std::integral<T> || std::floating_point<T>
    void add(const T& val) {
        add(&val, sizeof(T));
    }

    XXH128_hash_t hash() {
        return XXH3_128bits_digest(state);
    }

    std::string to_string() {
        auto res = hash();

        std::stringstream ss;
        // Format high and low 64-bit parts as 16-character padded hex strings
        ss << std::hex << std::setfill('0')
           << std::setw(16) << res.high64
           << std::setw(16) << res.low64;
        return ss.str();
    }

private:
    XXH3_state_t* const state = XXH3_createState();
};

std::string getBinaryHash() {
    std::ifstream file("/proc/self/exe", std::ios::binary);
    assert(file);

    Hash128 h;
    char buffer[65536];
    while (file.read(buffer, sizeof(buffer))) {
        h.add(buffer, file.gcount());
    }
    h.add(buffer, file.gcount());

    auto hash = h.to_string();
    return hash;
}

struct CacheEntry {
    double MMendL, MMbegL, MMmidL;
    double MMendK, MMbegK, MMmidK;
    double MMendKsimple, MMbegKsimple, MMmidKsimple;
    double MLendL, MLbegL, MLmidL;
    double MLendK, MLbegK, MLmidK;
    double MLendKsimple, MLbegKsimple, MLmidKsimple;
    double LMendL, LMbegL, LMmidL;
    double LMendK, LMbegK, LMmidK;

    template<class Archive>
    void serialize(Archive& archive) {
        archive(
            MMendL, MMbegL, MMmidL,
            MMendK, MMbegK, MMmidK,
            MMendKsimple, MMbegKsimple, MMmidKsimple,
            MLendL, MLbegL, MLmidL,
            MLendK, MLbegK, MLmidK,
            MLendKsimple, MLbegKsimple, MLmidKsimple,
            LMendL, LMbegL, LMmidL,
            LMendK, LMbegK, LMmidK
        );
    }
};

// vibe-coded cache
class ProfileCache {
public:
    ~ProfileCache() {
        save();
    }

    std::string computeCacheKey(const Profile &profile, const Float *letterFreqs,
                                    int sequenceLength, int border, int numOfSequences) {
        Hash128 h;
        h.add(binaryHash);
        h.add(profile.name);
        h.add(profile.width);
        h.add(profile.length);
        h.add(profile.values, profile.width * (profile.length + 1) * sizeof(Float));
        h.add(letterFreqs, (profile.width - nonLetterWidth) * sizeof(Float));
        h.add(sequenceLength);
        h.add(border);
        h.add(numOfSequences);
        h.add(INSERT1);
        h.add(INSERT2);
        h.add(DELETE1);
        h.add(DELETE2);
        h.add(BACKGROUND_FRAMESHIFT_RATE);
        h.add(BACKGROUND_FRAMESHIFT_RATE_2);
        h.add(STOP_CODON_PROB);
        h.add(BG_STOP_CODON_PROB);
        h.add(TANTAN_MASK_THRESHOLD);

        return h.to_string();
    }
private:
    std::unordered_map<std::string, CacheEntry> entries;
    //std::unordered_set<std::string> read_entries;
    std::filesystem::path cacheFilePath;
    std::string binaryHash = getBinaryHash();
    std::mutex cacheMutex;
    bool loaded = false;

    static std::filesystem::path getCacheFilePath() {
        std::filesystem::path cacheDir;

        if (const char* xdgCache = std::getenv("XDG_CACHE_HOME"); xdgCache && *xdgCache) {
            cacheDir = std::filesystem::path(xdgCache);
        } else if (const char* home = std::getenv("HOME")) {
            cacheDir = std::filesystem::path(home) / ".cache";
        } else {
            cacheDir = std::filesystem::current_path();
        }

        std::filesystem::path appCacheDir = cacheDir / "dummer";
        if (!std::filesystem::exists(appCacheDir)) {
            std::error_code ec;
            std::filesystem::create_directories(appCacheDir, ec);
            if (ec) {
                std::cerr << "# Warning: Could not create cache directory: " << ec.message() << "\n";
                appCacheDir = std::filesystem::current_path();
            }
        }

        return appCacheDir / "cache.bin";
    }

    void load() {
        if (loaded) return;

        cacheFilePath = getCacheFilePath();
        std::cout << "# Cache file: " << cacheFilePath << std::endl;

        for (int tries = 0; tries < 5; tries++) {
            if (std::filesystem::exists(cacheFilePath)) {
                try {
                    std::ifstream in(cacheFilePath, std::ios::binary);
                    if (in.is_open() && in.peek() != std::ifstream::traits_type::eof()) {
                        cereal::BinaryInputArchive archive(in);
                        archive(entries);
                        break;
                    }
                } catch (const std::exception& e) {
                    std::cerr << "# Cache Load Error: " << e.what() << ".\n";
                    std::error_code delete_ec;
                    std::filesystem::remove(cacheFilePath, delete_ec);
                    entries.clear();
                }
            }
        }

        loaded = true;
    }

public:
    bool lookup(const Profile &profile, const Float *letterFreqs, int sequenceLength,
                int border, int numOfSequences, CacheEntry &outEntry) {
        std::scoped_lock lock(cacheMutex);
        load();

        std::string key = computeCacheKey(profile, letterFreqs, sequenceLength, border, numOfSequences);
        if (auto it = entries.find(key); it != entries.end()) {
            outEntry = it->second;
            return true;
        }
        return false;
    }

    void store(const Profile &profile, const Float *letterFreqs, int sequenceLength,
              int border, int numOfSequences, const CacheEntry &entry) {
        std::scoped_lock lock(cacheMutex);
        load();

        std::string key = computeCacheKey(profile, letterFreqs, sequenceLength, border, numOfSequences);
        entries[key] = entry;
    }

    void save() {
        std::scoped_lock lock(cacheMutex);
        load();

        // std::erase_if(entries, [&](const auto& pair) {
        //     return read_entries.find(pair.first) == read_entries.end();
        // });

        try {
            std::ofstream out(cacheFilePath, std::ios::binary | std::ios::trunc);
            if (out) {
                cereal::BinaryOutputArchive archive(out);
                archive(entries);
            }
        } catch (const std::exception& e) {
            std::cerr << "# Cache Save Error: " << e.what() << std::endl;
        }
    }
} cache;

void estimateK(Profile &profile, const Float *letterFreqs, char *sequence, int sequenceLength,
               int border, int numOfSequences, int printVerbosity, ThreadPool &threadPool, std::vector<DPScratch> &threadScratches) {
    Float estimateK_lambdas = 0;
    Float estimateK_n = 0;

    CacheEntry entry;
    if (cache.lookup(profile, letterFreqs, sequenceLength, border, numOfSequences, entry)) {
        if (printVerbosity > 1) {
            std::cout << "# Warning: using cached results\n";
        }

        profile.gumbelKendAnchored = entry.MMendK;
        profile.gumbelKbegAnchored = entry.MMbegK;
        profile.gumbelKmidAnchored = entry.MMmidK;
        profile.lambda = entry.MMmidL;
    } else {
        int alphabetSize = profile.width - nonLetterWidth;


#ifdef ESTIMATOR_USE_RANDOM_CODONS

        std::vector<double> aaFreqs;
        Float sum = 0;
        if (alphabetSize > 4) {
            aaFreqs.resize(alphabetSize + 1);
            for (int k = 0; k < alphabetSize; ++k) {
                int n = aa2codons.at(getAlphabet(alphabetSize)[k]).size();
                aaFreqs[k] = letterFreqs[k] * n;
                sum += aaFreqs[k];
            }
            aaFreqs[alphabetSize] = BG_STOP_CODON_PROB;
            sum += aaFreqs[alphabetSize];
        } else {
            aaFreqs.assign(letterFreqs, letterFreqs + alphabetSize);
        }

        std::cout << "# sum is " << sum << std::endl;
        std::discrete_distribution<> dist(aaFreqs.begin(), aaFreqs.end());
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

        int numBatches = (numOfSequences + simdWidth - 1) / simdWidth;
        std::vector<std::array<std::vector<char>, simdWidth>> threadLocalSeqs(threadScratches.size());
        for (auto &localSeqs : threadLocalSeqs) {
            for (int lane = 0; lane < simdWidth; ++lane) {
                localSeqs[lane].resize(sequenceLength + border + 16);
            }
        }

        std::atomic<int> completedBatches(0);
        std::mutex mtx;
        std::condition_variable cv;

        for (int batchIdx = 0; batchIdx < numBatches; ++batchIdx) {
            threadPool.enqueue([&, batchIdx](int threadId) {
                DPScratch &threadScratch = threadScratches[threadId];
                auto &localSeqs = threadLocalSeqs[threadId];

                int start = batchIdx * simdWidth;
                int activeCount = std::min(static_cast<int>(simdWidth), numOfSequences - start);

                std::array<std::vector<uint8_t>, simdWidth> decoded;
                std::array<std::vector<uint8_t>*, simdWidth> decodedPtrs = {};
                std::array<Float, simdWidth> minProbRatio;
                minProbRatio.fill(-2.0f);

                for (int lane = 0; lane < activeCount; ++lane) {
                    int trialIdx = start + lane;
                    // Core-independent deterministic seeding based on trial index
                    std::mt19937_64 trialRandGen(5489 + trialIdx);
                    char *seqBuf = localSeqs[lane].data();

#ifdef ESTIMATOR_USE_RANDOM_CODONS
                    const char bases[] = {'A', 'C', 'G', 'T'};
                    std::uniform_int_distribution<int> distDNA(0, 3);
                    std::uniform_int_distribution<int> distOffset(0, 2);

                    int offset = distOffset(trialRandGen);
                    for (int j = 0; j < offset; j++) {
                        seqBuf[j] = charToNumber[bases[distDNA(trialRandGen)]];
                    }
                    for (int j = offset; j <= sequenceLength; j += 3) {
                        double r = std::generate_canonical<double, 10>(trialRandGen);
                        if (r < BACKGROUND_FRAMESHIFT_RATE) {
                            if (j <= sequenceLength) {
                                seqBuf[j] = charToNumber[bases[distDNA(trialRandGen)]];
                            }
                            j -= 2;
                            continue;
                        } else if (r < BACKGROUND_FRAMESHIFT_RATE + BACKGROUND_FRAMESHIFT_RATE_2) {
                            if (j <= sequenceLength) {
                                seqBuf[j] = charToNumber[bases[distDNA(trialRandGen)]];
                            }
                            if (j + 1 <= sequenceLength) {
                                seqBuf[j + 1] = charToNumber[bases[distDNA(trialRandGen)]];
                            }
                            j -= 1;
                            continue;
                        }
                        int x = dist(trialRandGen);
                        char aa = (x < alphabetSize) ? alphabet[x] : '*';
                        const auto &codons = aa2codons.at(aa);
                        std::uniform_int_distribution<int> dist2(0, (int)codons.size() - 1);
                        const auto &xx = codons[dist2(trialRandGen)];
                        for (int k = 0; k < 3; k++) {
                            if (j + k <= sequenceLength) {
                                seqBuf[j + k] = charToNumber[xx[k]];
                            }
                        }
                    }
#else
                    for (int j = 0; j <= sequenceLength; ++j)
                        seqBuf[j] = dist(trialRandGen);
#endif

                    for (int j = 0; j < border; ++j)
                        seqBuf[sequenceLength + j] = seqBuf[j];

                    decoded[lane] = decodeSequence(seqBuf, sequenceLength + border, alphabet, charToNumber);
                    decodedPtrs[lane] = &decoded[lane];
                }

                std::array<std::vector<AlignedSimilarity>, simdWidth> simsSIMD;
                for (int lane = 0; lane < activeCount; ++lane) {
                    findSimilaritiesBlockDP(simsSIMD[lane], profile, *decodedPtrs[lane], minProbRatio[lane], threadScratch.blockScratch);
                }

                for (int lane = 0; lane < activeCount; ++lane) {
                    int trialIdx = start + lane;
                    const auto &sims = simsSIMD[lane];
                    endScores[trialIdx] = log(sims[0].probRatio);
                    begScores[trialIdx] = log(sims[1].probRatio);
                    midScores[trialIdx] = log(sims[2].probRatio);

                    if (printVerbosity > 1) {
                        std::lock_guard<std::mutex> lock(g_cout_mutex);
                        std::cout << (trialIdx + 1) << "\t" << sims[0].anchor1 << "\t" << sims[0].anchor2 << "\t"
                                  << log2(sims[0].probRatio) + shift << "\t" << sims[1].anchor1 << "\t"
                                  << sims[1].anchor2 << "\t" << log2(sims[1].probRatio) + shift << "\t"
                                  << sims[2].anchor1 << "\t" << sims[2].anchor2 << "\t"
                                  << log2(sims[2].probRatio) + shift << std::endl;
                    }
                }

                if (++completedBatches == numBatches) {
                    std::lock_guard<std::mutex> lock(mtx);
                    cv.notify_one();
                }
            });
        }

        if (numBatches > 0) {
            std::unique_lock<std::mutex> lock(mtx);
            cv.wait(lock, [&]{ return completedBatches == numBatches; });
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

        profile.gumbelKendAnchored = MMendK;
        profile.gumbelKbegAnchored = MMbegK;
        profile.gumbelKmidAnchored = MMmidK;
        profile.lambda = MMmidL;

        entry = {
            MMendL, MMbegL, MMmidL,
            MMendK, MMbegK, MMmidK,
            MMendKsimple, MMbegKsimple, MMmidKsimple,
            MLendL, MLbegL, MLmidL,
            MLendK, MLbegK, MLmidK,
            MLendKsimple, MLbegKsimple, MLmidKsimple,
            LMendL, LMbegL, LMmidL,
            LMendK, LMbegK, LMmidK
        };
        cache.store(profile, letterFreqs, sequenceLength, border, numOfSequences, entry);
    }

    double s = scale;

    if (printVerbosity > 1) {
        std::cout << "#\tend-\tstart-\tmid-anchored\n";

        std::cout << "#lamMM\t" << entry.MMendL << "\t" << entry.MMbegL << "\t" << entry.MMmidL << "\n"

                  << "#kMM\t" << entry.MMendK / pow(s, entry.MMendL) << "\t" << entry.MMbegK / pow(s, entry.MMbegL) << "\t"
                  << entry.MMmidK / pow(s, entry.MMmidL) << "\n"

                  << "#kMM1\t" << entry.MMendKsimple / scale << "\t" << entry.MMbegKsimple / scale << "\t"
                  << entry.MMmidKsimple / scale << "\n";

        std::cout << "#lamML\t" << entry.MLendL << "\t" << entry.MLbegL << "\t" << entry.MLmidL << "\n"

                  << "#kML\t" << entry.MLendK / pow(s, entry.MLendL) << "\t" << entry.MLbegK / pow(s, entry.MLbegL) << "\t"
                  << entry.MLmidK / pow(s, entry.MLmidL) << "\n"

                  << "#kML1\t" << entry.MLendKsimple / scale << "\t" << entry.MLbegKsimple / scale << "\t"
                  << entry.MLmidKsimple / scale << "\n";

        std::cout << "#lamLM\t" << entry.LMendL << "\t" << entry.LMbegL << "\t" << entry.LMmidL << "\n"

                  << "#kLM\t" << entry.LMendK / pow(s, entry.LMendL) << "\t" << entry.LMbegK / pow(s, entry.LMbegL) << "\t"
                  << entry.LMmidK / pow(s, entry.LMmidL) << "\n";
    } else if (printVerbosity > 0) {
        std::cout << "# K: " << entry.MMendKsimple / scale << " " << entry.MMbegKsimple / scale << " "
                  << entry.MMmidKsimple / scale << "\n";
    } else {
        std::cout << "# K: " << entry.MMmidKsimple / scale << "\n";
    }

    std::cout << "# Lambda: " << entry.MMmidL << "\n";
    estimateK_n++, estimateK_lambdas += entry.MMmidL;
    std::cout << "# Avg Lambda: " << (estimateK_lambdas / estimateK_n) << "\n";
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
        if (tantanProbs[i] >= TANTAN_MASK_THRESHOLD)
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

        double alphaFS1 = INSERT1;
        double alphaFS2 = INSERT2;
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
        double deltaFS1 = DELETE1; // simulate delete
        double deltaFS2 = DELETE2;
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
            if (tantanProbs[i] >= TANTAN_MASK_THRESHOLD) {
                probs[k] = end[k];
            } else {
                double p = probs[k];
                int codonCount = (alphabetSize > 4) ? aa2codons.at(alphabet[k - 4]).size() : 1;
                probs[k] = ((1 - STOP_CODON_PROB) /* minus stop codon */ * p / codonCount);
            }
        }
        if (alphabetSize == 20) {
            probs[4 + 20] = probs[4 + 1]; // selenocysteine = cysteine
            probs[4 + 21] = probs[4 + 8]; // pyrrolysine = lysine
        }
        probs[4 + alphabetSize + 2] = 1.0 / 64.0; // for masked sequence letters
        probs[4 + alphabetSize + 3] = STOP_CODON_PROB / 3.0;
        probs[4 + alphabetSize + 4] = 0.0; // zero_idx padding
        if (tantanProbs[i] >= TANTAN_MASK_THRESHOLD)
            consensusSequence[i] |= 32;
    }

    // extra padding for wavefront DP
    for (int pad = 0; pad < 40; ++pad) {
        p.values_v2.push_back({0});
    }

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
                profile.name = word;
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
                values.insert(values.end(), 40 * profile.width, 0.0); // 40 padded rows
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
        v += p.width * (p.length + 1 + 40);
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
    int numThreadsOpt = std::thread::hardware_concurrency();
    if (numThreadsOpt == 0) {
        numThreadsOpt = 1;
    }

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
  -T N, --threads N number of threads to use (default: CPU cores)\n\
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

    const char sOpts[] = "hVve:s:m:d:D:t:l:b:T:";

    static struct option lOpts[] = {{"help", no_argument, 0, 'h'},
                                    {"version", no_argument, 0, 'V'},
                                    {"verbose", no_argument, 0, 'v'},
                                    {"threads", required_argument, 0, 'T'},
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
        case 'T':
            numThreadsOpt = intFromText(optarg);
            if (numThreadsOpt < 1)
                return badOpt();
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

    ThreadPool threadPool(numThreadsOpt);
    std::vector<DPScratch> threadScratches(numThreadsOpt);

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
        estimateK(p, bgProbs, &charVec[seqIdx], randomSeqLen, border, randomSeqNum, printVerbosity, threadPool, threadScratches);
#else
        NucDist dist = *reinterpret_cast<NucDist *>(p.debug);
        Float bgProbsDNA[256] = {0};
        bgProbsDNA[charToNumber['A']] = dist.overall['A'];
        bgProbsDNA[charToNumber['C']] = dist.overall['C'];
        bgProbsDNA[charToNumber['G']] = dist.overall['G'];
        bgProbsDNA[charToNumber['T']] = dist.overall['T'];
        estimateK(p, bgProbsDNA, &charVec[seqIdx], randomSeqLen, border, randomSeqNum,
                  printVerbosity, threadPool, threadScratches);
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
                                            1.0 / 1.0 /* p.lambda */))
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

    findFinalSimilaritiesBatched(similarities, allRequests, profiles, charVec.data(), threadPool, threadScratches);

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
