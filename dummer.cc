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
#include <utility>

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
#define OPT_t 20
#define OPT_l 5000
#define OPT_b 100
#define OPT_x 1e-5 // 0 to enable full DP mode
const int XDROP_MIN_I = 32;
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
using simd_t = Kokkos::Experimental::simd<Float>;
constexpr auto simdWidth = simd_t::size();

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
    Float delta_prime[3];
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
    bool has_pipeline_fields = false;
    std::vector<std::array<int, 4>> seeds;
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
    std::vector<std::array<int, 4>> seeds;

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
        sequence.seeds.clear();
        std::string line, word;
        getline(in, line);
        std::istringstream iss(line);
        if (!(iss >> word))
            return fail(in, "bad sequence data: no name");
#ifdef PIPELINE_MODE
        size_t slash = word.find('/');
        if (slash != std::string::npos) {
            std::string chr = word.substr(0, slash);
            std::string range = word.substr(slash + 1);
            size_t dash = range.find('-');
            if (dash != std::string::npos) {
                sequence.w_start = std::stoi(range.substr(0, dash));
                sequence.w_end = std::stoi(range.substr(dash + 1));
            }
            std::string length, profile, strand;
            if (iss >> length >> profile >> strand) {
                sequence.true_length = stoll(length.substr(length.find('=') + 1));
                sequence.target_profile = profile.substr(profile.find('=') + 1);
                sequence.is_plus = strand.find("plus_strand") != std::string::npos;
            }
            std::string seed_token;
            if (iss >> seed_token) {
                if (seed_token.rfind("seed=", 0) == 0) {
                    auto s = seed_token.substr(5);
                    size_t pos = 0;
                    while (pos < s.size()) {
                        size_t semi = s.find(';', pos);
                        std::string seed_str = (semi == std::string::npos) ? s.substr(pos) : s.substr(pos, semi - pos);
                        std::array<int, 4> seed = {};
                        size_t c = 0, st = 0;
                        for (size_t i = 0; i <= seed_str.size(); ++i) {
                            if (i == seed_str.size() || seed_str[i] == ',') {
                                if (c < 4) seed[c] = std::stoi(seed_str.substr(st, i - st));
                                ++c;
                                st = i + 1;
                            }
                        }
                        if (c == 4) sequence.seeds.push_back(seed);
                        if (semi == std::string::npos) break;
                        pos = semi + 1;
                    }
                }
            }
            sequence.has_pipeline_fields = true;
            word = chr;
        }
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
    char strand = "+-"[!s.is_plus];
    const char *seq = sim.alignedSequences.data();
    int length = sim.alignedSequences.size() / 2;
    int span1 = length - std::count(seq, seq + length, '-');
    int span2 = length - std::count(seq + length, seq + length * 2, '-');
    int start2 = strandPosition(sim.strandNum, s.length, sim.start2);
    if (s.has_pipeline_fields) {
        if (s.is_plus) {
            start2 = s.w_start - 1 + start2;
        } else {
            start2 = s.true_length - s.w_end + start2;
        }
    }
    int reportSeqLength = s.has_pipeline_fields ? s.true_length : s.length;
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

using DP_Cell = Float;
template <typename T, bool Rolling = false, int PrePad = 0, int PostPad = 0>
class FlatMatrix {
    std::vector<T> data;
    size_t logical_cols;
    size_t logical_rows;
    int offset_ = 0;

    size_t phys_cols() const { return PrePad + logical_cols + PostPad; }

    int phys_index(int j) const { return (j - offset_) + PrePad; }

public:
    FlatMatrix() : logical_cols(0), logical_rows(0) {}

    int offset() const { return offset_; }
    void set_offset(int off) { offset_ = off; }

    T get(size_t i, int j) const {
        int p = phys_index(j);
        if (p < 0 || p >= (int)phys_cols()) return T{};
        if constexpr (Rolling)
            return data[(i & 1) * phys_cols() + p];
        else
            return data[i * phys_cols() + p];
    }

    void set(size_t i, int j, T val) {
        int p = phys_index(j);
        if (p < 0 || p >= (int)phys_cols()) return;
        if constexpr (Rolling)
            data[(i & 1) * phys_cols() + p] = val;
        else
            data[i * phys_cols() + p] = val;
    }

    void resize(size_t r, size_t c) {
        offset_ = 0;
        logical_rows = r;
        logical_cols = c;
        size_t n = Rolling ? 2 : r;
        size_t total = n * phys_cols();
        data.resize(total);
    }

    void assign(size_t r, size_t c, T init = T()) {
        offset_ = 0;
        logical_rows = r;
        logical_cols = c;
        size_t n = Rolling ? 2 : r;
        size_t total = n * phys_cols();
        data.assign(total, init);
        for (size_t i = 0; i < n; ++i) {
            T *row = data.data() + i * phys_cols();
            std::fill_n(row, PrePad, T{});
            std::fill_n(row + PrePad + logical_cols, PostPad, T{});
        }
    }

    inline T &operator()(size_t i, int j) {
        if constexpr (Rolling) {
            return data[(i & 1) * phys_cols() + phys_index(j)];
        } else {
            return data[i * phys_cols() + phys_index(j)];
        }
    }

    inline const T &operator()(size_t i, int j) const {
        if constexpr (Rolling) {
            return data[(i & 1) * phys_cols() + phys_index(j)];
        } else {
            return data[i * phys_cols() + phys_index(j)];
        }
    }

    inline void clear_row(size_t i, T init_val = T()) {
        size_t actual_row = Rolling ? (i & 1) : i;
        auto row_start = data.begin() + actual_row * phys_cols() + PrePad;
        std::fill(row_start, row_start + logical_cols, init_val);
    }

    inline T *row_ptr(size_t i) {
        if constexpr (Rolling) {
            return data.data() + (i & 1) * phys_cols() + PrePad;
        } else {
            return data.data() + i * phys_cols() + PrePad;
        }
    }

    inline const T *row_ptr(size_t i) const {
        if constexpr (Rolling) {
            return data.data() + (i & 1) * phys_cols() + PrePad;
        } else {
            return data.data() + i * phys_cols() + PrePad;
        }
    }

    size_t cols() const { return logical_cols; }

    void set_cols(size_t c) {
        offset_ = 0;
        logical_cols = c;
        size_t n = Rolling ? 2 : logical_rows;
        size_t total = n * phys_cols();
        data.resize(total);
    }

    size_t rows() const { return Rolling ? 2 : logical_rows; }
};

template <typename T, int PrePad = 0, int PostPad = 0>
class PaddedVec {
    std::vector<T> data_;
    int logical_size_ = 0;
public:
    PaddedVec() = default;

    T &operator[](int i) { return data_[PrePad + i]; }

    const T &operator[](int i) const { return data_[PrePad + i]; }

    T *data() { return data_.data() + PrePad; }

    const T *data() const { return data_.data() + PrePad; }

    int size() const { return logical_size_; }
    bool empty() const { return logical_size_ == 0; }

    auto begin() { return data(); }
    auto end() { return data() + logical_size_; }
    auto begin() const { return data(); }
    auto end() const { return data() + logical_size_; }

    void resize(int n, T val = T{}) {
        logical_size_ = n;
        data_.assign(n + PrePad + PostPad, val);
        std::fill_n(data_.data(), PrePad, T{});
        std::fill_n(data_.data() + PrePad + n, PostPad, T{});
    }

    void assign(int n, T val = T{}) { resize(n, val); }

    void fill_logical(T val) {
        std::fill_n(data_.data() + PrePad, logical_size_, val);
    }

    void set_logical_size(int n) {
        logical_size_ = n;
        data_.resize(n + PrePad + PostPad, T{});
        std::fill_n(data_.data(), PrePad, T{});
        std::fill_n(data_.data() + PrePad + n, PostPad, T{});
    }

    int logical_size() const { return logical_size_; }

    void shiftExpand(const std::array<int64_t, simdWidth>& delta,
                     int new_logical, int activeCount,
                     std::vector<Float>& buf) {
        int old_phys = logical_size_ + PrePad + PostPad;
        int new_phys = new_logical + PrePad + PostPad;
        int max_cols = std::max(old_phys, new_phys);
        int needed = max_cols * simdWidth * 2;
        if ((int)buf.size() < needed) buf.resize(needed);

        // Extract old data into first half of buf
        for (int j = 0; j < max_cols; j++) {
            if (j < old_phys) {
                simd_unchecked_store(data_[j], &buf[j * simdWidth],
                                     Kokkos::Experimental::simd_flag_default);
            } else {
                std::fill_n(&buf[j * simdWidth], simdWidth, Float(0));
            }
        }

        // Shift into second half: shifted[j][k] = buf[(j - delta[k])][k]
        // Data at old position p → new position p + delta[k]
        Float* shifted = buf.data() + max_cols * simdWidth;
        for (int j = 0; j < max_cols; j++) {
            Float* dst = &shifted[j * simdWidth];
            std::fill_n(dst, simdWidth, Float(0));
            for (int k = 0; k < activeCount; k++) {
                int src = j - delta[k];
                if (src >= 0)
                    dst[k] = buf[src * simdWidth + k];
            }
        }

        logical_size_ = new_logical;
        int total = new_logical + PrePad + PostPad;
        data_.resize(total);
        for (int j = 0; j < total; j++)
            data_[j] = Kokkos::Experimental::simd_unchecked_load<simd_t>(&shifted[j * simdWidth]);
    }
};

struct DPScratch {
    // Per-sequence scalar FlatMatrices (each sequence has its own width)
    std::array<FlatMatrix<Float>, simdWidth> W1;
    std::array<FlatMatrix<Float, true>, simdWidth> W1_rolling;
    std::array<FlatMatrix<Float>, simdWidth> X;
    std::array<FlatMatrix<Float, false, 3, 0>, simdWidth> X_pfx;

    // Reusable temporary buffers for findSimilarities
    PaddedVec<simd_t, 0, 0> W0_curr, W0_next;
    PaddedVec<simd_t, 0, 0> null_probs_prefix, null_probs_suffix;
    PaddedVec<simd_t, 0, 0> Y0_next, Y0_curr;
    PaddedVec<simd_t, 0, 0> null_model_prefix, null_model_suffix;
    PaddedVec<simd_t, 3, 0> right_side;
    PaddedVec<simd_t, 0, 3> left_side;
    std::array<std::vector<AlignedSimilarity>, simdWidth> opt_profile_position;
    std::array<std::vector<bool>, simdWidth> aligned;
    std::vector<uint8_t> transposed_decoded;
    PaddedVec<simd_t, 4, 4> bg_codon_probs;

    // Per-lane anchor tracking — indexed by logical sequence position (j - seq_offset)
    std::array<std::vector<Float>, simdWidth> best_wMid;
    std::array<std::vector<Float>, simdWidth> best_wEnd;
    std::array<std::vector<int>, simdWidth> best_i;

    // Per-row per-lane band bounds for X-drop tracking (physical j coords)
    std::array<std::vector<int>, simdWidth> fwd_band_lo, fwd_band_hi;
    std::array<std::vector<int>, simdWidth> rev_band_lo, rev_band_hi;

    // Per-lane cumulative sequence offset for band harmonization
    std::array<int64_t, simdWidth> seq_offset = {};
    // Reusable temp buffer for column shifting
    std::vector<Float> shift_buf;
    // Current DP column count (grows on reoffset)
    int active_dp_width = 0;
};
struct DP_Cell_v2 {
    Float metric;
    int i, j;
    bool emit = false;

    constexpr bool operator< (const DP_Cell_v2 &other) const {
        return metric < other.metric;
    }
};

void addForwardAlignment(int idx, size_t profileLength, std::vector<SegmentPair> &alignment, int iBeg, int jBeg,
                         DPScratch &scratch) {
    int off = scratch.seq_offset[idx];
    int cols = (int)scratch.W1[idx].cols();
    int dp_bound = off + cols - 4;

    int i = iBeg, abs_pos = jBeg;
    while (i <= profileLength && abs_pos < dp_bound) {
        auto choice = std::max({
            DP_Cell_v2{.metric=(abs_pos + 3 < off + cols ? scratch.X[idx].get(i, abs_pos + 3) + scratch.W1[idx].get(i + 1, abs_pos + 3) : -INFINITY), .i=i + 1, .j=abs_pos + 3, .emit=true},
            DP_Cell_v2{.metric=scratch.W1[idx].get(i + 1, abs_pos), .i=i + 1, .j=abs_pos, .emit=false},
            DP_Cell_v2{.metric=scratch.W1[idx].get(i, abs_pos + 1), .i=i, .j=abs_pos + 1, .emit=false},
            DP_Cell_v2{.metric=scratch.left_side[abs_pos][idx], .i=INT_MAX, .j=INT_MAX, .emit=false},
        });

        if (choice.emit) {
            addForwardMatch(alignment, i, abs_pos + 1 - off);
        }
        i = choice.i, abs_pos = choice.j;
    }
}

void addReverseAlignment(int idx, std::vector<SegmentPair> &alignment, int iEnd, int jEnd,
                         DPScratch &scratch) {
    int off = scratch.seq_offset[idx];
    int i = iEnd - 1, abs_pos = jEnd;
    while (i >= 0 && abs_pos >= off) {
        DP_Cell opt_succ = 0;
        if (i >= 1 && abs_pos >= off + 3) {
            opt_succ = scratch.X_pfx[idx].get(i - 1, abs_pos - 3);
        }
        auto choice = std::max({
            DP_Cell_v2{.metric=(abs_pos >= off + 3 ? scratch.X[idx].get(i, abs_pos) + opt_succ : -INFINITY), .i=i - 1, .j=abs_pos - 3, .emit=true},
            DP_Cell_v2{.metric=(i >= 1 ? scratch.X_pfx[idx].get(i - 1, abs_pos) : -INFINITY), .i=i - 1, .j=abs_pos, .emit=false},
            DP_Cell_v2{.metric=(abs_pos > off ? scratch.X_pfx[idx].get(i, abs_pos - 1) : -INFINITY), .i=i, .j=abs_pos - 1, .emit=false},
            DP_Cell_v2{.metric=scratch.right_side[abs_pos][idx], .i=-1, .j=-1, .emit=false},
        });

        if (choice.emit) {
            addReverseMatch(alignment, i, abs_pos - 2 - off);
        }
        i = choice.i, abs_pos = choice.j;
    }
}

void addMidAnchored(int idx, size_t profileLength, std::vector<AlignedSimilarity> &similarities, int anchor1, int anchor2,
                    Float wBegAnchored, Float wEndAnchored, DPScratch &scratch) {
    Float wMidAnchored = wEndAnchored * wBegAnchored;
    AlignedSimilarity s = {wMidAnchored / scale, anchor1, anchor2, wEndAnchored};
#ifdef ALIGN
    addForwardAlignment(idx, profileLength, s.alignment, anchor1, anchor2, scratch);
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

void setCharToNumber(char *charToNumber, const char *alphabet) {
    for (int i = 0; alphabet[i]; ++i) {
        int c = alphabet[i];
        charToNumber[toupper(c)] = charToNumber[tolower(c)] = i;
    }
}

std::unordered_map<char, std::vector<std::string>> aa2codons;
std::unordered_map<char, std::vector<std::string>> &build_standard_genetic_code() {
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
void reoffset(DPScratch& scratch, int activeCount,
              const std::array<std::vector<uint8_t>*, simdWidth>& decoded,
              const Profile& profile, const Float* bg_probs_ptr,
              Float* first_arr, Float* last_arr,
              simd_t row_best_j, simd_t row_best_metric, bool capture_bands) {
    std::array<int64_t, simdWidth> deltas = {};

    Float best_reference = 0;
    int64_t reference = -1;
    for (int k = 0; k < activeCount; k++) {
        Float cur_val = row_best_metric[k];
        if (first_arr[k] <= last_arr[k] && cur_val > best_reference) {
            best_reference = cur_val;
            reference = row_best_j[k];
        }
    }

    if (reference == -1) {
        return;
    }

    for (int k = 0; k < activeCount; k++) {
        if (first_arr[k] <= last_arr[k]) {
            deltas[k] = reference - (int64_t)row_best_j[k];
        }
    }

    int64_t extra_offset_to_0_idx_left = INT64_MIN;
    for (int k = 0; k < activeCount; k++) {
        extra_offset_to_0_idx_left = std::max(extra_offset_to_0_idx_left, -((int64_t)first_arr[k] + deltas[k]));
    }

    for (int k = 0; k < activeCount; k++) {
        deltas[k] += extra_offset_to_0_idx_left;
    }

    // perform actual offsetting using deltas
    int64_t fixup_prevent_negative_offset = 0; // negative offset breaks left, right
    for (int k = 0; k < activeCount; k++) {
        scratch.seq_offset[k] += deltas[k];
        fixup_prevent_negative_offset = std::max(fixup_prevent_negative_offset, -scratch.seq_offset[k]);
    }

    scratch.active_dp_width = 0;
    for (int k = 0; k < activeCount; k++) {
        scratch.seq_offset[k] += fixup_prevent_negative_offset;
        deltas[k] += fixup_prevent_negative_offset;

        assert(scratch.seq_offset[k] >= 0);
        scratch.active_dp_width = std::max((int64_t)scratch.active_dp_width, (int64_t)decoded[k]->size() + scratch.seq_offset[k]);
    }

    int alphabetSize = profile.width - nonLetterWidth;
    int zero_idx = alphabetSize + 4;
    scratch.transposed_decoded.assign((scratch.active_dp_width + 8) * simdWidth, zero_idx);
    uint8_t* transposed_base = scratch.transposed_decoded.data() + 4 * simdWidth;
    for (int k = 0; k < activeCount; k++) {
        int offset = scratch.seq_offset[k];
        const auto& src = *decoded[k];
        for (int j = 0; j < scratch.active_dp_width; j++) {
            int src_pos = j - offset;
            transposed_base[j * simdWidth + k] =
                (0 <= src_pos && src_pos < (int)src.size()) ? src[src_pos] : zero_idx;
        }
    }

    scratch.bg_codon_probs.set_logical_size(scratch.active_dp_width);
    for (int j = -4; j < scratch.active_dp_width + 4; j++) {
        const char* indices = (const char*)&transposed_base[j * simdWidth];
        scratch.bg_codon_probs[j] = simd_t(simdLookup(bg_probs_ptr, indices));
    }
    scratch.null_probs_prefix.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    scratch.null_probs_suffix.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    scratch.W0_curr.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    scratch.W0_next.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    scratch.Y0_next.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    scratch.Y0_curr.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    scratch.null_model_prefix.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    scratch.null_model_suffix.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    if (!capture_bands) {
    scratch.right_side.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    scratch.left_side.shiftExpand(deltas, scratch.active_dp_width, activeCount, scratch.shift_buf);
    }

    // Shift per-lane band arrays into the new physical frame so that
    // fwd/rev bands (used by the out_lo/out_hi merge) stay consistent.
    {
        int n_rows = (int)profile.length + 2;
        for (int k = 0; k < activeCount; k++) {
            int64_t d = deltas[k];
            int lane_lo = (int)scratch.seq_offset[k];
            int lane_hi = (int)(scratch.seq_offset[k] + (int64_t)decoded[k]->size() - 1);
            auto shift_pair = [d, lane_lo, lane_hi](int lo, int hi) {
                if (lo == INT_MAX || hi == INT_MIN) return std::make_pair(INT_MAX, INT_MIN);
                int64_t nlo = (int64_t)lo + d;
                int64_t nhi = (int64_t)hi + d;
                if (nlo > lane_hi || nhi < lane_lo) return std::make_pair(INT_MAX, INT_MIN);
                nlo = std::max<int64_t>(nlo, lane_lo);
                nhi = std::min<int64_t>(nhi, lane_hi);
                return std::make_pair((int)nlo, (int)nhi);
            };
            for (int r = 0; r < n_rows; r++) {
                auto f = shift_pair(scratch.fwd_band_lo[k][r], scratch.fwd_band_hi[k][r]);
                scratch.fwd_band_lo[k][r] = f.first;
                scratch.fwd_band_hi[k][r] = f.second;
                auto rv = shift_pair(scratch.rev_band_lo[k][r], scratch.rev_band_hi[k][r]);
                scratch.rev_band_lo[k][r] = rv.first;
                scratch.rev_band_hi[k][r] = rv.second;
            }
        }
    }

    for (int k = 0; k < activeCount; k++) {
        scratch.W1_rolling[k].set_offset(scratch.seq_offset[k]);
        if (!capture_bands) {
            scratch.W1[k].set_offset(scratch.seq_offset[k]);
            scratch.X[k].set_offset(scratch.seq_offset[k]);
            scratch.X_pfx[k].set_offset(scratch.seq_offset[k]);
        }
    }

    for (int k = 0; k < activeCount; k++) {
        first_arr[k] += (Float)deltas[k];
        last_arr[k] += (Float)deltas[k];
    }
}

template <typename MatrixType>
simd_t simdGather(const MatrixType* matrices, Float* buf, int activeCount, int i_row, int col) {
    for (int k = 0; k < activeCount; k++)
        buf[k] = matrices[k].get(i_row, col);
    for (int k = activeCount; k < simdWidth; k++)
        buf[k] = Float(0);
    return Kokkos::Experimental::simd_unchecked_load<simd_t>(buf);
}

template <typename MatrixType>
void simdScatter(simd_t val, MatrixType* matrices, Float* buf, int activeCount, int i_row, int col) {
    simd_unchecked_store(val, buf, Kokkos::Experimental::simd_flag_default);
    for (int k = 0; k < activeCount; k++)
        matrices[k].set(i_row, col, buf[k]);
}

uint8_t* buildTransposedDecoded(DPScratch& scratch, int active_dp_width, int activeCount,
                                 int zero_idx, const std::array<std::vector<uint8_t>*, simdWidth>& decoded) {
    scratch.transposed_decoded.assign((active_dp_width + 8) * simdWidth, zero_idx);
    uint8_t* transposed_base = scratch.transposed_decoded.data() + 4 * simdWidth;
    for (int idx = 0; idx < activeCount; idx++) {
        for (size_t j = 0; j < decoded[idx]->size(); j++) {
            transposed_base[j * simdWidth + idx] = (*decoded[idx])[j];
        }
    }
    return transposed_base;
}

struct DPArgs {
    const std::array<std::vector<std::vector<int>>, simdWidth>* seed_j = nullptr;
    const std::array<std::vector<int>, simdWidth>* in_lo = nullptr;
    const std::array<std::vector<int>, simdWidth>* in_hi = nullptr;
    std::array<std::vector<int>, simdWidth>* out_lo = nullptr;
    std::array<std::vector<int>, simdWidth>* out_hi = nullptr;
    int seed_margin = 0;
};

#ifdef DEBUG_PRINT_ROWS
static void logBandStats(const char* label, int i, int activeCount,
                         Float* first_arr, Float* last_arr, simd_t best_metric,
                         const std::array<std::vector<uint8_t>*, simdWidth>& decoded) {
    std::ostringstream o;
    o << label << " row " << i << " band:";
    for (int k = 0; k < activeCount; k++) {
        int len = (int)decoded[k]->size();
        if (first_arr[k] <= last_arr[k]) {
            double pct = (last_arr[k] - first_arr[k] + 1.0) / len * 100.0;
            o << " " << k << ":" << (int)pct << "% (" << best_metric[k] << ") ";
        }
    }
    o << std::endl;
    std::cerr << o.str();
}
#endif

struct SeedGateResult {
    bool has_seeds;
    std::vector<simd_t> gate;
};

static SeedGateResult buildSeedGatedNullGate(
    int seed_j_shift, int row, int activeCount, int dp_width, int seed_margin,
    const std::array<std::vector<std::vector<int>>, simdWidth>& seed_j,
    DPScratch& scratch,
    const std::array<std::vector<uint8_t>*, simdWidth>& decoded) {
    SeedGateResult result{false, {}};
    for (int k = 0; k < activeCount && !result.has_seeds; k++)
        result.has_seeds = !seed_j[k][row].empty();
    if (result.has_seeds) {
        result.gate.assign(dp_width + 4, simd_t(0));
        for (int k = 0; k < activeCount; k++) {
            const auto& seeds = seed_j[k][row];
            if (seeds.empty()) continue;
            alignas(64) Float one_f[simdWidth] = {};
            one_f[k] = 1.0f;
            simd_t lane_one = Kokkos::Experimental::simd_unchecked_load<simd_t>(one_f);
            int lane_lo = (int)scratch.seq_offset[k];
            int lane_hi = (int)(scratch.seq_offset[k] + (int64_t)decoded[k]->size() - 1);
            for (int sj : seeds) {
                // seed_j stores logical (sequence) coords; convert to physical j.
                int sj_phys = sj + (int)scratch.seq_offset[k] + seed_j_shift;
                int lo = std::max(lane_lo, sj_phys - seed_margin);
                int hi = std::min(lane_hi, sj_phys + seed_margin);
                for (int gj = lo; gj <= hi; gj++)
                    result.gate[gj] = Kokkos::max(result.gate[gj], lane_one);
            }
        }
    }
    return result;
}

void findSimilarities(std::array<std::vector<AlignedSimilarity>, simdWidth> &similarities, const Profile &profile,
             const std::array<std::vector<uint8_t>*, simdWidth> &decoded, std::array<Float, simdWidth> minProbRatio,
             DPScratch &scratch, int activeCount, bool useXDrop, const DPArgs& args = DPArgs{}) {
    assert(0 < activeCount && activeCount <= simdWidth);

    const bool capture_bands = (args.out_lo != nullptr);
    const bool seed_gating = (args.seed_j != nullptr);

    const simd_t simd_OPT_x(Float(OPT_x));

    int maxSequenceLength = 0;
    for (int idx = 0; idx < activeCount; idx++) {
        maxSequenceLength = std::max(maxSequenceLength, (int)decoded[idx]->size());
        scratch.seq_offset[idx] = 0;
    }

    scratch.active_dp_width = maxSequenceLength;

    int alphabetSize = profile.width - nonLetterWidth;
    int zero_idx = alphabetSize + 4; // new padded index that maps to 0.0

    uint8_t* transposed_base = buildTransposedDecoded(scratch, scratch.active_dp_width, activeCount, zero_idx, decoded);

    alignas(64) Float actual_sequence_length[simdWidth] = {};
    for (int idx = 0; idx < activeCount; idx++) {
        actual_sequence_length[idx] = (Float)decoded[idx]->size();
    }
    simd_t actual_seq_len = Kokkos::Experimental::simd_unchecked_load<simd_t>(actual_sequence_length);

    auto &null_probs_prefix = scratch.null_probs_prefix;
    auto &null_probs_suffix = scratch.null_probs_suffix;
    null_probs_prefix.resize(scratch.active_dp_width + 4); null_probs_suffix.resize(scratch.active_dp_width + 4);

    null_probs_suffix[scratch.active_dp_width] = 0;
    for (int i = scratch.active_dp_width - 1; i >= 0; i--) {
        const Float *bg_probs_ptr = profile.log2_bg_probs.data() + 4;
        const char* indices = (const char*)&transposed_base[i * simdWidth];
        SimdFloat bg_raw = simdLookup(bg_probs_ptr, indices);
        simd_t bg_codon_emit_probs(bg_raw);

        Kokkos::Experimental::simd_mask<Float> msk = i + 2 < actual_seq_len;
        Kokkos::Experimental::simd_mask<Float> msk_fs2 = i + 1 < actual_seq_len;
        Kokkos::Experimental::simd_mask<Float> msk_fs1 = i < actual_seq_len;
        simd_t full_codon = (Float)log2(1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2) + bg_codon_emit_probs + null_probs_suffix[i + 3];
        simd_t partial_codon = (Float)log2(1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2) + (Float)log2(0.25) * (actual_seq_len - (Float)i);
        simd_t t1 = Kokkos::Experimental::condition(msk, full_codon, partial_codon);

        simd_t fs1 = (Float)log2(BACKGROUND_FRAMESHIFT_RATE * 0.25) + null_probs_suffix[i + 1];
        simd_t fs2 = (Float)log2(BACKGROUND_FRAMESHIFT_RATE_2 * 0.0625) + null_probs_suffix[i + 2];
        simd_t neg_inf_vec((Float)-INFINITY);
        simd_t t_fs1 = Kokkos::Experimental::condition(msk_fs1, fs1, neg_inf_vec);
        simd_t t_fs2 = Kokkos::Experimental::condition(msk_fs2, fs2, neg_inf_vec);

        null_probs_suffix[i] = log2_sum_exp(t1, log2_sum_exp(t_fs1, t_fs2));
    }
    null_probs_prefix[scratch.active_dp_width] = 0;
    const Float *log2_bg_probs_ptr = profile.log2_bg_probs.data() + 4;
    const Float log2_1_bg_fs = log2(1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2);
    const Float log2_bg_fs_025 = log2(BACKGROUND_FRAMESHIFT_RATE * 0.25);
    const Float log2_bg_fs2_00625 = log2(BACKGROUND_FRAMESHIFT_RATE_2 * 0.0625);
    const Float log2_025 = log2(0.25);
    const simd_t simd_log2_1_bg_fs(log2_1_bg_fs);
    const simd_t simd_log2_bg_fs_025(log2_bg_fs_025);
    const simd_t simd_log2_bg_fs2_00625(log2_bg_fs2_00625);
    const simd_t simd_log2_025(log2_025);
    const simd_t simd_neg_inf(-INFINITY);

    for (int i = 0; i < scratch.active_dp_width; i++) {
        // Vectorized bg emission lookup via transposed_base
        const char* indices = (const char*)&transposed_base[(i - 2) * simdWidth];
        SimdFloat bg_raw = simdLookup(log2_bg_probs_ptr, indices);
        simd_t bg_codon_emit_probs(bg_raw);

        // Mask for variable-length sequences (per-lane)
        Kokkos::Experimental::simd_mask<Float> msk_valid = (Float)i < actual_seq_len;

        // t1: codon branch — i is scalar, so use plain if/else
        simd_t t1;
        if (i >= 3) {
            t1 = simd_log2_1_bg_fs + bg_codon_emit_probs + null_probs_prefix[i - 3];
        } else if (i == 2) {
            t1 = simd_log2_1_bg_fs + bg_codon_emit_probs;
        } else {
            t1 = simd_log2_1_bg_fs + simd_log2_025 * (Float)(i + 1);
        }

        // t2: 1-bp frameshift branch
        simd_t t2 = simd_log2_bg_fs_025 + (i > 0 ? null_probs_prefix[i - 1] : simd_t(0));

        // t3: 2-bp frameshift branch
        simd_t t3;
        if (i >= 2) {
            t3 = simd_log2_bg_fs2_00625 + null_probs_prefix[i - 2];
        } else if (i == 1) {
            t3 = simd_log2_bg_fs2_00625 + simd_t(0);
        } else {
            t3 = simd_neg_inf;
        }

        // log2_sum_exp and mask out-of-bounds lanes
        simd_t result = log2_sum_exp(t1, log2_sum_exp(t2, t3));
        null_probs_prefix[i] = Kokkos::Experimental::condition(msk_valid, result, simd_neg_inf);
    }

    alignas(64) Float null_emit_tmp[simdWidth] = {0}, null_seq_log_prob[simdWidth] = {0};

    for (int idx = 0; idx < activeCount; idx++) {
        if (actual_sequence_length[idx] >= 3) {
            null_seq_log_prob[idx] =
            log2_sum_exp(log2_sum_exp(null_probs_prefix[actual_sequence_length[idx] - 1][idx],null_probs_prefix[actual_sequence_length[idx] - 2][idx]),
                null_probs_prefix[actual_sequence_length[idx] - 3][idx]
            );
        } else {
            // For very short sequences, treat as non‑alignable (e.g., -INF or 0)
            null_seq_log_prob[idx] = -INFINITY;
        }

        Float null_emit_raw = exp2(-(null_seq_log_prob[idx] / actual_sequence_length[idx]));
        null_emit_tmp[idx] = null_emit_raw;
    }

    auto null_emit_1 = Kokkos::Experimental::simd_unchecked_load<simd_t>(null_emit_tmp);
    auto null_emit_2 = null_emit_1 * null_emit_1;
    auto null_emit_3 = null_emit_2 * null_emit_1;

    null_emit_2 *= (Float)(0.25 * 0.25);
    null_emit_1 *= (Float)0.25;

    auto null_seq_log_prob_simd = Kokkos::Experimental::simd_unchecked_load<simd_t>(null_seq_log_prob);
    simd_t inv_actual_seq_len = (Float)1.0 / actual_seq_len;
    simd_t null_prob_per_pos = null_seq_log_prob_simd * inv_actual_seq_len;

    scratch.W0_curr.resize(scratch.active_dp_width + 4, simd_t(0.0));
    scratch.W0_next.resize(scratch.active_dp_width + 4, simd_t(0.0));

    for (int k = 0; k < activeCount; k++) {
        if (capture_bands)
            scratch.W1_rolling[k].assign(profile.length + 2, scratch.active_dp_width + 4);
        else
            scratch.W1[k].assign(profile.length + 2, scratch.active_dp_width + 4);
    }

    if (args.in_lo && args.in_hi) {
        for (int k = 0; k < activeCount; k++) {
            scratch.fwd_band_lo[k] = (*args.in_lo)[k];
            scratch.fwd_band_hi[k] = (*args.in_hi)[k];
            scratch.rev_band_lo[k] = (*args.in_lo)[k];
            scratch.rev_band_hi[k] = (*args.in_hi)[k];
        }
    } else {
        for (int k = 0; k < activeCount; k++) {
            int lane_lo = (int)scratch.seq_offset[k];
            int lane_hi = (int)(scratch.seq_offset[k] + (int64_t)decoded[k]->size() - 1);
            scratch.fwd_band_lo[k].assign(profile.length + 2, seed_gating ? INT_MAX : lane_lo);
            scratch.fwd_band_hi[k].assign(profile.length + 2, seed_gating ? INT_MIN : lane_hi);
            scratch.rev_band_lo[k].assign(profile.length + 2, seed_gating ? INT_MAX : lane_lo);
            scratch.rev_band_hi[k].assign(profile.length + 2, seed_gating ? INT_MIN : lane_hi);
        }
    }

    const size_t padded_seq_len = scratch.active_dp_width + 4;
    auto &Y0_next = scratch.Y0_next; Y0_next.resize(padded_seq_len, 0.0);
    auto &Y0_curr = scratch.Y0_curr; Y0_curr.resize(padded_seq_len, 0.0);

    auto &null_model_prefix = scratch.null_model_prefix; null_model_prefix.resize(padded_seq_len, 0.0);

    for (int j = 0; j < scratch.active_dp_width; j++) {
        simd_t exponent = -null_prob_per_pos * (actual_seq_len - (Float)1.0 - (Float)j) + null_probs_suffix[j + 1];
        null_model_prefix[j] = Kokkos::exp2(exponent);
    }
    auto &null_model_suffix = scratch.null_model_suffix; null_model_suffix = null_model_prefix;
    auto &left_side = scratch.left_side;
    auto &right_side = scratch.right_side;
#ifdef ALIGN
    left_side.assign(null_model_prefix.size(), 0.0);
    right_side.assign(null_model_prefix.size(), 0.0);
    if (!capture_bands) {
        for (int k = 0; k < activeCount; k++) {
            scratch.X[k].assign(profile.length + 2, scratch.active_dp_width + 4);
            scratch.X_pfx[k].assign(profile.length + 2, scratch.active_dp_width + 4);
        }
    }
#endif
        const Float *bg_probs_ptr = profile.bg_probs.data() + 4;
        scratch.bg_codon_probs.resize(scratch.active_dp_width);
        simd_t* bg_codon_probs_base = scratch.bg_codon_probs.data();
        for (int j = -4; j < scratch.active_dp_width + 4; j++) {
            const char* indices = (const char*)&transposed_base[j * simdWidth];
            bg_codon_probs_base[j] = simd_t(simdLookup(bg_probs_ptr, indices));
        }

        simd_t global_best_backward = simd_t(0.0);
        {
        alignas(64) Float gather_buf[simdWidth];
        auto gather_w1 = [&](int i_row, int col) -> simd_t {
            return capture_bands
                ? simdGather(scratch.W1_rolling.data(), gather_buf, activeCount, i_row, col)
                : simdGather(scratch.W1.data(), gather_buf, activeCount, i_row, col);
        };
        auto scatter_w1 = [&](simd_t val, int i_row, int col) {
            capture_bands
                ? simdScatter(val, scratch.W1_rolling.data(), gather_buf, activeCount, i_row, col)
                : simdScatter(val, scratch.W1.data(), gather_buf, activeCount, i_row, col);
        };

        for (int i = profile.length; i >= 0; i--) {
            const Params &params_cur = profile.values_v2[i];
            const Float *params_emission_probabilities = profile.values + (i)*profile.width + 4;

            const simd_t C_enter = params_cur.enter_match_probability * null_emit_3;
            const simd_t C_delta0 = params_cur.delta_prime[0];
            const simd_t C_delta1 = params_cur.delta_prime[1] * null_emit_2;
            const simd_t C_delta2 = params_cur.delta_prime[2] * null_emit_1;
            const simd_t C_alpha0 = params_cur.alpha_prime[0] * null_emit_3;
            const simd_t C_alpha1 = params_cur.alpha_prime[1] * null_emit_1;
            const simd_t C_alpha2 = params_cur.alpha_prime[2] * null_emit_2;
            const simd_t C_beta0 = params_cur.beta_prime[0] * null_emit_3;
            simd_t Z0_ring[4] = {0, 0, 0, 0};

            const simd_t C_eps0 = params_cur.epsilon_prime;
            const simd_t C_scale = scale;

            // Seed band expansion: revive/extend bands at seed rows
            if (seed_gating) {
                for (int k = 0; k < activeCount; k++) {
                    const auto& seeds = (*args.seed_j)[k][i];
                    if (seeds.empty()) continue;
#ifdef DEBUG_PRINT_ROWS
                    {
                        std::ostringstream o;
                        o << "BWD row=" << i << " lane=" << k << " seeds_j:";
                        for (int sj : seeds) o << " " << sj;
                        o << "\n";
                        std::cerr << o.str();
                    }
#endif

                    int lane_lo = (int)scratch.seq_offset[k];
                    int lane_hi = (int)(scratch.seq_offset[k] + (int64_t)decoded[k]->size() - 1);
                    int slo = std::clamp(*std::min_element(seeds.begin(), seeds.end()) - 1 - args.seed_margin + (int)scratch.seq_offset[k], lane_lo, lane_hi);
                    int shi = std::clamp(*std::max_element(seeds.begin(), seeds.end()) - 1 + args.seed_margin + (int)scratch.seq_offset[k], lane_lo, lane_hi);
                    scratch.rev_band_lo[k][i] = std::min(scratch.rev_band_lo[k][i], slo);
                    scratch.rev_band_hi[k][i] = std::max(scratch.rev_band_hi[k][i], shi);
                    scratch.fwd_band_lo[k][i] = std::min(scratch.fwd_band_lo[k][i], slo);
                    scratch.fwd_band_hi[k][i] = std::max(scratch.fwd_band_hi[k][i], shi);
                }
            }

            auto [row_has_seeds, null_gate] = seed_gating
                ? buildSeedGatedNullGate(-1, i, activeCount, scratch.active_dp_width, args.seed_margin, *args.seed_j, scratch, decoded)
                : SeedGateResult{false, {}};

            bool use_rev_band = seed_gating || (args.in_lo != nullptr) || useXDrop;
            int rev_lo, rev_hi;
            if (use_rev_band) {
                rev_lo = INT_MAX; rev_hi = INT_MIN;
                for (int k = 0; k < activeCount; k++) {
                    if (scratch.rev_band_lo[k][i] <= scratch.rev_band_hi[k][i]) {
                        rev_lo = std::min(rev_lo, scratch.rev_band_lo[k][i]);
                        rev_hi = std::max(rev_hi, scratch.rev_band_hi[k][i]);
                    }
                }
            } else {
                rev_lo = 0; rev_hi = scratch.active_dp_width - 1;
            }

            alignas(64) Float rev_lane_lo_f[simdWidth], rev_lane_hi_f[simdWidth];
            for (int k = 0; k < activeCount; k++) {
                rev_lane_lo_f[k] = (Float)(scratch.rev_band_lo[k][i]);
                rev_lane_hi_f[k] = (Float)(scratch.rev_band_hi[k][i]);
            }
            simd_t rev_lane_lo_vec = Kokkos::Experimental::simd_unchecked_load<simd_t>(rev_lane_lo_f);
            simd_t rev_lane_hi_vec = Kokkos::Experimental::simd_unchecked_load<simd_t>(rev_lane_hi_f);

            simd_t bwd_lane_first_active = simd_t((Float)INT_MAX);
            simd_t bwd_lane_last_active = simd_t((Float)INT_MIN);

            for (int j = rev_hi; j >= rev_lo; j--) {
                int r_0 = j & 3;
                int r_1 = (j + 1) & 3;
                int r_2 = (j + 2) & 3;
                int r_3 = (j + 3) & 3;

                const char* indices = (const char*)&transposed_base[(j + 1) * simdWidth];
                SimdFloat codon_raw = simdLookup(params_emission_probabilities, indices);
                simd_t codon_emit_probs(codon_raw);
                simd_t bg_codon_emit_probs = bg_codon_probs_base[j + 1];

                simd_t w_val =
                    gather_w1(i + 1, j + 3) * codon_emit_probs * C_enter +
                    Y0_next[j + 0] * C_delta0 +
                    gather_w1(i + 1, j + 2) * C_delta1 +
                    gather_w1(i + 1, j + 1) * C_delta2 +
                    Z0_ring[r_3] * bg_codon_emit_probs * C_alpha0 +
                    gather_w1(i, j + 1) * C_alpha1 +
                    gather_w1(i, j + 2) * C_alpha2
                 + null_model_prefix[j] * (row_has_seeds ?  null_gate[j] : simd_t(1)) * C_scale;

                // Apply per-lane band mask first so X-drop threshold is correct
                if (use_rev_band) {
                    simd_t j_vec((Float)j);
                    simd_t rev_in_band = Kokkos::Experimental::condition(
                        (j_vec >= rev_lane_lo_vec) && (j_vec <= rev_lane_hi_vec), simd_t(1), simd_t(0));
                    w_val *= rev_in_band;
                }
                // Backward X-drop: score = W1[i][j] * null_model_prefix
                simd_t bwd_score = w_val * scratch.null_model_prefix[j];
                global_best_backward = Kokkos::max(global_best_backward, bwd_score);
                simd_t bwd_active = useXDrop
                    ? Kokkos::Experimental::condition(bwd_score >= global_best_backward * simd_OPT_x && bwd_score > simd_t(0), simd_t(1), simd_t(0))
                    : simd_t(1);
                w_val *= bwd_active;
                scatter_w1(w_val, i, j);
#ifdef ALIGN
                if (!capture_bands)
                    right_side[j] += w_val;
#endif

                Y0_curr[j] = Kokkos::fma(C_eps0, Y0_next[j], w_val);
                simd_t z0_future = Z0_ring[r_3] * bg_codon_emit_probs;
                Z0_ring[r_0] = Kokkos::fma(C_beta0, z0_future, w_val);

                bwd_lane_first_active = Kokkos::Experimental::condition(
                    bwd_active > 0 && (Float)j < bwd_lane_first_active,
                    simd_t((Float)j), bwd_lane_first_active);
                bwd_lane_last_active = Kokkos::Experimental::condition(
                    bwd_active > 0 && (Float)j > bwd_lane_last_active,
                    simd_t((Float)j), bwd_lane_last_active);
            }

            if (useXDrop && i > 0) {
                alignas(64) Float first_arr[simdWidth], last_arr[simdWidth];
                simd_unchecked_store(bwd_lane_first_active, first_arr, Kokkos::Experimental::simd_flag_default);
                simd_unchecked_store(bwd_lane_last_active, last_arr, Kokkos::Experimental::simd_flag_default);

                for (int k = 0; k < activeCount; k++) {
                    int first = (int)first_arr[k], last = (int)last_arr[k];
                    if (first >= 0 && last < scratch.active_dp_width && first <= last) {
                        int lane_lo = (int)scratch.seq_offset[k];
                        int lane_hi = (int)(scratch.seq_offset[k] + (int64_t)decoded[k]->size() - 1);
                        int seq_lo = std::max(first - 3, lane_lo);
                        int seq_hi = std::min(last, lane_hi);
                        if (seq_lo <= seq_hi) {
                            scratch.rev_band_lo[k][i - 1] = std::clamp(seq_lo, lane_lo, lane_hi);
                            scratch.rev_band_hi[k][i - 1] = std::clamp(seq_hi, lane_lo, lane_hi);
                        } else {
                            scratch.rev_band_lo[k][i - 1] = INT_MAX;
                            scratch.rev_band_hi[k][i - 1] = INT_MIN;
                        }
                    } else {
                        scratch.rev_band_lo[k][i - 1] = INT_MAX;
                        scratch.rev_band_hi[k][i - 1] = INT_MIN;
                    }
                }

            if (useXDrop) {
#ifdef DEBUG_PRINT_ROWS
                logBandStats("# bwd", i, activeCount, first_arr, last_arr, global_best_backward, decoded);
#endif
            }
            }

            std::swap(Y0_curr, Y0_next);
        }
        }

    std::fill(Y0_next.begin(), Y0_next.begin() + padded_seq_len, simd_t(0.0));

    if (!capture_bands) {
    for (int k = 0; k < activeCount; k++) {
        int len = (int)decoded[k]->size();
        scratch.best_wMid[k].assign(len, Float(-INFINITY));
        scratch.best_wEnd[k].assign(len, Float(0.0));
        scratch.best_i[k].assign(len, -1);
    }
    }
    for (int j = 0; j < scratch.active_dp_width; j++) {
        simd_t exponent = -null_prob_per_pos * (Float)(j + 1) + null_probs_prefix[j];

        null_model_prefix[j] = Kokkos::exp2(exponent);
    }

#ifdef ALIGN
    if (!capture_bands) {
    // seems to be a clean way of getting expected value of null-sided junctions
    // TODO: verify logic
    for (int j = scratch.active_dp_width - 1; j >= 0; j--) {
        const char* indices = (const char*)&transposed_base[(j - 2) * simdWidth];
        SimdFloat bg_raw = simdLookup(bg_probs_ptr, indices);
        simd_t bg_codon_emit_probs(bg_raw);

        right_side[j - 3] += (Float)(1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2) * bg_codon_emit_probs *
                             null_emit_3 * right_side[j];
        right_side[j - 1] += (Float)(BACKGROUND_FRAMESHIFT_RATE * 0.25) * null_emit_1 * right_side[j];
        right_side[j - 2] += (Float)(BACKGROUND_FRAMESHIFT_RATE_2 * 0.0625) * null_emit_2 * right_side[j];

        right_side[j] *= null_model_prefix[j]; // TODO: this one specifically (might be off by 1 idk)
        right_side[j] = Kokkos::max(right_side[j], Float(0.0));
    }
    for (int j = 1; j < scratch.active_dp_width; j++) {
        right_side[j] += right_side[j - 1];
    }
    }
#endif

    size_t active_ij = 0, total_ij = 0;
    size_t active_cells = 0, total_cells = 0;
    size_t band_cells = 0;
    simd_t global_best = simd_t(0.0);
    for (int i = 0; i <= profile.length; i++) {
        const Params &params_cur = profile.values_v2[i];
        const Float *params_emission_probabilities = profile.values + (i)*profile.width + 4;

        // Pre-calculate constants
        const simd_t C_enter = params_cur.enter_match_probability * null_emit_3;
        const simd_t C_alpha0 = params_cur.alpha_prime[0];
        const simd_t C_alpha1 = params_cur.alpha_prime[1] * null_emit_1;
        const simd_t C_alpha2 = params_cur.alpha_prime[2] * null_emit_2;
        const simd_t C_beta0 = params_cur.beta_prime[0];
        const simd_t C_delta0 = params_cur.delta_prime[0];
        const simd_t C_delta1 = params_cur.delta_prime[1] * null_emit_2;
        const simd_t C_delta2 = params_cur.delta_prime[2] * null_emit_1;
        const simd_t C_eps0 = params_cur.epsilon_prime;
        const simd_t C_scale = scale;

        simd_t Z0_ring[4] = {0, 0, 0, 0};
        simd_t Z1_ring[4] = {0, 0, 0, 0};
        simd_t Z2_ring[4] = {0, 0, 0, 0};

        // Raw row pointers — avoid repeated i*cols in inner loop
        simd_t *__restrict__ w0_row_i = scratch.W0_curr.data();
        simd_t *__restrict__ w0_row_ip1 = (i + 1 <= profile.length) ? scratch.W0_next.data() : nullptr;

        alignas(64) Float gather_buf[simdWidth];
        auto gather_w1 = [&](int i_row, int col) -> simd_t {
            return capture_bands
                ? simdGather(scratch.W1_rolling.data(), gather_buf, activeCount, i_row, col)
                : simdGather(scratch.W1.data(), gather_buf, activeCount, i_row, col);
        };

        bool w1_ip1_avail = (i + 1 <= profile.length);
        auto gather_w1_ip1 = [&](int col) -> simd_t {
            if (!w1_ip1_avail || capture_bands) return simd_t(0);
            return gather_w1(i + 1, col);
        };
        auto gather_w1_i = [&](int col) -> simd_t {
            if (capture_bands) return simd_t(0);
            return gather_w1(i, col);
        };

#ifdef ALIGN
        auto gather_x = [&](int i_row, int col) -> simd_t {
            return simdGather(scratch.X.data(), gather_buf, activeCount, i_row, col);
        };
        auto gather_xpfx = [&](int i_row, int col) -> simd_t {
            return simdGather(scratch.X_pfx.data(), gather_buf, activeCount, i_row, col);
        };
        auto scatter_x = [&](simd_t val, int i_row, int col) {
            simdScatter(val, scratch.X.data(), gather_buf, activeCount, i_row, col);
        };
        auto scatter_xpfx = [&](simd_t val, int i_row, int col) {
            simdScatter(val, scratch.X_pfx.data(), gather_buf, activeCount, i_row, col);
        };
#endif

        // Seed band expansion: revive/extend bands at seed rows
        if (seed_gating) {
            for (int k = 0; k < activeCount; k++) {
                const auto& seeds = (*args.seed_j)[k][i];
                if (seeds.empty()) continue;
#ifdef DEBUG_PRINT_ROWS
                {
                    std::ostringstream o;
                    o << "FWD row=" << i << " lane=" << k << " seeds_j:";
                    for (int sj : seeds) o << " " << sj;
                    o << "\n";
                    std::cerr << o.str();
                }
#endif
                int lane_lo = (int)scratch.seq_offset[k];
                int lane_hi = (int)(scratch.seq_offset[k] + (int64_t)decoded[k]->size() - 1);
                int slo = std::clamp(*std::min_element(seeds.begin(), seeds.end()) + (int)scratch.seq_offset[k], lane_lo, lane_hi);
                int shi = std::clamp(*std::max_element(seeds.begin(), seeds.end()) + (int)scratch.seq_offset[k], lane_lo, lane_hi);
                scratch.fwd_band_lo[k][i] = std::min(scratch.fwd_band_lo[k][i], slo);
                scratch.fwd_band_hi[k][i] = std::max(scratch.fwd_band_hi[k][i], shi);
            }
        }

        auto [fwd_row_has_seeds, fwd_null_gate] = seed_gating
            ? buildSeedGatedNullGate(0, i, activeCount, scratch.active_dp_width, args.seed_margin, *args.seed_j, scratch, decoded)
            : SeedGateResult{false, {}};

        bool use_fwd_band = seed_gating || (args.in_lo != nullptr) || useXDrop;
        int lo, hi;
        if (use_fwd_band) {
            lo = INT_MAX, hi = INT_MIN;
            for (int k = 0; k < activeCount; k++) {
                if (scratch.fwd_band_lo[k][i] <= scratch.fwd_band_hi[k][i]) {
                    lo = std::min(lo, scratch.fwd_band_lo[k][i]);
                    hi = std::max(hi, scratch.fwd_band_hi[k][i]);
                }
            }
        } else {
            lo = 0; hi = scratch.active_dp_width - 1;
        }

        // Precompute per-lane band mask vectors for this row
        alignas(64) Float fwd_lane_lo_f[simdWidth], fwd_lane_hi_f[simdWidth];
        for (int k = 0; k < activeCount; k++) {
            fwd_lane_lo_f[k] = (Float)(scratch.fwd_band_lo[k][i]);
            fwd_lane_hi_f[k] = (Float)(scratch.fwd_band_hi[k][i]);
        }
        simd_t fwd_lane_lo_vec = Kokkos::Experimental::simd_unchecked_load<simd_t>(fwd_lane_lo_f);
        simd_t fwd_lane_hi_vec = Kokkos::Experimental::simd_unchecked_load<simd_t>(fwd_lane_hi_f);

        int seq_start = lo;
        band_cells += std::max(hi - seq_start + 1, 0);

        // Shift register for w[1..3] — avoids 3 matrix reads per iteration
        simd_t w_shift[3] = {0,0, 0}; // w_shift[0]=w0(i,j-1), [1]=w0(i,j-2), [2]=w0(i,j-3)
        simd_t pfx_prev = simd_t(0.0); // X_pfx(i, j-1) rolling value

        const simd_t simd_invScale(invScale);

        simd_t lane_first_active = simd_t((Float)scratch.active_dp_width);
        simd_t lane_last_active = simd_t(0.0);
        simd_t row_active_count = simd_t(0);
        simd_t row_best_metric = simd_t(-INFINITY);
        simd_t row_best_j = simd_t(0.0);

        alignas(64) Float seq_end_f[simdWidth] = {};
        for (int k = 0; k < activeCount; k++)
            seq_end_f[k] = (Float)(scratch.seq_offset[k] + (int)decoded[k]->size());
        simd_t seq_end = Kokkos::Experimental::simd_unchecked_load<simd_t>(seq_end_f);

        active_ij += hi - seq_start + 1;
        total_ij += maxSequenceLength;
        for (int j = seq_start; j <= hi; j++) {
            bool in_band = (j >= lo);

            simd_t w1 = w_shift[0], w2 = w_shift[1], w3 = w_shift[2];

            int r_0 = j & 3;
            int r_3 = (j - 3) & 3;

            const char* indices = (const char*)&transposed_base[(j - 2) * simdWidth];
            SimdFloat codon_raw = simdLookup(params_emission_probabilities, indices);

            simd_t codon_emit_probs(codon_raw);
            simd_t bg_codon_emit_probs = bg_codon_probs_base[j - 2];

            simd_t X_ij = C_enter * codon_emit_probs * w3;

#ifdef ALIGN
            if (!capture_bands) {
                simd_t X_ij_EV = X_ij * gather_w1_ip1(j) * simd_invScale;
                scatter_x(X_ij_EV, i, j);

                simd_t opt_succ = (i - 1 >= 0) ? gather_xpfx(i - 1, j - 3) : simd_t(0);

                simd_t pfx_mx = (i - 1 >= 0) ? gather_xpfx(i - 1, j) : simd_t(0);
                pfx_mx = Kokkos::max(pfx_mx, pfx_prev);

                pfx_mx = Kokkos::max(pfx_mx, X_ij_EV + opt_succ);
                pfx_mx = Kokkos::max(pfx_mx, right_side[j]);
                scatter_xpfx(pfx_mx, i, j);
                pfx_prev = pfx_mx;
            }
#endif
            Z0_ring[r_0] =
                bg_codon_emit_probs * null_emit_3 *
                (C_alpha0 * w3 + C_beta0 * Z0_ring[r_3]);
            Z1_ring[r_0] = C_alpha1 * w1;
            Z2_ring[r_0] = C_alpha2 * w2;

            simd_t w0 = w0_row_i[j];
            w0 += Z0_ring[r_0] + Z1_ring[r_0] + Z2_ring[r_0] + null_model_prefix[j] * (fwd_row_has_seeds ?  fwd_null_gate[j] : simd_t(1)) * C_scale;
#ifdef ALIGN
            if (!capture_bands)
                left_side[j] += w0;
#endif

            // Apply per-lane band mask first so X-drop threshold is correct
            simd_t j_vec((Float)j);
            if (use_fwd_band) {
                simd_t fwd_in_band = Kokkos::Experimental::condition(
                    (j_vec >= fwd_lane_lo_vec) && (j_vec <= fwd_lane_hi_vec), simd_t(1), simd_t(0));
                w0 *= fwd_in_band;
            }

            simd_t wMid = w0 * gather_w1_i(j) * simd_invScale;

            // X-drop: track running max and determine active cells
            simd_t xdrop_score = (capture_bands)
                ? w0 * scratch.null_model_suffix[j]
                : wMid;
            if (in_band) {
                global_best = Kokkos::max(global_best, xdrop_score);
            }
            simd_t is_active = (in_band && useXDrop)
                ? Kokkos::Experimental::condition(xdrop_score >= global_best * simd_OPT_x && xdrop_score > 0 /* disallow 0 probability regardless to avoid infinite extension */, simd_t(1), simd_t(0))
                : Kokkos::Experimental::condition(fwd_lane_lo_vec <= j_vec && j_vec <= fwd_lane_hi_vec, simd_t(1), simd_t(0));
            row_active_count += is_active;
            w0 *= is_active;
            w0_row_i[j] = w0;

            Z0_ring[r_0] *= is_active;
            Z1_ring[r_0] *= is_active;
            Z2_ring[r_0] *= is_active;

            w_shift[2] = w_shift[1];
            w_shift[1] = w_shift[0];
            w_shift[0] = w0;

            Y0_curr[j] = C_delta0 * w0 + C_eps0 * Y0_next[j];
            Y0_curr[j] *= is_active;
            if (w0_row_ip1)
                w0_row_ip1[j] += X_ij + Y0_curr[j] + C_delta1 * w2 + C_delta2 * w1;

            // Per-lane anchor tracking — update best values per logical position
            if (in_band) {
                lane_first_active = Kokkos::Experimental::condition(
                    is_active > 0 && (Float)j < lane_first_active,
                    simd_t((Float)j), lane_first_active);
                lane_last_active = Kokkos::Experimental::condition(
                    is_active > 0 && (Float)j > lane_last_active,
                    simd_t((Float)j), lane_last_active);

                if (!capture_bands) {
                alignas(64) Float wMid_arr[simdWidth], w0_arr[simdWidth];
                simd_unchecked_store(wMid, wMid_arr, Kokkos::Experimental::simd_flag_default);
                simd_unchecked_store(w0, w0_arr, Kokkos::Experimental::simd_flag_default);

                for (int k = 0; k < activeCount; k++) {
                    int logical_j = j - scratch.seq_offset[k];
                    if (logical_j >= 0 && logical_j < (int)decoded[k]->size()) {
                        if (wMid_arr[k] > scratch.best_wMid[k][logical_j]) {
                            scratch.best_wMid[k][logical_j] = wMid_arr[k];
                            scratch.best_wEnd[k][logical_j] = w0_arr[k];
                            scratch.best_i[k][logical_j] = i;
                        }
                    }
                }
                }

                Kokkos::Experimental::simd_mask<Float> row_mask = (xdrop_score > row_best_metric) && ((Float)j < seq_end);
                row_best_metric = Kokkos::Experimental::condition(row_mask, xdrop_score, row_best_metric);
                row_best_j = Kokkos::Experimental::condition(row_mask, simd_t((Float)j), row_best_j);
            }
        }

        int64_t row_active_this;
        {
            Float row_active_arr[simdWidth];
            simd_unchecked_store(row_active_count, row_active_arr, Kokkos::Experimental::simd_flag_default);
            row_active_this = 0;
            for (int k = 0; k < activeCount; k++) {
                row_active_this += (int64_t)row_active_arr[k];
                active_cells += (int64_t)row_active_arr[k];
            }
            total_cells += activeCount * (hi - seq_start + 1);
        }

        Float first_arr[simdWidth], last_arr[simdWidth];
        simd_unchecked_store(lane_first_active, first_arr, Kokkos::Experimental::simd_flag_default);
        simd_unchecked_store(lane_last_active, last_arr, Kokkos::Experimental::simd_flag_default);

        int64_t row_total = (int64_t)activeCount * (hi - seq_start + 1);
        if (row_active_this > row_total) {
            std::ostringstream o;
            o << "BUG at i=" << i << ": row_active=" << row_active_this << " row_total=" << row_total
              << " lo=" << lo << " hi=" << hi << " seq_start=" << seq_start
              << " active_dp_width=" << scratch.active_dp_width << "\n";
            for (int k = 0; k < activeCount; k++)
                o << "  lane " << k << ": first=" << (int)first_arr[k] << " last=" << (int)last_arr[k] << "\n";
            std::cerr << o.str();
        }
        if (useXDrop) {
            // Re-offset: re-center the DP frame when the active band drifts
            // too far. Throttled to avoid rebuilding transposed_decoded every row.
            int bw = hi - lo;
            bool do_reoffset = (i >= XDROP_MIN_I * 2) && (i % XDROP_MIN_I == 0) &&
                               (i < profile.length) && (bw > scratch.active_dp_width / 2);
            if (do_reoffset) {
                reoffset(scratch, activeCount, decoded, profile, bg_probs_ptr,
                         first_arr, last_arr, row_best_j, row_best_metric, capture_bands);
                transposed_base = scratch.transposed_decoded.data() + 4 * simdWidth;
                bg_codon_probs_base = scratch.bg_codon_probs.data();
            }

            for (int k = 0; k < activeCount; k++) {
                int first = first_arr[k];
                int last  = last_arr[k];
                if (first <= last) {
                    int lane_lo = (int)scratch.seq_offset[k];
                    int lane_hi = (int)(scratch.seq_offset[k] + (int64_t)decoded[k]->size() - 1);
                    int seq_lo = std::max(first, lane_lo);
                    int seq_hi = std::min(last + 3, lane_hi);
                    scratch.fwd_band_lo[k][i + 1] = std::clamp(seq_lo, lane_lo, lane_hi);
                    scratch.fwd_band_hi[k][i + 1] = std::clamp(seq_hi, lane_lo, lane_hi);
                } else {
                    scratch.fwd_band_lo[k][i + 1] = INT_MAX;
                    scratch.fwd_band_hi[k][i + 1] = INT_MIN;
                }
            }
        }

#ifdef DEBUG_PRINT_ROWS
        logBandStats("#", i, activeCount, first_arr, last_arr, row_best_metric, decoded);
#endif

        std::swap(Y0_curr, Y0_next);
        std::fill(Y0_curr.begin(), Y0_curr.end(), simd_t(0.0));

        std::swap(scratch.W0_curr, scratch.W0_next);
        std::fill(scratch.W0_next.begin(), scratch.W0_next.end(), simd_t(0.0));
    }

    if (true) {
        std::ostringstream ss;
        ss << "# X-drop enabled: " << useXDrop << " (i, j) used " <<  active_ij << "/" << total_ij
           << " (" << (100.0 * active_ij / total_ij) << "%)"
           << "  SIMD lanes active: " << active_cells << "/" << total_cells
           << " (" << (100.0 * active_cells / total_cells) << "%)" << std::endl;
        std::cerr << ss.str();
    }

#ifdef ALIGN
    if (!capture_bands) {
    for (int j = 0; j < scratch.active_dp_width; j++) {
        const char* indices = (const char*)&transposed_base[(j + 1) * simdWidth];
        SimdFloat bg_raw = simdLookup(bg_probs_ptr, indices);

        simd_t bg_codon_emit_probs(bg_raw);

        left_side[j + 3] += (Float)(1 - BACKGROUND_FRAMESHIFT_RATE - BACKGROUND_FRAMESHIFT_RATE_2) * bg_codon_emit_probs * null_emit_3 * left_side[j];
        left_side[j + 1] += (Float)(BACKGROUND_FRAMESHIFT_RATE * 0.25) * null_emit_1 * left_side[j];
        left_side[j + 2] += (Float)(BACKGROUND_FRAMESHIFT_RATE_2 * 0.0625) * null_emit_2 * left_side[j];

        left_side[j] *= null_model_suffix[j];
        left_side[j] = Kokkos::max(left_side[j], Float(0.0));
    }
    for (int j = scratch.active_dp_width - 2; j >= 0; j--) {
        left_side[j] += left_side[j + 1];
    }

    {
    alignas(64) Float gather_buf[simdWidth];
    auto gather_w1 = [&](int i_row, int col) -> simd_t {
        return simdGather(scratch.W1.data(), gather_buf, activeCount, i_row, col);
    };
    auto gather_x = [&](int i_row, int col) -> simd_t {
        return simdGather(scratch.X.data(), gather_buf, activeCount, i_row, col);
    };
    auto scatter_w1 = [&](simd_t val, int i_row, int col) {
        simdScatter(val, scratch.W1.data(), gather_buf, activeCount, i_row, col);
    };

    for (int i = profile.length; i >= 0; i--) {
        simd_t opt_right_rolling = simd_t(0.0);

        bool ip1_avail = (i + 1 <= profile.length);
        for (int j = scratch.active_dp_width - 1; j >= 0; j--) {
            simd_t opt_succ = ip1_avail ? gather_w1(i + 1, j + 3) : simd_t(0);
            simd_t opt_down = ip1_avail ? gather_w1(i + 1, j) : simd_t(0);

            auto opt = Kokkos::max(opt_down, opt_right_rolling);
            opt = Kokkos::max(opt, left_side[j]);
            opt = Kokkos::max(opt, gather_x(i, j) + opt_succ);
            scatter_w1(opt, i, j);
            opt_right_rolling = opt;
        }
    }
    }
    }
#endif

    if (capture_bands) {
        if (args.out_lo && args.out_hi) {
            // Translate merged bands from pass-1's (possibly reoffset-shifted)
            // physical frame back to sequence coordinates for pass 2, which
            // resets seq_offset=0 and rebuilds transposed_decoded at offset 0.
            for (int i = 0; i <= profile.length + 1; ++i) {
                for (int k = 0; k < activeCount; k++) {
                    int mlo = std::min(scratch.rev_band_lo[k][i], scratch.fwd_band_lo[k][i]);
                    int mhi = std::max(scratch.rev_band_hi[k][i], scratch.fwd_band_hi[k][i]);
                    if (mlo == INT_MAX || mhi == INT_MIN || mlo > mhi) {
                        (*args.out_lo)[k][i] = INT_MAX;
                        (*args.out_hi)[k][i] = INT_MIN;
                        continue;
                    }
                    int seq_lo = mlo - (int)scratch.seq_offset[k];
                    int seq_hi = mhi - (int)scratch.seq_offset[k];
                    int len = (int)decoded[k]->size();
                    if (seq_hi < 0 || seq_lo > len - 1) {
                        (*args.out_lo)[k][i] = INT_MAX;
                        (*args.out_hi)[k][i] = INT_MIN;
                    } else {
                        (*args.out_lo)[k][i] = std::clamp(seq_lo, 0, len - 1);
                        (*args.out_hi)[k][i] = std::clamp(seq_hi, 0, len - 1);
                    }
                }
            }
        }
        return;
    }

    for (int idx = 0; idx < activeCount; idx++) {
        int len = (int)decoded[idx]->size();
        int offset = scratch.seq_offset[idx];

        scratch.opt_profile_position[idx].assign(len, AlignedSimilarity(-INFINITY));

        Float best_overall = -INFINITY;
        for (int j = 0; j < len; j++) {
            int phys_j = j + offset;
            Float best_prob = scratch.best_wMid[idx][j];
            if (best_prob > -INFINITY) {
                best_overall = std::max(best_overall, best_prob);
                scratch.opt_profile_position[idx][j] = {
                    best_prob,
                    scratch.best_i[idx][j],
                    phys_j,
                    (Float)scratch.best_wEnd[idx][j]
                };
            }
        }

        if (minProbRatio[idx] >= 0) {
            if (best_overall >= minProbRatio[idx]) {
                std::cerr << "has" << std::endl;
                std::ranges::sort(scratch.opt_profile_position[idx], std::greater<>());
                auto &aligned = scratch.aligned;
                aligned[idx].assign(decoded[idx]->size() + 0, false);

                for (auto &aligned_similarity : scratch.opt_profile_position[idx]) {
                    int logical_j = aligned_similarity.anchor2 - offset;

                    if (aligned_similarity.probRatio >= minProbRatio[idx] &&
                        !aligned[idx][logical_j]) {

                        addMidAnchored(idx, profile.length, similarities[idx], aligned_similarity.anchor1,
                                       aligned_similarity.anchor2,
                                       aligned_similarity.probRatio * scale /
                                           aligned_similarity.wEndAnchored,
                                       aligned_similarity.wEndAnchored, scratch);
                        auto &x = similarities[idx].back();
                        finishMidAnchored(idx, x, scratch);
                        x.anchor2 -= offset;

                        int startIdx = std::max(logical_j - 12 * profile.length, 0);
                        int endIdx = std::min(logical_j + 12 * profile.length, (int)decoded[idx]->size() + 0);
                        std::fill(aligned[idx].begin() + startIdx, aligned[idx].begin() + endIdx, true);
                        }
                }
            }
        } else {
            auto sel = *std::max_element(scratch.opt_profile_position[idx].begin(), scratch.opt_profile_position[idx].end());
            sel.anchor2 -= offset;
            AlignedSimilarity b = sel;
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
    bool hasSeeds = false;
    for (int i = 0; i < activeCount; i++) {
        if (!req[i].seeds.empty()) { hasSeeds = true; break; }
    }
    if (hasSeeds) {
        std::array<std::vector<std::vector<int>>, simdWidth> seed_j;
        for (int k = 0; k < activeCount; k++) {
            seed_j[k].resize(profile.length + 2);
            for (const auto& s : req[k].seeds) {
                int off = scratch.seq_offset[k];
                if (s[0] >= 0 && s[0] <= profile.length)
                    seed_j[k][s[0]].push_back(s[1] + off);
                if (s[2] >= 0 && s[2] <= profile.length)
                    seed_j[k][s[2]].push_back(s[3] + off);
            }
        }

        std::array<std::vector<int>, simdWidth> lane_lo, lane_hi;
        for (int k = 0; k < activeCount; k++) {
            lane_lo[k].resize(profile.length + 2);
            lane_hi[k].resize(profile.length + 2);
        }

        // Pass 1: seed-gated band discovery
        findSimilarities(sims, profile, decoded, minProbRatio, scratch, activeCount, true,
            DPArgs{.seed_j = &seed_j, .out_lo = &lane_lo, .out_hi = &lane_hi});

        // Pass 2: gated DP within discovered bands
        findSimilarities(sims, profile, decoded, minProbRatio, scratch, activeCount, false,
            DPArgs{.in_lo = &lane_lo, .in_hi = &lane_hi});
    } else {
        findSimilarities(sims, profile, decoded, minProbRatio, scratch, activeCount, false);
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

    int batchStride = simdWidth;
    static const char* lanes_env = getenv("DUMMER_MAX_LANES");
    if (lanes_env) batchStride = std::clamp(std::atoi(lanes_env), 1, simdWidth);
    std::vector<BatchJob> jobs;
    for (size_t i = 0; i < profiles.size(); ++i) {
        const auto &requests = allRequests[i];
        size_t n = requests.size();
        for (size_t start = 0; start < n; start += batchStride) {
            int count = std::min(static_cast<size_t>(batchStride), n - start);
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
    int iters = 0;
    do {
        if (++iters > 1000) return 1.0;
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

        int batchStride = simdWidth;
        static const char* lanes_env = getenv("DUMMER_MAX_LANES");
        if (lanes_env) batchStride = std::clamp(std::atoi(lanes_env), 1, simdWidth);
        int numBatches = (numOfSequences + batchStride - 1) / batchStride;
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

                int start = batchIdx * batchStride;
                int activeCount = std::min(static_cast<int>(batchStride), numOfSequences - start);

                std::array<std::vector<uint8_t>, simdWidth> decoded;
                std::array<std::vector<uint8_t>*, simdWidth> decodedPtrs = {};
                std::array<Float, simdWidth> minProbRatio;
                minProbRatio.fill(-2.0f);

                for (int lane = 0; lane < activeCount; ++lane) {
                    int trialIdx = start + lane;
                    // Core-independent deterministic seeding based on trial index
                    std::mt19937_64 trialRandGen(5489 + trialIdx);
                    char *seqBuf = localSeqs[lane].data();

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

                    for (int j = 0; j < border; ++j)
                        seqBuf[sequenceLength + j] = seqBuf[j];

                    decoded[lane] = decodeSequence(seqBuf, sequenceLength + border, alphabet, charToNumber);
                    decodedPtrs[lane] = &decoded[lane];
                }

                std::array<std::vector<AlignedSimilarity>, simdWidth> simsSIMD;
                findSimilarities(simsSIMD, profile, decodedPtrs, minProbRatio, threadScratch, activeCount, false);

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
        double alpha = probs[0] * (1 - INSERT1 - INSERT2 - DELETE1 - DELETE2);
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

        double delta = probs[2] * (1 - INSERT1 - INSERT2 - DELETE1 - DELETE2);;
        double epsilon = probs[3];
        if (i == p.length)
            break;

        double delta1 = probs[p.width + 2];
        double epsilon1 = probs[p.width + 3];
        double deltaFS1 = DELETE1; // simulate delete
        double deltaFS2 = DELETE2;
        p.values_v2.rbegin()->delta_prime[0] = delta * (1 - epsilon1);
        p.values_v2.rbegin()->delta_prime[1] = deltaFS1 * (1 - epsilon1);
        p.values_v2.rbegin()->delta_prime[2] = deltaFS2 * (1 - epsilon1);
        p.values_v2.rbegin()->epsilon_prime = epsilon * (1 - epsilon1) / (1 - epsilon);
        p.values_v2.rbegin()->enter_match_probability = (1 - alpha - alphaFS1 - alphaFS2 - delta - deltaFS1 - deltaFS2);

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
            val = 0;
            break;
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
        estimateK(p, bgProbs, &charVec[seqIdx], randomSeqLen, border, randomSeqNum, printVerbosity, threadPool, threadScratches);
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
                    if (sequence.target_profile.empty() || !strcmp(&charVec[p.nameIdx], sequence.target_profile.c_str())) {
#endif
                        Float minProbRatio =
                            (evalueOpt > 0)
                                ? (std::pow(p.gumbelKmidAnchored * totSequenceLength / evalueOpt,
                                            1.0 / 1.0 /* p.lambda */))
                                : -1;
                        if (verbosity > 1)
                            std::cerr << "Profile: " << &charVec[p.nameIdx] << "\n";

                        allRequests[j].push_back({sd, minProbRatio, sequence.seeds});
#ifdef PIPELINE_MODE
                    }
#endif
                }
            }
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
