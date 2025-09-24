// Author: Martin C. Frith 2025
// SPDX-License-Identifier: BSD-3-Clause

// Look, it's 2025 and we still can't have good, easy, portable SIMD?
// So write hardware-specific SIMD, like an animal

#if defined __SSSE3__
#include <immintrin.h>
#elif defined __ARM_NEON
#include <arm_neon.h>
#endif

#if defined __AVX512F__

typedef __m512  SimdFlt;
typedef __m512d SimdDbl;

const int simdFltLen = 16;
const int simdDblLen = 8;

SimdFlt simdFill(float  x) { return _mm512_set1_ps(x); }
SimdDbl simdFill(double x) { return _mm512_set1_pd(x); }

SimdFlt simdLoad(const float  *p) { return _mm512_loadu_ps(p); }
SimdDbl simdLoad(const double *p) { return _mm512_loadu_pd(p); }

void simdStore(float  *p, SimdFlt x) { _mm512_storeu_ps(p, x); }
void simdStore(double *p, SimdDbl x) { _mm512_storeu_pd(p, x); }

SimdFlt simdMax(SimdFlt x, SimdFlt y) { return _mm512_max_ps(x, y); }
SimdDbl simdMax(SimdDbl x, SimdDbl y) { return _mm512_max_pd(x, y); }

SimdFlt simdAdd(SimdFlt x, SimdFlt y) { return _mm512_add_ps(x, y); }
SimdDbl simdAdd(SimdDbl x, SimdDbl y) { return _mm512_add_pd(x, y); }

SimdFlt simdMul(SimdFlt x, SimdFlt y) { return _mm512_mul_ps(x, y); }
SimdDbl simdMul(SimdDbl x, SimdDbl y) { return _mm512_mul_pd(x, y); }

__mmask16 simdGe(SimdFlt x, SimdFlt y) {  // greater or equal: x >= y
  return _mm512_cmp_ps_mask(x, y, 29);
}
__mmask8 simdGe(SimdDbl x, SimdDbl y) {  // greater or equal: x >= y
  return _mm512_cmp_pd_mask(x, y, 29);
}

SimdFlt simdLowItem(SimdFlt x) {
  return _mm512_broadcastss_ps(_mm512_castps512_ps128(x));
}
SimdDbl simdLowItem(SimdDbl x) {
  return _mm512_broadcastsd_pd(_mm512_castpd512_pd128(x));
}

SimdFlt simdHighItem(SimdFlt x) {
  return _mm512_permutexvar_ps(_mm512_set1_epi32(-1), x);
}
SimdDbl simdHighItem(SimdDbl x) {
  return _mm512_permutexvar_pd(_mm512_set1_epi32(-1), x);
}

SimdFlt simdLookup(const float  *v, const char *i) {
  return _mm512_set_ps(v[i[15]], v[i[14]], v[i[13]], v[i[12]],
		       v[i[11]], v[i[10]], v[i[9]], v[i[8]],
		       v[i[7]], v[i[6]], v[i[5]], v[i[4]],
		       v[i[3]], v[i[2]], v[i[1]], v[i[0]]);
}
SimdDbl simdLookup(const double *v, const char *i) {
  return _mm512_set_pd(v[i[7]], v[i[6]], v[i[5]], v[i[4]],
		       v[i[3]], v[i[2]], v[i[1]], v[i[0]]);
}

SimdFlt simdLookup(SimdFlt v, const char *i) {
  const __m128i *j = (const __m128i *)i;
  return _mm512_permutevar_ps(v, _mm512_cvtepi8_epi32(_mm_loadu_si128(j)));
}
SimdDbl simdLookup(SimdDbl v, const char *i) {
  return _mm512_permutexvar_pd(_mm512_cvtepi8_epi64(_mm_loadu_si64(i)), v);
}

SimdFlt simdLookup(SimdFlt a, SimdFlt b, const char *i) {
  __m512i k = _mm512_cvtepi8_epi32(_mm_loadu_si128((const __m128i *)i));
  return _mm512_permutex2var_ps(a, k, b);
}
SimdDbl simdLookup(SimdDbl a, SimdDbl b, const char *i) { return a; }  // unimp

// shift high item of zOld into low position of zNew
SimdFlt simdShiftFwd(SimdFlt zOld, SimdFlt zNew) {
  return _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_castps_si512(zNew),
						 _mm512_castps_si512(zOld),
						 15));
}
SimdDbl simdShiftFwd(SimdDbl zOld, SimdDbl zNew) {
  return _mm512_castsi512_pd(_mm512_alignr_epi64(_mm512_castpd_si512(zNew),
						 _mm512_castpd_si512(zOld),
						 7));
}

// shift low item of zOld into high position of zNew
SimdFlt simdShiftRev(SimdFlt zOld, SimdFlt zNew) {
  return _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_castps_si512(zOld),
						 _mm512_castps_si512(zNew),
						 1));
}
SimdDbl simdShiftRev(SimdDbl zOld, SimdDbl zNew) {
  return _mm512_castsi512_pd(_mm512_alignr_epi64(_mm512_castpd_si512(zOld),
						 _mm512_castpd_si512(zNew),
						 1));
}

SimdFlt simdCumulateFwd(SimdFlt v, SimdFlt c) {
  __m512i z = _mm512_setzero_si512();
  SimdFlt a;
  //a = _mm512_castsi512_ps(_mm512_maskz_alignr_epi32(0xFFFE, v, v, 15);
  a = _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_castps_si512(v), z, 15));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_castps_si512(v), z, 14));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_castps_si512(v), z, 12));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_castps_si512(v), z, 8));
  return simdAdd(v, simdMul(c, a));
}
SimdDbl simdCumulateFwd(SimdDbl v, SimdDbl c) {
  __m512i z = _mm512_setzero_si512();
  SimdDbl a;
  a = _mm512_castsi512_pd(_mm512_alignr_epi64(_mm512_castpd_si512(v), z, 7));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_pd(_mm512_alignr_epi64(_mm512_castpd_si512(v), z, 6));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_pd(_mm512_alignr_epi64(_mm512_castpd_si512(v), z, 4));
  return simdAdd(v, simdMul(c, a));
}

SimdFlt simdCumulateRev(SimdFlt v, SimdFlt c) {
  __m512i z = _mm512_setzero_si512();
  SimdFlt a;
  a = _mm512_castsi512_ps(_mm512_alignr_epi32(z, _mm512_castps_si512(v), 1));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_ps(_mm512_alignr_epi32(z, _mm512_castps_si512(v), 2));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_ps(_mm512_alignr_epi32(z, _mm512_castps_si512(v), 4));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_ps(_mm512_alignr_epi32(z, _mm512_castps_si512(v), 8));
  return simdAdd(v, simdMul(c, a));
}
SimdDbl simdCumulateRev(SimdDbl v, SimdDbl c) {
  __m512i z = _mm512_setzero_si512();
  SimdDbl a;
  a = _mm512_castsi512_pd(_mm512_alignr_epi64(z, _mm512_castpd_si512(v), 1));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_pd(_mm512_alignr_epi64(z, _mm512_castpd_si512(v), 2));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm512_castsi512_pd(_mm512_alignr_epi64(z, _mm512_castpd_si512(v), 4));
  return simdAdd(v, simdMul(c, a));
}

#elif defined __AVX2__

typedef __m256  SimdFlt;
typedef __m256d SimdDbl;

const int simdFltLen = 8;
const int simdDblLen = 4;

SimdFlt simdFill(float  x) { return _mm256_set1_ps(x); }
SimdDbl simdFill(double x) { return _mm256_set1_pd(x); }

SimdFlt simdLoad(const float  *p) { return _mm256_loadu_ps(p); }
SimdDbl simdLoad(const double *p) { return _mm256_loadu_pd(p); }

void simdStore(float  *p, SimdFlt x) { _mm256_storeu_ps(p, x); }
void simdStore(double *p, SimdDbl x) { _mm256_storeu_pd(p, x); }

SimdFlt simdMax(SimdFlt x, SimdFlt y) { return _mm256_max_ps(x, y); }
SimdDbl simdMax(SimdDbl x, SimdDbl y) { return _mm256_max_pd(x, y); }

SimdFlt simdAdd(SimdFlt x, SimdFlt y) { return _mm256_add_ps(x, y); }
SimdDbl simdAdd(SimdDbl x, SimdDbl y) { return _mm256_add_pd(x, y); }

SimdFlt simdMul(SimdFlt x, SimdFlt y) { return _mm256_mul_ps(x, y); }
SimdDbl simdMul(SimdDbl x, SimdDbl y) { return _mm256_mul_pd(x, y); }

int simdGe(SimdFlt x, SimdFlt y) {  // greater or equal: x >= y
  return _mm256_movemask_ps(_mm256_cmp_ps(x, y, 29));
}
int simdGe(SimdDbl x, SimdDbl y) {  // greater or equal: x >= y
  return _mm256_movemask_pd(_mm256_cmp_pd(x, y, 29));
}

SimdFlt simdLowItem(SimdFlt x) {
  return _mm256_broadcastss_ps(_mm256_castps256_ps128(x));
}
SimdDbl simdLowItem(SimdDbl x) {
  return _mm256_broadcastsd_pd(_mm256_castpd256_pd128(x));
}

SimdFlt simdHighItem(SimdFlt x) {
  return _mm256_permutevar8x32_ps(x, _mm256_set1_epi32(-1));
}
SimdDbl simdHighItem(SimdDbl x) { return _mm256_permute4x64_pd(x, 255); }

SimdFlt simdLookup(const float  *v, const char *i) {
  return _mm256_set_ps(v[i[7]], v[i[6]], v[i[5]], v[i[4]],
		       v[i[3]], v[i[2]], v[i[1]], v[i[0]]);
}
SimdDbl simdLookup(const double *v, const char *i) {
  return _mm256_set_pd(v[i[3]], v[i[2]], v[i[1]], v[i[0]]);
}

SimdFlt simdLookup(SimdFlt v, const char *i) {
  return _mm256_permutevar_ps(v, _mm256_cvtepi8_epi32(_mm_loadu_si64(i)));
}
SimdDbl simdLookup(SimdDbl v, const char *i) { return v; }  // not implemented
SimdFlt simdLookup(SimdFlt a, SimdFlt b, const char *i) { return a; }
SimdDbl simdLookup(SimdDbl a, SimdDbl b, const char *i) { return a; }

// shift high item of zOld into low position of zNew
SimdFlt simdShiftFwd(SimdFlt zOld, SimdFlt zNew) {
  SimdFlt a = _mm256_permute2f128_ps(zNew, zOld, 3);
  return _mm256_castsi256_ps(_mm256_alignr_epi8(_mm256_castps_si256(zNew),
						_mm256_castps_si256(a), 12));
}
SimdDbl simdShiftFwd(SimdDbl zOld, SimdDbl zNew) {
  SimdDbl a = _mm256_permute2f128_pd(zNew, zOld, 3);
  return _mm256_castsi256_pd(_mm256_alignr_epi8(_mm256_castpd_si256(zNew),
						_mm256_castpd_si256(a), 8));
}

// shift low item of zOld into high position of zNew
SimdFlt simdShiftRev(SimdFlt zOld, SimdFlt zNew) {
  SimdFlt a = _mm256_permute2f128_ps(zOld, zNew, 3);
  return _mm256_castsi256_ps(_mm256_alignr_epi8(_mm256_castps_si256(a),
						_mm256_castps_si256(zNew), 4));
}
SimdDbl simdShiftRev(SimdDbl zOld, SimdDbl zNew) {
  SimdDbl a = _mm256_permute2f128_pd(zOld, zNew, 3);
  return _mm256_castsi256_pd(_mm256_alignr_epi8(_mm256_castpd_si256(a),
						_mm256_castpd_si256(zNew), 8));
}

SimdFlt simdCumulateFwd(SimdFlt v, SimdFlt c) {
  SimdFlt a = _mm256_permute2f128_ps(v, v, 8);
  a = _mm256_castsi256_ps(_mm256_alignr_epi8(_mm256_castps_si256(v),
					     _mm256_castps_si256(a), 12));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm256_permute2f128_ps(v, v, 8);
  a = _mm256_castsi256_ps(_mm256_alignr_epi8(_mm256_castps_si256(v),
					     _mm256_castps_si256(a), 8));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm256_permute2f128_ps(v, v, 8);
  return simdAdd(v, simdMul(c, a));
}
SimdDbl simdCumulateFwd(SimdDbl v, SimdDbl c) {
  SimdDbl a = _mm256_permute2f128_pd(v, v, 8);
  a = _mm256_castsi256_pd(_mm256_alignr_epi8(_mm256_castpd_si256(v),
					     _mm256_castpd_si256(a), 8));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm256_permute2f128_pd(v, v, 8);
  return simdAdd(v, simdMul(c, a));
}

SimdFlt simdCumulateRev(SimdFlt v, SimdFlt c) {
  SimdFlt a = _mm256_permute2f128_ps(v, v, 0x83);
  a = _mm256_castsi256_ps(_mm256_alignr_epi8(_mm256_castps_si256(a),
					     _mm256_castps_si256(v), 4));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm256_permute2f128_ps(v, v, 0x83);
  a = _mm256_castsi256_ps(_mm256_alignr_epi8(_mm256_castps_si256(a),
					     _mm256_castps_si256(v), 8));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm256_permute2f128_ps(v, v, 0x83);
  return simdAdd(v, simdMul(c, a));
}
SimdDbl simdCumulateRev(SimdDbl v, SimdDbl c) {
  SimdDbl a = _mm256_permute2f128_pd(v, v, 0x83);
  a = _mm256_castsi256_pd(_mm256_alignr_epi8(_mm256_castpd_si256(a),
					     _mm256_castpd_si256(v), 8));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm256_permute2f128_pd(v, v, 0x83);
  return simdAdd(v, simdMul(c, a));
}

#elif defined __SSSE3__

typedef __m128  SimdFlt;
typedef __m128d SimdDbl;

const int simdFltLen = 4;
const int simdDblLen = 2;

SimdFlt simdFill(float  x) { return _mm_set1_ps(x); }
SimdDbl simdFill(double x) { return _mm_set1_pd(x); }

SimdFlt simdLoad(const float  *p) { return _mm_loadu_ps(p); }
SimdDbl simdLoad(const double *p) { return _mm_loadu_pd(p); }

void simdStore(float  *p, SimdFlt x) { _mm_storeu_ps(p, x); }
void simdStore(double *p, SimdDbl x) { _mm_storeu_pd(p, x); }

SimdFlt simdMax(SimdFlt x, SimdFlt y) { return _mm_max_ps(x, y); }
SimdDbl simdMax(SimdDbl x, SimdDbl y) { return _mm_max_pd(x, y); }

SimdFlt simdAdd(SimdFlt x, SimdFlt y) { return _mm_add_ps(x, y); }
SimdDbl simdAdd(SimdDbl x, SimdDbl y) { return _mm_add_pd(x, y); }

SimdFlt simdMul(SimdFlt x, SimdFlt y) { return _mm_mul_ps(x, y); }
SimdDbl simdMul(SimdDbl x, SimdDbl y) { return _mm_mul_pd(x, y); }

int simdGe(SimdFlt x, SimdFlt y) {  // greater or equal: x >= y
  return _mm_movemask_ps(_mm_cmpge_ps(x, y));
}
int simdGe(SimdDbl x, SimdDbl y) {  // greater or equal: x >= y
  return _mm_movemask_pd(_mm_cmpge_pd(x, y));
}

SimdFlt simdLowItem(SimdFlt x) { return _mm_shuffle_ps(x, x, 0); }
SimdDbl simdLowItem(SimdDbl x) { return _mm_shuffle_pd(x, x, 0); }

SimdFlt simdHighItem(SimdFlt x) { return _mm_shuffle_ps(x, x, 255); }
SimdDbl simdHighItem(SimdDbl x) { return _mm_shuffle_pd(x, x, 3); }

SimdFlt simdLookup(const float  *v, const char *i) {
  return _mm_set_ps(v[i[3]], v[i[2]], v[i[1]], v[i[0]]);
}
SimdDbl simdLookup(const double *v, const char *i) {
  return _mm_set_pd(v[i[1]], v[i[0]]);
}

SimdFlt simdLookup(SimdFlt v, const char *i) { return v; }  // not implemented
SimdDbl simdLookup(SimdDbl v, const char *i) { return v; }  // not implemented
SimdFlt simdLookup(SimdFlt a, SimdFlt b, const char *i) { return a; }
SimdDbl simdLookup(SimdDbl a, SimdDbl b, const char *i) { return a; }

// shift high item of zOld into low position of zNew
SimdFlt simdShiftFwd(SimdFlt zOld, SimdFlt zNew) {  // SSSE3
  return _mm_castsi128_ps(_mm_alignr_epi8(_mm_castps_si128(zNew),
					  _mm_castps_si128(zOld), 12));
}
SimdDbl simdShiftFwd(SimdDbl zOld, SimdDbl zNew) {  // SSSE3
  return _mm_castsi128_pd(_mm_alignr_epi8(_mm_castpd_si128(zNew),
					  _mm_castpd_si128(zOld), 8));
}

// shift low item of zOld into high position of zNew
SimdFlt simdShiftRev(SimdFlt zOld, SimdFlt zNew) {  // SSSE3
  return _mm_castsi128_ps(_mm_alignr_epi8(_mm_castps_si128(zOld),
					  _mm_castps_si128(zNew), 4));
}
SimdDbl simdShiftRev(SimdDbl zOld, SimdDbl zNew) {  // SSSE3
  return _mm_castsi128_pd(_mm_alignr_epi8(_mm_castpd_si128(zOld),
					  _mm_castpd_si128(zNew), 8));
}

SimdFlt simdCumulateFwd(SimdFlt v, SimdFlt c) {
  SimdFlt a = _mm_castsi128_ps(_mm_slli_si128(_mm_castps_si128(v), 4));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm_castsi128_ps(_mm_slli_si128(_mm_castps_si128(v), 8));
  return simdAdd(v, simdMul(c, a));
}
SimdDbl simdCumulateFwd(SimdDbl v, SimdDbl c) {
  SimdDbl a = _mm_castsi128_pd(_mm_slli_si128(_mm_castpd_si128(v), 8));
  return simdAdd(v, simdMul(c, a));
}

SimdFlt simdCumulateRev(SimdFlt v, SimdFlt c) {
  SimdFlt a = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v), 4));
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(v), 8));
  return simdAdd(v, simdMul(c, a));
}
SimdDbl simdCumulateRev(SimdDbl v, SimdDbl c) {
  SimdDbl a = _mm_castsi128_pd(_mm_srli_si128(_mm_castpd_si128(v), 8));
  return simdAdd(v, simdMul(c, a));
}

#elif defined __ARM_NEON

typedef float32x4_t SimdFlt;
typedef float64x2_t SimdDbl;

const int simdFltLen = 4;
const int simdDblLen = 2;

SimdFlt simdFill(float  x) { return vdupq_n_f32(x); }
SimdDbl simdFill(double x) { return vdupq_n_f64(x); }

SimdFlt simdLoad(const float  *p) { return vld1q_f32(p); }
SimdDbl simdLoad(const double *p) { return vld1q_f64(p); }

void simdStore(float  *p, SimdFlt x) { vst1q_f32(p, x); }
void simdStore(double *p, SimdDbl x) { vst1q_f64(p, x); }

SimdFlt simdMax(SimdFlt x, SimdFlt y) { return vmaxq_f32(x, y); }
SimdDbl simdMax(SimdDbl x, SimdDbl y) { return vmaxq_f64(x, y); }

SimdFlt simdAdd(SimdFlt x, SimdFlt y) { return vaddq_f32(x, y); }
SimdDbl simdAdd(SimdDbl x, SimdDbl y) { return vaddq_f64(x, y); }

SimdFlt simdMul(SimdFlt x, SimdFlt y) { return vmulq_f32(x, y); }
SimdDbl simdMul(SimdDbl x, SimdDbl y) { return vmulq_f64(x, y); }

uint64_t simdGe(SimdFlt x, SimdFlt y) {  // greater or equal: x >= y
  return vget_lane_u64(vreinterpret_u64_u16(vmovn_u32(vcgeq_f32(x, y))), 0);
}
uint64_t simdGe(SimdDbl x, SimdDbl y) {  // greater or equal: x >= y
  return vget_lane_u64(vreinterpret_u64_u32(vmovn_u64(vcgeq_f64(x, y))), 0);
}

SimdFlt simdLowItem(SimdFlt x) { return vdupq_laneq_f32(x, 0); }
SimdDbl simdLowItem(SimdDbl x) { return vdupq_laneq_f64(x, 0); }

SimdFlt simdHighItem(SimdFlt x) { return vdupq_laneq_f32(x, 3); }
SimdDbl simdHighItem(SimdDbl x) { return vdupq_laneq_f64(x, 1); }

SimdFlt simdLookup(const float *v, const char *i) {
  float a[] = { v[i[0]], v[i[1]], v[i[2]], v[i[3]] };
  return simdLoad(a);
}
SimdDbl simdLookup(const double *v, const char *i) {
  double a[] = { v[i[0]], v[i[1]] };
  return simdLoad(a);
}

SimdFlt simdLookup(SimdFlt v, const char *i) { return v; }  // not implemented
SimdDbl simdLookup(SimdDbl v, const char *i) { return v; }  // not implemented
SimdFlt simdLookup(SimdFlt a, SimdFlt b, const char *i) { return a; }
SimdDbl simdLookup(SimdDbl a, SimdDbl b, const char *i) { return a; }

// shift high item of zOld into low position of zNew
SimdFlt simdShiftFwd(SimdFlt zOld, SimdFlt zNew) {
  return vextq_f32(zOld, zNew, 3);
}
SimdDbl simdShiftFwd(SimdDbl zOld, SimdDbl zNew) {
  return vextq_f64(zOld, zNew, 1);
}

// shift low item of zOld into high position of zNew
SimdFlt simdShiftRev(SimdFlt zOld, SimdFlt zNew) {
  return vextq_f32(zNew, zOld, 1);
}
SimdDbl simdShiftRev(SimdDbl zOld, SimdDbl zNew) {
  return vextq_f64(zNew, zOld, 1);
}

SimdFlt simdCumulateFwd(SimdFlt v, SimdFlt c) {
  SimdFlt z = vdupq_n_f32(0);
  SimdFlt a = vextq_f32(z, v, 3);
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = vextq_f32(z, v, 2);
  return simdAdd(v, simdMul(c, a));
}
SimdDbl simdCumulateFwd(SimdDbl v, SimdDbl c) {
  SimdDbl z = vdupq_n_f64(0);
  SimdDbl a = vextq_f64(z, v, 1);
  return simdAdd(v, simdMul(c, a));
}

SimdFlt simdCumulateRev(SimdFlt v, SimdFlt c) {
  SimdFlt z = vdupq_n_f32(0);
  SimdFlt a = vextq_f32(v, z, 1);
  v = simdAdd(v, simdMul(c, a));
  c = simdMul(c, c);
  a = vextq_f32(v, z, 2);
  return simdAdd(v, simdMul(c, a));
}
SimdDbl simdCumulateRev(SimdDbl v, SimdDbl c) {
  SimdDbl z = vdupq_n_f64(0);
  SimdDbl a = vextq_f64(v, z, 1);
  return simdAdd(v, simdMul(c, a));
}

#else

typedef float  SimdFlt;
typedef double SimdDbl;

const int simdFltLen = 1;
const int simdDblLen = 1;

SimdFlt simdFill(float  x) { return x; }
SimdDbl simdFill(double x) { return x; }

SimdFlt simdLoad(const float  *p) { return *p; }
SimdDbl simdLoad(const double *p) { return *p; }

void simdStore(float  *p, SimdFlt x) { *p = x; }
void simdStore(double *p, SimdDbl x) { *p = x; }

SimdFlt simdMax(SimdFlt x, SimdFlt y) { return x > y ? x : y; }
SimdDbl simdMax(SimdDbl x, SimdDbl y) { return x > y ? x : y; }

SimdFlt simdAdd(SimdFlt x, SimdFlt y) { return x + y; }
SimdDbl simdAdd(SimdDbl x, SimdDbl y) { return x + y; }

SimdFlt simdMul(SimdFlt x, SimdFlt y) { return x * y; }
SimdDbl simdMul(SimdDbl x, SimdDbl y) { return x * y; }

bool simdGe(SimdFlt x, SimdFlt y) { return x >= y; }
bool simdGe(SimdDbl x, SimdDbl y) { return x >= y; }

SimdFlt simdLowItem(SimdFlt x) { return x; }
SimdDbl simdLowItem(SimdDbl x) { return x; }

SimdFlt simdHighItem(SimdFlt x) { return x; }
SimdDbl simdHighItem(SimdDbl x) { return x; }

SimdFlt simdLookup(const float  *v, const char *i) { return v[*i]; }
SimdDbl simdLookup(const double *v, const char *i) { return v[*i]; }

SimdFlt simdLookup(SimdFlt v, const char *i) { return v; }  // not implemented
SimdDbl simdLookup(SimdDbl v, const char *i) { return v; }  // not implemented
SimdFlt simdLookup(SimdFlt a, SimdFlt b, const char *i) { return a; }
SimdDbl simdLookup(SimdDbl a, SimdDbl b, const char *i) { return a; }

SimdFlt simdShiftFwd(SimdFlt zOld, SimdFlt zNew) { return zOld; }
SimdDbl simdShiftFwd(SimdDbl zOld, SimdDbl zNew) { return zOld; }

SimdFlt simdShiftRev(SimdFlt zOld, SimdFlt zNew) { return zOld; }
SimdDbl simdShiftRev(SimdDbl zOld, SimdDbl zNew) { return zOld; }

SimdFlt simdCumulateFwd(SimdFlt v, SimdFlt c) { return v; }
SimdDbl simdCumulateFwd(SimdDbl v, SimdDbl c) { return v; }

SimdFlt simdCumulateRev(SimdFlt v, SimdFlt c) { return v; }
SimdDbl simdCumulateRev(SimdDbl v, SimdDbl c) { return v; }

#endif
