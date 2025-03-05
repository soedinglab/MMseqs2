// Author: Martin C. Frith 2019
// SPDX-License-Identifier: MPL-2.0

#ifndef MCF_SIMD_HH
#define MCF_SIMD_HH

#if defined __SSE4_1__
#include <immintrin.h>
#elif defined __ARM_NEON
#include <arm_neon.h>
#endif

#include <stddef.h>  // size_t

namespace mcf {

#if defined __AVX2__

typedef __m256i SimdInt;
typedef __m256i SimdUint1;
typedef __m256d SimdDbl;

const int simdBytes = 32;

static inline SimdInt simdZero() {
  return _mm256_setzero_si256();
}

static inline SimdInt simdZero1() {
  return _mm256_setzero_si256();
}

static inline SimdDbl simdZeroDbl() {
  return _mm256_setzero_pd();
}

static inline SimdInt simdOnes1() {
  return _mm256_set1_epi32(-1);
}

static inline SimdInt simdLoad(const void *p) {
  return _mm256_loadu_si256((const SimdInt *)p);
}

static inline SimdInt simdLoad1(const void *p) {
  return _mm256_loadu_si256((const SimdInt *)p);
}

static inline SimdDbl simdLoadDbl(const double *p) {
  return _mm256_loadu_pd(p);
}

static inline void simdStore(void *p, SimdInt x) {
  _mm256_storeu_si256((SimdInt *)p, x);
}

static inline void simdStore1(void *p, SimdInt x) {
  _mm256_storeu_si256((SimdInt *)p, x);
}

static inline void simdStoreDbl(double *p, SimdDbl x) {
  _mm256_storeu_pd(p, x);
}

static inline SimdInt simdOr1(SimdInt x, SimdInt y) {
  return _mm256_or_si256(x, y);
}

static inline SimdInt simdBlend(SimdInt x, SimdInt y, SimdInt mask) {
  return _mm256_blendv_epi8(x, y, mask);
}

const int simdLen = 8;
const int simdDblLen = 4;

static inline SimdInt simdSet(int i7, int i6, int i5, int i4,
			      int i3, int i2, int i1, int i0) {
  return _mm256_set_epi32(i7, i6, i5, i4, i3, i2, i1, i0);
}

static inline SimdInt simdSet1(char jF, char jE, char jD, char jC,
			       char jB, char jA, char j9, char j8,
			       char j7, char j6, char j5, char j4,
			       char j3, char j2, char j1, char j0,
			       char iF, char iE, char iD, char iC,
			       char iB, char iA, char i9, char i8,
			       char i7, char i6, char i5, char i4,
			       char i3, char i2, char i1, char i0) {
  return _mm256_set_epi8(jF, jE, jD, jC, jB, jA, j9, j8,
			 j7, j6, j5, j4, j3, j2, j1, j0,
			 iF, iE, iD, iC, iB, iA, i9, i8,
			 i7, i6, i5, i4, i3, i2, i1, i0);
}

static inline SimdDbl simdSetDbl(double i3, double i2, double i1, double i0) {
  return _mm256_set_pd(i3, i2, i1, i0);
}

static inline SimdInt simdFill(int x) {
  return _mm256_set1_epi32(x);
}

static inline SimdInt simdFill1(char x) {
  return _mm256_set1_epi8(x);
}

static inline SimdDbl simdFillDbl(double x) {
  return _mm256_set1_pd(x);
}

static inline SimdInt simdGt(SimdInt x, SimdInt y) {
  return _mm256_cmpgt_epi32(x, y);
}

static inline SimdInt simdGe1(SimdInt x, SimdInt y) {
  return _mm256_cmpeq_epi8(_mm256_min_epu8(x, y), y);
}

static inline SimdInt simdAdd(SimdInt x, SimdInt y) {
  return _mm256_add_epi32(x, y);
}

static inline SimdInt simdAdd1(SimdInt x, SimdInt y) {
  return _mm256_add_epi8(x, y);
}

static inline SimdInt simdAdds1(SimdInt x, SimdInt y) {
  return _mm256_adds_epu8(x, y);
}

static inline SimdDbl simdAddDbl(SimdDbl x, SimdDbl y) {
  return _mm256_add_pd(x, y);
}

static inline SimdInt simdSub(SimdInt x, SimdInt y) {
  return _mm256_sub_epi32(x, y);
}

static inline SimdInt simdSub1(SimdInt x, SimdInt y) {
  return _mm256_sub_epi8(x, y);
}

static inline SimdDbl simdMulDbl(SimdDbl x, SimdDbl y) {
  return _mm256_mul_pd(x, y);
}

static inline SimdInt simdQuadruple1(SimdInt x) {
  return _mm256_slli_epi32(x, 2);
}

static inline SimdInt simdMax(SimdInt x, SimdInt y) {
  return _mm256_max_epi32(x, y);
}

static inline SimdInt simdMin1(SimdInt x, SimdInt y) {
  return _mm256_min_epu8(x, y);
}

static inline int simdHorizontalMax(SimdInt x) {
  __m128i z = _mm256_castsi256_si128(x);
  z = _mm_max_epi32(z, _mm256_extracti128_si256(x, 1));
  z = _mm_max_epi32(z, _mm_shuffle_epi32(z, 0x4E));
  z = _mm_max_epi32(z, _mm_shuffle_epi32(z, 0xB1));
  return _mm_cvtsi128_si32(z);
}

static inline int simdHorizontalMin1(SimdInt x) {
  __m128i z = _mm256_castsi256_si128(x);
  z = _mm_min_epu8(z, _mm256_extracti128_si256(x, 1));
  z = _mm_min_epu8(z, _mm_srli_epi16(z, 8));
  z = _mm_minpos_epu16(z);
  return _mm_extract_epi16(z, 0);
}

static inline double simdHorizontalAddDbl(SimdDbl x) {
  __m128d z = _mm256_castpd256_pd128(x);
  z = _mm_add_pd(z, _mm256_extractf128_pd(x, 1));
  return _mm_cvtsd_f64(_mm_hadd_pd(z, z));
}

static inline SimdInt simdChoose1(SimdInt items, SimdInt choices) {
  return _mm256_shuffle_epi8(items, choices);
}

#elif defined __SSE4_1__

typedef __m128i SimdInt;
typedef __m128i SimdUint1;
typedef __m128d SimdDbl;

const int simdBytes = 16;

static inline SimdInt simdZero() {
  return _mm_setzero_si128();
}

static inline SimdInt simdZero1() {
  return _mm_setzero_si128();
}

static inline SimdDbl simdZeroDbl() {
  return _mm_setzero_pd();
}

static inline SimdInt simdOnes1() {
  return _mm_set1_epi32(-1);
}

static inline SimdInt simdLoad(const void *p) {
  return _mm_loadu_si128((const SimdInt *)p);
}

static inline SimdInt simdLoad1(const void *p) {
  return _mm_loadu_si128((const SimdInt *)p);
}

static inline SimdDbl simdLoadDbl(const double *p) {
  return _mm_loadu_pd(p);
}

static inline void simdStore(void *p, SimdInt x) {
  _mm_storeu_si128((SimdInt *)p, x);
}

static inline void simdStore1(void *p, SimdInt x) {
  _mm_storeu_si128((SimdInt *)p, x);
}

static inline void simdStoreDbl(double *p, SimdDbl x) {
  _mm_storeu_pd(p, x);
}

static inline SimdInt simdOr1(SimdInt x, SimdInt y) {
  return _mm_or_si128(x, y);
}

static inline SimdInt simdBlend(SimdInt x, SimdInt y, SimdInt mask) {
  return _mm_blendv_epi8(x, y, mask);  // SSE4.1
}

const int simdLen = 4;
const int simdDblLen = 2;

static inline SimdInt simdSet(int i3, int i2, int i1, int i0) {
  return _mm_set_epi32(i3, i2, i1, i0);
}

static inline SimdInt simdSet1(char iF, char iE, char iD, char iC,
			       char iB, char iA, char i9, char i8,
			       char i7, char i6, char i5, char i4,
			       char i3, char i2, char i1, char i0) {
  return _mm_set_epi8(iF, iE, iD, iC, iB, iA, i9, i8,
		      i7, i6, i5, i4, i3, i2, i1, i0);
}

static inline SimdDbl simdSetDbl(double i1, double i0) {
  return _mm_set_pd(i1, i0);
}

static inline SimdInt simdFill(int x) {
  return _mm_set1_epi32(x);
}

static inline SimdInt simdFill1(char x) {
  return _mm_set1_epi8(x);
}

static inline SimdDbl simdFillDbl(double x) {
  return _mm_set1_pd(x);
}

static inline SimdInt simdGt(SimdInt x, SimdInt y) {
  return _mm_cmpgt_epi32(x, y);
}

static inline SimdInt simdGe1(SimdInt x, SimdInt y) {
  return _mm_cmpeq_epi8(_mm_min_epu8(x, y), y);
}

static inline SimdInt simdAdd(SimdInt x, SimdInt y) {
  return _mm_add_epi32(x, y);
}

static inline SimdInt simdAdd1(SimdInt x, SimdInt y) {
  return _mm_add_epi8(x, y);
}

static inline SimdInt simdAdds1(SimdInt x, SimdInt y) {
  return _mm_adds_epu8(x, y);
}

static inline SimdDbl simdAddDbl(SimdDbl x, SimdDbl y) {
  return _mm_add_pd(x, y);
}

static inline SimdInt simdSub(SimdInt x, SimdInt y) {
  return _mm_sub_epi32(x, y);
}

static inline SimdInt simdSub1(SimdInt x, SimdInt y) {
  return _mm_sub_epi8(x, y);
}

static inline SimdDbl simdMulDbl(SimdDbl x, SimdDbl y) {
  return _mm_mul_pd(x, y);
}

static inline SimdInt simdQuadruple1(SimdInt x) {
  return _mm_slli_epi32(x, 2);
}

static inline SimdInt simdMax(SimdInt x, SimdInt y) {
  return _mm_max_epi32(x, y);  // SSE4.1
}

static inline SimdInt simdMin1(SimdInt x, SimdInt y) {
  return _mm_min_epu8(x, y);
}

static inline int simdHorizontalMax(SimdInt x) {
  x = simdMax(x, _mm_shuffle_epi32(x, 0x4E));
  x = simdMax(x, _mm_shuffle_epi32(x, 0xB1));
  return _mm_cvtsi128_si32(x);
}

static inline int simdHorizontalMin1(SimdInt x) {
  x = _mm_min_epu8(x, _mm_srli_epi16(x, 8));
  x = _mm_minpos_epu16(x);  // SSE4.1
  return _mm_extract_epi16(x, 0);
}

static inline double simdHorizontalAddDbl(SimdDbl x) {
  return _mm_cvtsd_f64(_mm_hadd_pd(x, x));
}

static inline SimdInt simdChoose1(SimdInt items, SimdInt choices) {
  return _mm_shuffle_epi8(items, choices);  // SSSE3
}

#elif defined __ARM_NEON

typedef int32x4_t SimdInt;
typedef uint32x4_t SimdUint;
typedef uint8x16_t SimdUint1;
typedef float64x2_t SimdDbl;

const int simdBytes = 16;

static inline SimdInt simdZero() {
  return vdupq_n_s32(0);
}

static inline SimdUint1 simdZero1() {
  return vdupq_n_u8(0);
}

static inline SimdDbl simdZeroDbl() {
  return vdupq_n_f64(0);
}

static inline SimdUint1 simdOnes1() {
  return vdupq_n_u8(-1);
}

static inline SimdInt simdLoad(const int *p) {
  return vld1q_s32(p);
}

static inline SimdUint1 simdLoad1(const unsigned char *p) {
  return vld1q_u8(p);
}

static inline SimdDbl simdLoadDbl(const double *p) {
  return vld1q_f64(p);
}

static inline void simdStore(int *p, SimdInt x) {
  vst1q_s32(p, x);
}

static inline void simdStore1(unsigned char *p, SimdUint1 x) {
  vst1q_u8(p, x);
}

static inline void simdStoreDbl(double *p, SimdDbl x) {
  vst1q_f64(p, x);
}

static inline SimdUint1 simdOr1(SimdUint1 x, SimdUint1 y) {
  return vorrq_u8(x, y);
}

static inline SimdInt simdBlend(SimdInt x, SimdInt y, SimdUint mask) {
  return vbslq_s32(mask, y, x);
}

const int simdLen = 4;
const int simdDblLen = 2;

static inline SimdInt simdSet(unsigned i3, unsigned i2,
                              unsigned i1, unsigned i0) {
  size_t lo = i1;
  size_t hi = i3;
  return
    vcombine_s32(vcreate_s32((lo << 32) | i0), vcreate_s32((hi << 32) | i2));
}

static inline SimdUint1 simdSet1(unsigned char iF, unsigned char iE,
				 unsigned char iD, unsigned char iC,
				 unsigned char iB, unsigned char iA,
				 unsigned char i9, unsigned char i8,
				 unsigned char i7, unsigned char i6,
				 unsigned char i5, unsigned char i4,
				 unsigned char i3, unsigned char i2,
				 unsigned char i1, unsigned char i0) {
  size_t lo =
    (size_t)i0       | (size_t)i1 <<  8 | (size_t)i2 << 16 | (size_t)i3 << 24 |
    (size_t)i4 << 32 | (size_t)i5 << 40 | (size_t)i6 << 48 | (size_t)i7 << 56;

  size_t hi =
    (size_t)i8       | (size_t)i9 <<  8 | (size_t)iA << 16 | (size_t)iB << 24 |
    (size_t)iC << 32 | (size_t)iD << 40 | (size_t)iE << 48 | (size_t)iF << 56;

  return vcombine_u8(vcreate_u8(lo), vcreate_u8(hi));
}

static inline SimdDbl simdSetDbl(double i1, double i0) {
  return vcombine_f64(vdup_n_f64(i0), vdup_n_f64(i1));
}

static inline SimdInt simdFill(int x) {
  return vdupq_n_s32(x);
}

static inline SimdUint1 simdFill1(unsigned char x) {
  return vdupq_n_u8(x);
}

static inline SimdDbl simdFillDbl(double x) {
  return vdupq_n_f64(x);
}

static inline SimdUint simdGt(SimdInt x, SimdInt y) {
  return vcgtq_s32(x, y);
}

static inline SimdUint1 simdGe1(SimdUint1 x, SimdUint1 y) {
  return vcgeq_u8(x, y);
}

static inline SimdInt simdAdd(SimdInt x, SimdInt y) {
  return vaddq_s32(x, y);
}

static inline SimdUint1 simdAdd1(SimdUint1 x, SimdUint1 y) {
  return vaddq_u8(x, y);
}

static inline SimdUint1 simdAdds1(SimdUint1 x, SimdUint1 y) {
  return vqaddq_u8(x, y);
}

static inline SimdDbl simdAddDbl(SimdDbl x, SimdDbl y) {
  return vaddq_f64(x, y);
}

static inline SimdInt simdSub(SimdInt x, SimdInt y) {
  return vsubq_s32(x, y);
}

static inline SimdUint1 simdSub1(SimdUint1 x, SimdUint1 y) {
  return vsubq_u8(x, y);
}

static inline SimdDbl simdMulDbl(SimdDbl x, SimdDbl y) {
  return vmulq_f64(x, y);
}

static inline SimdUint1 simdQuadruple1(SimdUint1 x) {
  return vshlq_n_u8(x, 2);
}

static inline SimdInt simdMax(SimdInt x, SimdInt y) {
  return vmaxq_s32(x, y);
}

static inline SimdUint1 simdMin1(SimdUint1 x, SimdUint1 y) {
  return vminq_u8(x, y);
}

static inline int simdHorizontalMax(SimdInt x) {
  return vmaxvq_s32(x);
}

static inline int simdHorizontalMin1(SimdUint1 x) {
  return vminvq_u8(x);
}

static inline double simdHorizontalAddDbl(SimdDbl x) {
  return vaddvq_f64(x);
}

static inline SimdUint1 simdChoose1(SimdUint1 items, SimdUint1 choices) {
  return vqtbl1q_u8(items, choices);
}

#else

typedef int SimdInt;
typedef double SimdDbl;
const int simdBytes = 1;
const int simdLen = 1;
const int simdDblLen = 1;
static inline int simdZero() { return 0; }
static inline double simdZeroDbl() { return 0; }
static inline int simdSet(int x) { return x; }
static inline double simdSetDbl(double x) { return x; }
static inline int simdFill(int x) { return x; }
static inline int simdLoad(const int *p) { return *p; }
static inline double simdLoadDbl(const double *p) { return *p; }
static inline void simdStore(int *p, int x) { *p = x; }
static inline void simdStoreDbl(double *p, double x) { *p = x; }
static inline double simdFillDbl(double x) { return x; }
static inline int simdGt(int x, int y) { return x > y; }
static inline int simdAdd(int x, int y) { return x + y; }
static inline double simdAddDbl(double x, double y) { return x + y; }
static inline int simdSub(int x, int y) { return x - y; }
static inline double simdMulDbl(double x, double y) { return x * y; }
static inline int simdMax(int x, int y) { return x > y ? x : y; }
static inline int simdBlend(int x, int y, int mask) { return mask ? y : x; }
static inline int simdHorizontalMax(int a) { return a; }
static inline double simdHorizontalAddDbl(double x) { return x; }

#endif

}

#endif
