// SIMD helper
// optimze based on technolegy double, float and integer (32) SIMD instructions
// writen by Martin Steinegger

// #ifdef SIMD_H
// #error "simd.h included multiple times in the same translation unit"
// #endif

#ifndef SIMD_H
#define SIMD_H
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <iostream>
#include <float.h>

#define AVX512_ALIGN_DOUBLE		64
#define AVX512_VECSIZE_DOUBLE	8
#define AVX512_ALIGN_FLOAT		64
#define AVX512_VECSIZE_FLOAT	16
#define AVX512_ALIGN_INT		64
#define AVX512_VECSIZE_INT		16

#define AVX_ALIGN_DOUBLE		32
#define AVX_VECSIZE_DOUBLE		4
#define AVX_ALIGN_FLOAT			32
#define AVX_VECSIZE_FLOAT		8
#define AVX2_ALIGN_INT			32
#define AVX2_VECSIZE_INT		8

#define SSE_ALIGN_DOUBLE		16
#define SSE_VECSIZE_DOUBLE		2
#define SSE_ALIGN_FLOAT			16
#define SSE_VECSIZE_FLOAT		4
#define SSE_ALIGN_INT			16
#define SSE_VECSIZE_INT			4

#define MAX_ALIGN_DOUBLE	AVX512_ALIGN_DOUBLE
#define MAX_VECSIZE_DOUBLE	AVX512_VECSIZE_DOUBLE
#define MAX_ALIGN_FLOAT		AVX512_ALIGN_FLOAT
#define MAX_VECSIZE_FLOAT	AVX512_VECSIZE_FLOAT
#define MAX_ALIGN_INT		AVX512_ALIGN_INT
#define MAX_VECSIZE_INT		AVX512_VECSIZE_INT

#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/simde-features.h>

#if defined(SIMDE_X86_AVX512F_NATIVE) && defined(SIMDE_X86_AVX512BW_NATIVE)
#define AVX512
#endif

#if defined(SIMDE_X86_AVX2_NATIVE) && !defined(AVX512)
#define AVX2
#endif

#ifdef AVX512
// FIXME: Remove after updating SIMDe, headers are buggy in this versions
#define SIMDE_X86_AVX512_CMPLE_H
#define SIMDE_X86_AVX512_CMPGE_H
#include <simde/x86/avx512.h>

inline float simdf32_hmax_avx512(const __m512 buffer) {
    const __m512 shuffle1 = _mm512_shuffle_ps(buffer, buffer,(_MM_PERM_ENUM)_MM_SHUFFLE(1, 0, 3, 2));
    const __m512 max1 = _mm512_max_ps(buffer, shuffle1);
    const __m512 shuffle2 = _mm512_shuffle_ps(max1, max1, (_MM_PERM_ENUM)_MM_SHUFFLE(2, 3, 0, 1));
    const __m512 max2 = _mm512_max_ps(max1, shuffle2);
    const __m512 shuffle3 = _mm512_shuffle_f32x4(max2, max2, _MM_SHUFFLE(1, 0, 3, 2));
    const __m512 max3 = _mm512_max_ps(max2, shuffle3);
    const __m512 shuffle4 = _mm512_shuffle_f32x4(max3, max3, _MM_SHUFFLE(2, 3, 0, 1));
    const __m512 max4 = _mm512_max_ps(max3, shuffle4);
    // const __m128 max128 = _mm512_castps512_ps128(max4);
    return _mm512_cvtss_f32(max4);
}

inline uint32_t simdi32_hmax_avx512(const __m512i buffer) {
    const __m512i shuffle1 = _mm512_shuffle_epi32(buffer, (_MM_PERM_ENUM)_MM_SHUFFLE(1, 0, 3, 2));
    const __m512i max1 = _mm512_max_epi32(buffer, shuffle1);
    const __m512i shuffle2 = _mm512_shuffle_epi32(max1, (_MM_PERM_ENUM)_MM_SHUFFLE(2, 3, 0, 1));
    const __m512i max2 = _mm512_max_epi32(max1, shuffle2);
    const __m512i shuffle3 = _mm512_shuffle_i32x4(max2, max2, _MM_SHUFFLE(1, 0, 3, 2));
    const __m512i max3 = _mm512_max_epi32(max2, shuffle3);
    const __m512i shuffle4 = _mm512_shuffle_i32x4(max3, max3, _MM_SHUFFLE(2, 3, 0, 1));
    const __m512i max4 = _mm512_max_epi32(max3, shuffle4);
    return (uint32_t)_mm512_cvtsi512_si32(max4);
}

uint16_t simd_hmax16_sse(const __m128i buffer);

inline uint16_t simdi16_hmax_avx512(const __m512i buffer) {
    const __m128i abcd = _mm512_castsi512_si128(buffer);
    const uint16_t first = simd_hmax16_sse(abcd);
    const __m128i efgh = _mm512_extracti32x4_epi32(buffer, 1);
    const uint16_t second = simd_hmax16_sse(efgh);
    const uint16_t lower = std::max(first, second);

    const __m128i hijk = _mm512_extracti32x4_epi32(buffer, 2);
    const uint16_t third = simd_hmax16_sse(hijk);
    const __m128i lmno = _mm512_extracti32x4_epi32(buffer, 3);
    const uint16_t forth = simd_hmax16_sse(lmno);
    const uint16_t upper = std::max(third, forth);

    return std::max(lower, upper);
}

inline uint8_t simdi8_hmax_avx512(const __m512i buffer) {
    // https://github.com/EddyRivasLab/easel/blob/07ca83ba9ef0414dba9ce0a9331d465b5eb58f2b/esl_avx512.h#L35-L63
    // Use AVX instructions for this because AVX-512 can't extract 8-bit quantities
    // Intel has stated that there will be no performance penalty for switching between AVX-512 and AVX
    __m256i b = _mm256_max_epu8(_mm512_extracti64x4_epi64(buffer, 0), _mm512_extracti64x4_epi64(buffer, 1)); //changed from extract i32x8
    b = _mm256_max_epu8(b, _mm256_permute2x128_si256(b, b, 0x01));
    b = _mm256_max_epu8(b, _mm256_shuffle_epi32     (b,    0x4e));
    b = _mm256_max_epu8(b, _mm256_shuffle_epi32     (b,    0xb1));
    b = _mm256_max_epu8(b, _mm256_shufflelo_epi16   (b,    0xb1));
    b = _mm256_max_epu8(b, _mm256_srli_si256        (b,    1));
    return _mm256_extract_epi8(b, 0);  // epi8 is fine here. gets cast properly to uint8_t on return.
}


inline float simdf32_hadd(const __m512 buffer) {
    const __m512 shuffle1 = _mm512_shuffle_ps(buffer,buffer, (_MM_PERM_ENUM)_MM_SHUFFLE(1, 0, 3, 2));
    const __m512 max1 = _mm512_add_ps(buffer, shuffle1);
    const __m512 shuffle2 = _mm512_shuffle_ps(max1, max1, (_MM_PERM_ENUM)_MM_SHUFFLE(2, 3, 0, 1));
    const __m512 max2 = _mm512_add_ps(max1, shuffle2);
    const __m512 shuffle3 = _mm512_shuffle_f32x4(max2, max2, _MM_SHUFFLE(1, 0, 3, 2));
    const __m512 max3 = _mm512_add_ps(max2, shuffle3);
    const __m512 shuffle4 = _mm512_shuffle_f32x4(max3, max3, _MM_SHUFFLE(2, 3, 0, 1));
    const __m512 max4 = _mm512_add_ps(max3, shuffle4);
    return _mm512_cvtss_f32(max4);
}

template  <unsigned int N>
inline __m512i _mm512_shift_left(__m512i a) {
    //Only work for N <= 16
    __m512i mask = _mm512_permutex2var_epi64(_mm512_setzero_si512(), _mm512_set_epi64(13, 12, 11, 10, 9, 8, 7, 6), a);
    return _mm512_alignr_epi8(a,mask,16-N);
}

inline __m512i _mm512_shift_left4(__m512i a) {
    //Only work for N == 4
    return _mm512_alignr_epi32(a, _mm512_setzero_si512(), (64 - 4)/ 4);
}

inline __m512i _mm512_shift_left8(__m512i a) {
    //Only work for N == 8
    return _mm512_alignr_epi32(a, _mm512_setzero_si512(), (64 - 8)/ 4);
}

inline __m512i _mm512_shift_left16(__m512i a) {
    //Only work for N == 16
    return _mm512_alignr_epi32(a, _mm512_setzero_si512(), (64 - 16)/ 4);
}

inline __m512i _mm512_shift_left32(__m512i a) {
    //Only work for N == 32
    return _mm512_alignr_epi32(a, _mm512_setzero_si512(), (64 - 32)/ 4);
}

// Materialize a compare mask (k-register) into an all-ones/all-zeros vector so the
// result matches the AVX2/SSE contract (compares return vectors that callers and/or/blendv).
// _mm512_maskz_set1_* is the idiomatic single-instruction spelling: it lowers to one
// VPMOVM2{B,W,D,Q} on any AVX512(+DQ) compiler, whereas the older
// _mm512_mask_mov(setzero, mask, set1(-1)) idiom relies on the compiler recognizing a
// 3-input pattern (GCC does not always fold it, emitting an extra broadcast + merge-move).
inline __m512i simdi32_gt_avx512(__m512i a, __m512i b) {
    __mmask16 mask = _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_NLE);
    return _mm512_maskz_set1_epi32(mask, -1);
}

inline __m512i simdi16_gt_avx512(__m512i a, __m512i b) {
    __mmask32 mask = _mm512_cmp_epi16_mask(a, b, _MM_CMPINT_NLE);
    return _mm512_maskz_set1_epi16(mask, -1);
}

inline __m512i simdi8_gt_avx512(__m512i a, __m512i b) {
    __mmask64 mask = _mm512_cmp_epi8_mask(a, b, _MM_CMPINT_NLE);
    return _mm512_maskz_set1_epi8(mask, -1);
}

inline __m512i simdi32_eq_avx512(__m512i a, __m512i b) {
    __mmask16 mask = _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ);
    return _mm512_maskz_set1_epi32(mask, -1);
}

inline __m512i simdi16_eq_avx512(__m512i a, __m512i b) {
    __mmask32 mask = _mm512_cmp_epi16_mask(a, b, _MM_CMPINT_EQ);
    return _mm512_maskz_set1_epi16(mask, -1);
}

inline __m512i simdi8_eq_avx512(__m512i a, __m512i b) {
    __mmask64 mask = _mm512_cmp_epi8_mask(a, b, _MM_CMPINT_EQ);
    return _mm512_maskz_set1_epi8(mask, -1);
}

inline __m512 simdf32_gt_avx512(__m512 a, __m512 b) {
    __mmask16 mask = _mm512_cmp_ps_mask(a, b, _CMP_GT_OS);
    return _mm512_castsi512_ps(_mm512_maskz_set1_epi32(mask, -1));
}

inline __m512d simdf64_gt_avx512(__m512d a, __m512d b) {
    __mmask8 mask = _mm512_cmp_pd_mask(a, b, _CMP_GT_OS);
    return _mm512_castsi512_pd(_mm512_maskz_set1_epi64(mask, -1));
}

inline __m512d simdf64_lt_avx512(__m512d a, __m512d b) {
    __mmask8 mask = _mm512_cmp_pd_mask(a, b, _CMP_LT_OS);
    return _mm512_castsi512_pd(_mm512_maskz_set1_epi64(mask, -1));
}

inline __m512 simdf32_eq_avx512(__m512 a, __m512 b) {
    __mmask16 mask = _mm512_cmp_ps_mask(a, b, _CMP_EQ_OS);
    return _mm512_castsi512_ps(_mm512_maskz_set1_epi32(mask, -1));
}

inline __m512 simdf32_le_avx512(__m512 a, __m512 b) {
    __mmask16 mask = _mm512_cmp_ps_mask(a, b, _CMP_LE_OQ);
    return _mm512_castsi512_ps(_mm512_maskz_set1_epi32(mask, -1));
}

// AI Idea to get around z in _mm512_mask_blend_ps(z,x,y) beeing __mmask16 
inline __mmask16 simdf32_mask_from_ps_avx512(__m512 m) {
    return _mm512_movepi32_mask(_mm512_castps_si512(m));
}

inline bool simd_any_avx512(const __m512i buffer) {
    const __m512i vZero = _mm512_setzero_si512();
    const uint64_t mask = (uint64_t) _mm512_cmpeq_epi8_mask(buffer, vZero);
    return (mask != 0xFFFFFFFFFFFFFFFFULL);
}

inline bool simd_eq_all_avx512(const __m512i a, const __m512i b) {
    const uint64_t mask = (uint64_t)_mm512_cmpeq_epi8_mask(a, b);
    return (mask == 0xFFFFFFFFFFFFFFFFULL);
}

inline __m512 simdf32_reverse_sse(const __m512 buffer) {
    return _mm512_shuffle_ps(buffer, buffer, _MM_SHUFFLE(0, 1, 2, 3));
}


// float simdf32_hmax_sse(const __m128 buffer);
// inline float simdf32_hmax_avx512(const __m512 buffer) {
//     const __m128 abcd = _mm256_extractf128_ps(buffer, 3); // highest 128 bits
//     const __m128 efgh = _mm256_extractf128_ps(buffer, 2); // second highest 128 bits
//     const __m128 ijkl = _mm256_extractf128_ps(buffer, 1); // second lowest 128 bits
//     const __m128 mnop = _mm256_extractf128_ps(buffer, 0); // lowest 128 bits
//     const float abcd_max = simdf32_hmax_sse(abcd);
//     const float efgh_max = simdf32_hmax_sse(efgh);
//     const float ijkl_max = simdf32_hmax_sse(ijkl);
//     const float mnop_max = simdf32_hmax_sse(mnop);
//     return std::max(abcd_max, efgh_max, ijkl_max, mnop_max);
// }

// double support
#ifndef SIMD_DOUBLE
#define SIMD_DOUBLE
#define ALIGN_DOUBLE        AVX512_ALIGN_DOUBLE
#define VECSIZE_DOUBLE      AVX512_VECSIZE_DOUBLE
typedef __m512d simd_double;
#define simdf64_add(x,y)    _mm512_add_pd(x,y)
#define simdf64_sub(x,y)    _mm512_sub_pd(x,y)
#define simdf64_mul(x,y)    _mm512_mul_pd(x,y)
#define simdf64_div(x,y)    _mm512_div_pd(x,y)
#define simdf64_max(x,y)    _mm512_max_pd(x,y)
#define simdf64_load(x)     _mm512_load_pd(x)
#define simdf64_store(x,y)  _mm512_store_pd(x,y)
#define simdf64_set(x)      _mm512_set1_pd(x)
#define simdf64_setzero(x)  _mm512_setzero_pd()
#define simdf64_gt(x,y)     simdf64_gt_avx512(x,y)
#define simdf64_lt(x,y)     simdf64_lt_avx512(x,y)
#define simdf64_or(x,y)     _mm512_or_pd(x,y)
#define simdf64_and(x,y)    _mm512_and_pd (x,y)
#define simdf64_andnot(x,y) _mm512_andnot_pd(x,y)
#define simdf64_xor(x,y)    _mm512_xor_pd(x,y)
#endif //SIMD_DOUBLE
// float support
#ifndef SIMD_FLOAT
#define SIMD_FLOAT
#define ALIGN_FLOAT         AVX512_ALIGN_FLOAT
#define VECSIZE_FLOAT       AVX512_VECSIZE_FLOAT
typedef __m512  simd_float;
#define simdf32_add(x,y)    _mm512_add_ps(x,y)
#define simdf32_sub(x,y)    _mm512_sub_ps(x,y)
#define simdf32_mul(x,y)    _mm512_mul_ps(x,y)
#define simdf32_div(x,y)    _mm512_div_ps(x,y)
#define simdf32_sqrt(x)     _mm512_sqrt_ps(x)
#define simdf32_rcp(x)      _mm512_rcp14_ps(x)
#define simdf32_max(x,y)    _mm512_max_ps(x,y)
#define simdf32_min(x,y)    _mm512_min_ps(x,y)
#define simdf32_load(x)     _mm512_load_ps(x)
#define simdf32_loadu(x)    _mm512_loadu_ps(x)
#define simdf32_store(x,y)  _mm512_store_ps(x,y)
#define simdf32_storeu(x,y) _mm512_storeu_ps(x,y)
#define simdf32_set(x)      _mm512_set1_ps(x)
#define simdf32_setzero(x)  _mm512_setzero_ps()
#define simdf32_gt(x,y)     simdf32_gt_avx512(x,y)
#define simdf32_eq(x,y)     simdf32_eq_avx512(x,y)
#define simdf32_lt(x,y)     simdf32_gt_avx512(y,x)
#define simdf32_le(x,y)     simdf32_le_avx512(x,y)
#define simdf32_cmp(x,y,z)  _mm512_castsi512_ps(_mm512_maskz_set1_epi32(_mm512_cmp_ps_mask(x, y, z), -1))
#define simdf32_or(x,y)     _mm512_or_ps(x,y)
#define simdf32_and(x,y)    _mm512_and_ps(x,y)
#define simdf32_andnot(x,y) _mm512_andnot_ps(x,y)
#define simdf32_xor(x,y)    _mm512_xor_ps(x,y)
#define simdi32_i2f(x) 	    _mm512_cvtepi32_ps(x)  // convert integer to s.p. float
#define simdi_i2fcast(x)    _mm512_castsi512_ps(x)
#define simdf32_round(x)    _mm512_roundscale_ps(x, SIMDE_MM_FROUND_TO_NEAREST_INT | SIMDE_MM_FROUND_NO_EXC)
#define simdf32_blendv_ps(x,y,z) _mm512_mask_blend_ps(simdf32_mask_from_ps_avx512(z), x, y) // AI Idea, prob broken
#define simdf32_reverse(x)  _mm512_permutexvar_ps(_mm512_setr_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0), x) 
#define simdf32_fmadd(x,y,z) _mm512_fmadd_ps(x,y,z)
#define simdf32_hmax(x)     simdf32_hmax_avx512(x)
#define simdf32_f2i(x) 	    _mm512_cvtps_epi32(x)  // convert s.p. float to integer
#define simdf_f2icast(x)    _mm512_castps_si512(x)
#endif //SIMD_FLOAT

// integer support 
#ifndef SIMD_INT
#define SIMD_INT
#define ALIGN_INT           AVX512_ALIGN_INT
#define VECSIZE_INT         AVX512_VECSIZE_INT
typedef __m512i simd_int;
#define simdi32_add(x,y)    _mm512_add_epi32(x,y)
#define simdi16_add(x,y)    _mm512_add_epi16(x,y)
#define simdi16_adds(x,y)   _mm512_adds_epi16(x,y)
#define simdui8_adds(x,y)   _mm512_adds_epu8(x,y)
#define simdi32_sub(x,y)    _mm512_sub_epi32(x,y)
#define simdi16_sub(x,y)    _mm512_sub_epi16(x,y)
#define simdui8_subs(x,y)   _mm512_subs_epu8(x,y)
#define simdui16_subs(x,y)  _mm512_subs_epu16(x,y)
#define simdui32_subs(x,y)  _mm512_max_epi32(_mm512_sub_epi32(x,y), _mm512_setzero_si512())
#define simdi32_mul(x,y)    _mm512_mullo_epi32(x,y)
#define simdui8_max(x,y)    _mm512_max_epu8(x,y)

#define simdi32_max(x,y)    _mm512_max_epi32(x, y)
#define simdi32_hmax(x)     simdi32_hmax_avx512(x)
#define simdi16_max(x,y)    _mm512_max_epi16(x, y)
#define simdi16_hmax(x)     simdi16_hmax_avx512(x)
#define simdi8_max(x,y)     _mm512_max_epi8(x, y)
#define simdi8_hmax(x)      simdi8_hmax_avx512(x)

#define simdi16_min(x,y)    _mm512_min_epi16(x,y)
#define simdui8_avg(x,y)    _mm512_avg_epu8(x,y)
#define simdui16_avg(x,y)   _mm512_avg_epu16(x,y)
#define simdi_load(x)       _mm512_load_si512(x)
#define simdi_loadu(x)      _mm512_loadu_si512(x)
#define simdi_streamload(x) _mm512_stream_load_si512(x)
#define simdi_store(x,y)    _mm512_store_si512(x,y)
#define simdi_storeu(x,y)   _mm512_storeu_si512(x,y)
#define simdi32_set(x)      _mm512_set1_epi32(x)
#define simdi16_set(x)      _mm512_set1_epi16(x)
#define simdi8_set(x)       _mm512_set1_epi8(x)
#define simdi32_shuffle(x,y) _mm512_shuffle_epi32(x,y)
#define simdi16_shuffle(x,y) NOT_YET_IMP(x,y) //_mm512_shuffle_epi16(x,y) does not exist
#define simdi8_shuffle(x,y) _mm512_shuffle_epi8(x,y)
#define simdi_setzero()     _mm512_setzero_si512()
#define simdi32_gt(x,y)     simdi32_gt_avx512(x, y)
#define simdi16_gt(x,y)     simdi16_gt_avx512(x, y)
#define simdi8_gt(x,y)      simdi8_gt_avx512(x, y)
#define simdi32_eq(x,y)     simdi32_eq_avx512(x, y)
#define simdi16_eq(x,y)     simdi16_eq_avx512(x, y)
#define simdi8_eq(x,y)      simdi8_eq_avx512(x, y)
#define simdi32_lt(x,y)     simdi32_gt_avx512(y, x)
#define simdi16_lt(x,y)     simdi16_gt_avx512(y, x)
#define simdi8_lt(x,y)      simdi8_gt_avx512(y, x)

#define SIMD_MOVEMASK_MAX   0xffffffffffffffff  // not sure if correct
typedef uint64_t            movemask_max_t;
#define popcount(x)         __builtin_popcountll(x)
#define simd_any(x)         simd_any_avx512(x)
#define simd_eq_all(x,y)    simd_eq_all_avx512(x,y)
#define simdi_or(x,y)       _mm512_or_si512(x,y)
#define simdi_and(x,y)      _mm512_and_si512(x,y)
#define simdi_andnot(x,y)   _mm512_andnot_si512(x,y)
#define simdi_xor(x,y)      _mm512_xor_si512(x,y)
#define simdi8_shiftl(x,y)  _mm512_shift_left<y>(x)
//#define simdi8_shiftr(x,y)  simdi8_shift_right<y>(x)
#define simdi8_movemask(x)  _mm512_movepi8_mask(x)
#define simdi16_extract(x,y) NOT_YET_IMP()         // no 16 bit version available
#define simdi16_slli(x,y)	_mm512_slli_epi16(x,y) // shift integers in a left by y
#define simdi16_srli(x,y)	_mm512_srli_epi16(x,y) // shift integers in a right by y
#define simdi32_slli(x,y)	_mm512_slli_epi32(x,y) // shift integers in a left by y
#define simdi32_srli(x,y)	_mm512_srli_epi32(x,y) // shift integers in a right by y
#define simdi32_srai(x,y)	_mm512_srai_epi32(x,y)
#define simdi16_pack(x,y)   _mm512_packs_epi16(x,y)
#define simdi32_pack(x,y)   _mm512_packs_epi32(x,y)
#define simdi8_blend(x,y,z)   _mm512_mask_blend_epi8(_mm512_movepi8_mask(z),x,y)
#define simdi32_insert(x,y,z) _mm512_mask_set1_epi32(x, 1UL<<z, y) // https://stackoverflow.com/questions/58303958/how-to-implement-16-and-32-bit-integer-insert-and-extract-operations-with-avx-51

#endif //SIMD_INT
#endif //AVX512_SUPPORT


#ifdef AVX2
#include <simde/x86/avx2.h>
#include <simde/x86/fma.h>
// integer support  (usable with AVX2)
#ifndef SIMD_INT
#define SIMD_INT
#define ALIGN_INT           AVX2_ALIGN_INT
#define VECSIZE_INT         AVX2_VECSIZE_INT
uint32_t simd_hmax32_sse(const __m128i buffer);
uint16_t simd_hmax16_sse(const __m128i buffer);
uint8_t simd_hmax8_sse(const __m128i buffer);
bool simd_any_sse(const __m128i buffer);

inline uint32_t simd_hmax32_avx(const __m256i buffer) {
    const __m128i abcd = _mm256_castsi256_si128(buffer);
    const uint32_t first = simd_hmax32_sse(abcd);
    const __m128i efgh = _mm256_extracti128_si256(buffer, 1);
    const uint32_t second = simd_hmax32_sse(efgh);
    return std::max(first, second);
}

inline uint16_t simd_hmax16_avx(const __m256i buffer) {
    const __m128i abcd = _mm256_castsi256_si128(buffer);
    const uint16_t first = simd_hmax16_sse(abcd);
    const __m128i efgh = _mm256_extracti128_si256(buffer, 1);
    const uint16_t second = simd_hmax16_sse(efgh);
    return std::max(first, second);
}

inline uint8_t simd_hmax8_avx(const __m256i buffer) {
    const __m128i abcd = _mm256_castsi256_si128(buffer);
    const uint8_t first = simd_hmax8_sse(abcd);
    const __m128i efgh = _mm256_extracti128_si256(buffer, 1);
    const uint8_t second = simd_hmax8_sse(efgh);
    return std::max(first, second);
}

inline bool simd_any_avx(const __m256i buffer) {
#if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
    const __m128i lo = _mm256_castsi256_si128(buffer);
    const __m128i hi = _mm256_extracti128_si256(buffer, 1);
    const __m128i lo_or_hi = _mm_or_si128(lo, hi);
    return simd_any_sse(lo_or_hi);
#else
    const __m256i vZero = _mm256_set1_epi32(0);
    const __m256i vTemp = _mm256_cmpeq_epi8(buffer, vZero);
    const uint32_t mask = _mm256_movemask_epi8(vTemp);
    return (mask != 0xffffffff);
#endif
}

inline bool simd_eq_all_avx(const __m256i a, const __m256i b) {
    const __m256i vector_mask = _mm256_cmpeq_epi8(a, b);
#if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
    const __m128i lo = _mm256_castsi256_si128(vector_mask);
    const __m128i hi = _mm256_extracti128_si256(vector_mask, 1);
    const __m128i lo_and_hi = _mm_and_si128(lo, hi);
    const uint32_t min_dword = vminvq_u32(vreinterpretq_u32_s64(lo_and_hi));
    return (min_dword == 0xffffffff);
#else
    const uint32_t bit_mask = _mm256_movemask_epi8(vector_mask);
    return (bit_mask == 0xffffffff);
#endif
}

float simdf32_hmax_sse(const __m128 buffer);

inline float simdf32_hmax_avx(const __m256 buffer) {
    const __m128 lower = _mm256_castps256_ps128(buffer); // Lower 128 bits
    const __m128 upper = _mm256_extractf128_ps(buffer, 1); // Upper 128 bits
    const float lower_max = simdf32_hmax_sse(lower);
    const float upper_max = simdf32_hmax_sse(upper);
    return std::max(lower_max, upper_max);
}
__m128 simdf32_reverse_sse(const __m128 buffer);

inline __m256 simdf32_reverse_avx256(const __m256 buffer) {
    const __m128 lower = _mm256_castps256_ps128(buffer); // Lower 128 bits
    const __m128 upper = _mm256_extractf128_ps(buffer, 1); // Upper 128 bits
    const __m128 lower_rev = simdf32_reverse_sse(lower);
    const __m128 upper_rev = simdf32_reverse_sse(upper);
    return _mm256_insertf128_ps(_mm256_castps128_ps256(upper_rev), lower_rev, 0x1);
}


inline float simdf32_hadd(const __m256 x) {
    /* ( x3+x7, x2+x6, x1+x5, x0+x4 ) */
    const __m128 x128 = _mm_add_ps(_mm256_extractf128_ps(x, 1), _mm256_castps256_ps128(x));
    /* ( -, -, x1+x3+x5+x7, x0+x2+x4+x6 ) */
    const __m128 x64 = _mm_add_ps(x128, _mm_movehl_ps(x128, x128));
    /* ( -, -, -, x0+x1+x2+x3+x4+x5+x6+x7 ) */
    const __m128 x32 = _mm_add_ss(x64, _mm_shuffle_ps(x64, x64, 0x55));
    /* Conversion to float is a no-op on x86-64 */
    return _mm_cvtss_f32(x32);
}

template  <unsigned int N>
inline __m256i _mm256_shift_left(__m256i a) {
    __m256i mask = _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0) );
    return _mm256_alignr_epi8(a,mask,16-N);
}

inline unsigned short extract_epi16(__m256i v, int pos) {
    switch(pos){
        case 0: return _mm256_extract_epi16(v, 0);
        case 1: return _mm256_extract_epi16(v, 1);
        case 2: return _mm256_extract_epi16(v, 2);
        case 3: return _mm256_extract_epi16(v, 3);
        case 4: return _mm256_extract_epi16(v, 4);
        case 5: return _mm256_extract_epi16(v, 5);
        case 6: return _mm256_extract_epi16(v, 6);
        case 7: return _mm256_extract_epi16(v, 7);
        case 8: return _mm256_extract_epi16(v, 8);
        case 9: return _mm256_extract_epi16(v, 9);
        case 10: return _mm256_extract_epi16(v, 10);
        case 11: return _mm256_extract_epi16(v, 11);
        case 12: return _mm256_extract_epi16(v, 12);
        case 13: return _mm256_extract_epi16(v, 13);
        case 14: return _mm256_extract_epi16(v, 14);
        case 15: return _mm256_extract_epi16(v, 15);
    }
    return 0;
}

typedef __m256i simd_int;
#define simdi32_add(x,y)    _mm256_add_epi32(x,y)
#define simdi16_add(x,y)    _mm256_add_epi16(x,y)
#define simdi16_adds(x,y)   _mm256_adds_epi16(x,y)
#define simdi16_sub(x,y)    _mm256_sub_epi16(x,y)
#define simdui8_adds(x,y)   _mm256_adds_epu8(x,y)
#define simdi32_sub(x,y)    _mm256_sub_epi32(x,y)
#define simdui32_subs(x,y)  _mm256_max_epi32(_mm256_sub_epi32(x,y), _mm256_setzero_si256())
#define simdui16_subs(x,y)  _mm256_subs_epu16(x,y)
#define simdui8_subs(x,y)   _mm256_subs_epu8(x,y)
#define simdi32_mul(x,y)    _mm256_mullo_epi32(x,y)
#define simdi32_max(x,y)    _mm256_max_epi32(x,y)
#define simdi16_max(x,y)    _mm256_max_epi16(x,y)
#define simdi16_min(x,y)    _mm256_min_epi16(x,y)
#define simdi32_insert(x,y,z) _mm256_insert_epi32(x,y,z)
#define simdi32_extract(x,y) _mm256_extract_epi32(x,y)
#define simdi32_hmax(x)     simd_hmax32_avx(x)
#define simdi16_hmax(x)     simd_hmax16_avx(x)
#define simdui8_max(x,y)    _mm256_max_epu8(x,y)
#define simdi8_hmax(x)      simd_hmax8_avx(x)
#define simd_any(x)         simd_any_avx(x)
#define simd_eq_all(x,y)    simd_eq_all_avx(x,y)
#define simdi_load(x)       _mm256_load_si256(x)
#define simdf_load(x)       _mm256_load_ps(x)
#define simdui8_avg(x,y)    _mm256_avg_epu8(x,y)
#define simdui16_avg(x,y)   _mm256_avg_epu16(x,y)
#define simdi_loadu(x)       _mm256_loadu_si256(x)
#define simdi_streamload(x) _mm256_stream_load_si256(x)
#define simdi_store(x,y)    _mm256_store_si256(x,y)
#define simdi_storeu(x,y)   _mm256_storeu_si256(x,y)
#define simdi32_set(x)      _mm256_set1_epi32(x)
#define simdi16_set(x)      _mm256_set1_epi16(x)
#define simdi8_set(x)       _mm256_set1_epi8(x)
#define simdi32_shuffle(x,y) _mm256_shuffle_epi32(x,y)
#define simdi16_shuffle(x,y) NOT_YET_IMP() //_mm256_shuffle_epi16(x,y) does not exist
#define simdi8_shuffle(x,y)  _mm256_shuffle_epi8(x,y)
#define simdi_setzero()     _mm256_setzero_si256()
#define simdi8_blend(x,y,z) _mm256_blendv_epi8(x,y,z)
#define simdi32_gt(x,y)     _mm256_cmpgt_epi32(x,y)
#define simdi8_gt(x,y)      _mm256_cmpgt_epi8(x,y)
#define simdi16_gt(x,y)     _mm256_cmpgt_epi16(x,y)
#define simdi8_eq(x,y)      _mm256_cmpeq_epi8(x,y)
#define simdi16_eq(x,y)     _mm256_cmpeq_epi16(x,y)
#define simdi32_eq(x,y)     _mm256_cmpeq_epi32(x,y)
#define simdi32_lt(x,y)     _mm256_cmpgt_epi32(y,x) // inverse
#define simdi16_lt(x,y)     _mm256_cmpgt_epi16(y,x) // inverse
#define simdi8_lt(x,y)      _mm256_cmpgt_epi8(y,x)
#define simdi_or(x,y)       _mm256_or_si256(x,y)
#define simdi_and(x,y)      _mm256_and_si256(x,y)
#define simdi_andnot(x,y)   _mm256_andnot_si256(x,y)
#define simdi_xor(x,y)      _mm256_xor_si256(x,y)
#define simdi8_shiftl(x,y)  _mm256_shift_left<y>(x)
//TODO fix like shift_left
#define simdi8_shiftr(x,y)  _mm256_srli_si256(x,y)
#define SIMD_MOVEMASK_MAX   0xffffffff
typedef uint32_t            movemask_max_t;
#define popcount(x)         __builtin_popcount(x)
#define simdi8_movemask(x)  _mm256_movemask_epi8(x)
#define simdi16_extract(x,y) extract_epi16(x,y)
#define simdi32_pack(x,y)   _mm256_packs_epi32(x,y)
#define simdi16_pack(x,y)   _mm256_packs_epi16(x,y)
#define simdi16_slli(x,y)	_mm256_slli_epi16(x,y) // shift integers in a left by y
#define simdi16_srli(x,y)	_mm256_srli_epi16(x,y) // shift integers in a right by y
#define simdi32_slli(x,y)   _mm256_slli_epi32(x,y) // shift integers in a left by y
#define simdi32_srli(x,y)   _mm256_srli_epi32(x,y) // shift integers in a right by y
#define simdi32_srai(x,y)   _mm256_srai_epi32(x,y)
#define simdi32_i2f(x) 	    _mm256_cvtepi32_ps(x)  // convert integer to s.p. float
#define simdi_i2fcast(x)    _mm256_castsi256_ps(x)
#endif

#include <simde/x86/avx.h>
// double support (usable with AVX1)
#ifndef SIMD_DOUBLE
#define SIMD_DOUBLE
#define ALIGN_DOUBLE        AVX_ALIGN_DOUBLE
#define VECSIZE_DOUBLE      AVX_VECSIZE_DOUBLE
typedef __m256d simd_double;
#define simdf64_add(x,y)    _mm256_add_pd(x,y)
#define simdf64_sub(x,y)    _mm256_sub_pd(x,y)
#define simdf64_mul(x,y)    _mm256_mul_pd(x,y)
#define simdf64_div(x,y)    _mm256_div_pd(x,y)
#define simdf64_max(x,y)    _mm256_max_pd(x,y)
#define simdf64_load(x)     _mm256_load_pd(x)
#define simdf64_store(x,y)  _mm256_store_pd(x,y)
#define simdf64_set(x)      _mm256_set1_pd(x)
#define simdf64_setzero(x)  _mm256_setzero_pd()
#define simdf64_gt(x,y)     _mm256_cmp_pd(x,y,_CMP_GT_OS)
#define simdf64_lt(x,y)     _mm256_cmp_pd(x,y,_CMP_LT_OS)
#define simdf64_or(x,y)     _mm256_or_pd(x,y)
#define simdf64_and(x,y)    _mm256_and_pd(x,y)
#define simdf64_andnot(x,y) _mm256_andnot_pd(x,y)
#define simdf64_xor(x,y)    _mm256_xor_pd(x,y)
#endif //SIMD_DOUBLE
// float support (usable with AVX1)
#ifndef SIMD_FLOAT
#define SIMD_FLOAT
#define ALIGN_FLOAT         AVX_ALIGN_FLOAT
#define VECSIZE_FLOAT       AVX_VECSIZE_FLOAT
typedef __m256 simd_float;
#define simdf32_add(x,y)    _mm256_add_ps(x,y)
#define simdf32_sub(x,y)    _mm256_sub_ps(x,y)
#define simdf32_mul(x,y)    _mm256_mul_ps(x,y)
#define simdf32_div(x,y)    _mm256_div_ps(x,y)
#define simdf32_sqrt(x)     _mm256_sqrt_ps(x)
#define simdf32_rcp(x)      _mm256_rcp_ps(x)
#define simdf32_max(x,y)    _mm256_max_ps(x,y)
#define simdf32_min(x,y)    _mm256_min_ps(x,y)
#define simdf32_load(x)     _mm256_load_ps(x)
#define simdf32_store(x,y)  _mm256_store_ps(x,y)
#define simdf32_loadu(x)    _mm256_loadu_ps(x)
#define simdf32_storeu(x,y) _mm256_storeu_ps(x,y)
#define simdf32_set(x)      _mm256_set1_ps(x)
#define simdf32_setzero(x)  _mm256_setzero_ps()
#define simdf32_gt(x,y)     _mm256_cmp_ps(x,y,_CMP_GT_OS)
#define simdf32_eq(x,y)     _mm256_cmp_ps(x,y,_CMP_EQ_OS)
#define simdf32_lt(x,y)     _mm256_cmp_ps(x,y,_CMP_LT_OS)
#define simdf32_le(x,y)     _mm256_cmp_ps(x,y,_CMP_LE_OQ) // OQ or OS?
#define simdf32_cmp(x,y,z) _mm256_cmp_ps(x,y,z)
#define simdf32_or(x,y)     _mm256_or_ps(x,y)
#define simdf32_and(x,y)    _mm256_and_ps(x,y)
#define simdf32_andnot(x,y) _mm256_andnot_ps(x,y)
#define simdf32_xor(x,y)    _mm256_xor_ps(x,y)
#define simdf32_f2i(x) 	    _mm256_cvtps_epi32(x)  // convert s.p. float to integer
#define simdf_f2icast(x)    _mm256_castps_si256(x) // compile time cast
#define simdf32_round(x)    _mm256_round_ps(x, SIMDE_MM_FROUND_TO_NEAREST_INT | SIMDE_MM_FROUND_NO_EXC)
#define simdf32_blendv_ps(x,y,z) _mm256_blendv_ps(x,y,z)
#define simdf32_movemask_ps(x) _mm256_movemask_ps(x)
#define simdf32_hmax(x)     simdf32_hmax_avx(x)
#define simdf32_reverse(x)  simdf32_reverse_avx256(x)
#define simdf32_fmadd(x,y,z) _mm256_fmadd_ps(x,y,z)
#endif
#endif

#include <simde/x86/sse4.1.h>
// see https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction
inline uint32_t simd_hmax32_sse(const __m128i buffer) {
#if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
    return vmaxvq_u32(vreinterpretq_u32_s64(buffer));
#else
    __m128i hi64  = _mm_shuffle_epi32(buffer, _MM_SHUFFLE(1, 0, 3, 2));
    __m128i max64 = _mm_max_epi32(hi64, buffer);
    __m128i hi32  = _mm_shufflelo_epi16(max64, _MM_SHUFFLE(1, 0, 3, 2)); // Swap the low two elements
    __m128i max32 = _mm_max_epi32(max64, hi32);
    return _mm_cvtsi128_si32(max32); // SSE2 movd
#endif
}

inline uint16_t simd_hmax16_sse(const __m128i buffer) {
#if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
    return vmaxvq_u16(vreinterpretq_u16_s64(buffer));
#else
    __m128i tmp1 = _mm_subs_epu16(_mm_set1_epi16((short)65535), buffer);
    __m128i tmp3 = _mm_minpos_epu16(tmp1);
    return (65535 - _mm_cvtsi128_si32(tmp3));
#endif
}

inline uint8_t simd_hmax8_sse(const __m128i buffer) {
#if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
    return vmaxvq_u8(vreinterpretq_u8_s64(buffer));
#else
    __m128i tmp1 = _mm_subs_epu8(_mm_set1_epi8((char)255), buffer);
    __m128i tmp2 = _mm_min_epu8(tmp1, _mm_srli_epi16(tmp1, 8));
    __m128i tmp3 = _mm_minpos_epu16(tmp2);
    return (int8_t)(255 -(int8_t) _mm_cvtsi128_si32(tmp3));
#endif
}

inline float simdf32_hmax_sse(const __m128 buffer) {
    __m128 temp = _mm_shuffle_ps(buffer, buffer, _MM_SHUFFLE(2, 3, 0, 1)); // Swap adjacent elements
    __m128 max1 = _mm_max_ps(buffer, temp);

    temp = _mm_shuffle_ps(max1, max1, _MM_SHUFFLE(1, 0, 3, 2)); // Compare across the remaining pairs
    __m128 max2 = _mm_max_ps(max1, temp);

    return _mm_cvtss_f32(max2);
}

inline bool simd_any_sse(const __m128i buffer) {
#if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
    const uint32_t non_zero = vmaxvq_u32(vreinterpretq_u32_s64(buffer));
    return static_cast<bool>(non_zero);
#else
    const __m128i vZero = _mm_set1_epi32(0);
    const __m128i vTemp = _mm_cmpeq_epi8(buffer, vZero);
    const uint32_t mask = _mm_movemask_epi8(vTemp);
    return (mask != 0xffff);
#endif
}

inline bool simd_eq_all_sse(const __m128i a, const __m128i b) {
    const __m128i vector_mask = _mm_cmpeq_epi8(a, b);
#if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
    const uint32_t min_dword = vminvq_u32(vreinterpretq_u32_s64(vector_mask));
    return (min_dword == 0xffffffff);
#else
    const uint32_t bit_mask = _mm_movemask_epi8(vector_mask);
    return (bit_mask == 0xffff);
#endif
}

inline __m128 simdf32_reverse_sse(const __m128 buffer) {
    return _mm_shuffle_ps(buffer, buffer, _MM_SHUFFLE(0, 1, 2, 3));
}

// see https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction
inline float simdf32_hadd(const __m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = _mm_add_ss(sums, shuf);
    return        _mm_cvtss_f32(sums);
}

// double support
#ifndef SIMD_DOUBLE
#define SIMD_DOUBLE
#define ALIGN_DOUBLE        SSE_ALIGN_DOUBLE
#define VECSIZE_DOUBLE      SSE_VECSIZE_DOUBLE
typedef __m128d simd_double;
#define simdf64_add(x,y)    _mm_add_pd(x,y)
#define simdf64_sub(x,y)    _mm_sub_pd(x,y)
#define simdf64_mul(x,y)    _mm_mul_pd(x,y)
#define simdf64_div(x,y)    _mm_div_pd(x,y)
#define simdf64_max(x,y)    _mm_max_pd(x,y)
#define simdf64_load(x)     _mm_load_pd(x)
#define simdf64_store(x,y)  _mm_store_pd(x,y)
#define simdf64_set(x)      _mm_set1_pd(x)
#define simdf64_setzero(x)  _mm_setzero_pd()
#define simdf64_gt(x,y)     _mm_cmpgt_pd(x,y)
#define simdf64_lt(x,y)     _mm_cmplt_pd(x,y)
#define simdf64_or(x,y)     _mm_or_pd(x,y)
#define simdf64_and(x,y)    _mm_and_pd(x,y)
#define simdf64_andnot(x,y) _mm_andnot_pd(x,y)
#define simdf64_xor(x,y)    _mm_xor_pd(x,y)
#endif //SIMD_DOUBLE

// float support
#ifndef SIMD_FLOAT
#define SIMD_FLOAT
#define ALIGN_FLOAT         SSE_ALIGN_FLOAT
#define VECSIZE_FLOAT       SSE_VECSIZE_FLOAT
typedef __m128  simd_float;
#define simdf32_add(x,y)    _mm_add_ps(x,y)
#define simdf32_sub(x,y)    _mm_sub_ps(x,y)
#define simdf32_mul(x,y)    _mm_mul_ps(x,y)
#define simdf32_div(x,y)    _mm_div_ps(x,y)
#define simdf32_rcp(x)      _mm_rcp_ps(x)
#define simdf32_max(x,y)    _mm_max_ps(x,y)
#define simdf32_min(x,y)    _mm_min_ps(x,y)
#define simdf32_load(x)     _mm_load_ps(x)
#define simdf32_store(x,y)  _mm_store_ps(x,y)
#define simdf32_loadu(x)    _mm_loadu_ps(x)
#define simdf32_storeu(x,y) _mm_storeu_ps(x,y)
#define simdf32_set(x)      _mm_set1_ps(x)
#define simdf32_setzero(x)  _mm_setzero_ps()
#define simdf32_gt(x,y)     _mm_cmpgt_ps(x,y)
#define simdf32_eq(x,y)     _mm_cmpeq_ps(x,y)
#define simdf32_lt(x,y)     _mm_cmplt_ps(x,y)
#define simdf32_le(x,y)     _mm_cmple_ps(x,y)
#define simdf32_cmp(x,y,z)  _mm_cmpunord_ps(x,y)
#define simdf32_or(x,y)     _mm_or_ps(x,y)
#define simdf32_and(x,y)    _mm_and_ps(x,y)
#define simdf32_andnot(x,y) _mm_andnot_ps(x,y)
#define simdf32_xor(x,y)    _mm_xor_ps(x,y)
#define simdf32_f2i(x) 	    _mm_cvtps_epi32(x)  // convert s.p. float to integer
#define simdf_f2icast(x)    _mm_castps_si128(x) // compile time cast
#define simdf32_round(x)    _mm_round_ps(x, SIMDE_MM_FROUND_TO_NEAREST_INT | SIMDE_MM_FROUND_NO_EXC) // SSE4.1
#define simdf32_blendv_ps(x,y,z) _mm_blendv_ps(x,y,z)
#define simdf32_movemask_ps(x) _mm_movemask_ps(x)
#define simdf32_hmax(x)    simdf32_hmax_sse(x)
#define simdf32_reverse(x) simdf32_reverse_sse(x) //_mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 1, 2, 3))
#define simdf32_fmadd(x,y,z) _mm_add_ps(_mm_mul_ps(x,y),z)
#endif //SIMD_FLOAT

// integer support
#ifndef SIMD_INT
#define SIMD_INT
inline unsigned short extract_epi16(__m128i v, int pos) {
    switch(pos){
        case 0: return _mm_extract_epi16(v, 0);
        case 1: return _mm_extract_epi16(v, 1);
        case 2: return _mm_extract_epi16(v, 2);
        case 3: return _mm_extract_epi16(v, 3);
        case 4: return _mm_extract_epi16(v, 4);
        case 5: return _mm_extract_epi16(v, 5);
        case 6: return _mm_extract_epi16(v, 6);
        case 7: return _mm_extract_epi16(v, 7);
    }
    return 0;
}
#define ALIGN_INT           SSE_ALIGN_INT
#define VECSIZE_INT         SSE_VECSIZE_INT
typedef __m128i simd_int;
#define simdi32_add(x,y)    _mm_add_epi32(x,y)
#define simdi16_add(x,y)    _mm_add_epi16(x,y)
#define simdi16_sub(x,y)    _mm_sub_epi16(x,y)

#define simdi16_adds(x,y)   _mm_adds_epi16(x,y)
#define simdui8_adds(x,y)   _mm_adds_epu8(x,y)
#define simdi32_sub(x,y)    _mm_sub_epi32(x,y)
#define simdui32_subs(x,y)  _mm_max_epi32(_mm_sub_epi32(x,y), _mm_setzero_si128())
#define simdui16_subs(x,y)  _mm_subs_epu16(x,y)
#define simdui8_subs(x,y)   _mm_subs_epu8(x,y)
#define simdi32_mul(x,y)    _mm_mullo_epi32(x,y) // SSE4.1
#define simdi32_max(x,y)    _mm_max_epi32(x,y) // SSE4.1
#define simdi16_max(x,y)    _mm_max_epi16(x,y)
#define simdi16_min(x,y)    _mm_min_epi16(x,y)
#define simdi32_insert(x,y,z) _mm_insert_epi32(x,y,z)
#define simdi32_extract(x,y) _mm_extract_epi32(x,y)
#define simdi32_hmax(x)     simd_hmax32_sse(x)
#define simdi16_hmax(x)     simd_hmax16_sse(x)
#define simdui8_max(x,y)    _mm_max_epu8(x,y)
#define simdi8_hmax(x)      simd_hmax8_sse(x)
#define simd_any(x)         simd_any_sse(x)
#define simd_eq_all(x,y)    simd_eq_all_sse(x,y)
#define simdui8_avg(x,y)    _mm_avg_epu8(x,y)
#define simdui16_avg(x,y)   _mm_avg_epu16(x,y)
#define simdi_load(x)       _mm_load_si128(x)
#define simdi_loadu(x)      _mm_loadu_si128(x)
#define simdi_streamload(x) _mm_stream_load_si128(x)
#define simdi_storeu(x,y)   _mm_storeu_si128(x,y)
#define simdi_store(x,y)    _mm_store_si128(x,y)
#define simdi32_set(x)      _mm_set1_epi32(x)
#define simdi16_set(x)      _mm_set1_epi16(x)
#define simdi8_set(x)       _mm_set1_epi8(x)
#define simdi32_shuffle(x,y) _mm_shuffle_epi32(x,y)
#define simdi16_shuffle(x,y) _mm_shuffle_epi16(x,y)
#define simdi8_shuffle(x,y)  _mm_shuffle_epi8(x,y)
#define simdi_setzero()     _mm_setzero_si128()
#define simdi8_blend(x,y,z) _mm_blendv_epi8(x,y,z)
#define simdi32_gt(x,y)     _mm_cmpgt_epi32(x,y)
#define simdi8_gt(x,y)      _mm_cmpgt_epi8(x,y)
#define simdi32_eq(x,y)     _mm_cmpeq_epi32(x,y)
#define simdi16_eq(x,y)     _mm_cmpeq_epi16(x,y)
#define simdi8_eq(x,y)      _mm_cmpeq_epi8(x,y)
#define simdi32_lt(x,y)     _mm_cmplt_epi32(x,y)
#define simdi16_lt(x,y)     _mm_cmplt_epi16(x,y)
#define simdi8_lt(x,y)      _mm_cmplt_epi8(x,y)
#define simdi16_gt(x,y)     _mm_cmpgt_epi16(x,y)
#define simdi_or(x,y)       _mm_or_si128(x,y)
#define simdi_and(x,y)      _mm_and_si128(x,y)
#define simdi_andnot(x,y)   _mm_andnot_si128(x,y)
#define simdi_xor(x,y)      _mm_xor_si128(x,y)
#define simdi8_shiftl(x,y)  _mm_slli_si128(x,y)
#define simdi8_shiftr(x,y)  _mm_srli_si128(x,y)
#define SIMD_MOVEMASK_MAX   0xffff
typedef uint32_t            movemask_max_t;
#define popcount(x)         __builtin_popcount(x)
#define simdi8_movemask(x)  _mm_movemask_epi8(x)
#define simdi16_extract(x,y) extract_epi16(x,y)
#define simdi32_pack(x,y)   _mm_packs_epi32(x,y)
#define simdi16_pack(x,y)   _mm_packs_epi16(x,y)
#define simdi16_slli(x,y)	_mm_slli_epi16(x,y) // shift integers in a left by y
#define simdi16_srli(x,y)	_mm_srli_epi16(x,y) // shift integers in a right by y
#define simdi32_slli(x,y)	_mm_slli_epi32(x,y) // shift integers in a left by y
#define simdi32_srli(x,y)	_mm_srli_epi32(x,y) // shift integers in a right by y
#define simdi32_srai(x,y)   _mm_srai_epi32(x,y)
#define simdi32_i2f(x) 	    _mm_cvtepi32_ps(x)  // convert integer to s.p. float
#define simdi_i2fcast(x)    _mm_castsi128_ps(x)
#endif //SIMD_INT

// ---------------------------------------------------------------------------
// Fused compare + reduce-to-bool helpers.
// These express "is any/all lane comparison true" as a single boolean, used by
// the Striped-Smith-Waterman lazy-F early-exit and per-column change tests.
// On AVX512 the comparison stays in a k-register (compare -> scalar test on the
// mask), avoiding the compare -> all-ones vector -> movemask round-trip that the
// generic simd_any(simdiX_gt(..)) / simdi8_movemask(simdiX_gt(..)) idiom incurs.
// On AVX2/SSE they expand to exactly the previous idiom (unchanged codegen).
// ---------------------------------------------------------------------------
static inline bool simdi16_gt_any(const simd_int a, const simd_int b) {
#ifdef AVX512
    return _mm512_cmp_epi16_mask(a, b, _MM_CMPINT_NLE) != 0;
#else
    return simd_any(simdi16_gt(a, b));
#endif
}

static inline bool simdi32_gt_any(const simd_int a, const simd_int b) {
#ifdef AVX512
    return _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_NLE) != 0;
#else
    return simdi8_movemask(simdi32_gt(a, b)) != 0;
#endif
}

static inline bool simdi32_eq_all(const simd_int a, const simd_int b) {
#ifdef AVX512
    return _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ) == 0xFFFF;
#else
    return (movemask_max_t) simdi8_movemask(simdi32_eq(a, b)) == SIMD_MOVEMASK_MAX;
#endif
}

inline void *mem_align(size_t boundary, size_t size) {
    void *pointer;
    if (posix_memalign(&pointer, boundary, size) != 0) {
#define MEM_ALIGN_ERROR "mem_align could not allocate memory.\n"
        fwrite(MEM_ALIGN_ERROR, sizeof(MEM_ALIGN_ERROR), 1, stderr);
#undef MEM_ALIGN_ERROR
        exit(3);
    }
    return pointer;
}
#ifdef SIMD_FLOAT
inline simd_float * malloc_simd_float(const size_t size) {
    return (simd_float *) mem_align(ALIGN_FLOAT, size);
}
#endif
#ifdef SIMD_DOUBLE
inline simd_double * malloc_simd_double(const size_t size) {
    return (simd_double *) mem_align(ALIGN_DOUBLE, size);
}
#endif
#ifdef SIMD_INT
inline simd_int * malloc_simd_int(const size_t size) {
    return (simd_int *) mem_align(ALIGN_INT, size);
}
#endif

template <typename T>
T** malloc_matrix(int dim1, int dim2) {
#define ICEIL(x_int, fac_int) ((x_int + fac_int - 1) / fac_int) * fac_int
    // Compute mem sizes rounded up to nearest multiple of ALIGN_FLOAT
    size_t size_pointer_array = ICEIL(dim1*sizeof(T*), ALIGN_FLOAT);
    size_t dim2_padded = ICEIL(dim2*sizeof(T), ALIGN_FLOAT)/sizeof(T);

    T** matrix = (T**) mem_align( ALIGN_FLOAT, size_pointer_array + dim1*dim2_padded*sizeof(T) );
    if (matrix == NULL)
        return matrix;

    T* ptr = (T*) (matrix + (size_pointer_array/sizeof(T*)) );
    for (int i=0; i<dim1; ++i) {
        matrix[i] = ptr;
        ptr += dim2_padded;
    }
#undef ICEIL
    return matrix;
}


inline simd_float simdf32_fpow2(simd_float X) {

    simd_int* xPtr = (simd_int*) &X;    // store address of float as pointer to int

    const simd_float CONST32_05f       = simdf32_set(0.5f); // Initialize a vector (4x32) with 0.5f
    // (3 << 22) --> Initialize a large integer vector (shift left)
    const simd_int CONST32_3i          = simdi32_set(3);
    const simd_int CONST32_3shift22    = simdi32_slli(CONST32_3i, 22);
    const simd_float CONST32_1f        = simdf32_set(1.0f);
    const simd_float CONST32_FLTMAXEXP = simdf32_set(FLT_MAX_EXP);
    const simd_float CONST32_FLTMAX    = simdf32_set(FLT_MAX);
    const simd_float CONST32_FLTMINEXP = simdf32_set(FLT_MIN_EXP);
    // fifth order
    const simd_float CONST32_A = simdf32_set(0.00187682f);
    const simd_float CONST32_B = simdf32_set(0.00898898f);
    const simd_float CONST32_C = simdf32_set(0.0558282f);
    const simd_float CONST32_D = simdf32_set(0.240153f);
    const simd_float CONST32_E = simdf32_set(0.693153f);

    simd_float tx;
    simd_int lx;
    simd_float dx;
    simd_float result    = simdf32_set(0.0f);
    simd_float maskedMax = simdf32_set(0.0f);
    simd_float maskedMin = simdf32_set(0.0f);

    // Check wheter one of the values is bigger or smaller than FLT_MIN_EXP or FLT_MAX_EXP
    // The correct FLT_MAX_EXP value is written to the right place
    maskedMax = simdf32_gt(X, CONST32_FLTMAXEXP);
    maskedMin = simdf32_gt(X, CONST32_FLTMINEXP);
    maskedMin = simdf32_xor(maskedMin, maskedMax);
    // If a value is bigger than FLT_MAX_EXP --> replace the later result with FLTMAX
    maskedMax = simdf32_and(CONST32_FLTMAX, simdf32_gt(X, CONST32_FLTMAXEXP));

    tx = simdf32_add((simd_float ) CONST32_3shift22, simdf32_sub(X, CONST32_05f)); // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
    // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127),
    // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)

    lx = simdf32_f2i(tx);                                       // integer value of x

    dx = simdf32_sub(X, simdi32_i2f(lx));                       // float remainder of x

    //   x = 1.0f + dx*(0.693153f             // polynomial apporoximation of 2^x for x in the range [0, 1]
    //            + dx*(0.240153f             // Gives relative deviation < 2.3E-7
    //            + dx*(0.0558282f            // Speed: 2.3E-8s
    //            + dx*(0.00898898f
    //            + dx* 0.00187682f ))));
    X = simdf32_mul(dx, CONST32_A);
    X = simdf32_add(CONST32_B, X);  // add constant B
    X = simdf32_mul(dx, X);
    X = simdf32_add(CONST32_C, X);  // add constant C
    X = simdf32_mul(dx, X);
    X = simdf32_add(CONST32_D, X);  // add constant D
    X = simdf32_mul(dx, X);
    X = simdf32_add(CONST32_E, X);  // add constant E
    X = simdf32_mul(dx, X);
    X = simdf32_add(X, CONST32_1f); // add 1.0f

    simd_int lxExp = simdi32_slli(lx, 23); // add integer power of 2 to exponent

    *xPtr = simdi32_add(*xPtr, lxExp); // add integer power of 2 to exponent

    // Add all Values that are greater than min and less than max
    result = simdf32_and(maskedMin, X);
    // Add MAX_FLT values where entry values were > FLT_MAX_EXP
    result = simdf32_or(result, maskedMax);

    return result;
}

static inline simd_float polynomial_5 (simd_float x, simd_float c0, simd_float c1, simd_float c2, simd_float c3, simd_float c4, simd_float c5) {
    simd_float x2 = simdf32_mul(x, x);
    simd_float x4 = simdf32_mul(x2, x2);
    return simdf32_fmadd(simdf32_fmadd(c3, x, c2), x2, simdf32_fmadd(simdf32_fmadd(c5, x, c4), x4, simdf32_fmadd(c1, x, c0)));
}

static inline simd_float polynomial_8(simd_float x, simd_float c0, simd_float c1, simd_float c2, simd_float c3, simd_float c4, simd_float c5, simd_float c6, simd_float c7, simd_float c8) {
    simd_float x2 = simdf32_mul(x, x);
    simd_float x4 = simdf32_mul(x2, x2);
    simd_float x8 = simdf32_mul(x4, x4);

    return simdf32_fmadd(simdf32_fmadd(simdf32_fmadd(c7, x, c6), x2, simdf32_fmadd(c5, x, c4)), x4,
        simdf32_fmadd(simdf32_fmadd(c3, x, c2), x2, simdf32_fmadd(c1, x, c0) + simdf32_mul(x8, c8)));
}

static inline simd_float simdf32_pow2n(simd_float n) {
    simd_float a = simdf32_add(n, simdf32_set(127.0f + 8388608.0f));
    simd_int b = simdf_f2icast(a);
    simd_int c = simdi32_slli(b, 23);
    return simdi_i2fcast(c);
}

// Robot idea, seems to make FwBw faster on AVX512
static inline simd_float simdf32_abs(simd_float a) {
  #ifdef AVX512
      const __m512i mask = _mm512_set1_epi32(0x7fffffff);
      return _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(a), mask));
  #else
      const simd_float mask = simdi_i2fcast(simdi32_set(0x7FFFFFFF));
      return simdf32_and(a, mask);
  #endif
}

static inline simd_int simdi32_is_finite(simd_float a) {
    simd_int t1 = simdf_f2icast(a);
    simd_int t2 = simdi32_slli(t1, 1);

    simd_int exponent_mask = simdi32_set(0xFF000000);
    simd_int t3 = simdi_and(t2, exponent_mask);

    simd_int all_ones = simdi32_set(0xFF000000);
    simd_int is_not_inf_or_nan = simdi32_eq(t3, all_ones);
    simd_int result = simdi_xor(is_not_inf_or_nan, simdi32_set(-1)); 

    return result;
}

static inline simd_float simdf32_exp(simd_float x_init) {
    const simd_float P0 = simdf32_set(1.0f / 2.0f);
    const simd_float P1 = simdf32_set(1.0f / 6.0f);
    const simd_float P2 = simdf32_set(1.0f / 24.0f);
    const simd_float P3 = simdf32_set(1.0f / 120.0f);
    const simd_float P4 = simdf32_set(1.0f / 720.0f);
    const simd_float P5 = simdf32_set(1.0f / 5040.0f);

    const simd_float negLN2_HI = simdf32_set(-0.693359375f);
    const simd_float negLN2_LO = simdf32_set(2.12194440e-4f);
    const simd_float VM_LOG2E = simdf32_set(1.44269504088896340736);

    simd_float x = x_init;
    simd_float r = simdf32_round(simdf32_mul(x_init, VM_LOG2E));
    x = simdf32_fmadd(r, negLN2_HI, x);
    x = simdf32_fmadd(r, negLN2_LO, x);
    simd_float x2 = simdf32_mul(x, x);
    simd_float z = polynomial_5(x, P0, P1, P2, P3, P4, P5);
    z = simdf32_fmadd(z, x2, x);
    simd_float n2 = simdf32_pow2n(r);
    z = simdf32_fmadd(z, n2, n2);

    #ifdef AVX512
        //special cases
        const simd_float MAX_X = simdf32_set(87.3f);
        __mmask16 inrange = _mm512_cmp_ps_mask(simdf32_abs(x_init), MAX_X, _CMP_LT_OS);
        const simd_float InfVec = simdf32_set(std::numeric_limits<float>::infinity());
        const __mmask16 signBit = _mm512_cmp_epi32_mask(_mm512_castps_si512(x_init), _mm512_setzero_si512(), _MM_CMPINT_LT);
        __mmask16 isNan = _mm512_cmp_ps_mask(x_init, x_init, 3);
        r = _mm512_mask_blend_ps(signBit, InfVec, simdf32_set(0.0f)); // value in case of -
        z = _mm512_mask_blend_ps(inrange, r, z);     // +/- underflow
        z = _mm512_mask_blend_ps(isNan, z, x_init); // NAN goes through
    #else
        //special cases
        const simd_float MAX_X = simdf32_set(87.3f);
        simd_float inrange = simdf32_lt(simdf32_abs(x_init), MAX_X);
        const simd_float InfVec = simdf32_set(std::numeric_limits<float>::infinity());
        simd_float signBit = simdi_i2fcast(simdi32_srai(simdf_f2icast(x_init), 31));
        simd_float isNan = simdf32_cmp(x_init, x_init, 3);
        r = simdf32_blendv_ps(InfVec, simdf32_set(0.0f), signBit); // value in case of -
        z = simdf32_blendv_ps(r, z, inrange);     // +/- underflow
        z = simdf32_blendv_ps(z, x_init, isNan); // NAN goes through
    #endif
    return z;
}

#ifdef AVX512
static inline simd_float simdf32_log(simd_float x_init) {
    // Same as AVX2 but uses masked to avoid having to convert mask to vector back and forth 
    // Constants
    const simd_float LN2f_HI = simdf32_set(0.693359375f);
    const simd_float LN2f_LO = simdf32_set(-2.12194440e-4f);
    const simd_float P0LOGF  = simdf32_set(3.3333331174E-1f);
    const simd_float P1LOGF  = simdf32_set(-2.4999993993E-1f);
    const simd_float P2LOGF  = simdf32_set(2.0000714765E-1f);
    const simd_float P3LOGF  = simdf32_set(-1.6668057665E-1f);
    const simd_float P4LOGF  = simdf32_set(1.4249322787E-1f);
    const simd_float P5LOGF  = simdf32_set(-1.2420140846E-1f);
    const simd_float P6LOGF  = simdf32_set(1.1676998740E-1f);
    const simd_float P7LOGF  = simdf32_set(-1.1514610310E-1f);
    const simd_float P8LOGF  = simdf32_set(7.0376836292E-2f);
    const simd_float SQRT2_threshold   = simdf32_set(1.41421356237309504880*0.5);
    const simd_float one = simdf32_set(1.0f);

    // separate mantissa from exponent
    simd_int xi = simdf_f2icast(x_init);
    simd_int mi = simdi_or(simdi_and(xi, simdi32_set(0x007FFFFF)), simdi32_set(0x3F000000));
    simd_float m = simdi_i2fcast(mi);

    simd_int ei = simdi32_sub(simdi32_srli(simdi32_slli(xi, 1), 24), simdi32_set(0x7F));
    simd_float e = simdi32_i2f(ei);
    __mmask16 blend_mask = _mm512_cmp_ps_mask(m, SQRT2_threshold, _CMP_GT_OS);
    __mmask16 not_blend_mask = _mm512_cmp_ps_mask(m, SQRT2_threshold, _CMP_LE_OQ);
    simd_float m_2 = simdf32_add(m, m);

    m = _mm512_mask_blend_ps(not_blend_mask, m, m_2);
    m = simdf32_sub(m, one);

    simd_float e_1 = simdf32_add(e, one);
    e = _mm512_mask_blend_ps(blend_mask, e, e_1);

    simd_float res = polynomial_8(m, P0LOGF, P1LOGF, P2LOGF, P3LOGF, P4LOGF, P5LOGF, P6LOGF, P7LOGF, P8LOGF);
    simd_float m2 = simdf32_mul(m, m);
    res = simdf32_mul(res, simdf32_mul(m2, m));

    res = simdf32_fmadd(e, LN2f_LO, res);
    res = simdf32_add(res, simdf32_sub(m, simdf32_mul(m2, simdf32_set(0.5f))));
    res = simdf32_fmadd(e, LN2f_HI, res);

    // Special cases
    const simd_float VM_SMALLEST_NORMALF = simdf32_set(1.17549435e-38f);
    __mmask16 overflow =_mm512_cmp_epi32_mask(simdi32_is_finite(x_init), _mm512_setzero_si512(), _MM_CMPINT_EQ);
    __mmask16 underflow = _mm512_cmp_ps_mask(x_init, VM_SMALLEST_NORMALF, _CMP_LT_OS);

    const simd_float negNanVec = simdf32_set(-std::numeric_limits<float>::quiet_NaN());
    const simd_float negInfVec = simdf32_set(-std::numeric_limits<float>::infinity());
    // if overflow(+- INF or NaN) gives x_init
    res = _mm512_mask_blend_ps(overflow, res, x_init);
    // if underflow(<1.17549435e-38f) gives -NAN
    res = _mm512_mask_blend_ps(underflow, res, negNanVec);
    simd_int x_exponent = simdi_and(xi, simdi32_set(0x7F800000));
    // if x == 0 or subnormal gives -INF
    __mmask16 maskZeroOrSubnormal = _mm512_cmp_epi32_mask(x_exponent, _mm512_setzero_si512(), _MM_CMPINT_EQ); // x == 0 or subnormal
    res = _mm512_mask_blend_ps(maskZeroOrSubnormal, res, negInfVec );
    return res;
}
#else
static inline simd_float simdf32_log(simd_float x_init) {
    // Constants
    const simd_float LN2f_HI = simdf32_set(0.693359375f);
    const simd_float LN2f_LO = simdf32_set(-2.12194440e-4f);
    const simd_float P0LOGF  = simdf32_set(3.3333331174E-1f);
    const simd_float P1LOGF  = simdf32_set(-2.4999993993E-1f);
    const simd_float P2LOGF  = simdf32_set(2.0000714765E-1f);
    const simd_float P3LOGF  = simdf32_set(-1.6668057665E-1f);
    const simd_float P4LOGF  = simdf32_set(1.4249322787E-1f);
    const simd_float P5LOGF  = simdf32_set(-1.2420140846E-1f);
    const simd_float P6LOGF  = simdf32_set(1.1676998740E-1f);
    const simd_float P7LOGF  = simdf32_set(-1.1514610310E-1f);
    const simd_float P8LOGF  = simdf32_set(7.0376836292E-2f);
    const simd_float SQRT2_threshold   = simdf32_set(1.41421356237309504880*0.5);
    const simd_float one = simdf32_set(1.0f);

    // separate mantissa from exponent
    simd_int xi = simdf_f2icast(x_init);
    simd_int mi = simdi_or(simdi_and(xi, simdi32_set(0x007FFFFF)), simdi32_set(0x3F000000));
    simd_float m = simdi_i2fcast(mi);

    simd_int ei = simdi32_sub(simdi32_srli(simdi32_slli(xi, 1), 24), simdi32_set(0x7F));
    simd_float e = simdi32_i2f(ei);
    simd_float blend_mask = simdf32_gt(m, SQRT2_threshold);
    simd_float not_blend_mask = simdf32_le(m, SQRT2_threshold);
    simd_float m_2 = simdf32_add(m, m);

    m = simdf32_blendv_ps(m, m_2, not_blend_mask);
    m = simdf32_sub(m, one);

    simd_float e_1 = simdf32_add(e, one);
    e = simdf32_blendv_ps(e, e_1, blend_mask);

    simd_float res = polynomial_8(m, P0LOGF, P1LOGF, P2LOGF, P3LOGF, P4LOGF, P5LOGF, P6LOGF, P7LOGF, P8LOGF);
    simd_float m2 = simdf32_mul(m, m);
    res = simdf32_mul(res, simdf32_mul(m2, m));

    res = simdf32_fmadd(e, LN2f_LO, res);
    res = simdf32_add(res, simdf32_sub(m, simdf32_mul(m2, simdf32_set(0.5f))));
    res = simdf32_fmadd(e, LN2f_HI, res);

    // Special cases
    const simd_float VM_SMALLEST_NORMALF = simdf32_set(1.17549435e-38f);
    simd_float overflow = simdi_i2fcast(simdi_xor(simdi32_is_finite(x_init), simdi32_set(-1)));
    simd_float underflow = simdf32_lt(x_init, VM_SMALLEST_NORMALF);

    const simd_float negNanVec = simdf32_set(-std::numeric_limits<float>::quiet_NaN());
    const simd_float negInfVec = simdf32_set(-std::numeric_limits<float>::infinity());
    // if overflow(+- INF or NaN) gives x_init
    res = simdf32_blendv_ps(res, x_init, overflow);
    // if underflow(<1.17549435e-38f) gives -NAN
    res = simdf32_blendv_ps(res, negNanVec, underflow);
    simd_int x_exponent = simdi_and(xi, simdi32_set(0x7F800000));
    // if x == 0 or subnormal gives -INF
    simd_float maskZeroOrSubnormal = simdi_i2fcast(simdi32_eq(x_exponent, simdi32_set(0))); // x == 0 or subnormal
    res = simdf32_blendv_ps(res, negInfVec, maskZeroOrSubnormal);
    return res;
}
#endif

inline float ScalarProd20(const float* qi, const float* tj) {
//#ifdef AVX
//  float __attribute__((aligned(ALIGN_FLOAT))) res;
//  __m256 P; // query 128bit SSE2 register holding 4 floats
//  __m256 S; // aux register
//  __m256 R; // result
//  __m256* Qi = (__m256*) qi;
//  __m256* Tj = (__m256*) tj;
//
//  R = _mm256_mul_ps(*(Qi++),*(Tj++));
//  P = _mm256_mul_ps(*(Qi++),*(Tj++));
//  S = _mm256_mul_ps(*Qi,*Tj); // floats A, B, C, D, ?, ?, ? ,?
//  R = _mm256_add_ps(R,P);     // floats 0, 1, 2, 3, 4, 5, 6, 7
//  P = _mm256_permute2f128_ps(R, R, 0x01); // swap hi and lo 128 bits: 4, 5, 6, 7, 0, 1, 2, 3
//  R = _mm256_add_ps(R,P);     // 0+4, 1+5, 2+6, 3+7, 0+4, 1+5, 2+6, 3+7
//  R = _mm256_add_ps(R,S);     // 0+4+A, 1+5+B, 2+6+C, 3+7+D, ?, ?, ? ,?
//  R = _mm256_hadd_ps(R,R);    // 04A15B, 26C37D, ?, ?, 04A15B, 26C37D, ?, ?
//  R = _mm256_hadd_ps(R,R);    // 01234567ABCD, ?, 01234567ABCD, ?, 01234567ABCD, ?, 01234567ABCD, ?
//  _mm256_store_ps(&res, R);
//  return res;
//#else
//
//
//TODO fix this
    float __attribute__((aligned(16))) res;
    __m128 P; // query 128bit SSE2 register holding 4 floats
    __m128 R;// result
    __m128* Qi = (__m128*) qi;
    __m128* Tj = (__m128*) tj;

    __m128 P1 = _mm_mul_ps(*(Qi),*(Tj));
    __m128 P2 = _mm_mul_ps(*(Qi+1),*(Tj+1));
    __m128 R1 = _mm_add_ps(P1, P2);

    __m128 P3 = _mm_mul_ps(*(Qi + 2), *(Tj + 2));
    __m128 P4 = _mm_mul_ps(*(Qi + 3), *(Tj + 3));
    __m128 R2 = _mm_add_ps(P3, P4);
    __m128 P5 = _mm_mul_ps(*(Qi+4), *(Tj+4));

    R = _mm_add_ps(R1, R2);
    R = _mm_add_ps(R,P5);

//    R = _mm_hadd_ps(R,R);
//    R = _mm_hadd_ps(R,R);
    P = _mm_shuffle_ps(R, R, _MM_SHUFFLE(2,0,2,0));
    R = _mm_shuffle_ps(R, R, _MM_SHUFFLE(3,1,3,1));
    R = _mm_add_ps(R,P);
    P = _mm_shuffle_ps(R, R, _MM_SHUFFLE(2,0,2,0));
    R = _mm_shuffle_ps(R, R, _MM_SHUFFLE(3,1,3,1));
    R = _mm_add_ps(R,P);
    _mm_store_ss(&res, R);
    return res;
//#endif
//    return tj[0] * qi[0] + tj[1] * qi[1] + tj[2] * qi[2] + tj[3] * qi[3]
//            + tj[4] * qi[4] + tj[5] * qi[5] + tj[6] * qi[6] + tj[7] * qi[7]
//            + tj[8] * qi[8] + tj[9] * qi[9] + tj[10] * qi[10] + tj[11] * qi[11]
//            + tj[12] * qi[12] + tj[13] * qi[13] + tj[14] * qi[14]
//            + tj[15] * qi[15] + tj[16] * qi[16] + tj[17] * qi[17]
//            + tj[18] * qi[18] + tj[19] * qi[19];
}

#endif //SIMD_H
