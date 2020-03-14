// sse2altivec is still very incomplete
// licensed under GPLv3 see LICENCE file
#ifndef SSE2ALTIVEC
#define SSE2ALTIVEC

// ignore all warnings
#pragma GCC system_header

#include <altivec.h>
#define SSE 1

typedef __vector double __m128d;
typedef __vector float __m128;
typedef __vector int __m128i;

typedef __vector   signed char simd_s8;
typedef __vector unsigned char simd_u8;
typedef __vector   signed short simd_s16;
typedef __vector unsigned short simd_u16;
typedef __vector int64_t simd_s64;
typedef __vector uint64_t simd_u64;

#define _mm_add_ps(x,y)          (__m128)((__m128)(x) + (__m128)(y))
#define _mm_sub_ps(x,y)          (__m128)((__m128)(x) - (__m128)(y))
#define _mm_mul_ps(x,y)          (__m128)((__m128)(x) * (__m128)(y))
#define _mm_div_ps(x,y)          (__m128)((__m128)(x) / (__m128)(y))
#define _mm_rcp_ps(x)            (__m128)vec_re((__m128)(x))
#define _mm_max_ps(x,y)          (__m128)vec_max((__m128)(x),(__m128)(y))
#define _mm_min_ps(x,y)          (__m128)vec_min((__m128)(x),(__m128)(y))
#define _mm_load_ps(x)           (__m128)vec_vsx_ld(0, (__m128 const*)(x))
#define _mm_store_ps(x,y)        vec_vsx_st((__m128)(y),0,(__m128*)(x))
#define _mm_store_ss(x,y)        vec_vsx_st((__m128)(y),0,(__m128*)(x))
#define _mm_set1_ps(x)           (__m128)vec_splats((float)(x))
#define _mm_setzero_ps(x)        (__m128)vec_splats((float)0)
#define _mm_cmpgt_ps(x,y)        (__m128)vec_cmpgt((__m128)(x),(__m128)(y))
#define _mm_cmpeq_ps(x,y)        (__m128)vec_cmpeq((__m128)(x),(__m128)(y))
#define _mm_cmplt_ps(x,y)        (__m128)vec_cmplt((__m128)(x),(__m128)(y))
#define _mm_or_ps(x,y)           (__m128)vec_or((__m128)(x),(__m128)(y))
#define _mm_and_ps(x,y)          (__m128)vec_and((__m128)(x),(__m128)(y))
#define _mm_andnot_ps(x,y)       (__m128)vec_andc((__m128)(x),(__m128)(y))
#define _mm_xor_ps(x,y)          (__m128)vec_xor((__m128)(x),(__m128)(y))
#define _mm_cvtps_epi32(x)       (__m128i)vec_cts((x),0)
#define _mm_castps_si128(x)      (__m128i)(x)
#define _mm_add_epi32(x,y)       (__m128i)vec_add((__m128i)(x),(__m128i)(y))
#define _mm_add_epi16(x,y)       (__m128i)vec_add((simd_s16)(x),(simd_s16)(y))
#define _mm_add_epi8(x,y)        (__m128i)vec_add((simd_s8)(x),(simd_s8)(y))
#define _mm_adds_epi16(x,y)      (__m128i)vec_adds((simd_s16)(x),(simd_s16)(y))
#define _mm_adds_epu8(x,y)       (__m128i)vec_adds((simd_u8)(x),(simd_u8)(y))
#define _mm_sub_epi32(x,y)       (__m128i)vec_sub((__m128i)(x),(__m128i)(y))
#define _mm_sub_epi16(x,y)       (__m128i)vec_sub((simd_s16)(x),(simd_s16)(y))
#define _mm_sub_epi8(x,y)        (__m128i)vec_sub((simd_s8)(x),(simd_s8)(y))
#define _mm_subs_epu16(x,y)      (__m128i)vec_subs((simd_u16)(x),(simd_u16)(y))
#define _mm_subs_epu8(x,y)       (__m128i)vec_subs((simd_u8)(x),(simd_u8)(y))
#define _mm_mullo_epi32(x,y)     (__m128i)vec_mul((__m128i)(x),(__m128i)(y))
#define _mm_max_epi32(x,y)       (__m128i)vec_max((__m128i)(x),(__m128i)(y))
#define _mm_max_epi16(x,y)       (__m128i)vec_max((simd_s16)(x),(simd_s16)(y))
#define _mm_max_epu8(x,y)        (__m128i)vec_max((simd_u8)(x),(simd_u8)(y))
#define _mm_min_epu8(x,y)        (__m128i)vec_min((simd_u8)(x),(simd_u8)(y))
#define _mm_load_si128(x)        (__m128i)vec_vsx_ld(0,(__m128i const*)(x))
#define _mm_loadu_si128(x)       (__m128i)vec_vsx_ld(0,(__m128i const*)(x))
#define _mm_storeu_si128(x,y)    vec_vsx_st((__m128i)(y),0,(__m128i*)(x))
#define _mm_store_si128(x,y)     vec_vsx_st((__m128i)(y),0,(__m128i*)(x))
#define _mm_set1_epi32(x)        (__m128i)vec_splats((signed int)(x))
#define _mm_set1_epi16(x)        (__m128i)vec_splats((signed short)(x))
#define _mm_set1_epi8(x)         (__m128i)vec_splats((signed char)(x))
#define _mm_setzero_si128(x)     (__m128i)vec_splats(0)
#define _mm_cmpgt_epi32(x,y)     (__m128i)vec_cmpgt((__m128i)(x),(__m128i)(y))
#define _mm_cmpgt_epi16(x,y)     (__m128i)vec_cmpgt((simd_s16)(x),(simd_s16)(y))
#define _mm_cmpgt_epi8(x,y)      (__m128i)vec_cmpgt((simd_s8)(x),(simd_s8)(y))
#define _mm_cmpeq_epi32(x,y)     (__m128i)vec_cmpeq((__m128i)(x),(__m128i)(y))
#define _mm_cmpeq_epi16(x,y)     (__m128i)vec_cmpeq((simd_s16)(x),(simd_s16)(y))
#define _mm_cmpeq_epi8(x,y)      (__m128i)vec_cmpeq((simd_s8)(x),(simd_s8)(y))
#define _mm_cmplt_epi32(x,y)     (__m128i)vec_cmplt((__m128i)(x),(__m128i)(y))
#define _mm_cmplt_epi16(x,y)     (__m128i)vec_cmplt((simd_s16)(x),(simd_s16)(y))
#define _mm_cmplt_epi8(x,y)      (__m128i)vec_cmplt((simd_s8)(x),(simd_s8)(y))
#define _mm_or_si128(x,y)        (__m128i)vec_or((__m128i)(x),(__m128i)(y))
#define _mm_and_si128(x,y)       (__m128i)vec_and((__m128i)(x),(__m128i)(y))
#define _mm_andnot_si128(x,y)    (__m128i)vec_andc((__m128i)(x),(__m128i)(y))
#define _mm_xor_si128(x,y)       (__m128i)vec_xor((__m128i)(x),(__m128i)(y))
#define _mm_extract_epi16(x,imm) (int16_t)vec_extract((simd_s16)(x),(imm))
#define _mm_extract_epi8(x,imm)  (int8_t)vec_extract((simd_s8)(x),(imm))
#define _mm_slli_epi16(x,y)      (simd_s16)vec_sl((simd_s16)(x),vec_splats((unsigned short)(y)))
#define _mm_srli_epi16(x,y)      (simd_s16)vec_sr((simd_s16)(x),vec_splats((unsigned short)(y)))
#define _mm_slli_epi32(x,y)      (__m128i)vec_sl((__m128i)(x),vec_splats((unsigned int)(y)))
#define _mm_srli_epi32(x,y)      (__m128i)vec_sr((__m128i)(x),vec_splats((unsigned int)(y)))
#define _mm_cvtepi32_ps(x)       (__m128)vec_ctf((__m128i)(x),0)
#define _mm_castsi128_ps(x)      (__m128)(x)
#define _mm_slli_si128(x,y)      (__m128i)vec_slo((simd_u8)(x),(simd_u8)vec_splats((char)(y << 3)))
#define _mm_srli_si128(x,y)      (__m128i)vec_sro((simd_u8)(x),(simd_u8)vec_splats((char)(y << 3)))
#define _mm_cvtsi128_si64(a)     (int64_t)vec_extract((simd_s64)(a),0)
#define _mm_cvtsi128_si32(a)     (int32_t)vec_extract((__m128i)(a),0)
#define _mm_cvtsi64_si128(a)     (__m128i)((simd_s64){(int64_t)(a),0})
#define _mm_cvtsi32_si128(a)     (__m128i){(int)(a),0,0,0}
#define _mm_packs_epi32(x,y)     (simd_s16)vec_packs((__m128i)(x), (__m128i)(y))
#define _mm_packus_epi16(x,y)    (simd_u8)vec_packsu((simd_s16)(x), (simd_s16)(y))
#define _mm_set_epi32(e3,e2,e1,e0)  (__m128i){(e0),(e1),(e2),(e3)}
#define _mm_setr_epi32(e3,e2,e1,e0) (__m128i){(e3),(e2),(e1),(e0)}
#define _mm_set_epi16(e7,e6,e5,e4,e3,e2,e1,e0) \
                                    (__m128i)((simd_s16){(e0),(e1),(e2),(e3),(e4),(e5),(e6),(e7)})
#define _mm_setr_epi16(e7,e6,e5,e4,e3,e2,e1,e0) \
                                    (__m128i)((simd_s16){(e7),(e6),(e5),(e4),(e3),(e2),(e1),(e0)})
#define _mm_set_epi8(e15,e14,e13,e12,e11,e10,e9,e8,e7,e6,e5,e4,e3,e2,e1,e0) \
                                    (__m128i)((simd_s8){(e0),(e1),(e2),(e3),(e4),(e5),(e6),(e7),(e8),(e9),(e10),(e11),(e12),(e13),(e14),(e15)})
#define _mm_setr_epi8(e15,e14,e13,e12,e11,e10,e9,e8,e7,e6,e5,e4,e3,e2,e1,e0) \
                                    (__m128i)((simd_s8){(e15),(e14),(e13),(e12),(e11),(e10),(e9),(e8),(e7),(e6),(e5),(e4),(e3),(e2),(e1),(e0)})

static inline void _mm_storel_epi64(__m128i* mem_addr, __m128i a) {
    *((int64_t*)mem_addr) = (int64_t)vec_extract((simd_s64)(a), 0);
}

// From OpenCV
// https://github.com/opencv/opencv/pull/15235
// 3-Clause BSD License
static inline unsigned short _mm_movemask_epi8(__m128i value) {
    static const simd_u8 perm = {120, 112, 104, 96, 88, 80, 72, 64, 56, 48, 40, 32, 24, 16, 8, 0};
    return vec_extract((__m128i)vec_vbpermq((simd_u8)value, perm), 2);
}

// From reedsolomon
// https://github.com/NicolasT/reedsolomon/blob/master/cbits/reedsolomon.c
// MIT License
static inline __m128i _mm_shuffle_epi8(__m128i a, __m128i b) {
  const __m128i zero = (__m128i)vec_splats((unsigned char)0);
  return (__m128i)vec_perm((simd_u8)a, (simd_u8)zero, (simd_u8)b);
}

#endif
