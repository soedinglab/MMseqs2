#ifndef LIBMARV_PTX_WRAPPERS_CUH
#define LIBMARV_PTX_WRAPPERS_CUH

#if !defined(__HIPCC__)
#include <cuda_fp16.h>
#else
#include <hip/hip_fp16.h>
#endif

#include "namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN


#if (defined(__CUDA_ARCH_FAMILY_SPECIFIC__) && (__CUDA_ARCH_FAMILY_SPECIFIC__ >= 1200))
    #define IS_BLACKWELL_120F_PTX
#endif

#if (__CUDACC_VER_MAJOR__ > 13) || (__CUDACC_VER_MAJOR__ == 13 && __CUDACC_VER_MINOR__ >= 2)
    #define IS_AT_LEAST_CUDA_13_2
#endif

#if defined(IS_BLACKWELL_120F_PTX) && defined(IS_AT_LEAST_CUDA_13_2) 
    #define HAS_BLACKWELL_INT8_PTX
#endif    

extern __device__ void blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();


__device__ __forceinline__
unsigned int ptx_add_u8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("add.u8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

__device__ __forceinline__
unsigned int ptx_add_sat_u8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("add.sat.u8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

__device__ __forceinline__
unsigned int ptx_sub_u8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("sub.u8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

__device__ __forceinline__
unsigned int ptx_sub_sat_u8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("sub.sat.u8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

__device__ __forceinline__
unsigned int ptx_max_u8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("max.u8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

__device__ __forceinline__
unsigned int ptx_min_u8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("min.u8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}


__device__ __forceinline__
unsigned int ptx_add_s8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("add.s8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

__device__ __forceinline__
unsigned int ptx_add_sat_s8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("add.sat.s8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

__device__ __forceinline__
unsigned int ptx_sub_s8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("sub.s8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

__device__ __forceinline__
unsigned int ptx_sub_sat_s8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("sub.sat.s8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

__device__ __forceinline__
unsigned int ptx_max_s8x4(unsigned int a, unsigned int b){
    #ifdef HAS_BLACKWELL_INT8_PTX
    unsigned int result;
    asm("max.s8x4 %0, %1, %2;\n\t" : "=r"(result) : "r"(a), "r"(b));
    return result;
    #else
    (void)a;
    (void)b;
    blackwells_int8_ptx_instructions_are_only_available_for_sm120f_or_higher_since_cuda13_2();
    return 0;
    #endif
}

extern __device__ void ptx_add_f32_f16_is_only_available_for_sm100_or_higher();
__device__ __forceinline__
float ptx_add_f32_f16(float b, half a){
    #if __CUDA_ARCH__ >= 1000
    static_assert(sizeof(half) == sizeof(unsigned short));
    unsigned short as;
    memcpy(&as, &a, sizeof(unsigned short));
    float result;
    asm("add.f32.f16  %0, %1, %2;\n\t" : "=f"(result) : "h"(as), "f"(b));
    return result;
    #else
    (void)a;
    (void)b;
    ptx_add_f32_f16_is_only_available_for_sm100_or_higher();
    return 0;
    #endif
}


LIBMARV_NAMESPACE_WITH_NESTING_END


#endif