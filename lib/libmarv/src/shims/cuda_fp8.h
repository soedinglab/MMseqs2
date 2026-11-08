// Placeholder fp8 types (named in type traits only).
#ifndef CUDA_FP8_H_METAL_SHIM
#define CUDA_FP8_H_METAL_SHIM

#include <cstdint>

struct __nv_fp8_e4m3   { std::uint8_t  __x; };
struct __nv_fp8_e5m2   { std::uint8_t  __x; };
struct __nv_fp8x2_e4m3 { std::uint16_t __x; };

#endif
