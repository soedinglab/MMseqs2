// Placeholder half types (named in type traits only).
#ifndef CUDA_FP16_H_METAL_SHIM
#define CUDA_FP16_H_METAL_SHIM

#include <cstdint>

using __half = _Float16;

struct __half2 {
    _Float16 x;
    _Float16 y;
};

using half  = __half;
using half2 = __half2;

static_assert(sizeof(__half) == 2, "half must be 2 bytes");
static_assert(sizeof(__half2) == 4, "half2 must be 4 bytes");

inline __half  __float2half(float f)             { return static_cast<_Float16>(f); }
inline float   __half2float(__half h)            { return static_cast<float>(h); }
inline __half2 __halves2half2(__half a, __half b) { return __half2{a, b}; }
inline __half  __low2half(__half2 h)             { return h.x; }
inline __half  __high2half(__half2 h)            { return h.y; }
inline __half2 make_half2(__half a, __half b)    { return __half2{a, b}; }

#endif
