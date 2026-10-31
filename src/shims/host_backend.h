// Force-included into the Metal and SIMD builds: lets libmarv's CUDA host code compile as plain C++ on the shims in src/shims.
#ifndef MARV_HOST_BACKEND_H
#define MARV_HOST_BACKEND_H

#define __global__
#define __device__
#define __host__
#define __forceinline__ inline

#ifndef NO_NVTOOLSEXT
#define NO_NVTOOLSEXT 1
#endif

#include <iostream>
#include <cstdlib>
// cuda_helpers.cuh defines CUERR only for CUDA/HIP compilers
#define CUERR {                                                              \
    cudaError_t err;                                                         \
    if ((err = cudaGetLastError()) != cudaSuccess) {                         \
        std::cout << "libmarv error: " << cudaGetErrorString(err) << " : "     \
                  << __FILE__ << ", line " << __LINE__ << std::endl;         \
        exit(1);                                                             \
    }                                                                        \
}

#include <sys/mman.h>
#ifndef MAP_POPULATE
#define MAP_POPULATE 0
#endif

#include <cstdint>
#include <cstddef>

// referenced by kernel templates that are parsed but never instantiated
struct dim3 {
    unsigned int x, y, z;
    constexpr dim3(unsigned int x_ = 1, unsigned int y_ = 1, unsigned int z_ = 1) : x(x_), y(y_), z(z_) {}
};
inline dim3 threadIdx{0, 0, 0};
inline dim3 blockIdx{0, 0, 0};
inline dim3 blockDim{1, 1, 1};
inline dim3 gridDim{1, 1, 1};

struct char2  { signed char x, y; };
struct char4  { signed char x, y, z, w; };
struct short2 { short x, y; };
struct short4 { short x, y, z, w; };
struct alignas(8)  int2   { int x, y; };
struct alignas(16) int4   { int x, y, z, w; };
struct alignas(8)  float2 { float x, y; };
struct alignas(16) float4 { float x, y, z, w; };

#include <cuda_fp16.h>
#include <cuda_fp8.h>
#include <cuda_runtime.h>
#include <cuda/std/span>
#include <cuda/std/mdspan>
#include <cuda/std/tuple>
#include <cuda/std/array>
#include <cuda/std/iterator>
#include <cuda/cmath>
#include <cuda/functional>
#include <thrust/thrust_impl.h>
#include <cub/cub.cuh>

#endif // MARV_HOST_BACKEND_H
