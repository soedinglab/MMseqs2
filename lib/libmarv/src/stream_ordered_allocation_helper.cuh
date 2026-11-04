#ifndef LIBMARV_STREAM_ORDERED_ALLOCATION_HELPER
#define LIBMARV_STREAM_ORDERED_ALLOCATION_HELPER

#include "cuda_hip_rename.h"

#include "namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

inline
cudaError_t helper_gpuMallocAsync(void** devPtr, size_t size, cudaStream_t hStream){
    #if defined(__CUDACC__)
        return cudaMallocAsync(devPtr, size, hStream);
    #else
        //AMD's stream-order allocation API is beta and does not always work correctly, use synchronous function
        return cudaMalloc(devPtr, size);
    #endif
}

inline
cudaError_t helper_gpuFreeAsync(void* devPtr, cudaStream_t hStream){
    #if defined(__CUDACC__)
        return cudaFreeAsync(devPtr, hStream);
    #else
        //AMD's stream-order allocation API is beta and does not always work correctly, use synchronous function
        return cudaFree(devPtr);
    #endif
}

LIBMARV_NAMESPACE_END

#endif