// Metal services of the runtime facade: kernel dispatch and pointer resolution.
#ifndef CUDA_METAL_LAUNCH_H
#define CUDA_METAL_LAUNCH_H

#include <cuda_runtime.h>
#include <cstddef>

struct MetalKernelArg {
    const void* ptr;
    size_t size;
    bool isBuffer;
    static MetalKernelArg buffer(const void* p) { return MetalKernelArg{p, 0, true}; }
    template<class T>
    static MetalKernelArg bytes(const T& v) { return MetalKernelArg{&v, sizeof(T), false}; }
};

cudaError_t metalLaunchKernel(const char* kernelName, const MetalKernelArg* args, int numArgs,
                              unsigned int gridThreadgroups, unsigned int threadsPerThreadgroup,
                              size_t threadgroupBytes, cudaStream_t stream);

size_t metalKernelMaxThreads(const char* kernelName);


bool metalResolvePointer(const void* ptr, size_t* byteOffsetOut);


size_t metalMaxThreadgroupMemory(void);

#endif
