// CUDA runtime subset used by libmarv's host code; implemented by each backend runtime (metal/cuda_metal_runtime.mm, simd/marv_simd_runtime.cpp).
#ifndef CUDA_RUNTIME_H_METAL_SHIM
#define CUDA_RUNTIME_H_METAL_SHIM

#include <cstddef>
#include <cstdint>
#include <cstring>

typedef struct CUstream_metal*  cudaStream_t;
typedef struct CUevent_metal*   cudaEvent_t;
typedef struct CUmemPool_metal* cudaMemPool_t;
typedef void (*cudaHostFn_t)(void* userData);

typedef enum cudaError {
    cudaSuccess = 0,
    cudaErrorInvalidValue = 1,
    cudaErrorMemoryAllocation = 2,
    cudaErrorInvalidDevice = 101,
    cudaErrorPeerAccessAlreadyEnabled = 704,
    cudaErrorPeerAccessNotEnabled = 705,
    cudaErrorNotSupported = 801,
    cudaErrorUnknown = 999
} cudaError_t;

typedef enum cudaMemcpyKind {
    cudaMemcpyHostToHost = 0,
    cudaMemcpyHostToDevice = 1,
    cudaMemcpyDeviceToHost = 2,
    cudaMemcpyDeviceToDevice = 3,
    cudaMemcpyDefault = 4
} cudaMemcpyKind;

typedef enum cudaDeviceAttr {
    cudaDevAttrPageableMemoryAccessUsesHostPageTables = 100
} cudaDeviceAttr;

typedef enum cudaMemPoolAttr {
    cudaMemPoolReuseFollowEventDependencies = 1,
    cudaMemPoolReuseAllowOpportunistic = 2,
    cudaMemPoolReuseAllowInternalDependencies = 3,
    cudaMemPoolAttrReleaseThreshold = 4,
    cudaMemPoolAttrReservedMemCurrent = 5,
    cudaMemPoolAttrReservedMemHigh = 6,
    cudaMemPoolAttrUsedMemCurrent = 7,
    cudaMemPoolAttrUsedMemHigh = 8
} cudaMemPoolAttr;

enum {
    cudaEventDefault = 0,
    cudaEventBlockingSync = 1,
    cudaEventDisableTiming = 2,
    cudaStreamDefault = 0,
    cudaStreamNonBlocking = 1,
    cudaHostRegisterDefault = 0,
    cudaIpcMemLazyEnablePeerAccess = 1};

// implicit streams
#define cudaStreamLegacy    ((cudaStream_t) 0x1)
#define cudaStreamPerThread ((cudaStream_t) 0x2)

typedef struct cudaIpcMemHandle_st { char reserved[64]; } cudaIpcMemHandle_t;

cudaError_t cudaMalloc(void** devPtr, size_t size);
cudaError_t cudaFree(void* devPtr);
cudaError_t cudaMallocHost(void** ptr, size_t size);
cudaError_t cudaFreeHost(void* ptr);
cudaError_t cudaMemcpy(void* dst, const void* src, size_t count, cudaMemcpyKind kind);
cudaError_t cudaMemset(void* devPtr, int value, size_t count);
cudaError_t cudaMemGetInfo(size_t* freeBytes, size_t* totalBytes);

cudaError_t cudaGetDevice(int* device);
cudaError_t cudaSetDevice(int device);
cudaError_t cudaGetDeviceCount(int* count);
cudaError_t cudaDeviceSynchronize(void);
cudaError_t cudaGetLastError(void);
cudaError_t cudaPeekAtLastError(void);
const char* cudaGetErrorString(cudaError_t error);

cudaError_t cudaStreamCreate(cudaStream_t* stream);
cudaError_t cudaStreamCreateWithFlags(cudaStream_t* stream, unsigned int flags);
cudaError_t cudaStreamDestroy(cudaStream_t stream);
cudaError_t cudaStreamSynchronize(cudaStream_t stream);
cudaError_t cudaStreamQuery(cudaStream_t stream);
cudaError_t cudaStreamWaitEvent(cudaStream_t stream, cudaEvent_t event, unsigned int flags = 0);

cudaError_t cudaEventCreate(cudaEvent_t* event);
cudaError_t cudaEventCreateWithFlags(cudaEvent_t* event, unsigned int flags);
cudaError_t cudaEventDestroy(cudaEvent_t event);
cudaError_t cudaEventRecord(cudaEvent_t event, cudaStream_t stream = 0);
cudaError_t cudaEventSynchronize(cudaEvent_t event);
cudaError_t cudaEventQuery(cudaEvent_t event);
cudaError_t cudaEventElapsedTime(float* ms, cudaEvent_t start, cudaEvent_t end);

cudaError_t cudaLaunchHostFunc(cudaStream_t stream, cudaHostFn_t fn, void* userData);

#include "../backend_hooks.h"

template<class T>
inline cudaError_t cudaMalloc(T** devPtr, size_t size) { return cudaMalloc(reinterpret_cast<void**>(devPtr), size); }
template<class T>
inline cudaError_t cudaMallocHost(T** ptr, size_t size) { return cudaMallocHost(reinterpret_cast<void**>(ptr), size); }

inline cudaError_t cudaMallocAsync(void** devPtr, size_t size, cudaStream_t = 0) { return cudaMalloc(devPtr, size); }
template<class T>
inline cudaError_t cudaMallocAsync(T** devPtr, size_t size, cudaStream_t = 0) { return cudaMalloc(reinterpret_cast<void**>(devPtr), size); }
inline cudaError_t cudaFreeAsync(void* devPtr, cudaStream_t = 0) { return cudaFree(devPtr); }

inline cudaError_t cudaMemcpyAsync(void* dst, const void* src, size_t count, cudaMemcpyKind kind, cudaStream_t = 0) {
    return cudaMemcpy(dst, src, count, kind);
}
inline cudaError_t cudaMemsetAsync(void* devPtr, int value, size_t count, cudaStream_t = 0) {
    return cudaMemset(devPtr, value, count);
}

inline cudaError_t cudaHostRegister(void*, size_t, unsigned int = 0) { return cudaSuccess; }

inline cudaError_t cudaDeviceGetAttribute(int* value, cudaDeviceAttr attr, int) {
    if (value == nullptr) return cudaErrorInvalidValue;
    // unified memory: lets GpuDatabase use the host database in place
    *value = (attr == cudaDevAttrPageableMemoryAccessUsesHostPageTables) ? 1 : 0;
    return cudaSuccess;
}

inline cudaError_t cudaDeviceGetDefaultMemPool(cudaMemPool_t* pool, int) { if (pool) *pool = reinterpret_cast<cudaMemPool_t>(1); return cudaSuccess; }
inline cudaError_t cudaMemPoolSetAttribute(cudaMemPool_t, cudaMemPoolAttr, void*) { return cudaSuccess; }
inline cudaError_t cudaMemPoolGetAttribute(cudaMemPool_t, cudaMemPoolAttr, void* value) { if (value) *static_cast<std::uint64_t*>(value) = 0; return cudaSuccess; }

inline cudaError_t cudaDeviceCanAccessPeer(int* canAccess, int, int) { if (canAccess) *canAccess = 0; return cudaSuccess; }
inline cudaError_t cudaDeviceEnablePeerAccess(int, unsigned int = 0) { return cudaSuccess; }
inline cudaError_t cudaDeviceDisablePeerAccess(int) { return cudaSuccess; }

inline cudaError_t cudaIpcGetMemHandle(cudaIpcMemHandle_t* handle, void*) { if (handle) memset(handle, 0, sizeof(*handle)); return cudaErrorNotSupported; }
inline cudaError_t cudaIpcOpenMemHandle(void** devPtr, cudaIpcMemHandle_t, unsigned int) { if (devPtr) *devPtr = nullptr; return cudaErrorNotSupported; }

#endif
