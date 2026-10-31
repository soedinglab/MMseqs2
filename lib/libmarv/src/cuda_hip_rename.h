#if defined(__HIPCC__)

#ifndef __grid_constant__
#define __grid_constant__
#endif

#include <hip/hip_fp8.h>
using __nv_fp8_e4m3 = __hip_fp8_e4m3;
using __nv_fp8_e5m2 = __hip_fp8_e5m2;
using __nv_fp8x2_e4m3 = __hip_fp8x2_e4m3;

#define cudaStream_t hipStream_t
#define cudaEvent_t hipEvent_t
#define cudaError_t hipError_t
#define cudaMemcpyKind hipMemcpyKind
#define cudaMemPool_t hipMemPool_t
#define cudaIpcMemHandle_t hipIpcMemHandle_t
#define cudaDeviceProp hipDeviceProp_t

#define cudaSuccess hipSuccess
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyHostToHost hipMemcpyHostToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemcpyDefault hipMemcpyDefault
#define cudaEventDisableTiming hipEventDisableTiming
#define cudaStreamLegacy hipStreamLegacy
#define cudaStreamPerThread hipStreamPerThread
#define cudaStreamNonBlocking hipStreamNonBlocking
#define cudaErrorPeerAccessAlreadyEnabled hipErrorPeerAccessAlreadyEnabled
#define cudaErrorPeerAccessNotEnabled hipErrorPeerAccessNotEnabled
#define cudaIpcMemLazyEnablePeerAccess hipIpcMemLazyEnablePeerAccess
#define cudaDevAttrComputeCapabilityMajor hipDeviceAttributeComputeCapabilityMajor
#define cudaDevAttrComputeCapabilityMinor hipDeviceAttributeComputeCapabilityMinor
#define cudaDevAttrMultiProcessorCount hipDeviceAttributeMultiprocessorCount
#define cudaDevAttrMaxSharedMemoryPerBlockOptin hipDeviceAttributeSharedMemPerBlockOptin
#define cudaDevAttrPageableMemoryAccessUsesHostPageTables hipDeviceAttributePageableMemoryAccessUsesHostPageTables
#define cudaDeviceAttributeComputeCapabilityMajor hipDeviceAttributeComputeCapabilityMajor
#define cudaDeviceAttributeComputeCapabilityMinor hipDeviceAttributeComputeCapabilityMinor
#define cudaDeviceAttributeMultiprocessorCount hipDeviceAttributeMultiprocessorCount
#define cudaFuncAttributeMaxDynamicSharedMemorySize hipFuncAttributeMaxDynamicSharedMemorySize
#define cudaFuncCachePreferShared hipFuncCachePreferShared
#define cudaMemPoolAttrReleaseThreshold hipMemPoolAttrReleaseThreshold
#define cudaMemPoolAttrReservedMemCurrent hipMemPoolAttrReservedMemCurrent
#define cudaMemPoolAttrUsedMemCurrent hipMemPoolAttrUsedMemCurrent



#define cudaGetLastError hipGetLastError
#define cudaPeekAtLastError hipPeekAtLastError
#define cudaGetErrorString hipGetErrorString
#define cudaGetDevice hipGetDevice
#define cudaSetDevice hipSetDevice
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaDeviceGetAttribute hipDeviceGetAttribute
#define cudaDeviceCanAccessPeer hipDeviceCanAccessPeer
#define cudaDeviceEnablePeerAccess hipDeviceEnablePeerAccess
#define cudaDeviceDisablePeerAccess hipDeviceDisablePeerAccess
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaDeviceSetCacheConfig hipDeviceSetCacheConfig
#define cudaDeviceGetDefaultMemPool hipDeviceGetDefaultMemPool
#define cudaMemPoolSetAttribute hipMemPoolSetAttribute
#define cudaMemPoolTrimTo hipMemPoolTrimTo
#define cudaMemPoolGetAttribute hipMemPoolGetAttribute
#define cudaMemGetInfo hipMemGetInfo
#define cudaMalloc hipMalloc
#define cudaMallocAsync hipMallocAsync
#define cudaFree hipFree
#define cudaFreeAsync hipFreeAsync
#define cudaMallocHost hipHostMalloc
#define cudaFreeHost hipHostFree
#define cudaMemcpy hipMemcpy
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemsetAsync hipMemsetAsync
#define cudaLaunchHostFunc hipLaunchHostFunc
#define cudaFuncSetAttribute hipFuncSetAttribute
#define cudaOccupancyMaxActiveBlocksPerMultiprocessor hipOccupancyMaxActiveBlocksPerMultiprocessor
#define cudaEventCreate hipEventCreate
#define cudaEventCreateWithFlags hipEventCreateWithFlags
#define cudaEventDestroy hipEventDestroy
#define cudaEventQuery hipEventQuery
#define cudaEventRecord hipEventRecord
#define cudaEventSynchronize hipEventSynchronize
#define cudaEventElapsedTime hipEventElapsedTime
#define cudaStreamCreate hipStreamCreate
#define cudaStreamCreateWithFlags hipStreamCreateWithFlags
#define cudaStreamDestroy hipStreamDestroy
#define cudaStreamQuery hipStreamQuery
#define cudaStreamSynchronize hipStreamSynchronize
#define cudaStreamWaitEvent hipStreamWaitEvent
#define cudaIpcGetMemHandle hipIpcGetMemHandle
#define cudaIpcOpenMemHandle hipIpcOpenMemHandle

#include <thrust/system/hip/execution_policy.h>
namespace thrust {
namespace cuda = hip;
}

#ifndef LIBMARV_HIP_THRUST_NO_INIT
#define LIBMARV_HIP_THRUST_NO_INIT
namespace thrust {
    struct MarvNoInitTag {
        template<class T>
        constexpr operator T() const { return T{}; }
    };
    static constexpr MarvNoInitTag no_init{};
}
#endif

#ifndef LIBMARV_HIP_CG_ASYNC_SHIM
#define LIBMARV_HIP_CG_ASYNC_SHIM
namespace cuda {
    template<size_t Align>
    struct aligned_size_t {
        size_t value;
        __host__ __device__ constexpr aligned_size_t(size_t v) : value(v) {}
        __host__ __device__ constexpr operator size_t() const { return value; }
    };
}
namespace cooperative_groups {
    template<class Group, class DstT, class DstSize, class SrcT, class SrcSize>
    __device__ inline void memcpy_async(Group group, DstT dst, DstSize dstSize, SrcT src, SrcSize srcSize) {
        const size_t numBytes = size_t(dstSize) < size_t(srcSize) ? size_t(dstSize) : size_t(srcSize);
        char* d = reinterpret_cast<char*>(dst);
        const char* s = reinterpret_cast<const char*>(src);
        for (size_t i = group.thread_rank(); i < numBytes; i += group.num_threads()) {
            d[i] = s[i];
        }
    }
    template<class Group>
    __device__ inline void wait(Group group) { group.sync(); }
    template<int Stage, class Group>
    __device__ inline void wait_prior(Group group) { group.sync(); }
}
#endif

#endif