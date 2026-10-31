// CUDA runtime facade for Metal: device memory is page-aligned host memory in shared MTLBuffers, one command queue,
// host access to device memory drains the pending command buffers, events are timestamps.
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <cuda_runtime.h>
#include "host_memory.h"
#include "cuda_metal_launch.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <map>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>
#include <unistd.h>

extern "C" {
extern unsigned char marv_metallib[];
extern unsigned int  marv_metallib_len;
}

namespace {

struct Region {
    uintptr_t base;
    size_t length;
    void* raw;
    id<MTLBuffer> buffer;
    size_t references = 1;
};

std::mutex g_mutex;
// owned allocations are disjoint, host registrations may overlap after page rounding
std::map<uintptr_t, Region> g_regions;
std::map<uintptr_t, Region> g_hostRegions;
id<MTLDevice> g_device = nil;
id<MTLCommandQueue> g_queue = nil;
id<MTLLibrary> g_library = nil;
bool g_libraryTried = false;
std::unordered_map<std::string, void*> g_pipelines;
std::vector<void*> g_pending;
std::mutex g_encodeMutex;
id<MTLCommandBuffer> g_recording = nil;
unsigned g_recordedDispatches = 0;

// called with g_encodeMutex held
void submitRecording() {
    if (g_recording == nil) return;
    @autoreleasepool {
        [g_recording commit];
        std::lock_guard<std::mutex> lk(g_mutex);
        g_pending.push_back((__bridge_retained void*) g_recording);
        g_recording = nil;
        g_recordedDispatches = 0;
    }
}
cudaError_t g_lastError = cudaSuccess;

struct EventRecord { std::chrono::steady_clock::time_point time; bool recorded = false; };

size_t pageSize() { return (size_t) getpagesize(); }
size_t pageRoundUp(size_t n) { const size_t pg = pageSize(); return ((n + pg - 1) / pg) * pg; }

void ensureDevice() {
    if (g_device != nil) return;
    g_device = MTLCreateSystemDefaultDevice();
    if (g_device == nil) {
        fprintf(stderr, "[marv-metal] no Metal device available\n");
        abort();
    }
    g_queue = [g_device newCommandQueue];
}

void ensureLibrary() {
    if (g_libraryTried) return;
    g_libraryTried = true;
    ensureDevice();
    NSError* err = nil;
    dispatch_data_t data = dispatch_data_create(marv_metallib, marv_metallib_len,
                                                dispatch_get_global_queue(QOS_CLASS_DEFAULT, 0),
                                                DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    g_library = [g_device newLibraryWithData:data error:&err];
    if (g_library == nil) {
        fprintf(stderr, "[marv-metal] failed to load the embedded kernel library: %s\n",
                err ? err.localizedDescription.UTF8String : "unknown error");
    }
}

Region* findRegion(const void* p) {
    const uintptr_t u = (uintptr_t) p;
    auto it = g_regions.upper_bound(u);
    if (it != g_regions.begin()) {
        Region& r = std::prev(it)->second;
        if (u - r.base < r.length) return &r;
    }
    // owned allocations first; a shorter overlapping registration must not hide an enclosing one
    auto host = g_hostRegions.upper_bound(u);
    while (host != g_hostRegions.begin()) {
        Region& r = (--host)->second;
        if (u >= r.base && u - r.base < r.length) return &r;
    }
    return nullptr;
}

cudaError_t allocRegistered(void** outPtr, size_t size) {
    if (outPtr == nullptr) return (g_lastError = cudaErrorInvalidValue);
    ensureDevice();
    if (size == 0) size = 1;
    if (size > size_t(g_device.maxBufferLength)) return (g_lastError = cudaErrorMemoryAllocation);
    const size_t bufLen = pageRoundUp(size);
    if (bufLen > size_t(g_device.maxBufferLength)) return (g_lastError = cudaErrorMemoryAllocation);
    void* raw = nullptr;
    if (posix_memalign(&raw, pageSize(), bufLen) != 0 || raw == nullptr) return (g_lastError = cudaErrorMemoryAllocation);
    id<MTLBuffer> buffer = [g_device newBufferWithBytesNoCopy:raw length:bufLen options:MTLResourceStorageModeShared deallocator:nil];
    if (buffer == nil) { free(raw); return (g_lastError = cudaErrorMemoryAllocation); }
    std::lock_guard<std::mutex> lk(g_mutex);
    g_regions[(uintptr_t) raw] = Region{(uintptr_t) raw, bufLen, raw, buffer};
    *outPtr = raw;
    return cudaSuccess;
}

void* pipelineForName(const char* name) {
    ensureLibrary();
    if (g_library == nil) return nullptr;
    std::string key(name);
    {
        std::lock_guard<std::mutex> lk(g_mutex);
        auto it = g_pipelines.find(key);
        if (it != g_pipelines.end()) return it->second;
    }
    id<MTLFunction> fn = [g_library newFunctionWithName:[NSString stringWithUTF8String:name]];
    if (fn == nil) return nullptr;
    NSError* err = nil;
    id<MTLComputePipelineState> pso = [g_device newComputePipelineStateWithFunction:fn error:&err];
    if (pso == nil) {
        fprintf(stderr, "[marv-metal] failed to build pipeline '%s': %s\n", name, err ? err.localizedDescription.UTF8String : "unknown");
        return nullptr;
    }
    void* retained = (__bridge_retained void*) pso;
    std::lock_guard<std::mutex> lk(g_mutex);
    g_pipelines[key] = retained;
    return retained;
}

}

cudaError_t cudaMalloc(void** devPtr, size_t size) { return allocRegistered(devPtr, size); }
cudaError_t cudaMallocHost(void** ptr, size_t size) { return allocRegistered(ptr, size); }

cudaError_t cudaFree(void* devPtr) {
    if (devPtr == nullptr) return cudaSuccess;
    marvBackendSync();
    std::lock_guard<std::mutex> lk(g_mutex);
    auto it = g_regions.find((uintptr_t) devPtr);
    if (it == g_regions.end() || it->second.raw == nullptr) return (g_lastError = cudaErrorInvalidValue);
    void* raw = it->second.raw;
    g_regions.erase(it); // release the no-copy buffer before its backing storage
    free(raw);
    return cudaSuccess;
}
cudaError_t cudaFreeHost(void* ptr) { return cudaFree(ptr); }

// every host-side access drains the queue first: this gives the stream ordering the host code relies on
cudaError_t cudaMemcpy(void* dst, const void* src, size_t count, cudaMemcpyKind) {
    marvBackendSync();
    if (count) memmove(dst, src, count);
    return cudaSuccess;
}
cudaError_t cudaMemset(void* devPtr, int value, size_t count) {
    marvBackendSync();
    if (count) memset(devPtr, value, count);
    return cudaSuccess;
}
cudaError_t cudaMemGetInfo(size_t* freeBytes, size_t* totalBytes) {
    ensureDevice();
    size_t total = (size_t) [g_device recommendedMaxWorkingSetSize];
    if (total == 0) total = (size_t) 8 << 30;
    if (totalBytes) *totalBytes = total;
    const size_t allocated = size_t(g_device.currentAllocatedSize);
    if (freeBytes) *freeBytes = std::min(marvAvailableHostMemory(total),
                                       allocated < total ? total - allocated : size_t(0));
    return cudaSuccess;
}

cudaError_t cudaGetDevice(int* device) { if (device) *device = 0; return cudaSuccess; }
cudaError_t cudaSetDevice(int device) { return device == 0 ? cudaSuccess : (g_lastError = cudaErrorInvalidDevice); }
cudaError_t cudaGetDeviceCount(int* count) {
    if (count == nullptr) return cudaErrorInvalidValue;
    @autoreleasepool { *count = (MTLCreateSystemDefaultDevice() != nil) ? 1 : 0; }
    return cudaSuccess;
}
cudaError_t cudaDeviceSynchronize(void) { marvBackendSync(); return cudaSuccess; }
cudaError_t cudaGetLastError(void) { cudaError_t e = g_lastError; g_lastError = cudaSuccess; return e; }
cudaError_t cudaPeekAtLastError(void) { return g_lastError; }
const char* cudaGetErrorString(cudaError_t error) {
    switch (error) {
        case cudaSuccess: return "no error";
        case cudaErrorInvalidValue: return "invalid argument";
        case cudaErrorMemoryAllocation: return "out of memory";
        case cudaErrorInvalidDevice: return "invalid device ordinal";
        case cudaErrorPeerAccessAlreadyEnabled: return "peer access is already enabled";
        case cudaErrorPeerAccessNotEnabled: return "peer access has not been enabled";
        case cudaErrorNotSupported: return "operation not supported on the Metal backend";
        default: return "unknown error";
    }
}

cudaError_t cudaStreamCreate(cudaStream_t* stream) { if (stream) *stream = reinterpret_cast<cudaStream_t>(new int(0)); return cudaSuccess; }
cudaError_t cudaStreamCreateWithFlags(cudaStream_t* stream, unsigned int) { return cudaStreamCreate(stream); }
cudaError_t cudaStreamDestroy(cudaStream_t stream) {
    if (stream != nullptr && stream != cudaStreamLegacy && stream != cudaStreamPerThread) delete reinterpret_cast<int*>(stream);
    return cudaSuccess;
}
cudaError_t cudaStreamSynchronize(cudaStream_t) { marvBackendSync(); return cudaSuccess; }
cudaError_t cudaStreamQuery(cudaStream_t) { marvBackendSync(); return cudaSuccess; }
cudaError_t cudaStreamWaitEvent(cudaStream_t, cudaEvent_t, unsigned int) { marvBackendSync(); return cudaSuccess; }

cudaError_t cudaEventCreate(cudaEvent_t* event) { if (event) *event = reinterpret_cast<cudaEvent_t>(new EventRecord()); return cudaSuccess; }
cudaError_t cudaEventCreateWithFlags(cudaEvent_t* event, unsigned int) { return cudaEventCreate(event); }
cudaError_t cudaEventDestroy(cudaEvent_t event) { delete reinterpret_cast<EventRecord*>(event); return cudaSuccess; }
cudaError_t cudaEventRecord(cudaEvent_t event, cudaStream_t) {
    if (event == nullptr) return cudaErrorInvalidValue;
    marvBackendSync();
    EventRecord* e = reinterpret_cast<EventRecord*>(event);
    e->time = std::chrono::steady_clock::now(); e->recorded = true;
    return cudaSuccess;
}
cudaError_t cudaEventSynchronize(cudaEvent_t) { marvBackendSync(); return cudaSuccess; }
cudaError_t cudaEventQuery(cudaEvent_t) { return cudaSuccess; }
cudaError_t cudaEventElapsedTime(float* ms, cudaEvent_t start, cudaEvent_t end) {
    if (ms == nullptr || start == nullptr || end == nullptr) return cudaErrorInvalidValue;
    EventRecord* a = reinterpret_cast<EventRecord*>(start);
    EventRecord* b = reinterpret_cast<EventRecord*>(end);
    if (!a->recorded || !b->recorded) return cudaErrorInvalidValue;
    *ms = (float) std::chrono::duration<double, std::milli>(b->time - a->time).count();
    return cudaSuccess;
}

cudaError_t cudaLaunchHostFunc(cudaStream_t, cudaHostFn_t fn, void* userData) {
    marvBackendSync();
    if (fn) fn(userData);
    return cudaSuccess;
}

size_t marvBackendMaxBufferLength() { ensureDevice(); return size_t(g_device.maxBufferLength); }
size_t marvBackendAllocationSize(size_t bytes) { return bytes ? pageRoundUp(bytes) : 0; }

bool marvBackendRegisterHostMemory(const void* ptr, size_t bytes) {
    if (bytes == 0) return true;
    if (ptr == nullptr) return false;
    ensureDevice();
    const size_t pg = pageSize();
    const uintptr_t u = (uintptr_t) ptr;
    const uintptr_t pageStart = u & ~(uintptr_t)(pg - 1);
    const size_t maximum = size_t(g_device.maxBufferLength);
    if (bytes > maximum || u - pageStart > maximum - bytes) return false;
    const size_t len = pageRoundUp((u - pageStart) + bytes);
    if (len > maximum) return false;
    std::lock_guard<std::mutex> lk(g_mutex);
    auto existing = g_hostRegions.find(u);
    if (existing != g_hostRegions.end() && len <= existing->second.length) {
        ++existing->second.references;
        return true;
    }
    id<MTLBuffer> buffer = [g_device newBufferWithBytesNoCopy:(void*) pageStart length:len options:MTLResourceStorageModeShared deallocator:nil];
    if (buffer == nil) return false;
    const size_t references = existing == g_hostRegions.end() ? 1 : existing->second.references + 1;
    g_hostRegions[u] = Region{pageStart, len, nullptr, buffer, references};
    return true;
}

void marvBackendUnregisterHostMemory(const void* ptr) {
    if (ptr == nullptr) return;
    marvBackendSync();
    std::lock_guard<std::mutex> lk(g_mutex);
    auto it = g_hostRegions.find((uintptr_t) ptr);
    if (it != g_hostRegions.end() && --it->second.references == 0) g_hostRegions.erase(it);
}

bool metalResolvePointer(const void* ptr, size_t* byteOffsetOut) {
    std::lock_guard<std::mutex> lk(g_mutex);
    Region* r = findRegion(ptr);
    if (r == nullptr) return false;
    if (byteOffsetOut) *byteOffsetOut = (uintptr_t) ptr - r->base;
    return true;
}

size_t metalKernelMaxThreads(const char* name) {
    void* p = pipelineForName(name);
    if (p == nullptr) return 0;
    id<MTLComputePipelineState> pso = (__bridge id<MTLComputePipelineState>) p;
    return (size_t) [pso maxTotalThreadsPerThreadgroup];
}

size_t metalMaxThreadgroupMemory(void) {
    ensureDevice();
    return (size_t) [g_device maxThreadgroupMemoryLength];
}

cudaError_t metalLaunchKernel(const char* kernelName, const MetalKernelArg* args, int numArgs,
                              unsigned int gridThreadgroups, unsigned int threadsPerThreadgroup,
                              size_t threadgroupBytes, cudaStream_t) {
    ensureDevice();
    void* psoPtr = pipelineForName(kernelName);
    if (psoPtr == nullptr) {
        fprintf(stderr, "[marv-metal] kernel '%s' not found in the embedded library\n", kernelName);
        return (g_lastError = cudaErrorNotSupported);
    }
    id<MTLComputePipelineState> pso = (__bridge id<MTLComputePipelineState>) psoPtr;
    @autoreleasepool {
        std::lock_guard<std::mutex> encodeLock(g_encodeMutex);
        if (g_recording == nil) g_recording = [g_queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [g_recording computeCommandEncoder];
        [enc setComputePipelineState:pso];
        for (int i = 0; i < numArgs; ++i) {
            if (args[i].isBuffer) {
                id<MTLBuffer> buffer = nil;
                {
                    std::lock_guard<std::mutex> lk(g_mutex);
                    Region* r = findRegion(args[i].ptr);
                    if (r) buffer = r->buffer;
                }
                if (buffer == nil) {
                    [enc endEncoding];
                    fprintf(stderr, "[marv-metal] launch '%s': argument %d (%p) is not device or registered memory\n", kernelName, i, args[i].ptr);
                    return (g_lastError = cudaErrorInvalidValue);
                }
                [enc setBuffer:buffer offset:0 atIndex:i];
            } else {
                [enc setBytes:args[i].ptr length:args[i].size atIndex:i];
            }
        }
        if (threadgroupBytes > 0) [enc setThreadgroupMemoryLength:threadgroupBytes atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(gridThreadgroups, 1, 1) threadsPerThreadgroup:MTLSizeMake(threadsPerThreadgroup, 1, 1)];
        [enc endEncoding];
        if (++g_recordedDispatches >= 16) submitRecording();
    }
    return cudaSuccess;
}

void marvBackendSync(void) {
    std::vector<void*> pending;
    {
        std::lock_guard<std::mutex> encodeLock(g_encodeMutex);
        submitRecording();
        std::lock_guard<std::mutex> lk(g_mutex);
        pending.swap(g_pending);
    }
    for (void* p : pending) {
        id<MTLCommandBuffer> cb = (__bridge_transfer id<MTLCommandBuffer>) p;
        [cb waitUntilCompleted];
        if (cb.status == MTLCommandBufferStatusError) {
            fprintf(stderr, "[marv-metal] command buffer error: %s\n", cb.error ? cb.error.localizedDescription.UTF8String : "unknown");
            g_lastError = cudaErrorUnknown;
        }
    }
}
