// CUDA runtime facade for the SIMD (CPU) backend: device memory is host memory, streams and events are trivial.
#include <cuda_runtime.h>

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

namespace {
struct EventRecord { std::chrono::steady_clock::time_point time; bool recorded = false; };
cudaError_t g_lastError = cudaSuccess;
}

cudaError_t cudaMalloc(void** devPtr, size_t size) {
    if (devPtr == nullptr) return (g_lastError = cudaErrorInvalidValue);
    void* p = nullptr;
    if (posix_memalign(&p, 64, size ? size : 1) != 0) return (g_lastError = cudaErrorMemoryAllocation);
    *devPtr = p;
    return cudaSuccess;
}
cudaError_t cudaFree(void* devPtr) { free(devPtr); return cudaSuccess; }
cudaError_t cudaMallocHost(void** ptr, size_t size) { return cudaMalloc(ptr, size); }
cudaError_t cudaFreeHost(void* ptr) { return cudaFree(ptr); }
cudaError_t cudaMemcpy(void* dst, const void* src, size_t count, cudaMemcpyKind) { if (count) memmove(dst, src, count); return cudaSuccess; }
cudaError_t cudaMemset(void* devPtr, int value, size_t count) { if (count) memset(devPtr, value, count); return cudaSuccess; }
// CPU kernels read any host memory (also mapped files larger than RAM) in place: report an unbounded device
cudaError_t cudaMemGetInfo(size_t* freeBytes, size_t* totalBytes) {
    if (totalBytes) *totalBytes = SIZE_MAX / 2;
    if (freeBytes) *freeBytes = SIZE_MAX / 2;
    return cudaSuccess;
}

cudaError_t cudaGetDevice(int* device) { if (device) *device = 0; return cudaSuccess; }
cudaError_t cudaSetDevice(int device) { return device == 0 ? cudaSuccess : (g_lastError = cudaErrorInvalidDevice); }
cudaError_t cudaGetDeviceCount(int* count) { if (count) *count = 1; return cudaSuccess; }
cudaError_t cudaDeviceSynchronize(void) { return cudaSuccess; }
cudaError_t cudaGetLastError(void) { cudaError_t e = g_lastError; g_lastError = cudaSuccess; return e; }
cudaError_t cudaPeekAtLastError(void) { return g_lastError; }
const char* cudaGetErrorString(cudaError_t error) {
    switch (error) {
        case cudaSuccess: return "no error";
        case cudaErrorInvalidValue: return "invalid argument";
        case cudaErrorMemoryAllocation: return "out of memory";
        case cudaErrorInvalidDevice: return "invalid device ordinal";
        case cudaErrorNotSupported: return "operation not supported on the SIMD backend";
        default: return "unknown error";
    }
}

cudaError_t cudaStreamCreate(cudaStream_t* stream) { if (stream) *stream = reinterpret_cast<cudaStream_t>(new int(0)); return cudaSuccess; }
cudaError_t cudaStreamCreateWithFlags(cudaStream_t* stream, unsigned int) { return cudaStreamCreate(stream); }
cudaError_t cudaStreamDestroy(cudaStream_t stream) {
    if (stream != nullptr && stream != cudaStreamLegacy && stream != cudaStreamPerThread) delete reinterpret_cast<int*>(stream);
    return cudaSuccess;
}
cudaError_t cudaStreamSynchronize(cudaStream_t) { return cudaSuccess; }
cudaError_t cudaStreamQuery(cudaStream_t) { return cudaSuccess; }
cudaError_t cudaStreamWaitEvent(cudaStream_t, cudaEvent_t, unsigned int) { return cudaSuccess; }

cudaError_t cudaEventCreate(cudaEvent_t* event) { if (event) *event = reinterpret_cast<cudaEvent_t>(new EventRecord()); return cudaSuccess; }
cudaError_t cudaEventCreateWithFlags(cudaEvent_t* event, unsigned int) { return cudaEventCreate(event); }
cudaError_t cudaEventDestroy(cudaEvent_t event) { delete reinterpret_cast<EventRecord*>(event); return cudaSuccess; }
cudaError_t cudaEventRecord(cudaEvent_t event, cudaStream_t) {
    if (event == nullptr) return cudaErrorInvalidValue;
    EventRecord* e = reinterpret_cast<EventRecord*>(event);
    e->time = std::chrono::steady_clock::now(); e->recorded = true;
    return cudaSuccess;
}
cudaError_t cudaEventSynchronize(cudaEvent_t) { return cudaSuccess; }
cudaError_t cudaEventQuery(cudaEvent_t) { return cudaSuccess; }
cudaError_t cudaEventElapsedTime(float* ms, cudaEvent_t start, cudaEvent_t end) {
    if (ms == nullptr || start == nullptr || end == nullptr) return cudaErrorInvalidValue;
    EventRecord* a = reinterpret_cast<EventRecord*>(start);
    EventRecord* b = reinterpret_cast<EventRecord*>(end);
    if (!a->recorded || !b->recorded) return cudaErrorInvalidValue;
    *ms = (float) std::chrono::duration<double, std::milli>(b->time - a->time).count();
    return cudaSuccess;
}
cudaError_t cudaLaunchHostFunc(cudaStream_t, cudaHostFn_t fn, void* userData) { if (fn) fn(userData); return cudaSuccess; }

// backend hooks: nothing is asynchronous and all memory is host memory
void marvBackendSync(void) {}
bool marvBackendRegisterHostMemory(const void*, size_t) { return true; }
size_t marvBackendMaxBufferLength() { return SIZE_MAX; }
size_t marvBackendAllocationSize(size_t bytes) { return bytes; }
void marvBackendUnregisterHostMemory(const void*) {}
