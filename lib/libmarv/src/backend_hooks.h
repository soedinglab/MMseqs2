#ifndef LIBMARV_BACKEND_HOOKS_H
#define LIBMARV_BACKEND_HOOKS_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

// Host backends (Metal, SIMD) implement these in their runtime; CUDA/HIP need none of them.
#if defined(__CUDACC__) || defined(__HIPCC__)
inline void marvBackendSync(void) {}
inline bool marvBackendRegisterHostMemory(const void*, size_t) { return true; }
inline void marvBackendUnregisterHostMemory(const void*) {}
inline size_t marvBackendMaxBufferLength() { return SIZE_MAX; }
inline size_t marvBackendAllocationSize(size_t bytes) { return bytes; }
#else
void marvBackendSync(void);
bool marvBackendRegisterHostMemory(const void* ptr, size_t bytes);
void marvBackendUnregisterHostMemory(const void* ptr);
size_t marvBackendMaxBufferLength();
size_t marvBackendAllocationSize(size_t bytes);
#endif

namespace libmarv {
// Registers the host database arrays with the backend as one unit; empty on failure.
template<class Data>
std::shared_ptr<void> registerHostDatabase(const Data& data) {
    auto registrations = std::shared_ptr<std::vector<const void*>>(
        new std::vector<const void*>, [](auto* pointers) {
            for (const void* p : *pointers) marvBackendUnregisterHostMemory(p);
            delete pointers;
        });
    const std::pair<const void*, size_t> arrays[] = {
        {data.chars(), data.numChars()},
        {data.lengths(), sizeof(*data.lengths()) * data.numSequences()},
        {data.offsets(), sizeof(*data.offsets()) * (data.numSequences() + 1)}
    };
    registrations->reserve(3);
    for (const auto& array : arrays) {
        if (!marvBackendRegisterHostMemory(array.first, array.second)) return {};
        if (array.second) registrations->push_back(array.first);
    }
    return registrations;
}
} // namespace libmarv

#endif
