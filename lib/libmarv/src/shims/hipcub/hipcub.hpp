// CPU cub::DeviceSelect::If. Namespace hipcub because the sources alias cub = hipcub on non-CUDA compilers.
#ifndef HIPCUB_HPP_METAL_SHIM
#define HIPCUB_HPP_METAL_SHIM

#include <cuda_runtime.h>
#include <cstddef>
#include <cstring>

namespace hipcub {

struct DeviceSelect {
    template<class InIt, class OutIt, class NumSelectedIt, class Pred>
    static cudaError_t If(void* d_temp_storage, size_t& temp_storage_bytes, InIt in, OutIt out, NumSelectedIt num_selected,
                          int num_items, Pred pred, cudaStream_t = 0) {
        // size query: report a non-zero size so that the second call is the execution
        if (d_temp_storage == nullptr) { temp_storage_bytes = 16; return cudaSuccess; }
        marvBackendSync();
        int count = 0;
        for (int i = 0; i < num_items; ++i) {
            auto v = in[i];
            if (pred(v)) { *out = v; ++out; ++count; }
        }
        *num_selected = count;
        return cudaSuccess;
    }
};

}

#endif
