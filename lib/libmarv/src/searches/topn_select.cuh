#ifndef LIBMARV_TOPN_SELECT_CUH
#define LIBMARV_TOPN_SELECT_CUH

#include <thrust/sort.h>
#include <thrust/execution_policy.h>

#include "../util.cuh"
#include "../backend_hooks.h"
#if defined(__METAL_BACKEND__) || defined(__SIMD_BACKEND__)
#include "../shims/host_topn.h"
#endif

#include "../namespace.hpp"

LIBMARV_NAMESPACE_WITH_NESTING_BEGIN

// Orders the first count (score, id) pairs so that the best topN come first: descending score, ascending id.
template<class Id>
void selectTopN(int* scores, Id* ids, int count, int topN, cudaStream_t stream) {
#if defined(__METAL_BACKEND__) || defined(__SIMD_BACKEND__)
    if (topN <= 4096) {
        marvBackendSync();
        selectHostTopN(scores, ids, count, topN);
        return;
    }
#endif
    thrust::stable_sort_by_key(
        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
        scores, scores + count, ids, thrust::greater<int>()
    );
}

LIBMARV_NAMESPACE_WITH_NESTING_END

#endif
