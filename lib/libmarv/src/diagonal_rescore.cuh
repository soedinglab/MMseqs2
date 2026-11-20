#ifndef LIBMARV_DIAGONAL_RESCORE_CUH
#define LIBMARV_DIAGONAL_RESCORE_CUH

#include <cuda/cmath>

#include <cstdint>
#if defined(__CUDACC__)
#include <cuda_fp16.h>
#else
#include <hip/hip_fp16.h>
#endif

#if defined(__CUDACC__)
#include <cooperative_groups.h>
#else
#include <hip/hip_cooperative_groups.h>
#endif

#include "mathops.cuh"
#include "cuda_errorcheck.cuh"
#include "namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN


__global__
void diagonal_rescore_simple_kernel(
    __grid_constant__ int* const __restrict__ updatedScoreOutput, //result[i] will be added to output[i]
    __grid_constant__ const std::int8_t* const __restrict__ subjectData, // must be in [0,11]
    __grid_constant__ const int* const __restrict__ subjectLengths,
    __grid_constant__ const size_t* const __restrict__ subjectOffsets,
    __grid_constant__ const int numSubjects,
    __grid_constant__ const int* const __restrict__ subjectEndPositionsExcl,
    __grid_constant__ const int* const __restrict__ queryEndPositionsExcl,
    __grid_constant__ const half* const pssm_aa12, //12+1 rows, queryLength columns. last row is filled with 0
    __grid_constant__ const int queryLength
){
    //1 warp per subject
    constexpr int warpsize = 32;
    auto gridgroup = cooperative_groups::this_grid();
    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<warpsize>(block);

    const int warpIdInGrid = gridgroup.thread_rank() / warp.num_threads();
    const int numWarpsInGrid = gridgroup.num_threads() / warp.num_threads();

    for(int alignmentId = warpIdInGrid; alignmentId < numSubjects; alignmentId += numWarpsInGrid){
        const auto* subjectBasePtr = subjectData + subjectOffsets[alignmentId];
        const int subjectLength = subjectLengths[alignmentId];
        const int subjectEndExcl = subjectEndPositionsExcl[alignmentId];
        const int queryEndExcl = queryEndPositionsExcl[alignmentId];
        //query oriented along x-axis, subject along y-axis
        const int bestDiagonal = queryEndExcl - subjectEndExcl;
        const int currentAlignmentScore = updatedScoreOutput[alignmentId];

        const int start_Q = max(0, bestDiagonal);
        const int start_S = max(0, -bestDiagonal);
        const int max_end_Q = (bestDiagonal >= 0) ? queryLength : min(queryLength, start_Q + (subjectLength - start_S));
        const int max_end_S = (bestDiagonal < 0) ? subjectLength : min(subjectLength, start_S + (queryLength - start_Q));
        const int cellsOnDiagonal = min(max_end_Q - start_Q, max_end_S - start_S);
        // const int end_Q = start_Q + cellsOnDiagonal;
        const int end_S = start_S + cellsOnDiagonal;

        float seg_score = 0;
        float max12st_score = 0;
        for(int s = start_S, q = start_Q; s < end_S; s += warpsize){
            const int chunkSize = cuda::std::min(warpsize, end_S - s);
            int element = 12; //oob
            if(warp.thread_rank() < chunkSize){
                element = subjectBasePtr[s + warp.thread_rank()];
            }

            for(int i = 0; i < chunkSize; i++){
                int currentElement = warp.shfl(element, i);

                const half substscore = pssm_aa12[currentElement * queryLength + q];
                seg_score = MathOps<float>::add_relu(seg_score, substscore);
                max12st_score = max(max12st_score, seg_score);

                q++;
            }
        }

        if(warp.thread_rank() == 0){
            updatedScoreOutput[alignmentId] = currentAlignmentScore + max12st_score;
        }
    }

}

void call_diagonal_rescore_simple_kernel(
    int* const d_updatedScoreOutput, //result[i] will be added to output[i]
    const std::int8_t* const d_subjectData, // must be in [0,11]
    const int* const d_subjectLengths,
    const size_t* const d_subjectOffsets,
    const int numSubjects,
    const int* const d_subjectEndPositionsExcl,
    const int* const d_queryEndPositionsExcl,
    const half* const d_pssm_aa12, //12+1 rows, queryLength columns. last row is filled with 0
    const int queryLength,
    cudaStream_t stream
){
    const int blocksize = 128;
    const int numblocks = cuda::ceil_div(numSubjects, blocksize / 32);

    diagonal_rescore_simple_kernel<<<numblocks, blocksize, 0, stream>>>(
        d_updatedScoreOutput, //result[i] will be added to output[i]
        d_subjectData, // must be in [0,11]
        d_subjectLengths,
        d_subjectOffsets,
        numSubjects,
        d_subjectEndPositionsExcl,
        d_queryEndPositionsExcl,
        d_pssm_aa12, //12+1 rows, queryLength columns. last row is filled with 0
        queryLength
    );
    CUDACHECKASYNC;
}






LIBMARV_NAMESPACE_WITH_NESTING_END




#endif

