#ifndef LIBMARV_SEQUENCE_MASKING_CUH
#define LIBMARV_SEQUENCE_MASKING_CUH

//#define WITH_TIMING

#include "hpc_helpers/all_helpers.cuh"

#if defined(__CUDACC__)
#include <cub/cub.cuh>
#else
#include <hipcub/hipcub.hpp>
namespace cub = hipcub;
#endif

#include <thrust/device_vector.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/gather.h>
#include <thrust/binary_search.h>

#include <cuda/functional>
#include <cuda/cmath>
#include <cuda/std/type_traits>

#if defined(__CUDACC__)
#include <cuda/barrier>
#endif



#include "cuda_errorcheck.cuh"
// #include "util.cuh"
// #include "mathops.cuh"

#if defined(__CUDACC__)
#include <cooperative_groups.h>
#include <cooperative_groups/scan.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/memcpy_async.h>
#else
#include <hip/hip_cooperative_groups.h>
#endif

#include "cuda_hip_compatibility.cuh"

#include "namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN


namespace maskingsinglethread{

    template<typename T>
    constexpr bool is_power_of_two(T x){
        return (x > 0) && ((x & (x - 1)) == 0);
    }


    template<int numWarpsInBlock, int maxUintsPerThread>
    struct SingleTileKernel1Smem{
        uint32_t elements[numWarpsInBlock][maxUintsPerThread][32+1];        
        int segmentLengths[numWarpsInBlock][32];
        size_t segmentOffsets[numWarpsInBlock][32];
    };
    
    template<int blocksize, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    __global__
    void singleTileKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        __builtin_assume(blockDim.x == blocksize);
        static_assert(maxSegmentLength % 4 == 0);
        
        auto block = cooperative_groups::this_thread_block();
        auto warp = cooperative_groups::tiled_partition<32>(block);
        
        constexpr int numWarpsInBlock = blocksize / 32;
        const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
        const int warpIdInBlock = warp.meta_group_rank();
        const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;
        
        constexpr int maxUintsPerThread = maxSegmentLength / 4;
    
        using Smem = SingleTileKernel1Smem<numWarpsInBlock, maxUintsPerThread>;
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);

        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };
    
    
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("numSegments %d\n", numSegments);
        // }
    
        //each warp processes 32 sequences
        const int numBatches = cuda::ceil_div(numSegments, 32);
        for(int batch32 = warpIdInGrid; batch32 < numBatches; batch32 += numWarpsInGrid){
            const int workSegmentId = batch32 * 32 + warp.thread_rank();
            int workSegmentLength = 0;
            size_t workSegmentOffset = 0;
            if(workSegmentId < numSegments){
                workSegmentLength = d_lengths[workSegmentId];
                workSegmentOffset = d_offsets[workSegmentId];
            }
            const int workSegmentUints = cuda::ceil_div(workSegmentLength, 4);
            
            smem.segmentLengths[warpIdInBlock][warp.thread_rank()] = workSegmentLength;
            smem.segmentOffsets[warpIdInBlock][warp.thread_rank()] = workSegmentOffset;
            // if(batch32 == 0){
            //     printf("workSegmentId %d, workSegmentLength %d, workSegmentOffset %lu\n", 
            //         workSegmentId, workSegmentLength, workSegmentOffset);
            // }
            warp.sync();
    
                    
            {
               
    
                if constexpr(!is_power_of_two(maxUintsPerThread)){
                    uint32_t temploaded[32];

                    if(batch32 != numBatches-1){
                        for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                                const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                                const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                if(loadLane < segmentUints){
                                    temploaded[localLoadSegmentId] = segmentPtr[loadLane];
                                }
                            }
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] 
                                    = applyInputTransform(temploaded[localLoadSegmentId]);
                            }
                        }
                    }else{
                        const int numSegmentsInBatch = cuda::std::min(32, numSegments - batch32 * 32);
                        for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                if(localLoadSegmentId < numSegmentsInBatch){
                                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                    const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                    if(loadLane < segmentUints){
                                        temploaded[localLoadSegmentId] = segmentPtr[loadLane];
                                    }
                                }
                            }
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] 
                                    = applyInputTransform(temploaded[localLoadSegmentId]);
                            }
                        }
                    }
                }else{
                    uint32_t temploaded[maxUintsPerThread];
                    // #pragma unroll 32
                    for(int ui = 0; ui < maxUintsPerThread; ui++){
                        const int k = ui * warp.num_threads() + warp.thread_rank();
                        const int localLoadSegmentId = k / maxUintsPerThread;
                        const int loadLane = k % maxUintsPerThread;
        
                        // if(localLoadSegmentId < 32){
                            const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                            const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                            const int segmentUints = cuda::ceil_div(segmentLength, 4);
                            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                            if(loadLane < segmentUints){
                                temploaded[ui] = segmentPtr[loadLane];
                            }
                        // }
                    }
                    // #pragma unroll 32
                    for(int ui = 0; ui < maxUintsPerThread; ui++){
                        const int k = ui * warp.num_threads() + warp.thread_rank();
                        const int localLoadSegmentId = k / maxUintsPerThread;
                        const int loadLane = k % maxUintsPerThread;
        
                        // if(localLoadSegmentId < 32){
                            smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] 
                                = applyInputTransform(temploaded[ui]);
                        // }
                    }
                }
            }
    
            warp.sync();
    
    
            #if 1
            if(workSegmentId < numSegments){
                char prev = 0;
                int count = 0;
                int startOfRun = 0;
                int numRuns = 0;
                
    
                for(int uintindex = 0; uintindex < workSegmentUints; uintindex++){
                    uint32_t current4 = smem.elements[warpIdInBlock][uintindex][warp.thread_rank()];
                    alignas(4) char current4chars[4];
                    memcpy(&current4chars, &current4, 4);
    
                    // if(workSegmentId == 444){
                    //     printf("uintindex %d, load %2d %2d %2d %2d from %p\n", 
                    //         uintindex,
                    //         current4chars[0], current4chars[1], current4chars[2], current4chars[3],
                    //         (void*)&smem.elements[warpIdInBlock][uintindex][warp.thread_rank()]
                    //     );
                    // }
                    
                    if(uintindex == 0){
                        prev = current4chars[0];
                        //smem.runBeginPositions[warpIdInBlock][numRuns][warp.thread_rank()] = startOfRun;
                        numRuns++;
                    }
                    for(int posIn4 = 0; posIn4 < 4; posIn4++){
                        const int posInSegment = uintindex * 4 + posIn4;
                        if(posInSegment < workSegmentLength){
                            char cur = current4chars[posIn4];
                            // if(workSegmentId == 444){
                            //     printf("posInSegment %d, prev %d, prevcount %d, cur %d\n", posInSegment, prev, count, cur);
                            // }
                            if(cur == prev){
                                count++;
                            }else{
                                if(count >= maskingThreshold){
                                        // if(workSegmentId == 444){
                                        //     printf("masking loop %d - %d\n", startOfRun, posInSegment);
                                        // }
    
                                        #if 1
                                        for(int i = startOfRun; i < posInSegment; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
                                        #else
    
                                        const int nextborder = cuda::round_up(startOfRun, 4);
                                        const int lastborder = cuda::round_down(posInSegment, 4);
    
                                        for(int i = startOfRun; i < nextborder; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
    
                                        for(int i = nextborder; i < lastborder; i += 4){
                                            const int uindex2 = i / 4;
                                            ((char4*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[0] 
                                                = make_char4(maskingLetter, maskingLetter, maskingLetter, maskingLetter);
                                        }
    
                                        for(int i = lastborder; i < posInSegment; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
    
                                        #endif
                                }
                                prev = cur;
                                count = 1;
                                startOfRun = posInSegment;
    
                                //smem.runBeginPositions[warpIdInBlock][numRuns][warp.thread_rank()] = startOfRun;
                                numRuns++;
                            }
                        }
                    }
                }
                //handle last run
                if(count >= maskingThreshold){
    
                    // if(workSegmentId == 444){
                    //     printf("mask last run, startOfRun %d, workSegmentLength %d\n", startOfRun, workSegmentLength);
                    // }
                    for(int i = startOfRun; i < workSegmentLength; i++){
                        const int uindex2 = i / 4;
                        const int posui2 = i % 4;
                        ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                    }
                }
        
                if(applyExtraPaddingLetter){
                    for(int i = workSegmentLength; i < cuda::round_up(workSegmentLength,4); i++){
                        const int uindex2 = i / 4;
                        const int posui2 = i % 4;
                        ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = paddingLetter;
                    }
                }
            }
            #endif
    
            //wait until masking is done in smem
            warp.sync();
    
            //write smem back to gmem (inverse operation of loading )

            #if 1
            for(int ui = 0; ui < maxUintsPerThread; ui++){
                const int k = ui * warp.num_threads() + warp.thread_rank();
                const int localLoadSegmentId = k / maxUintsPerThread;
                const int loadLane = k % maxUintsPerThread;
    
                // if(localLoadSegmentId < 32){
                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);
                    if(loadLane < segmentUints){
                        segmentPtr[loadLane] = smem.elements[warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }
            #else

            const int numSegmentsInBatch = cuda::std::min(32, numSegments - batch32 * 32);
            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                // if(localLoadSegmentId < numSegmentsInBatch){
                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);

                    for(int loadLane = warp.thread_rank(); loadLane < segmentUints; loadLane += warp.num_threads()){
                        segmentPtr[loadLane] = smem.elements[warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }

            #endif
    
            //wait until smem can be re-used
            warp.sync();
    
        }
    }
    
    
    template<int blocksize, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    void callSingleTileKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));

        int maxBlocksPerSM = 0;
    
        constexpr int maxUintsPerThread = maxSegmentLength / 4;
        constexpr int numWarpsInBlock = blocksize / 32;
        auto kernel = singleTileKernel<blocksize, maxSegmentLength, InputOffsetIterator, InputTransformOp>;
        using SmemType = SingleTileKernel1Smem<numWarpsInBlock, maxUintsPerThread>;
    
        int smembytes = sizeof(SmemType);
    
        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }
    
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));
        
        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";
    
        #ifdef WITH_TIMING
        helpers::GpuTimer timer1("maskingsinglethread::singleTileKernel<" + std::to_string(blocksize) + "," + std::to_string(maxSegmentLength) + ">");
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );

        #ifdef WITH_TIMING
        timer1.print();
        #endif
    }


    template<int numWarpsInBlock, int maxUintsPerThread>
    struct SingleTileWithPlainFlagsKernel1Smem{
        uint32_t elements[numWarpsInBlock][maxUintsPerThread][32+1];        
        int segmentLengths[numWarpsInBlock][32];
        size_t segmentOffsets[numWarpsInBlock][32];
    };

    template<int blocksize, int maxMaskingThreshold, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    __global__
    void singleTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        __builtin_assume(blockDim.x == blocksize);
        static_assert(maxSegmentLength % 4 == 0);
        
        auto block = cooperative_groups::this_thread_block();
        auto warp = cooperative_groups::tiled_partition<32>(block);
        
        constexpr int numWarpsInBlock = blocksize / 32;
        const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
        const int warpIdInBlock = warp.meta_group_rank();
        const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;
        
        constexpr int maxUintsPerThread = maxSegmentLength / 4;
        constexpr int maxMaskingHaloUints = cuda::ceil_div(maxMaskingThreshold, 4);

        const int maskingHaloUints = cuda::ceil_div(maskingThreshold, 4);
        __builtin_assume(maskingThreshold <= maxMaskingThreshold);
        __builtin_assume(maskingHaloUints <= maxMaskingHaloUints);
    
        using Smem = SingleTileWithPlainFlagsKernel1Smem<numWarpsInBlock, maxUintsPerThread>;
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);


        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };
    
    
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("numSegments %d\n", numSegments);
        // }
    
        //each warp processes 32 sequences
        const int numBatches = cuda::ceil_div(numSegments, 32);
        for(int batch32 = warpIdInGrid; batch32 < numBatches; batch32 += numWarpsInGrid){
            const int workSegmentId = batch32 * 32 + warp.thread_rank();
            int workSegmentLength = 0;
            size_t workSegmentOffset = 0;
            if(workSegmentId < numSegments){
                workSegmentLength = d_lengths[workSegmentId];
                workSegmentOffset = d_offsets[workSegmentId];
            }
            const int workSegmentUints = cuda::ceil_div(workSegmentLength, 4);
            
            smem.segmentLengths[warpIdInBlock][warp.thread_rank()] = workSegmentLength;
            smem.segmentOffsets[warpIdInBlock][warp.thread_rank()] = workSegmentOffset;
            // if(batch32 == 0){
            //     printf("workSegmentId %d, workSegmentLength %d, workSegmentOffset %lu\n", 
            //         workSegmentId, workSegmentLength, workSegmentOffset);
            // }
            warp.sync();
    
                    
            {
               
    
                if constexpr(!is_power_of_two(maxUintsPerThread)){
                    uint32_t temploaded[32];

                    if(batch32 != numBatches-1){
                        for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                                const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                                const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                if(loadLane < segmentUints){
                                    temploaded[localLoadSegmentId] = segmentPtr[loadLane];
                                }
                            }
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] 
                                    = applyInputTransform(temploaded[localLoadSegmentId]);
                            }
                        }
                    }else{
                        const int numSegmentsInBatch = cuda::std::min(32, numSegments - batch32 * 32);
                        for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                if(localLoadSegmentId < numSegmentsInBatch){
                                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                    const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                    if(loadLane < segmentUints){
                                        temploaded[localLoadSegmentId] = segmentPtr[loadLane];
                                    }
                                }
                            }
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] 
                                    = applyInputTransform(temploaded[localLoadSegmentId]);
                            }
                        }
                    }
                }else{
                    uint32_t temploaded[maxUintsPerThread];
                    // #pragma unroll 32
                    for(int ui = 0; ui < maxUintsPerThread; ui++){
                        const int k = ui * warp.num_threads() + warp.thread_rank();
                        const int localLoadSegmentId = k / maxUintsPerThread;
                        const int loadLane = k % maxUintsPerThread;
        
                        // if(localLoadSegmentId < 32){
                            const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                            const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                            const int segmentUints = cuda::ceil_div(segmentLength, 4);
                            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                            if(loadLane < segmentUints){
                                temploaded[ui] = segmentPtr[loadLane];
                            }
                        // }
                    }
                    // #pragma unroll 32
                    for(int ui = 0; ui < maxUintsPerThread; ui++){
                        const int k = ui * warp.num_threads() + warp.thread_rank();
                        const int localLoadSegmentId = k / maxUintsPerThread;
                        const int loadLane = k % maxUintsPerThread;
        
                        // if(localLoadSegmentId < 32){
                            smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] 
                                = applyInputTransform(temploaded[ui]);
                        // }
                    }
                }
            }
    
            warp.sync();
    
    
            #if 1
            if(workSegmentId < numSegments){
                
                alignas(4) char flags[(maxUintsPerThread+maxMaskingHaloUints) * 4];
                alignas(4) char charElements[maxUintsPerThread * 4];

                //load elements
                #pragma unroll
                for(int i = 0; i < maxUintsPerThread; i++){
                    if(i < workSegmentUints){
                        uint32_t tmp = smem.elements[warpIdInBlock][i][warp.thread_rank()];
                        memcpy(&charElements[4*i], &tmp, 4);
                    }
                }

                //compute flags
                {
                    #pragma unroll
                    for(int i = 0; i < (maxUintsPerThread+maxMaskingHaloUints)*4; i++){
                        if(i == 0){
                            flags[i] = 1;
                        }else{
                            if(i < workSegmentLength){
                                flags[i] = charElements[i-1] != charElements[i];
                            }else{
                                flags[i] = 1;
                            }
                        }
                    }
                }

                //counting
                int count = 0;
                int prevcount = 0;
                #pragma unroll
                for(int i = 0; i < maxUintsPerThread*4; i++){
                    if(flags[i] == 0){
                        count = prevcount;
                    }else{
                        //new segment begins at i, count to the right
                        count = 1;
                        #pragma unroll
                        for(int k = 0; k < maxMaskingThreshold-1; k++){
                            if(flags[i+1+k] == 1){
                                break;
                            }
                            count++;
                        }
                    }
                    if(count >= maskingThreshold){
                        charElements[i] = maskingLetter;
                    }
                    prevcount = count;
                }

                if(applyExtraPaddingLetter){
                    #pragma unroll
                    for(int i = 0; i < maxUintsPerThread; i++){
                        if(i >= workSegmentLength){
                            charElements[i] = paddingLetter;
                        }
                    }
                }

                //store updated elements to smem
                #pragma unroll
                for(int i = 0; i < maxUintsPerThread; i++){
                    if(i < workSegmentUints){
                        uint32_t tmp;
                        memcpy(&tmp, &charElements[4*i], 4);
                        smem.elements[warpIdInBlock][i][warp.thread_rank()] = tmp;
                    }
                }
            }

            #endif
    
            //wait until masking is done in smem
            warp.sync();
    
            //write smem back to gmem (inverse operation of loading )

            #if 1
            for(int ui = 0; ui < maxUintsPerThread; ui++){
                const int k = ui * warp.num_threads() + warp.thread_rank();
                const int localLoadSegmentId = k / maxUintsPerThread;
                const int loadLane = k % maxUintsPerThread;
    
                // if(localLoadSegmentId < 32){
                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);
                    if(loadLane < segmentUints){
                        segmentPtr[loadLane] = smem.elements[warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }
            #else

            const int numSegmentsInBatch = cuda::std::min(32, numSegments - batch32 * 32);
            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                // if(localLoadSegmentId < numSegmentsInBatch){
                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);

                    for(int loadLane = warp.thread_rank(); loadLane < segmentUints; loadLane += warp.num_threads()){
                        segmentPtr[loadLane] = smem.elements[warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }

            #endif
    
            //wait until smem can be re-used
            warp.sync();
    
        }
    }
    
    
    template<int blocksize, int maxMaskingThreshold, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    void callSingleTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        if(maskingThreshold > maxMaskingThreshold){
            throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
        }
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));

        int maxBlocksPerSM = 0;
    
        constexpr int maxUintsPerThread = maxSegmentLength / 4;
        constexpr int numWarpsInBlock = blocksize / 32;
        auto kernel = singleTileWithPlainFlagsKernel<blocksize, maxMaskingThreshold, maxSegmentLength, InputOffsetIterator, InputTransformOp>;
        using SmemType = SingleTileWithPlainFlagsKernel1Smem<numWarpsInBlock, maxUintsPerThread>;
    
        int smembytes = sizeof(SmemType);
    
        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }
    
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));
        
        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";
    
        #ifdef WITH_TIMING
        helpers::GpuTimer timer1("maskingsinglethread::singleTileWithPlainFlagsKernel<" + std::to_string(blocksize) + "," + std::to_string(maxSegmentLength) + ">");
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );
        #ifdef WITH_TIMING
        timer1.print();
        #endif
    }

#if defined(__CUDACC__)
    template<int blocksize, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    __global__
    void singleTileAsyncLoadKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.y == 1);
        __builtin_assume(blockDim.z == 1);

        static_assert(maxSegmentLength % 4 == 0);
        
        auto block = cooperative_groups::this_thread_block();
        auto warp = cooperative_groups::tiled_partition<32>(block);
        
        constexpr int numWarpsInBlock = blocksize / 32;
        const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
        const int warpIdInBlock = warp.meta_group_rank();
        const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;
        
        constexpr int maxUintsPerThread = maxSegmentLength / 4;

        #pragma nv_diag_suppress static_var_with_dynamic_init
        __shared__ cuda::barrier<cuda::thread_scope_block> barriers[numWarpsInBlock];
        auto& warpBarrier = barriers[warp.meta_group_rank()];
        if(warp.thread_rank() == 0){
            init(&warpBarrier, warp.num_threads());
        }
        warp.sync();

        using Smem = SingleTileKernel1Smem<numWarpsInBlock, maxUintsPerThread>;
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);


        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };
    

    
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("numSegments %d\n", numSegments);
        // }
    
        //each warp processes 32 sequences
        const int numBatches = cuda::ceil_div(numSegments, 32);
        for(int batch32 = warpIdInGrid; batch32 < numBatches; batch32 += numWarpsInGrid){
            const int workSegmentId = batch32 * 32 + warp.thread_rank();
            int workSegmentLength = 0;
            size_t workSegmentOffset = 0;
            if(workSegmentId < numSegments){
                workSegmentLength = d_lengths[workSegmentId];
                workSegmentOffset = d_offsets[workSegmentId];
            }
            const int workSegmentUints = cuda::ceil_div(workSegmentLength, 4);
            
            smem.segmentLengths[warpIdInBlock][warp.thread_rank()] = workSegmentLength;
            smem.segmentOffsets[warpIdInBlock][warp.thread_rank()] = workSegmentOffset;
            // if(batch32 == 0){
            //     printf("workSegmentId %d, workSegmentLength %d, workSegmentOffset %lu\n", 
            //         workSegmentId, workSegmentLength, workSegmentOffset);
            // }
            warp.sync();
    
                    
            {
               
    
                if constexpr(false && !is_power_of_two(maxUintsPerThread)){
                    if(batch32 != numBatches-1){
                        for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                                const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                                const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                if(loadLane < segmentUints){
                                    //smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] = segmentPtr[loadLane];
                                    cuda::memcpy_async(
                                        &smem.elements[warpIdInBlock][loadLane][localLoadSegmentId],
                                        &segmentPtr[loadLane],
                                        cuda::aligned_size_t<4>(sizeof(uint32_t)),
                                        warpBarrier
                                    );
                                }
                            }
                        }
                    }else{
                        const int numSegmentsInBatch = cuda::std::min(32, numSegments - batch32 * 32);
                        for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                if(localLoadSegmentId < numSegmentsInBatch){
                                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                    const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                    if(loadLane < segmentUints){
                                        //smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] = segmentPtr[loadLane];
                                        cuda::memcpy_async(
                                            &smem.elements[warpIdInBlock][loadLane][localLoadSegmentId],
                                            &segmentPtr[loadLane],
                                            cuda::aligned_size_t<4>(sizeof(uint32_t)),
                                            warpBarrier
                                        );
                                    }
                                }
                            }
                        }
                    }
                }else{
                    // #pragma unroll 32
                    for(int ui = 0; ui < maxUintsPerThread; ui++){
                        const int k = ui * warp.num_threads() + warp.thread_rank();
                        const int localLoadSegmentId = k / maxUintsPerThread;
                        const int loadLane = k % maxUintsPerThread;
        
                        // if(localLoadSegmentId < 32){
                            const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                            const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                            const int segmentUints = cuda::ceil_div(segmentLength, 4);
                            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                            if(loadLane < segmentUints){
                                //smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] = segmentPtr[loadLane];
                                cuda::memcpy_async(
                                    &smem.elements[warpIdInBlock][loadLane][localLoadSegmentId],
                                    &segmentPtr[loadLane],
                                    cuda::aligned_size_t<4>(sizeof(uint32_t)),
                                    warpBarrier
                                );
                            }
                        // }
                    }
                }
            }
    
            // warp.sync();
            warpBarrier.arrive_and_wait();
    
    
            #if 1
            if(workSegmentId < numSegments){
                char prev = 0;
                int count = 0;
                int startOfRun = 0;
                int numRuns = 0;
                
    
                for(int uintindex = 0; uintindex < workSegmentUints; uintindex++){
                    uint32_t current4 = smem.elements[warpIdInBlock][uintindex][warp.thread_rank()];
                    alignas(4) char current4chars[4];
                    memcpy(&current4chars, &current4, 4);
    
                    // if(workSegmentId == 444){
                    //     printf("uintindex %d, load %2d %2d %2d %2d from %p\n", 
                    //         uintindex,
                    //         current4chars[0], current4chars[1], current4chars[2], current4chars[3],
                    //         (void*)&smem.elements[warpIdInBlock][uintindex][warp.thread_rank()]
                    //     );
                    // }
                    
                    if(uintindex == 0){
                        prev = current4chars[0];
                        //smem.runBeginPositions[warpIdInBlock][numRuns][warp.thread_rank()] = startOfRun;
                        numRuns++;
                    }
                    for(int posIn4 = 0; posIn4 < 4; posIn4++){
                        const int posInSegment = uintindex * 4 + posIn4;
                        if(posInSegment < workSegmentLength){
                            char cur = current4chars[posIn4];
                            // if(workSegmentId == 444){
                            //     printf("posInSegment %d, prev %d, prevcount %d, cur %d\n", posInSegment, prev, count, cur);
                            // }
                            if(cur == prev){
                                count++;
                            }else{
                                if(count >= maskingThreshold){
                                        // if(workSegmentId == 444){
                                        //     printf("masking loop %d - %d\n", startOfRun, posInSegment);
                                        // }
    
                                        #if 1
                                        for(int i = startOfRun; i < posInSegment; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
                                        #else
    
                                        const int nextborder = cuda::round_up(startOfRun, 4);
                                        const int lastborder = cuda::round_down(posInSegment, 4);
    
                                        for(int i = startOfRun; i < nextborder; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
    
                                        for(int i = nextborder; i < lastborder; i += 4){
                                            const int uindex2 = i / 4;
                                            ((char4*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[0] 
                                                = make_char4(maskingLetter, maskingLetter, maskingLetter, maskingLetter);
                                        }
    
                                        for(int i = lastborder; i < posInSegment; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
    
                                        #endif
                                }
                                prev = cur;
                                count = 1;
                                startOfRun = posInSegment;
    
                                //smem.runBeginPositions[warpIdInBlock][numRuns][warp.thread_rank()] = startOfRun;
                                numRuns++;
                            }
                        }
                    }
                }
                //handle last run
                if(count >= maskingThreshold){
    
                    // if(workSegmentId == 444){
                    //     printf("mask last run, startOfRun %d, workSegmentLength %d\n", startOfRun, workSegmentLength);
                    // }
                    for(int i = startOfRun; i < workSegmentLength; i++){
                        const int uindex2 = i / 4;
                        const int posui2 = i % 4;
                        ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                    }
                }
        
                if(applyExtraPaddingLetter){
                    for(int i = workSegmentLength; i < cuda::round_up(workSegmentLength,4); i++){
                        const int uindex2 = i / 4;
                        const int posui2 = i % 4;
                        ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = paddingLetter;
                    }
                }
            }
            #endif
    
            //wait until masking is done in smem
            warp.sync();
    
            //write smem back to gmem (inverse operation of loading )

            #if 1
            for(int ui = 0; ui < maxUintsPerThread; ui++){
                const int k = ui * warp.num_threads() + warp.thread_rank();
                const int localLoadSegmentId = k / maxUintsPerThread;
                const int loadLane = k % maxUintsPerThread;
    
                // if(localLoadSegmentId < 32){
                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);
                    if(loadLane < segmentUints){
                        segmentPtr[loadLane] = smem.elements[warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }
            #else

            const int numSegmentsInBatch = cuda::std::min(32, numSegments - batch32 * 32);
            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                // if(localLoadSegmentId < numSegmentsInBatch){
                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);

                    for(int loadLane = warp.thread_rank(); loadLane < segmentUints; loadLane += warp.num_threads()){
                        segmentPtr[loadLane] = smem.elements[warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }

            #endif
    
            //wait until smem can be re-used
            warp.sync();
    
        }
    }
    
    
    template<int blocksize, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    void callSingleTileAsyncLoadKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));

        int maxBlocksPerSM = 0;
    
        constexpr int maxUintsPerThread = maxSegmentLength / 4;
        constexpr int numWarpsInBlock = blocksize / 32;
        auto kernel = singleTileAsyncLoadKernel<blocksize, maxSegmentLength, InputOffsetIterator, InputTransformOp>;
        using SmemType = SingleTileKernel1Smem<numWarpsInBlock, maxUintsPerThread>;
    
        int smembytes = sizeof(SmemType);
    
        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }
    
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));
        
        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";

        #ifdef WITH_TIMING
        helpers::GpuTimer timer1("maskingsinglethread::singleTileAsyncLoadKernel<" + std::to_string(blocksize) + "," + std::to_string(maxSegmentLength) + ">");
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );
        #ifdef WITH_TIMING
        timer1.print();
        #endif
    }
#endif // defined(__CUDACC__)




    template<int numWarpsInBlock, int smemMaxUintsPerTile, int maxMaskingThreshold>
    struct MultiTileKernelSmem{
        static constexpr int extraMaskingUints = cuda::ceil_div(maxMaskingThreshold, 4);

        uint32_t elements[numWarpsInBlock][smemMaxUintsPerTile + extraMaskingUints][32+1];
        int segmentLengths[numWarpsInBlock][32];
        size_t segmentOffsets[numWarpsInBlock][32];
    };

    template<int blocksize, int smemMaxPositionsPerTile, int maxMaskingThreshold, class InputOffsetIterator, class InputTransformOp>
    __global__
    void multiTileKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        __builtin_assume(blockDim.x == blocksize);
        static_assert(smemMaxPositionsPerTile % 4 == 0);
        
        auto block = cooperative_groups::this_thread_block();
        auto warp = cooperative_groups::tiled_partition<32>(block);
        
        constexpr int numWarpsInBlock = blocksize / 32;
        const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
        const int warpIdInBlock = warp.meta_group_rank();
        const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;
        
        static_assert(smemMaxPositionsPerTile % 4 == 0);
        constexpr int smemMaxUintsPerTile = smemMaxPositionsPerTile / 4;
        constexpr int extraMaskingUints = cuda::ceil_div(maxMaskingThreshold, 4);

        using Smem = MultiTileKernelSmem<numWarpsInBlock, smemMaxUintsPerTile, maxMaskingThreshold>;
        static_assert(Smem::extraMaskingUints == extraMaskingUints);
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);

        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };

        //each warp processes 32 sequences
        const int numBatches = cuda::ceil_div(numSegments, 32);
        for(int batch32 = warpIdInGrid; batch32 < numBatches; batch32 += numWarpsInGrid){
            const int workSegmentId = batch32 * 32 + warp.thread_rank();
            int workSegmentLength = 0;
            size_t workSegmentOffset = 0;
            if(workSegmentId < numSegments){
                workSegmentLength = d_lengths[workSegmentId];
                workSegmentOffset = d_offsets[workSegmentId];
            }
            const int workSegmentUints = cuda::ceil_div(workSegmentLength, 4);
            const int workSegmentNumTiles = cuda::ceil_div(workSegmentUints, smemMaxUintsPerTile);
            
            smem.segmentLengths[warpIdInBlock][warp.thread_rank()] = workSegmentLength;
            smem.segmentOffsets[warpIdInBlock][warp.thread_rank()] = workSegmentOffset;

            //wait for smem lengths and offsets
            warp.sync();

            int maximumSegmentLengthInBatch = compatibility::group_reduce_max(warp, workSegmentLength);

            //each segment/sequence is split into tiles of size smemMaxPositionsPerTile
            const int maxNumTiles = cuda::ceil_div(maximumSegmentLengthInBatch, smemMaxPositionsPerTile);
            const int numSegmentsInBatch = cuda::std::min(32, numSegments - batch32 * 32);

            //the current masking state
            char prev = 0;
            int count = 0;
            int startOfRunInSmem = 0;
            int numRuns = 0;
            int numMaskedRuns = 0;

            for(int tileNr = 0; tileNr < maxNumTiles; tileNr++){

                //for each sequence in batch, load current tile into smem. tiles are loaded collaboratively and stored transposed to smem
                {                        
                    #if 0
                        uint32_t tmp[32];                            
                        if(batch32 != numBatches-1){
                            for(int smemLoadLane = warp.thread_rank(); smemLoadLane < (extraMaskingUints + smemMaxUintsPerTile); smemLoadLane += warp.num_threads()){
                                const int globalLoadLane = tileNr * (smemMaxUintsPerTile) + smemLoadLane;
                                for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                    const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                    if(globalLoadLane < segmentUints){
                                        tmp[localLoadSegmentId] = segmentPtr[globalLoadLane];
                                    }
                                }
                                for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                    smem.elements[warpIdInBlock][smemLoadLane][localLoadSegmentId] 
                                        = applyInputTransform(tmp[localLoadSegmentId]);
                                }
                            }
                        }else{
                            for(int smemLoadLane = warp.thread_rank(); smemLoadLane < (extraMaskingUints + smemMaxUintsPerTile); smemLoadLane += warp.num_threads()){
                                const int globalLoadLane = tileNr * (smemMaxUintsPerTile) + smemLoadLane;
                                for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                    if(localLoadSegmentId < numSegmentsInBatch){
                                        const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                                        const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                                        const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                        const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                        if(globalLoadLane < segmentUints){
                                            tmp[localLoadSegmentId] = segmentPtr[globalLoadLane];
                                        }
                                    }
                                }
                                for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                    smem.elements[warpIdInBlock][smemLoadLane][localLoadSegmentId] 
                                        = applyInputTransform(tmp[localLoadSegmentId]);
                                }
                            }
                        }
                    #else
                        uint32_t tmp[(extraMaskingUints + smemMaxUintsPerTile)];
                        
                        for(int loadNr = 0; loadNr < (extraMaskingUints + smemMaxUintsPerTile); loadNr++){
                            const int k = loadNr * warp.num_threads() + warp.thread_rank();

                            const int localLoadSegmentId = k / (extraMaskingUints + smemMaxUintsPerTile);
                            const int smemLoadLane = (k % (extraMaskingUints + smemMaxUintsPerTile));
                            const int globalLoadLane = tileNr * (smemMaxUintsPerTile) + smemLoadLane;

                            const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                            const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                            const int segmentUints = cuda::ceil_div(segmentLength, 4);
                            const uint32_t* segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);

                            if(globalLoadLane < segmentUints){
                                tmp[loadNr] = segmentPtr[globalLoadLane];
                            }
                        }
                        for(int loadNr = 0; loadNr < (extraMaskingUints + smemMaxUintsPerTile); loadNr++){
                            const int k = loadNr * warp.num_threads() + warp.thread_rank();
                            const int localLoadSegmentId = k / (extraMaskingUints + smemMaxUintsPerTile);
                            const int smemLoadLane = (k % (extraMaskingUints + smemMaxUintsPerTile));

                            smem.elements[warpIdInBlock][smemLoadLane][localLoadSegmentId] 
                                = applyInputTransform(tmp[loadNr]);
                        }
                    #endif
                }

                //wait until tiles are loaded to smem
                warp.sync();

                #if 1
                if(workSegmentId < numSegments && tileNr < workSegmentNumTiles){
                    const int alreadyProcessedUints = tileNr * smemMaxUintsPerTile;
                    const int numUintsToProcess = cuda::std::min(smemMaxUintsPerTile, workSegmentUints - alreadyProcessedUints);

                    
                    for(int uintindex = 0; uintindex < numUintsToProcess; uintindex++){
                        uint32_t current4 = smem.elements[warpIdInBlock][uintindex][warp.thread_rank()];
                        alignas(4) char current4chars[4];
                        memcpy(&current4chars, &current4, 4);

                        if(tileNr == 0 && uintindex == 0){
                            prev = current4chars[0];
                            numRuns++;
                        }
                        for(int posIn4 = 0; posIn4 < 4; posIn4++){
                            const int posInSegmentSmem = uintindex * 4 + posIn4;
                            const int posInSegment = posInSegmentSmem + tileNr * smemMaxPositionsPerTile;
                            if(posInSegment < workSegmentLength){
                                char cur = current4chars[posIn4];
                                if(cur == prev){
                                    count++;
                                }else{
                                    if(count >= maskingThreshold){
                                        for(int i = startOfRunInSmem; i < posInSegmentSmem; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
                                        numMaskedRuns++;
                                    }
                                    prev = cur;
                                    count = 1;
                                    startOfRunInSmem = posInSegmentSmem;
        
                                    numRuns++;
                                }
                            }
                        }
                    }

                    // if we processed a tile which is not the last tile, 
                    // process the extra masking uints to handle overlap with next tile

                    if(tileNr < workSegmentNumTiles-1){
                        // if(workSegmentId == 842){
                        //     printf("tileNr %d, prev %d, count %d, startOfRunInSmem %d\n",
                        //         tileNr, prev, count, startOfRunInSmem
                        //     );
                        // }
                        const int count_backup = count;
                        bool done = false;
                        //check if the current run continues with a total length of at least maskingThreshold
                        for(int uintindex = numUintsToProcess; uintindex < numUintsToProcess + extraMaskingUints; uintindex++){
                            if(!done){
                                uint32_t current4 = smem.elements[warpIdInBlock][uintindex][warp.thread_rank()];
                                alignas(4) char current4chars[4];
                                memcpy(&current4chars, &current4, 4);

                                for(int posIn4 = 0; posIn4 < 4; posIn4++){
                                    // const int posInSegmentSmem = uintindex * 4 + posIn4;
                                    char cur = current4chars[posIn4];
                                    if(cur == prev){
                                        count++;
                                    }else{
                                        done = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if(count >= maskingThreshold){
                            // the run continued long enough, do mask
                            for(int i = startOfRunInSmem; i < 4*numUintsToProcess; i++){
                                const int uindex2 = i / 4;
                                const int posui2 = i % 4;
                                ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                            }
                            numMaskedRuns++;
                        }
                        //for next tile, run begins at begin of smem tile
                        startOfRunInSmem = 0;
                        count = count_backup;
                    }else{
                        //we processed the last tile of sequence.
                        const int numPositionsInLastTile = workSegmentLength - alreadyProcessedUints * 4;
                        //handle last run
                        if(count >= maskingThreshold){
                            
                            for(int i = startOfRunInSmem; i < numPositionsInLastTile; i++){
                                const int uindex2 = i / 4;
                                const int posui2 = i % 4;
                                ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                            }
                            numMaskedRuns++;
                        }
                        if(applyExtraPaddingLetter){
                            for(int i = numPositionsInLastTile; i < cuda::round_up(numPositionsInLastTile,4); i++){
                                const int uindex2 = i / 4;
                                const int posui2 = i % 4;
                                // assert(uindex2 < (extraMaskingUints + smemMaxUintsPerTile));
                                // assert(i < smemMaxPositionsPerTile);
                                // printf("uindex2 %d, posui2 %d, numPositionsInLastTile %d %d\n", 
                                //     uindex2, posui2, numPositionsInLastTile, cuda::round_up(numPositionsInLastTile,4));
                                ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = paddingLetter;
                            }
                        }
                        // //debug
                        // ((char*)&smem.elements[warpIdInBlock][numMaskedRuns / 4][warp.thread_rank()])[0] = paddingLetter;
                    }
                }
                
                //wait until tiles are processed in smem
                warp.sync();
                #endif

                //write smem tiles back to gmem
                {
                    #if 1
                        const int maxPositionsToLoadPerSequenceWithoutExtra = cuda::std::min(
                            smemMaxPositionsPerTile, 
                            maximumSegmentLengthInBatch - tileNr * smemMaxPositionsPerTile
                        );
                        const int maxUintsToLoadPerSequenceWithoutExtra = cuda::ceil_div(maxPositionsToLoadPerSequenceWithoutExtra, 4);
                        const int maxUintsToLoadPerSequence = extraMaskingUints + maxUintsToLoadPerSequenceWithoutExtra;
                        const int maxUintsToLoad = warp.num_threads() * maxUintsToLoadPerSequence;
                        
                        for(int k = warp.thread_rank(); k < maxUintsToLoad; k += warp.num_threads()){
                            const int localLoadSegmentId = k / maxUintsToLoadPerSequence;
                            const int smemLoadLane = (k % maxUintsToLoadPerSequence);
                            const int globalLoadLane = tileNr * smemMaxUintsPerTile + smemLoadLane;
                            // const int loadSegmentId = batch32 * 32 + localLoadSegmentId;
                            

                            const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                            const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                            const int segmentUints = cuda::ceil_div(segmentLength, 4);
                            uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);

                            if(globalLoadLane < segmentUints){
                                segmentPtr[globalLoadLane] = smem.elements[warpIdInBlock][smemLoadLane][localLoadSegmentId];
                            }
                        }
                    #else
                        for(int ui = 0; ui < smemMaxUintsPerTile; ui++){
                            const int k = ui * warp.num_threads() + warp.thread_rank();
                            const int localLoadSegmentId = k / smemMaxUintsPerTile;
                            const int smemLoadLane = k % smemMaxUintsPerTile;
                            const int globalLoadLane = tileNr * smemMaxUintsPerTile + smemLoadLane;
                
                            // if(localLoadSegmentId < 32){
                                const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                                const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                                const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);
                                if(globalLoadLane < segmentUints){
                                    segmentPtr[globalLoadLane] = smem.elements[warpIdInBlock][smemLoadLane][localLoadSegmentId];
                                }
                            // }
                        }
                    #endif
                }
                

                //wait until smem can be re-used for next chunk
                warp.sync();
            }
        }


    }

    template<int blocksize, int smemMaxPositionsPerTile, int maxMaskingThreshold, class InputOffsetIterator, class InputTransformOp>
    void callMultiTileKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        if(maskingThreshold > maxMaskingThreshold){
            throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
        }
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));

        int maxBlocksPerSM = 0;

        constexpr int smemMaxUintsPerTile = smemMaxPositionsPerTile / 4;
        constexpr int numWarpsInBlock = blocksize / 32;
        auto kernel = multiTileKernel<blocksize, smemMaxPositionsPerTile, maxMaskingThreshold, InputOffsetIterator, InputTransformOp>;
        using SmemType = MultiTileKernelSmem<numWarpsInBlock, smemMaxUintsPerTile, maxMaskingThreshold>;

        int smembytes = sizeof(SmemType);

        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }

        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));   

        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";

        #ifdef WITH_TIMING
        helpers::GpuTimer timer1("maskingsinglethread::multiTileKernel<" + std::to_string(blocksize) + "," + std::to_string(smemMaxPositionsPerTile) + ">");
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );
        #ifdef WITH_TIMING
        timer1.print();
        #endif
    }





    template<int numWarpsInBlock, int maxUintsPerThread, int numStages>
    struct SingleTileStagedKernel1Smem{
        uint32_t elements[1][numWarpsInBlock][maxUintsPerThread][32+1];
        // uint32_t elements[numWarpsInBlock][maxUintsPerThread][32+1];
        int segmentLengths[numStages][numWarpsInBlock][32];
        size_t segmentOffsets[numStages][numWarpsInBlock][32];
    };


    template<int blocksize, int maxSegmentLength, int numStages, class InputOffsetIterator, class InputTransformOp>
    __global__
    void singleTileStagedKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        
        __builtin_assume(blockDim.x == blocksize);
        static_assert(maxSegmentLength % 4 == 0);
        
        auto block = cooperative_groups::this_thread_block();
        auto warp = cooperative_groups::tiled_partition<32>(block);
        
        constexpr int numWarpsInBlock = blocksize / 32;
        const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
        const int warpIdInBlock = warp.meta_group_rank();
        const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;
        
        constexpr int maxUintsPerThread = maxSegmentLength / 4;
        //static_assert(maxUintsPerThread <= 32);
        if constexpr(maxUintsPerThread > 32){
            return;
        }
    
        using Smem = SingleTileStagedKernel1Smem<numWarpsInBlock, maxUintsPerThread, numStages>;
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);

        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };
    
    
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("numSegments %d\n", numSegments);
        // }
    
        //each warp processes 32 sequences
        const int numBatches = cuda::ceil_div(numSegments, 32);
        int stage = 0;

        
        uint32_t temploaded_current[32];
        uint32_t temploaded_next[32];

        constexpr int loadApproach = 1; //false && !is_power_of_two(maxUintsPerThread)

        //load data of first batch next data
        if(warpIdInGrid < numBatches){
            const int fetchbatch32 = warpIdInGrid;
            const int fetchSegmentId = fetchbatch32 * 32 + warp.thread_rank();
            const int fetchstage = stage;
            int fetchSegmentLength = 0;
            size_t fetchSegmentOffset = 0;
            if(fetchSegmentId < numSegments){
                fetchSegmentLength = d_lengths[fetchSegmentId];
                fetchSegmentOffset = d_offsets[fetchSegmentId];
            }
            smem.segmentLengths[fetchstage][warpIdInBlock][warp.thread_rank()] = fetchSegmentLength;
            smem.segmentOffsets[fetchstage][warpIdInBlock][warp.thread_rank()] = fetchSegmentOffset;
            warp.sync();

            if constexpr(loadApproach == 0){
                if(fetchbatch32 != numBatches-1){
                    for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                        for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                            const int segmentLength = smem.segmentLengths[fetchstage][warpIdInBlock][localLoadSegmentId];
                            const size_t segmentOffset = smem.segmentOffsets[fetchstage][warpIdInBlock][localLoadSegmentId];
                            const int segmentUints = cuda::ceil_div(segmentLength, 4);
                            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                            if(loadLane < segmentUints){
                                temploaded_next[localLoadSegmentId] = segmentPtr[loadLane];
                            }
                        }
                    }
                }else{
                    const int numSegmentsInBatch = cuda::std::min(32, numSegments - fetchbatch32 * 32);
                    for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                        for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                            if(localLoadSegmentId < numSegmentsInBatch){
                                const int segmentLength = smem.segmentLengths[fetchstage][warpIdInBlock][localLoadSegmentId];
                                const size_t segmentOffset = smem.segmentOffsets[fetchstage][warpIdInBlock][localLoadSegmentId];
                                const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                if(loadLane < segmentUints){
                                    temploaded_next[localLoadSegmentId] = segmentPtr[loadLane];
                                }
                            }
                        }
                    }
                }
            }else{
                for(int ui = 0; ui < maxUintsPerThread; ui++){
                    const int k = ui * warp.num_threads() + warp.thread_rank();
                    const int localLoadSegmentId = k / maxUintsPerThread;
                    const int loadLane = k % maxUintsPerThread;
    
                    // if(localLoadSegmentId < 32){
                        const int segmentLength = smem.segmentLengths[fetchstage][warpIdInBlock][localLoadSegmentId];
                        const size_t segmentOffset = smem.segmentOffsets[fetchstage][warpIdInBlock][localLoadSegmentId];
                        const int segmentUints = cuda::ceil_div(segmentLength, 4);
                        const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                        if(loadLane < segmentUints){
                            temploaded_next[ui] = segmentPtr[loadLane];
                        }
                    // }
                }
            }
        }


        for(int batch32 = warpIdInGrid; batch32 < numBatches; batch32 += numWarpsInGrid, stage = (stage + 1) % numStages){

            //make prefetched data the current data
            {
                for(int i = 0; i < 32; i++){
                    temploaded_current[i] = temploaded_next[i];
                }
            }

            //prefetch next data
            if(batch32 + numWarpsInGrid < numBatches){
                const int fetchbatch32 = batch32 + numWarpsInGrid;
                const int fetchSegmentId = fetchbatch32 * 32 + warp.thread_rank();
                const int fetchstage = (stage + 1) % 2;
                int fetchSegmentLength = 0;
                size_t fetchSegmentOffset = 0;
                if(fetchSegmentId < numSegments){
                    fetchSegmentLength = d_lengths[fetchSegmentId];
                    fetchSegmentOffset = d_offsets[fetchSegmentId];
                }
                smem.segmentLengths[fetchstage][warpIdInBlock][warp.thread_rank()] = fetchSegmentLength;
                smem.segmentOffsets[fetchstage][warpIdInBlock][warp.thread_rank()] = fetchSegmentOffset;
                warp.sync();

                if constexpr(loadApproach == 0){
                    if(fetchbatch32 != numBatches-1){
                        for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                const int segmentLength = smem.segmentLengths[fetchstage][warpIdInBlock][localLoadSegmentId];
                                const size_t segmentOffset = smem.segmentOffsets[fetchstage][warpIdInBlock][localLoadSegmentId];
                                const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                if(loadLane < segmentUints){
                                    temploaded_next[localLoadSegmentId] = segmentPtr[loadLane];
                                }
                            }
                        }
                    }else{
                        const int numSegmentsInBatch = cuda::std::min(32, numSegments - fetchbatch32 * 32);
                        for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                                if(localLoadSegmentId < numSegmentsInBatch){
                                    const int segmentLength = smem.segmentLengths[fetchstage][warpIdInBlock][localLoadSegmentId];
                                    const size_t segmentOffset = smem.segmentOffsets[fetchstage][warpIdInBlock][localLoadSegmentId];
                                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                                    const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                                    if(loadLane < segmentUints){
                                        temploaded_next[localLoadSegmentId] = segmentPtr[loadLane];
                                    }
                                }
                            }
                        }
                    }
                }else{
                    for(int ui = 0; ui < maxUintsPerThread; ui++){
                        const int k = ui * warp.num_threads() + warp.thread_rank();
                        const int localLoadSegmentId = k / maxUintsPerThread;
                        const int loadLane = k % maxUintsPerThread;
        
                        // if(localLoadSegmentId < 32){
                            const int segmentLength = smem.segmentLengths[fetchstage][warpIdInBlock][localLoadSegmentId];
                            const size_t segmentOffset = smem.segmentOffsets[fetchstage][warpIdInBlock][localLoadSegmentId];
                            const int segmentUints = cuda::ceil_div(segmentLength, 4);
                            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                            if(loadLane < segmentUints){
                                temploaded_next[ui] = segmentPtr[loadLane];
                            }
                        // }
                    }
                }
            }

            if constexpr(loadApproach == 0){
                for(int loadLane = warp.thread_rank(); loadLane < maxUintsPerThread; loadLane += warp.num_threads()){
                    for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                        smem.elements[0][warpIdInBlock][loadLane][localLoadSegmentId] 
                            = applyInputTransform(temploaded_current[localLoadSegmentId]);
                    }
                }
            }else{
                for(int ui = 0; ui < maxUintsPerThread; ui++){
                    const int k = ui * warp.num_threads() + warp.thread_rank();
                    const int localLoadSegmentId = k / maxUintsPerThread;
                    const int loadLane = k % maxUintsPerThread;
    
                    // if(localLoadSegmentId < 32){
                        smem.elements[0][warpIdInBlock][loadLane][localLoadSegmentId] 
                            = applyInputTransform(temploaded_current[ui]);
                    // }
                }
            }


            warp.sync();

            const int workSegmentId = batch32 * 32 + warp.thread_rank();
            int workSegmentLength = smem.segmentLengths[stage][warpIdInBlock][warp.thread_rank()];
            size_t workSegmentOffset = smem.segmentOffsets[stage][warpIdInBlock][warp.thread_rank()];
            const int workSegmentUints = cuda::ceil_div(workSegmentLength, 4);

    
            #if 1
            if(workSegmentId < numSegments){
                char prev = 0;
                int count = 0;
                int startOfRun = 0;
                int numRuns = 0;
                
    
                for(int uintindex = 0; uintindex < workSegmentUints; uintindex++){
                    uint32_t current4 = smem.elements[0][warpIdInBlock][uintindex][warp.thread_rank()];
                    alignas(4) char current4chars[4];
                    memcpy(&current4chars, &current4, 4);
    
                    // if(workSegmentId == 444){
                    //     printf("uintindex %d, load %2d %2d %2d %2d from %p\n", 
                    //         uintindex,
                    //         current4chars[0], current4chars[1], current4chars[2], current4chars[3],
                    //         (void*)&smem.elements[0][warpIdInBlock][uintindex][warp.thread_rank()]
                    //     );
                    // }
                    
                    if(uintindex == 0){
                        prev = current4chars[0];
                        numRuns++;
                    }
                    for(int posIn4 = 0; posIn4 < 4; posIn4++){
                        const int posInSegment = uintindex * 4 + posIn4;
                        if(posInSegment < workSegmentLength){
                            char cur = current4chars[posIn4];
                            // if(workSegmentId == 444){
                            //     printf("posInSegment %d, prev %d, prevcount %d, cur %d\n", posInSegment, prev, count, cur);
                            // }
                            if(cur == prev){
                                count++;
                            }else{
                                if(count >= maskingThreshold){
                                        // if(workSegmentId == 444){
                                        //     printf("masking loop %d - %d\n", startOfRun, posInSegment);
                                        // }
    
                                        #if 1
                                        for(int i = startOfRun; i < posInSegment; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[0][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
                                        #else
    
                                        const int nextborder = cuda::round_up(startOfRun, 4);
                                        const int lastborder = cuda::round_down(posInSegment, 4);
    
                                        for(int i = startOfRun; i < nextborder; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[0][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
    
                                        for(int i = nextborder; i < lastborder; i += 4){
                                            const int uindex2 = i / 4;
                                            ((char4*)&smem.elements[0][warpIdInBlock][uindex2][warp.thread_rank()])[0] 
                                                = make_char4(maskingLetter, maskingLetter, maskingLetter, maskingLetter);
                                        }
    
                                        for(int i = lastborder; i < posInSegment; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[0][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
    
                                        #endif
                                }
                                prev = cur;
                                count = 1;
                                startOfRun = posInSegment;
    
                                numRuns++;
                            }
                        }
                    }
                }
                //handle last run
                if(count >= maskingThreshold){
    
                    // if(workSegmentId == 444){
                    //     printf("mask last run, startOfRun %d, workSegmentLength %d\n", startOfRun, workSegmentLength);
                    // }
                    for(int i = startOfRun; i < workSegmentLength; i++){
                        const int uindex2 = i / 4;
                        const int posui2 = i % 4;
                        ((char*)&smem.elements[0][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                    }
                }
        
                if(applyExtraPaddingLetter){
                    for(int i = workSegmentLength; i < cuda::round_up(workSegmentLength,4); i++){
                        const int uindex2 = i / 4;
                        const int posui2 = i % 4;
                        ((char*)&smem.elements[0][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = paddingLetter;
                    }
                }
            }
            #endif
    
            //wait until masking is done in smem
            warp.sync();
    
            //write smem back to gmem (inverse operation of loading )

            #if 1
            for(int ui = 0; ui < maxUintsPerThread; ui++){
                const int k = ui * warp.num_threads() + warp.thread_rank();
                const int localLoadSegmentId = k / maxUintsPerThread;
                const int loadLane = k % maxUintsPerThread;
    
                // if(localLoadSegmentId < 32){
                    const int segmentLength = smem.segmentLengths[stage][warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[stage][warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);
                    if(loadLane < segmentUints){
                        segmentPtr[loadLane] = smem.elements[0][warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }
            #else

            const int numSegmentsInBatch = cuda::std::min(32, numSegments - batch32 * 32);
            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                // if(localLoadSegmentId < numSegmentsInBatch){
                    const int segmentLength = smem.segmentLengths[stage][warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[stage][warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);

                    for(int loadLane = warp.thread_rank(); loadLane < segmentUints; loadLane += warp.num_threads()){
                        segmentPtr[loadLane] = smem.elements[0][warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }

            #endif
    
            //wait until smem can be re-used
            warp.sync();
    
        }
    }
    
    
    template<int blocksize, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    void callSingleTileStagedKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));

        int maxBlocksPerSM = 0;
    
        constexpr int maxUintsPerThread = maxSegmentLength / 4;
        constexpr int numWarpsInBlock = blocksize / 32;
        constexpr int numStages = 2;
        auto kernel = singleTileStagedKernel<blocksize, maxSegmentLength, numStages, InputOffsetIterator, InputTransformOp>;
        using SmemType = SingleTileStagedKernel1Smem<numWarpsInBlock, maxUintsPerThread, numStages>;
    
        int smembytes = sizeof(SmemType);
    
        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }
    
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));

        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";

        #ifdef WITH_TIMING
        helpers::GpuTimer timer1("maskingsinglethread::singleTileStagedKernel<" + std::to_string(blocksize) + "," + std::to_string(maxSegmentLength) + ">");
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );
        #ifdef WITH_TIMING
        timer1.print();
        #endif
    }

#if defined(__CUDACC__)
    template<int numWarpsInBlock, int maxUintsPerThread, int numStages>
    struct SingleTileStagedAsyncCopyKernel1Smem{
        uint32_t elements[2][numWarpsInBlock][maxUintsPerThread][32+1];
        int segmentLengths[numStages][numWarpsInBlock][32];
        size_t segmentOffsets[numStages][numWarpsInBlock][32];
    };


    template<int blocksize, int maxSegmentLength, int numStages, class InputOffsetIterator, class InputTransformOp>
    __global__
    void singleTileStagedAsyncCopyKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        __builtin_assume(blockDim.x == blocksize);
        static_assert(maxSegmentLength % 4 == 0);
        
        auto block = cooperative_groups::this_thread_block();
        auto warp = cooperative_groups::tiled_partition<32>(block);
        
        constexpr int numWarpsInBlock = blocksize / 32;
        const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
        const int warpIdInBlock = warp.meta_group_rank();
        const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;
        
        constexpr int maxUintsPerThread = maxSegmentLength / 4;
        // static_assert(maxUintsPerThread <= 32);
        if constexpr(maxUintsPerThread > 32){
            return;
        }

        static_assert(numStages <= 2);

        using Barrier = cuda::barrier<cuda::thread_scope_block>;
        #pragma nv_diag_suppress static_var_with_dynamic_init
        __shared__ Barrier barriers[numWarpsInBlock*numStages];
        auto* warpBarriers = &barriers[2*warp.meta_group_rank()];
        if(warp.thread_rank() < 2){
            init(warpBarriers + warp.thread_rank(), warp.num_threads());
        }
        warp.sync();
    
        using Smem = SingleTileStagedAsyncCopyKernel1Smem<numWarpsInBlock, maxUintsPerThread, numStages>;
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);


        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };
    
    
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("numSegments %d\n", numSegments);
        // }
    
        //each warp processes 32 sequences
        const int numBatches = cuda::ceil_div(numSegments, 32);
        int stage = 0;


        //load data of first batch
        Barrier::arrival_token arrivalTokens[numStages];

        if(warpIdInGrid < numBatches){
            const int fetchbatch32 = warpIdInGrid;
            const int fetchSegmentId = fetchbatch32 * 32 + warp.thread_rank();
            const int fetchstage = stage;
            int fetchSegmentLength = 0;
            size_t fetchSegmentOffset = 0;
            if(fetchSegmentId < numSegments){
                fetchSegmentLength = d_lengths[fetchSegmentId];
                fetchSegmentOffset = d_offsets[fetchSegmentId];
            }
            smem.segmentLengths[fetchstage][warpIdInBlock][warp.thread_rank()] = fetchSegmentLength;
            smem.segmentOffsets[fetchstage][warpIdInBlock][warp.thread_rank()] = fetchSegmentOffset;
            warp.sync();


            for(int ui = 0; ui < maxUintsPerThread; ui++){
                const int k = ui * warp.num_threads() + warp.thread_rank();
                const int localLoadSegmentId = k / maxUintsPerThread;
                const int loadLane = k % maxUintsPerThread;

                // if(localLoadSegmentId < 32){
                    const int segmentLength = smem.segmentLengths[fetchstage][warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[fetchstage][warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                    if(loadLane < segmentUints){
                        //smem.elements[fetchstage][warpIdInBlock][loadLane][localLoadSegmentId] = segmentPtr[loadLane];
                        cuda::memcpy_async(
                            &smem.elements[fetchstage][warpIdInBlock][loadLane][localLoadSegmentId],
                            &segmentPtr[loadLane],
                            cuda::aligned_size_t<4>(sizeof(uint32_t)),
                            warpBarriers[fetchstage]
                        );
                    }
                // }
            }
            if(fetchstage == 0){
                arrivalTokens[0] = warpBarriers[fetchstage].arrive();
            }else{
                arrivalTokens[1] = warpBarriers[fetchstage].arrive();
            }
        }


        for(int batch32 = warpIdInGrid; batch32 < numBatches; batch32 += numWarpsInGrid, stage = (stage + 1) % numStages){


            //prefetch next data
            if(batch32 + numWarpsInGrid < numBatches){
                const int fetchbatch32 = batch32 + numWarpsInGrid;
                const int fetchSegmentId = fetchbatch32 * 32 + warp.thread_rank();
                const int fetchstage = (stage + 1) % 2;
                int fetchSegmentLength = 0;
                size_t fetchSegmentOffset = 0;
                if(fetchSegmentId < numSegments){
                    fetchSegmentLength = d_lengths[fetchSegmentId];
                    fetchSegmentOffset = d_offsets[fetchSegmentId];
                }
                smem.segmentLengths[fetchstage][warpIdInBlock][warp.thread_rank()] = fetchSegmentLength;
                smem.segmentOffsets[fetchstage][warpIdInBlock][warp.thread_rank()] = fetchSegmentOffset;
                warp.sync();

                for(int ui = 0; ui < maxUintsPerThread; ui++){
                    const int k = ui * warp.num_threads() + warp.thread_rank();
                    const int localLoadSegmentId = k / maxUintsPerThread;
                    const int loadLane = k % maxUintsPerThread;
    
                    // if(localLoadSegmentId < 32){
                        const int segmentLength = smem.segmentLengths[fetchstage][warpIdInBlock][localLoadSegmentId];
                        const size_t segmentOffset = smem.segmentOffsets[fetchstage][warpIdInBlock][localLoadSegmentId];
                        const int segmentUints = cuda::ceil_div(segmentLength, 4);
                        const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
                        if(loadLane < segmentUints){
                            //smem.elements[fetchstage][warpIdInBlock][loadLane][localLoadSegmentId] = segmentPtr[loadLane];
                            cuda::memcpy_async(
                                &smem.elements[fetchstage][warpIdInBlock][loadLane][localLoadSegmentId],
                                &segmentPtr[loadLane],
                                cuda::aligned_size_t<4>(sizeof(uint32_t)),
                                warpBarriers[fetchstage]
                            );
                        }
                    // }
                }
                if(fetchstage == 0){
                    arrivalTokens[0] = warpBarriers[fetchstage].arrive();
                }else{
                    arrivalTokens[1] = warpBarriers[fetchstage].arrive();
                }
            }




            
            const int workSegmentId = batch32 * 32 + warp.thread_rank();
            int workSegmentLength = smem.segmentLengths[stage][warpIdInBlock][warp.thread_rank()];
            size_t workSegmentOffset = smem.segmentOffsets[stage][warpIdInBlock][warp.thread_rank()];
            const int workSegmentUints = cuda::ceil_div(workSegmentLength, 4);
            
            //wait for current elements to be ready
            if(stage == 0){
                warpBarriers[stage].wait(std::move(arrivalTokens[0]));
            }else{
                warpBarriers[stage].wait(std::move(arrivalTokens[1]));
            }
            
    
            #if 1
            if(workSegmentId < numSegments){
                char prev = 0;
                int count = 0;
                int startOfRun = 0;
                int numRuns = 0;
                
    
                for(int uintindex = 0; uintindex < workSegmentUints; uintindex++){
                    uint32_t current4 = smem.elements[stage][warpIdInBlock][uintindex][warp.thread_rank()];
                    current4 = applyInputTransform(current4);
                    smem.elements[stage][warpIdInBlock][uintindex][warp.thread_rank()] = current4;
                    alignas(4) char current4chars[4];
                    memcpy(&current4chars, &current4, 4);
    
                    // if(workSegmentId == 444){
                    //     printf("uintindex %d, load %2d %2d %2d %2d from %p\n", 
                    //         uintindex,
                    //         current4chars[0], current4chars[1], current4chars[2], current4chars[3],
                    //         (void*)&smem.elements[stage][warpIdInBlock][uintindex][warp.thread_rank()]
                    //     );
                    // }
                    
                    if(uintindex == 0){
                        prev = current4chars[0];
                        numRuns++;
                    }
                    for(int posIn4 = 0; posIn4 < 4; posIn4++){
                        const int posInSegment = uintindex * 4 + posIn4;
                        if(posInSegment < workSegmentLength){
                            char cur = current4chars[posIn4];
                            // if(workSegmentId == 444){
                            //     printf("posInSegment %d, prev %d, prevcount %d, cur %d\n", posInSegment, prev, count, cur);
                            // }
                            if(cur == prev){
                                count++;
                            }else{
                                if(count >= maskingThreshold){
                                        // if(workSegmentId == 444){
                                        //     printf("masking loop %d - %d\n", startOfRun, posInSegment);
                                        // }
    
                                        #if 1
                                        for(int i = startOfRun; i < posInSegment; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[stage][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
                                        #else
    
                                        const int nextborder = cuda::round_up(startOfRun, 4);
                                        const int lastborder = cuda::round_down(posInSegment, 4);
    
                                        for(int i = startOfRun; i < nextborder; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[stage][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
    
                                        for(int i = nextborder; i < lastborder; i += 4){
                                            const int uindex2 = i / 4;
                                            ((char4*)&smem.elements[stage][warpIdInBlock][uindex2][warp.thread_rank()])[0] 
                                                = make_char4(maskingLetter, maskingLetter, maskingLetter, maskingLetter);
                                        }
    
                                        for(int i = lastborder; i < posInSegment; i++){
                                            const int uindex2 = i / 4;
                                            const int posui2 = i % 4;
                                            ((char*)&smem.elements[stage][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                        }
    
                                        #endif
                                }
                                prev = cur;
                                count = 1;
                                startOfRun = posInSegment;
    
                                numRuns++;
                            }
                        }
                    }
                }
                //handle last run
                if(count >= maskingThreshold){
    
                    // if(workSegmentId == 444){
                    //     printf("mask last run, startOfRun %d, workSegmentLength %d\n", startOfRun, workSegmentLength);
                    // }
                    for(int i = startOfRun; i < workSegmentLength; i++){
                        const int uindex2 = i / 4;
                        const int posui2 = i % 4;
                        ((char*)&smem.elements[stage][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                    }
                }
        
                if(applyExtraPaddingLetter){
                    for(int i = workSegmentLength; i < cuda::round_up(workSegmentLength,4); i++){
                        const int uindex2 = i / 4;
                        const int posui2 = i % 4;
                        ((char*)&smem.elements[stage][warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = paddingLetter;
                    }
                }
            }
            #endif
    
            //wait until masking is done in smem
            warp.sync();
    
            //write smem back to gmem (inverse operation of loading )

            #if 1
            for(int ui = 0; ui < maxUintsPerThread; ui++){
                const int k = ui * warp.num_threads() + warp.thread_rank();
                const int localLoadSegmentId = k / maxUintsPerThread;
                const int loadLane = k % maxUintsPerThread;
    
                // if(localLoadSegmentId < 32){
                    const int segmentLength = smem.segmentLengths[stage][warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[stage][warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);
                    if(loadLane < segmentUints){
                        segmentPtr[loadLane] = smem.elements[stage][warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }
            #else

            const int numSegmentsInBatch = cuda::std::min(32, numSegments - batch32 * 32);
            for(int localLoadSegmentId = 0; localLoadSegmentId < 32; localLoadSegmentId++){
                // if(localLoadSegmentId < numSegmentsInBatch){
                    const int segmentLength = smem.segmentLengths[stage][warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[stage][warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);

                    for(int loadLane = warp.thread_rank(); loadLane < segmentUints; loadLane += warp.num_threads()){
                        segmentPtr[loadLane] = smem.elements[stage][warpIdInBlock][loadLane][localLoadSegmentId];
                    }
                // }
            }

            #endif
    
            //wait until smem can be re-used
            warp.sync();
    
        }
    }
    
    
    template<int blocksize, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    void callSingleTileStagedAsyncCopyKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));

        int maxBlocksPerSM = 0;
    
        constexpr int maxUintsPerThread = maxSegmentLength / 4;
        constexpr int numWarpsInBlock = blocksize / 32;
        constexpr int numStages = 2;
        auto kernel = singleTileStagedAsyncCopyKernel<blocksize, maxSegmentLength, numStages, InputOffsetIterator, InputTransformOp>;
        using SmemType = SingleTileStagedAsyncCopyKernel1Smem<numWarpsInBlock, maxUintsPerThread, numStages>;
    
        int smembytes = sizeof(SmemType);
    
        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }
    
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));

        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";
    
        #ifdef WITH_TIMING
        helpers::GpuTimer timer1("maskingsinglethread::singleTileStagedAsyncCopyKernel<" + std::to_string(blocksize) + "," + std::to_string(maxSegmentLength) + ">");
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );
        #ifdef WITH_TIMING
        timer1.print();
        #endif
    }
#endif // defined(__CUDACC__)

}



namespace maskingsinglegroup{

    template<class InputOffsetIterator>
    __global__
    void applyPaddingLetterKernel(
        char* __restrict__ d_data,
        const int* __restrict__ d_inputLengths,
        InputOffsetIterator d_inputOffsets,
        int numSegments,
        char paddingLetter,
        int padTo
    ){
        const int tid = threadIdx.x + blockIdx.x * blockDim.x;
        const int segmentId = tid / padTo;
        const int tidInSegment = tid % padTo;

        if(segmentId < numSegments){

            const size_t segmentBeginOffset = d_inputOffsets[segmentId];
            const int segmentLength = d_inputLengths[segmentId];
            const int numPaddingBytes = (segmentLength % padTo == 0) ? 0 : padTo - (segmentLength % padTo);
            if(tidInSegment < numPaddingBytes){
                d_data[segmentBeginOffset + segmentLength + tidInSegment] = paddingLetter;
            }
        }
    }


    template<int numGroupsInBlock, int maxUintsPerSequence, int maxMaskingHaloUints>
    struct SingleTileWithPlainFlagsKernel1Smem{
        // uint32_t elements[numGroupsInBlock][maxMaskingHaloUints + maxUintsPerSequence + maxMaskingHaloUints];
        uint32_t byteFlags[numGroupsInBlock][maxMaskingHaloUints + maxUintsPerSequence + maxMaskingHaloUints];
    };
    
    template<int blocksize, int groupsize, int maxMaskingThreshold, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    __global__
    void singleTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        __builtin_assume(blockDim.x == blocksize);
        static_assert(maxSegmentLength % 4 == 0);
        
        // static_assert(groupsize == 32);
        auto block = cooperative_groups::this_thread_block();
        auto warp = cooperative_groups::tiled_partition<32>(block);
        auto group = cooperative_groups::tiled_partition<groupsize>(block);
        
        // constexpr int numWarpsInBlock = blocksize / 32;
        // const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
        // const int warpIdInBlock = warp.meta_group_rank();
        // const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;

        constexpr int numGroupsInBlock = blocksize / groupsize;
        const int numGroupsInGrid = numGroupsInBlock * gridDim.x;
        const int groupIdInBlock = group.meta_group_rank();
        const int groupIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        
        constexpr int maxUintsPerSegment = maxSegmentLength / 4;
        constexpr int maxUintsPerThread = cuda::ceil_div(maxUintsPerSegment, groupsize);
        constexpr int maxMaskingHaloUints = cuda::ceil_div(maxMaskingThreshold, 4);


        const int maskingHaloUints = cuda::ceil_div(maskingThreshold, 4);
        __builtin_assume(maskingThreshold <= maxMaskingThreshold);
        __builtin_assume(maskingHaloUints <= maxMaskingHaloUints);
    
        using Smem = SingleTileWithPlainFlagsKernel1Smem<numGroupsInBlock, maxUintsPerSegment, maxMaskingHaloUints>;
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);

        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };
    

    
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("numSegments %d\n", numSegments);
        // }
    
        //each group processes 1 sequence

        for(int segmentId = groupIdInGrid; segmentId < numSegments; segmentId += numGroupsInGrid){
            const int segmentLength = d_lengths[segmentId];
            const size_t segmentOffset = d_offsets[segmentId];
            const int segmentUints = cuda::ceil_div(segmentLength, 4);
            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
            uint32_t* const segmentOutputPtr = (uint32_t*)(d_outputData + segmentOffset);
            __builtin_assume(segmentLength <= maxSegmentLength);
            __builtin_assume(segmentUints <= maxUintsPerSegment);

            uint32_t myTmpElements[maxUintsPerThread];
            for(int i = 0; i < maxUintsPerThread; i++){
                const int index = i * group.num_threads() + group.thread_rank();
                if(index < segmentUints){
                    myTmpElements[i] = segmentPtr[index];
                }
            }
            for(int i = 0; i < maxUintsPerThread; i++){
                const int index = i * group.num_threads() + group.thread_rank();
                myTmpElements[i] = applyInputTransform(myTmpElements[i]);

                const int threadIdOfLeftNeighborElements = group.thread_rank() > 0 ? group.thread_rank() - 1 : group.num_threads() - 1;
                uint32_t shflelem = myTmpElements[i];
                if(i > 0 && group.thread_rank() == group.num_threads() - 1){
                    shflelem = myTmpElements[i-1];
                }
                uint32_t neighborElements = group.shfl(shflelem, threadIdOfLeftNeighborElements);
                if(i == 0 && group.thread_rank() == 0){
                    neighborElements = 0;
                }

                alignas(4) char charElements[8];
                memcpy(&charElements[0], &neighborElements, 4);
                memcpy(&charElements[4], &myTmpElements[i], 4);

                alignas(4) char currentflags[4];
                for(int k = 4; k < 8; k++){
                    currentflags[k-4] = (charElements[k-1] != charElements[k]);
                }
                uint32_t flagsvec;
                memcpy(&flagsvec, &currentflags[0], 4);

                if(index < segmentUints){
                    smem.byteFlags[groupIdInBlock][maxMaskingHaloUints+index] = flagsvec;
                }
            }

            //if(isFirstChunk)
            if(group.thread_rank() == 0){
                //mark begin of whole segment
                ((char*)&smem.byteFlags[groupIdInBlock][maxMaskingHaloUints])[0] = 1; 
            }
            //if(isLastChunk)
            if(group.thread_rank() == segmentUints % group.num_threads()){
                //mark first oob element
                ((char*)&smem.byteFlags[groupIdInBlock][maxMaskingHaloUints])[segmentLength] = 1; 
                
            }
            group.sync();

            //counting
            for(int x = 0; x < maxUintsPerThread; x++){
                const int index = x * group.num_threads() + group.thread_rank();
                
                alignas(4) char charflagarray[4*(maxMaskingHaloUints + 1 + maxMaskingHaloUints)];
                for(int i = 0; i < (maxMaskingHaloUints + 1 + maxMaskingHaloUints); i++){
                    uint32_t tmp = 0x01010101;
                    memcpy(&charflagarray[4*i], &tmp, 4);
                }

                if(index < segmentUints){
                    uint32_t center = smem.byteFlags[groupIdInBlock][maxMaskingHaloUints + index];
                    memcpy(&charflagarray[4*maxMaskingHaloUints], &center, 4);

                    #pragma unroll
                    for(int i = 0; i < maxMaskingHaloUints; i++){
                        if(i < maskingHaloUints){
                            uint32_t l = smem.byteFlags[groupIdInBlock][maxMaskingHaloUints + index-1-i];
                            uint32_t r = smem.byteFlags[groupIdInBlock][maxMaskingHaloUints + index+1+i];
                            memcpy(&charflagarray[4*(maxMaskingHaloUints-1-i)], &l, 4);
                            memcpy(&charflagarray[4*(maxMaskingHaloUints+1+i)], &r, 4);

                            // left4[maxMaskingHaloUints-1-i] = smem.byteFlags[groupIdInBlock][maxMaskingHaloUints + index-1-i];
                            // right4[i] = smem.byteFlags[groupIdInBlock][maxMaskingHaloUints + index+1+i];
                        }
                    }
                }
                
                int counts[4];

                //for first element, count left
                counts[0] = 0;
                for(int i = 0; i < maxMaskingThreshold; i++){
                    counts[0]++;
                    if(charflagarray[4*maxMaskingHaloUints - i] == 1){
                        break;
                    }
                }
                //for first element, count right
                for(int i = 0; i < maxMaskingThreshold-1; i++){
                    if(charflagarray[4*maxMaskingHaloUints + 1 + i] == 1){
                        break;
                    }
                    counts[0]++;
                }
                for(int k = 1; k < 4; k++){
                    if(charflagarray[4*maxMaskingHaloUints+k] == 0){
                        counts[k] = counts[k-1]; //same segment, same masking criterion
                    }else{
                        //new segment begins at k, count to the right
                        counts[k] = 1;
                        for(int i = 0; i < maxMaskingThreshold-1; i++){
                            if(charflagarray[4*maxMaskingHaloUints + 1 + k + i] == 1){
                                break;
                            }
                            counts[k]++;
                        }
                    }
                }
                
                alignas(4) char charElements[4];
                memcpy(&charElements[0], &myTmpElements[x], 4);
                for(int k = 0; k < 4; k++){
                    if(counts[k] >= maskingThreshold){
                        charElements[k] = maskingLetter;
                    }
                }
                if(applyExtraPaddingLetter){
                    if(index == segmentUints-1){
                        for(int k = 0; k < 4; k++){
                            if(index * 4 + k >= segmentLength){
                                charElements[k] = paddingLetter;
                            }
                        }
                    }
                }
                memcpy(&myTmpElements[x], &charElements[0], 4);

                if(index < segmentUints){
                    segmentOutputPtr[index] = myTmpElements[x];
                }
            }
            group.sync();

        }
    }
    
    
    template<int blocksize, int groupsize, int maxMaskingThreshold, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    void callSingleTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        if(maskingThreshold > maxMaskingThreshold){
            throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
        }
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    
        int maxBlocksPerSM = 0;
    
        constexpr int maxUintsPerSegment = maxSegmentLength / 4;
        constexpr int maxMaskingHaloUints = cuda::ceil_div(maxMaskingThreshold, 4);
        constexpr int numGroupsInBlock = blocksize / groupsize;
        auto kernel = singleTileWithPlainFlagsKernel<blocksize, groupsize, maxMaskingThreshold, 
            maxSegmentLength, InputOffsetIterator, InputTransformOp>;
        using SmemType = SingleTileWithPlainFlagsKernel1Smem<numGroupsInBlock, maxUintsPerSegment, maxMaskingHaloUints>;
    
        int smembytes = sizeof(SmemType);
    
        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }
    
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));
        
        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";

        #ifdef WITH_TIMING
        std::string timername = "maskingsinglegroup::singleTileWithPlainFlagsKernel<";
        timername += std::to_string(blocksize) + ",";
        timername += std::to_string(groupsize) + ",";
        timername += std::to_string(maxMaskingThreshold) + ",";
        timername += std::to_string(maxSegmentLength) + ">";
        helpers::GpuTimer timer1(timername);
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );
        #ifdef WITH_TIMING
        timer1.print();
        #endif

    }








    template<int numGroupsInBlock, int maxUintsPerTile, int maxMaskingHaloUints>
    struct MultiTileWithPlainFlagsKernel1Smem{
        // uint32_t elements[numGroupsInBlock][maxMaskingHaloUints + maxUintsPerTile + maxMaskingHaloUints];
        uint32_t byteFlags[numGroupsInBlock][maxMaskingHaloUints + maxUintsPerTile + maxMaskingHaloUints];
    };
    
    template<int blocksize, int groupsize, int maxMaskingThreshold, int maxTilesize, class InputOffsetIterator, class InputTransformOp>
    __global__
    void multiTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        __builtin_assume(blockDim.x == blocksize);
        static_assert(maxTilesize % 4 == 0);
        
        // static_assert(groupsize == 32);
        auto block = cooperative_groups::this_thread_block();
        auto warp = cooperative_groups::tiled_partition<32>(block);
        auto group = cooperative_groups::tiled_partition<groupsize>(block);
        
        // constexpr int numWarpsInBlock = blocksize / 32;
        // const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
        // const int warpIdInBlock = warp.meta_group_rank();
        // const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;

        constexpr int numGroupsInBlock = blocksize / groupsize;
        const int numGroupsInGrid = numGroupsInBlock * gridDim.x;
        const int groupIdInBlock = group.meta_group_rank();
        const int groupIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        
        constexpr int maxUintsPerTile = maxTilesize / 4;
        constexpr int maxMaskingHaloUints = cuda::ceil_div(maxMaskingThreshold, 4);
        constexpr int maxUintsPerThread = cuda::ceil_div(maxUintsPerTile+maxMaskingHaloUints, groupsize);

        const int maskingHaloUints = cuda::ceil_div(maskingThreshold, 4);
        __builtin_assume(maskingThreshold <= maxMaskingThreshold);
        __builtin_assume(maskingHaloUints <= maxMaskingHaloUints);
    
        using Smem = MultiTileWithPlainFlagsKernel1Smem<numGroupsInBlock, maxUintsPerTile, maxMaskingHaloUints>;
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);

        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };
    
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("numSegments %d\n", numSegments);
        // }
    
        //each group processes 1 sequence

        for(int segmentId = groupIdInGrid; segmentId < numSegments; segmentId += numGroupsInGrid){
            const int segmentLength = d_lengths[segmentId];
            const size_t segmentOffset = d_offsets[segmentId];
            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
            uint32_t* const segmentOutputPtr = (uint32_t*)(d_outputData + segmentOffset);
            const int segmentUints = cuda::ceil_div(segmentLength, 4);
            
            const int numTiles = cuda::ceil_div(segmentUints, maxUintsPerTile);

            for(int tileNr = 0; tileNr < numTiles; tileNr++){
                const bool isFirstTile = tileNr == 0;
                const bool isLastTile = tileNr == numTiles-1;

                const int tileBegin = tileNr * maxUintsPerTile;
                const int tileEnd = cuda::std::min(segmentUints, tileBegin + maxUintsPerTile);
                const int tileUints = tileEnd - tileBegin;
                // const int tileSize = tileEnd - tileBegin;
                
                //const int tileLoadBegin = isFirstTile ? tileBegin : tileBegin - maskingHaloUints;
                const int tileLoadBegin = tileBegin - maskingHaloUints;
                const int tileLoadEnd = (!isLastTile) ? tileEnd + maskingHaloUints : tileEnd;
                const int tileWithHaloUints = tileLoadEnd - tileLoadBegin;


                uint32_t myTmpElements[maxUintsPerThread];
                for(int i = 0; i < maxUintsPerThread; i++){
                    const int index = i * group.num_threads() + group.thread_rank();
                    if(index < tileWithHaloUints){
                        if(tileBegin - maskingHaloUints + index >= 0){
                            myTmpElements[i] = segmentPtr[tileBegin - maskingHaloUints + index];
                            // printf("thread %d, i %d, load from gmem %d\n", group.thread_rank(), i, tileBegin - maskingHaloUints + index);
                        }
                    }
                }

                // if(group.thread_rank() == 0){
                //     for(int i = 0; i < segmentLength + maxMaskingThreshold; i++){
                //         ((char*)&smem.byteFlags[groupIdInBlock][0])[i] = 42;
                //     }
                // }
                // group.sync();
                // if(group.thread_rank() == 0){
                //     printf("flags\n");
                //     for(int i = 0; i < segmentLength + maxMaskingThreshold; i++){
                //         printf("%d ", ((char*)&smem.byteFlags[groupIdInBlock][0])[i]);
                //     }
                //     printf("\n");
                // }
                // group.sync();

                for(int i = 0; i < maxUintsPerThread; i++){
                    const int index = i * group.num_threads() + group.thread_rank();
                    myTmpElements[i] = applyInputTransform(myTmpElements[i]);

                    const int threadIdOfLeftNeighborElements = group.thread_rank() > 0 ? group.thread_rank() - 1 : group.num_threads() - 1;
                    uint32_t shflelem = myTmpElements[i];
                    if(i > 0 && group.thread_rank() == group.num_threads() - 1){
                        shflelem = myTmpElements[i-1];
                    }
                    uint32_t neighborElements = group.shfl(shflelem, threadIdOfLeftNeighborElements);
                    if(i == 0 && group.thread_rank() == 0){
                        neighborElements = 0;
                    }

                    alignas(4) char charElements[8];
                    memcpy(&charElements[0], &neighborElements, 4);
                    memcpy(&charElements[4], &myTmpElements[i], 4);

                    alignas(4) char currentflags[4];
                    for(int k = 4; k < 8; k++){
                        currentflags[k-4] = (charElements[k-1] != charElements[k]);
                    }
                    uint32_t flagsvec;
                    memcpy(&flagsvec, &currentflags[0], 4);

                    if(index < tileWithHaloUints){
                        smem.byteFlags[groupIdInBlock][index] = flagsvec;
                    }
                }

                if(isFirstTile){
                    if(group.thread_rank() == 0){
                        //mark begin of whole segment
                        ((char*)&smem.byteFlags[groupIdInBlock][maskingHaloUints])[0] = 1; 
                    }
                }
                if(isLastTile){
                    const int segmentLengthInLastTile = segmentLength - tileNr * maxTilesize;
                    if(group.thread_rank() == ((maskingHaloUints + tileUints) % group.num_threads())){
                        //mark first oob element
                        ((char*)&smem.byteFlags[groupIdInBlock][maskingHaloUints])[segmentLengthInLastTile] = 1;                         
                    }
                }
                group.sync();


                //counting
                for(int x = 0; x < maxUintsPerThread; x++){
                    const int index = x * group.num_threads() + group.thread_rank();
                    
                    alignas(4) char charflagarray[4*(maxMaskingHaloUints + 1 + maxMaskingHaloUints)];
                    for(int i = 0; i < (maxMaskingHaloUints + 1 + maxMaskingHaloUints); i++){
                        uint32_t tmp = 0x01010101;
                        memcpy(&charflagarray[4*i], &tmp, 4);
                    }

                    if(maskingHaloUints <= index && index < maskingHaloUints + tileUints){
                        uint32_t center = smem.byteFlags[groupIdInBlock][index];
                        memcpy(&charflagarray[4*maxMaskingHaloUints], &center, 4);

                        #pragma unroll
                        for(int i = 0; i < maxMaskingHaloUints; i++){
                            if(i < maskingHaloUints){
                                uint32_t l = smem.byteFlags[groupIdInBlock][index-1-i];
                                uint32_t r = smem.byteFlags[groupIdInBlock][index+1+i];
                                memcpy(&charflagarray[4*(maxMaskingHaloUints-1-i)], &l, 4);
                                memcpy(&charflagarray[4*(maxMaskingHaloUints+1+i)], &r, 4);

                                // left4[maxMaskingHaloUints-1-i] = smem.byteFlags[groupIdInBlock][maxMaskingHaloUints + index-1-i];
                                // right4[i] = smem.byteFlags[groupIdInBlock][maxMaskingHaloUints + index+1+i];
                            }
                        }
                    }
                    
                    int counts[4];

                    //for first element, count left
                    counts[0] = 0;
                    for(int i = 0; i < maxMaskingThreshold; i++){
                        counts[0]++;
                        if(charflagarray[4*maxMaskingHaloUints - i] == 1){
                            break;
                        }
                    }
                    //for first element, count right
                    for(int i = 0; i < maxMaskingThreshold-1; i++){
                        if(charflagarray[4*maxMaskingHaloUints + 1 + i] == 1){
                            break;
                        }
                        counts[0]++;
                    }
                    for(int k = 1; k < 4; k++){
                        if(charflagarray[4*maxMaskingHaloUints+k] == 0){
                            counts[k] = counts[k-1]; //same segment, same masking criterion
                        }else{
                            //new segment begins at k, count to the right
                            counts[k] = 1;
                            for(int i = 0; i < maxMaskingThreshold-1; i++){
                                if(charflagarray[4*maxMaskingHaloUints + 1 + k + i] == 1){
                                    break;
                                }
                                counts[k]++;
                            }
                        }
                    }
                    
                    alignas(4) char charElements[4];
                    memcpy(&charElements[0], &myTmpElements[x], 4);
                    for(int k = 0; k < 4; k++){
                        if(counts[k] >= maskingThreshold){
                            charElements[k] = maskingLetter;
                        }
                    }
                    if(applyExtraPaddingLetter){
                        if((index-maskingHaloUints) == tileUints-1){
                            for(int k = 0; k < 4; k++){
                                if((tileBegin + (index-maskingHaloUints)) * 4 + k >= segmentLength){
                                    charElements[k] = paddingLetter;
                                }
                            }
                        }
                    }
                    memcpy(&myTmpElements[x], &charElements[0], 4);

                    if(maskingHaloUints <= index && index < maskingHaloUints + tileUints){
                        segmentOutputPtr[tileBegin - maskingHaloUints + index] = myTmpElements[x];
                        // printf("thread %d, i %d, write to gmem %d, tileUints %d, maxUintsPerThread %d, tileWithHaloUints %d\n", 
                        //     group.thread_rank(), x, tileBegin - maskingHaloUints + index, tileUints, maxUintsPerThread, tileWithHaloUints);
                    }
                }
                group.sync();
            }
        }

    }
    
    
    template<int blocksize, int groupsize, int maxMaskingThreshold, int maxTilesize, class InputOffsetIterator, class InputTransformOp>
    void callMultiTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        if(maskingThreshold > maxMaskingThreshold){
            throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
        }
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    
        int maxBlocksPerSM = 0;
    
        constexpr int maxUintsPerTile = maxTilesize / 4;
        constexpr int maxMaskingHaloUints = cuda::ceil_div(maxMaskingThreshold, 4);
        constexpr int numGroupsInBlock = blocksize / groupsize;
        auto kernel = multiTileWithPlainFlagsKernel<blocksize, groupsize, maxMaskingThreshold, 
                maxTilesize, InputOffsetIterator, InputTransformOp>;
        using SmemType = MultiTileWithPlainFlagsKernel1Smem<numGroupsInBlock, maxUintsPerTile, maxMaskingHaloUints>;
    
        int smembytes = sizeof(SmemType);
    
        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }
    
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));
        
        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";

        #ifdef WITH_TIMING
        std::string timername = "maskingsinglegroup::multiTileWithPlainFlagsKernel<";
        timername += std::to_string(blocksize) + ",";
        timername += std::to_string(groupsize) + ",";
        timername += std::to_string(maxMaskingThreshold) + ",";
        timername += std::to_string(maxTilesize) + ">";
        helpers::GpuTimer timer1(timername);
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );
        #ifdef WITH_TIMING
        timer1.print();
        #endif

    }

}




namespace maskingsingleblock{

    template<class InputOffsetIterator>
    __global__
    void applyPaddingLetterKernel(
        char* __restrict__ d_data,
        const int* __restrict__ d_inputLengths,
        InputOffsetIterator d_inputOffsets,
        int numSegments,
        char paddingLetter,
        int padTo
    ){
        const int tid = threadIdx.x + blockIdx.x * blockDim.x;
        const int segmentId = tid / padTo;
        const int tidInSegment = tid % padTo;

        if(segmentId < numSegments){

            const size_t segmentBeginOffset = d_inputOffsets[segmentId];
            const int segmentLength = d_inputLengths[segmentId];
            const int numPaddingBytes = (segmentLength % padTo == 0) ? 0 : padTo - (segmentLength % padTo);
            if(tidInSegment < numPaddingBytes){
                d_data[segmentBeginOffset + segmentLength + tidInSegment] = paddingLetter;
            }
        }
    }


    template<int blocksize, int maxUintsPerSequence, int maxMaskingHaloUints>
    struct SingleTileWithPlainFlagsKernel1Smem{
        uint32_t communication_elements[blocksize];
        uint32_t byteFlags[maxMaskingHaloUints + maxUintsPerSequence + maxMaskingHaloUints];
    };
    
    template<int blocksize, int maxMaskingThreshold, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    __global__
    void singleTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.y == 1);
        __builtin_assume(blockDim.z == 1);
        static_assert(maxSegmentLength % 4 == 0);
        
        auto block = cooperative_groups::this_thread_block();
        
        constexpr int maxUintsPerSegment = maxSegmentLength / 4;
        constexpr int maxUintsPerThread = cuda::ceil_div(maxUintsPerSegment, blocksize);
        constexpr int maxMaskingHaloUints = cuda::ceil_div(maxMaskingThreshold, 4);


        const int maskingHaloUints = cuda::ceil_div(maskingThreshold, 4);
        __builtin_assume(maskingThreshold <= maxMaskingThreshold);
        __builtin_assume(maskingHaloUints <= maxMaskingHaloUints);
    
        using Smem = SingleTileWithPlainFlagsKernel1Smem<blocksize, maxUintsPerSegment, maxMaskingHaloUints>;
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);

        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };

    
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("numSegments %d\n", numSegments);
        // }
    
        //each group processes 1 sequence

        for(int segmentId = blockIdx.x; segmentId < numSegments; segmentId += gridDim.x){
            const int segmentLength = d_lengths[segmentId];
            const size_t segmentOffset = d_offsets[segmentId];
            const int segmentUints = cuda::ceil_div(segmentLength, 4);
            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
            uint32_t* const segmentOutputPtr = (uint32_t*)(d_outputData + segmentOffset);
            __builtin_assume(segmentLength <= maxSegmentLength);
            __builtin_assume(segmentUints <= maxUintsPerSegment);

            uint32_t myTmpElements[maxUintsPerThread];
            for(int i = 0; i < maxUintsPerThread; i++){
                const int index = i * block.num_threads() + block.thread_rank();
                if(index < segmentUints){
                    myTmpElements[i] = segmentPtr[index];
                }
            }
            for(int i = 0; i < maxUintsPerThread; i++){
                const int index = i * block.num_threads() + block.thread_rank();
                myTmpElements[i] = applyInputTransform(myTmpElements[i]);

                const int threadIdOfLeftNeighborElements = block.thread_rank() > 0 ? block.thread_rank() - 1 : block.num_threads() - 1;
                uint32_t shflelem = myTmpElements[i];
                if(i > 0 && block.thread_rank() == block.num_threads() - 1){
                    shflelem = myTmpElements[i-1];
                }
                smem.communication_elements[block.thread_rank()] = shflelem;
                block.sync();
                uint32_t neighborElements = smem.communication_elements[threadIdOfLeftNeighborElements];
                block.sync();
                if(i == 0 && block.thread_rank() == 0){
                    neighborElements = 0;
                }

                alignas(4) char charElements[8];
                memcpy(&charElements[0], &neighborElements, 4);
                memcpy(&charElements[4], &myTmpElements[i], 4);

                // if(segmentId == 41){
                //     if(index >= 191 && index <= 194){
                //         printf("index %d, %d %d %d %d %d %d %d %d\n",
                //             index, 
                //             charElements[0], charElements[1], charElements[2], charElements[3],
                //             charElements[4], charElements[5], charElements[6], charElements[7]
                //         );
                //     }
                // }

                alignas(4) char currentflags[4];
                for(int k = 4; k < 8; k++){
                    currentflags[k-4] = (charElements[k-1] != charElements[k]);
                }
                uint32_t flagsvec;
                memcpy(&flagsvec, &currentflags[0], 4);

                if(index < segmentUints){
                    smem.byteFlags[maxMaskingHaloUints+index] = flagsvec;
                }
            }

            //if(isFirstChunk)
            if(block.thread_rank() == 0){
                //mark begin of whole segment
                ((char*)&smem.byteFlags[maxMaskingHaloUints])[0] = 1; 
            }
            //if(isLastChunk)
            if(block.thread_rank() == segmentUints % block.num_threads()){
                //mark first oob element
                ((char*)&smem.byteFlags[maxMaskingHaloUints])[segmentLength] = 1; 
                
            }
            block.sync();

            //counting
            for(int x = 0; x < maxUintsPerThread; x++){
                const int index = x * block.num_threads() + block.thread_rank();
                
                alignas(4) char charflagarray[4*(maxMaskingHaloUints + 1 + maxMaskingHaloUints)];
                for(int i = 0; i < (maxMaskingHaloUints + 1 + maxMaskingHaloUints); i++){
                    uint32_t tmp = 0x01010101;
                    memcpy(&charflagarray[4*i], &tmp, 4);
                }

                if(index < segmentUints){
                    uint32_t center = smem.byteFlags[maxMaskingHaloUints + index];
                    memcpy(&charflagarray[4*maxMaskingHaloUints], &center, 4);

                    #pragma unroll
                    for(int i = 0; i < maxMaskingHaloUints; i++){
                        if(i < maskingHaloUints){
                            uint32_t l = smem.byteFlags[maxMaskingHaloUints + index-1-i];
                            uint32_t r = smem.byteFlags[maxMaskingHaloUints + index+1+i];
                            memcpy(&charflagarray[4*(maxMaskingHaloUints-1-i)], &l, 4);
                            memcpy(&charflagarray[4*(maxMaskingHaloUints+1+i)], &r, 4);

                            // left4[maxMaskingHaloUints-1-i] = smem.byteFlags[maxMaskingHaloUints + index-1-i];
                            // right4[i] = smem.byteFlags[maxMaskingHaloUints + index+1+i];
                        }
                    }
                }
                
                int counts[4];

                //for first element, count left
                counts[0] = 0;
                for(int i = 0; i < maxMaskingThreshold; i++){
                    counts[0]++;
                    if(charflagarray[4*maxMaskingHaloUints - i] == 1){
                        break;
                    }
                }
                //for first element, count right
                for(int i = 0; i < maxMaskingThreshold-1; i++){
                    if(charflagarray[4*maxMaskingHaloUints + 1 + i] == 1){
                        break;
                    }
                    counts[0]++;
                }
                for(int k = 1; k < 4; k++){
                    if(charflagarray[4*maxMaskingHaloUints+k] == 0){
                        counts[k] = counts[k-1]; //same segment, same masking criterion
                    }else{
                        //new segment begins at k, count to the right
                        counts[k] = 1;
                        for(int i = 0; i < maxMaskingThreshold-1; i++){
                            if(charflagarray[4*maxMaskingHaloUints + 1 + k + i] == 1){
                                break;
                            }
                            counts[k]++;
                        }
                    }
                }
                
                alignas(4) char charElements[4];
                memcpy(&charElements[0], &myTmpElements[x], 4);
                for(int k = 0; k < 4; k++){
                    if(counts[k] >= maskingThreshold){
                        charElements[k] = maskingLetter;
                    }
                }
                if(applyExtraPaddingLetter){
                    if((index) == segmentUints-1){
                        for(int k = 0; k < 4; k++){
                            if(((index)) * 4 + k >= segmentLength){
                                charElements[k] = paddingLetter;
                            }
                        }
                    }
                }
                memcpy(&myTmpElements[x], &charElements[0], 4);

                if(index < segmentUints){
                    segmentOutputPtr[index] = myTmpElements[x];
                }
            }
            block.sync();

        }
    }
    
    
    template<int blocksize, int maxMaskingThreshold, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
    void callSingleTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        if(maskingThreshold > maxMaskingThreshold){
            throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
        }
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    
        int maxBlocksPerSM = 0;
    
        constexpr int maxUintsPerSegment = maxSegmentLength / 4;
        constexpr int maxMaskingHaloUints = cuda::ceil_div(maxMaskingThreshold, 4);
        auto kernel = singleTileWithPlainFlagsKernel<blocksize, maxMaskingThreshold, 
            maxSegmentLength, InputOffsetIterator, InputTransformOp>;
        using SmemType = SingleTileWithPlainFlagsKernel1Smem<blocksize, maxUintsPerSegment, maxMaskingHaloUints>;
    
        int smembytes = sizeof(SmemType);
    
        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }
    
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));
        
        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";

        #ifdef WITH_TIMING
        std::string timername = "maskingsingleblock::multiTileWithPlainFlagsKernel<";
        timername += std::to_string(blocksize) + ",";
        timername += std::to_string(maxMaskingThreshold) + ",";
        timername += std::to_string(maxSegmentLength) + ">";
        helpers::GpuTimer timer1(timername);
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );
        #ifdef WITH_TIMING
        timer1.print();
        #endif

    }






    template<int blocksize, int maxUintsPerTile, int maxMaskingHaloUints>
    struct MultiTileWithPlainFlagsKernel1Smem{
        uint32_t communication_elements[blocksize];
        uint32_t byteFlags[maxMaskingHaloUints + maxUintsPerTile + maxMaskingHaloUints];
    };
    
    template<int blocksize, int maxMaskingThreshold, int maxTilesize, class InputOffsetIterator, class InputTransformOp>
    __global__
    void multiTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp
    ){
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.y == 1);
        __builtin_assume(blockDim.z == 1);
        static_assert(maxTilesize % 4 == 0);
        
        auto block = cooperative_groups::this_thread_block();

        constexpr int maxUintsPerTile = maxTilesize / 4;
        constexpr int maxMaskingHaloUints = cuda::ceil_div(maxMaskingThreshold, 4);
        constexpr int maxUintsPerThread = cuda::ceil_div(maxUintsPerTile+maxMaskingHaloUints, blocksize);

        const int maskingHaloUints = cuda::ceil_div(maskingThreshold, 4);
        __builtin_assume(maskingThreshold <= maxMaskingThreshold);
        __builtin_assume(maskingHaloUints <= maxMaskingHaloUints);
    
        using Smem = MultiTileWithPlainFlagsKernel1Smem<blocksize, maxUintsPerTile, maxMaskingHaloUints>;
        extern __shared__ char externsmem[];
        Smem& smem = *((Smem*)externsmem);

        auto applyInputTransform = [&](uint32_t data){
            alignas(4) char elements[4];
            memcpy(&elements[0], &data, 4);
            for(int i = 0; i < 4; i++){
                elements[i] = inputTransformOp(elements[i]);
            }
            memcpy(&data, &elements[0], 4);
            return data;
        };
    
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("numSegments %d\n", numSegments);
        // }
    
        //each group processes 1 sequence

        for(int segmentId = blockIdx.x; segmentId < numSegments; segmentId += gridDim.x){
            const int segmentLength = d_lengths[segmentId];
            const size_t segmentOffset = d_offsets[segmentId];
            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);
            uint32_t* const segmentOutputPtr = (uint32_t*)(d_outputData + segmentOffset);
            const int segmentUints = cuda::ceil_div(segmentLength, 4);
            
            const int numTiles = cuda::ceil_div(segmentUints, maxUintsPerTile);

            for(int tileNr = 0; tileNr < numTiles; tileNr++){
                const bool isFirstTile = tileNr == 0;
                const bool isLastTile = tileNr == numTiles-1;

                const int tileBegin = tileNr * maxUintsPerTile;
                const int tileEnd = cuda::std::min(segmentUints, tileBegin + maxUintsPerTile);
                const int tileUints = tileEnd - tileBegin;
                // const int tileSize = tileEnd - tileBegin;
                
                //const int tileLoadBegin = isFirstTile ? tileBegin : tileBegin - maskingHaloUints;
                const int tileLoadBegin = tileBegin - maskingHaloUints;
                const int tileLoadEnd = (!isLastTile) ? tileEnd + maskingHaloUints : tileEnd;
                const int tileWithHaloUints = tileLoadEnd - tileLoadBegin;


                uint32_t myTmpElements[maxUintsPerThread];
                for(int i = 0; i < maxUintsPerThread; i++){
                    const int index = i * block.num_threads() + block.thread_rank();
                    if(index < tileWithHaloUints){
                        if(tileBegin - maskingHaloUints + index >= 0){
                            myTmpElements[i] = segmentPtr[tileBegin - maskingHaloUints + index];
                            // printf("thread %d, i %d, load from gmem %d\n", group.thread_rank(), i, tileBegin - maskingHaloUints + index);
                        }
                    }
                }

                // if(group.thread_rank() == 0){
                //     for(int i = 0; i < segmentLength + maxMaskingThreshold; i++){
                //         ((char*)&smem.byteFlags[groupIdInBlock][0])[i] = 42;
                //     }
                // }
                // group.sync();
                // if(group.thread_rank() == 0){
                //     printf("flags\n");
                //     for(int i = 0; i < segmentLength + maxMaskingThreshold; i++){
                //         printf("%d ", ((char*)&smem.byteFlags[groupIdInBlock][0])[i]);
                //     }
                //     printf("\n");
                // }
                // group.sync();

                for(int i = 0; i < maxUintsPerThread; i++){
                    const int index = i * block.num_threads() + block.thread_rank();
                    myTmpElements[i] = applyInputTransform(myTmpElements[i]);

                    const int threadIdOfLeftNeighborElements = block.thread_rank() > 0 ? block.thread_rank() - 1 : block.num_threads() - 1;
                    uint32_t shflelem = myTmpElements[i];
                    if(i > 0 && block.thread_rank() == block.num_threads() - 1){
                        shflelem = myTmpElements[i-1];
                    }
                    smem.communication_elements[block.thread_rank()] = shflelem;
                    block.sync();
                    uint32_t neighborElements = smem.communication_elements[threadIdOfLeftNeighborElements];
                    block.sync();
                    if(i == 0 && block.thread_rank() == 0){
                        neighborElements = 0;
                    }

                    alignas(4) char charElements[8];
                    memcpy(&charElements[0], &neighborElements, 4);
                    memcpy(&charElements[4], &myTmpElements[i], 4);

                    alignas(4) char currentflags[4];
                    for(int k = 4; k < 8; k++){
                        currentflags[k-4] = (charElements[k-1] != charElements[k]);
                    }
                    uint32_t flagsvec;
                    memcpy(&flagsvec, &currentflags[0], 4);

                    if(index < tileWithHaloUints){
                        smem.byteFlags[index] = flagsvec;
                    }
                }

                if(isFirstTile){
                    if(block.thread_rank() == 0){
                        //mark begin of whole segment
                        ((char*)&smem.byteFlags[maskingHaloUints])[0] = 1; 
                    }
                }
                if(isLastTile){
                    const int segmentLengthInLastTile = segmentLength - tileNr * maxTilesize;
                    if(block.thread_rank() == ((maskingHaloUints + tileUints) % block.num_threads())){
                        //mark first oob element
                        ((char*)&smem.byteFlags[maskingHaloUints])[segmentLengthInLastTile] = 1;                         
                    }
                }
                block.sync();


                //counting
                for(int x = 0; x < maxUintsPerThread; x++){
                    const int index = x * block.num_threads() + block.thread_rank();
                    
                    alignas(4) char charflagarray[4*(maxMaskingHaloUints + 1 + maxMaskingHaloUints)];
                    for(int i = 0; i < (maxMaskingHaloUints + 1 + maxMaskingHaloUints); i++){
                        uint32_t tmp = 0x01010101;
                        memcpy(&charflagarray[4*i], &tmp, 4);
                    }

                    if(maskingHaloUints <= index && index < maskingHaloUints + tileUints){
                        uint32_t center = smem.byteFlags[index];
                        memcpy(&charflagarray[4*maxMaskingHaloUints], &center, 4);

                        #pragma unroll
                        for(int i = 0; i < maxMaskingHaloUints; i++){
                            if(i < maskingHaloUints){
                                uint32_t l = smem.byteFlags[index-1-i];
                                uint32_t r = smem.byteFlags[index+1+i];
                                memcpy(&charflagarray[4*(maxMaskingHaloUints-1-i)], &l, 4);
                                memcpy(&charflagarray[4*(maxMaskingHaloUints+1+i)], &r, 4);
                            }
                        }
                    }
                    
                    int counts[4];

                    //for first element, count left
                    counts[0] = 0;
                    for(int i = 0; i < maxMaskingThreshold; i++){
                        counts[0]++;
                        if(charflagarray[4*maxMaskingHaloUints - i] == 1){
                            break;
                        }
                    }
                    //for first element, count right
                    for(int i = 0; i < maxMaskingThreshold-1; i++){
                        if(charflagarray[4*maxMaskingHaloUints + 1 + i] == 1){
                            break;
                        }
                        counts[0]++;
                    }
                    for(int k = 1; k < 4; k++){
                        if(charflagarray[4*maxMaskingHaloUints+k] == 0){
                            counts[k] = counts[k-1]; //same segment, same masking criterion
                        }else{
                            //new segment begins at k, count to the right
                            counts[k] = 1;
                            for(int i = 0; i < maxMaskingThreshold-1; i++){
                                if(charflagarray[4*maxMaskingHaloUints + 1 + k + i] == 1){
                                    break;
                                }
                                counts[k]++;
                            }
                        }
                    }
                    
                    alignas(4) char charElements[4];
                    memcpy(&charElements[0], &myTmpElements[x], 4);
                    for(int k = 0; k < 4; k++){
                        if(counts[k] >= maskingThreshold){
                            charElements[k] = maskingLetter;
                        }
                    }
                    if(applyExtraPaddingLetter){
                        if((index-maskingHaloUints) == tileUints-1){
                            for(int k = 0; k < 4; k++){
                                if((tileBegin + (index-maskingHaloUints)) * 4 + k >= segmentLength){
                                    charElements[k] = paddingLetter;
                                }
                            }
                        }
                    }
                    memcpy(&myTmpElements[x], &charElements[0], 4);

                    if(maskingHaloUints <= index && index < maskingHaloUints + tileUints){
                        segmentOutputPtr[tileBegin - maskingHaloUints + index] = myTmpElements[x];
                    }
                }
                block.sync();
            }
        }

    }
    
    
    template<int blocksize, int maxMaskingThreshold, int maxTilesize, class InputOffsetIterator, class InputTransformOp>
    void callMultiTileWithPlainFlagsKernel(
        char* __restrict__ d_outputData,
        const int* __restrict__ d_lengths,
        InputOffsetIterator d_offsets,
        const char* __restrict__ d_inputData,
        int numSegments,
        int maskingThreshold,
        char maskingLetter,
        bool applyExtraPaddingLetter,
        char paddingLetter,
        InputTransformOp inputTransformOp,
        cudaStream_t stream
    ){
        if(maskingThreshold > maxMaskingThreshold){
            throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
        }
        int deviceId;
        int numSMs;
        CUDACHECK(cudaGetDevice(&deviceId));
        CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    
        int maxBlocksPerSM = 0;
    
        constexpr int maxUintsPerTile = maxTilesize / 4;
        constexpr int maxMaskingHaloUints = cuda::ceil_div(maxMaskingThreshold, 4);
        auto kernel = multiTileWithPlainFlagsKernel<blocksize, maxMaskingThreshold, 
                maxTilesize, InputOffsetIterator, InputTransformOp>;
        using SmemType = MultiTileWithPlainFlagsKernel1Smem<blocksize, maxUintsPerTile, maxMaskingHaloUints>;
    
        int smembytes = sizeof(SmemType);
    
        if(smembytes > 48*1024){
            cudaError_t status = cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes);
            if(status != cudaSuccess){
                cudaGetLastError();
                throw std::runtime_error("Unable to set cudaFuncAttributeMaxDynamicSharedMemorySize to " + std::to_string(smembytes));
            }
        }
    
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            smembytes
        ));
        
        const int numBlocks = numSMs * maxBlocksPerSM;
        
        // std::cout << "numSegments "  << numSegments << ", num blocks " << numBlocks << "\n";

        #ifdef WITH_TIMING
        std::string timername = "maskingsingleblock::multiTileWithPlainFlagsKernel<";
        timername += std::to_string(blocksize) + ",";
        timername += std::to_string(maxMaskingThreshold) + ",";
        timername += std::to_string(maxTilesize) + ">";
        helpers::GpuTimer timer1(timername);
        #endif
        kernel<<<numBlocks, blocksize, smembytes, stream>>>(
            d_outputData,
            d_lengths,
            d_offsets,
            d_inputData,
            numSegments,
            maskingThreshold,
            maskingLetter,
            applyExtraPaddingLetter,
            paddingLetter,
            inputTransformOp
        );
        #ifdef WITH_TIMING
        timer1.print();
        #endif

    }

}





LIBMARV_NAMESPACE_WITH_NESTING_END

#endif