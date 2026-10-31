#ifndef LIBMARV_SEQUENCE_MASKING_OLD_CUH
#define LIBMARV_SEQUENCE_MASKING_OLD_CUH

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
#include <cuda_fp16.h>
#include <cuda_fp8.h>
#else
#include <hip/hip_fp16.h>
#include <hip/hip_fp8.h>
#endif

#include "cuda_errorcheck.cuh"
#include "util.cuh"
#include "mathops.cuh"

#if defined(__CUDACC__)
#include <cooperative_groups.h>
#include <cooperative_groups/scan.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/memcpy_async.h>
#else
#include <hip/hip_cooperative_groups.h>
#endif

// #include "cuda_hip_compatibility.cuh"

#include "namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN



template<int blocksize>
__global__
void setPaddingBytesKernel(
    char* d_data,
    const int* d_lengths,
    const size_t* d_offsets,
    int numSegments,
    char paddingLetter
){
    auto group = cooperative_groups::tiled_partition<4>(cooperative_groups::this_thread_block());
    constexpr int numGroupsInBlock = blocksize / 4;
    const int numGroupsInGrid = numGroupsInBlock * gridDim.x;
    const int groupIdInBlock = group.meta_group_rank();
    const int groupIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 4;

    for(int s = groupIdInGrid; s < numSegments; s += numGroupsInGrid){
        const size_t offset = d_offsets[s];
        const int length = d_lengths[s];
        const int paddedLength = cuda::ceil_div(length, 4) * 4;

        if(length + group.thread_rank() < paddedLength){
            d_data[offset + length + group.thread_rank()] = paddingLetter;
        }
    }
}



/*
    Requirements:

*/
template<class InputDataIterator, class InputOffsetIterator>
void maskSequences_rle_rld_impl(
    char* d_outputData,
    InputDataIterator d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter,
    cudaStream_t stream
){
    using LetterIdTuple = cuda::std::tuple<char, uint16_t>;
    auto inputDataWithSegmentId = thrust::make_transform_iterator(
        thrust::make_counting_iterator<size_t>(0),
        cuda::proclaim_return_type<LetterIdTuple>([=] __device__ (size_t index){
            char letter = d_inputData[index];

            auto upIter = thrust::upper_bound(thrust::seq, d_inputOffsets, d_inputOffsets + numInputSequences, index);
            int segmentId = numInputSequences-1;
            if(upIter != d_inputOffsets + numInputSequences){
                segmentId = cuda::std::distance(d_inputOffsets, upIter) - 1;
            }
            return LetterIdTuple(letter, segmentId % 65536);
        })
    );

    // thrust::device_vector<LetterIdTuple> foo(inputDataWithSegmentId, inputDataWithSegmentId + totalInputElements);
    // std::cout << "foo\n";
    // for(int i = 0; i< totalInputElements; i++){
    //     LetterIdTuple t = foo[i];
    //     std::cout<< "(" << cuda::std::get<0>(t) << "," << cuda::std::get<1>(t) << ")" << " ";
    // }
    // std::cout << "\n";

    int numRuns;


    // thrust::device_vector<int, thrust_async_allocator<int>> d_numRuns(1, thrust_async_allocator<int>(stream));
    // //thrust::device_vector<char, thrust_async_allocator<char>> d_runValues(totalInputElements, thrust_async_allocator<char>(stream));

    // thrust::device_vector<LetterIdTuple, thrust_async_allocator<LetterIdTuple>> d_runValues_tuple(totalInputElements, thrust_async_allocator<LetterIdTuple>(stream));
    // thrust::device_vector<int, thrust_async_allocator<int>> d_runLengths(totalInputElements, thrust_async_allocator<int>(stream));


    int* d_numRuns_ptr; 
    CUDACHECK(cudaMallocAsync(&d_numRuns_ptr, sizeof(int) * 1, stream));
    // LetterIdTuple* d_runValues_tuple_ptr; 
    // CUDACHECK(cudaMallocAsync(&d_runValues_tuple_ptr, sizeof(LetterIdTuple) * totalInputElements, stream));
    char* d_runValues_ptr; 
    CUDACHECK(cudaMallocAsync(&d_runValues_ptr, sizeof(char) * totalInputElements, stream));
    uint16_t* d_runValues_segmentIds; 
    CUDACHECK(cudaMallocAsync(&d_runValues_segmentIds, sizeof(uint16_t) * totalInputElements, stream));
    int* d_runLengths_ptr; 
    CUDACHECK(cudaMallocAsync(&d_runLengths_ptr, sizeof(int) * totalInputElements, stream));


    auto rle_output = thrust::make_zip_iterator(
        d_runValues_ptr,
        d_runValues_segmentIds
    );


    // auto d_nonTupleRunValues = thrust::make_transform_output_iterator(
    //     d_runValues.begin(),
    //     cuda::proclaim_return_type<char>([] __device__ (LetterIdTuple tup){
    //         return cuda::std::get<0>(tup);
    //     })
    // );

    size_t cubBytes = 0;
    CUDACHECK(cub::DeviceRunLengthEncode::Encode(
        nullptr,
        cubBytes,
        inputDataWithSegmentId,
        rle_output,
        d_runLengths_ptr,
        d_numRuns_ptr,
        totalInputElements,
        stream
    ));

    char* d_temp; cudaMallocAsync(&d_temp, cubBytes, stream);

    CUDACHECK(cub::DeviceRunLengthEncode::Encode(
        d_temp,
        cubBytes,
        inputDataWithSegmentId,
        rle_output,
        d_runLengths_ptr,
        d_numRuns_ptr,
        totalInputElements,
        stream
    ));
    CUDACHECK(cudaFreeAsync(d_temp, stream));

    CUDACHECK(cudaMemcpyAsync(&numRuns, d_numRuns_ptr, sizeof(int), cudaMemcpyDeviceToHost, stream));
    CUDACHECK(cudaStreamSynchronize(stream));
    std::cout << "maskSequences: totalInputElements " << totalInputElements << ", numRuns " << numRuns << "\n";

    // std::cout << "d_runLengths\n";
    // for(int i = 0; i < numRuns; i++){
    //     std::cout << d_runLengths[i] << " ";
    // }
    // std::cout << "\n";

    #if 1

    //run-length decode

    thrust::inclusive_scan(
        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
        d_runLengths_ptr,
        d_runLengths_ptr + numRuns,
        d_runLengths_ptr
    );

    thrust::transform(
        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
        thrust::make_counting_iterator<size_t>(1),
        thrust::make_counting_iterator<size_t>(1) + totalInputElements,
        d_outputData,
        cuda::proclaim_copyable_arguments(
            [
                runLengthInclPS = d_runLengths_ptr, 
                d_runValues_ptr,
                numRuns,
                maskingThreshold,
                maskingLetter
            ] __device__ (size_t index){
                auto it = thrust::lower_bound(thrust::seq, runLengthInclPS, runLengthInclPS + numRuns, index);
                const int rle_segment = cuda::std::distance(runLengthInclPS, it);
                const int rle_length = rle_segment == 0 ? runLengthInclPS[0] : runLengthInclPS[rle_segment] - runLengthInclPS[rle_segment-1];
                
                if(rle_length >= maskingThreshold){
                    return maskingLetter;
                }else{
                    return d_runValues_ptr[rle_segment];
                }
            }
        )
    );


    #else

        //apply masking letter to runs >= maskingThreshold
        thrust::scatter_if(
            thrust::cuda::par_nosync.on(stream),
            thrust::make_constant_iterator(maskingLetter),
            thrust::make_constant_iterator(maskingLetter) + numRuns,
            thrust::make_counting_iterator(0),
            d_runLengths_ptr,
            d_runValues_ptr,
            // thrust::make_transform_output_iterator(
            //     d_runValues_tuple_ptr,
            //     cuda::proclaim_return_type<LetterIdTuple>([] __device__ (char letter){
            //         return LetterIdTuple(letter, 0);
            //     })
            // ),
            [=] __device__ (int runlength){
                return runlength >= maskingThreshold;
            }
        );


        //run-length decode

        thrust::inclusive_scan(
            thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
            d_runLengths_ptr,
            d_runLengths_ptr + numRuns,
            d_runLengths_ptr
        );

        auto srcLocations = thrust::make_transform_iterator(
            thrust::counting_iterator<int>(1),
            cuda::proclaim_return_type<int>([runLengthInclPS = d_runLengths_ptr, numRuns] __device__ (int index){
                auto it = thrust::lower_bound(thrust::seq, runLengthInclPS, runLengthInclPS + numRuns, index);
                return int(cuda::std::distance(runLengthInclPS, it));
            })
        );

        thrust::gather(
            thrust::cuda::par_nosync.on(stream),
            srcLocations,
            srcLocations + totalInputElements,
            d_runValues_ptr,
            // thrust::make_transform_iterator(
            //     d_runValues_tuple_ptr,
            //     cuda::proclaim_return_type<char>([] __device__ (LetterIdTuple tup){
            //         return cuda::std::get<0>(tup);
            //     })
            // ),    
            d_outputData
        );

    #endif

setPaddingBytesKernel<128><<<cuda::ceil_div(numInputSequences, 128 / 4), 128, 0, stream>>>(
        d_outputData,
        d_inputLengths,
        d_inputOffsets,
        numInputSequences,
        paddingLetter
    );
    CUDACHECKASYNC;


    CUDACHECK(cudaFreeAsync(d_runLengths_ptr, stream));
    // CUDACHECK(cudaFreeAsync(d_runValues_tuple_ptr, stream));
    CUDACHECK(cudaFreeAsync(d_runValues_ptr, stream));
    CUDACHECK(cudaFreeAsync(d_runValues_segmentIds, stream));
    CUDACHECK(cudaFreeAsync(d_numRuns_ptr, stream));
}




/*
    Requirements:

*/
template<class InputDataIterator, class InputOffsetIterator>
void maskSequences_simple_for_each_impl(
    char* d_outputData,
    InputDataIterator d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter,
    cudaStream_t stream
){
    using LetterIdTuple = cuda::std::tuple<char, int>;
    auto inputDataWithSegmentId = thrust::make_transform_iterator(
        thrust::make_counting_iterator<size_t>(0),
        cuda::proclaim_return_type<LetterIdTuple>([=] __device__ (size_t index){
            char letter = d_inputData[index];

            auto upIter = thrust::upper_bound(thrust::seq, d_inputOffsets, d_inputOffsets + numInputSequences, index);
            int segmentId = numInputSequences-1;
            if(upIter != d_inputOffsets + numInputSequences){
                segmentId = cuda::std::distance(d_inputOffsets, upIter) - 1;
            }
            return LetterIdTuple(letter, segmentId);
        })
    );

    thrust::for_each(
        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
        thrust::make_counting_iterator<int64_t>(0),
        thrust::make_counting_iterator<int64_t>(0) + totalInputElements,
        [
            maskingThreshold,
            maskingLetter,
            inputDataWithSegmentId,
            totalInputElements,
            d_outputData
        ]__device__ (int64_t index){
            const LetterIdTuple myData = inputDataWithSegmentId[index];

            int countLeft = 0;
            int countRight = 0;

            //check to the left
            for(int64_t i = index - 1; i > index - maskingThreshold; i--){
                if(i < 0){
                    break;
                }
                auto other = inputDataWithSegmentId[i];
                if(other == myData){
                    countLeft++;
                }else{
                    break;
                }
            }

            //check to the right
            for(int64_t i = index + 1; i < index + maskingThreshold; i++){
                if(i >= totalInputElements){
                    break;
                }
                auto other = inputDataWithSegmentId[i];
                if(other == myData){
                    countRight++;
                }else{
                    break;
                }
            }

            int count = countLeft + 1 + countRight;
            if(count >= maskingThreshold){
                d_outputData[index] = maskingLetter;
            }else{
                d_outputData[index] = cuda::std::get<0>(myData);
            }
        }
    );


    setPaddingBytesKernel<128><<<cuda::ceil_div(numInputSequences, 128 / 4), 128, 0, stream>>>(
        d_outputData,
        d_inputLengths,
        d_inputOffsets,
        numInputSequences,
        paddingLetter
    );
    CUDACHECKASYNC;

}



template<int maximumMaskingThreshold, class InputDataIterator, class SegmentIdIterator>
__global__
void maskingKernel(
    char* d_outputData,
    InputDataIterator d_inputData,
    SegmentIdIterator d_segmentIds,
    int64_t totalInputElements,
    int maskingThreshold,
    char maskingLetter
){

    using LetterIdTuple = cuda::std::tuple<char, int>;
    constexpr int64_t maxChunkSize = 1024;
    constexpr int smemElements = maxChunkSize + 2 * (maximumMaskingThreshold-1);
    __shared__ char smem_letters[smemElements];
    __shared__ int smem_segmentIds[smemElements];

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    for(int64_t chunkId = blockIdx.x; chunkId < numChunks; chunkId += gridDim.x){
        const int64_t chunkBegin = chunkId * maxChunkSize;
        const int64_t chunkEnd_excl = cuda::std::min(totalInputElements, (chunkId+1) * maxChunkSize);
        // const int chunkSize = chunkEnd_excl - chunkBegin;

        const int64_t loadBegin = cuda::std::max(int64_t(0), chunkBegin - (maximumMaskingThreshold-1));
        const int64_t loadEnd_excl = cuda::std::min(totalInputElements, chunkEnd_excl + (maximumMaskingThreshold-1));
        const int elementsToLoad = loadEnd_excl - loadBegin;

        // if(threadIdx.x == 0){
            
        //     if(elementsToLoad > smemElements){
        //         printf("loadBegin %ld, chunkBegin %lu, chunkEnd_excl %lu, loadEnd_excl %ld, elementsToLoad %d, smemElements %d\n", 
        //             loadBegin, chunkBegin, chunkEnd_excl, loadEnd_excl, elementsToLoad, smemElements);
        //     }
        //     assert(elementsToLoad <= smemElements);
        // }

        //load chunk to shared memory
        __syncthreads();
        for(int64_t i = loadBegin + threadIdx.x; i < loadEnd_excl; i += blockDim.x){
            smem_letters[i - loadBegin] = d_inputData[i];
            smem_segmentIds[i - loadBegin] = d_segmentIds[i];
        }
        __syncthreads();

        //process chunk
        for(int64_t element = chunkBegin + threadIdx.x; element < chunkEnd_excl; element += blockDim.x){
            const int indexInSmem = element - loadBegin;
            const char myLetter = smem_letters[indexInSmem];
            const int mySegmentId = smem_segmentIds[indexInSmem];

            int count = 1; //account for my letter
            //check to the left
            for(int i = indexInSmem - 1; i > indexInSmem - maskingThreshold; i--){
                if(i < 0){
                    break;
                }
                const char otherLetter = smem_letters[i];
                const int otherSegmentId = smem_segmentIds[i];
                if(myLetter == otherLetter && mySegmentId == otherSegmentId){
                    count++;
                }else{
                    break;
                }
            }

            //check to the right
            for(int64_t i = indexInSmem + 1; i < indexInSmem + maskingThreshold; i++){
                if(i >= elementsToLoad){
                    break;
                }
                const char otherLetter = smem_letters[i];
                const int otherSegmentId = smem_segmentIds[i];
                if(myLetter == otherLetter && mySegmentId == otherSegmentId){
                    count++;
                }else{
                    break;
                }
            }


            if(count >= maskingThreshold){
                d_outputData[element] = maskingLetter;
            }else{
                d_outputData[element] = myLetter;
            }
        }

    }
}

/*
    Requirements:

*/
template<int maximumMaskingThreshold, class InputDataIterator, class InputOffsetIterator>
void maskSequences_custom_kernel_impl(
    char* d_outputData,
    InputDataIterator d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter,
    cudaStream_t stream
){
    assert(maskingThreshold <= maximumMaskingThreshold);

    auto d_segmentIds = thrust::make_transform_iterator(
        thrust::make_counting_iterator<size_t>(0),
        cuda::proclaim_return_type<int>([=] __device__ (size_t index){
            auto upIter = thrust::upper_bound(thrust::seq, d_inputOffsets, d_inputOffsets + numInputSequences, index);
            int segmentId = numInputSequences-1;
            if(upIter != d_inputOffsets + numInputSequences){
                segmentId = cuda::std::distance(d_inputOffsets, upIter) - 1;
            }
            return segmentId;
        })
    );

    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 128;
    auto kernel = maskingKernel<maximumMaskingThreshold, InputDataIterator, decltype(d_segmentIds)>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    kernel<<<numSMs*maxBlocksPerSM, blocksize, 0, stream>>>(
        d_outputData,
        d_inputData,
        d_segmentIds,
        totalInputElements,
        maskingThreshold,
        maskingLetter
    );
    CUDACHECKASYNC;

    setPaddingBytesKernel<128><<<cuda::ceil_div(numInputSequences, 128 / 4), 128, 0, stream>>>(
        d_outputData,
        d_inputLengths,
        d_inputOffsets,
        numInputSequences,
        paddingLetter
    );
    CUDACHECKASYNC;

}




    
template<int maximumMaskingThreshold, int maximumChunkSize, class InputDataIterator, class SegmentIdIterator>
__global__
void maskingKernel2(
    char* d_outputData,
    InputDataIterator d_inputData,
    SegmentIdIterator d_segmentIds,
    int64_t totalInputElements,
    int maskingThreshold,
    char maskingLetter
){

    using LetterIdTuple = cuda::std::tuple<char, int>;
    constexpr int64_t maxChunkSize = maximumChunkSize;
    constexpr int smemElements = maxChunkSize + 2 * (maximumMaskingThreshold-1);

    struct Pair{
        char letter;
        uint16_t id;
    };
    __shared__ Pair smem_pairs[smemElements];

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    for(int64_t chunkId = blockIdx.x; chunkId < numChunks; chunkId += gridDim.x){
        const int64_t chunkBegin = chunkId * maxChunkSize;
        const int64_t chunkEnd_excl = cuda::std::min(totalInputElements, (chunkId+1) * maxChunkSize);
        // const int chunkSize = chunkEnd_excl - chunkBegin;

        const int64_t loadBegin = cuda::std::max(int64_t(0), chunkBegin - (maximumMaskingThreshold-1));
        const int64_t loadEnd_excl = cuda::std::min(totalInputElements, chunkEnd_excl + (maximumMaskingThreshold-1));
        const int elementsToLoad = loadEnd_excl - loadBegin;

        // if(threadIdx.x == 0){
            
        //     if(elementsToLoad > smemElements){
        //         printf("loadBegin %ld, chunkBegin %lu, chunkEnd_excl %lu, loadEnd_excl %ld, elementsToLoad %d, smemElements %d\n", 
        //             loadBegin, chunkBegin, chunkEnd_excl, loadEnd_excl, elementsToLoad, smemElements);
        //     }
        //     assert(elementsToLoad <= smemElements);
        // }

        //load chunk to shared memory
        __syncthreads();
        for(int64_t i = loadBegin + threadIdx.x; i < loadEnd_excl; i += blockDim.x){
            Pair p;
            p.letter = d_inputData[i];
            // p.id = d_segmentIds[i] % (65536);
            p.id = d_segmentIds[i];
            smem_pairs[i - loadBegin] = p;
        }
        __syncthreads();

        //process chunk
        for(int64_t element = chunkBegin + threadIdx.x; element < chunkEnd_excl; element += blockDim.x){
            const int indexInSmem = element - loadBegin;
            const Pair myPair = smem_pairs[indexInSmem];
            const char myLetter = myPair.letter;
            const int mySegmentId = myPair.id;

            int count = 1; //account for my letter
            //check to the left
            for(int i = indexInSmem - 1; i > indexInSmem - maskingThreshold; i--){
                if(i < 0){
                    break;
                }
                const Pair otherPair = smem_pairs[i];
                const char otherLetter = otherPair.letter;
                const int otherSegmentId = otherPair.id;
                if(myLetter == otherLetter && mySegmentId == otherSegmentId){
                    count++;
                }else{
                    break;
                }
            }

            //check to the right
            for(int64_t i = indexInSmem + 1; i < indexInSmem + maskingThreshold; i++){
                if(count >= maskingThreshold){
                    break;
                }
                if(i >= elementsToLoad){
                    break;
                }
                const Pair otherPair = smem_pairs[i];
                const char otherLetter = otherPair.letter;
                const int otherSegmentId = otherPair.id;
                if(myLetter == otherLetter && mySegmentId == otherSegmentId){
                    count++;
                }else{
                    break;
                }
            }

            // if(element >= 736120 && element < 736130){
            //     printf("element %ld, letter %d, count %d\n", element, int(myLetter), count);
            // }


            if(count >= maskingThreshold){
                d_outputData[element] = maskingLetter;
            }else{
                d_outputData[element] = myLetter;
            }
        }

    }
}




template<int maximumMaskingThreshold, class InputDataIterator, class InputOffsetIterator>
void maskSequences_custom_kernel_impl2(
    char* d_outputData,
    InputDataIterator d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter,
    cudaStream_t stream
){
    assert(maskingThreshold <= maximumMaskingThreshold);

    // auto d_segmentIds = thrust::make_transform_iterator(
    //     thrust::make_counting_iterator<size_t>(0),
    //     cuda::proclaim_return_type<int>([=] __device__ (size_t index){
    //         auto upIter = thrust::upper_bound(thrust::seq, d_inputOffsets, d_inputOffsets + numInputSequences, index);
    //         int segmentId = numInputSequences-1;
    //         if(upIter != d_inputOffsets + numInputSequences){
    //             segmentId = cuda::std::distance(d_inputOffsets, upIter) - 1;
    //         }
    //         return segmentId;
    //     })
    // );

    //we only need to distinguish segments in each subchunk in the kernel, ushort is sufficient and provides safe wrap-around for the scan.
    constexpr int maximumChunkSize = 1024;
    uint16_t* d_segmentIds; 
    CUDACHECK(cudaMallocAsync(&d_segmentIds, sizeof(uint16_t) * totalInputElements, stream));
    CUDACHECK(cudaMemsetAsync(d_segmentIds, 0, sizeof(uint16_t) * totalInputElements, stream));
    thrust::scatter(
        thrust::cuda::par_nosync.on(stream),
        thrust::make_constant_iterator<uint16_t>(1),
        thrust::make_constant_iterator<uint16_t>(1) + numInputSequences,
        d_inputOffsets,
        d_segmentIds
    );
    thrust::inclusive_scan(
        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
        d_segmentIds,
        d_segmentIds + totalInputElements,
        d_segmentIds
    );

    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 128;
    auto kernel = maskingKernel2<maximumMaskingThreshold, maximumChunkSize, InputDataIterator, decltype(d_segmentIds)>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    kernel<<<numSMs*maxBlocksPerSM, blocksize, 0, stream>>>(
        d_outputData,
        d_inputData,
        d_segmentIds,
        totalInputElements,
        maskingThreshold,
        maskingLetter
    );
    CUDACHECKASYNC;

    setPaddingBytesKernel<128><<<cuda::ceil_div(numInputSequences, 128 / 4), 128, 0, stream>>>(
        d_outputData,
        d_inputLengths,
        d_inputOffsets,
        numInputSequences,
        paddingLetter
    );
    CUDACHECKASYNC;

    CUDACHECK(cudaFreeAsync(d_segmentIds, stream));

}




template<int blocksize, class InputDataIterator, class SegmentIdIterator>
__global__
void maskingKernel3(
    char* d_outputData,
    InputDataIterator d_inputData,
    SegmentIdIterator d_segmentIds,
    int64_t totalInputElements,
    int maskingThreshold,
    char maskingLetter
){
    __builtin_assume(blocksize == blockDim.x);

    struct RunIdCountPair{
        int runId;
        int count;
    };
    struct RunIdCountPairScanOp{
        __device__
        RunIdCountPair operator()(const RunIdCountPair& left, const RunIdCountPair& right) const{
            RunIdCountPair result;
            result.runId = right.runId;
            if(left.runId == right.runId){
                result.count = right.count + left.count;
            }else{
                result.count = right.count;                
            }
            return result;
        }
    };

    constexpr uint32_t oob_element = 0xFFFFFFFF;



    constexpr int items_per_thread = 16;
    constexpr int items_per_block = blocksize * items_per_thread;

    using BlockExchange = cub::BlockExchange<uint32_t, blocksize, items_per_thread>;
    using BlockDiscontinuity = cub::BlockDiscontinuity<uint32_t, blocksize>;
    using BlockScanInt = cub::BlockScan<int, blocksize>;
    using BlockScanPair = cub::BlockScan<RunIdCountPair, blocksize>;
    

    __shared__ struct{
        typename BlockExchange::TempStorage temp_exchange;
        typename BlockDiscontinuity::TempStorage temp_discontinuity;
        typename BlockScanInt::TempStorage temp_scanInt;
        typename BlockScanPair::TempStorage temp_scanPair;
    } temp;

    __shared__ char smem_letters[items_per_block];
    __shared__ int smem_runLengths[items_per_block];



    /*
        Example for chunking: with 2048 items per block and maskingThreshold == 3, 2045 is the last index that could be 
        the start of a run to mask in the current chunk. index 2046 has to be checked by the next chunk. Thus, next chunk
        needs to start at index 2044 to check if 2046 is the end of a run to mask.
    */

    for(int64_t chunkBegin = blockIdx.x * (items_per_block - 2*(maskingThreshold - 1)) ; 
            chunkBegin < totalInputElements; 
            chunkBegin += gridDim.x * (items_per_block - 2*(maskingThreshold - 1))){

        const int64_t chunkEnd_excl = cuda::std::min(totalInputElements, chunkBegin + items_per_block);
        // const int64_t numInChunk = chunkEnd_excl - chunkBegin;

        uint32_t myElements[items_per_thread];        

        #pragma unroll
        for(int i = 0; i < items_per_thread; i++){
            const int loadIndex = chunkBegin + i * blockDim.x + threadIdx.x;
            if(loadIndex < chunkEnd_excl){
                uint32_t letter = d_inputData[loadIndex];
                //uint32_t id = d_segmentIds[loadIndex] % (65536);
                uint32_t id = d_segmentIds[loadIndex];
                uint32_t element = (letter << 16) | id;
                myElements[i] = element;
            }else{
                myElements[i] = oob_element;
            }
        }
        BlockExchange(temp.temp_exchange).StripedToBlocked(myElements, myElements);
        // __syncthreads();

        // if(threadIdx.x == 0){
        //     printf("elements\n");
        //     #pragma unroll
        //     for(int i = 0; i < items_per_thread; i++){
        //         printf("%u ", myElements[i]);
        //     }
        //     printf("\n");
        // }

        // const int printindex = 134088448;

        // #pragma unroll
        // for(int i = 0; i < items_per_thread; i++){
        //     if(chunkBegin + threadIdx.x * items_per_thread + i >= (printindex-5) && chunkBegin + threadIdx.x * items_per_thread + i <= (printindex+5)){
        //         printf("chunkId %ld, element %ld, %u\n", chunkId, chunkBegin + threadIdx.x * items_per_thread + i, myElements[i]);
        //     }
        // }


        int tailFlags[items_per_thread];
        BlockDiscontinuity(temp.temp_discontinuity).FlagTails(tailFlags, myElements, cuda::std::not_equal_to{});
        #pragma unroll
        for(int i = 0; i < items_per_thread; i++){
            if(chunkBegin + threadIdx.x * items_per_thread + i >= chunkEnd_excl){
                tailFlags[i] = 0;
            }
        }

        // if(threadIdx.x == 0){
        //     printf("tailFlags\n");
        //     #pragma unroll
        //     for(int i = 0; i < items_per_thread; i++){
        //         printf("%u ", tailFlags[i]);
        //     }
        //     printf("\n");
        // }

        // #pragma unroll
        // for(int i = 0; i < items_per_thread; i++){
        //     if(chunkBegin + threadIdx.x * items_per_thread + i >= (printindex-5) && chunkBegin + threadIdx.x * items_per_thread + i <= (printindex+5)){
        //         printf("chunkId %ld, tailFlags %ld, %u\n", chunkId, chunkBegin + threadIdx.x * items_per_thread + i, tailFlags[i]);
        //     }
        // }

        //exclusive sum to compute run ids marked by tail flags
        int prefixsum[items_per_thread];
        BlockScanInt(temp.temp_scanInt).ExclusiveSum(tailFlags, prefixsum);
        // __syncthreads();

        // if(threadIdx.x == 0){
        //     printf("prefixsum\n");
        //     #pragma unroll
        //     for(int i = 0; i < items_per_thread; i++){
        //         printf("%u ", prefixsum[i]);
        //     }
        //     printf("\n");
        // }

        // #pragma unroll
        // for(int i = 0; i < items_per_thread; i++){
        //     if(chunkBegin + threadIdx.x * items_per_thread + i >= (printindex-5) && chunkBegin + threadIdx.x * items_per_thread + i <= (printindex+5)){
        //         printf("chunkId %ld, prefixsum %ld, %u\n", chunkId, chunkBegin + threadIdx.x * items_per_thread + i, prefixsum[i]);
        //     }
        // }



        //inclusive scan to compute run lengths
        RunIdCountPair zippedprefixsum[items_per_thread];
        {
            RunIdCountPair zipped[items_per_thread];
            #pragma unroll
            for(int i = 0; i < items_per_thread; i++){
                zipped[i].count = 1;
                zipped[i].runId = prefixsum[i];
            }
            RunIdCountPairScanOp pairScanOp{};
            BlockScanPair(temp.temp_scanPair).InclusiveScan(zipped, zippedprefixsum, pairScanOp);
            // __syncthreads();

        }


        #pragma unroll
        for(int i = 0; i < items_per_thread; i++){
            if(tailFlags[i] == 1){
                smem_runLengths[prefixsum[i]] = zippedprefixsum[i].count;
            }
        }
        __syncthreads();

// if(warp.thread_rank() == 0) printf("E\n");
        int myRunLengths[items_per_thread];
        #pragma unroll
        for(int i = 0; i < items_per_thread; i++){
            myRunLengths[i] = smem_runLengths[zippedprefixsum[i].runId];
        }

        #pragma unroll
        for(int i = 0; i < items_per_thread; i++){
            if(chunkBegin + threadIdx.x * items_per_thread + i < chunkEnd_excl){
                char letter = (myElements[i] >> 16);
                if(myRunLengths[i] >= maskingThreshold){
                    letter = maskingLetter;
                }
                smem_letters[threadIdx.x * items_per_thread + i] = letter;
            }
        }
        __syncthreads();


// if(warp.thread_rank() == 0) printf("F\n");

        //write to output

        const bool isFirstChunk = chunkBegin == 0;
        const bool isLastChunk = chunkEnd_excl == totalInputElements;

        #pragma unroll
        for(int i = 0; i < items_per_thread; i++){
            const int index = i * blockDim.x + threadIdx.x;
            if(!isFirstChunk && !isLastChunk){
                if(index >= (maskingThreshold - 1)){
                    if(index < items_per_block - (maskingThreshold - 1) ){
                        d_outputData[chunkBegin + index] = smem_letters[index];

                        // if(chunkBegin + index == 2046){
                        //     printf("AAA chunkId %ld, write %c to %d. index %d, \n", chunkId, smem_letters[index], 2046, index);
                        // }
                    }
                }
            }else if(!isFirstChunk && isLastChunk){
                if(index >= (maskingThreshold - 1)){
                    if(chunkBegin + index < chunkEnd_excl){
                        d_outputData[chunkBegin + index] = smem_letters[index];

                        // if(chunkBegin + index == 2046){
                        //     printf("BBB chunkId %ld, write %c to %d. index %d, \n", chunkId, smem_letters[index], 2046, index);
                        // }
                    }
                }
            }else if(isFirstChunk && !isLastChunk){
                if(index < items_per_block - (maskingThreshold - 1) ){
                    d_outputData[chunkBegin + index] = smem_letters[index];

                    // if(chunkBegin + index == 2046){
                    //     printf("CCC chunkId %ld, write %c to %d. index %d, \n", chunkId, smem_letters[index], 2046, index);
                    // }
                }
            }else{
                //first and last
                if(chunkBegin + index < chunkEnd_excl){
                    d_outputData[chunkBegin + index] = smem_letters[index];

                    // if(chunkBegin + index == 2046){
                    //     printf("DDD chunkId %ld, write %c to %d. index %d, \n", chunkId, smem_letters[index], 2046, index);
                    // }
                }
            }
        }


        // #pragma unroll
        // for(int i = 0; i < items_per_thread; i++){
        //     const int index = i * blockDim.x + threadIdx.x;
        //     if(chunkBegin + index < chunkEnd_excl){
        //         if(chunkId == numChunks-1){
        //             if(chunkBegin + index == printindex){
        //                 printf("AAA write %c to %d\n", smem_letters[index], printindex);
        //                 }
        //             d_outputData[chunkBegin + index] = smem_letters[index];
        //         }else{
        //             if(index < items_per_block - (maskingThreshold - 1)){
        //                 if(chunkBegin + index == printindex){
        //                     printf("BBB chunkId %ld, write %c to %d. index %d, \n", chunkId, smem_letters[index], printindex, index);
        //                 }
        //                 d_outputData[chunkBegin + index] = smem_letters[index];
        //             }
        //         }
        //     }
        // }

        // #pragma unroll
        // for(int i = 0; i < items_per_thread; i++){
        //     const int index = i * blockDim.x + threadIdx.x;
        //     if(chunkBegin + index < chunkEnd_excl){
        //         if(chunkId == 0){
        //             if(chunkBegin + index == printindex){
        //                 printf("AAA write %c to %d\n", smem_letters[index], printindex);
        //                 }
        //             d_outputData[chunkBegin + index] = smem_letters[index];
        //         }else{
        //             //avoid overwriting results from previous overlapping chunk
        //             if(index >= (maskingThreshold - 1)){
        //                 if(chunkBegin + index == printindex){
        //                     printf("BBB chunkId %ld, write %c to %d. index %d, \n", chunkId, smem_letters[index], printindex, index);
        //                 }
        //                 d_outputData[chunkBegin + index] = smem_letters[index];
        //             }
        //         }
        //     }
        // }

// if(warp.thread_rank() == 0) printf("G\n");


        __syncthreads();
    }
}




template<class InputDataIterator, class InputOffsetIterator>
void maskSequences_custom_kernel_impl3(
    char* d_outputData,
    InputDataIterator d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter,
    cudaStream_t stream
){


    // auto d_segmentIds = thrust::make_transform_iterator(
    //     thrust::make_counting_iterator<size_t>(0),
    //     cuda::proclaim_return_type<int>([=] __device__ (size_t index){
    //         auto upIter = thrust::upper_bound(thrust::seq, d_inputOffsets, d_inputOffsets + numInputSequences, index);
    //         int segmentId = numInputSequences-1;
    //         if(upIter != d_inputOffsets + numInputSequences){
    //             segmentId = cuda::std::distance(d_inputOffsets, upIter) - 1;
    //         }
    //         return segmentId;
    //     })
    // );

    //we only need to distinguish segments in each subchunk in the kernel, ushort is sufficient and provides safe wrap-around for the scan.

    uint16_t* d_segmentIds; 
    CUDACHECK(cudaMallocAsync(&d_segmentIds, sizeof(uint16_t) * totalInputElements, stream));
    CUDACHECK(cudaMemsetAsync(d_segmentIds, 0, sizeof(uint16_t) * totalInputElements, stream));
    thrust::scatter(
        thrust::cuda::par_nosync.on(stream),
        thrust::make_constant_iterator<uint16_t>(1),
        thrust::make_constant_iterator<uint16_t>(1) + numInputSequences,
        d_inputOffsets,
        d_segmentIds
    );
    thrust::inclusive_scan(
        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
        d_segmentIds,
        d_segmentIds + totalInputElements,
        d_segmentIds
    );

    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 128;
    auto kernel = maskingKernel3<blocksize, InputDataIterator, decltype(d_segmentIds)>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    kernel<<<numSMs*maxBlocksPerSM, blocksize, 0, stream>>>(
        d_outputData,
        d_inputData,
        d_segmentIds,
        totalInputElements,
        maskingThreshold,
        maskingLetter
    );
    CUDACHECKASYNC;

    setPaddingBytesKernel<128><<<cuda::ceil_div(numInputSequences, 128 / 4), 128, 0, stream>>>(
        d_outputData,
        d_inputLengths,
        d_inputOffsets,
        numInputSequences,
        paddingLetter
    );
    CUDACHECKASYNC;

    CUDACHECK(cudaFreeAsync(d_segmentIds, stream));

}




template<int blocksize, int items_per_thread, class InputDataIterator, class InputOffsetIterator>
__global__
void maskingKernel4_smem_rle_rld_1warpPerSegment(
    char* d_outputData,
    InputDataIterator d_inputData,
    const int* d_lengths,
    InputOffsetIterator d_offsets,
    int64_t totalInputElements,
    int numSegments,
    int maskingThreshold,
    char maskingLetter,
    char paddingLetter
){

    constexpr int items_per_warp = items_per_thread * 32;

    auto warp = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());
    constexpr int numWarpsInBlock = blocksize / 32;
    const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
    const int warpIdInBlock = warp.meta_group_rank();
    const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;

    constexpr char oob_element = cuda::std::numeric_limits<char>::max();

    struct RunIdCountPair{
        int runId;
        int count;
    };
    struct RunIdCountPairScanOp{
        __device__
        RunIdCountPair operator()(const RunIdCountPair& left, const RunIdCountPair& right) const{
            RunIdCountPair result;
            result.runId = right.runId;
            if(left.runId == right.runId){
                result.count = right.count + left.count;
            }else{
                result.count = right.count;                
            }
            return result;
        }
    };

    using WarpLoadChar = cub::WarpLoad<char, items_per_thread, cub::WARP_LOAD_TRANSPOSE>;
    using WarpLoadChar4 = cub::WarpLoad<char4, items_per_thread / 4, cub::WARP_LOAD_TRANSPOSE>;
    using WarpScanInt = cub::WarpScan<int>;
    using WarpScanPair = cub::WarpScan<RunIdCountPair>;

    __shared__ struct Temp{
        struct U{
            typename WarpLoadChar::TempStorage temp_warploadChar[numWarpsInBlock];
            //typename WarpLoadChar4::TempStorage temp_warploadChar4[numWarpsInBlock];
            // char4 loaded_letters4[numWarpsInBlock][items_per_warp/4];
            char masked_letters[numWarpsInBlock][items_per_warp];
        } u;
        typename WarpScanInt::TempStorage temp_warpScanInt[numWarpsInBlock];
        typename WarpScanPair::TempStorage temp_warpScanPair[numWarpsInBlock];
        int runOffsets[numWarpsInBlock][items_per_warp];
        char runValues[numWarpsInBlock][items_per_warp];
    } temp;



    for(int segmentId = warpIdInGrid; segmentId < numSegments; segmentId += numWarpsInGrid){
        const size_t segmentOffset = d_offsets[segmentId];
        const int segmentSize = d_lengths[segmentId];
        const auto segmentInputData = d_inputData + segmentOffset;
        char* const segmentOutputData = d_outputData + segmentOffset;

        const int numChunks = cuda::ceil_div(segmentSize, items_per_warp);

        char lastElementOfPreviousChunk = oob_element;
        int lastCountOfPreviousChunk = 0;
        int writtenOutputElements = 0;

        for(int chunkId = 0; chunkId < numChunks; chunkId++){
            const int chunkBegin = chunkId * items_per_warp;
            const int chunkEnd_excl = cuda::std::min(segmentSize, chunkBegin + items_per_warp);
            const int numInChunk = chunkEnd_excl - chunkBegin;

            // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
            //     printf("segmentOffset %lu, chunkBegin %d, chunkEnd_excl %d, numInChunk %d\n", segmentOffset, chunkBegin, chunkEnd_excl, numInChunk);
            // }

            alignas(4) char items[items_per_thread];

            WarpLoadChar(temp.u.temp_warploadChar[warpIdInBlock]).Load(segmentInputData + chunkBegin, items, numInChunk, oob_element);

            // {
            //     for(int i = warp.thread_rank(); i < cuda::ceil_div(numInChunk,4); i += warp.num_threads()){
            //         temp.u.loaded_letters4[warpIdInBlock][i] = ((const char4*)(segmentInputData + chunkBegin))[i];
            //     }
            //     warp.sync();
            //     #pragma unroll
            //     for(int i = 0; i < items_per_thread / 4; i++){
            //         if(warp.thread_rank() * (items_per_thread/4) + i < cuda::ceil_div(numInChunk,4)){
            //             char4 tmp = temp.u.loaded_letters4[warpIdInBlock][warp.thread_rank() * (items_per_thread/4)  + i];
            //             memcpy(&items[4*i], &tmp, 4);
            //         }else{
            //             unsigned int tmp = 0x7F7F7F7F;
            //             memcpy(&items[4*i], &tmp, 4);
            //         }
            //     }
            //     #pragma unroll
            //     for(int i = 0; i < items_per_thread; i++){
            //         if(warp.thread_rank() * (items_per_thread) + i >= numInChunk){
            //             items[i] = oob_element;
            //         }
            //     }
            // }

            // char4 items4[items_per_thread / 4];
            // WarpLoadChar4(temp.u.temp_warploadChar4[warpIdInBlock]).Load((const char4*)(segmentInputData + chunkBegin), items4, cuda::ceil_div(numInChunk,4), char4{127,127,127,127});
            // memcpy(&items[0], &items4[0], sizeof(char) * items_per_thread);

            warp.sync(); //sync to ensure smem union u can be reused

            // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
            //     printf("items\n");
            //     #pragma unroll
            //     for(int i = 0; i < items_per_thread; i++){
            //         printf("%d ", int(items[i]));
            //     }
            //     printf("\n");
            // }

            int tailFlags[items_per_thread];
            {
                auto successoritem = warp.shfl_down(items[0], 1);
                if(warp.thread_rank() == warp.num_threads()-1){
                    if(chunkEnd_excl < segmentSize){
                        successoritem = segmentInputData[chunkEnd_excl];
                    }else{
                        successoritem = oob_element;
                    }
                }
                #pragma unroll
                for(int i = 0; i < items_per_thread-1; i++){
                    tailFlags[i] = items[i] != items[i+1];
                }
                tailFlags[items_per_thread-1] = items[items_per_thread-1] != successoritem;
            }

            // if(warpIdInBlock == 0 && warp.thread_rank() == 10){
            //     printf("tailFlags\n");
            //     #pragma unroll
            //     for(int i = 0; i < items_per_thread; i++){
            //         printf("%d ", tailFlags[i]);
            //     }
            //     printf("\n");
            // }

            //exclusive scan of tail flags to compute num runs
            int prefixsum[items_per_thread];
            int warpNumRuns = 0;
            {
                prefixsum[0] = 0;
                #pragma unroll
                for(int i = 1; i < items_per_thread; i++){
                    prefixsum[i] = prefixsum[i-1] + tailFlags[i-1];
                }
                int inclTotal = prefixsum[items_per_thread-1] + tailFlags[items_per_thread-1];
                int predTotal = 0;
                WarpScanInt(temp.temp_warpScanInt[warpIdInBlock]).ExclusiveSum(inclTotal, predTotal, warpNumRuns);
                for(int i = 0; i < items_per_thread; i++){
                    prefixsum[i] += predTotal;
                }
            }

            // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
            //     printf("warpNumRuns %d\n", warpNumRuns);
            // }

            // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
            //     printf("tailFlags PS\n");
            //     #pragma unroll
            //     for(int i = 0; i < items_per_thread; i++){
            //         printf("%d ", prefixsum[i]);
            //     }
            //     printf("\n");
            // }

            //inclusive scan to compute run lengths
            int run_lengths[items_per_thread];
            {
                RunIdCountPair zippedprefixsum[items_per_thread];
                RunIdCountPair zipped[items_per_thread];
                #pragma unroll
                for(int i = 0; i < items_per_thread; i++){
                    zipped[i].count = 1;
                    zipped[i].runId = prefixsum[i];
                }
                if(chunkId > 0 && warp.thread_rank() == 0){
                    if(lastElementOfPreviousChunk == items[0]){
                        zipped[0].count = 1 + lastCountOfPreviousChunk;
                    }
                }
                RunIdCountPairScanOp pairScanOp{};

                //inclusive scan
                zippedprefixsum[0] = zipped[0];
                #pragma unroll
                for(int i = 1; i < items_per_thread; i++){
                    zippedprefixsum[i] = pairScanOp(zippedprefixsum[i-1], zipped[i]);
                }
                RunIdCountPair inclTotal = zippedprefixsum[items_per_thread-1];
                RunIdCountPair predTotal;
                WarpScanPair(temp.temp_warpScanPair[warpIdInBlock]).ExclusiveScan(inclTotal, predTotal, RunIdCountPair{0, 0}, pairScanOp);
                for(int i = 0; i < items_per_thread; i++){
                    zippedprefixsum[i] = pairScanOp(predTotal, zippedprefixsum[i]);
                }

                // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
                //     printf("zippedprefixsum count\n");
                //     #pragma unroll
                //     for(int i = 0; i < items_per_thread; i++){
                //         printf("%d ", zippedprefixsum[i].count);
                //     }
                //     printf("\n");
                // }
                // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
                //     printf("zippedprefixsum runId\n");
                //     #pragma unroll
                //     for(int i = 0; i < items_per_thread; i++){
                //         printf("%d ", zippedprefixsum[i].runId);
                //     }
                //     printf("\n");
                // }

                lastElementOfPreviousChunk = warp.shfl(items[items_per_thread-1], warp.num_threads() - 1);
                lastCountOfPreviousChunk = warp.shfl(zippedprefixsum[items_per_thread-1].count, warp.num_threads() - 1);

                #pragma unroll
                for(int i = 0; i < items_per_thread; i++){
                    run_lengths[i] = tailFlags[i] == 1 ? zippedprefixsum[i].count : 0;
                }
            }

            // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
            //     printf("run_lengths\n");
            //     #pragma unroll
            //     for(int i = 0; i < items_per_thread; i++){
            //         printf("%d ", run_lengths[i]);
            //     }
            //     printf("\n");
            // }


            //exclusive scan of run lengths
            int run_lengths_ps[items_per_thread];
            int warpNumDecodedElements = 0;            
            {
                run_lengths_ps[0] = 0;
                #pragma unroll
                for(int i = 1; i < items_per_thread; i++){
                    run_lengths_ps[i] = run_lengths_ps[i-1] + run_lengths[i-1];
                }
                int inclTotal = run_lengths_ps[items_per_thread-1] + run_lengths[items_per_thread-1];
                int predTotal = 0;
                WarpScanInt(temp.temp_warpScanInt[warpIdInBlock]).ExclusiveSum(inclTotal, predTotal, warpNumDecodedElements);
                for(int i = 0; i < items_per_thread; i++){
                    run_lengths_ps[i] += predTotal;
                }
            }

            // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
            //     printf("run_lengths excl PS\n");
            //     #pragma unroll
            //     for(int i = 0; i < items_per_thread; i++){
            //         printf("%d ", run_lengths_ps[i]);
            //     }
            //     printf("\n");
            // }


            // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
            //     printf("warpNumDecodedElements %d\n", warpNumDecodedElements);
            // }

            // if(warpIdInBlock == 0){
            //     if(segmentId == 0){

            //         warp.sync();
            //         if(warp.thread_rank() == 0){
            //             printf("warpNumDecodedElements %d\n", warpNumDecodedElements);
            //         }
            //         warp.sync();
            //         for(int t = 0; t < warp.num_threads(); t++){
            //             if(t == warp.thread_rank()){
            //                 for(int i = 0; i < items_per_thread; i++){
            //                     printf("thread %d, letter,tailFlag,prefixsum,run_lengths %d, %d, %d, %d, %d\n", 
            //                         warp.thread_rank(), int(items[i]), tailFlags[i], prefixsum[i], run_lengths[i], run_lengths_ps[i]);
            //                 }
            //             }
            //         }
            //         warp.sync();
            //     }
            // }


            //write masked RLE result to smem
            for(int i = 0; i < items_per_thread; i++){ 
                if(tailFlags[i] == 1){
                    char letter = items[i];
                    const int count = run_lengths[i];
                    if(count >= maskingThreshold){
                        letter = maskingLetter;
                    }
                    const int smemOutputOffset = prefixsum[i];
                    const int runOffset = run_lengths_ps[i];

                    temp.runOffsets[warpIdInBlock][smemOutputOffset] = runOffset; 
                    temp.runValues[warpIdInBlock][smemOutputOffset] = letter;
                    // printf("thread %d, write (%c %d) to smem.\n", warp.thread_rank(), letter, runOffset);
                }
            }
            warp.sync(); //wait for RLE results

            // if(warpIdInBlock == 0 && warp.thread_rank() == 0){
            //     if(segmentId == 0){
            //         printf("runs %d. warpNumDecodedElements %d. (value, offset)\n", warpNumRuns, warpNumDecodedElements);
            //         for(int i = 0; i < warpNumRuns; i++){
            //             printf("(%d %d), ", int(temp.runValues[warpIdInBlock][i]), temp.runOffsets[warpIdInBlock][i]);
            //         }
            //         printf("\n");
            //     }
            // }
            

            //decode RLE results into smem and write back to global output
            {
                const int* const smemOffsetsBegin = &temp.runOffsets[warpIdInBlock][0];
                const int* const smemOffsetsEnd = smemOffsetsBegin + warpNumRuns;
                const int iterations = cuda::ceil_div(warpNumDecodedElements, items_per_warp);
                for(int iteration = 0; iteration < iterations; iteration++){
                    const int iterationBeginIndex = iteration * items_per_warp;
                    const int numInIteration = cuda::std::min(warpNumDecodedElements - iterationBeginIndex, items_per_warp);

                    for(int i = warp.thread_rank(); i < numInIteration; i += warp.num_threads()){ 
                        const int outputIndex = iterationBeginIndex + i;
                        if(outputIndex < warpNumDecodedElements){
                            auto it = thrust::lower_bound(thrust::seq, smemOffsetsBegin, smemOffsetsEnd, outputIndex+1);
                            const int runId = cuda::std::distance(smemOffsetsBegin, it) - 1;
                            // if(segmentId == 42 && iteration > 0){
                            //     printf("outputIndex %d, runId %d, iterationBeginIndex %d\n", outputIndex, runId, iterationBeginIndex);
                            // }
                            temp.u.masked_letters[warpIdInBlock][i] = temp.runValues[warpIdInBlock][runId];
                        }
                    }
                    warp.sync();

                    // if(segmentId == 0){
                    //     if(warp.thread_rank() == 0){
                    //         printf("smem masked letters. numInIteration %d.\n", numInIteration);
                    //         for(int i = 0; i < numInIteration; i++){
                    //             printf("%d ", int(temp.u.masked_letters[warpIdInBlock][i]));
                    //         }
                    //         printf("\n");
                    //     }
                    // }

                    // write smem segment to global output;
                    for(int i = warp.thread_rank(); i < numInIteration; i += warp.num_threads()){
                        segmentOutputData[writtenOutputElements + i] = temp.u.masked_letters[warpIdInBlock][i];
                    }
                    warp.sync();
                    writtenOutputElements += numInIteration;
                }
            }

        }

        // if(warpIdInBlock == 0 && warp.thread_rank() == 3){
        //     printf("segmentSize %d, writtenOutputElements %d\n", segmentSize, writtenOutputElements);
        // }
        // if(segmentSize != writtenOutputElements && segmentId == 13322){
        //     if(warp.thread_rank() == 0){
        //         printf("segmentId %d, segmentSize %d, writtenOutputElements %d\n", segmentId, segmentSize, writtenOutputElements);
        //         for(int i = 0; i < segmentSize; i++){
        //             printf("%d ", segmentInputData[i]);
        //         }
        //         printf("\n");
        //     }
        // }
        assert(segmentSize == writtenOutputElements);

        //apply 4-byte padding
        const int paddedLength = cuda::ceil_div(segmentSize, 4) * 4;
        if(segmentSize + warp.thread_rank() < paddedLength){
            segmentOutputData[segmentSize + warp.thread_rank()] = paddingLetter;
        }
    }
   
}

/*
    Requirements:

*/
template<class InputDataIterator, class InputOffsetIterator>
void maskSequences_custom_kernel_impl4(
    char* d_outputData,
    InputDataIterator d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter, //to pad sequences to multiple of 4 bytes
    cudaStream_t stream
){


    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 256;
    constexpr int items_per_thread = 4;
    auto kernel = maskingKernel4_smem_rle_rld_1warpPerSegment<blocksize, items_per_thread, InputDataIterator, InputOffsetIterator>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    const int requiredNumBlocks = cuda::ceil_div(numInputSequences, blocksize / 32);
    const int numBlocks = std::min(requiredNumBlocks, numSMs*maxBlocksPerSM);

    // std::cout << "numBlocks " << numBlocks << "\n";

    kernel<<<numBlocks, blocksize, 0, stream>>>(
        d_outputData,
        d_inputData,
        d_inputLengths,
        d_inputOffsets,
        totalInputElements,
        numInputSequences,
        maskingThreshold,
        maskingLetter,
        paddingLetter
    );
    CUDACHECKASYNC;
}






template<int blocksize, int maximumChunkSize, class InputTransformOp, class InputDataIterator, class InputOffsetIterator, class OutputOffsetIterator>
__global__
void maskingKernel5_smem_leftright_1warpPerSegment(
    char* d_outputData,
    OutputOffsetIterator d_outputOffsets,
    InputDataIterator d_inputData,
    const int* d_lengths,
    InputOffsetIterator d_offsets,
    int64_t totalInputElements,
    int numSegments,
    int maskingThreshold,
    char maskingLetter,
    char paddingLetter,
    InputTransformOp inputTransformOp
){
    __builtin_assume(blocksize == blockDim.x);

    constexpr int smemElements = maximumChunkSize;
    static_assert(smemElements % 4 == 0);
    static_assert(smemElements % 16 == 0);

    auto warp = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());
    constexpr int numWarpsInBlock = blocksize / 32;
    const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
    const int warpIdInBlock = warp.meta_group_rank();
    const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;

    constexpr bool isIdentityTransform = cuda::std::is_same_v<InputTransformOp, cuda::std::identity>;


    __shared__ alignas(4) struct Smem{
        char chunkElements[numWarpsInBlock][smemElements];
    } smem;


    for(int segmentId = warpIdInGrid; segmentId < numSegments; segmentId += numWarpsInGrid){
        const size_t segmentOffset = d_offsets[segmentId];
        const int segmentSize = d_lengths[segmentId];
        const auto segmentInputData = d_inputData + segmentOffset;
        char* const segmentOutputData = d_outputData + segmentOffset;
        

        //chunks are overlapping by maskingThreshold-1 positions
        // const int numChunks = cuda::ceil_div(segmentSize, (smemElements - maskingThreshold - 1));

        // for(int chunkId = 0; chunkId < numChunks; chunkId++){
        //     const int chunkBegin = chunkId * (smemElements - maskingThreshold - 1);
        //     const int chunkEnd_excl = cuda::std::min(segmentSize, chunkBegin + smemElements);
        //     const int numInChunk = chunkEnd_excl - chunkBegin;

        for(int chunkBegin = 0; chunkBegin < segmentSize; chunkBegin += smemElements - 2*(maskingThreshold - 1)){
            const int chunkEnd_excl = cuda::std::min(segmentSize, chunkBegin + smemElements);
            const int numInChunk = chunkEnd_excl - chunkBegin;

            //load chunk to shared memory
            warp.sync(); //wait before re-using smem
            // #if defined(__CUDACC__)
            if (isIdentityTransform && cuda::std::is_pointer_v<InputDataIterator> && chunkBegin % 4 == 0){
                assert(uintptr_t(segmentInputData) % 4 == 0);
                cooperative_groups::memcpy_async(
                    warp, 
                    &smem.chunkElements[warpIdInBlock][0], 
                    cuda::aligned_size_t<4>(smemElements), 
                    segmentInputData + chunkBegin, 
                    cuda::aligned_size_t<4>(cuda::ceil_div(numInChunk,4)*4) //segments are padded to 4 bytes so reading more is safe
                );
                cooperative_groups::wait(warp);
            }else{
                for(int i = chunkBegin + warp.thread_rank(); i < chunkEnd_excl; i += warp.num_threads()){
                    smem.chunkElements[warpIdInBlock][i - chunkBegin] = inputTransformOp(segmentInputData[i]);
                }
                warp.sync();
            }
            // #else
            // for(int i = chunkBegin + warp.thread_rank(); i < chunkEnd_excl; i += warp.num_threads()){
            //     smem.chunkElements[warpIdInBlock][i - chunkBegin] = inputTransformOp(segmentInputData[i]);
            // }
            // warp.sync();
            // #endif


            if constexpr(!isIdentityTransform){
                // uint4* vectorsmem = reinterpret_cast<uint4*>(&smem.chunkElements[warpIdInBlock][0]);
                // for(int i = warp.thread_rank(); i < cuda::ceil_div(numInChunk,16); i += warp.num_threads()){
                //     uint4 tmp = vectorsmem[i];
                //     alignas(16) char letters[16];
                //     memcpy(&letters[0], &tmp, 16);

                //     #pragma unroll
                //     for(int k = 0; k < 16; k++){
                //         letters[k] = inputTransformOp(letters[k]);
                //     }

                //     memcpy(&tmp, &letters[0], 16);

                //     vectorsmem[i] = tmp;
                // }
                // warp.sync();
            }
            




            int firstOutputElementIndex = chunkBegin;
            if(firstOutputElementIndex > 0){
                firstOutputElementIndex += (maskingThreshold-1);
            }
            int lastOutputElementIndex_excl = chunkEnd_excl;
            if(lastOutputElementIndex_excl < segmentSize){
                lastOutputElementIndex_excl -= (maskingThreshold-1);
            }
            // if(segmentId == 11){
            // if(warp.thread_rank() == 0){
            //     printf("chunkBegin %d, chunkEnd_excl %d, firstOutputElementIndex %d, lastOutputElementIndex_excl %d\n", 
            //         chunkBegin, chunkEnd_excl, firstOutputElementIndex, lastOutputElementIndex_excl);
            // }
            // }

            for(int elementIndex = firstOutputElementIndex + warp.thread_rank(); elementIndex < lastOutputElementIndex_excl; elementIndex += warp.num_threads()){
                const int indexInSmem = elementIndex - chunkBegin;
                //const char myLetter = inputTransformOp(smem.chunkElements[warpIdInBlock][indexInSmem]);
                const char myLetter = (smem.chunkElements[warpIdInBlock][indexInSmem]);

                int count = 1; //account for my letter
                //check to the left
                for(int i = indexInSmem - 1; i >= cuda::std::max(0, indexInSmem - maskingThreshold + 1); i--){
                    //const char otherLetter = inputTransformOp(smem.chunkElements[warpIdInBlock][i]);
                    const char otherLetter = (smem.chunkElements[warpIdInBlock][i]);
                    if(myLetter == otherLetter){
                        count++;
                    }else{
                        break;
                    }
                }

                //check to the right
                for(int64_t i = indexInSmem + 1; i < cuda::std::min(numInChunk, indexInSmem + maskingThreshold); i++){
                    //const char otherLetter = inputTransformOp(smem.chunkElements[warpIdInBlock][i]);
                    const char otherLetter = (smem.chunkElements[warpIdInBlock][i]);
                    if(myLetter == otherLetter){
                        count++;
                        if(count >= maskingThreshold){
                            break;
                        }
                    }else{
                        break;
                    }
                }

                // if(segmentId == 11){
                //     if(elementIndex >= 10630){
                //         printf("elementIndex %d, letter %d, count %d\n", elementIndex, int(myLetter), count);
                //     }
                // }

                if(count >= maskingThreshold){
                    segmentOutputData[elementIndex] = maskingLetter;
                }else{
                    segmentOutputData[elementIndex] = myLetter;
                }
            }
        }

        //apply 4-byte padding
        const int paddedLength = cuda::ceil_div(segmentSize, 4) * 4;
        if(segmentSize + warp.thread_rank() < paddedLength){
            segmentOutputData[segmentSize + warp.thread_rank()] = paddingLetter;
        }
    }
   
}



/*
    Requirements:
        sequence[i] has length inputLengths[i] and begins at input position inputOffsets[i]
        sequence[i] is padded to 4 byte boundary
        d_outputData must not overlap d_inputData
        d_outputData must be able to contain 4-byte padded sequences
*/
template<class InputTransformOp, class InputDataIterator, class InputOffsetIterator, class OutputOffsetIterator>
void maskSequences_custom_kernel_impl5(
    char* d_outputData,
    OutputOffsetIterator d_outputOffsets,
    InputDataIterator d_inputData, 
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter, //to pad output sequences to 4 byte boundary
    InputTransformOp inputTransformOp, // input letters are given by inputTransformOp(d_inputData[i])
    cudaStream_t stream
){
    constexpr int maximumChunkSize = 1024;

    if(2*maskingThreshold >= maximumChunkSize){
        throw std::runtime_error("maskingThreshold too large for hardcoded maximumChunkSize\n");
    }

    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 256;
    auto kernel = maskingKernel5_smem_leftright_1warpPerSegment<blocksize, maximumChunkSize, 
        InputTransformOp, InputDataIterator, InputOffsetIterator, OutputOffsetIterator>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    const int requiredNumBlocks = cuda::ceil_div(numInputSequences, blocksize / 32);
    const int numBlocks = std::min(requiredNumBlocks, numSMs*maxBlocksPerSM);

    // std::cout << "numBlocks " << numBlocks << "\n";

    kernel<<<numBlocks, blocksize, 0, stream>>>(
        d_outputData,
        d_outputOffsets,
        d_inputData,
        d_inputLengths,
        d_inputOffsets,
        totalInputElements,
        numInputSequences,
        maskingThreshold,
        maskingLetter,
        paddingLetter,
        inputTransformOp
    );
    CUDACHECKASYNC;
}






template<int items_per_thread, int maxMaskingThreshold, class InputDataIterator, class InputOffsetIterator>
__global__
void maskingKernel6_reg_leftright_perSegment(
    char* d_outputData,
    InputDataIterator d_inputData,
    const int* d_lengths,
    InputOffsetIterator d_offsets,
    int64_t totalInputElements,
    int numSegments,
    int maskingThreshold,
    char maskingLetter,
    char paddingLetter
){
    static_assert(maxMaskingThreshold <= 32); //this ensures that at most two item indices need to be accessed for each search iteration

    constexpr int items_per_warp = 32 * items_per_thread;

    auto warp = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());
    const int numWarpsInBlock = blockDim.x / 32;
    const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
    // const int warpIdInBlock = warp.meta_group_rank();
    const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;

    constexpr int oob_element = cuda::std::numeric_limits<int>::max(); //use masking letter instead?


    for(int segmentId = warpIdInGrid; segmentId < numSegments; segmentId += numWarpsInGrid){
        const size_t segmentOffset = d_offsets[segmentId];
        const int segmentSize = d_lengths[segmentId];
        const auto segmentInputData = d_inputData + segmentOffset;
        char* const segmentOutputData = d_outputData + segmentOffset;

        for(int chunkBegin = 0; chunkBegin < segmentSize; chunkBegin += items_per_warp - 2*(maskingThreshold - 1)){
            const int chunkEnd_excl = cuda::std::min(segmentSize, chunkBegin + items_per_warp);
            // const int numInChunk = chunkEnd_excl - chunkBegin;

            //load chunk to registers

            int items[items_per_thread];
            
            #pragma unroll
            for(int i = 0; i < items_per_thread; i++){
                //striped index
                const int index = chunkBegin + i * warp.num_threads() + warp.thread_rank();
                if(index < chunkEnd_excl){
                    items[i] = segmentInputData[index];
                }else{
                    items[i] = oob_element;
                }
            }

            int firstOutputElementIndex = chunkBegin;
            if(firstOutputElementIndex > 0){
                firstOutputElementIndex += (maskingThreshold-1);
            }
            int lastOutputElementIndex_excl = chunkEnd_excl;
            if(lastOutputElementIndex_excl < segmentSize){
                lastOutputElementIndex_excl -= (maskingThreshold-1);
            }

            //scan surrounding elements to check if run lenght meets threshold, then write to output
            #pragma unroll
            for(int currentIndex = 0; currentIndex < items_per_thread; currentIndex++){
                //striped index
                const int outputindex = chunkBegin + currentIndex * warp.num_threads() + warp.thread_rank();

                //check left and right of currentIndex to count occurences
                int myLetter = items[currentIndex];
                int count = 1;
                
                //check left
                int leftNeighborTid = warp.thread_rank() == 0 ? warp.num_threads() - 1 : warp.thread_rank() - 1;
                bool stopLeft = false;
                #pragma unroll
                for(int distance = 1; distance < maxMaskingThreshold; distance++){
                    if(distance < maskingThreshold){
                        const int neighborLetterFromPrevIndex = currentIndex == 0 ? oob_element : warp.shfl(items[currentIndex-1], leftNeighborTid);
                        const int neighborLetterFromSameIndex = warp.shfl(items[currentIndex], leftNeighborTid);
                        const int neighborLetter = leftNeighborTid >= warp.thread_rank() ? neighborLetterFromPrevIndex : neighborLetterFromSameIndex;
                        // if(segmentId == 1029){
                        //     if(outputindex == 55){
                        //         printf("left. chunkBegin %d, my thread id %d, currentIndex %d, distance %d, myLetter %d, neighborLetter %d, leftNeighborTid %d\n",
                        //         chunkBegin, warp.thread_rank(), currentIndex, distance, myLetter, neighborLetter, leftNeighborTid);
                        //     }
                        // }
                        if(neighborLetter == myLetter){
                            if(!stopLeft){
                                count++;
                            }
                        }else{
                            //run which includes myLetter is finished
                            stopLeft = true; //thread is done, but we cannot break out of the loop because shfl needs to be executed by full warp
                        }
                        leftNeighborTid = leftNeighborTid == 0 ? warp.num_threads() - 1 : leftNeighborTid - 1;
                    }
                }

                //check right
                int rightNeighborTid = (warp.thread_rank() + 1) % warp.num_threads();
                bool stopRight = false;
                for(int distance = 1; distance < maxMaskingThreshold; distance++){
                    if(distance < maskingThreshold){
                        const int neighborLetterFromNextIndex = currentIndex == items_per_thread-1 ? oob_element : warp.shfl(items[currentIndex+1], rightNeighborTid);
                        const int neighborLetterFromSameIndex = warp.shfl(items[currentIndex], rightNeighborTid);
                        const int neighborLetter = rightNeighborTid <= warp.thread_rank() ? neighborLetterFromNextIndex : neighborLetterFromSameIndex;
                        // if(segmentId == 1029){
                        //     if(outputindex == 55){
                        //         printf("right. chunkBegin %d, my thread id %d, currentIndex %d, distance %d, myLetter %d, neighborLetter %d, rightNeighborTid %d\n",
                        //         chunkBegin, warp.thread_rank(), currentIndex, distance, myLetter, neighborLetter, rightNeighborTid);
                        //     }
                        // }
                        if(neighborLetter == myLetter){
                            if(!stopRight){
                                count++;
                            }
                        }else{
                            //run which includes myLetter is finished
                            stopRight = true; //thread is done, but we cannot break out of the loop because shfl needs to be executed by full warp
                        }
                        rightNeighborTid = (rightNeighborTid + 1) % warp.num_threads();
                    }
                }


                if(firstOutputElementIndex <= outputindex && outputindex < lastOutputElementIndex_excl){
                    if(count >= maskingThreshold){
                        segmentOutputData[outputindex] = maskingLetter;
                    }else{
                        segmentOutputData[outputindex] = myLetter;
                    }
                }

            }
        }

        //apply 4-byte padding
        const int paddedLength = cuda::ceil_div(segmentSize, 4) * 4;
        if(segmentSize + warp.thread_rank() < paddedLength){
            segmentOutputData[segmentSize + warp.thread_rank()] = paddingLetter;
        }
    }
   
}

/*
    Requirements:

*/
template<int maxMaskingThreshold, class InputDataIterator, class InputOffsetIterator>
void maskSequences_custom_kernel_impl6(
    char* d_outputData,
    InputDataIterator d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter, //to pad sequences to multiple of 4 bytes
    cudaStream_t stream
){
    if(maskingThreshold > maxMaskingThreshold){
        throw std::runtime_error("maskingThreshold too large for hardcoded maxMaskingThreshold\n");
    }

    constexpr int items_per_thread = 4;

    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 256;
    auto kernel = maskingKernel6_reg_leftright_perSegment<items_per_thread, maxMaskingThreshold, InputDataIterator, InputOffsetIterator>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    const int requiredNumBlocks = cuda::ceil_div(numInputSequences, blocksize / 32);
    const int numBlocks = std::min(requiredNumBlocks, numSMs*maxBlocksPerSM);

    // std::cout << "numBlocks " << numBlocks << "\n";

    kernel<<<numBlocks, blocksize, 0, stream>>>(
        d_outputData,
        d_inputData,
        d_inputLengths,
        d_inputOffsets,
        totalInputElements,
        numInputSequences,
        maskingThreshold,
        maskingLetter,
        paddingLetter
    );
    CUDACHECKASYNC;
}



#if !defined(FOOBARAOS)


/*
    In constrast to maskingKernel5 which splits input sequences into chunks, maskingKernel7 splits output sequences into chunks
*/
template<int blocksize, int outputChunkSize, int maxMaskingThreshold, class InputTransformOp, class InputDataIterator, class InputOffsetIterator, class OutputOffsetIterator>
__global__
void maskingKernel7_smem_leftright_1warpPerSegment(
    char* d_outputData,
    OutputOffsetIterator d_outputOffsets,
    InputDataIterator d_inputData,
    const int* d_lengths,
    InputOffsetIterator d_offsets,
    int64_t totalInputElements,
    int numSegments,
    int maskingThreshold,
    char maskingLetter,
    char paddingLetter,
    InputTransformOp inputTransformOp
){
    __builtin_assume(blocksize == blockDim.x);
    __builtin_assume(maskingThreshold <= maxMaskingThreshold);

    static_assert(outputChunkSize % 16 == 0);
    constexpr int smemNumInputElements = outputChunkSize + 2 * (cuda::ceil_div((maxMaskingThreshold-1), 16) * 16);
    // constexpr int smemNumOutputElements = outputChunkSize;

    auto warp = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());
    constexpr int numWarpsInBlock = blocksize / 32;
    const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
    const int warpIdInBlock = warp.meta_group_rank();
    const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;

    constexpr bool isIdentityTransform = cuda::std::is_same_v<InputTransformOp, cuda::std::identity>;

    constexpr int numStages = 2;
    
    __shared__ alignas(16) struct Smem{
        char inputElements[numWarpsInBlock][numStages][smemNumInputElements];
        //char outputElements[numWarpsInBlock][smemNumOutputElements];
    } smem;


    for(int segmentId = warpIdInGrid; segmentId < numSegments; segmentId += numWarpsInGrid){

        const size_t segmentOffset = d_offsets[segmentId];
        const int segmentSize = d_lengths[segmentId];
        const auto segmentInputData = d_inputData + segmentOffset;
        char* const segmentOutputData = d_outputData + segmentOffset;
        const int numOutputChunks = cuda::ceil_div(segmentSize, outputChunkSize);

        for(int outputChunk = 0, fetchChunk = 0; outputChunk < numOutputChunks; outputChunk++){
            for (; fetchChunk < numOutputChunks && fetchChunk < (outputChunk + numStages); ++fetchChunk) {
                const int outputChunkBegin = fetchChunk * outputChunkSize;
                const int outputChunkEnd_excl = cuda::std::min(segmentSize, outputChunkBegin + outputChunkSize);
                // const int numInOutputChunk = outputChunkEnd_excl - outputChunkBegin;
    
                const int inputChunkBegin = cuda::std::max(0, outputChunkBegin - (maskingThreshold-1));
                const int inputChunkEnd_excl = cuda::std::min(segmentSize, outputChunkEnd_excl + (maskingThreshold-1));
    
                //align load begin to 4 bytes. valid since each segment begin is aligned to 4 bytes.
                const int inputLoadFrontPadding = inputChunkBegin % 4;
                const int inputChunkLoadBegin = inputChunkBegin - inputLoadFrontPadding;
                //align load end to 4 bytes. valid since all segments are padded to 4 bytes
                const int inputChunkLoadEnd_excl = cuda::ceil_div(inputChunkEnd_excl, 4) * 4;
                const int loadSize = inputChunkLoadEnd_excl - inputChunkLoadBegin;
                // assert(loadSize <= smemNumInputElements);

                // if(segmentId == 0){
                //     if(warp.thread_rank() == 0){
                //         printf("fetchChunk %d, load %d elements, %d - %d, smemNumInputElements %d\n", 
                //                 fetchChunk,loadSize, inputChunkLoadBegin,inputChunkLoadEnd_excl, smemNumInputElements);
                //         for(int i = 0; i < 16; i++){
                //             printf("%d ", int((segmentInputData + inputChunkLoadBegin)[i]));
                //         }
                //         printf("\n");
                //     }
                // }
                cooperative_groups::memcpy_async(
                    warp, 
                    &smem.inputElements[warpIdInBlock][fetchChunk % numStages][0], 
                    cuda::aligned_size_t<4>(smemNumInputElements), 
                    segmentInputData + inputChunkLoadBegin, 
                    cuda::aligned_size_t<4>(loadSize)
                );
            }
            const int outputChunkBegin = outputChunk * outputChunkSize;
            const int outputChunkEnd_excl = cuda::std::min(segmentSize, outputChunkBegin + outputChunkSize);
            const int numInOutputChunk = outputChunkEnd_excl - outputChunkBegin;

            const int inputChunkBegin = cuda::std::max(0, outputChunkBegin - (maskingThreshold-1));
            const int inputChunkEnd_excl = cuda::std::min(segmentSize, outputChunkEnd_excl + (maskingThreshold-1));

            //align load begin to 4 bytes. valid since each segment begin is aligned to 4 bytes.
            const int inputLoadFrontPadding = inputChunkBegin % 4;
            const int inputChunkLoadBegin = inputChunkBegin - inputLoadFrontPadding;
            //align load end to 4 bytes. valid since all segments are padded to 4 bytes
            const int inputChunkLoadEnd_excl = cuda::ceil_div(inputChunkEnd_excl, 4) * 4;
            const int loadSize = inputChunkLoadEnd_excl - inputChunkLoadBegin;




            if constexpr(!isIdentityTransform){
                for(int i = 4 * warp.thread_rank(); i < loadSize; i += 4*warp.num_threads()){
                    char4 four = *reinterpret_cast<const char4*>(&smem.inputElements[warpIdInBlock][outputChunk % numStages][i]);
                    four.x = inputTransformOp(four.x);
                    four.y = inputTransformOp(four.y);
                    four.z = inputTransformOp(four.z);
                    four.w = inputTransformOp(four.w);
                    *reinterpret_cast<char4*>(&smem.inputElements[warpIdInBlock][outputChunk % numStages][i]) = four;
                }
            }

            const int inputOutputSmemOffset = outputChunkBegin - inputChunkLoadBegin;
            const int inputSmemEndIndex_excl = inputChunkEnd_excl - inputChunkLoadBegin;
            for(int outputElementLocalIndex = warp.thread_rank(); outputElementLocalIndex < numInOutputChunk; outputElementLocalIndex += warp.num_threads()){

                const int inputIndexInSmem = inputOutputSmemOffset + outputElementLocalIndex;
                const char myLetter = (smem.inputElements[warpIdInBlock][outputChunk % numStages][inputIndexInSmem]);

                // if(segmentId == 0){
                //     if(outputChunk == 1){
                //         if(warp.thread_rank() == 0){
                //             printf("outputElementLocalIndex %d, inputIndexInSmem %d, myLetter %d\n", 
                //                 outputElementLocalIndex, inputIndexInSmem, myLetter);
                //         }
                //     }
                // }

                int count = 1; //account for my letter
                //check to the left
                for(int i = inputIndexInSmem - 1; i >= cuda::std::max(0, inputIndexInSmem - maskingThreshold + 1); i--){
                    //const char otherLetter = inputTransformOp(smem.inputElements[warpIdInBlock][i]);
                    const char otherLetter = (smem.inputElements[warpIdInBlock][outputChunk % numStages][i]);
                    if(myLetter == otherLetter){
                        count++;
                    }else{
                        break;
                    }
                }

                //check to the right
                const int loopLimit = cuda::std::min(inputSmemEndIndex_excl, inputIndexInSmem + maskingThreshold);
                for(int64_t i = inputIndexInSmem + 1; i < loopLimit; i++){                    
                    //const char otherLetter = inputTransformOp(smem.chunkElements[warpIdInBlock][i]);
                    const char otherLetter = (smem.inputElements[warpIdInBlock][outputChunk % numStages][i]);
                    if(myLetter == otherLetter){
                        count++;
                        if(count >= maskingThreshold){
                            break;
                        }
                    }else{
                        break;
                    }
                }

                // if(warp.thread_rank() == 0){
                //     if(segmentId == 0){
                //         printf("write to gmem pos %d, count %d\n", outputChunkBegin + outputElementLocalIndex, count);
                //     }
                // }

                const char outputLetter = (count >= maskingThreshold) ? maskingLetter : myLetter;
                segmentOutputData[outputChunkBegin + outputElementLocalIndex] = outputLetter;
            }
            
            warp.sync(); //done with shared memory of stage, can be reused for prefetching
        }
        //apply 4-byte padding
        const int paddedLength = cuda::ceil_div(segmentSize, 4) * 4;
        if(segmentSize + warp.thread_rank() < paddedLength){
            segmentOutputData[segmentSize + warp.thread_rank()] = paddingLetter;
        }

    }
   
}



/*
    Requirements:
        sequence[i] has length inputLengths[i] and begins at input position inputOffsets[i]
        sequence[i] is padded to 4 byte boundary
        d_outputData must not overlap d_inputData
        d_outputData must be able to contain 4-byte padded sequences
*/
template<int outputChunkSize, int maxMaskingThreshold, class InputTransformOp, class InputDataIterator, class InputOffsetIterator, class OutputOffsetIterator>
void maskSequences_custom_kernel_impl7(
    char* d_outputData,
    OutputOffsetIterator d_outputOffsets,
    InputDataIterator d_inputData, 
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter, //to pad output sequences to 4 byte boundary
    InputTransformOp inputTransformOp, // input letters are given by inputTransformOp(d_inputData[i])
    cudaStream_t stream
){
    if(maskingThreshold > maxMaskingThreshold){
        throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
    }
    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 256;
    auto kernel = maskingKernel7_smem_leftright_1warpPerSegment<blocksize, outputChunkSize, maxMaskingThreshold,  
        InputTransformOp, InputDataIterator, InputOffsetIterator, OutputOffsetIterator>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    const int requiredNumBlocks = cuda::ceil_div(numInputSequences, blocksize / 32);
    const int numBlocks = std::min(requiredNumBlocks, numSMs*maxBlocksPerSM);

    // std::cout << "numBlocks " << numBlocks << "\n";

    kernel<<<numBlocks, blocksize, 0, stream>>>(
        d_outputData,
        d_outputOffsets,
        d_inputData,
        d_inputLengths,
        d_inputOffsets,
        totalInputElements,
        numInputSequences,
        maskingThreshold,
        maskingLetter,
        paddingLetter,
        inputTransformOp
    );
    CUDACHECKASYNC;
}

#endif







template<int blocksize, int /*outputChunkSize*/, int maxMaskingThreshold, class InputTransformOp, class InputDataIterator>
__global__
void maskingKernel8_gmem_leftright_consecutivesegments(
    char* d_outputData,
    InputDataIterator d_inputData,
    const uint32_t* d_bitarray,
    int64_t totalInputElements,
    int numSegments,
    int maskingThreshold,
    char maskingLetter,
    char paddingLetter,
    InputTransformOp inputTransformOp
){
    __builtin_assume(blocksize == blockDim.x);
    __builtin_assume(maskingThreshold <= maxMaskingThreshold);

    constexpr bool isIdentityTransform = cuda::std::is_same_v<InputTransformOp, cuda::std::identity>;

    auto bitarray_isSet = [&](int64_t position){
        const int64_t dataIndex = position / (8 * sizeof(uint32_t));
        const int64_t bitIndex = position % (8 * sizeof(uint32_t));

        return (d_bitarray[dataIndex] >> bitIndex) & 1;
    };

    for(int64_t myElementIndex = threadIdx.x + blockIdx.x * blockDim.x; myElementIndex < totalInputElements; myElementIndex += int64_t(blockDim.x) * gridDim.x){
        const char myLetter = isIdentityTransform ? d_inputData[myElementIndex] : inputTransformOp(d_inputData[myElementIndex]);
        const bool myFlag = bitarray_isSet(myElementIndex);

        int count = 1; //account for my letter


        //check to the left
        //comparing to the left must not exceed a flagged position
        if(!myFlag){
            for(int i = myElementIndex - 1; i >= cuda::std::max(int64_t(0), myElementIndex - maskingThreshold + 1); i--){
                const char otherLetter = isIdentityTransform ? d_inputData[i] : inputTransformOp(d_inputData[i]);
                const bool otherFlag = bitarray_isSet(i);
                if(myLetter == otherLetter){
                    count++;
                    // do not continue to the left.
                    if(otherFlag){
                        break;
                    }
                }else{
                    break;
                }
            }
        }

        //check to the right
        //comparing to the right must stop before including a flagged position
        const int loopLimit = cuda::std::min(totalInputElements, myElementIndex + maskingThreshold);
        for(int64_t i = myElementIndex + 1; i < loopLimit; i++){
            const char otherLetter = isIdentityTransform ? d_inputData[i] : inputTransformOp(d_inputData[i]);
            const bool otherFlag = bitarray_isSet(i);
            if(myLetter == otherLetter && !otherFlag){
                count++;
                if(count >= maskingThreshold){
                    break;
                }
            }else{
                break;
            }
        }

        const char outputLetter = (count >= maskingThreshold) ? maskingLetter : myLetter;
        d_outputData[myElementIndex] = outputLetter;
    }

                

   
}





template<int outputChunkSize, int maxMaskingThreshold, class InputTransformOp, class InputDataIterator>
void maskSequences_custom_kernel_impl8(
    char* d_outputData,
    InputDataIterator d_inputData, 
    const uint32_t* d_bitarray,
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter, //to pad output sequences to 4 byte boundary
    InputTransformOp inputTransformOp, // input letters are given by inputTransformOp(d_inputData[i])
    cudaStream_t stream
){
    if(maskingThreshold > maxMaskingThreshold){
        throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
    }
    if(maskingThreshold < 4){
        throw std::runtime_error("masking Threshold must be greater than the maximum number of padding bytes per sequence");
    }
    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 256;
    auto kernel = maskingKernel8_gmem_leftright_consecutivesegments<blocksize, outputChunkSize, maxMaskingThreshold,  
        InputTransformOp, InputDataIterator>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    const int requiredNumBlocks = cuda::ceil_div(numInputSequences, blocksize / 32);
    const int numBlocks = std::min(requiredNumBlocks, numSMs*maxBlocksPerSM);

    // std::cout << "numBlocks " << numBlocks << "\n";

    kernel<<<numBlocks, blocksize, 0, stream>>>(
        d_outputData,
        d_inputData,
        d_bitarray,
        totalInputElements,
        numInputSequences,
        maskingThreshold,
        maskingLetter,
        paddingLetter,
        inputTransformOp
    );
    CUDACHECKASYNC;
}


#if !defined(FOOBARAOS)

template<int blocksize, int outputChunkSize, int maxMaskingThreshold, class InputTransformOp, class InputDataIterator>
__global__
void maskingKernel9_smem_leftright_consecutivesegments(
    char* d_outputData,
    InputDataIterator d_inputData,
    const uint32_t* d_bitarray,
    int64_t totalInputElements,
    int numSegments,
    int maskingThreshold,
    char maskingLetter,
    char paddingLetter,
    InputTransformOp inputTransformOp
){
    __builtin_assume(blocksize == blockDim.x);
    __builtin_assume(maskingThreshold <= maxMaskingThreshold);

    static_assert(outputChunkSize % 16 == 0);
    constexpr int smemNumInputElements = outputChunkSize + 2 * (cuda::ceil_div((maxMaskingThreshold-1), 16) * 16);
    constexpr int smemBitArrayElements = smemNumInputElements;
    // constexpr int smemNumOutputElements = outputChunkSize;

    auto warp = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());
    constexpr int numWarpsInBlock = blocksize / 32;
    const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
    const int warpIdInBlock = warp.meta_group_rank();
    const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;

    constexpr bool isIdentityTransform = cuda::std::is_same_v<InputTransformOp, cuda::std::identity>;

    constexpr int numStages = 2;
    
    __shared__ alignas(16) struct Smem{
        char inputElements[numWarpsInBlock][numStages][smemNumInputElements];
        char flagArray[numWarpsInBlock][numStages][smemBitArrayElements];
    } smem;

    auto bitarray_isSet = [&](int64_t position){
        const int64_t dataIndex = position / (8 * sizeof(uint32_t));
        const int64_t bitIndex = position % (8 * sizeof(uint32_t));

        return (d_bitarray[dataIndex] >> bitIndex) & 1;
    };

    const int totalInputElements_int = totalInputElements;

    const int numOutputChunks = cuda::ceil_div(totalInputElements_int, outputChunkSize);

    for(int outputChunk = warpIdInGrid, fetchChunk = warpIdInGrid, outputStage = 0, fetchStage = 0; outputChunk < numOutputChunks; outputChunk += numWarpsInGrid, outputStage++){
        for (; fetchChunk < numOutputChunks && fetchChunk < (outputChunk + numStages * numWarpsInGrid); fetchChunk += numWarpsInGrid, fetchStage++) {

            const int outputChunkBegin = fetchChunk * outputChunkSize;
            const int outputChunkEnd_excl = cuda::std::min(totalInputElements_int, outputChunkBegin + outputChunkSize);
            // const int numInOutputChunk = outputChunkEnd_excl - outputChunkBegin;

            const int inputChunkBegin = cuda::std::max(0, outputChunkBegin - (maskingThreshold-1));
            const int inputChunkEnd_excl = cuda::std::min(totalInputElements_int, outputChunkEnd_excl + (maskingThreshold-1));

            //align load begin to 4 bytes. valid since each segment begin is aligned to 4 bytes.
            const int inputLoadFrontPadding = inputChunkBegin % 4;
            const int inputChunkLoadBegin = inputChunkBegin - inputLoadFrontPadding;
            //align load end to 4 bytes. valid since all segments are padded to 4 bytes
            const int inputChunkLoadEnd_excl = cuda::ceil_div(inputChunkEnd_excl, 4) * 4;
            const int loadSize = inputChunkLoadEnd_excl - inputChunkLoadBegin;
            // if(warp.thread_rank()==0){
            //     printf("warpIdInGrid %d, outputChunk %d, outputStage %d, fetchChunk %d, fetchStage %d, dst ptr %p, smemNumInputElements %d, loadSize %d, inputChunkLoadBegin %d\n",
            //         warpIdInGrid, outputChunk, outputStage, fetchChunk, fetchStage, 
            //         &smem.inputElements[warpIdInBlock][fetchStage % numStages][0], smemNumInputElements, loadSize, inputChunkLoadBegin
            //     );
            // }

            cooperative_groups::memcpy_async(
                warp, 
                &smem.inputElements[warpIdInBlock][fetchStage % numStages][0], 
                cuda::aligned_size_t<4>(smemNumInputElements), 
                d_inputData + inputChunkLoadBegin, 
                cuda::aligned_size_t<4>(loadSize)
            );
        }

        const int outputChunkBegin = outputChunk * outputChunkSize;
        const int outputChunkEnd_excl = cuda::std::min(totalInputElements_int, outputChunkBegin + outputChunkSize);
        const int numInOutputChunk = outputChunkEnd_excl - outputChunkBegin;

        const int inputChunkBegin = cuda::std::max(0, outputChunkBegin - (maskingThreshold-1));
        const int inputChunkEnd_excl = cuda::std::min(totalInputElements_int, outputChunkEnd_excl + (maskingThreshold-1));

        //align load begin to 4 bytes. valid since each segment begin is aligned to 4 bytes.
        const int inputLoadFrontPadding = inputChunkBegin % 4;
        const int inputChunkLoadBegin = inputChunkBegin - inputLoadFrontPadding;
        //align load end to 4 bytes. valid since all segments are padded to 4 bytes
        const int inputChunkLoadEnd_excl = cuda::ceil_div(inputChunkEnd_excl, 4) * 4;
        const int loadSize = inputChunkLoadEnd_excl - inputChunkLoadBegin;

        for(int i = warp.thread_rank(); i < loadSize; i += warp.num_threads()){
            smem.flagArray[warpIdInBlock][outputStage % numStages][i] = bitarray_isSet(inputChunkLoadBegin + i);
        }

        //wait for async load from gmem to smem of our input chunk
        if(outputChunk < numOutputChunks - 1 - int(numWarpsInGrid)){
            cooperative_groups::wait_prior<1>(warp);
            // cooperative_groups::wait(warp);
        }else{
            cooperative_groups::wait(warp);
        }
        //and wait for smem.flagArray
        // warp.sync();

        if constexpr(!isIdentityTransform){
            for(int i = 4 * warp.thread_rank(); i < loadSize; i += 4*warp.num_threads()){
                char4 four = *reinterpret_cast<const char4*>(&smem.inputElements[warpIdInBlock][outputStage % numStages][i]);
                four.x = inputTransformOp(four.x);
                four.y = inputTransformOp(four.y);
                four.z = inputTransformOp(four.z);
                four.w = inputTransformOp(four.w);

                // half2 parts[2];
                // parts[0] = make_half2(four.x, four.y);
                // parts[1] = make_half2(four.z, four.w);
            
                // parts[0] = parts[0] * make_half2(1/12.f, 1/12.f);
                // parts[1] = parts[1] * make_half2(1/12.f, 1/12.f);

                // four.x = parts[0].x;
                // four.y = parts[0].y;
                // four.z = parts[1].x;
                // four.w = parts[1].y;

                *reinterpret_cast<char4*>(&smem.inputElements[warpIdInBlock][outputStage % numStages][i]) = four;
            }
            warp.sync();
        }
        auto maybeTransform = [&](auto x){
            return x;
            // return inputTransformOp(x);
        };

        const int inputOutputSmemOffset = outputChunkBegin - inputChunkLoadBegin;
        const int inputSmemEndIndex_excl = inputChunkEnd_excl - inputChunkLoadBegin;
        for(int outputElementLocalIndex = warp.thread_rank(); outputElementLocalIndex < numInOutputChunk; outputElementLocalIndex += warp.num_threads()){

            const int inputIndexInSmem = inputOutputSmemOffset + outputElementLocalIndex;
            // const char myLetter = (smem.inputElements[warpIdInBlock][outputChunk % numStages][inputIndexInSmem]);
            const char myLetter = maybeTransform(smem.inputElements[warpIdInBlock][outputStage % numStages][inputIndexInSmem]);
            const char myFlag = smem.flagArray[warpIdInBlock][outputStage % numStages][inputIndexInSmem];


            int count = 1; //account for my letter
            //check to the left
            //comparing to the left must not exceed a flagged position
            if(!myFlag){
                for(int i = inputIndexInSmem - 1; i >= cuda::std::max(0, inputIndexInSmem - maskingThreshold + 1); i--){
                    const char otherLetter = maybeTransform(smem.inputElements[warpIdInBlock][outputStage % numStages][i]);

                    const char otherFlag = smem.flagArray[warpIdInBlock][outputStage % numStages][i];
                    if(myLetter == otherLetter){
                        count++;
                        // do not continue to the left.
                        if(otherFlag == 1){
                            break;
                        }
                    }else{
                        break;
                    }
                }
            }

            //check to the right
            //comparing to the right must stop before including a flagged position
            const int loopLimit = cuda::std::min(inputSmemEndIndex_excl, inputIndexInSmem + maskingThreshold);
            for(int64_t i = inputIndexInSmem + 1; i < loopLimit; i++){                    
                const char otherLetter = maybeTransform(smem.inputElements[warpIdInBlock][outputStage % numStages][i]);

                const char otherFlag = smem.flagArray[warpIdInBlock][outputStage % numStages][i];
                if(myLetter == otherLetter && otherFlag == 0){
                    count++;
                    if(count >= maskingThreshold){
                        break;
                    }
                }else{
                    break;
                }
            }

            const char outputLetter = (count >= maskingThreshold) ? maskingLetter : myLetter;
            d_outputData[outputChunkBegin + outputElementLocalIndex] = outputLetter;
        }
        
        // if(warp.thread_rank() == 0){
        //     printf("warpIdInGrid %d, outputChunk %d, outputStage %d, chunk end sync\n",
        //                 warpIdInGrid, outputChunk, outputStage
        //             );
        // }
        warp.sync(); //done with shared memory of stage, can be reused for prefetching
    }
}



template<int outputChunkSize, int maxMaskingThreshold, class InputTransformOp, class InputDataIterator>
void maskSequences_custom_kernel_impl9(
    char* d_outputData,
    InputDataIterator d_inputData, 
    const uint32_t* d_bitarray,
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter, //to pad output sequences to 4 byte boundary
    InputTransformOp inputTransformOp, // input letters are given by inputTransformOp(d_inputData[i])
    cudaStream_t stream
){
    if(maskingThreshold > maxMaskingThreshold){
        throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
    }
    if(maskingThreshold < 4){
        throw std::runtime_error("masking Threshold must be greater than the maximum number of padding bytes per sequence");
    }
    if(totalInputElements > size_t(std::numeric_limits<int>::max())){
        throw std::runtime_error("> 2GB inputs are not supported for masking kernel 9");
    }

    // totalInputElements = 1024;

    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 128;
    auto kernel = maskingKernel9_smem_leftright_consecutivesegments<blocksize, outputChunkSize, maxMaskingThreshold,  
        InputTransformOp, InputDataIterator>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    const int requiredNumBlocks = cuda::ceil_div(numInputSequences, blocksize / 32);
    const int numBlocks = std::min(requiredNumBlocks, numSMs*maxBlocksPerSM);

    // std::cout << "numBlocks " << numBlocks << "\n";

    kernel<<<numBlocks, blocksize, 0, stream>>>(
        d_outputData,
        d_inputData,
        d_bitarray,
        totalInputElements,
        numInputSequences,
        maskingThreshold,
        maskingLetter,
        paddingLetter,
        inputTransformOp
    );
    CUDACHECKASYNC;
}


#endif



template<class InputTransformOp>
__global__
void maskingComputeHeadFlagsKernel1(
    //uint32_t* bitarrayFlags,
    char* __restrict__ flags,
    const char* __restrict__ d_inputData,
    int64_t totalInputElements,
    InputTransformOp inputTransformOp
){
    //const int64_t chunksOf16 = totalInputElements / 16;

    const int tid = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = blockDim.x * gridDim.x;

    const char* inputAlignedPtr = cuda::std::assume_aligned<16>(d_inputData);
    char* flagsAlignedPtr = cuda::std::assume_aligned<16>(flags);
    if(tid == 0){
        flagsAlignedPtr[0] = 1;
    }
    for(int64_t i = 1 + tid; i < totalInputElements; i += stride){
        const char left = inputAlignedPtr[i-1];
        const char right = inputAlignedPtr[i];
        const char leftTransformed = inputTransformOp(left);
        const char rightTransformed = inputTransformOp(right);

        if(leftTransformed != rightTransformed){
            flagsAlignedPtr[i] = 1;
        }else{
            flagsAlignedPtr[i] = 0;
        }
    }
}



template<class InputTransformOp>
__global__
void maskingComputeHeadFlagsKernel2(
    //uint32_t* bitarrayFlags,
    char* __restrict__ flags,
    const char* __restrict__ d_inputData,
    int64_t totalInputElements,
    InputTransformOp inputTransformOp
){
    const int64_t chunksOf16 = totalInputElements / 16;

    const int tid = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = blockDim.x * gridDim.x;

    const char* inputAlignedPtr = cuda::std::assume_aligned<16>(d_inputData);
    char* flagsAlignedPtr = cuda::std::assume_aligned<16>(flags);
    if(tid == 0){
        flagsAlignedPtr[0] = 1;
    }

    //first chunk of 16 (which does not have a left neighbor chunk)
    for(int64_t i = 1 + tid; i < cuda::std::min(int64_t(16), totalInputElements); i += stride){
        const char left = inputAlignedPtr[i-1];
        const char right = inputAlignedPtr[i];
        const char leftTransformed = inputTransformOp(left);
        const char rightTransformed = inputTransformOp(right);

        if(leftTransformed != rightTransformed){
            flagsAlignedPtr[i] = 1;
        }else{
            flagsAlignedPtr[i] = 0;
        }
    }

    for(int64_t c = 1 + tid; c < chunksOf16; c += stride){
        const float4 leftChunkTmp = ((const float4*)inputAlignedPtr)[c-1];
        const float4 rightChunkTmp = ((const float4*)inputAlignedPtr)[c];
        alignas(16) char leftChunk[16];
        memcpy(&leftChunk[0], &leftChunkTmp, 16);
        alignas(16) char rightChunk[16];
        memcpy(&rightChunk[0], &rightChunkTmp, 16);

        alignas(16) char flags[16];


        #pragma unroll
        for(int k = 0; k < 16; k++){
            const char leftTransformed = inputTransformOp(leftChunk[k]);
            const char rightTransformed = inputTransformOp(rightChunk[k]);

            if(leftTransformed != rightTransformed){
                flags[k] = 1;
            }else{
                flags[k] = 0;
            }
        }

        float4 flagsAs16; memcpy(&flagsAs16, &flags[0], 16);
        ((float4*)flagsAlignedPtr)[c] = flagsAs16;
    }

    //remaining elements
    for(int64_t i = chunksOf16*16; i < totalInputElements; i += stride){
        const char left = inputAlignedPtr[i-1];
        const char right = inputAlignedPtr[i];
        const char leftTransformed = inputTransformOp(left);
        const char rightTransformed = inputTransformOp(right);

        if(leftTransformed != rightTransformed){
            flagsAlignedPtr[i] = 1;
        }else{
            flagsAlignedPtr[i] = 0;
        }
    }
}

#if !defined(FOOBARAOS)

template<int64_t maxChunkSize, class InputTransformOp>
__global__
void maskingComputeHeadFlagsKernel3(
    //uint32_t* bitarrayFlags,
    char* __restrict__ flags,
    const char* __restrict__ d_inputData,
    int64_t totalInputElements,
    InputTransformOp inputTransformOp
){
    static_assert(maxChunkSize % 16 == 0);

    const int tid = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = blockDim.x * gridDim.x;

    auto block = cooperative_groups::this_thread_block();

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    char* flagsAlignedPtr = cuda::std::assume_aligned<16>(flags);

    constexpr int numStages = 2;

    __shared__ alignas(16) struct Smem{
        char inputElements[numStages][16 + maxChunkSize];
    } smem;

    int computeStage = 0;
    int fetchStage = 0;
    for(int chunk = blockIdx.x, fetchChunk = blockIdx.x; chunk < numChunks; chunk += gridDim.x, computeStage++){
        for (; fetchChunk < numChunks && fetchChunk < (chunk + numStages * gridDim.x); fetchChunk += gridDim.x) {
            const bool isFirstChunk = fetchChunk == 0;
            const bool isLastChunk = fetchChunk == numChunks-1;

            const int64_t chunkBegin = fetchChunk * maxChunkSize;
            const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);
            const int64_t chunkLoadBegin = (!isFirstChunk) ? chunkBegin - 16 : chunkBegin;
            const int64_t chunkLoadEnd = (!isLastChunk) ? cuda::ceil_div(chunkEnd, 16) * 16
                : cuda::ceil_div(chunkEnd, 4) * 4;

            // if(block.thread_rank() == 0){
            //     printf("memcpy_async: stage %d, chunk %d, isFirstChunk %d, isLastChunk %d, chunkBegin %ld, load from offset %ld until %ld\n",
            //     fetchStage, fetchChunk, isFirstChunk, isLastChunk, chunkBegin, chunkLoadBegin, chunkLoadEnd
            //     );
            // }
            char* dstsmemptr = (isFirstChunk) ? &smem.inputElements[fetchStage % numStages][16] : &smem.inputElements[fetchStage % numStages][0];
            if(isLastChunk){
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin)
                );
            }else{
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin)
                );
            }
            fetchStage++;
        }


        const bool isFirstChunk = chunk == 0;
        const bool isLastChunk = chunk == numChunks-1;

        const int64_t chunkBegin = chunk * maxChunkSize;
        const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);
        // const int64_t chunkLoadBegin = (!isFirstChunk) ? chunkBegin - 16 : chunkBegin;
        // const int64_t chunkLoadEnd = (!isLastChunk) ? cuda::ceil_div(chunkEnd, 16) * 16
        //     : cuda::ceil_div(chunkEnd, 4) * 4;

        //wait for async load from gmem to smem of our input chunk
        if(chunk < numChunks - 1 - int(gridDim.x)){
            cooperative_groups::wait_prior<1>(block);
        }else{
            cooperative_groups::wait(block);
        }
        
        if(isFirstChunk){
            //mark first element of input
            if(threadIdx.x == 0){
                flagsAlignedPtr[chunkBegin] = 1;
            }
            const int chunkSize = chunkEnd - chunkBegin;
            for(int i = 1 + threadIdx.x; i < chunkSize; i += blockDim.x){
                const char left = smem.inputElements[computeStage % numStages][16 + i-1];
                const char right = smem.inputElements[computeStage % numStages][16 + i];
                const char leftTransformed = inputTransformOp(left);
                const char rightTransformed = inputTransformOp(right);
        
                if(leftTransformed != rightTransformed){
                    flagsAlignedPtr[chunkBegin + i] = 1;
                }else{
                    flagsAlignedPtr[chunkBegin + i] = 0;
                }
            }
        }else{
            const int chunkSize = chunkEnd - chunkBegin;
            for(int i = threadIdx.x; i < chunkSize; i += blockDim.x){
                const char left = smem.inputElements[computeStage % numStages][16 + i-1];
                const char right = smem.inputElements[computeStage % numStages][16 + i];
                const char leftTransformed = inputTransformOp(left);
                const char rightTransformed = inputTransformOp(right);
        
                if(leftTransformed != rightTransformed){
                    flagsAlignedPtr[chunkBegin + i] = 1;
                }else{
                    flagsAlignedPtr[chunkBegin + i] = 0;
                }
            }
        }
        //chunk is done, can re-use smem
        block.sync();
    }
}

#endif

#if !defined(FOOBARAOS)
template<int64_t maxChunkSize, class InputTransformOp>
__global__
void maskingComputeHeadFlagsKernel4(
    uint32_t* __restrict__ d_bitarrayFlags,
    const char* __restrict__ d_inputData,
    int64_t totalInputElements,
    // InputOffsetIterator d_inputOffsets,
    // int numSegments,
    // const int* __restrict__ searchRangesBegin,
    InputTransformOp inputTransformOp
){
    static_assert(maxChunkSize % 16 == 0);

    const int tid = threadIdx.x + blockIdx.x * blockDim.x;
    const int stride = blockDim.x * gridDim.x;

    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<32>(block);

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    constexpr int numStages = 2;
    constexpr int smemBitArrayNumUint32 = cuda::ceil_div(maxChunkSize, 32);

    __shared__ alignas(16) struct Smem{
        char inputElements[numStages][16 + maxChunkSize];
        uint32_t bitarray[smemBitArrayNumUint32];
    } smem;

    int computeStage = 0;
    int fetchStage = 0;
    for(int chunk = blockIdx.x, fetchChunk = blockIdx.x; chunk < numChunks; chunk += gridDim.x, computeStage++){
        for (; fetchChunk < numChunks && fetchChunk < (chunk + numStages * gridDim.x); fetchChunk += gridDim.x) {
            const bool isFirstChunk = fetchChunk == 0;
            const bool isLastChunk = fetchChunk == numChunks-1;

            const int64_t chunkBegin = fetchChunk * maxChunkSize;
            const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);
            const int64_t chunkLoadBegin = (!isFirstChunk) ? chunkBegin - 16 : chunkBegin;
            const int64_t chunkLoadEnd = (!isLastChunk) ? cuda::ceil_div(chunkEnd, 16) * 16
                : cuda::ceil_div(chunkEnd, 4) * 4;

            // if(block.thread_rank() == 0){
            //     printf("memcpy_async: stage %d, chunk %d, isFirstChunk %d, isLastChunk %d, chunkBegin %ld, load from offset %ld until %ld\n",
            //     fetchStage, fetchChunk, isFirstChunk, isLastChunk, chunkBegin, chunkLoadBegin, chunkLoadEnd
            //     );
            // }
            char* dstsmemptr = (isFirstChunk) ? &smem.inputElements[fetchStage % numStages][16] : &smem.inputElements[fetchStage % numStages][0];
            if(isLastChunk){
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin)
                );
            }else{
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin)
                );
            }
            fetchStage++;
        }


        const bool isFirstChunk = chunk == 0;
        const bool isLastChunk = chunk == numChunks-1;

        const int64_t chunkBegin = chunk * maxChunkSize;
        const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);

        
        // const int64_t chunkLoadBegin = (!isFirstChunk) ? chunkBegin - 16 : chunkBegin;
        // const int64_t chunkLoadEnd = (!isLastChunk) ? cuda::ceil_div(chunkEnd, 16) * 16
        //     : cuda::ceil_div(chunkEnd, 4) * 4;

        //wait for async load from gmem to smem of our input chunk
        if(chunk < numChunks - 1 - int(gridDim.x)){
            cooperative_groups::wait_prior<1>(block);
        }else{
            cooperative_groups::wait(block);
        }        

        const int chunkSize = chunkEnd - chunkBegin;
        const int numBitArrayUint32 = cuda::ceil_div(chunkSize, 32);
        const int loopLimit = numBitArrayUint32 * 32;
        for(int i = threadIdx.x; i < loopLimit; i += blockDim.x){
            char left = 0;
            char right = 0;
            if(i < chunkSize){
                left = smem.inputElements[computeStage % numStages][16 + i-1];
                right = smem.inputElements[computeStage % numStages][16 + i];
            }
            const char leftTransformed = inputTransformOp(left);
            const char rightTransformed = inputTransformOp(right);
            uint32_t flag = leftTransformed != rightTransformed;
            if(i == 0 && isFirstChunk){
                flag = 1;
            }
            const int bitPositionIn32 = i % 32;
            const int indexOf32 = i / 32;
            flag <<= bitPositionIn32;
            const uint32_t flags = compatibility::cooperative_groups::reduce(warp, flag, compatibility::cooperative_groups::bit_or<uint32_t>{});
            smem.bitarray[indexOf32] = flags;
        }
        //wait for smem.bitarray
        block.sync();

        // const int searchRangeBegin = searchRangesBegin[chunk];
        // const int searchRangeEnd = chunk < numChunks ? searchRangesBegin[chunk] : numSegments;

        // for(int i = threadIdx.x; i < chunkSize; i += blockDim.x){
        //     const int searchValue = chunkBegin + i;

        //     auto it = thrust::lower_bound(thrust::seq, d_inputOffsets + searchRangeBegin, d_inputOffsets + searchRangeEnd, searchValue);
        //     if(it != d_inputOffsets + searchRangeEnd){
        //         if(*it == searchValue){
        //             const int bitPositionIn32 = i % 32;
        //             const int indexOf32 = i / 32;
        //             uint32_t* const address = &smem.bitarray[indexOf32];
        //             uint32_t old = *address;
        //             uint32_t assumed = 0;
        //             do {
        //                 assumed = old;
        //                 const uint32_t newVal = (old & ~(1U << bitPositionIn32)) | (1U << bitPositionIn32);
        //                 old = atomicCAS(address, assumed, newVal);
        //             } while (assumed != old);
        //         }
        //     }
        // }

        // //wait for smem.bitarray
        // block.sync();

        //copy bitarray from smem to gmem
        cooperative_groups::memcpy_async(
            block,
            (char*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)),
            cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32),
            (const char*)&smem.bitarray[0],
            cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32)
        );


        //chunk is done, can re-use smem
        block.sync();
    }
}

#endif


#if !defined(FOOBARAOS)

template<bool outputTransformedData, int64_t maxChunkSize, class InputTransformOp>
__global__
void masking10_compute_flags_kernel(
    char* __restrict__ d_transformedOutput, //only accessed if outputTransformedData == true
    uint32_t* __restrict__ d_bitarrayFlags,
    const char* __restrict__ d_inputData,
    int64_t totalInputElements,
    InputTransformOp inputTransformOp
){
    static_assert(maxChunkSize % 16 == 0);


    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<32>(block);

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    constexpr int numStages = 2;
    constexpr int numUint32ForMaxChunkSize = cuda::ceil_div(maxChunkSize, 32);
    constexpr int smemBitArrayNumUint32Padded = cuda::ceil_div(numUint32ForMaxChunkSize, 4) * 4; //pad to 16 bytes

    __shared__ alignas(16) struct Smem{
        char inputElements[numStages][16 + maxChunkSize];
        uint32_t inputBitarray[numStages][smemBitArrayNumUint32Padded];
        uint32_t bitarray[smemBitArrayNumUint32Padded];
        // char outputElements[outputTransformedData ? maxChunkSize : 16];
    } smem;

    int computeStage = 0;
    int fetchStage = 0;
    for(int chunk = blockIdx.x, fetchChunk = blockIdx.x; chunk < numChunks; chunk += gridDim.x, computeStage++){
        for (; fetchChunk < numChunks && fetchChunk < (chunk + numStages * gridDim.x); fetchChunk += gridDim.x) {
            const bool isFirstChunk = fetchChunk == 0;
            const bool isLastChunk = fetchChunk == numChunks-1;

            const int64_t chunkBegin = fetchChunk * maxChunkSize;
            const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);
            const int64_t chunkLoadBegin = (!isFirstChunk) ? chunkBegin - 16 : chunkBegin;
            const int64_t chunkLoadEnd = (!isLastChunk) ? cuda::ceil_div(chunkEnd, 16) * 16
                : cuda::ceil_div(chunkEnd, 4) * 4;

            // if(block.thread_rank() == 0){
            //     printf("memcpy_async: stage %d, chunk %d, isFirstChunk %d, isLastChunk %d, chunkBegin %ld, load from offset %ld until %ld\n",
            //     fetchStage, fetchChunk, isFirstChunk, isLastChunk, chunkBegin, chunkLoadBegin, chunkLoadEnd
            //     );
            // }
            char* dstsmemptr = (isFirstChunk) ? &smem.inputElements[fetchStage % numStages][16] : &smem.inputElements[fetchStage % numStages][0];
            if(isLastChunk){
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin)
                );
            }else{
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin)
                );
            }

            //load the input bitarray of chunk
            const int chunkSize = chunkEnd - chunkBegin;
            const int bitarrayBegin = chunkBegin / 32;

            cooperative_groups::memcpy_async(
                block,
                (char*)(&smem.inputBitarray[fetchStage % numStages][0]),
                cuda::aligned_size_t<4>(cuda::ceil_div(chunkSize, 32) * sizeof(uint32_t)),
                (const char*)(d_bitarrayFlags + bitarrayBegin),
                cuda::aligned_size_t<4>(cuda::ceil_div(chunkSize, 32) * sizeof(uint32_t))
            );
            

            fetchStage++;
        }


        const bool isFirstChunk = chunk == 0;
        // const bool isLastChunk = chunk == numChunks-1;

        const int64_t chunkBegin = chunk * maxChunkSize;
        const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);

        //wait for async load from gmem to smem of our input chunk
        if(chunk < numChunks - 1 - int(gridDim.x)){
            //each cg::memcpy_async commits immediately, so it is 2 outstanding calls per fetch
            cooperative_groups::wait_prior<2>(block);
        }else{
            cooperative_groups::wait(block);
        }        

        const int chunkSize = chunkEnd - chunkBegin;
        const int numBitArrayUint32 = cuda::ceil_div(chunkSize, 32);
        const int loopLimit1 = numBitArrayUint32 * 32;
        for(int i = threadIdx.x; i < loopLimit1; i += blockDim.x){
            char left = 0;
            char right = 0;
            if(i < chunkSize){
                left = smem.inputElements[computeStage % numStages][16 + i-1];
                right = smem.inputElements[computeStage % numStages][16 + i];
            }
            const char leftTransformed = inputTransformOp(left);
            const char rightTransformed = inputTransformOp(right);
            if constexpr(outputTransformedData){
                //smem.inputElements[computeStage % numStages][16 + i] = rightTransformed;
                if(i < chunkSize){
                    d_transformedOutput[chunkBegin + i] = rightTransformed;
                }
            }
            uint32_t flag = leftTransformed != rightTransformed;
            if(i == 0 && isFirstChunk){
                flag = 1;
            }
            const int bitPositionIn32 = 31 - (i % 32);
            const int indexOf32 = i / 32;
            flag <<= bitPositionIn32;
            const uint32_t flags = compatibility::cooperative_groups::reduce(warp, flag, compatibility::cooperative_groups::bit_or<uint32_t>{});
            //smem.bitarray[indexOf32] = flags | smem.inputBitarray[computeStage % numStages][indexOf32];
            (d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32))[indexOf32] = flags | smem.inputBitarray[computeStage % numStages][indexOf32];
        }

        // const int loopLimit2 = cuda::ceil_div(chunkSize, blockDim.x) * blockDim.x;
        // for(int i = threadIdx.x; i < loopLimit2; i += blockDim.x){
        //     char left = 0;
        //     char right = 0;
        //     if(i < chunkSize){
        //         left = smem.inputElements[computeStage % numStages][16 + i-1];
        //         right = smem.inputElements[computeStage % numStages][16 + i];
        //     }
        //     const char leftTransformed = inputTransformOp(left);
        //     const char rightTransformed = inputTransformOp(right);
        //     block.sync();
        //     if(i < loopLimit1){
        //         if constexpr(outputTransformedData){
        //             smem.inputElements[computeStage % numStages][16 + i] = rightTransformed;
        //         }
        //         uint32_t flag = leftTransformed != rightTransformed;
        //         if(i == 0 && isFirstChunk){
        //             flag = 1;
        //         }
        //         const int bitPositionIn32 = 31 - (i % 32);
        //         const int indexOf32 = i / 32;
        //         flag <<= bitPositionIn32;
        //         const uint32_t flags = cooperative_groups::reduce(warp, flag, cooperative_groups::bit_or<uint32_t>{});
        //         smem.bitarray[indexOf32] = flags;

        //         // if constexpr(outputTransformedData){
        //         //     smem.outputElements[i] = right;
        //         // }
        //     }
        // }
        //wait for smem.bitarray
        block.sync();

        // //apply input flags to computed flags
        // for(int i = threadIdx.x; i < numBitArrayUint32; i += blockDim.x){
        //     smem.bitarray[i] |= smem.inputBitarray[computeStage % numStages][i];
        // }

        // //wait for smem.bitarray
        // block.sync();

        //copy bitarray from smem to gmem
        // cooperative_groups::memcpy_async(
        //     block,
        //     (char*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)),
        //     (const char*)&smem.bitarray[0],
        //     cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32)
        // );

        //copy transformed from smem to gmem
        if constexpr(outputTransformedData){

            // for(int i = 4*threadIdx.x; i < cuda::ceil_div(chunkSize, 4); i += 4*blockDim.x){
            //     char4 elems = ((const char4*)(&smem.inputElements[computeStage % numStages][16]))[i/4];
            //     elems.x = inputTransformOp(elems.x);
            //     elems.y = inputTransformOp(elems.y);
            //     elems.z = inputTransformOp(elems.z);
            //     elems.w = inputTransformOp(elems.w);
            //     ((char4*)(d_transformedOutput + chunkBegin))[i/4] = elems;
            // }

            // cooperative_groups::memcpy_async(
            //     block,
            //     (char*)(d_transformedOutput + chunkBegin),
            //     //(const char*)&smem.outputElements[0],
            //     (const char*)&smem.inputElements[computeStage % numStages][16],
            //     cuda::aligned_size_t<4>(sizeof(char) * (cuda::ceil_div(chunkSize, 4) * 4))
            // );
        }


        //chunk is done, can re-use smem
        block.sync();
    }
}

#endif

#if 0

#if !defined(FOOBARAOS)

template<bool outputTransformedData, int64_t maxChunkSize, class InputTransformOp>
__global__
void masking10_compute_flags_kernel_2(
    char* __restrict__ d_transformedOutput, //only accessed if outputTransformedData == true
    uint32_t* __restrict__ d_bitarrayFlags,
    const char* __restrict__ d_inputData,
    int64_t totalInputElements,
    InputTransformOp inputTransformOp
){
    static_assert(maxChunkSize % 16 == 0);


    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<32>(block);
    auto subwarp = cooperative_groups::tiled_partition<8>(block);

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    constexpr int numStages = 2;
    constexpr int numUint32ForMaxChunkSize = cuda::ceil_div(maxChunkSize, 32);
    constexpr int smemBitArrayNumUint32Padded = cuda::ceil_div(numUint32ForMaxChunkSize, 4) * 4; //pad to 16 bytes

    __shared__ alignas(16) struct Smem{
        char inputElements[numStages][16 + maxChunkSize];
        // uint32_t inputBitarray[numStages][smemBitArrayNumUint32Padded];
        // uint32_t bitarray[smemBitArrayNumUint32Padded];
        // char outputElements[outputTransformedData ? maxChunkSize : 16];
    } smem;

    int computeStage = 0;
    int fetchStage = 0;
    for(int chunk = blockIdx.x, fetchChunk = blockIdx.x; chunk < numChunks; chunk += gridDim.x, computeStage++){
        for (; fetchChunk < numChunks && fetchChunk < (chunk + numStages * gridDim.x); fetchChunk += gridDim.x) {
            const bool isFirstChunk = fetchChunk == 0;
            const bool isLastChunk = fetchChunk == numChunks-1;

            const int64_t chunkBegin = fetchChunk * maxChunkSize;
            const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);
            const int64_t chunkLoadBegin = (!isFirstChunk) ? chunkBegin - 16 : chunkBegin;
            const int64_t chunkLoadEnd = (!isLastChunk) ? cuda::ceil_div(chunkEnd, 16) * 16
                : cuda::ceil_div(chunkEnd, 4) * 4;

            // if(block.thread_rank() == 0){
            //     printf("memcpy_async: stage %d, chunk %d, isFirstChunk %d, isLastChunk %d, chunkBegin %ld, load from offset %ld until %ld\n",
            //     fetchStage, fetchChunk, isFirstChunk, isLastChunk, chunkBegin, chunkLoadBegin, chunkLoadEnd
            //     );
            // }
            char* dstsmemptr = (isFirstChunk) ? &smem.inputElements[fetchStage % numStages][16] : &smem.inputElements[fetchStage % numStages][0];
            if(isLastChunk){
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin)
                );
            }else{
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin)
                );
            }

            //load the input bitarray of chunk
            // const int chunkSize = chunkEnd - chunkBegin;
            // const int bitarrayBegin = chunkBegin / 32;

            // cooperative_groups::memcpy_async(
            //     block,
            //     (char*)(&smem.inputBitarray[fetchStage % numStages][0]),
            //     (const char*)(d_bitarrayFlags + bitarrayBegin),
            //     cuda::aligned_size_t<4>(cuda::ceil_div(chunkSize, 32) * sizeof(uint32_t))
            // );
            

            fetchStage++;
        }


        const bool isFirstChunk = chunk == 0;
        const bool isLastChunk = chunk == numChunks-1;

        const int64_t chunkBegin = chunk * maxChunkSize;
        const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);

        //wait for async load from gmem to smem of our input chunk
        if(chunk < numChunks - 1 - int(gridDim.x)){

            cooperative_groups::wait_prior<1>(block);
        }else{
            cooperative_groups::wait(block);
        }
        
        //

        const int chunkSize = chunkEnd - chunkBegin;
        const int numBitArrayUint32 = cuda::ceil_div(chunkSize, 32);

        //each subwarp computes 8 results. each results requires 8 collaborating threads
        for(int outputUintChunk = subwarp.meta_group_rank(); outputUintChunk < cuda::ceil_div(numBitArrayUint32, subwarp.size()); outputUintChunk += subwarp.meta_group_size()){
            uint32_t myResult;

            const uint32_t* smeminputAsUint = (const uint32_t*)&smem.inputElements[computeStage % numStages][16];
            #pragma unroll
            for(int k = 0; k < 8; k++){
                if(outputUintChunk + k < numBitArrayUint32){
                    uint32_t loadedLeft = smeminputAsUint[8 * (outputUintChunk + k) + subwarp.thread_rank() - 1];
                    uint32_t loadedRight = smeminputAsUint[8 * (outputUintChunk + k) + subwarp.thread_rank()];
                    loadedLeft = __funnelshift_r(loadedRight, loadedLeft, 8);
                    //TODO transform
                    #ifdef HAS_BLACKWELL_INT8_PTX
                    //if difference is not 0, corresponding input bytes are not equal
                    uint32_t differences = ptx_sub_s8x4(loadedLeft,loadedRight);
                    uint32_t packedFlags = ((((differences >> 24) & 0xFF) > 0 ? 1u : 0u) << 3) 
                        | ((((differences >> 16) & 0xFF) > 0 ? 1u : 0u) << 2) 
                        | ((((differences >> 8) & 0xFF) > 0 ? 1u : 0u) << 1) 
                        | ((((differences) & 0xFF) > 0 ? 1u : 0u)) ;
                    #else
                    uint32_t compresult = __vcmpne4(loadedLeft, loadedRight);
                    uint32_t packedFlags = ((compresult >> 24) & 1) | ((compresult >> 16) & 1) | ((compresult >> 8) & 1) | (compresult & 1);
                    #endif
                    packedFlags <<= subwarp.thread_rank();
                    uint32_t flags = cooperative_groups::reduce(warp, packedFlags, cooperative_groups::bit_or<uint32_t>{});
                    if(k == subwarp.thread_rank()){
                        myResult = flags;
                    }
                }
            }

            if(outputUintChunk + subwarp.thread_rank() < numBitArrayUint32){
                (d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32))[outputUintChunk + subwarp.thread_rank()] = myResult;
                // smem.bitarray[outputUintChunk + subwarp.thread_rank()] = myResult;
            }
        }

        // block.sync();
        // //copy bitarray from smem to gmem
        // cooperative_groups::memcpy_async(
        //     block,
        //     (char*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)),
        //     (const char*)&smem.bitarray[0],
        //     cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32)
        // );

        //chunk is done, can re-use smem
        block.sync();
    }
}

#endif
#endif 



#if !defined(FOOBARAOS)

template<int blocksize, bool outputTransformedData, int64_t maxChunkSize, class InputTransformOp>
__global__
void masking10_compute_flags_kernel_3(
    char* __restrict__ d_transformedOutput, //only accessed if outputTransformedData == true
    uint32_t* __restrict__ d_bitarrayFlags,
    const char* __restrict__ d_inputData,
    int64_t totalInputElements,
    // InputOffsetIterator d_inputOffsets,
    // int numSegments,
    // const int* __restrict__ searchRangesBegin,
    InputTransformOp inputTransformOp
){
    static_assert(maxChunkSize % 16 == 0);
    __builtin_assume(blockDim.x == blocksize);
    __builtin_assume(blockDim.y == 1);
    __builtin_assume(blockDim.z == 1);

    // const int tid = threadIdx.x + blockIdx.x * blockDim.x;
    // const int stride = blockDim.x * gridDim.x;

    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<32>(block);
    auto groupOf2 = cooperative_groups::tiled_partition<2>(block);

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    constexpr int numStages = 2;
    constexpr int smemBitArrayNumUint32 = cuda::ceil_div(maxChunkSize, 32);

    __shared__ alignas(16) struct Smem{
        char inputElements[numStages][16 + maxChunkSize];
        uint32_t bitarray[smemBitArrayNumUint32];
    } smem;

    int computeStage = 0;
    int fetchStage = 0;
    for(int chunk = blockIdx.x, fetchChunk = blockIdx.x; chunk < numChunks; chunk += gridDim.x, computeStage++){
        for (; fetchChunk < numChunks && fetchChunk < (chunk + numStages * gridDim.x); fetchChunk += gridDim.x) {
            const bool isFirstChunk = fetchChunk == 0;
            const bool isLastChunk = fetchChunk == numChunks-1;

            const int64_t chunkBegin = fetchChunk * maxChunkSize;
            const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);
            const int64_t chunkLoadBegin = (!isFirstChunk) ? chunkBegin - 16 : chunkBegin;
            const int64_t chunkLoadEnd = (!isLastChunk) ? cuda::ceil_div(chunkEnd, 16) * 16
                : cuda::ceil_div(chunkEnd, 4) * 4;

            // if(block.thread_rank() == 0){
            //     printf("memcpy_async: stage %d, chunk %d, isFirstChunk %d, isLastChunk %d, chunkBegin %ld, load from offset %ld until %ld\n",
            //     fetchStage, fetchChunk, isFirstChunk, isLastChunk, chunkBegin, chunkLoadBegin, chunkLoadEnd
            //     );
            // }
            char* dstsmemptr = (isFirstChunk) ? &smem.inputElements[fetchStage % numStages][16] : &smem.inputElements[fetchStage % numStages][0];
            if(isLastChunk){
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin)
                );
            }else{
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin)
                );
            }
            fetchStage++;
        }


        // const bool isFirstChunk = chunk == 0;
        // const bool isLastChunk = chunk == numChunks-1;

        const int64_t chunkBegin = chunk * maxChunkSize;
        const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);
        const int chunkSize = chunkEnd - chunkBegin;
        
        // const int64_t chunkLoadBegin = (!isFirstChunk) ? chunkBegin - 16 : chunkBegin;
        // const int64_t chunkLoadEnd = (!isLastChunk) ? cuda::ceil_div(chunkEnd, 16) * 16
        //     : cuda::ceil_div(chunkEnd, 4) * 4;

        //wait for async load from gmem to smem of our input chunk
        if(chunk < numChunks - 1 - int(gridDim.x)){
            cooperative_groups::wait_prior<1>(block);
        }else{
            cooperative_groups::wait(block);
        }

        if constexpr(outputTransformedData || !std::is_same_v<InputTransformOp, cuda::std::identity>){
            using CopyType = uint4;
            const int fullCopyChunks = chunkSize / sizeof(CopyType);
            for(int i = threadIdx.x; i < fullCopyChunks; i += blockDim.x){
                CopyType* smemptr = ((CopyType*)&smem.inputElements[computeStage % numStages][16]) + i;
                CopyType tmp = *smemptr;
                alignas(sizeof(CopyType)) char elements[sizeof(CopyType)];
                memcpy(&elements[0], &tmp, sizeof(CopyType));                
                for(int k = 0; k < sizeof(CopyType); k++){
                    elements[k] = inputTransformOp(elements[k]);
                }
                memcpy(&tmp, &elements[0],  sizeof(CopyType));
                *smemptr = tmp;

                if constexpr(outputTransformedData){
                    ((CopyType*)(d_transformedOutput + chunkBegin))[i] = tmp;
                }
            }
            //handle remainder element-wise
            for(int i = sizeof(CopyType)*fullCopyChunks + threadIdx.x; i < chunkSize; i += blockDim.x){
                char element = smem.inputElements[computeStage % numStages][16 + i];
                element = inputTransformOp(element);
                smem.inputElements[computeStage % numStages][16 + i] = element;
                
                if constexpr(outputTransformedData){
                    (d_transformedOutput + chunkBegin)[i] = element;
                }
            }

            block.sync();
        }

        //variant e
        #if 1
            const int numBitArrayUint32 = cuda::ceil_div(chunkSize, 32);
            //use 2 threads for each uint, each thread computes 16 bits of flags
            for(int indexOf32 = groupOf2.meta_group_rank(); indexOf32 < numBitArrayUint32; indexOf32 += groupOf2.meta_group_size()){

                uint4 tmp1 = *((uint4*)&smem.inputElements[computeStage % numStages][16 + (32 * indexOf32 + 16 * groupOf2.thread_rank())]);

                char left = smem.inputElements[computeStage % numStages][16 + (32 * indexOf32 + 16 * groupOf2.thread_rank() - 1)];

                //shift 16 bytes up by one position.
                alignas(16) uint32_t uuu[4];
                memcpy(&uuu[0], &tmp1, 16);

                uuu[3] = (uuu[2] >> 24) | (uuu[3] << 8);
                uuu[2] = (uuu[1] >> 24) | (uuu[2] << 8);
                uuu[1] = (uuu[0] >> 24) | (uuu[1] << 8);
                uuu[0] =   left | (uuu[0] << 8);

                #ifdef HAS_BLACKWELL_INT8_PTX

                    alignas(16) uint32_t center[4];
                    memcpy(&center[0], &tmp1, 16);

                    uint32_t flags = 0;
                    for(int i = 0; i < 4; i++){
                        //if difference is not 0, corresponding input bytes are not equal
                        const uint32_t differences = ptx_sub_s8x4(uuu[i],center[i]);
                        //masks given in little-endian
                        flags = (flags << 1) | ((differences & 0x000000FF) > 0);
                        flags = (flags << 1) | ((differences & 0x0000FF00) > 0);
                        flags = (flags << 1) | ((differences & 0x00FF0000) > 0);
                        flags = (flags << 1) | ((differences & 0xFF000000) > 0);
                    }
               

                #else

                    alignas(16) __nv_fp8x2_e4m3 centerfp8x2[8];
                    memcpy(&centerfp8x2[0], &tmp1, 16);
                    alignas(16) __nv_fp8x2_e4m3 leftfp8x2[8];
                    memcpy(&leftfp8x2[0], &uuu[0], 16);
                    
                    uint32_t flags = 0;
                    for(int i = 0; i < 8; i++){
                        half2 l = static_cast<half2>(leftfp8x2[i]);
                        half2 c = static_cast<half2>(centerfp8x2[i]);

                        flags = (flags << 1) | (l.x != c.x);
                        flags = (flags << 1) | (l.y != c.y);
                    }
                #endif

                if(groupOf2.thread_rank() == 0){
                    flags <<= 16;
                }

                //combine the 16 flags of two collaborating threads
                
                //set oob flags to 1
                const int unoccupied = 32 - cuda::std::min(32, chunkSize - indexOf32 * 32);
                flags |= (uint32_t(1) << unoccupied)-1;

                if(indexOf32 < numBitArrayUint32){
                    alignas(4) uint16_t u16[2];
                    memcpy(&u16[0], &flags, 4);

                    //((uint16_t*)&smem.bitarray[indexOf32])[1-groupOf2.thread_rank()] = groupOf2.thread_rank() == 0 ? u16[1] : u16[0];
                    ((uint16_t*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32) + indexOf32))[1-groupOf2.thread_rank()] = groupOf2.thread_rank() == 0 ? u16[1] : u16[0];
                }
            }
            //wait for smem.bitarray
            // block.sync();

            // //copy bitarray from smem to gmem
            // cooperative_groups::memcpy_async(
            //     block,
            //     (char*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)),
            //     (const char*)&smem.bitarray[0],
            //     cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32)
            // );
        #endif

        //variant d
        #if 0
            const int numBitArrayUint32 = cuda::ceil_div(chunkSize, 32);
            //use 2 threads for each uint, each thread computes 16 bits of flags
            for(int indexOf32 = groupOf2.meta_group_rank(); indexOf32 < numBitArrayUint32; indexOf32 += groupOf2.meta_group_size()){

                uint4 tmp1 = *((uint4*)&smem.inputElements[computeStage % numStages][16 + (32 * indexOf32 + 16 * groupOf2.thread_rank())]);

                char left = smem.inputElements[computeStage % numStages][16 + (32 * indexOf32 + 16 * groupOf2.thread_rank() - 1)];

                // if(chunk == 18){
                //     if(indexOf32 == 12){
                //             alignas(16) char pp[16];
                //             memcpy(&pp[0], &tmp1, 16);

                //             if(groupOf2.thread_rank() == 0){
                //                 printf("input chars\n");
                //                 for(int i = 0; i < 16; i++){
                //                     printf("%d ", int(pp[i]));
                //                 }
                //                 // printf("\n");
                //             }
                //             groupOf2.sync();

                //             if(groupOf2.thread_rank() == 1){
                //                 for(int i = 0; i < 16; i++){
                //                     printf("%d ", int(pp[i]));
                //                 }
                //                 printf("\n");
                //             }
                //             groupOf2.sync();
                //     }
                // }

                //shift 16 bytes up by one position.
                alignas(16) uint32_t uuu[4];
                memcpy(&uuu[0], &tmp1, 16);

                uuu[3] = (uuu[2] >> 24) | (uuu[3] << 8);
                uuu[2] = (uuu[1] >> 24) | (uuu[2] << 8);
                uuu[1] = (uuu[0] >> 24) | (uuu[1] << 8);
                uuu[0] =   left | (uuu[0] << 8);

                #ifdef HAS_BLACKWELL_INT8_PTX

                alignas(16) uint32_t center[4];
                memcpy(&center[0], &tmp1, 16);

                uint32_t flags = 0;
                for(int i = 0; i < 4; i++){
                    //if difference is not 0, corresponding input bytes are not equal
                    const uint32_t differences = ptx_sub_s8x4(uuu[i],center[i]);
                    //masks given in little-endian
                    flags = (flags << 1) | ((differences & 0x000000FF) > 0);
                    flags = (flags << 1) | ((differences & 0x0000FF00) > 0);
                    flags = (flags << 1) | ((differences & 0x00FF0000) > 0);
                    flags = (flags << 1) | ((differences & 0xFF000000) > 0);
                }
               

                #else

                alignas(16) __nv_fp8x2_e4m3 centerfp8x2[8];
                memcpy(&centerfp8x2[0], &tmp1, 16);
                alignas(16) __nv_fp8x2_e4m3 leftfp8x2[8];
                memcpy(&leftfp8x2[0], &uuu[0], 16);
                
                uint32_t flags = 0;
                for(int i = 0; i < 8; i++){
                    half2 l = static_cast<half2>(leftfp8x2[i]);
                    half2 c = static_cast<half2>(centerfp8x2[i]);

                    flags = (flags << 1) | (l.x != c.x);
                    flags = (flags << 1) | (l.y != c.y);
                }
                #endif
                if(groupOf2.thread_rank() == 0){
                    flags <<= 16;
                }

                //combine the 16 flags of two collaborating threads
                flags = cooperative_groups::reduce(groupOf2, flags, cooperative_groups::bit_or<uint32_t>{});

                //set oob flags to 1
                const int unoccupied = 32 - cuda::std::min(32, chunkSize - indexOf32 * 32);
                flags |= (uint32_t(1) << unoccupied)-1;
                if(indexOf32 < numBitArrayUint32){
                    if(groupOf2.thread_rank() == 0){
                        smem.bitarray[indexOf32] = flags;
                    }
                }
            }
            //wait for smem.bitarray
            block.sync();

            //copy bitarray from smem to gmem
            cooperative_groups::memcpy_async(
                block,
                (char*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)),
                (const char*)&smem.bitarray[0],
                cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32)
            );
        #endif

        // variant c
        #if 0 
            const int numBitArrayUint32 = cuda::ceil_div(chunkSize, 32);
            const int smemLoadIterations = cuda::ceil_div(numBitArrayUint32, 32);
            for(int smemLoadIteration = warp.meta_group_rank(); smemLoadIteration < smemLoadIterations; smemLoadIteration += warp.meta_group_size()){
                const int indexOf32 = smemLoadIteration * warp.size() + warp.thread_rank();


                uint4 tmp1_a_foobar = make_uint4(0,0,0,0);
                const int loadOffset1 = smemLoadIteration * 32 * 32 + 16 * warp.thread_rank();
                if(loadOffset1 < maxChunkSize){
                    tmp1_a_foobar = *((uint4*)&smem.inputElements[computeStage % numStages][16 + loadOffset1]);
                }

                uint4 tmp2_a_foobar = make_uint4(0,0,0,0);
                const int loadOffset2 = smemLoadIteration * 32 * 32 + 16 * warp.size() + 16 * warp.thread_rank();
                if(loadOffset2 < maxChunkSize){
                    tmp2_a_foobar = *((uint4*)&smem.inputElements[computeStage % numStages][16 + loadOffset2]);
                }

                const uint4 tmp1_foobar_1 = warp.shfl(tmp1_a_foobar, (2 * warp.thread_rank()) % 32);
                const uint4 tmp1_foobar_2 = warp.shfl(tmp2_a_foobar, (2 * warp.thread_rank()) % 32);
                const uint4 tmp2_foobar_1 = warp.shfl(tmp1_a_foobar, (2 * warp.thread_rank() + 1) % 32);
                const uint4 tmp2_foobar_2 = warp.shfl(tmp2_a_foobar, (2 * warp.thread_rank() + 1) % 32);
                const uint4 tmp1_foobar = warp.thread_rank() < 16 ? tmp1_foobar_1 : tmp1_foobar_2;
                const uint4 tmp2_foobar = warp.thread_rank() < 16 ? tmp2_foobar_1 : tmp2_foobar_2;
                
                const uint4 tmp1 = tmp1_foobar;
                const uint4 tmp2 = tmp2_foobar;
                char left = smem.inputElements[computeStage % numStages][16 + (32 * indexOf32 - 1)];
                alignas(32) char myChars[32];
                memcpy(&myChars[0], &tmp1, 16);
                memcpy(&myChars[16], &tmp2, 16);
                uint32_t flags = 0;
                for(int i = 0; i < 32; i++){
                    uint32_t comp = left != myChars[i];
                    const int bitPositionIn32 = 31 - (i % 32);
                    flags |= (comp << bitPositionIn32);
                    left = myChars[i];
                }
                //set oob flags to 1
                const int unoccupied = 32 - cuda::std::min(32, chunkSize - indexOf32 * 32);
                flags |= (uint32_t(1) << unoccupied)-1;
                if(indexOf32 < numBitArrayUint32){
                    smem.bitarray[indexOf32] = flags;
                }
            }
            //wait for smem.bitarray
            block.sync();

            //copy bitarray from smem to gmem
            cooperative_groups::memcpy_async(
                block,
                (char*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)),
                (const char*)&smem.bitarray[0],
                cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32)
            );
        #endif

        //variant b
        #if 0
            const int numBitArrayUint32 = cuda::ceil_div(chunkSize, 32);
            const int loopLimit = cuda::ceil_div(numBitArrayUint32, 32) * 32;
            for(int indexOf32 = threadIdx.x; indexOf32 < loopLimit; indexOf32 += blockDim.x){
                //two adjacent threads collaborate on loading 32 byte from shared memory,
                const int indexOf32LoadIndex0 = 2*(indexOf32 / 2);
                const int lane = warp.thread_rank() % 2;
                uint4 tmp1_a_foobar = make_uint4(0,0,0,0);
                //lower thread loads first 16 byte of lower thread, upper thread loads second 16 byte of lower thread
                if((32 * indexOf32LoadIndex0) + 16 * lane < maxChunkSize){
                    tmp1_a_foobar = *((uint4*)&smem.inputElements[computeStage % numStages][16 + (32 * indexOf32LoadIndex0) + 16*lane]);
                }
                const int indexOf32LoadIndex1 = indexOf32LoadIndex0 + 1;
                uint4 tmp2_a_foobar = make_uint4(0,0,0,0);
                if(32 * indexOf32LoadIndex1 + 16 * lane < maxChunkSize){
                    tmp2_a_foobar = *((uint4*)&smem.inputElements[computeStage % numStages][16 + (32 * indexOf32LoadIndex1) + 16*lane]);
                }
                const uint4 shfl_u_foobar = warp.shfl_up(tmp2_a_foobar, 1);
                const uint4 shfl_d_foobar = warp.shfl_down(tmp1_a_foobar, 1);
                const uint4 tmp1_foobar = lane == 0 ? tmp1_a_foobar : shfl_u_foobar;
                const uint4 tmp2_foobar = lane == 0 ? shfl_d_foobar : tmp2_a_foobar;


                const uint4 tmp1 = tmp1_foobar;
                const uint4 tmp2 = tmp2_foobar;

                char left = smem.inputElements[computeStage % numStages][16 + (32 * indexOf32 - 1)];
                alignas(32) char myChars[32];
                memcpy(&myChars[0], &tmp1, 16);
                memcpy(&myChars[16], &tmp2, 16);
                uint32_t flags = 0;
                for(int i = 0; i < 32; i++){
                    uint32_t comp = left != myChars[i];
                    const int bitPositionIn32 = 31 - (i % 32);
                    flags |= (comp << bitPositionIn32);
                    left = myChars[i];
                }
                //set oob flags to 1
                const int unoccupied = 32 - cuda::std::min(32, chunkSize - indexOf32 * 32);
                flags |= (uint32_t(1) << unoccupied)-1;
                if(indexOf32 < numBitArrayUint32){
                    smem.bitarray[indexOf32] = flags;
                }
            }
            //wait for smem.bitarray
            block.sync();

            //copy bitarray from smem to gmem
            cooperative_groups::memcpy_async(
                block,
                (char*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)),
                (const char*)&smem.bitarray[0],
                cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32)
            );
        #endif

        //variant a
        #if 0
            //each thread computes 32-bit flag
            const int numBitArrayUint32 = cuda::ceil_div(chunkSize, 32);
            for(int indexOf32 = threadIdx.x; indexOf32 < numBitArrayUint32; indexOf32 += blockDim.x){
                uint4 tmp1 = *((uint4*)&smem.inputElements[computeStage % numStages][16 + (32 * indexOf32)]);
                uint4 tmp2 = make_uint4(0,0,0,0);
                if(32*indexOf32 + 16 < maxChunkSize){
                    tmp2 = *((uint4*)&smem.inputElements[computeStage % numStages][16 + (32 * indexOf32 + 16)]);
                }
                
                char left = smem.inputElements[computeStage % numStages][16 + (32 * indexOf32 - 1)];
                alignas(32) char myChars[32];
                memcpy(&myChars[0], &tmp1, 16);
                memcpy(&myChars[16], &tmp2, 16);
                uint32_t flags = 0;
                for(int i = 0; i < 32; i++){
                    uint32_t comp = left != myChars[i];
                    const int bitPositionIn32 = 31 - (i % 32);
                    flags |= (comp << bitPositionIn32);
                    left = myChars[i];
                }
                //set oob flags to 1
                const int unoccupied = 32 - cuda::std::min(32, chunkSize - indexOf32 * 32);
                flags |= (uint32_t(1) << unoccupied)-1;
                if(indexOf32 < numBitArrayUint32){
                    smem.bitarray[indexOf32] = flags;
                }
            }
            //wait for smem.bitarray
            block.sync();

            //copy bitarray from smem to gmem
            cooperative_groups::memcpy_async(
                block,
                (char*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)),
                (const char*)&smem.bitarray[0],
                cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32)
            );
        #endif

        //original impl
        #if 0
            const int numBitArrayUint32 = cuda::ceil_div(chunkSize, 32);
            const int loopLimit = numBitArrayUint32 * 32;
            for(int i = threadIdx.x; i < loopLimit; i += blockDim.x){
                char left = 0;
                char right = 0;
                if(i < chunkSize){
                    left = smem.inputElements[computeStage % numStages][16 + i-1];
                    right = smem.inputElements[computeStage % numStages][16 + i];
                }
                const char leftTransformed = left; //inputTransformOp(left);
                const char rightTransformed = right; //inputTransformOp(right);
                uint32_t flag = leftTransformed != rightTransformed;
                const int bitPositionIn32 = 31 - (i % 32);
                const int indexOf32 = i / 32;
                flag <<= bitPositionIn32;
                const uint32_t flags = cooperative_groups::reduce(warp, flag, cooperative_groups::bit_or<uint32_t>{});
                smem.bitarray[indexOf32] = flags;
            }
            //wait for smem.bitarray
            block.sync();

            //copy bitarray from smem to gmem
            cooperative_groups::memcpy_async(
                block,
                (char*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)),
                (const char*)&smem.bitarray[0],
                cuda::aligned_size_t<4>(sizeof(uint32_t) * numBitArrayUint32)
            );
        #endif


        //chunk is done, can re-use smem
        block.sync();
    }
}

#endif


template<int blocksize, bool outputTransformedData, int64_t maxChunkSize, class InputTransformOp>
__global__
void masking10_compute_flags_kernel_4(
    char* __restrict__ d_transformedOutput, //only accessed if outputTransformedData == true
    uint32_t* __restrict__ d_bitarrayFlags,
    const char* __restrict__ d_inputData,
    int64_t totalInputElements,
    InputTransformOp inputTransformOp
){
    static_assert(maxChunkSize % 16 == 0);
    __builtin_assume(blockDim.x == blocksize);
    __builtin_assume(blockDim.y == 1);
    __builtin_assume(blockDim.z == 1);

    // const int tid = threadIdx.x + blockIdx.x * blockDim.x;
    // const int stride = blockDim.x * gridDim.x;

    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<32>(block);
    auto groupOf2 = cooperative_groups::tiled_partition<2>(block);

    const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;
    const int numWarpsInGrid = (blockDim.x * gridDim.x) / 32;

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    //each warp processes maxChunkSize input elements

    for(int chunk = warpIdInGrid; chunk < numChunks; chunk += numWarpsInGrid){
        const bool isFirstChunk = chunk == 0;
        // const bool isLastChunk = chunk == numChunks-1;
        const int64_t chunkBegin = chunk * maxChunkSize;
        const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);
        const int chunkSize = chunkEnd - chunkBegin;

        //each thread in warp processes 16 input elements
        const int num16chunks = cuda::ceil_div(chunkSize, 16);
        const int loopLimit1 = cuda::ceil_div(num16chunks, warp.size()) * warp.size();

        alignas(16) char myElementsAsChars[16];
        char leftNeighborElement = 0;
        uint4 myElements_ui4_next = make_uint4(0,0,0,0);
        //load my first 16 elements
        if(warp.thread_rank() < num16chunks){
            const int c16 = warp.thread_rank();
            myElements_ui4_next = *reinterpret_cast<const uint4*>(d_inputData + chunkBegin + 16 * c16);
        }
        for(int c16 = warp.thread_rank(); c16 < loopLimit1; c16 += warp.size()){
            const char shuffledLeftNeighborElement_b = warp.shfl(myElementsAsChars[15], warp.size() - 1);

            uint4 myElements_ui4 = myElements_ui4_next;

            //fire up load for next iteration
            if(c16 + warp.size() < num16chunks){
                const int c16_next = c16 + warp.size();
                myElements_ui4_next = *reinterpret_cast<const uint4*>(d_inputData + chunkBegin + 16 * c16_next);
            }
            memcpy(&myElementsAsChars[0], &myElements_ui4, 16);

            // if(chunk == 451){
            //     if(c16 == 96){
            //         alignas(16) char pp[16];
            //         memcpy(&pp[0], &myElements_ui4, 16);
            //             printf("input chars\n");
            //             for(int i = 0; i < 16; i++){
            //                 printf("%d ", int(pp[i]));
            //             }
            //             printf("\n");
            //     }
            //     warp.sync();
            // }

            //for first iteration, load left neighbor from global memory
            if(c16 == 0){
                if(!isFirstChunk){
                    leftNeighborElement = *(d_inputData + chunkBegin + 16 * c16 - 1);
                }
            }

            //other iterations fetch left element from neighbor thread's elements
            const char shuffledLeftNeighborElement_a = warp.shfl_up(myElementsAsChars[15], 1);
            if(warp.thread_rank() > 0){
                leftNeighborElement = shuffledLeftNeighborElement_a;
            }else{
                if(c16 > 0){
                    leftNeighborElement = shuffledLeftNeighborElement_b;
                }
            }
            // if(chunk == 451){
            //     if(c16 == 96){
            //         printf("leftNeighborElement %d\n", int(leftNeighborElement));
            //     }
            // }

            //transform
            if constexpr(outputTransformedData || !std::is_same_v<InputTransformOp, cuda::std::identity>){
                for(int k = 0; k < 16; k++){
                    myElementsAsChars[k] = inputTransformOp(myElementsAsChars[k]);
                }
                memcpy(&myElements_ui4, &myElementsAsChars[0], 16);
                if constexpr(outputTransformedData){
                    if(c16 < num16chunks){
                        *reinterpret_cast<uint4*>(d_transformedOutput + chunkBegin + 16 * c16) = myElements_ui4;
                    }
                }
            }

            #if 1

            uint32_t flags = leftNeighborElement != myElementsAsChars[0];
            for(int i = 1; i < 16; i++){
                flags = (flags << 1) | (myElementsAsChars[i-1] != myElementsAsChars[i]);
            }

            #else 
                //some more elaborate approaches, but not really faster. memory-bound after all

                //#ifdef HAS_BLACKWELL_INT8_PTX
                #if 0
                    //shift elements up by one position to get a left neighbor element for each of my elements.
                    alignas(16) uint32_t left16[4];
                    memcpy(&left16[0], &myElements_ui4, 16);

                    left16[3] = (left16[2] >> 24) | (left16[3] << 8);
                    left16[2] = (left16[1] >> 24) | (left16[2] << 8);
                    left16[1] = (left16[0] >> 24) | (left16[1] << 8);
                    left16[0] = leftNeighborElement | (left16[0] << 8);
                    
                    alignas(16) uint32_t center[4];
                    memcpy(&center[0], &myElements_ui4, 16);

                    uint32_t flags = 0;
                    for(int i = 0; i < 4; i++){
                        //if difference is not 0, corresponding input bytes are not equal
                        const uint32_t differences = ptx_sub_s8x4(left16[i],center[i]);
                        //masks given in little-endian
                        flags = (flags << 1) | ((differences & 0x000000FF) > 0);
                        flags = (flags << 1) | ((differences & 0x0000FF00) > 0);
                        flags = (flags << 1) | ((differences & 0x00FF0000) > 0);
                        flags = (flags << 1) | ((differences & 0xFF000000) > 0);
                    }
                

                #else

                    //shift elements up by one position to get a left neighbor element for each of my elements.
                    alignas(16) uint32_t left16[4];
                    memcpy(&left16[0], &myElements_ui4, 16);
        
                    left16[3] = (left16[2] >> 24) | (left16[3] << 8);
                    left16[2] = (left16[1] >> 24) | (left16[2] << 8);
                    left16[1] = (left16[0] >> 24) | (left16[1] << 8);
                    left16[0] = leftNeighborElement | (left16[0] << 8);

                    //using half2 enables HSETP2 which can set comparison flag for 2 halfs

                    alignas(16) __nv_fp8x2_e4m3 centerfp8x2[8];
                    memcpy(&centerfp8x2[0], &myElements_ui4, 16);
                    alignas(16) __nv_fp8x2_e4m3 leftfp8x2[8];
                    memcpy(&leftfp8x2[0], &left16[0], 16);
                    
                    uint32_t flags = 0;
                    for(int i = 0; i < 8; i++){
                        half2 l = static_cast<half2>(leftfp8x2[i]);
                        half2 c = static_cast<half2>(centerfp8x2[i]);

                        flags = (flags << 1) | (l.x != c.x);
                        flags = (flags << 1) | (l.y != c.y);
                    }
                #endif
            #endif


            //set oob flags to 1
            const int unoccupied = 16 - cuda::std::min(16, chunkSize - c16 * 16);
            flags |= (uint32_t(1) << unoccupied)-1;

            if(c16 < num16chunks){
                const int outputindex = groupOf2.thread_rank() == 0 ? c16+1 : c16-1;
                ((uint16_t*)(d_bitarrayFlags + chunk * cuda::ceil_div(maxChunkSize, 32)))[outputindex] = flags;
            }
        }

    }
}

#if !defined(FOOBARAOS)

template<int64_t maxChunkSize>
__global__
void masking10_UseFlagsKernel(
    char* __restrict__ d_outputData,
    const uint32_t* __restrict__ d_bitarrayFlags,
    const char* __restrict__ d_inputData,
    int64_t totalInputElements,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter
){
    static_assert(maxChunkSize % 16 == 0);
    assert(maskingThreshold <= 32);


    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<32>(block);

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    constexpr int numStages = 2;
    constexpr int numUint32ForMaxChunkSize = 1 + cuda::ceil_div(maxChunkSize, 32) + 1;
    constexpr int smemBitArrayNumUint32Padded = cuda::ceil_div(numUint32ForMaxChunkSize, 4) * 4; //pad to 16 bytes

    __shared__ alignas(16) struct Smem{
        char inputElements[numStages][maxChunkSize];
        uint32_t inputBitarray[numStages][smemBitArrayNumUint32Padded];
    } smem;

    const int bitarrayTotalElements = cuda::ceil_div(totalInputElements, 32);
    int computeStage = 0;
    int fetchStage = 0;
    for(int chunk = blockIdx.x, fetchChunk = blockIdx.x; chunk < numChunks; chunk += gridDim.x, computeStage++){
        for (; fetchChunk < numChunks && fetchChunk < (chunk + numStages * gridDim.x); fetchChunk += gridDim.x) {
            const bool isFirstChunk = fetchChunk == 0;
            const bool isLastChunk = fetchChunk == numChunks-1;

            const int64_t chunkBegin = fetchChunk * maxChunkSize;
            const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);
            const int64_t chunkLoadBegin = chunkBegin;
            const int64_t chunkLoadEnd = (!isLastChunk) ? cuda::ceil_div(chunkEnd, 16) * 16
                : cuda::ceil_div(chunkEnd, 4) * 4;

            // if(block.thread_rank() == 0){
            //     printf("memcpy_async: stage %d, chunk %d, isFirstChunk %d, isLastChunk %d, chunkBegin %ld, load from offset %ld until %ld\n",
            //     fetchStage, fetchChunk, isFirstChunk, isLastChunk, chunkBegin, chunkLoadBegin, chunkLoadEnd
            //     );
            // }
            char* dstsmemptr = &smem.inputElements[fetchStage % numStages][0];
            if(isLastChunk){
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<4>(chunkLoadEnd - chunkLoadBegin)
                );
            }else{
                cooperative_groups::memcpy_async(
                    block,
                    dstsmemptr,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin),
                    d_inputData + chunkLoadBegin,
                    cuda::aligned_size_t<16>(chunkLoadEnd - chunkLoadBegin)
                );
            }

            //load the input bitarray of chunk
            const int chunkSize = chunkEnd - chunkBegin;
            const int bitarrayBegin = chunkBegin / 32;
            const int bitarrayEnd = bitarrayBegin + cuda::ceil_div(chunkSize, 32);
            const int bitarrayLoadBegin = isFirstChunk ? bitarrayBegin : bitarrayBegin - 1;
            const int bitarrayLoadEnd = cuda::std::min(bitarrayTotalElements, bitarrayEnd + 1);
            const int bitarrayLoadElements = bitarrayLoadEnd - bitarrayLoadBegin;
            if(isFirstChunk){
                cooperative_groups::memcpy_async(
                    block,
                    (char*)(&smem.inputBitarray[fetchStage % numStages][1]),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * bitarrayLoadElements),
                    (const char*)(d_bitarrayFlags + bitarrayLoadBegin),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * bitarrayLoadElements)
                );
                if(threadIdx.x == 0){
                    smem.inputBitarray[fetchStage % numStages][0] = 0xFFFFFFFF;
                }
            }else{
                cooperative_groups::memcpy_async(
                    block,
                    (char*)(&smem.inputBitarray[fetchStage % numStages][0]),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * bitarrayLoadElements),
                    (const char*)(d_bitarrayFlags + bitarrayLoadBegin),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * bitarrayLoadElements)
                );
            }
            if(isLastChunk){
                if(threadIdx.x == 0){
                    smem.inputBitarray[fetchStage % numStages][bitarrayLoadElements] = 0xFFFFFFFF;
                }
            }
            

            fetchStage++;
        }

        // const bool isFirstChunk = chunk == 0;
        // const bool isLastChunk = chunk == numChunks-1;

        const int64_t chunkBegin = chunk * maxChunkSize;
        const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);

        //wait for async load from gmem to smem of our input chunk
        if(chunk < numChunks - 1 - int(gridDim.x)){
            //each cg::memcpy_async commits immediately, so it is 2 outstanding calls per fetch
            cooperative_groups::wait_prior<2>(block);
        }else{
            cooperative_groups::wait(block);
        }

        /*
            Use bit-level operations to count run length of run containing current position.
            In bitarray, runs begin with a 1 (inclusive) and end with a 1 (exclusive)

            Example: 0b0000'0000'0000'0000'0000'0100'0101'1111
                                                   ^
                                                   |
                                           computing run length of this position

            ffs gives the length of portion left of current position (including current position)
            ffs(????'????'0000'0000'0000'0000'0000'0100) = 3

            clz gives the length of portion right of current position
            clz(0101'1111'1111'1111'1111'1111'1111'1111) = 1
            
            3+1 = 4
        */

        // if(chunk == 0){
        //     if(threadIdx.x == 0){
        //         printf("in kernel. flags: %u %u %u\n",
        //             smem.inputBitarray[computeStage % numStages][0],
        //             smem.inputBitarray[computeStage % numStages][1],
        //             smem.inputBitarray[computeStage % numStages][2]
        //         );
        //     }
        // }

        //

        const int chunkSize = chunkEnd - chunkBegin;
        for(int i = threadIdx.x; i < chunkSize; i += blockDim.x){
            
            const int bitPositionIn32 = 31-(i % 32);
            const int indexOf32 = i / 32;
            
            const unsigned int centerUint = smem.inputBitarray[computeStage % numStages][1+indexOf32];
            const unsigned int leftUint = smem.inputBitarray[computeStage % numStages][indexOf32];
            const unsigned int rightUint = smem.inputBitarray[computeStage % numStages][1+indexOf32+1];
            
            const unsigned int ffsinput = __funnelshift_r(centerUint, leftUint, bitPositionIn32);
            const int leftAndSelf = ffsinput == 0 ? 32 : __ffs(ffsinput);

            const unsigned int clzinput = __funnelshift_lc(rightUint, centerUint, 32-bitPositionIn32);
            const int right = __clz(clzinput);
            const int count = leftAndSelf + right;

            // if(chunk == 0){
            //     if(i < 4){
            //         printf("indexOf32 %d, bitPositionIn32 %d, leftAndSelf %d, right %d, centerUint %u, rightUint %u, clzinput %u\n", 
            //             indexOf32, bitPositionIn32, leftAndSelf, right, centerUint, rightUint, clzinput);
            //     }
            // }

                // if(chunkBegin + i == 134217683){
                //     printf("global position %ld, letter %d, indexOf32 %d, bitPositionIn32 %d, position flag is %u, leftAndSelf %d, right %d, "
                //         "leftUint %u, centerUint %u, rightUint %u, ffsinput %u, clzinput %u\n", 
                //         chunkBegin + i,
                //         int(smem.inputElements[computeStage % numStages][i]),
                //         indexOf32, bitPositionIn32, (centerUint >> bitPositionIn32) & 1, 
                //         leftAndSelf, right, leftUint, centerUint, rightUint, ffsinput, clzinput);
                // }
            
            if(count >= maskingThreshold){
                smem.inputElements[computeStage % numStages][i] = maskingLetter;
            }           
        }
        //wait for smem.inputElements
        block.sync();

        //copy masked input elements from smem to gmem
        cooperative_groups::memcpy_async(
            block,
            (char*)(d_outputData + chunkBegin),
            cuda::aligned_size_t<4>(sizeof(char) * chunkSize),
            (const char*)&smem.inputElements[computeStage % numStages][0],
            cuda::aligned_size_t<4>(sizeof(char) * chunkSize)
        );


        //chunk is done, can re-use smem
        block.sync();
    }
}

#endif


#if !defined(FOOBARAOS)

template<int64_t maxChunkSize>
__global__
void masking10_UseFlagsKernel2(
    char* __restrict__ d_outputData,
    const uint32_t* __restrict__ d_bitarrayFlags,
    int64_t totalInputElements,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter
){
    static_assert(maxChunkSize % 16 == 0);
    assert(maskingThreshold <= 32);


    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<32>(block);

    const int64_t numChunks = cuda::ceil_div(totalInputElements, maxChunkSize);

    constexpr int numStages = 2;
    constexpr int numUint32ForMaxChunkSize = 1 + cuda::ceil_div(maxChunkSize, 32) + 1;
    constexpr int smemBitArrayNumUint32Padded = cuda::ceil_div(numUint32ForMaxChunkSize, 4) * 4; //pad to 16 bytes

    __shared__ alignas(16) struct Smem{
        uint32_t inputBitarray[numStages][smemBitArrayNumUint32Padded];
    } smem;

    const int bitarrayTotalElements = cuda::ceil_div(totalInputElements, 32);
    int computeStage = 0;
    int fetchStage = 0;
    for(int chunk = blockIdx.x, fetchChunk = blockIdx.x; chunk < numChunks; chunk += gridDim.x, computeStage++){
        for (; fetchChunk < numChunks && fetchChunk < (chunk + numStages * gridDim.x); fetchChunk += gridDim.x) {
            const bool isFirstChunk = fetchChunk == 0;
            const bool isLastChunk = fetchChunk == numChunks-1;

            const int64_t chunkBegin = fetchChunk * maxChunkSize;
            const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);


            //load the input bitarray of chunk
            const int chunkSize = chunkEnd - chunkBegin;
            const int bitarrayBegin = chunkBegin / 32;
            const int bitarrayEnd = bitarrayBegin + cuda::ceil_div(chunkSize, 32);
            const int bitarrayLoadBegin = isFirstChunk ? bitarrayBegin : bitarrayBegin - 1;
            const int bitarrayLoadEnd = cuda::std::min(bitarrayTotalElements, bitarrayEnd + 1);
            const int bitarrayLoadElements = bitarrayLoadEnd - bitarrayLoadBegin;
            if(isFirstChunk){
                cooperative_groups::memcpy_async(
                    block,
                    (char*)(&smem.inputBitarray[fetchStage % numStages][1]),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * bitarrayLoadElements),
                    (const char*)(d_bitarrayFlags + bitarrayLoadBegin),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * bitarrayLoadElements)
                );
                if(threadIdx.x == 0){
                    smem.inputBitarray[fetchStage % numStages][0] = 0xFFFFFFFF;
                }
            }else{
                cooperative_groups::memcpy_async(
                    block,
                    (char*)(&smem.inputBitarray[fetchStage % numStages][0]),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * bitarrayLoadElements),
                    (const char*)(d_bitarrayFlags + bitarrayLoadBegin),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * bitarrayLoadElements)
                );
            }
            if(isLastChunk){
                if(threadIdx.x == 0){
                    smem.inputBitarray[fetchStage % numStages][bitarrayLoadElements] = 0xFFFFFFFF;
                }
            }
            

            fetchStage++;
        }

        // const bool isFirstChunk = chunk == 0;
        // const bool isLastChunk = chunk == numChunks-1;

        const int64_t chunkBegin = chunk * maxChunkSize;
        const int64_t chunkEnd = cuda::std::min(totalInputElements, chunkBegin + maxChunkSize);

        //wait for async load from gmem to smem of our input chunk
        if(chunk < numChunks - 1 - int(gridDim.x)){
            cooperative_groups::wait_prior<1>(block);
        }else{
            cooperative_groups::wait(block);
        }

        /*
            Use bit-level operations to count run length of run containing current position.
            In bitarray, runs begin with a 1 (inclusive) and end with a 1 (exclusive)

            Example: 0b0000'0000'0000'0000'0000'0100'0101'1111
                                                   ^
                                                   |
                                           computing run length of this position

            ffs gives the length of portion left of current position (including current position)
            ffs(????'????'0000'0000'0000'0000'0000'0100) = 3

            clz gives the length of portion right of current position
            clz(0101'1111'1111'1111'1111'1111'1111'1111) = 1
            
            3+1 = 4
        */

        // if(chunk == 0){
        //     if(threadIdx.x == 0){
        //         printf("in kernel. flags: %u %u %u\n",
        //             smem.inputBitarray[computeStage % numStages][0],
        //             smem.inputBitarray[computeStage % numStages][1],
        //             smem.inputBitarray[computeStage % numStages][2]
        //         );
        //     }
        // }

        //

        const int chunkSize = chunkEnd - chunkBegin;
        for(int i = threadIdx.x; i < chunkSize; i += blockDim.x){
            
            const int bitPositionIn32 = 31-(i % 32);
            const int indexOf32 = i / 32;
            
            const unsigned int centerUint = smem.inputBitarray[computeStage % numStages][1+indexOf32];
            const unsigned int leftUint = smem.inputBitarray[computeStage % numStages][indexOf32];
            const unsigned int rightUint = smem.inputBitarray[computeStage % numStages][1+indexOf32+1];
            
            const unsigned int ffsinput = __funnelshift_r(centerUint, leftUint, bitPositionIn32);
            const int leftAndSelf = ffsinput == 0 ? 32 : __ffs(ffsinput);

            const unsigned int clzinput = __funnelshift_lc(rightUint, centerUint, 32-bitPositionIn32);
            const int right = __clz(clzinput);
            const int count = leftAndSelf + right;

            // if(chunk == 0){
            //     if(i < 4){
            //         printf("indexOf32 %d, bitPositionIn32 %d, leftAndSelf %d, right %d, centerUint %u, rightUint %u, clzinput %u\n", 
            //             indexOf32, bitPositionIn32, leftAndSelf, right, centerUint, rightUint, clzinput);
            //     }
            // }

                // if(chunkBegin + i == 134217683){
                //     printf("global position %ld, letter %d, indexOf32 %d, bitPositionIn32 %d, position flag is %u, leftAndSelf %d, right %d, "
                //         "leftUint %u, centerUint %u, rightUint %u, ffsinput %u, clzinput %u\n", 
                //         chunkBegin + i,
                //         int(smem.inputElements[computeStage % numStages][i]),
                //         indexOf32, bitPositionIn32, (centerUint >> bitPositionIn32) & 1, 
                //         leftAndSelf, right, leftUint, centerUint, rightUint, ffsinput, clzinput);
                // }
            
            if(count >= maskingThreshold){
                (d_outputData + chunkBegin)[i] = maskingLetter;
            }           
        }

        //chunk is done, can re-use smem
        block.sync();
    }
}

#endif

template<int blocksize>
__global__
void masking10_UseFlagsKernel3_billo(
    char* __restrict__ d_outputData,
    const uint32_t* __restrict__ d_bitarrayFlags,
    int64_t totalInputElements,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter
){
    assert(maskingThreshold <= 32);
    __builtin_assume(blockDim.x == blocksize);

    const size_t bitarrayDataElements = cuda::ceil_div(totalInputElements, 32);

    for(int64_t i = threadIdx.x + blockIdx.x * blockDim.x; i < totalInputElements; i += blockDim.x * gridDim.x){
        const int bitPositionIn32 = 31-(i % 32);
        const int indexOf32 = i / 32;
        unsigned int leftUint = 0;
        unsigned int rightUint = 0xFFFFFFFF;
        const unsigned int centerUint = d_bitarrayFlags[indexOf32];
        if(indexOf32 > 0){
            leftUint = d_bitarrayFlags[indexOf32-1];
        }
        if(indexOf32 < bitarrayDataElements-1){
            rightUint = d_bitarrayFlags[indexOf32+1];
        }
        
        const unsigned int ffsinput = __funnelshift_r(centerUint, leftUint, bitPositionIn32);
        const int leftAndSelf = ffsinput == 0 ? 32 : __ffs(ffsinput);

        const unsigned int clzinput = __funnelshift_lc(rightUint, centerUint, 32-bitPositionIn32);
        const int right = __clz(clzinput);
        const int count = leftAndSelf + right;
        
        if(count >= maskingThreshold){
            d_outputData[i] = maskingLetter;
        }  
    }
}

template<int blocksize, int smemElements>
__global__
void masking10_UseFlagsKernel4(
    char* __restrict__ d_outputData,
    const uint32_t* __restrict__ d_bitarrayFlags,
    int bitarrayDataElements,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter
){
    static_assert(smemElements % 4 == 0); //assert 16 byte alignment
    assert(maskingThreshold <= 32);

    __builtin_assume(blockDim.x == blocksize);
    __builtin_assume(blockDim.y == 1);
    __builtin_assume(blockDim.z == 1);

    constexpr int warpsize = 32;
    constexpr int numWarpsInBlock = blocksize / warpsize;
    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<warpsize>(block);

    constexpr int maxChunkSize = smemElements - 2 * 4;
    const int64_t numChunks = cuda::ceil_div(bitarrayDataElements, maxChunkSize);

    constexpr int numStages = 2;

    // using WarpExchange = cub::WarpExchange<int, 32, warpsize>;

    __shared__ alignas(16) struct Smem{
        uint32_t inputBitarray[numStages][smemElements];
        //typename WarpExchange::TempStorage warpExchange[numWarpsInBlock];
        int communicate[numWarpsInBlock][warpsize * 32];
    } smem;

    int computeStage = 0;
    int fetchStage = 0;
    for(int chunk = blockIdx.x, fetchChunk = blockIdx.x; chunk < numChunks; chunk += gridDim.x, computeStage++){
        for (; fetchChunk < numChunks && fetchChunk < (chunk + numStages * gridDim.x); fetchChunk += gridDim.x) {
            const bool isFirstChunk = fetchChunk == 0;
            const bool isLastChunk = fetchChunk == numChunks-1;

            const int chunkBegin = fetchChunk * maxChunkSize;
            const int chunkEnd = cuda::std::min(bitarrayDataElements, chunkBegin + maxChunkSize);

            const int chunkLoadBegin = isFirstChunk ? chunkBegin : chunkBegin - 1;
            const int chunkLoadEnd = cuda::std::min(bitarrayDataElements, chunkEnd + 1);
            const int chunkLoadElements = chunkLoadEnd - chunkLoadBegin;
            if(isFirstChunk){
                cooperative_groups::memcpy_async(
                    block,
                    (char*)(&smem.inputBitarray[fetchStage % numStages][1]),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * chunkLoadElements),
                    (const char*)(d_bitarrayFlags + chunkLoadBegin),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * chunkLoadElements)
                );
                if(threadIdx.x == 0){
                    smem.inputBitarray[fetchStage % numStages][0] = 0xFFFFFFFF;
                }
            }else{
                cooperative_groups::memcpy_async(
                    block,
                    (char*)(&smem.inputBitarray[fetchStage % numStages][0]),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * chunkLoadElements),
                    (const char*)(d_bitarrayFlags + chunkLoadBegin),
                    cuda::aligned_size_t<4>(sizeof(uint32_t) * chunkLoadElements)
                );
            }
            if(isLastChunk){
                if(threadIdx.x == 0){
                    smem.inputBitarray[fetchStage % numStages][chunkLoadElements] = 0xFFFFFFFF;
                }
            }
            

            fetchStage++;
        }

        // const bool isFirstChunk = chunk == 0;
        // const bool isLastChunk = chunk == numChunks-1;

        const int chunkBegin = chunk * maxChunkSize;
        const int chunkEnd = cuda::std::min(bitarrayDataElements, chunkBegin + maxChunkSize);

        //wait for async load from gmem to smem of our input chunk
        if(chunk < numChunks - 1 - int(gridDim.x)){
            cooperative_groups::wait_prior<1>(block);
        }else{
            cooperative_groups::wait(block);
        }

        /*
            Use bit-level operations to count run length of run containing current position.
            In bitarray, runs begin with a 1 (inclusive) and end with a 1 (exclusive)

            Example: 0b0000'0000'0000'0000'0000'0100'0101'1111
                                                   ^
                                                   |
                                           computing run length of this position

            ffs gives the length of portion left of current position (including current position)
            ffs(????'????'0000'0000'0000'0000'0000'0100) = 3

            clz gives the length of portion right of current position
            clz(0101'1111'1111'1111'1111'1111'1111'1111) = 1
            
            3+1 = 4
        */

        const int chunkSize = chunkEnd - chunkBegin;
        char* const chunkOutputData = d_outputData + size_t(32) * chunkBegin;
        //iterate over uint32_t elements of bitarray
        const int loopLimit = cuda::ceil_div(chunkSize, warp.size()) * warp.size();
        for(int indexOf32 = threadIdx.x; indexOf32 < loopLimit; indexOf32 += blockDim.x){
            const unsigned int leftUint = (indexOf32 < chunkSize) ? smem.inputBitarray[computeStage % numStages][indexOf32] : 0;
            const unsigned int centerUint = (indexOf32 < chunkSize) ? smem.inputBitarray[computeStage % numStages][1+indexOf32] : 0;
            const unsigned int rightUint = (indexOf32 < chunkSize) ? smem.inputBitarray[computeStage % numStages][1+indexOf32+1] : 0;

            int count = 0;

            int numMasked = 0;
            int maskedPositions[32];
            if(indexOf32 < chunkSize){
                #pragma unroll
                for(int i = 0; i < 32; i++){

                    const int bitPositionIn32 = 31-i;
                    if(i == 0){
                        const unsigned int ffsinput = __funnelshift_r(centerUint, leftUint, bitPositionIn32);
                        const int leftAndSelf = ffsinput == 0 ? 32 : __ffs(ffsinput);

                        const unsigned int clzinput = __funnelshift_lc(rightUint, centerUint, 32-bitPositionIn32);
                        const int right = __clz(clzinput);
                        count = leftAndSelf + right;
                    }else{
                        uint32_t bit = (centerUint >> bitPositionIn32) & 1;
                        //if bit is 0, it is the same segment as previous bit -> same count
                        if(bit == 1){
                            const unsigned int ffsinput = __funnelshift_r(centerUint, leftUint, bitPositionIn32);
                            const int leftAndSelf = ffsinput == 0 ? 32 : __ffs(ffsinput);

                            const unsigned int clzinput = __funnelshift_lc(rightUint, centerUint, 32-bitPositionIn32);
                            const int right = __clz(clzinput);
                            count = leftAndSelf + right;
                        }
                    }

                    // if((chunkBegin + indexOf32) * 32 + i == 32052){
                    //     printf("global position %ld, indexOf32 %d, i %d, bitPositionIn32 %d, position flag is %u, leftAndSelf %d, right %d, "
                    //         "leftUint %u, centerUint %u, rightUint %u, ffsinput %u, clzinput %u\n", 
                    //         int64_t(chunkBegin + indexOf32) * 32 + i,
                    //         indexOf32, i, bitPositionIn32, (centerUint >> bitPositionIn32) & 1, 
                    //         leftAndSelf, right, leftUint, centerUint, rightUint, ffsinput, clzinput);
                    // }
                
                    if(count >= maskingThreshold){
                        //chunkOutputData[indexOf32 * 32 + i] = maskingLetter;
                        maskedPositions[i] = indexOf32 * 32 + i;
                        numMasked++;
                    }else{
                        maskedPositions[i] = -1;
                    }
                }
            }else{
                #pragma unroll
                for(int i = 0; i < 32; i++){
                    maskedPositions[i] = -1;
                }
            }

            const int exclPs = compatibility::cooperative_groups::exclusive_scan<32, blocksize>(warp, numMasked);
            const int numMaskedInWarp = compatibility::cooperative_groups::reduce(warp, numMasked, compatibility::cooperative_groups::plus<int>{});
            for(int i = 0, k = 0; i < 32; i++){
                if(maskedPositions[i] != -1){
                    smem.communicate[warp.meta_group_rank()][exclPs + k] = maskedPositions[i];
                    k++;
                }
            }
            warp.sync();
            for(int m = warp.thread_rank(); m < numMaskedInWarp; m += warp.size()){
                const int pos = smem.communicate[warp.meta_group_rank()][m];
                chunkOutputData[pos] = maskingLetter;
            }
            
            //WarpExchange(smem.warpExchange).BlockedToStriped(maskedPositions, maskedPositions);

            // for(int i = 0; i < 32; i++){
            //     if(maskedPositions[i] != -1){
            //         chunkOutputData[indexOf32 * 32 + i] = maskingLetter;
            //     }
            // }
        }

        //chunk is done, can re-use smem
        block.sync();
    }
}


template<class InputOffsetIterator>
__global__
void applyPaddingLetterKernel(
    char* __restrict__ d_data,
    const int* __restrict__ d_inputLengths,
    InputOffsetIterator d_inputOffsets,
    int numSegments,
    char paddingLetter
){
    const int tid = threadIdx.x + blockIdx.x * blockDim.x;
    const int segmentId = tid / 4;
    const int tidInSegment = tid % 4;

    if(segmentId < numSegments){

        const size_t segmentBeginOffset = d_inputOffsets[segmentId];
        const int segmentLength = d_inputLengths[segmentId];
        const int numPaddingBytes = (segmentLength % 4 == 0) ? 0 : 4 - (segmentLength % 4);
        if(tidInSegment < numPaddingBytes){
            d_data[segmentBeginOffset + segmentLength + tidInSegment] = paddingLetter;
        }
    }
}

/*
    d_outputData and d_inputData must be aligned to 16 bytes, and padded to 16 byte boundary

    All segments must be located one after each other without gaps. 
    Each segment must begin at 4 byte boundary. 
    padding bytes between two segments will be set to paddingLetter if applyExtraPaddingLetter == true
*/
template<class InputTransformOp, class InputOffsetIterator>
void maskSequences_consecutive_segments(
    bool inputAndOutArePadded,
    char* d_outputData,
    const char* d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    bool applyExtraPaddingLetter,
    char paddingLetter,
    InputTransformOp inputTransformOp,
    cudaStream_t stream
){
    if(!inputAndOutArePadded){
        throw std::runtime_error("inputAndOutArePadded padding is required for maskSequences_consecutive_segments");
    }
    if(maskingThreshold > 32){
        throw std::runtime_error("Cannot use maskSequences_consecutive_segments with maskingTreshold > 32");
    }

    const size_t bitarrayDataElements = cuda::ceil_div(totalInputElements, 32);
    thrust::device_vector<char, thrust_async_allocator<char>> d_bitarray_as_char 
        = allocateTempBuffer(sizeof(uint32_t) * bitarrayDataElements, stream);
    uint32_t* d_bitarray = (uint32_t*)d_bitarray_as_char.data().get();

    int deviceId;
    int numSMs;
    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));

    // {
    //     int maxBlocksPerSM;
    //     constexpr int blocksize = 128;
    //     auto kernel = masking10_compute_flags_kernel_3<blocksize, true, 4096, InputTransformOp>;
    //     CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
    //         &maxBlocksPerSM,
    //         kernel,
    //         blocksize, 
    //         0
    //     ));
    //     helpers::GpuTimer timer1("masking10_compute_flags_kernel_3");
    //     kernel<<<numSMs * maxBlocksPerSM, blocksize>>>(
    //         d_outputData,
    //         d_bitarray,
    //         d_inputData,
    //         totalInputElements,
    //         inputTransformOp
    //     );
    //     timer1.print();

    // }
    {
        int maxBlocksPerSM;
        constexpr int blocksize = 128;
        auto kernel = masking10_compute_flags_kernel_4<blocksize, true, 4096, InputTransformOp>;
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            0
        ));
        // helpers::GpuTimer timer1("masking10_compute_flags_kernel_4");
        kernel<<<numSMs * maxBlocksPerSM, blocksize, 0, stream>>>(
            d_outputData,
            d_bitarray,
            d_inputData,
            totalInputElements,
            inputTransformOp
        );
        // timer1.print();

    }

    {
        // helpers::GpuTimer timer0("set begin and oob flags");
        thrust::for_each(
            thrust::cuda::par_nosync.on(stream),
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(0) + numInputSequences,
            [
                d_inputOffsets,
                d_inputLengths,
                d_bitarray,
                numInputSequences,
                bitarrayDataElements,
                totalInputElements
            ] __device__ (int segmentId){
                const size_t segmentBeginOffset = d_inputOffsets[segmentId];
                const int segmentLength = d_inputLengths[segmentId];

                auto setBit = [&](size_t position){
                    const size_t dataIndex = position / 32;
                    const size_t bitIndex = 31 - (position % 32);
                    uint32_t value = 1u << bitIndex;
                    atomicOr(&d_bitarray[dataIndex], value);
                };

                //flag segment begin
                setBit(segmentBeginOffset);

                //flag padding bytes
                {
                    const int numPaddingBytes = (segmentLength % 4 == 0) ? 0 : 4 - (segmentLength % 4);
                    for(int b = 0; b < numPaddingBytes; b++){
                        setBit(segmentBeginOffset + segmentLength + b);
                    }
                }

                //set the unused bits at the end to 1
                if(segmentId == numInputSequences-1){
                    const size_t bitsAvailable = bitarrayDataElements * 32;
                    const size_t bitsUsed = totalInputElements;
                    const int bitsUnused = bitsAvailable - bitsUsed;
                    const uint32_t oobflags = (1u << bitsUnused) - 1u;
                    atomicOr(&d_bitarray[bitarrayDataElements-1], oobflags);
                }
            }
        );

        // timer0.print();
    }

    {
        #if 0
        {
        int maxBlocksPerSM;
        constexpr int blocksize = 128;
        auto kernel = masking10_UseFlagsKernel2<16384>;
        // auto kernel = masking10_UseFlagsKernel3_billo<blocksize>;
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            0
        ));
        // helpers::GpuTimer timer2("masking10_UseFlagsKernel2");
        kernel<<<numSMs * maxBlocksPerSM, blocksize, 0, stream>>>(
            d_outputData,
            d_bitarray,
            totalInputElements,
            maskingThreshold,
            maskingLetter
        );
        // timer2.print();
        }
        #endif

        #if 1
        {
            #if !defined(FOOBARAOS)
        int maxBlocksPerSM;
        constexpr int blocksize = 256;
        auto kernel = masking10_UseFlagsKernel4<blocksize, 2048>;
        // auto kernel = masking10_UseFlagsKernel3_billo<blocksize>;
        CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &maxBlocksPerSM,
            kernel,
            blocksize, 
            0
        ));
        // helpers::GpuTimer timer2("masking10_UseFlagsKernel4");
        kernel<<<numSMs * maxBlocksPerSM, blocksize, 0, stream>>>(
            d_outputData,
            d_bitarray,
            bitarrayDataElements,
            maskingThreshold,
            maskingLetter
        );
        #else
        assert(false); 
        std::cout << "ERROR MISSING KERNEL AMD\n";
        #endif
        // timer2.print();
        }
        #endif
    }

    if(applyExtraPaddingLetter){
        applyPaddingLetterKernel<<<cuda::ceil_div(numInputSequences, 128/4), 128, 0, stream>>>(
            d_outputData,
            d_inputLengths,
            d_inputOffsets,
            numInputSequences,
            paddingLetter
        );
    }
        
}















template<int blocksize, int maximumChunkSize, int maximumMaskingThreshold, class InputTransformOp, class InputOffsetIterator, class OutputOffsetIterator>
__global__
void maskingKernel5555_1warpPerSegment(
    char* __restrict__ d_outputData,
    OutputOffsetIterator d_outputOffsets,
    const char* __restrict__ d_inputData,
    const int* __restrict__ d_lengths,
    InputOffsetIterator d_offsets,
    int64_t totalInputElements,
    int numSegments,
    int maskingThreshold,
    char maskingLetter,
    char paddingLetter,
    InputTransformOp inputTransformOp
){
    __builtin_assume(blocksize == blockDim.x);
    __builtin_assume(maskingThreshold <= maximumMaskingThreshold);

    static_assert(maximumMaskingThreshold <= 16);
    static_assert(maximumChunkSize % 16 == 0);

    auto warp = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());
    constexpr int numWarpsInBlock = blocksize / 32;
    const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
    const int warpIdInBlock = warp.meta_group_rank();
    const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;

    // constexpr bool isIdentityTransform = cuda::std::is_same_v<InputTransformOp, cuda::std::identity>;

    constexpr int vecLoadPadding = 16;
    constexpr int maskingHaloSize = 16;
    const int maxChunkWithHaloSize = maskingHaloSize+maximumChunkSize+maskingHaloSize;

    __shared__ alignas(16) struct Smem{
        char chunkElements[numWarpsInBlock][maxChunkWithHaloSize + vecLoadPadding];
        char chunkFlags[numWarpsInBlock][maxChunkWithHaloSize];
        //uint32_t bitarrayFlags[numWarpsInBlock][1+cuda::ceil_div(maximumChunkSize,32)+1];
    } smem;


    for(int segmentId = warpIdInGrid; segmentId < numSegments; segmentId += numWarpsInGrid){
        const size_t segmentInputOffset = d_offsets[segmentId];
        const size_t segmentOutputOffset = d_outputOffsets[segmentId];
        const int segmentSize = d_lengths[segmentId];
        const char* const segmentInputData = cuda::std::assume_aligned<4>(d_inputData + segmentInputOffset);
        char* const segmentOutputData = cuda::std::assume_aligned<4>(d_outputData + segmentOutputOffset);
        // const int smemVecLoadOffset = uintptr_t(segmentInputData) % 16 == 0 ? 0 : 16 - uintptr_t(segmentInputData) % 16;

        //divide output range into chunks
        const int numChunks = cuda::ceil_div(segmentSize, maximumChunkSize);

        constexpr int numChar4toStoreMaxChunkWithHaloSize = cuda::ceil_div(maxChunkWithHaloSize, 4);
        constexpr int numChar4PerThreadtoStoreMaxChunkWithHaloSize = cuda::ceil_div(numChar4toStoreMaxChunkWithHaloSize, warp.size());

        uint32_t nextLoadedInputs[2][numChar4PerThreadtoStoreMaxChunkWithHaloSize];
        for(int prefetchDistance = 0; prefetchDistance < 2; prefetchDistance++){
            const int next_chunk = prefetchDistance;
            if(next_chunk < numChunks){
                const bool next_isFirstChunk = next_chunk == 0;
                const bool next_isLastChunk = next_chunk == numChunks-1;
    
                const int next_chunkBegin = next_chunk * maximumChunkSize;
                const int next_chunkEnd = cuda::std::min(segmentSize, next_chunkBegin + maximumChunkSize);
                // const int next_chunkSize = next_chunkEnd - next_chunkBegin;
                
                const int next_chunkLoadBegin = next_isFirstChunk ? next_chunkBegin : next_chunkBegin - maskingHaloSize;
                const int next_chunkLoadEnd = (!next_isLastChunk) ? next_chunkEnd + maskingHaloSize : cuda::round_up(next_chunkEnd, 4);
                const int next_chunkLoadSize = next_chunkLoadEnd - next_chunkLoadBegin;
                const int limit = cuda::ceil_div(next_chunkLoadSize, 4);
                for(int i = 0; i < numChar4PerThreadtoStoreMaxChunkWithHaloSize; i++){
                    const int index = i * warp.num_threads() + warp.thread_rank();
                    if(index < limit){
                        // if(segmentId == 0){
                        //     printf("prefetch next_chunk %d, thread %d, index %d, load from %d\n", 
                        //         next_chunk, warp.thread_rank(), index, next_chunkLoadBegin + 4*index
                        //     );
                        // }
                        nextLoadedInputs[prefetchDistance][i] = *((const uint32_t*)&segmentInputData[next_chunkLoadBegin + 4*index]);
                    }
                }
            }
        }

        for(int chunk = 0; chunk < numChunks; chunk++){
            const bool isFirstChunk = chunk == 0;
            const bool isLastChunk = chunk == numChunks-1;

            const int chunkBegin = chunk * maximumChunkSize;
            const int chunkEnd = cuda::std::min(segmentSize, chunkBegin + maximumChunkSize);
            const int chunkSize = chunkEnd - chunkBegin;
            
            const int chunkLoadBegin = isFirstChunk ? chunkBegin : chunkBegin - maskingHaloSize;
            const int chunkLoadEnd = (!isLastChunk) ? chunkEnd + maskingHaloSize : cuda::round_up(chunkEnd, 4);
            const int chunkLoadSize = chunkLoadEnd - chunkLoadBegin;

            //wait before re-using smem
            warp.sync();

            #if 1

                uint32_t currentLoadedInputs[numChar4PerThreadtoStoreMaxChunkWithHaloSize];
                if(chunk % 2 == 0){
                    for(int i = 0; i < numChar4PerThreadtoStoreMaxChunkWithHaloSize; i++){
                        currentLoadedInputs[i] = nextLoadedInputs[0][i];
                    }
                }else{
                    for(int i = 0; i < numChar4PerThreadtoStoreMaxChunkWithHaloSize; i++){
                        currentLoadedInputs[i] = nextLoadedInputs[1][i];
                    }
                }
                //fire up loading of next chunk
                {
                    const int next_chunk = chunk+2;
                    if(next_chunk < numChunks){
                        const bool next_isFirstChunk = next_chunk == 0;
                        const bool next_isLastChunk = next_chunk == numChunks-1;
            
                        const int next_chunkBegin = next_chunk * maximumChunkSize;
                        const int next_chunkEnd = cuda::std::min(segmentSize, next_chunkBegin + maximumChunkSize);
                        // const int next_chunkSize = next_chunkEnd - next_chunkBegin;
                        
                        const int next_chunkLoadBegin = next_isFirstChunk ? next_chunkBegin : next_chunkBegin - maskingHaloSize;
                        const int next_chunkLoadEnd = (!next_isLastChunk) ? next_chunkEnd + maskingHaloSize : cuda::round_up(next_chunkEnd, 4);
                        const int next_chunkLoadSize = next_chunkLoadEnd - next_chunkLoadBegin;
                        const int limit = cuda::ceil_div(next_chunkLoadSize, 4);
                        if(next_chunk % 2 == 0){
                            for(int i = 0; i < numChar4PerThreadtoStoreMaxChunkWithHaloSize; i++){
                                const int index = i * warp.num_threads() + warp.thread_rank();
                                if(index < limit){
                                    nextLoadedInputs[0][i] = *((const uint32_t*)&segmentInputData[next_chunkLoadBegin + 4*index]);
                                }
                            }
                        }else{
                            for(int i = 0; i < numChar4PerThreadtoStoreMaxChunkWithHaloSize; i++){
                                const int index = i * warp.num_threads() + warp.thread_rank();
                                if(index < limit){
                                    nextLoadedInputs[1][i] = *((const uint32_t*)&segmentInputData[next_chunkLoadBegin + 4*index]);
                                }
                            }
                        }
                    }
                }


                for(int i = 0; i < numChar4PerThreadtoStoreMaxChunkWithHaloSize; i++){
                    const int index = i * warp.num_threads() + warp.thread_rank();
                    const int limit = cuda::ceil_div(chunkLoadSize, 4);
                    if(index < limit){
                        char4 tmp = *((const char4*)&currentLoadedInputs[i]);
                        tmp.x = inputTransformOp(tmp.x);
                        tmp.y = inputTransformOp(tmp.y);
                        tmp.z = inputTransformOp(tmp.z);
                        tmp.w = inputTransformOp(tmp.w);
                        *((char4*)&smem.chunkElements[warpIdInBlock][maskingHaloSize + (chunkLoadBegin+4*index - chunkBegin)]) = tmp;
                    }
                }

                // //transform and store into smem
                // for(int i = chunkLoadBegin + 4*warp.thread_rank(); i < chunkLoadEnd; i += 4*warp.num_threads()){
                //     if(segmentId == 0){
                //         printf("chunk %d normal. thread %d, load from %d, store to %d\n", 
                //             chunk, warp.thread_rank(), i, maskingHaloSize + (i - chunkBegin)
                //         );
                //     }
                //     char4 tmp = *((const char4*)&segmentInputData[i]);
                //     tmp.x = inputTransformOp(tmp.x);
                //     tmp.y = inputTransformOp(tmp.y);
                //     tmp.z = inputTransformOp(tmp.z);
                //     tmp.w = inputTransformOp(tmp.w);
                //     *((char4*)&smem.chunkElements[warpIdInBlock][maskingHaloSize + (i - chunkBegin)]) = tmp;
                // }


                
                warp.sync();

                //compute chunk flags
                int vecchunks = maxChunkWithHaloSize / 4;
                using VecT = uint;
                for(int c = warp.thread_rank(); c < vecchunks; c += warp.num_threads()){
                    VecT current = ((const VecT*)&smem.chunkElements[warpIdInBlock][0])[c];
                    alignas(sizeof(VecT)) char myElementsAsChars[sizeof(VecT)];
                    memcpy(&myElementsAsChars[0], &current, sizeof(VecT));

                    char leftNeighborElement = 0;
                    if(c > 0){
                        leftNeighborElement = smem.chunkElements[warpIdInBlock][(sizeof(VecT) * c) - 1];
                    }

                    alignas(sizeof(VecT)) char currentflags[sizeof(VecT)];

                    currentflags[0] = leftNeighborElement != myElementsAsChars[0];
                    for(int i = 1; i < sizeof(VecT); i++){
                        currentflags[i] = (myElementsAsChars[i-1] != myElementsAsChars[i]);
                    }
                    VecT flagsvec;
                    memcpy(&flagsvec, &currentflags[0], sizeof(VecT));
                    ((VecT*)&smem.chunkFlags[warpIdInBlock][0])[c] = flagsvec;      
                }
                warp.sync();
                if(isFirstChunk){
                    smem.chunkFlags[warpIdInBlock][maskingHaloSize + (chunkBegin - chunkLoadBegin)] = 1; //mark begin of whole segment
                }
                if(isLastChunk){
                    if(maskingHaloSize + chunkSize + warp.thread_rank() < maxChunkWithHaloSize){
                        smem.chunkFlags[warpIdInBlock][maskingHaloSize + chunkSize + warp.thread_rank()] = 1; //mark oob elements;
                    }
                }
                warp.sync();


                //counting
                for(int c4 = warp.thread_rank(); c4 < cuda::ceil_div(chunkSize,4); c4 += warp.num_threads()){
                    uint center = ((uint*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4];
                    alignas(16) uint left4[4];
                    // if(c4 % 4 == 0){
                    //     *((uint4*)&left4) = ((uint4*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4-4];
                    // }else{
                        left4[0] = ((uint*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4-4];
                        left4[1] = ((uint*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4-3];
                        left4[2] = ((uint*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4-2];
                        left4[3] = ((uint*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4-1];
                    // }
                    alignas(16) uint right4[4];
                    // if(c4 % 4 == 3){
                    //     *((uint4*)&right4) = ((uint4*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4+1];
                    // }else{
                        right4[0] = ((uint*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4+1];
                        right4[1] = ((uint*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4+2];
                        right4[2] = ((uint*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4+3];
                        right4[3] = ((uint*)&smem.chunkFlags[warpIdInBlock][maskingHaloSize])[c4+4];
                    // }
                    
                    alignas(16) char chararray[36];
                    memcpy(&chararray[0], &left4[0], 16);
                    memcpy(&chararray[16], &center, 4);
                    memcpy(&chararray[20], &right4[0], 16);

                    int counts[4];

                    //for first element, count left
                    counts[0] = 0;
                    for(int i = 0; i < maximumMaskingThreshold; i++){
                        counts[0]++;
                        if(chararray[16 - i] == 1){
                            break;
                        }
                    }
                    //for first element, count right
                    for(int i = 0; i < maximumMaskingThreshold-1; i++){
                        if(chararray[16 + 1 + i] == 1){
                            break;
                        }
                        counts[0]++;
                    }
                    for(int k = 1; k < 4; k++){
                        if(chararray[16+k] == 0){
                            counts[k] = counts[k-1]; //same segment, same masking criterion
                        }else{
                            //new segment begins at k, count to the right
                            counts[k] = 1;
                            for(int i = 0; i < maximumMaskingThreshold-1; i++){
                                if(chararray[16 + 1 + k + i] == 1){
                                    break;
                                }
                                counts[k]++;
                            }
                        }
                    }
                    
                    for(int k = 0; k < 4; k++){
                        if(counts[k] >= maskingThreshold){
                            smem.chunkElements[warpIdInBlock][maskingHaloSize + c4 * 4 + k] = maskingLetter;
                        }
                    }
                }
                warp.sync();


                const int full = chunkSize / 4;
                const int remainder = chunkSize % 4;
                for(int i = warp.thread_rank(); i < full; i += warp.num_threads()){
                    ((uint*)&segmentOutputData[chunkBegin])[i] = ((const uint*)&smem.chunkElements[warpIdInBlock][maskingHaloSize])[i];
                }
                if(warp.thread_rank() < remainder){                    
                    segmentOutputData[chunkBegin + 4 * full + warp.thread_rank()] = smem.chunkElements[warpIdInBlock][maskingHaloSize + 4 * full + warp.thread_rank()];
                }
                // for(int outputIndex = warp.thread_rank(); outputIndex < chunkSize; outputIndex += warp.num_threads()){
                //     const int indexInSmem = 16 + outputIndex;
                //     segmentOutputData[chunkBegin + outputIndex] = smem.chunkElements[warpIdInBlock][indexInSmem];
                // }
                warp.sync();

            #endif

            #if 0

                //transform and store into smem
                for(int i = chunkLoadBegin + 4*warp.thread_rank(); i < chunkLoadEnd; i += 4*warp.num_threads()){
                    char4 tmp = *((const char4*)&segmentInputData[i]);
                    tmp.x = inputTransformOp(tmp.x);
                    tmp.y = inputTransformOp(tmp.y);
                    tmp.z = inputTransformOp(tmp.z);
                    tmp.w = inputTransformOp(tmp.w);
                    *((char4*)&smem.chunkElements[warpIdInBlock][maskingHaloSize + (i - chunkBegin)]) = tmp;
                }

                
                warp.sync();

                //compute chunk flags
                int vecchunks = maxChunkWithHaloSize / 4;
                using VecT = uint;
                for(int c = warp.thread_rank(); c < vecchunks; c += warp.num_threads()){
                    VecT current = ((const VecT*)&smem.chunkElements[warpIdInBlock][0])[c];
                    alignas(sizeof(VecT)) char myElementsAsChars[sizeof(VecT)];
                    memcpy(&myElementsAsChars[0], &current, sizeof(VecT));

                    char leftNeighborElement = 0;
                    if(c > 0){
                        leftNeighborElement = smem.chunkElements[warpIdInBlock][(sizeof(VecT) * c) - 1];
                    }

                    alignas(sizeof(VecT)) char currentflags[sizeof(VecT)];

                    currentflags[0] = leftNeighborElement != myElementsAsChars[0];
                    for(int i = 1; i < sizeof(VecT); i++){
                        currentflags[i] = (myElementsAsChars[i-1] != myElementsAsChars[i]);
                    }
                    VecT flagsvec;
                    memcpy(&flagsvec, &currentflags[0], sizeof(VecT));
                    ((VecT*)&smem.chunkFlags[warpIdInBlock][0])[c] = flagsvec;      
                }
                warp.sync();
                if(isFirstChunk){
                    smem.chunkFlags[warpIdInBlock][maskingHaloSize + (chunkBegin - chunkLoadBegin)] = 1; //mark begin of whole segment
                }
                if(isLastChunk){
                    if(maskingHaloSize + chunkSize + warp.thread_rank() < maxChunkWithHaloSize){
                        smem.chunkFlags[warpIdInBlock][maskingHaloSize + chunkSize + warp.thread_rank()] = 1; //mark oob elements;
                    }
                }

                //counting
                for(int outputIndex = warp.thread_rank(); outputIndex < chunkSize; outputIndex += warp.num_threads()){
                    const int indexInSmem = maskingHaloSize + outputIndex;

                    // const char myFlag = (smem.chunkFlags[warpIdInBlock][indexInSmem]);



                    //check to the left and self
                    int indexOfFirst1LeftCenter = -1;
                    int loopLimit1 = isFirstChunk ? maskingHaloSize : cuda::std::max(0, indexInSmem - maskingThreshold + 1);
                    for(int i = indexInSmem; i >= loopLimit1; i--){
                        //const char otherLetter = inputTransformOp(smem.chunkElements[warpIdInBlock][i]);
                        if(indexOfFirst1LeftCenter == -1){
                            const char flag = (smem.chunkFlags[warpIdInBlock][i]);
                            if(flag == 1){
                                indexOfFirst1LeftCenter = i;
                            }
                        }
                    }

                    int count = indexOfFirst1LeftCenter == -1 ? loopLimit1 : indexInSmem - indexOfFirst1LeftCenter + 1;

                    //check to the right
                    int loopLimit2 = isLastChunk ? indexInSmem + maskingThreshold : indexInSmem + maskingThreshold;
                    loopLimit2 = cuda::std::min(loopLimit2, indexInSmem + 1 + (maskingThreshold - count));
                    for(int i = indexInSmem + 1; i < loopLimit2; i++){
                        //const char otherLetter = inputTransformOp(smem.chunkElements[warpIdInBlock][i]);
                        const char flag = (smem.chunkFlags[warpIdInBlock][i]);
                        if(flag == 0){
                            count++;
                        }else{
                            break;
                        }
                    }

                    if(count >= maskingThreshold){
                        smem.chunkElements[warpIdInBlock][indexInSmem] = maskingLetter;
                    }
                }



                // cooperative_groups::memcpy_async(
                //     warp,
                //     (char*)(&segmentOutputData[chunkBegin]),
                //     (const char*)(&smem.chunkFlags[warpIdInBlock][16]),
                //     sizeof(char) * chunkSize
                // );
                const int full = chunkSize / 4;
                const int remainder = chunkSize % 4;
                for(int i = warp.thread_rank(); i < full; i += warp.num_threads()){
                    ((uint*)&segmentOutputData[chunkBegin])[i] = ((const uint*)&smem.chunkElements[warpIdInBlock][maskingHaloSize])[i];
                }
                if(warp.thread_rank() < remainder){                    
                    segmentOutputData[chunkBegin + 4 * full + warp.thread_rank()] = smem.chunkElements[warpIdInBlock][maskingHaloSize + 4 * full + warp.thread_rank()];
                }
                // for(int outputIndex = warp.thread_rank(); outputIndex < chunkSize; outputIndex += warp.num_threads()){
                //     const int indexInSmem = 16 + outputIndex;
                //     segmentOutputData[chunkBegin + outputIndex] = smem.chunkFlags[warpIdInBlock][indexInSmem];
                // }
                warp.sync();

            #endif

            #if 0
                for(int i = chunkLoadBegin + 4*warp.thread_rank(); i < chunkLoadEnd; i += 4*warp.num_threads()){
                    char4 tmp = *((const char4*)&segmentInputData[i]);
                    tmp.x = inputTransformOp(tmp.x);
                    tmp.y = inputTransformOp(tmp.y);
                    tmp.z = inputTransformOp(tmp.z);
                    tmp.w = inputTransformOp(tmp.w);
                    *((char4*)&smem.chunkElements[warpIdInBlock][16 + (i - chunkBegin)]) = tmp;
                }
                warp.sync();
                for(int outputIndex = warp.thread_rank(); outputIndex < chunkSize; outputIndex += warp.num_threads()){
                    const int indexInSmem = 16 + outputIndex;

                    const char myLetter = (smem.chunkElements[warpIdInBlock][indexInSmem]);

                    int count = 1; //account for my letter

                    //check to the left
                    int loopLimit1 = isFirstChunk ? 16 : cuda::std::max(0, indexInSmem - maskingThreshold + 1);
                    for(int i = indexInSmem - 1; i >= loopLimit1; i--){
                        //const char otherLetter = inputTransformOp(smem.chunkElements[warpIdInBlock][i]);
                        const char otherLetter = (smem.chunkElements[warpIdInBlock][i]);
                        if(myLetter == otherLetter){
                            count++;
                        }else{
                            break;
                        }
                    }

                    //check to the right
                    int loopLimit2 = isLastChunk ? 16 + (chunkLoadEnd - chunkBegin) : indexInSmem + maskingThreshold;
                    loopLimit2 = cuda::std::min(loopLimit2, indexInSmem + 1 + (maskingThreshold - count));
                    for(int i = indexInSmem + 1; i < loopLimit2; i++){
                        //const char otherLetter = inputTransformOp(smem.chunkElements[warpIdInBlock][i]);
                        const char otherLetter = (smem.chunkElements[warpIdInBlock][i]);
                        if(myLetter == otherLetter){
                            count++;
                            if(count >= maskingThreshold){
                                break;
                            }
                        }else{
                            break;
                        }
                    }


                    if(count >= maskingThreshold){
                        segmentOutputData[chunkBegin + outputIndex] = maskingLetter;
                    }else{
                        segmentOutputData[chunkBegin + outputIndex] = myLetter;
                    }
                }
                if(isLastChunk){
                    //apply 4-byte padding
                    const int paddedLength = cuda::ceil_div(segmentSize, 4) * 4;
                    if(segmentSize + warp.thread_rank() < paddedLength){
                        segmentOutputData[segmentSize + warp.thread_rank()] = paddingLetter;
                    }
                }
            #endif
            
        }



    }
   
}



/*
    Requirements:
        sequence[i] has length inputLengths[i] and begins at input position inputOffsets[i]
        sequence[i] is padded to 4 byte boundary
        d_outputData must not overlap d_inputData
        d_outputData must be able to contain 4-byte padded sequences
*/
template<int maximumChunkSize, int maximumMaskingThreshold, class InputTransformOp, class InputOffsetIterator, class OutputOffsetIterator>
void call_maskingKernel5555_1warpPerSegment(
    char* d_outputData,
    OutputOffsetIterator d_outputOffsets,
    const char* d_inputData, 
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    bool applyExtraPaddingLetter,
    char paddingLetter, //to pad output sequences to 4 byte boundary if applyExtraPaddingLetter == 1
    InputTransformOp inputTransformOp, // input letters are given by inputTransformOp(d_inputData[i])
    cudaStream_t stream
){
    static_assert(maximumMaskingThreshold <= 16);

    if(maskingThreshold > 16){
        throw std::runtime_error("maskingThreshold too large, limit 16\n");
    }

    int deviceId;
    int numSMs;
    int maxBlocksPerSM = 0;
    constexpr int blocksize = 256;
    auto kernel = maskingKernel5555_1warpPerSegment<blocksize, maximumChunkSize, maximumMaskingThreshold,
        InputTransformOp, InputOffsetIterator, OutputOffsetIterator>;

    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));
    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        0
    ));

    const int requiredNumBlocks = cuda::ceil_div(numInputSequences, blocksize / 32);
    const int numBlocks = std::min(requiredNumBlocks, numSMs*maxBlocksPerSM);

    // std::cout << "numBlocks " << numBlocks << "\n";

    kernel<<<numBlocks, blocksize, 0, stream>>>(
        d_outputData,
        d_outputOffsets,
        d_inputData,
        d_inputLengths,
        d_inputOffsets,
        totalInputElements,
        numInputSequences,
        maskingThreshold,
        maskingLetter,
        paddingLetter,
        inputTransformOp
    );
    CUDACHECKASYNC;

    if(applyExtraPaddingLetter){
        applyPaddingLetterKernel<<<cuda::ceil_div(numInputSequences, 128/4), 128, 0, stream>>>(
            d_outputData,
            d_inputLengths,
            d_outputOffsets,
            numInputSequences,
            paddingLetter
        );
    }
}


template<class InputTransformOp, class InputOffsetIterator, class OutputOffsetIterator>
void maskSequences5555(
    char* d_outputData,
    OutputOffsetIterator d_outputOffsets,
    const char* d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    bool applyExtraPaddingLetter,
    char paddingLetter, //to pad output sequences to 4 byte boundary if applyExtraPaddingLetter == 1
    InputTransformOp inputTransformOp,
    cudaStream_t stream
){
    #define CALLKERNEL(chunksize, threshold){ \
        call_maskingKernel5555_1warpPerSegment<chunksize, threshold>( \
            d_outputData, \
            d_outputOffsets, \
            d_inputData,  \
            d_inputLengths,  \
            d_inputOffsets, \
            totalInputElements, \
            numInputSequences, \
            maskingThreshold, \
            maskingLetter, \
            applyExtraPaddingLetter, \
            paddingLetter, \
            inputTransformOp, \
            stream \
        ); \
    }
    if(maskingThreshold == 6){
        CALLKERNEL(128,6);
    }else{
        CALLKERNEL(128,16);
    }
    #undef CALLKERNEL
}











template<class InputTransformOp, class InputDataIterator, class InputOffsetIterator, class OutputOffsetIterator>
void maskSequences(
    char* d_outputData,
    OutputOffsetIterator d_outputOffsets,
    InputDataIterator d_inputData, //returns the letter (char) to consider for masking. each letter must be in range [0, 127]
    const int* d_inputLengths, //size numInputSequences
    InputOffsetIterator d_inputOffsets, //size numInputSequences
    size_t totalInputElements,
    int numInputSequences,
    int maskingThreshold, // replace letters in all repeats >= maskingThreshold by masking maskingLetter
    char maskingLetter,
    char paddingLetter,
    InputTransformOp inputTransformOp,
    cudaStream_t stream
){

    maskSequences_custom_kernel_impl5(
        d_outputData,
        d_outputOffsets,
        d_inputData,
        d_inputLengths,
        d_inputOffsets,
        totalInputElements,
        numInputSequences,
        maskingThreshold,
        maskingLetter,
        paddingLetter,
        inputTransformOp,
        stream
    );
}








namespace maskingsinglethread{

template<int numWarpsInBlock, int maxUintsPerThread>
struct Kernel1Smem{
    uint32_t elements[numWarpsInBlock][maxUintsPerThread][32+1];
    // uint32_t outputelements[numWarpsInBlock][uintsPerThread][32+1];
    // int runBeginPositions[numWarpsInBlock][maxSegmentLength][32+1];
    int segmentLengths[numWarpsInBlock][32];
    size_t segmentOffsets[numWarpsInBlock][32];
};

template<int blocksize, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
__global__
void kernel1(
    char* __restrict__ d_outputData,
    const int* __restrict__ d_lengths,
    InputOffsetIterator d_offsets,
    const char* __restrict__ d_inputData,
    int numSegments,
    int maskingThreshold,
    char maskingLetter,
    char paddingLetter,
    InputTransformOp inputTransformOp
){
    __builtin_assume(blockDim.x == blocksize);
    
    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<32>(block);
    
    constexpr int numWarpsInBlock = blocksize / 32;
    const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
    const int warpIdInBlock = warp.meta_group_rank();
    const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;
    
    constexpr int maxUintsPerThread = maxSegmentLength / 4;

    extern __shared__ Kernel1Smem<numWarpsInBlock, maxUintsPerThread> externsmem[];

    auto& smem = externsmem[0];

    // if(threadIdx.x + blockIdx.x == 0){
    //     printf("numSegments %d\n", numSegments);
    // }

    //each warp processes 32 sequences
    for(int batch32 = warpIdInGrid; batch32 < cuda::ceil_div(numSegments, 32); batch32 += numWarpsInGrid){
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
        // __syncthreads();

        // if(batch32 == 0){
        //     for(int i = 0; i < 32; i++){
        //         printf("thread %d, i %d, segmentLengths %d\n",
        //             warp.thread_rank(), i, smem.segmentLengths[warpIdInBlock][i]);
        //     }
        // }

        const int maximumSegmentLengthInBatch = cooperative_groups::reduce(warp, workSegmentLength, cooperative_groups::greater<int>{});
        const int maximumUintsPerSegmentInBatch = cuda::ceil_div(maximumSegmentLengthInBatch, 4);

        // //need to load uintsPerThread 4byte elements for each of the 32 sequences
        const int totalElementsToLoad = maximumUintsPerSegmentInBatch * 32;

        // const int totalElementsToLoad = maxUintsPerThread * 32;
        
        for(int k = warp.thread_rank(); k < totalElementsToLoad; k += warp.num_threads()){
            // const int localLoadSegmentId = k / maxUintsPerThread;
            // const int loadLane = k % maxUintsPerThread;
            const int localLoadSegmentId = k / maximumUintsPerSegmentInBatch;
            const int loadLane = k % maximumUintsPerSegmentInBatch;
            // const int loadSegmentId = batch32 * 32 + localLoadSegmentId;
            

            const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
            const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
            const int segmentUints = cuda::ceil_div(segmentLength, 4);
            const uint32_t* const segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);

            // if(batch32 == 0){
            //     printf("k %d, localLoadSegmentId %d, loadLane %d, segmentLength %d, segmentOffset %lu, segmentUints %d\n", 
            //         k, localLoadSegmentId, loadLane, segmentLength, segmentOffset, segmentUints);
            // }

            if(loadLane < segmentUints){
                // if(loadSegmentId == 444){
                    smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] = segmentPtr[loadLane];
                // }else{
                //     smem.elements[warpIdInBlock][loadLane][localLoadSegmentId] = 0x7F7F7F7F;
                // }

                
                // if(loadSegmentId == 444){
                //     auto current4 = segmentPtr[loadLane];
                //     alignas(4) char current4chars[4];
                //     memcpy(&current4chars, &current4, 4);
                    
                //     printf("thread %d, loadLane %d, loadSegmentId %d, write %2d %2d %2d %2d to %p \n", 
                //         warp.thread_rank(), loadLane, loadSegmentId,
                //         current4chars[0], current4chars[1], current4chars[2], current4chars[3],
                //         (void*)&smem.elements[warpIdInBlock][loadLane][localLoadSegmentId]
                //     );
                // }

            }
        }
        warp.sync();



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
                                    for(int i = startOfRun; i < posInSegment; i++){
                                        const int uindex2 = i / 4;
                                        const int posui2 = i % 4;
                                        ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                                    }
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


            // int runBegin = 0;
            // for(int r = 0; r < numRuns-1; r++){
            //     const int runEnd = smem.runBeginPositions[warpIdInBlock][r+1][warp.thread_rank()];
            //     const int runLength = runEnd - runBegin;
            //     if(runLength >= maskingThreshold){
            //         // for(int i = runBegin; i < runEnd; i++){
            //         //     const int uindex2 = i / 4;
            //         //     const int posui2 = i % 4;
            //         //     ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
            //         // }
            //         ((char*)&smem.elements[warpIdInBlock][runBegin/4][warp.thread_rank()])[0] = maskingLetter;
            //     }
            //     runBegin = runEnd;
            // }
            // //handle last run
            // const int runEnd = workSegmentLength;
            // const int runLength = runEnd - runBegin;
            // if(runLength >= maskingThreshold){
            //     // for(int i = runBegin; i < runEnd; i++){
            //     //     const int uindex2 = i / 4;
            //     //     const int posui2 = i % 4;
            //     //     ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
            //     // }
            //     ((char*)&smem.elements[warpIdInBlock][runBegin/4][warp.thread_rank()])[0] = maskingLetter;
            // }

            // ((char*)&smem.elements[warpIdInBlock][numRuns/4][warp.thread_rank()])[0] = maskingLetter;

            // if(workSegmentId == 0){
            //     printf("\n");
            // }
                


            // if(workSegmentId == 444){
                for(int i = workSegmentLength; i < cuda::round_up(workSegmentLength,4); i++){
                    const int uindex2 = i / 4;
                    const int posui2 = i % 4;
                    ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = paddingLetter;
                }
            // }
            
            // for(int i = workSegmentLength; i < cuda::round_up(workSegmentLength,4); i++){
            //     const int uindex2 = i / 4;
            //     const int posui2 = i % 4;
            //     ((char*)&smem.outputelements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = paddingLetter;
            // }
            // warp.sync();

            // uint32_t* segmentOutputPtr = (uint32_t*)(d_outputData + workSegmentOffset);
            // for(int k = 0; k < cuda::ceil_div(workSegmentLength, 4); k++){
            //     //segmentOutputPtr[k] = smem.outputelements[warpIdInBlock][k][warp.thread_rank()];
            //     segmentOutputPtr[k] = smem.elements[warpIdInBlock][k][warp.thread_rank()];
            // }
        }

        //wait until masking is done in smem
        warp.sync();

        //write smem back to gmem (inverse operation of loading )
        for(int k = warp.thread_rank(); k < totalElementsToLoad; k += warp.num_threads()){
            // const int localLoadSegmentId = k / maxUintsPerThread;
            // const int loadLane = k % maxUintsPerThread;
            const int localLoadSegmentId = k / maximumUintsPerSegmentInBatch;
            const int loadLane = k % maximumUintsPerSegmentInBatch;
            // const int loadSegmentId = batch32 * 32 + localLoadSegmentId;            

            const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
            const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
            const int segmentUints = cuda::ceil_div(segmentLength, 4);
            uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);

            if(loadLane < segmentUints){
                segmentPtr[loadLane] = smem.elements[warpIdInBlock][loadLane][localLoadSegmentId];
            }
        }

        //wait until smem can be re-used
        warp.sync();

    }
}


template<int blocksize, int maxSegmentLength, class InputOffsetIterator, class InputTransformOp>
void callKernel1(
    char* __restrict__ d_outputData,
    const int* __restrict__ d_lengths,
    InputOffsetIterator d_offsets,
    const char* __restrict__ d_inputData,
    int numSegments,
    int maskingThreshold,
    char maskingLetter,
    char paddingLetter,
    InputTransformOp inputTransformOp
){
    int deviceId;
    int numSMs;
    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));

    int maxBlocksPerSM = 0;

    constexpr int maxUintsPerThread = maxSegmentLength / 4;
    constexpr int numWarpsInBlock = blocksize / 32;
    auto kernel = kernel1<blocksize, maxSegmentLength, size_t*, InputTransformOp>;
    using SmemType = Kernel1Smem<numWarpsInBlock, maxUintsPerThread>;

    int smembytes = sizeof(SmemType);

    CUDACHECK(cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes));

    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        smembytes
    ));   

    helpers::GpuTimer timer1("kernel1<" + std::to_string(blocksize) + "," + std::to_string(maxSegmentLength) + ">");
    kernel<<<cuda::ceil_div(numSegments, blocksize), blocksize, smembytes>>>(
        d_outputData,
        d_lengths,
        d_offsets,
        d_inputData,
        numSegments,
        maskingThreshold,
        maskingLetter,
        paddingLetter,
        inputTransformOp
    );
    timer1.print();
}




template<int numWarpsInBlock, int smemMaxUintsPerChunk, int maxMaskingThreshold>
struct Kernel2Smem{
    static constexpr int extraMaskingUints = cuda::ceil_div(maxMaskingThreshold, 4);
    uint32_t elements[numWarpsInBlock][smemMaxUintsPerChunk + cuda::ceil_div(maxMaskingThreshold, 4)][32+1];
    int segmentLengths[numWarpsInBlock][32];
    size_t segmentOffsets[numWarpsInBlock][32];
};

template<int blocksize, int smemMaxPositionsPerChunk, int maxMaskingThreshold, class InputOffsetIterator, class InputTransformOp>
__global__
void kernel2(
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
    
    auto block = cooperative_groups::this_thread_block();
    auto warp = cooperative_groups::tiled_partition<32>(block);
    
    constexpr int numWarpsInBlock = blocksize / 32;
    const int numWarpsInGrid = numWarpsInBlock * gridDim.x;
    const int warpIdInBlock = warp.meta_group_rank();
    const int warpIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / 32;
    
    static_assert(smemMaxPositionsPerChunk % 4 == 0);
    constexpr int smemMaxUintsPerChunk = smemMaxPositionsPerChunk / 4;
    constexpr int extraMaskingUints = cuda::ceil_div(maxMaskingThreshold, 4);

    using Smem = Kernel2Smem<numWarpsInBlock, smemMaxUintsPerChunk, maxMaskingThreshold>;
    static_assert(Smem::extraMaskingUints == extraMaskingUints);
    extern __shared__ Smem externsmem2[];
    auto& smem = externsmem2[0];

    //each warp processes 32 sequences
    for(int batch32 = warpIdInGrid; batch32 < cuda::ceil_div(numSegments, 32); batch32 += numWarpsInGrid){
        const int workSegmentId = batch32 * 32 + warp.thread_rank();
        int workSegmentLength = 0;
        size_t workSegmentOffset = 0;
        if(workSegmentId < numSegments){
            workSegmentLength = d_lengths[workSegmentId];
            workSegmentOffset = d_offsets[workSegmentId];
        }
        const int workSegmentUints = cuda::ceil_div(workSegmentLength, 4);
        const int workSegmentNumChunks = cuda::ceil_div(workSegmentUints, smemMaxUintsPerChunk);
        
        smem.segmentLengths[warpIdInBlock][warp.thread_rank()] = workSegmentLength;
        smem.segmentOffsets[warpIdInBlock][warp.thread_rank()] = workSegmentOffset;

        //wait for smem lengths and offsets
        warp.sync();

        int maximumSegmentLengthInBatch = cooperative_groups::reduce(warp, workSegmentLength, cooperative_groups::greater<int>{});

        //each segment/sequence is split into chunks of size smemMaxPositionsPerChunk
        const int numChunkIterations = cuda::ceil_div(maximumSegmentLengthInBatch, smemMaxPositionsPerChunk);

        //the current masking state
        char prev = 0;
        int count = 0;
        int startOfRunInSmem = 0;
        int numRuns = 0;
        int numMaskedRuns = 0;

        for(int chunkIter = 0; chunkIter < numChunkIterations; chunkIter++){

            //for each sequence in batch, load current chunk into smem. chunks are loaded collaboratively and stored transposed to smem
            {
                #if 0
                    const int maxPositionsToLoadPerSequenceWithoutExtra = cuda::std::min(
                        smemMaxPositionsPerChunk, 
                        maximumSegmentLengthInBatch - chunkIter * smemMaxPositionsPerChunk
                    );
                    const int maxUintsToLoadPerSequenceWithoutExtra = cuda::ceil_div(maxPositionsToLoadPerSequenceWithoutExtra, 4);
                    const int maxUintsToLoadPerSequence = extraMaskingUints + maxUintsToLoadPerSequenceWithoutExtra;
                    const int maxUintsToLoad = warp.size() * maxUintsToLoadPerSequence;
                    
                    for(int k = warp.thread_rank(); k < maxUintsToLoad; k += warp.num_threads()){
                        const int localLoadSegmentId = k / maxUintsToLoadPerSequence;
                        const int smemLoadLane = (k % maxUintsToLoadPerSequence);
                        const int globalLoadLane = chunkIter * smemMaxUintsPerChunk + smemLoadLane;
                        // const int loadSegmentId = batch32 * 32 + localLoadSegmentId;
                        

                        const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                        const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                        const int segmentUints = cuda::ceil_div(segmentLength, 4);
                        const uint32_t* segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);

                        if(globalLoadLane < segmentUints){
                            smem.elements[warpIdInBlock][smemLoadLane][localLoadSegmentId] = segmentPtr[globalLoadLane];
                        }
                    }
                #else
                    // for(int loadNr = 0; loadNr < (extraMaskingUints + smemMaxUintsPerChunk); loadNr++){
                    //     const int k = loadNr * warp.size() + warp.thread_rank();

                    //     const int localLoadSegmentId = k / (extraMaskingUints + smemMaxUintsPerChunk);
                    //     const int smemLoadLane = (k % (extraMaskingUints + smemMaxUintsPerChunk));
                    //     const int globalLoadLane = chunkIter * (smemMaxUintsPerChunk) + smemLoadLane;


                    //     const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                    //     const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                    //     const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    //     const uint32_t* segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);

                    //     if(globalLoadLane < segmentUints){
                    //         smem.elements[warpIdInBlock][smemLoadLane][localLoadSegmentId] = segmentPtr[globalLoadLane];
                    //     }
                    // }

                    uint32_t tmp[(extraMaskingUints + smemMaxUintsPerChunk)];
                    int smemLoadLane_[(extraMaskingUints + smemMaxUintsPerChunk)];
                    int localLoadSegmentId_[(extraMaskingUints + smemMaxUintsPerChunk)];
                    
                    for(int loadNr = 0; loadNr < (extraMaskingUints + smemMaxUintsPerChunk); loadNr++){
                        const int k = loadNr * warp.size() + warp.thread_rank();

                        const int localLoadSegmentId = k / (extraMaskingUints + smemMaxUintsPerChunk);
                        const int smemLoadLane = (k % (extraMaskingUints + smemMaxUintsPerChunk));
                        const int globalLoadLane = chunkIter * (smemMaxUintsPerChunk) + smemLoadLane;

                        smemLoadLane_[loadNr] = smemLoadLane;
                        localLoadSegmentId_[loadNr] = localLoadSegmentId;


                        const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                        const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                        const int segmentUints = cuda::ceil_div(segmentLength, 4);
                        const uint32_t* segmentPtr = (const uint32_t*)(d_inputData + segmentOffset);

                        if(globalLoadLane < segmentUints){
                            tmp[loadNr] = segmentPtr[globalLoadLane];
                        }else{
                            smemLoadLane_[loadNr] = -42;
                        }
                    }
                    for(int loadNr = 0; loadNr < (extraMaskingUints + smemMaxUintsPerChunk); loadNr++){
                        if(smemLoadLane_[loadNr] != -42){
                            smem.elements[warpIdInBlock][smemLoadLane_[loadNr]][localLoadSegmentId_[loadNr]] = tmp[loadNr];
                        }
                    }

                #endif
            }

            //wait until chunks are loaded to smem
            warp.sync();

            #if 1
            if(workSegmentId < numSegments && chunkIter < workSegmentNumChunks){
                const int alreadyProcessedUints = chunkIter * smemMaxUintsPerChunk;
                const int numUintsToProcess = cuda::std::min(smemMaxUintsPerChunk, workSegmentUints - alreadyProcessedUints);
                const int numPositionsToProcess = workSegmentLength - alreadyProcessedUints * 4;
                
                for(int uintindex = 0; uintindex < numUintsToProcess; uintindex++){
                    uint32_t current4 = smem.elements[warpIdInBlock][uintindex][warp.thread_rank()];
                    alignas(4) char current4chars[4];
                    memcpy(&current4chars, &current4, 4);

                    if(chunkIter == 0 && uintindex == 0){
                        prev = current4chars[0];
                        numRuns++;
                    }
                    for(int posIn4 = 0; posIn4 < 4; posIn4++){
                        const int posInSegmentSmem = uintindex * 4 + posIn4;
                        const int posInSegment = posInSegmentSmem + chunkIter * smemMaxPositionsPerChunk;
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

                // if we processed a chunk which is not the last chunk, 
                // process the extra masking uints to handle overlap with next chunk
                // if(alreadyProcessedUints + numUintsToProcess < workSegmentUints){
                if(chunkIter < workSegmentNumChunks-1){
                    // if(workSegmentId == 842){
                    //     printf("chunkIter %d, prev %d, count %d, startOfRunInSmem %d\n",
                    //         chunkIter, prev, count, startOfRunInSmem
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
                    //for next chunk, run begins at begin of smem chunk
                    startOfRunInSmem = 0;
                    count = count_backup;
                }else{
                    //we processed the last chunk of sequence.

                    //handle last run
                    if(count >= maskingThreshold){
                        
                        for(int i = startOfRunInSmem; i < numPositionsToProcess; i++){
                            const int uindex2 = i / 4;
                            const int posui2 = i % 4;
                            ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = maskingLetter;
                        }
                        numMaskedRuns++;
                    }
                    if(applyExtraPaddingLetter){
                        for(int i = numPositionsToProcess; i < cuda::round_up(numPositionsToProcess,4); i++){
                            const int uindex2 = i / 4;
                            const int posui2 = i % 4;
                            // assert(uindex2 < (extraMaskingUints + smemMaxUintsPerChunk));
                            // assert(i < smemMaxPositionsPerChunk);
                            // printf("uindex2 %d, posui2 %d, numPositionsToProcess %d %d\n", 
                            //     uindex2, posui2, numPositionsToProcess, cuda::round_up(numPositionsToProcess,4));
                            ((char*)&smem.elements[warpIdInBlock][uindex2][warp.thread_rank()])[posui2] = paddingLetter;
                        }
                    }
                    // //debug
                    // ((char*)&smem.elements[warpIdInBlock][numMaskedRuns / 4][warp.thread_rank()])[0] = paddingLetter;
                }
            }
            
            //wait until chunks are processed in smem
            warp.sync();
            #endif

            //write smem chunks back to gmem
            {
                const int maxPositionsToLoadPerSequenceWithoutExtra = cuda::std::min(
                    smemMaxPositionsPerChunk, 
                    maximumSegmentLengthInBatch - chunkIter * smemMaxPositionsPerChunk
                );
                const int maxUintsToLoadPerSequenceWithoutExtra = cuda::ceil_div(maxPositionsToLoadPerSequenceWithoutExtra, 4);
                const int maxUintsToLoadPerSequence = extraMaskingUints + maxUintsToLoadPerSequenceWithoutExtra;
                const int maxUintsToLoad = warp.size() * maxUintsToLoadPerSequence;
                
                for(int k = warp.thread_rank(); k < maxUintsToLoad; k += warp.num_threads()){
                    const int localLoadSegmentId = k / maxUintsToLoadPerSequence;
                    const int smemLoadLane = (k % maxUintsToLoadPerSequence);
                    const int globalLoadLane = chunkIter * smemMaxUintsPerChunk + smemLoadLane;
                    // const int loadSegmentId = batch32 * 32 + localLoadSegmentId;
                    

                    const int segmentLength = smem.segmentLengths[warpIdInBlock][localLoadSegmentId];
                    const size_t segmentOffset = smem.segmentOffsets[warpIdInBlock][localLoadSegmentId];
                    const int segmentUints = cuda::ceil_div(segmentLength, 4);
                    uint32_t* const segmentPtr = (uint32_t*)(d_outputData + segmentOffset);

                    if(globalLoadLane < segmentUints){
                        segmentPtr[globalLoadLane] = smem.elements[warpIdInBlock][smemLoadLane][localLoadSegmentId];
                    }
                }

            }
            

            //wait until smem can be re-used for next chunk
            warp.sync();
        }
    }


}

template<int blocksize, int smemMaxPositionsPerChunk, int maxMaskingThreshold, class InputOffsetIterator, class InputTransformOp>
void callKernel2(
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
    if(maskingThreshold > maxMaskingThreshold){
        throw std::runtime_error("maskingThreshold > maxMaskingThreshold");
    }
    int deviceId;
    int numSMs;
    CUDACHECK(cudaGetDevice(&deviceId));
    CUDACHECK(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId));

    int maxBlocksPerSM = 0;

    constexpr int smemMaxUintsPerChunk = smemMaxPositionsPerChunk / 4;
    constexpr int numWarpsInBlock = blocksize / 32;
    auto kernel = kernel2<blocksize, smemMaxPositionsPerChunk, maxMaskingThreshold, size_t*, InputTransformOp>;
    using SmemType = Kernel2Smem<numWarpsInBlock, smemMaxUintsPerChunk, maxMaskingThreshold>;

    int smembytes = sizeof(SmemType);

    CUDACHECK(cudaFuncSetAttribute(reinterpret_cast<const void*>(kernel), cudaFuncAttributeMaxDynamicSharedMemorySize, smembytes));

    CUDACHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &maxBlocksPerSM,
        kernel,
        blocksize, 
        smembytes
    ));   

    helpers::GpuTimer timer1("kernel2<" + std::to_string(blocksize) + "," + std::to_string(smemMaxPositionsPerChunk) + ">");
    kernel<<<cuda::ceil_div(numSegments, blocksize), blocksize, smembytes>>>(
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
    timer1.print();
}

}






LIBMARV_NAMESPACE_WITH_NESTING_END

#endif