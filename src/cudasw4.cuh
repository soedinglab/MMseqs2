#ifndef CUDASW4_CUH
#define CUDASW4_CUH

#include "hpc_helpers/cuda_raiiwrappers.cuh"
#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/nvtx_markers.cuh"
#include "hpc_helpers/simple_allocation.cuh"

#include "config.hpp"
#include "dbdata.hpp"
#include "length_partitions.hpp"
#include "util.cuh"
#include "kernels.cuh"
#include "blosum.hpp"
#include "types.hpp"
#include "dbbatching.cuh"
#include "convert.cuh"
#include "target_subject_ids.cuh"
#include "gapless_kernel_config.cuh"
#include "smithwaterman_kernel_config.cuh"
#include "types.hpp"

#include "pssm.cuh"
#include "pssmkernels_gapless.cuh"
#include "pssmkernels_smithwaterman.cuh"

#include "gpudatabaseallocation.cuh"

#include <thrust/binary_search.h>
#include <thrust/sort.h>
#include <thrust/equal.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/merge.h>

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <string_view>
#include <optional>
#include <sstream>

namespace cudasw4{

    template<class T, T factor>
    struct RoundToNextMultiple{
        __host__ __device__ 
        T operator()(const T& value){
            return SDIV(value, factor) * factor;
        }        
    };

    struct CompareScoresDescendingRefIdsAscending{
        template<class Tuple1, class Tuple2>
        __host__ __device__
        bool operator()(const Tuple1& lhs, const Tuple2& rhs) const{
            const auto scoreL = thrust::get<0>(lhs);
            const auto refIdL = thrust::get<1>(lhs);
            const auto scoreR = thrust::get<0>(rhs);
            const auto refIdR = thrust::get<1>(rhs);
            if(scoreL < scoreR) return false;
            if(scoreL > scoreR) return true;
            //scores are equal
            return refIdL < refIdR;
        }
    };

    __global__
    void addKernel(int* output, const int* input1, const int* input2){
        *output = *input1 + *input2;
    }


    __global__
    void sumNumOverflowsKernel(int* output, const int* input, int numGpus){
        int sum = input[0];
        for(int gpu = 1; gpu < numGpus; gpu++){
            sum += input[gpu];
        }
        output[0] = sum;
    }

    template<class PartitionOffsets, class Indices>
    __global__
    void transformLocalSequenceIndicesToGlobalIndices(
        int gpu,
        int N,
        PartitionOffsets partitionOffsets,
        Indices maxReduceArrayIndices
    ){
        const int tid = threadIdx.x + blockIdx.x * blockDim.x;

        if(tid < N){
            maxReduceArrayIndices[tid] = partitionOffsets.getGlobalIndex(gpu, maxReduceArrayIndices[tid]);
        }
    }

    template<bool isEncoded_>
    struct QueryView{
        static constexpr bool isEncoded = isEncoded_;

        QueryView(const char* ptr_, SequenceLengthT length_) : ptr(ptr_), length(length_){
            if(ptr == nullptr) throw std::runtime_error("QueryView constructed from nullptr");
        }
        QueryView(const QueryView&) = default;
        QueryView& operator=(const QueryView&) = default;

        const char* ptr{};
        SequenceLengthT length{};
    };
    using DecodedQueryView = QueryView<false>;
    using EncodedQueryView = QueryView<true>;

    struct KernelConfigFilenames{
        std::optional<std::string> gapless;
        std::optional<std::string> sw;
    };

    struct BenchmarkStats{
        int numOverflows{};
        double seconds{};
        double gcups{};
    };

    struct ScanResult{
        std::vector<int> scores{};
        std::vector<ReferenceIdT> referenceIds{};
        std::vector<AlignmentEndPosition> endPositions{};
        BenchmarkStats stats{};
    };

    struct MemoryConfig{
        size_t maxBatchBytes = 128ull * 1024ull * 1024ull;
        size_t maxBatchSequences = 10'000'000;
        size_t maxTempBytes = 4ull * 1024ull * 1024ull * 1024ull;
        size_t maxGpuMem = std::numeric_limits<size_t>::max();
    };

    struct HostGpuPartitionOffsets{
        int numGpus;
        int numLengthPartitions;
        std::vector<size_t> partitionSizes;
        std::vector<size_t> horizontalPS;
        std::vector<size_t> verticalPS;
        std::vector<size_t> totalPerLengthPartitionPS;

        HostGpuPartitionOffsets() = default;

        HostGpuPartitionOffsets(int numGpus_, int numLengthpartitions_, std::vector<size_t> partitionSizes_)
            : numGpus(numGpus_), 
            numLengthPartitions(numLengthpartitions_), 
            partitionSizes(std::move(partitionSizes_)),
            horizontalPS(numGpus * numLengthPartitions, 0),
            verticalPS(numGpus * numLengthPartitions, 0),
            totalPerLengthPartitionPS(numLengthPartitions, 0)
        {
            assert(partitionSizes.size() == numGpus * numLengthPartitions);

            for(int gpu = 0; gpu < numGpus; gpu++){
                for(int l = 1; l < numLengthPartitions; l++){
                    horizontalPS[gpu * numLengthPartitions + l] = horizontalPS[gpu * numLengthPartitions + l-1] + partitionSizes[gpu * numLengthPartitions + l-1];
                }
            }
            for(int l = 0; l < numLengthPartitions; l++){
                for(int gpu = 1; gpu < numGpus; gpu++){
                    verticalPS[gpu * numLengthPartitions + l] = verticalPS[(gpu-1) * numLengthPartitions + l] + partitionSizes[(gpu-1) * numLengthPartitions + l];
                }
            }
            for(int l = 1; l < numLengthPartitions; l++){
                totalPerLengthPartitionPS[l] = totalPerLengthPartitionPS[l-1] 
                    + (verticalPS[(numGpus - 1) * numLengthPartitions + (l-1)] + partitionSizes[(numGpus-1) * numLengthPartitions + (l-1)]);
            }
        }

        size_t getGlobalIndex(int gpu, size_t localIndex) const {
            const size_t* const myHorizontalPS = horizontalPS.data() + gpu * numLengthPartitions;
            const auto it = std::lower_bound(myHorizontalPS, myHorizontalPS + numLengthPartitions, localIndex+1);
            const int whichPartition = std::distance(myHorizontalPS, it) - 1;
            const size_t occurenceInPartition = localIndex - myHorizontalPS[whichPartition];
            const size_t globalPartitionBegin = totalPerLengthPartitionPS[whichPartition];
            const size_t elementsOfOtherPreviousGpusInPartition = verticalPS[gpu * numLengthPartitions + whichPartition];
            //std::cout << "whichPartition " << whichPartition << ", occurenceInPartition " << occurenceInPartition 
            //    << ", globalPartitionBegin " << globalPartitionBegin << ", elementsOfOtherPreviousGpusInPartition " << elementsOfOtherPreviousGpusInPartition << "\n";
            return globalPartitionBegin + elementsOfOtherPreviousGpusInPartition + occurenceInPartition;
        };

        void print(std::ostream& os){
            os << "numGpus " << numGpus << "\n";
            os << "numLengthPartitions " << numLengthPartitions << "\n";
            os << "partitionSizes\n";
            for(int gpu = 0; gpu < numGpus; gpu++){
                for(int l = 0; l < numLengthPartitions; l++){
                    os << partitionSizes[gpu * numLengthPartitions + l] << " ";
                }
                os << "\n";
            }
            os << "\n";
            os << "horizontalPS\n";
            for(int gpu = 0; gpu < numGpus; gpu++){
                for(int l = 0; l < numLengthPartitions; l++){
                    os << horizontalPS[gpu * numLengthPartitions + l] << " ";
                }
                os << "\n";
            }
            os << "\n";
            os << "verticalPS\n";
            for(int gpu = 0; gpu < numGpus; gpu++){
                for(int l = 0; l < numLengthPartitions; l++){
                    os << verticalPS[gpu * numLengthPartitions + l] << " ";
                }
                os << "\n";
            }
            os << "\n";
            os << "totalPerLengthPartitionPS\n";
            for(int l = 0; l < numLengthPartitions; l++){
                os << totalPerLengthPartitionPS[l] << " ";
            }
            os << "\n";
        }
    };

    struct DeviceGpuPartitionOffsets{
        template<class T>
        using MyDeviceBuffer = helpers::SimpleAllocationDevice<T, 0>;

        int numGpus;
        int numLengthPartitions;
        MyDeviceBuffer<size_t> partitionSizes;
        MyDeviceBuffer<size_t> horizontalPS;
        MyDeviceBuffer<size_t> verticalPS;
        MyDeviceBuffer<size_t> totalPerLengthPartitionPS;

        struct View{
            int numGpus;
            int numLengthPartitions;
            const size_t* partitionSizes;
            const size_t* horizontalPS;
            const size_t* verticalPS;
            const size_t* totalPerLengthPartitionPS;

            __device__
            size_t getGlobalIndex(int gpu, size_t localIndex) const {
                const size_t* const myHorizontalPS = horizontalPS + gpu * numLengthPartitions;
                const auto it = thrust::lower_bound(thrust::seq, myHorizontalPS, myHorizontalPS + numLengthPartitions, localIndex+1);
                const int whichPartition = thrust::distance(myHorizontalPS, it) - 1;
                const size_t occurenceInPartition = localIndex - myHorizontalPS[whichPartition];
                const size_t globalPartitionBegin = totalPerLengthPartitionPS[whichPartition];
                const size_t elementsOfOtherPreviousGpusInPartition = verticalPS[gpu * numLengthPartitions + whichPartition];
                return globalPartitionBegin + elementsOfOtherPreviousGpusInPartition + occurenceInPartition;
            };
        };

        DeviceGpuPartitionOffsets() = default;
        DeviceGpuPartitionOffsets(const HostGpuPartitionOffsets& hostData)
            : numGpus(hostData.numGpus),
            numLengthPartitions(hostData.numLengthPartitions),
            partitionSizes(numGpus * numLengthPartitions),
            horizontalPS(numGpus * numLengthPartitions),
            verticalPS(numGpus * numLengthPartitions),
            totalPerLengthPartitionPS(numLengthPartitions)
        {
            cudaMemcpyAsync(partitionSizes.data(), hostData.partitionSizes.data(), sizeof(size_t) * numGpus * numLengthPartitions, H2D, cudaStreamLegacy); CUERR;
            cudaMemcpyAsync(horizontalPS.data(), hostData.horizontalPS.data(), sizeof(size_t) * numGpus * numLengthPartitions, H2D, cudaStreamLegacy); CUERR;
            cudaMemcpyAsync(verticalPS.data(), hostData.verticalPS.data(), sizeof(size_t) * numGpus * numLengthPartitions, H2D, cudaStreamLegacy); CUERR;
            cudaMemcpyAsync(totalPerLengthPartitionPS.data(), hostData.totalPerLengthPartitionPS.data(), sizeof(size_t) * numLengthPartitions, H2D, cudaStreamLegacy); CUERR;
        }

        View getDeviceView() const{
            View view;
            view.numGpus = numGpus;
            view.numLengthPartitions = numLengthPartitions;
            view.partitionSizes = partitionSizes.data();
            view.horizontalPS = horizontalPS.data();
            view.verticalPS = verticalPS.data();
            view.totalPerLengthPartitionPS = totalPerLengthPartitionPS.data();
            return view;
        }
    };


    class CudaSW4{
    public:
        template<class T>
        using MyPinnedBuffer = helpers::SimpleAllocationPinnedHost<T, 0>;
        template<class T>
        using MyDeviceBuffer = helpers::SimpleAllocationDevice<T, 0>;

        struct GpuWorkingSet{

            //using MaxReduceArray = TopNMaximaArray<maxReduceArraySize>;
            using MaxReduceArray = TopNMaximaArray;
            using MaxReduceArrayWithEndPositions = TopNMaximaArrayWithExtra<AlignmentEndPosition>;

            GpuWorkingSet(
                size_t gpumemlimit,
                size_t maxBatchBytes,
                size_t maxBatchSequences,
                size_t maxTempBytes,
                const std::vector<DBdataView>& dbPartitions,
                const std::vector<DeviceBatchCopyToPinnedPlan>& dbBatches,
                bool needsPinnedStagingBuffers,
                int maxReduceArraySize_ = 512 * 1024
            ) : maxReduceArraySize(maxReduceArraySize_)
            {
                cudaGetDevice(&deviceId);

                size_t numSubjects = 0;
                size_t numSubjectBytes = 0;
                for(const auto& p : dbPartitions){
                    numSubjects += p.numSequences();
                    numSubjectBytes += p.numChars();
                }
        
                d_query.resize(1024*1024); CUERR
                gpuFullQueryPSSM.resize(10000, 21);

                numTempBytes = std::min(maxTempBytes, gpumemlimit);
                d_tempStorageHE.resize(numTempBytes);
                d_maxReduceArrayScores.resize(maxReduceArraySize);
                d_maxReduceArrayIndices.resize(maxReduceArraySize);
                d_maxReduceArrayExtras.resize(maxReduceArraySize);                

                size_t usedGpuMem = 0;
                usedGpuMem += numTempBytes;
                usedGpuMem += sizeof(float) * maxReduceArraySize; // d_maxReduceArrayScores
                usedGpuMem += sizeof(ReferenceIdT) * maxReduceArraySize; // d_maxReduceArrayIndices
                usedGpuMem += sizeof(AlignmentEndPosition) * maxReduceArraySize; //d_maxReduceArrayExtras

                if(usedGpuMem > gpumemlimit){
                    throw std::runtime_error("Out of memory working set");
                }
                
        
                //devAlignmentScoresFloat.resize(numSubjects);
        
                forkStreamEvent = CudaEvent{cudaEventDisableTiming}; CUERR;
                numWorkStreamsWithoutTemp = 1;
                workstreamIndex = 0;
                workStreamsWithoutTemp.resize(numWorkStreamsWithoutTemp);

                numCopyBuffers = 2;
        
                h_chardata_vec.resize(numCopyBuffers);
                h_lengthdata_vec.resize(numCopyBuffers);
                h_offsetdata_vec.resize(numCopyBuffers);
                d_chardata_vec.resize(numCopyBuffers);
                d_lengthdata_vec.resize(numCopyBuffers);
                d_offsetdata_vec.resize(numCopyBuffers);
                copyStreams.resize(numCopyBuffers);
                pinnedBufferEvents.resize(numCopyBuffers);
                deviceBufferEvents.resize(numCopyBuffers);
                //d_total_overflow_number.resize(1);
                //d_overflow_number.resize(numCopyBuffers);
                //h_overflow_number.resize(numCopyBuffers);
                //d_overflow_positions_vec.resize(numCopyBuffers);

                size_t memoryRequiredForFullDB = 0;
                memoryRequiredForFullDB += numSubjectBytes; // d_fulldb_chardata
                memoryRequiredForFullDB += sizeof(SequenceLengthT) * numSubjects; //d_fulldb_lengthdata
                memoryRequiredForFullDB += sizeof(size_t) * (numSubjects+1); //d_fulldb_offsetdata
                //memoryRequiredForFullDB += sizeof(ReferenceIdT) * numSubjects * 2; //d_overflow_positions_vec
        
               

                if(usedGpuMem + memoryRequiredForFullDB <= gpumemlimit){
                    numBatchesInCachedDB = dbBatches.size();
                    charsOfBatches = numSubjectBytes;
                    subjectsOfBatches = numSubjects;
                    d_cacheddb = std::make_shared<GpuDatabaseAllocation>(numSubjectBytes, numSubjects);

                    for(int i = 0; i < numCopyBuffers; i++){
                        if(needsPinnedStagingBuffers){
                            h_chardata_vec[i].resize(maxBatchBytes);
                            h_lengthdata_vec[i].resize(maxBatchSequences);
                            h_offsetdata_vec[i].resize(maxBatchSequences+1);
                        }
                        pinnedBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        deviceBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        //d_overflow_positions_vec[i].resize(numSubjects);
                    }
                }else{
                    //allocate a double buffer for batch transfering
                    size_t memoryRequiredForBatchedProcessing = 0;
                    memoryRequiredForBatchedProcessing += maxBatchBytes * numCopyBuffers; // d_chardata_vec
                    memoryRequiredForBatchedProcessing += sizeof(SequenceLengthT) * maxBatchSequences * numCopyBuffers; //d_lengthdata_vec
                    memoryRequiredForBatchedProcessing += sizeof(size_t) * (maxBatchSequences+1) * numCopyBuffers; //d_offsetdata_vec
                    usedGpuMem += memoryRequiredForBatchedProcessing;
                    if(usedGpuMem > gpumemlimit){
                        throw std::runtime_error("Out of memory working set");
                    }
                    
                    for(int i = 0; i < numCopyBuffers; i++){
                        if(needsPinnedStagingBuffers){
                            h_chardata_vec[i].resize(maxBatchBytes);
                            h_lengthdata_vec[i].resize(maxBatchSequences);
                            h_offsetdata_vec[i].resize(maxBatchSequences+1);
                        }
                        d_chardata_vec[i].resize(maxBatchBytes);
                        d_lengthdata_vec[i].resize(maxBatchSequences);
                        d_offsetdata_vec[i].resize(maxBatchSequences+1);
                        pinnedBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        deviceBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        //d_overflow_positions_vec[i].resize(numSubjects);
                    }

                    //count how many batches fit into remaining gpu memory

                    numBatchesInCachedDB = 0;
                    charsOfBatches = 0;
                    subjectsOfBatches = 0;
                    size_t totalRequiredMemForBatches = sizeof(size_t);
                    for(; numBatchesInCachedDB < dbBatches.size(); numBatchesInCachedDB++){
                        const auto& batch = dbBatches[numBatchesInCachedDB];
                        const size_t requiredMemForBatch = batch.usedSeq * sizeof(SequenceLengthT) + batch.usedSeq * sizeof(size_t) + batch.usedBytes;
                        if(usedGpuMem + totalRequiredMemForBatches + requiredMemForBatch <= gpumemlimit){
                            //ok, fits
                            totalRequiredMemForBatches += requiredMemForBatch;
                            charsOfBatches += batch.usedBytes;
                            subjectsOfBatches += batch.usedSeq;
                        }else{
                            //does not fit
                            break;
                        }
                    }
                    assert(numBatchesInCachedDB < dbBatches.size());

                    //std::cout << "numBatchesInCachedDB " << numBatchesInCachedDB << ", charsOfBatches " << charsOfBatches << ", subjectsOfBatches " << subjectsOfBatches << "\n";

                    assert(usedGpuMem + totalRequiredMemForBatches <= gpumemlimit);
                    d_cacheddb = std::make_shared<GpuDatabaseAllocation>(charsOfBatches, subjectsOfBatches);
                }
            }

            GpuWorkingSet(
                size_t gpumemlimit,
                size_t maxBatchBytes,
                size_t maxBatchSequences,
                size_t maxTempBytes,
                const std::vector<DBdataView>& dbPartitions,
                const std::vector<DeviceBatchCopyToPinnedPlan>& dbBatches,
                std::shared_ptr<GpuDatabaseAllocationBase> existingGpuDBAllocation,
                bool needsPinnedStagingBuffers,
                int maxReduceArraySize_ = 512 * 1024
            ) : maxReduceArraySize(maxReduceArraySize_)
            {
                cudaGetDevice(&deviceId);

                assert(existingGpuDBAllocation != nullptr);

                size_t numSubjects = 0;
                size_t numSubjectBytes = 0;
                for(const auto& p : dbPartitions){
                    numSubjects += p.numSequences();
                    numSubjectBytes += p.numChars();
                }
        
                d_query.resize(1024*1024); CUERR
                gpuFullQueryPSSM.resize(10000, 21);

                numTempBytes = std::min(maxTempBytes, gpumemlimit);
                d_tempStorageHE.resize(numTempBytes);
                d_maxReduceArrayScores.resize(maxReduceArraySize);
                d_maxReduceArrayIndices.resize(maxReduceArraySize);
                d_maxReduceArrayExtras.resize(maxReduceArraySize);

                size_t usedGpuMem = 0;
                usedGpuMem += numTempBytes;
                usedGpuMem += sizeof(float) * maxReduceArraySize; // d_maxReduceArrayScores
                usedGpuMem += sizeof(ReferenceIdT) * maxReduceArraySize; // d_maxReduceArrayIndices
                usedGpuMem += sizeof(AlignmentEndPosition) * maxReduceArraySize; //d_maxReduceArrayExtras

                if(usedGpuMem > gpumemlimit){
                    throw std::runtime_error("Out of memory working set");
                }
                
        
                //devAlignmentScoresFloat.resize(numSubjects);
        
                forkStreamEvent = CudaEvent{cudaEventDisableTiming}; CUERR;
                numWorkStreamsWithoutTemp = 10;
                workstreamIndex = 0;
                workStreamsWithoutTemp.resize(numWorkStreamsWithoutTemp);

                numCopyBuffers = 2;
        
                h_chardata_vec.resize(numCopyBuffers);
                h_lengthdata_vec.resize(numCopyBuffers);
                h_offsetdata_vec.resize(numCopyBuffers);
                d_chardata_vec.resize(numCopyBuffers);
                d_lengthdata_vec.resize(numCopyBuffers);
                d_offsetdata_vec.resize(numCopyBuffers);
                copyStreams.resize(numCopyBuffers);
                pinnedBufferEvents.resize(numCopyBuffers);
                deviceBufferEvents.resize(numCopyBuffers);
                //d_total_overflow_number.resize(1);
                //d_overflow_number.resize(numCopyBuffers);
                //h_overflow_number.resize(numCopyBuffers);
                //d_overflow_positions_vec.resize(numCopyBuffers);
        
                d_cacheddb = existingGpuDBAllocation;

                if(d_cacheddb->getNumChars() >= numSubjectBytes && d_cacheddb->getNumSubjects() >= numSubjects){
                    numBatchesInCachedDB = dbBatches.size();
                    charsOfBatches = numSubjectBytes;
                    subjectsOfBatches = numSubjects;

                    for(int i = 0; i < numCopyBuffers; i++){
                        if(needsPinnedStagingBuffers){
                            h_chardata_vec[i].resize(maxBatchBytes);
                            h_lengthdata_vec[i].resize(maxBatchSequences);
                            h_offsetdata_vec[i].resize(maxBatchSequences+1);
                        }
                        pinnedBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        deviceBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        //d_overflow_positions_vec[i].resize(numSubjects);
                    }
                }else{
                    //allocate a double buffer for batch transfering
                    size_t memoryRequiredForBatchedProcessing = 0;
                    memoryRequiredForBatchedProcessing += maxBatchBytes * numCopyBuffers; // d_chardata_vec
                    memoryRequiredForBatchedProcessing += sizeof(SequenceLengthT) * maxBatchSequences * numCopyBuffers; //d_lengthdata_vec
                    memoryRequiredForBatchedProcessing += sizeof(size_t) * (maxBatchSequences+1) * numCopyBuffers; //d_offsetdata_vec
                    usedGpuMem += memoryRequiredForBatchedProcessing;

                    //std::cout << "usedGpuMem " << usedGpuMem << ", gpumemlimit " << gpumemlimit << "\n";

                    //cached db is already accounted for because gpumemlimit was obtained after cached db was allocated

                    // if(usedGpuMem > gpumemlimit){
                    //     throw std::runtime_error("Out of memory working set");
                    // }
                    
                    for(int i = 0; i < numCopyBuffers; i++){
                        if(needsPinnedStagingBuffers){
                            h_chardata_vec[i].resize(maxBatchBytes);
                            h_lengthdata_vec[i].resize(maxBatchSequences);
                            h_offsetdata_vec[i].resize(maxBatchSequences+1);
                        }
                        d_chardata_vec[i].resize(maxBatchBytes);
                        d_lengthdata_vec[i].resize(maxBatchSequences);
                        d_offsetdata_vec[i].resize(maxBatchSequences+1);
                        pinnedBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        deviceBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        //d_overflow_positions_vec[i].resize(numSubjects);
                    }

                    //count how many batches fit into d_cacheddb

                    numBatchesInCachedDB = 0;
                    charsOfBatches = 0;
                    subjectsOfBatches = 0;
                    for(; numBatchesInCachedDB < dbBatches.size(); numBatchesInCachedDB++){
                        const auto& batch = dbBatches[numBatchesInCachedDB];
                        if(subjectsOfBatches + batch.usedSeq <= d_cacheddb->getNumSubjects() && charsOfBatches + batch.usedBytes <= d_cacheddb->getNumChars()){
                            //ok, fits
                            charsOfBatches += batch.usedBytes;
                            subjectsOfBatches += batch.usedSeq;
                        }else{
                            //does not fit
                            break;
                        }
                    }
                    assert(charsOfBatches <= d_cacheddb->getNumChars());
                    assert(subjectsOfBatches <= d_cacheddb->getNumSubjects());
                    assert(numBatchesInCachedDB < dbBatches.size());


                }
            }
        
            MaxReduceArray getMaxReduceArray(size_t offset){
                return MaxReduceArray(
                    d_maxReduceArrayScores.data(), 
                    d_maxReduceArrayIndices.data(), 
                    offset,
                    maxReduceArraySize
                );
            }
        
            void resetMaxReduceArray(cudaStream_t stream){
                thrust::fill(thrust::cuda::par_nosync.on(stream),
                    d_maxReduceArrayScores.data(),
                    d_maxReduceArrayScores.data() + maxReduceArraySize,
                    -1.f
                );
                cudaMemsetAsync(d_maxReduceArrayIndices.data(), 0, sizeof(ReferenceIdT) * maxReduceArraySize, stream); CUERR;
            }

            MaxReduceArrayWithEndPositions getMaxReduceArrayWithEndPositions(size_t offset){
                return MaxReduceArrayWithEndPositions(
                    d_maxReduceArrayScores.data(), 
                    d_maxReduceArrayIndices.data(), 
                    d_maxReduceArrayExtras.data(),
                    offset,
                    maxReduceArraySize
                );
            }
        
            void resetMaxReduceArrayWithEndPositions(cudaStream_t stream){
                resetMaxReduceArray(stream);
                cudaMemsetAsync(d_maxReduceArrayExtras.data(), 0, sizeof(AlignmentEndPosition) * maxReduceArraySize, stream); CUERR;
            }

            void resetTopNArrays(cudaStream_t stream){
                thrust::fill(thrust::cuda::par_nosync.on(stream),
                    d_topN_scores.data(),
                    d_topN_scores.data() + d_topN_scores.size(),
                    -1.f
                );
                cudaMemsetAsync(d_topN_refIds.data(), 0, sizeof(ReferenceIdT) * d_topN_refIds.size(), stream); CUERR;
                cudaMemsetAsync(d_topN_alignmentEndPositions.data(), 0, sizeof(AlignmentEndPosition) * d_topN_alignmentEndPositions.size(), stream); CUERR;
            }
        
            void setPartitionOffsets(const HostGpuPartitionOffsets& offsets){
                deviceGpuPartitionOffsets = DeviceGpuPartitionOffsets(offsets);
            }
            
            size_t getNumCharsInCachedDB() const{
                return charsOfBatches;
            }

            size_t getNumSequencesInCachedDB() const{
                return subjectsOfBatches;
            }

            size_t getNumBatchesInCachedDB() const{
                return numBatchesInCachedDB;
            }

            void setTopNSize(size_t topN){
                d_topN_scores.resize(2*topN);
                d_topN_refIds.resize(2*topN);
                d_topN_alignmentEndPositions.resize(2*topN);
                d_topN_scores_tmp.resize(2*topN);
                d_topN_refIds_tmp.resize(2*topN);
                d_topN_alignmentEndPositions_tmp.resize(2*topN);
            }
        

            int deviceId;
            int numCopyBuffers;
            int numWorkStreamsWithoutTemp = 1;
            int workstreamIndex;
            int copyBufferIndex = 0;
            int maxReduceArraySize = 512 * 1024;
            size_t numTempBytes;
            size_t numBatchesInCachedDB = 0;
            size_t charsOfBatches = 0;
            size_t subjectsOfBatches = 0;
        
            MyDeviceBuffer<float> d_maxReduceArrayScores;
            MyDeviceBuffer<ReferenceIdT> d_maxReduceArrayIndices;
            MyDeviceBuffer<AlignmentEndPosition> d_maxReduceArrayExtras;

            MyDeviceBuffer<float> d_topN_scores;
            MyDeviceBuffer<ReferenceIdT> d_topN_refIds;
            MyDeviceBuffer<float> d_topN_scores_tmp;
            MyDeviceBuffer<ReferenceIdT> d_topN_refIds_tmp;
            MyDeviceBuffer<AlignmentEndPosition> d_topN_alignmentEndPositions;
            MyDeviceBuffer<AlignmentEndPosition> d_topN_alignmentEndPositions_tmp;
        
            MyDeviceBuffer<char> d_query;
            MyDeviceBuffer<char> d_tempStorageHE;
            //MyDeviceBuffer<float> devAlignmentScoresFloat;
            // MyDeviceBuffer<size_t> d_selectedPositions;
            //MyDeviceBuffer<int> d_total_overflow_number;
            //MyDeviceBuffer<int> d_overflow_number;
            //MyPinnedBuffer<int> h_overflow_number;
            CudaStream hostFuncStream;
            CudaStream workStreamForTempUsage;
            CudaEvent forkStreamEvent;

            size_t maxNumBatchesInCachedDB = 0;
            std::shared_ptr<GpuDatabaseAllocationBase> d_cacheddb;

            
            std::vector<MyPinnedBuffer<char>> h_chardata_vec;
            std::vector<MyPinnedBuffer<SequenceLengthT>> h_lengthdata_vec;
            std::vector<MyPinnedBuffer<size_t>> h_offsetdata_vec;
            std::vector<MyDeviceBuffer<char>> d_chardata_vec;
            std::vector<MyDeviceBuffer<SequenceLengthT>> d_lengthdata_vec;
            std::vector<MyDeviceBuffer<size_t>> d_offsetdata_vec;
            std::vector<CudaStream> copyStreams;
            std::vector<CudaEvent> pinnedBufferEvents;
            std::vector<CudaEvent> deviceBufferEvents;
            std::vector<CudaStream> workStreamsWithoutTemp;
            //std::vector<MyDeviceBuffer<ReferenceIdT>> d_overflow_positions_vec;
        
            DeviceGpuPartitionOffsets deviceGpuPartitionOffsets;

            GpuPSSM gpuFullQueryPSSM;
            GpuPermutedPSSMforGapless gpuPermutedPSSMforGapless;
            GpuPermutedPSSMforSW gpuPermutedPSSMforSW;
        };

        struct SequenceLengthStatistics{
            SequenceLengthT max_length = 0;
            SequenceLengthT min_length = std::numeric_limits<SequenceLengthT>::max();
            size_t sumOfLengths = 0;
        };
    private:
        struct BatchDstInfo{
            bool isUploaded{};
            char* charsPtr{};
            SequenceLengthT* lengthsPtr{};
            size_t* offsetsPtr{};
        };

    public:

        CudaSW4(
            std::vector<int> deviceIds_, 
            int numTop,
            BlosumType blosumType,
            const MemoryConfig& memoryConfig,
            bool verbose_,
            const KernelConfigFilenames& kernelConfigFilenames
        ) : deviceIds(std::move(deviceIds_)), verbose(verbose_)
        {
            #ifdef CUDASW_DEBUG_CHECK_CORRECTNESS
                blosumType = BlosumType::BLOSUM62_20;
            #endif
            if(deviceIds.size() == 0){ 
                throw std::runtime_error("No device selected");
            
            }
            RevertDeviceId rdi{};

            initializeGpus();

            //resultNumOverflows.resize(1);

            const int numGpus = deviceIds.size();
            cudaSetDevice(deviceIds[0]);
            
            //d_resultNumOverflows.resize(numGpus);
            scanTimer = std::make_unique<helpers::GpuTimer>("Scan");
            totalTimer = std::make_unique<helpers::GpuTimer>("Total");

            setBlosum(blosumType);
            setNumTop(numTop);
            setMemoryConfig(memoryConfig);

            initializeListOfAvailableKernelConfigs(kernelConfigFilenames);

            dbIsReady = false;
        }

        CudaSW4() = delete;
        CudaSW4(const CudaSW4&) = delete;
        CudaSW4(CudaSW4&&) = default;
        CudaSW4& operator=(const CudaSW4&) = delete;
        CudaSW4& operator=(CudaSW4&&) = default;

        void setGapOpenScore(int score){
            if(verbose && score >= 0){
                std::cout << "Warning, gap open score set to non-negative value. Is this intended?\n";
            }
            gop = score;
        }
        void setGapExtendScore(int score){
            if(verbose && score >= 0){
                std::cout << "Warning, gap extend score set to non-negative value. Is this intended?\n";
            }
            gex = score;
        }

        void setDatabase(std::shared_ptr<DB> dbPtr){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReady();
        }

        void setDatabase(std::shared_ptr<DBWithVectors> dbPtr){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReady();
        }

        void setDatabase(std::shared_ptr<PseudoDB> dbPtr){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReady();
        }

        void setDatabase(std::shared_ptr<MMseqsDB> dbPtr){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReady();
        }

        void setDatabase(std::shared_ptr<ExternalDB> dbPtr){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReady();
        }

        void setDatabase(std::shared_ptr<DB> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReadyWithExistingFullGpuDB(existingFullGpuDBAllocations);
        }

        void setDatabase(std::shared_ptr<DBWithVectors> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReadyWithExistingFullGpuDB(existingFullGpuDBAllocations);
        }

        void setDatabase(std::shared_ptr<PseudoDB> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReadyWithExistingFullGpuDB(existingFullGpuDBAllocations);
        }

        void setDatabase(std::shared_ptr<MMseqsDB> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReadyWithExistingFullGpuDB(existingFullGpuDBAllocations);
        }

        void setDatabase(std::shared_ptr<ExternalDB> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReadyWithExistingFullGpuDB(existingFullGpuDBAllocations);
        }

        void setBlosum(BlosumType blosumType){
            setProgramWideBlosum(blosumType, deviceIds);
        }

        void setNumTop(int value){
            if(value > MaxNumberOfResults::value()){
                throw std::runtime_error("setNumTop: value too large");
            }
            if(value >= 0){
                numTop = value;
                updateNumResultsPerQuery();

                cub::SwitchDevice sd(deviceIds[0]);
                const int numGpus = deviceIds.size();           

                h_finalAlignmentScores.resize(results_per_query);
                h_finalReferenceIds.resize(results_per_query);
                h_finalEndPositions.resize(results_per_query);
                d_finalAlignmentScores_allGpus.resize(results_per_query * numGpus);
                d_finalReferenceIds_allGpus.resize(results_per_query * numGpus);  
                d_finalEndPositions_allGpus.resize(results_per_query * numGpus);              
            }
        }

        void setMemoryConfig(const MemoryConfig& val){
            memoryConfig = val;
        }

        std::vector<std::shared_ptr<GpuDatabaseAllocationBase>> getFullGpuDBAllocations(){
            if(!dbIsReady) return {};            

            prefetchDBToGpus();

            std::vector<std::shared_ptr<GpuDatabaseAllocationBase>> result;

            const int numGpus = deviceIds.size();
            for(int gpu = 0; gpu < numGpus; gpu++){
                auto& ws = *workingSets[gpu];
                
                result.push_back(ws.d_cacheddb);
            }

            return result;
        }

        std::string_view getReferenceHeader(ReferenceIdT referenceId) const{
            const auto& data = fullDB.getData();
            const char* const headerBegin = data.headers() + data.headerOffsets()[referenceId];
            const char* const headerEnd = data.headers() + data.headerOffsets()[referenceId+1];
            return std::string_view(headerBegin, std::distance(headerBegin, headerEnd));
        }

        int getReferenceLength(ReferenceIdT referenceId) const{
            const auto& data = fullDB.getData();
            return data.lengths()[referenceId];
        }

        std::string getReferenceSequence(ReferenceIdT referenceId) const{
            const auto& data = fullDB.getData();
            const char* const begin = data.chars() + data.offsets()[referenceId];
            const char* const end = begin + getReferenceLength(referenceId);

            std::string sequence(end - begin, '\0');
            std::transform(
                begin, 
                end,
                sequence.begin(),
                InverseConvertAA_20{}
            );

            return sequence;
        }

        void markCachedDBBatchesAsUploaded(int gpu){
            auto& ws = *workingSets[gpu];
            if(ws.getNumBatchesInCachedDB() > 0){
                batchPlansDstInfoVec_cachedDB[gpu][0].isUploaded = true;
                for(size_t i = 0; i < ws.getNumBatchesInCachedDB(); i++){
                    batchPlansDstInfoVec[gpu][i].isUploaded = true;
                }
            }
        }

        void prefetchDBToGpus(){
            nvtx::ScopedRange sr("prefetchDBToGpus", 1);
            RevertDeviceId rdi{};

            const int numGpus = deviceIds.size();
            std::vector<int> copyIds;

            helpers::CpuTimer copyTimer("transfer DB to GPUs");
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]);
                auto& ws = *workingSets[gpu];

                if(ws.getNumBatchesInCachedDB() > 0 && !batchPlansDstInfoVec_cachedDB[gpu][0].isUploaded){
                    const auto& plan = batchPlans_cachedDB[gpu][0];
                    const int currentBuffer = 0;
                    cudaStream_t H2DcopyStream = ws.copyStreams[currentBuffer];

                    executeCopyPlanH2DDirect(
                        plan,
                        ws.d_cacheddb->getCharData(),
                        ws.d_cacheddb->getLengthData(),
                        ws.d_cacheddb->getOffsetData(),
                        subPartitionsForGpus[gpu],
                        H2DcopyStream
                    );
                    
                    copyIds.push_back(gpu);

                    markCachedDBBatchesAsUploaded(gpu);
                }
            }
            for(int gpu : copyIds){
                cudaSetDevice(deviceIds[gpu]);
                cudaDeviceSynchronize(); CUERR;
            }
            copyTimer.stop();
            if(copyIds.size() > 0){
                if(verbose){
                    std::cout << "Transferred DB data in advance to GPU(s) ";
                    for(auto x : copyIds){
                        std::cout << x << " ";
                    }
                    std::cout << "\n";
                    copyTimer.print();
                }
            }
        }

        template<class QueryView>
        ScanResult scan(QueryView queryView, std::optional<const int8_t*> precomputedPssmOpt){
            nvtx::ScopedRange sr("scan", 6);
            if(!dbIsReady){
                throw std::runtime_error("DB not set correctly");
            }
            RevertDeviceId rdi{};

            // if(queryView.length <= getMaxSingleTileQueryLength_Gapless()){
            //     return ScanResult{};
            // }

            const int masterDeviceId = deviceIds[0];
            cudaSetDevice(masterDeviceId);

            scanTimer->reset();
            scanTimer->start();

            setQuery(queryView, precomputedPssmOpt);

            if(verbose && scanType == ScanType::GaplessPlusSW_Endpos && results_per_query == 0){
                std::cout << "Warning. Gapless+SW_Endpos selected, but results_per_query == 0\n";
            }

            switch(scanType){
                case ScanType::GaplessPlusSW_Endpos: 
                {
                    auto scanTypeOld = scanType;

                    scanType = ScanType::Gapless;
                    scanDatabaseForQuery_gapless();

                    std::vector<ReferenceIdT> vec(h_finalReferenceIds.begin(), h_finalReferenceIds.begin() + results_per_query);
                    // for(auto x : vec){ 
                    //     std::cout << x << " : " << getReferenceSequence(x) << "\n";
                    // }
                    // std::cout << "\n";

                    auto topIds = std::make_shared<TargetSubjectIds>(
                        std::move(vec)
                    );
                    setTargetSubjectIds(topIds);

                    scanType = ScanType::SW_Endpos;
                    scanDatabaseForQuery_sw_endpos();
                    setTargetSubjectIds(nullptr);

                    scanType = scanTypeOld;
                }; break;
                case ScanType::SW_Endpos: scanDatabaseForQuery_sw_endpos(); break;
                case ScanType::Gapless: //fall-through
                default: scanDatabaseForQuery_gapless(); break;
            }
            

            scanTimer->stop();

            totalProcessedQueryLengths += queryView.length;
            //totalNumOverflows += resultNumOverflows[0]; 

            const auto& sequenceLengthStatistics = getSequenceLengthStatistics();

            ScanResult result;
            size_t computedCells = 0;
            if(targetSubjectIds){
                for(const auto& x : targetSubjectIds->subjectIds){
                    computedCells += getReferenceLength(x);
                }
                computedCells *= queryView.length;
            }else{
                computedCells = sequenceLengthStatistics.sumOfLengths * queryView.length;
            }
            result.stats = makeBenchmarkStats(
                scanTimer->elapsed() / 1000, 
                computedCells, 
                0 //resultNumOverflows[0]
            );

            #ifdef CUDASW_DEBUG_CHECK_CORRECTNESS

            if(true || queryView.length > getMaxSingleTileQueryLength_Gapless()){
                std::vector<char> convertedQuery;            
                const char* query = queryView.ptr;
                if constexpr(QueryView::isEncoded){
                    convertedQuery.resize(queryView.length);
                    std::transform(
                        queryView.ptr,
                        queryView.ptr + queryView.length,
                        convertedQuery.data(),
                        InverseConvertAA_20{}
                    );
                    query = convertedQuery.data();
                }

                std::vector<int> cpuScores = computeAllScoresCPU_gaplessfilter_blosum62(query, queryView.length);
                //std::vector<int> cpuScoresExact = computeAllScoresCPU_exact_affine_blosum62(query, queryView.length);
                size_t numToCheck = cpuScores.size();

                //auto boundaries = getLengthPartitionBoundaries();
                
                //bool checkOk = true;
                int numErrors = 0;
                int numGreater2048 = 0;
                for(size_t i = 0; i < numToCheck; i++){
                    const auto refId = h_finalReferenceIds[i];
                    int gpu = h_finalAlignmentScores[i];
                    const int cpu = cpuScores[refId];


                    //if(getReferenceLength(refId) <= boundaries[boundaries.size() - 2]){
                        //gpu scores for all but last length partition are computed in half precision. don't report errors caused by rounding.
                        if(gpu > 2048){
                            gpu = cpu;
                            numGreater2048++;
                        }
                    //}
                    if(cpu != gpu){
                        if(numErrors == 0){
                            std::cout << "error. i " << i << ", sequence id " << refId 
                                << ", cpu score " << cpu << ", gpu score " << gpu 
                                //<< ", exact cpu score " << cpuScoresExact[refId] 
                                << ".";
                            std::cout << "Query:\n";
                            std::copy(query, query + queryView.length, std::ostream_iterator<char>(std::cout, ""));
                            std::cout << "\n";
                            std::cout << "db sequence:\n";
                            std::cout << getReferenceSequence(refId) << "\n";

                            // if(refId > 5 && refId < cpuScores.size() - 5){
                            //     for(int x = -5; x < 5; x++){
                            //         std::cout << getReferenceHeader(refId + x) << "\n";
                            //         std::cout << getReferenceSequence(refId + x) << "\n";
                            //     }
                            // }
                        }
                        numErrors++;
                    }
                }
                if(numErrors == 0){
                    std::cout << "Check ok, cpu and gpu produced same results. > 2048: " << numGreater2048 << " / " << numToCheck << "\n";
                }else{
                    std::cout << "Check not ok!!! " << numErrors << " sequences produced different results\n";
                }
            }

            //#endif
            #else

            result.scores.insert(result.scores.end(), h_finalAlignmentScores.begin(), h_finalAlignmentScores.begin() + results_per_query);
            result.referenceIds.insert(result.referenceIds.end(), h_finalReferenceIds.begin(), h_finalReferenceIds.begin() + results_per_query);
            result.endPositions.insert(result.endPositions.end(), h_finalEndPositions.begin(), h_finalEndPositions.begin() + results_per_query);
            
            #endif

            return result;
        }

        std::vector<int> computeAllScoresCPU_exact_affine_blosum62(const char* query, SequenceLengthT queryLength){
            const auto& view = fullDB.getData();
            size_t numSequences = view.numSequences();
            std::vector<int> result(numSequences);

            std::vector<char> convertedQuery(queryLength);
            auto convertQuery = [&](auto c){
                ConvertAA_20 charToMMseqs;
                ConvertAA_20_mmseqs_to_ncbi mmseqsToNcbi;
                return mmseqsToNcbi(charToMMseqs(c));
            };
            std::transform(
                query,
                query + queryLength,
                convertedQuery.data(),
                convertQuery
            );
            #pragma omp parallel for
            for(size_t i = 0; i < numSequences; i++){
                size_t offset = view.offsets()[i];
                int length = view.lengths()[i];
                const char* seq = view.chars() + offset;

                #if 1
                std::vector<char> convertedSubject(length);
                std::transform(
                    seq,
                    seq + length,
                    convertedSubject.data(),
                    ConvertAA_20_mmseqs_to_ncbi{}
                );
                seq = convertedSubject.data();
                #endif

                int score = affine_local_DP_host_protein_blosum62_converted(
                    convertedQuery.data(),
                    seq,
                    queryLength,
                    length,
                    gop,
                    gex
                );
                result[i] = score;
            }
            return result;
        }

        std::vector<int> computeAllScoresCPU_gaplessfilter_blosum62(const char* query, SequenceLengthT queryLength){
            const auto& view = fullDB.getData();
            size_t numSequences = view.numSequences();
            std::vector<int> result(numSequences);

            std::vector<char> convertedQuery(queryLength);
            auto convertQuery = [&](auto c){
                ConvertAA_20 charToMMseqs;
                ConvertAA_20_mmseqs_to_ncbi mmseqsToNcbi;
                return mmseqsToNcbi(charToMMseqs(c));
            };
            std::transform(
                query,
                query + queryLength,
                convertedQuery.data(),
                convertQuery
            );
            #pragma omp parallel for
            for(size_t i = 0; i < numSequences; i++){
                size_t offset = view.offsets()[i];
                int length = view.lengths()[i];
                const char* seq = view.chars() + offset;

                #if 1
                std::vector<char> convertedSubject(length);
                std::transform(
                    seq,
                    seq + length,
                    convertedSubject.data(),
                    ConvertAA_20_mmseqs_to_ncbi{}
                );
                seq = convertedSubject.data();
                #endif

                int score = GaplessFilter_host_protein_converted_blosum62(
                    convertedQuery.data(),
                    seq,
                    queryLength,
                    length
                );
                result[i] = score;
            }
            return result;
        }


        void printDBInfo() const{
            nvtx::ScopedRange sr("printDBInfo", 0);
            const size_t numSequences = fullDB.getData().numSequences();
            std::cout << numSequences << " sequences, " << fullDB.getData().numChars() << " characters\n";

            SequenceLengthStatistics stats = getSequenceLengthStatistics();

            std::cout << "Min length " << stats.min_length << ", max length " << stats.max_length 
                << ", avg length " << stats.sumOfLengths / numSequences << "\n";
        }

        void printDBLengthPartitions() const{
            auto lengthBoundaries = getLengthPartitionBoundaries();
            const int numLengthPartitions = getLengthPartitionBoundaries().size();

            for(int i = 0; i < numLengthPartitions; i++){
                std::cout << "<= " << lengthBoundaries[i] << ": " << fullDB_numSequencesPerLengthPartition[i] << "\n";
            }
        }

        void totalTimerStart(){
            RevertDeviceId rdi{};
            cudaSetDevice(deviceIds[0]);
            totalProcessedQueryLengths = 0;
            //totalNumOverflows = 0;
            totalTimer->start();
        }

        BenchmarkStats totalTimerStop(){
            RevertDeviceId rdi{};
            cudaSetDevice(deviceIds[0]);
            totalTimer->stop();

            const auto& sequenceLengthStatistics = getSequenceLengthStatistics();
            BenchmarkStats stats = makeBenchmarkStats(
                totalTimer->elapsed() / 1000,
                totalProcessedQueryLengths * sequenceLengthStatistics.sumOfLengths,
                0 //totalNumOverflows
            );

            return stats;
        }

        void setScanType(ScanType type){
            if(verbose){
                std::cout << "Set scan type to " << to_string(type) << "\n";
            }
            scanType = type;
        }

        void setTargetSubjectIds(std::shared_ptr<TargetSubjectIds> ptr){
            targetSubjectIds = ptr;

            if(targetSubjectIds){
                if(dbIsReady){
                    targetSubjectIds->removeOutOfBoundsTargets(fullDB.getData().numSequences());
                }
            }
        }

    private:
        void initializeGpus(){
            const int numGpus = deviceIds.size();

            for(int i = 0; i < numGpus; i++){
                cudaSetDevice(deviceIds[i]); CUERR
                helpers::init_cuda_context(); CUERR
                cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);CUERR
        
                cudaMemPool_t mempool;
                cudaDeviceGetDefaultMemPool(&mempool, deviceIds[i]); CUERR
                uint64_t threshold = UINT64_MAX;
                cudaMemPoolSetAttribute(mempool, cudaMemPoolAttrReleaseThreshold, &threshold);CUERR
            
                gpuStreams.emplace_back();
                gpuEvents.emplace_back(cudaEventDisableTiming);
            }
        }

        void makeReady(){
            nvtx::ScopedRange sr("makeReady", 0);
            #ifdef CUDASW_DEBUG_CHECK_CORRECTNESS
            const auto& dbData = fullDB.getData();
            size_t numDBSequences = dbData.numSequences();
            if(numDBSequences > size_t(std::numeric_limits<int>::max()))
                throw std::runtime_error("cannot check correctness for this db size");

            maxReduceArraySize = numDBSequences;
            results_per_query = maxReduceArraySize;
            setNumTopNoCheck(maxReduceArraySize);
            #endif

            dbSequenceLengthStatistics = nullptr;

            computeTotalNumSequencePerLengthPartition();
            partitionDBAmongstGpus();

            createDBBatchesForGpus();
            allocateGpuWorkingSets();
            assignBatchesToGpuMem();
            
            
            dbIsReady = true;
            updateNumResultsPerQuery();

            if(targetSubjectIds){
                targetSubjectIds->removeOutOfBoundsTargets(fullDB.getData().numSequences());
            }
        }

        void makeReadyWithExistingFullGpuDB(const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations){
            nvtx::ScopedRange sr("makeReadyWithExistingFullGpuDB", 0);
            #ifdef CUDASW_DEBUG_CHECK_CORRECTNESS
            const auto& dbData = fullDB.getData();
            size_t numDBSequences = dbData.numSequences();
            if(numDBSequences > size_t(std::numeric_limits<int>::max()))
                throw std::runtime_error("cannot check correctness for this db size");

            maxReduceArraySize = numDBSequences;
            results_per_query = maxReduceArraySize;
            setNumTopNoCheck(maxReduceArraySize);
            #endif

            dbSequenceLengthStatistics = nullptr;

            computeTotalNumSequencePerLengthPartition();
            partitionDBAmongstGpus();

            createDBBatchesForGpus();
            allocateGpuWorkingSetsWithExistingFullGpuDB(existingFullGpuDBAllocations);
            assignBatchesToGpuMem();

            const int numGpus = deviceIds.size();
            for(int gpu = 0; gpu < numGpus; gpu++){
                markCachedDBBatchesAsUploaded(gpu);
            }
            
            
            dbIsReady = true;
            updateNumResultsPerQuery();

            if(targetSubjectIds){
                targetSubjectIds->removeOutOfBoundsTargets(fullDB.getData().numSequences());
            }


            // const auto& data = fullDB.getData();

            // int pageableMemoryAccessUsesHostPageTables = 0;
            // int readOnlyHostRegisterSupported = 0;
            // cudaDeviceGetAttribute(&pageableMemoryAccessUsesHostPageTables, cudaDevAttrPageableMemoryAccessUsesHostPageTables, 0); CUERR;
            // std::cout << "pageableMemoryAccessUsesHostPageTables " << pageableMemoryAccessUsesHostPageTables << "\n";
            // //cudaDeviceGetAttribute(&readOnlyHostRegisterSupported, cudaDeviceAttrReadOnlyHostRegisterSupported, 0); CUERR;
            // // std::cout << "readOnlyHostRegisterSupported " << readOnlyHostRegisterSupported << "\n";

            // cudaDeviceProp prop;
            // cudaGetDeviceProperties(&prop, 0); CUERR;

            // std::cout << "prop.pageableMemoryAccess " << prop.pageableMemoryAccess << "\n";
            // std::cout << "prop.pageableMemoryAccessUsesHostPageTables " << prop.pageableMemoryAccessUsesHostPageTables << "\n";
            // std::cout << "prop.hostRegisterReadOnlySupported " << prop.hostRegisterReadOnlySupported << "\n";
            // std::cout << "prop.hostRegisterSupported " << prop.hostRegisterSupported << "\n";


            // cudaHostRegister((void*)data.chars(), sizeof(char) * data.numChars(), cudaHostRegisterDefault); CUERR;
            // cudaHostRegister((void*)data.lengths(), sizeof(SequenceLengthT) * data.numSequences(), cudaHostRegisterDefault); CUERR;
            // cudaHostRegister((void*)data.offsets(), sizeof(size_t) * data.numSequences(), cudaHostRegisterDefault); CUERR;
        }

        void computeTotalNumSequencePerLengthPartition(){
            nvtx::ScopedRange sr("computeTotalNumSequencePerLengthPartition", 1);
            auto lengthBoundaries = getLengthPartitionBoundaries();
            const int numLengthPartitions = getLengthPartitionBoundaries().size();

            fullDB_numSequencesPerLengthPartition.resize(numLengthPartitions);

            const auto& dbData = fullDB.getData();
            auto partitionBegin = dbData.lengths();
            for(int i = 0; i < numLengthPartitions; i++){
                //length k is in partition i if boundaries[i-1] < k <= boundaries[i]
                SequenceLengthT searchFor = lengthBoundaries[i];
                if(searchFor < std::numeric_limits<SequenceLengthT>::max()){
                    searchFor += 1;
                }
                auto partitionEnd = std::lower_bound(
                    partitionBegin, 
                    dbData.lengths() + dbData.numSequences(), 
                    searchFor
                );
                fullDB_numSequencesPerLengthPartition[i] = std::distance(partitionBegin, partitionEnd);
                partitionBegin = partitionEnd;
            }
        }

        void partitionDBAmongstGpus(){
            nvtx::ScopedRange sr("partitionDBAmongstGpus", 2);
            const int numGpus = deviceIds.size();
            const int numLengthPartitions = getLengthPartitionBoundaries().size();

            numSequencesPerLengthPartitionPrefixSum.clear();
            dbPartitionsByLengthPartitioning.clear();
            subPartitionsForGpus.clear();
            lengthPartitionIdsForGpus.clear();
            numSequencesPerGpu.clear();
            numSequencesPerGpuPrefixSum.clear();

            const auto& data = fullDB.getData();
    
            subPartitionsForGpus.resize(numGpus);
            lengthPartitionIdsForGpus.resize(numGpus);
            numSequencesPerGpu.resize(numGpus, 0);
            numSequencesPerGpuPrefixSum.resize(numGpus, 0);
    
            numSequencesPerLengthPartitionPrefixSum.resize(numLengthPartitions, 0);
            for(int i = 0; i < numLengthPartitions-1; i++){
                numSequencesPerLengthPartitionPrefixSum[i+1] = numSequencesPerLengthPartitionPrefixSum[i] + fullDB_numSequencesPerLengthPartition[i];
            }
    
            for(int i = 0; i < numLengthPartitions; i++){
                size_t begin = numSequencesPerLengthPartitionPrefixSum[i];
                size_t end = begin + fullDB_numSequencesPerLengthPartition[i];
                dbPartitionsByLengthPartitioning.emplace_back(data, begin, end);        
            }
    
            for(int lengthPartitionId = 0; lengthPartitionId < numLengthPartitions; lengthPartitionId++){
                const auto& lengthPartition = dbPartitionsByLengthPartitioning[lengthPartitionId];        
                const auto partitionedByGpu = partitionDBdata_by_numberOfChars(lengthPartition, lengthPartition.numChars() / numGpus);
        
                assert(int(partitionedByGpu.size()) <= numGpus);
                for(int gpu = 0; gpu < numGpus; gpu++){
                    if(gpu < int(partitionedByGpu.size())){
                        subPartitionsForGpus[gpu].push_back(partitionedByGpu[gpu]);
                        lengthPartitionIdsForGpus[gpu].push_back(lengthPartitionId);
                    }else{
                        //add empty partition
                        subPartitionsForGpus[gpu].push_back(DBdataView(data, 0, 0));
                        lengthPartitionIdsForGpus[gpu].push_back(0);
                    }
                }
            }
        
            for(int i = 0; i < numGpus; i++){
                for(const auto& p : subPartitionsForGpus[i]){
                    numSequencesPerGpu[i] += p.numSequences();
                }
            }
            for(int i = 0; i < numGpus-1; i++){
                numSequencesPerGpuPrefixSum[i+1] = numSequencesPerGpuPrefixSum[i] + numSequencesPerGpu[i];
            }
        
            numSequencesPerGpu_total.resize(numGpus);
            numSequencesPerGpuPrefixSum_total.resize(numGpus);
            numSequencesPerGpuPrefixSum_total[0] = 0;

        
            for(int i = 0; i < numGpus; i++){
                size_t num = numSequencesPerGpu[i];
                numSequencesPerGpu_total[i] = num;
                if(i < numGpus - 1){
                    numSequencesPerGpuPrefixSum_total[i+1] = numSequencesPerGpuPrefixSum_total[i] + num;
                }
            }

            std::vector<size_t> sequencesInPartitions(numGpus * numLengthPartitions);
            for(int gpu = 0; gpu < numGpus; gpu++){
                assert(subPartitionsForGpus[gpu].size() == numLengthPartitions);
                for(int i = 0; i < numLengthPartitions; i++){
                    sequencesInPartitions[gpu * numLengthPartitions + i] = subPartitionsForGpus[gpu][i].numSequences();
                }
            }
            hostGpuPartitionOffsets = HostGpuPartitionOffsets(numGpus, numLengthPartitions, std::move(sequencesInPartitions));
        }

        void allocateGpuWorkingSets(){
            nvtx::ScopedRange sr("allocateGpuWorkingSets", 3);
            const int numGpus = deviceIds.size();
            workingSets.clear();
            workingSets.resize(numGpus);

            if(verbose){
                std::cout << "Allocate Memory: \n";
            }
            //nvtx::push_range("ALLOC_MEM", 0);
            helpers::CpuTimer allocTimer("ALLOC_MEM");

            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]);

                size_t freeMem, totalMem;
                cudaMemGetInfo(&freeMem, &totalMem);
                constexpr size_t safety = 256*1024*1024;
                size_t memlimit = std::min(freeMem, memoryConfig.maxGpuMem);
                if(memlimit > safety){
                    memlimit -= safety;
                }

                if(verbose){
                    std::cout << "gpu " << gpu << " may use " << memlimit << " bytes. ";
                }

                const bool needsPinnedStagingBuffers = numGpus > 1;

                workingSets[gpu] = std::make_unique<GpuWorkingSet>(
                    memlimit,
                    memoryConfig.maxBatchBytes,
                    memoryConfig.maxBatchSequences,
                    memoryConfig.maxTempBytes,
                    subPartitionsForGpus[gpu],
                    batchPlans[gpu],
                    needsPinnedStagingBuffers,
                    maxReduceArraySize
                );

                if(verbose){
                    std::cout << "Using " << workingSets[gpu]->numTempBytes << " temp bytes. ";
                }
                if(verbose){
                    std::cout << workingSets[gpu]->getNumBatchesInCachedDB() << " out of " << batchPlans[gpu].size() << " DB batches will be cached in gpu memory\n";
                }

                //set gpu partition table
                workingSets[gpu]->setPartitionOffsets(hostGpuPartitionOffsets);

                const bool usesCallbackThread = numGpus > 1;

                if(usesCallbackThread){
                    //spin up the host callback thread
                    auto noop = [](void*){};
                    cudaLaunchHostFunc(
                        gpuStreams[gpu], 
                        noop, 
                        nullptr
                    ); CUERR
                }

                workingSets[gpu]->setTopNSize(results_per_query);
            }    

            if(verbose){
                allocTimer.print();
            }
        }

        void allocateGpuWorkingSetsWithExistingFullGpuDB(const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations){
            nvtx::ScopedRange sr("allocateGpuWorkingSetsWithExistingFullGpuDB", 3);
            const int numGpus = deviceIds.size();
            workingSets.clear();
            workingSets.resize(numGpus);

            if(verbose){
                std::cout << "Allocate Memory: \n";
            }
            //nvtx::push_range("ALLOC_MEM", 0);
            helpers::CpuTimer allocTimer("ALLOC_MEM");

            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]);

                size_t freeMem, totalMem;
                cudaMemGetInfo(&freeMem, &totalMem);
                constexpr size_t safety = 256*1024*1024;
                size_t memlimit = std::min(freeMem, memoryConfig.maxGpuMem);
                if(memlimit > safety){
                    memlimit -= safety;
                }

                if(verbose){
                    std::cout << "gpu " << gpu << " may use " << memlimit << " bytes. ";
                }

                const bool needsPinnedStagingBuffers = numGpus > 1;

                workingSets[gpu] = std::make_unique<GpuWorkingSet>(
                    memlimit,
                    memoryConfig.maxBatchBytes,
                    memoryConfig.maxBatchSequences,
                    memoryConfig.maxTempBytes,
                    subPartitionsForGpus[gpu],
                    batchPlans[gpu],
                    existingFullGpuDBAllocations[gpu],
                    needsPinnedStagingBuffers,
                    maxReduceArraySize
                );

                if(verbose){
                    std::cout << "Using " << workingSets[gpu]->numTempBytes << " temp bytes. ";
                }
                if(verbose){
                    std::cout << workingSets[gpu]->getNumBatchesInCachedDB() << " out of " << batchPlans[gpu].size() << " DB batches will be cached in gpu memory\n";
                }

                //set gpu partition table
                workingSets[gpu]->setPartitionOffsets(hostGpuPartitionOffsets);

                const bool usesCallbackThread = numGpus > 1;

                if(usesCallbackThread){
                    //spin up the host callback thread
                    auto noop = [](void*){};
                    cudaLaunchHostFunc(
                        gpuStreams[gpu], 
                        noop, 
                        nullptr
                    ); CUERR
                }

                workingSets[gpu]->setTopNSize(results_per_query);
            }    

            if(verbose){
                allocTimer.print();
            }
        }

        

        void createDBBatchesForGpus(){
            nvtx::ScopedRange sr("createDBBatchesForGpus", 4);
            const int numGpus = deviceIds.size();

            batchPlans.clear();
            batchPlans.resize(numGpus);
            batchPlans_cachedDB.clear();
            batchPlans_cachedDB.resize(numGpus);
    
            for(int gpu = 0; gpu < numGpus; gpu++){
                batchPlans[gpu] = computeDbCopyPlan(
                    subPartitionsForGpus[gpu],
                    lengthPartitionIdsForGpus[gpu],
                    memoryConfig.maxBatchBytes,
                    memoryConfig.maxBatchSequences
                );
                if(verbose){
                    std::cout << "Batch plan gpu " << gpu << ": " << batchPlans[gpu].size() << " batches\n";
                }
            }
        }

        void assignBatchesToGpuMem(){
            nvtx::ScopedRange sr("createDBBatchesForGpus", 5);
            const int numGpus = deviceIds.size();
            batchPlansDstInfoVec.clear();
            batchPlansDstInfoVec.resize(numGpus);
            batchPlansDstInfoVec_cachedDB.clear();
            batchPlansDstInfoVec_cachedDB.resize(numGpus);
    
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]);
                auto& ws = *workingSets[gpu];
                if(ws.getNumBatchesInCachedDB() > 0){
                    //can cache full db in gpu mem

                    auto plansForCachedDB = computeDbCopyPlan(
                        subPartitionsForGpus[gpu],
                        lengthPartitionIdsForGpus[gpu],
                        sizeof(char) * ws.getNumCharsInCachedDB(),
                        ws.getNumSequencesInCachedDB()
                    );
                    assert(plansForCachedDB.size() >= 1);
                    plansForCachedDB.erase(plansForCachedDB.begin() + 1, plansForCachedDB.end());
                    batchPlans_cachedDB[gpu] = plansForCachedDB;
                    // if(verbose){
                    //     std::cout << "Cached db single batch plan " << plansForCachedDB[0] << "\n";
                    // }

                    BatchDstInfo dstInfo;
                    dstInfo.isUploaded = false;
                    dstInfo.charsPtr = ws.d_cacheddb->getCharData();
                    dstInfo.lengthsPtr = ws.d_cacheddb->getLengthData();
                    dstInfo.offsetsPtr = ws.d_cacheddb->getOffsetData();
                    batchPlansDstInfoVec_cachedDB[gpu].push_back(dstInfo);
                }

                {
                    BatchDstInfo dstInfo;
                    dstInfo.isUploaded = false;
                    dstInfo.charsPtr = ws.d_cacheddb->getCharData();
                    dstInfo.lengthsPtr = ws.d_cacheddb->getLengthData();
                    dstInfo.offsetsPtr = ws.d_cacheddb->getOffsetData();

                    for(size_t i = 0; i < ws.getNumBatchesInCachedDB(); i++){
                        batchPlansDstInfoVec[gpu].push_back(dstInfo);
                        const auto& plan = batchPlans[gpu][i];
                        dstInfo.charsPtr += plan.usedBytes;
                        dstInfo.lengthsPtr += plan.usedSeq;
                        dstInfo.offsetsPtr += plan.usedSeq;
                    }

                    for(size_t i = ws.getNumBatchesInCachedDB(), buf = 0; i < batchPlans[gpu].size(); i++, buf = (buf+1)%ws.numCopyBuffers){
                        dstInfo.charsPtr = ws.d_chardata_vec[buf].data();
                        dstInfo.lengthsPtr = ws.d_lengthdata_vec[buf].data();
                        dstInfo.offsetsPtr = ws.d_offsetdata_vec[buf].data();
                        batchPlansDstInfoVec[gpu].push_back(dstInfo);
                    }
                }
            }
        }

        void printDBDataView(const DBdataView& view) const{
            std::cout << "Sequences: " << view.numSequences() << "\n";
            std::cout << "Chars: " << view.offsets()[0] << " - " << view.offsets()[view.numSequences()] << " (" << (view.offsets()[view.numSequences()] - view.offsets()[0]) << ")"
                << " " << view.numChars() << "\n";
        }

        void printDBDataViews(const std::vector<DBdataView>& views) const {
            size_t numViews = views.size();
            for(size_t p = 0; p < numViews; p++){
                const DBdataView& view = views[p];
        
                std::cout << "View " << p << "\n";
                printDBDataView(view);
            }
        }

        SequenceLengthStatistics getSequenceLengthStatistics() const{
            if(dbSequenceLengthStatistics == nullptr){
                dbSequenceLengthStatistics = std::make_unique<SequenceLengthStatistics>();
                const auto& data = fullDB.getData();
                size_t numSeq = data.numSequences();

                for (size_t i=0; i < numSeq; i++) {
                    if (data.lengths()[i] > dbSequenceLengthStatistics->max_length) dbSequenceLengthStatistics->max_length = data.lengths()[i];
                    if (data.lengths()[i] < dbSequenceLengthStatistics->min_length) dbSequenceLengthStatistics->min_length = data.lengths()[i];
                    dbSequenceLengthStatistics->sumOfLengths += data.lengths()[i];
                }
            }
            return *dbSequenceLengthStatistics;
        }

        std::vector<DeviceBatchCopyToPinnedPlan> computeDbCopyPlan(
            const std::vector<DBdataView>& dbPartitions,
            const std::vector<int>& lengthPartitionIds,
            size_t MAX_CHARDATA_BYTES,
            size_t MAX_SEQ
        ) const {
            std::vector<DeviceBatchCopyToPinnedPlan> result;
        
            size_t currentCopyPartition = 0;
            size_t currentCopySeqInPartition = 0;
        
            //size_t processedSequences = 0;
            while(currentCopyPartition < dbPartitions.size()){
                
                size_t usedBytes = 0;
                size_t usedSeq = 0;
        
                DeviceBatchCopyToPinnedPlan plan;
        
                while(currentCopyPartition < dbPartitions.size()){
                    if(dbPartitions[currentCopyPartition].numSequences() == 0){
                        currentCopyPartition++;
                        continue;
                    }
        
                    //figure out how many sequences to copy to pinned
                    size_t remainingBytes = MAX_CHARDATA_BYTES - usedBytes;
                    
                    auto dboffsetsBegin = dbPartitions[currentCopyPartition].offsets() + currentCopySeqInPartition;
                    auto dboffsetsEnd = dbPartitions[currentCopyPartition].offsets() + dbPartitions[currentCopyPartition].numSequences() + 1;
                    
                    auto searchFor = dbPartitions[currentCopyPartition].offsets()[currentCopySeqInPartition] + remainingBytes + 1; // +1 because remainingBytes is inclusive
                    auto it = std::lower_bound(
                        dboffsetsBegin,
                        dboffsetsEnd,
                        searchFor
                    );
        
                    size_t numToCopyByBytes = 0;
                    if(it != dboffsetsBegin){
                        numToCopyByBytes = std::distance(dboffsetsBegin, it) - 1;
                    }
                    if(numToCopyByBytes == 0 && currentCopySeqInPartition == 0){
                        std::cout << "Warning. copy buffer size too small. skipped a db portion\n";
                        break;
                    }
                    
                    size_t remainingSeq = MAX_SEQ - usedSeq;            
                    size_t numToCopyBySeq = std::min(dbPartitions[currentCopyPartition].numSequences() - currentCopySeqInPartition, remainingSeq);
                    size_t numToCopy = std::min(numToCopyByBytes,numToCopyBySeq);
        
                    if(numToCopy > 0){
                        DeviceBatchCopyToPinnedPlan::CopyRange copyRange;
                        copyRange.lengthPartitionId = lengthPartitionIds[currentCopyPartition];
                        copyRange.currentCopyPartition = currentCopyPartition;
                        copyRange.currentCopySeqInPartition = currentCopySeqInPartition;
                        copyRange.numToCopy = numToCopy;
                        plan.copyRanges.push_back(copyRange);
        
                        if(usedSeq == 0){
                            plan.h_partitionIds.push_back(lengthPartitionIds[currentCopyPartition]);
                            plan.h_numPerPartition.push_back(numToCopy);
                        }else{
                            //if is same length partition as previous copy 
                            if(plan.h_partitionIds.back() == lengthPartitionIds[currentCopyPartition]){
                                plan.h_numPerPartition.back() += numToCopy;
                            }else{
                                //new length partition
                                plan.h_partitionIds.push_back(lengthPartitionIds[currentCopyPartition]);
                                plan.h_numPerPartition.push_back(numToCopy);
                            }
                        }
                        usedBytes += (dbPartitions[currentCopyPartition].offsets()[currentCopySeqInPartition+numToCopy] 
                            - dbPartitions[currentCopyPartition].offsets()[currentCopySeqInPartition]);
                        usedSeq += numToCopy;
        
                        currentCopySeqInPartition += numToCopy;
                        if(currentCopySeqInPartition == dbPartitions[currentCopyPartition].numSequences()){
                            currentCopySeqInPartition = 0;
                            currentCopyPartition++;
                        }
                    }else{
                        break;
                    }
                }
        
                plan.usedBytes = usedBytes;
                plan.usedSeq = usedSeq;    
                
                if(usedSeq == 0 && currentCopyPartition < dbPartitions.size() && dbPartitions[currentCopyPartition].numSequences() > 0){
                    std::cout << "Warning. copy buffer size too small. skipped a db portion. stop\n";
                    break;
                }
        
                if(plan.usedSeq > 0){
                    result.push_back(plan);
                }
            }
        
            return result;
        }

        template<class QueryView>
        void setQuery(QueryView queryView, std::optional<const int8_t*> precomputedPssmOpt){
            nvtx::ScopedRange sr("setQuery", 0);

            const int queryLength = queryView.length;
            // const char* query = queryView.ptr;

            if(queryLength > MaxSequenceLength::value()){
                std::string msg = "Query length is " + std::to_string(queryLength) 
                    + ", but config allows only lengths <= " + std::to_string(MaxSequenceLength::value());
                throw std::runtime_error(msg);
            }
            
            currentQueryLength = queryLength;
            //pad query to multiple of 4 for char4 access
            //add sizeof(char4) * warpsize for unguarded accesses outside of the DP matrix
            currentQueryLengthWithPadding = SDIV(queryLength, 4) * 4 + sizeof(char4) * 32;

            PSSM hostFullQueryPSSM = [&](){
                if(precomputedPssmOpt.has_value()){
                    const int8_t* precomputedPssm = precomputedPssmOpt.value();
                    if(precomputedPssm == nullptr) throw std::runtime_error("setQuery pssm is nullptr");
                    return PSSM::fromPSSM(queryView.ptr, queryView.length, precomputedPssm, 21);
                }else{
                    if constexpr(QueryView::isEncoded){
                        return PSSM::fromBlosum(blosumType, queryView.ptr, queryView.length);
                    }else{
                        std::vector<char> currentQueryEncodedHost(queryView.length);
                        std::transform(
                            queryView.ptr,
                            queryView.ptr + queryView.length,
                            currentQueryEncodedHost.begin(),
                            ConvertAA_20{}
                        );
                        return PSSM::fromBlosum(blosumType, currentQueryEncodedHost.data(), currentQueryEncodedHost.size());
                    }
                }                
            }();

            // std::cout << "hostFullQueryPSSM\n";
            // for(int r = 0; r < 21; r++){
            //     for(int l = 0; l < queryLength; l++){
            //         std::cout << hostFullQueryPSSM[r][l] << " ";
            //     }
            //     std::cout << "\n";
            // }

            const int numGpus = deviceIds.size();
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                auto& ws = *workingSets[gpu];
                ws.d_query.resize(currentQueryLengthWithPadding);
                cudaMemsetAsync(ws.d_query.data() + currentQueryLength, 20, currentQueryLengthWithPadding - currentQueryLength, gpuStreams[gpu]);
                cudaMemcpyAsync(ws.d_query.data(), queryView.ptr, currentQueryLength, cudaMemcpyDefault, gpuStreams[gpu]); CUERR

                if constexpr(!QueryView::isEncoded){
                    thrust::transform(
                        thrust::cuda::par_nosync.on(gpuStreams[gpu]),
                        ws.d_query.data(),
                        ws.d_query.data() + currentQueryLength,
                        ws.d_query.data(),
                        ConvertAA_20{}
                    );
                }

                // std::vector<char> tmpvec(currentQueryLength);
                // cudaMemcpy(tmpvec.data(), ws.d_query.data(), sizeof(char) * currentQueryLength, cudaMemcpyDeviceToHost);
                // std::transform(
                //     tmpvec.data(),
                //     tmpvec.data() + queryView.length,
                //     tmpvec.data(),
                //     ConvertAA_20_mmseqs_to_ncbi{}
                // );
                // std::cout << "ws.d_query: ";
                // for(auto x : tmpvec){
                //     std::cout << int(x) << " ";
                // }
                // std::cout << "\n";

                ws.gpuFullQueryPSSM.upload(hostFullQueryPSSM, gpuStreams[gpu]);

                auto makeGaplessPSSM = [&](){
                    if(currentQueryLength <= getMaxSingleTileQueryLength_Gapless()){
                        auto config = getSingleTileGroupRegConfigForPSSM_Gapless(currentQueryLength);
                        if(verbose){
                            std::cout << "Query length " << currentQueryLength << ". Set up PSSM for single-tile processing. "
                                << "Tilesize " << (config.groupsize * config.numRegs * 2) << " = " << config.groupsize << " * " << config.numRegs << " * 2" 
                                ", dpx: " << config.dpx << ", approach: " << to_string(config.approach) << "\n";
                        }
                        constexpr int accessSize = 16; //kernel uses float4 for pssm access
                        if(!config.dpx){
                            ws.gpuPermutedPSSMforGapless.template fromGpuPSSMView<half, accessSize>(ws.gpuFullQueryPSSM.makeView(), config.groupsize, config.numRegs, gpuStreams[gpu]);
                        }else{
                            ws.gpuPermutedPSSMforGapless.template fromGpuPSSMView<short, accessSize>(ws.gpuFullQueryPSSM.makeView(), config.groupsize, config.numRegs, gpuStreams[gpu]);
                        }
                    }else{
                        auto config = getMultiTileGroupRegConfigForPSSM_Gapless(currentQueryLength);
                        if(verbose){
                            std::cout << "Query length " << currentQueryLength << ". Set up PSSM for multi-tile processing. "
                                << "Tilesize " << (config.groupsize * config.numRegs * 2) << " = " << config.groupsize << " * " << config.numRegs << " * 2" 
                                ", dpx: " << config.dpx << ", approach: " << to_string(config.approach) << "\n";
                        }
                        constexpr int accessSize = 16; //kernel uses float4 for pssm access
                        if(!config.dpx){
                            ws.gpuPermutedPSSMforGapless.template fromGpuPSSMView<half, accessSize>(ws.gpuFullQueryPSSM.makeView(), config.groupsize, config.numRegs, gpuStreams[gpu]);
                        }else{
                            ws.gpuPermutedPSSMforGapless.template fromGpuPSSMView<short, accessSize>(ws.gpuFullQueryPSSM.makeView(), config.groupsize, config.numRegs, gpuStreams[gpu]);
                        }
                    }
                };

                auto makeSWPSSM = [&](){
                    if(currentQueryLength <= getMaxSingleTileQueryLength_SW()){
                        auto config = getSingleTileGroupRegConfigForPSSM_SW(currentQueryLength);
                        if(verbose){
                            std::cout << "Query length " << currentQueryLength << ". Set up PSSM for single-tile processing. "
                                << "Tilesize " << (config.groupsize * config.numRegs) << " = " << config.groupsize << " * " << config.numRegs 
                                << ", dpx: " << config.dpx << ", approach: " << to_string(config.approach) << "\n";
                        }
                        constexpr int accessSize = 16; //kernel uses float4 for pssm access                        
                        if(!config.dpx){
                            ws.gpuPermutedPSSMforSW.template fromGpuPSSMView<accessSize, float>(ws.gpuFullQueryPSSM.makeView(), config.groupsize, config.numRegs, gpuStreams[gpu]);
                        }else{
                            ws.gpuPermutedPSSMforSW.template fromGpuPSSMView<accessSize, int>(ws.gpuFullQueryPSSM.makeView(), config.groupsize, config.numRegs, gpuStreams[gpu]);
                        }
                    }else{
                        auto config = getMultiTileGroupRegConfigForPSSM_SW(currentQueryLength);
                        if(verbose){
                            std::cout << "Query length " << currentQueryLength << ". Set up PSSM for single-tile processing. "
                                << "Tilesize " << (config.groupsize * config.numRegs) << " = " << config.groupsize << " * " << config.numRegs 
                                << ", dpx: " << config.dpx << ", approach: " << to_string(config.approach) << "\n";
                        }
                        constexpr int accessSize = 16; //kernel uses float4 for pssm access                        
                        if(!config.dpx){
                            ws.gpuPermutedPSSMforSW.template fromGpuPSSMView<accessSize, float>(ws.gpuFullQueryPSSM.makeView(), config.groupsize, config.numRegs, gpuStreams[gpu]);
                        }else{
                            ws.gpuPermutedPSSMforSW.template fromGpuPSSMView<accessSize, int>(ws.gpuFullQueryPSSM.makeView(), config.groupsize, config.numRegs, gpuStreams[gpu]);
                        }
                    }
                };

                if(scanType == ScanType::Gapless){
                    makeGaplessPSSM();
                }else if(scanType == ScanType::SW_Endpos){
                    makeSWPSSM();
                }else if(scanType == ScanType::GaplessPlusSW_Endpos){
                    makeGaplessPSSM();
                    makeSWPSSM();
                }

                //THIS cudaMemcpyToSymbolAsync IS ONLY REQUIRED FOR THE NON-PSSM KERNELS

                //TODO leave query in gmem, dont use cmem ???
                // cudaMemcpyToSymbolAsync(constantQuery4, ws.d_query.data(), currentQueryLength, 0, cudaMemcpyDeviceToDevice, gpuStreams[gpu]); CUERR

            }


            // for(int subjectLetter = 0; subjectLetter < hostFullQueryPSSM.alphabetSize; subjectLetter++){
            //     for(int queryLetter = 0; queryLetter < hostFullQueryPSSM.queryLength; queryLetter++){
            //         std::cout << hostFullQueryPSSM[subjectLetter][queryLetter] << " ";
            //     }
            //     std::cout << "\n";
            // }
            // std::cout << "\n";
            // std::exit(0);
        }

        void scanDatabaseForQuery_gapless(){
            nvtx::ScopedRange sr("scanDatabaseForQuery_gapless", 0);
            const int numGpus = deviceIds.size();
            const int masterDeviceId = deviceIds[0];
            const auto& masterStream1 = gpuStreams[0];
            auto& masterevent1 = gpuEvents[0];

            cudaSetDevice(masterDeviceId);
            // scanTimer->reset();
            // scanTimer->start();

            thrust::fill(
                thrust::cuda::par_nosync.on(masterStream1),
                d_finalAlignmentScores_allGpus.begin(),
                d_finalAlignmentScores_allGpus.end(),
                0
            );

            cudaSetDevice(masterDeviceId);           

            cudaEventRecord(masterevent1, masterStream1); CUERR;

            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                cudaStreamWaitEvent(gpuStreams[gpu], masterevent1, 0); CUERR;
            }

            if(!targetSubjectIds){
                processQueryOnGpus();
            }else{
                processQueryOnGpusWithTargetSubjectIds();
            }

            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                auto& ws = *workingSets[gpu];

                if(numGpus > 1){
                    //transform per gpu local sequence indices into global sequence indices
                    if(results_per_query > 0){
                        transformLocalSequenceIndicesToGlobalIndices<<<SDIV(results_per_query, 128), 128, 0, gpuStreams[gpu]>>>(
                            gpu,
                            results_per_query,
                            ws.deviceGpuPartitionOffsets.getDeviceView(),
                            ws.d_topN_refIds.data()
                        ); CUERR;
                    }
                }

                cudaMemcpyAsync(
                    d_finalAlignmentScores_allGpus.data() + results_per_query*gpu,
                    ws.d_topN_scores.data(),
                    sizeof(float) * results_per_query,
                    cudaMemcpyDeviceToDevice,
                    gpuStreams[gpu]
                ); CUERR;
                cudaMemcpyAsync(
                    d_finalReferenceIds_allGpus.data() + results_per_query*gpu,
                    ws.d_topN_refIds.data(),
                    sizeof(ReferenceIdT) * results_per_query,
                    cudaMemcpyDeviceToDevice,
                    gpuStreams[gpu]
                ); CUERR;                
                // cudaMemcpyAsync(
                //     d_resultNumOverflows.data() + gpu,
                //     ws.d_total_overflow_number.data(),
                //     sizeof(int),
                //     cudaMemcpyDeviceToDevice,
                //     gpuStreams[gpu]
                // ); CUERR;                

                cudaEventRecord(ws.forkStreamEvent, gpuStreams[gpu]); CUERR;

                cudaSetDevice(masterDeviceId);
                cudaStreamWaitEvent(masterStream1, ws.forkStreamEvent, 0); CUERR;
            }

            cudaSetDevice(masterDeviceId);

            if(numGpus > 1){
                //sort per-gpu top results to find overall top results
                auto sortInput = thrust::make_zip_iterator(
                    d_finalAlignmentScores_allGpus.begin(),
                    d_finalReferenceIds_allGpus.begin()
                );
                thrust::sort(
                    thrust::cuda::par_nosync(thrust_async_allocator<char>(masterStream1)).on(masterStream1),
                    sortInput,
                    sortInput + results_per_query * numGpus,
                    CompareScoresDescendingRefIdsAscending{}
                );

                // thrust::sort_by_key(
                //     thrust::cuda::par_nosync(thrust_async_allocator<char>(masterStream1)).on(masterStream1),
                //     d_finalAlignmentScores_allGpus.begin(),
                //     d_finalAlignmentScores_allGpus.begin() + results_per_query * numGpus,
                //     d_finalReferenceIds_allGpus.begin(),
                //     thrust::greater<float>()
                // );


                //sum the overflows per gpu
                //sumNumOverflowsKernel<<<1,1,0,masterStream1>>>(d_resultNumOverflows.data(), d_resultNumOverflows.data(), numGpus); CUERR;                
            }

            cudaMemcpyAsync(
                h_finalAlignmentScores.data(), 
                d_finalAlignmentScores_allGpus.data(), 
                sizeof(float) * results_per_query, 
                cudaMemcpyDeviceToHost, 
                masterStream1
            );  CUERR
            cudaMemcpyAsync(
                h_finalReferenceIds.data(), 
                d_finalReferenceIds_allGpus.data(), 
                sizeof(ReferenceIdT) * results_per_query, 
                cudaMemcpyDeviceToHost, 
                masterStream1
            );  CUERR
            // cudaMemcpyAsync(
            //     resultNumOverflows.data(), 
            //     d_resultNumOverflows.data(), 
            //     sizeof(int), 
            //     cudaMemcpyDeviceToHost, 
            //     masterStream1
            // );  CUERR

            cudaStreamSynchronize(masterStream1); CUERR;

            if(targetSubjectIds){
                //h_finalReferenceIds will contain numbers from 0 to num target subject ids. convert to proper target subject ids
                for(int i = 0; i < results_per_query; i++){
                    h_finalReferenceIds[i] = targetSubjectIds->subjectIds[h_finalReferenceIds[i]];
                }
            }
        }


        void scanDatabaseForQuery_sw_endpos(){
            nvtx::ScopedRange sr("scanDatabaseForQuery_sw_endpos", 0);
            const int numGpus = deviceIds.size();
            const int masterDeviceId = deviceIds[0];
            const auto& masterStream1 = gpuStreams[0];
            auto& masterevent1 = gpuEvents[0];

            cudaSetDevice(masterDeviceId);
            // scanTimer->reset();
            // scanTimer->start();

            thrust::fill(
                thrust::cuda::par_nosync.on(masterStream1),
                d_finalAlignmentScores_allGpus.begin(),
                d_finalAlignmentScores_allGpus.end(),
                0
            );

            cudaSetDevice(masterDeviceId);           

            cudaEventRecord(masterevent1, masterStream1); CUERR;

            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                cudaStreamWaitEvent(gpuStreams[gpu], masterevent1, 0); CUERR;
            }

            if(!targetSubjectIds){
                processQueryOnGpus();
            }else{
                processQueryOnGpusWithTargetSubjectIds();
            }

            if(!targetSubjectIds){


                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];

                    if(numGpus > 1){
                        //transform per gpu local sequence indices into global sequence indices
                        if(results_per_query > 0){
                            transformLocalSequenceIndicesToGlobalIndices<<<SDIV(results_per_query, 128), 128, 0, gpuStreams[gpu]>>>(
                                gpu,
                                results_per_query,
                                ws.deviceGpuPartitionOffsets.getDeviceView(),
                                ws.d_topN_refIds.data()
                            ); CUERR;
                        }
                    }

                    cudaMemcpyAsync(
                        d_finalAlignmentScores_allGpus.data() + results_per_query*gpu,
                        ws.d_topN_scores.data(),
                        sizeof(float) * results_per_query,
                        cudaMemcpyDeviceToDevice,
                        gpuStreams[gpu]
                    ); CUERR;
                    cudaMemcpyAsync(
                        d_finalReferenceIds_allGpus.data() + results_per_query*gpu,
                        ws.d_topN_refIds.data(),
                        sizeof(ReferenceIdT) * results_per_query,
                        cudaMemcpyDeviceToDevice,
                        gpuStreams[gpu]
                    ); CUERR;
                    cudaMemcpyAsync(
                        d_finalEndPositions_allGpus.data() + results_per_query*gpu,
                        ws.d_topN_alignmentEndPositions.data(),
                        sizeof(AlignmentEndPosition) * results_per_query,
                        cudaMemcpyDeviceToDevice,
                        gpuStreams[gpu]
                    ); CUERR;  
                    // cudaMemcpyAsync(
                    //     d_resultNumOverflows.data() + gpu,
                    //     ws.d_total_overflow_number.data(),
                    //     sizeof(int),
                    //     cudaMemcpyDeviceToDevice,
                    //     gpuStreams[gpu]
                    // ); CUERR;                

                    cudaEventRecord(ws.forkStreamEvent, gpuStreams[gpu]); CUERR;

                    cudaSetDevice(masterDeviceId);
                    cudaStreamWaitEvent(masterStream1, ws.forkStreamEvent, 0); CUERR;
                }
            }else{
                //processQueryOnGpusWithTargetSubjectIds currently does not utilize multiple gpus

                for(int gpu = 0; gpu < 1; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];

                    cudaMemcpyAsync(
                        d_finalAlignmentScores_allGpus.data() + results_per_query*gpu,
                        ws.d_topN_scores.data(),
                        sizeof(float) * results_per_query,
                        cudaMemcpyDeviceToDevice,
                        gpuStreams[gpu]
                    ); CUERR;
                    cudaMemcpyAsync(
                        d_finalReferenceIds_allGpus.data() + results_per_query*gpu,
                        ws.d_topN_refIds.data(),
                        sizeof(ReferenceIdT) * results_per_query,
                        cudaMemcpyDeviceToDevice,
                        gpuStreams[gpu]
                    ); CUERR;
                    cudaMemcpyAsync(
                        d_finalEndPositions_allGpus.data() + results_per_query*gpu,
                        ws.d_topN_alignmentEndPositions.data(),
                        sizeof(AlignmentEndPosition) * results_per_query,
                        cudaMemcpyDeviceToDevice,
                        gpuStreams[gpu]
                    ); CUERR;  
                    // cudaMemcpyAsync(
                    //     d_resultNumOverflows.data() + gpu,
                    //     ws.d_total_overflow_number.data(),
                    //     sizeof(int),
                    //     cudaMemcpyDeviceToDevice,
                    //     gpuStreams[gpu]
                    // ); CUERR;                

                    cudaEventRecord(ws.forkStreamEvent, gpuStreams[gpu]); CUERR;

                    cudaSetDevice(masterDeviceId);
                    cudaStreamWaitEvent(masterStream1, ws.forkStreamEvent, 0); CUERR;
                }
            }

            cudaSetDevice(masterDeviceId);

            if(!targetSubjectIds){
                if(numGpus > 1){
                    //sort per-gpu top results to find overall top results
                    auto sortInputKeys = thrust::make_zip_iterator(
                        d_finalAlignmentScores_allGpus.begin(),
                        d_finalReferenceIds_allGpus.begin()
                    );
                    thrust::sort_by_key(
                        thrust::cuda::par_nosync(thrust_async_allocator<char>(masterStream1)).on(masterStream1),
                        sortInputKeys,
                        sortInputKeys + results_per_query * numGpus,
                        d_finalEndPositions_allGpus.begin(),
                        CompareScoresDescendingRefIdsAscending{}
                    );

                    // thrust::sort_by_key(
                    //     thrust::cuda::par_nosync(thrust_async_allocator<char>(masterStream1)).on(masterStream1),
                    //     d_finalAlignmentScores_allGpus.begin(),
                    //     d_finalAlignmentScores_allGpus.begin() + results_per_query * numGpus,
                    //     thrust::make_zip_iterator(
                    //         d_finalReferenceIds_allGpus.begin(),
                    //         d_finalEndPositions_allGpus.begin()
                    //     ),
                    //     thrust::greater<float>()
                    // );


                    //sum the overflows per gpu
                    //sumNumOverflowsKernel<<<1,1,0,masterStream1>>>(d_resultNumOverflows.data(), d_resultNumOverflows.data(), numGpus); CUERR;                
                }
            }

            cudaMemcpyAsync(
                h_finalAlignmentScores.data(), 
                d_finalAlignmentScores_allGpus.data(), 
                sizeof(float) * results_per_query, 
                cudaMemcpyDeviceToHost, 
                masterStream1
            );  CUERR
            cudaMemcpyAsync(
                h_finalReferenceIds.data(), 
                d_finalReferenceIds_allGpus.data(), 
                sizeof(ReferenceIdT) * results_per_query, 
                cudaMemcpyDeviceToHost, 
                masterStream1
            );  CUERR
            cudaMemcpyAsync(
                h_finalEndPositions.data(), 
                d_finalEndPositions_allGpus.data(), 
                sizeof(AlignmentEndPosition) * results_per_query, 
                cudaMemcpyDeviceToHost, 
                masterStream1
            );  CUERR
            // cudaMemcpyAsync(
            //     resultNumOverflows.data(), 
            //     d_resultNumOverflows.data(), 
            //     sizeof(int), 
            //     cudaMemcpyDeviceToHost, 
            //     masterStream1
            // );  CUERR

            cudaStreamSynchronize(masterStream1); CUERR;

            if(targetSubjectIds){
                //h_finalReferenceIds will contain numbers from 0 to num target subject ids. convert to proper target subject ids
                for(int i = 0; i < results_per_query; i++){
                    h_finalReferenceIds[i] = targetSubjectIds->subjectIds[h_finalReferenceIds[i]];
                }
            }
        }

        void processQueryOnGpus(){

            // std::cout << "ProcessQueryOnGpus: dstinfos isUploaded\n";
            // for(size_t i = 0; i < batchPlans[0].size(); i++){
            //     std::cout << batchPlansDstInfoVec[0][i].isUploaded << " ";
            // }
            // std::cout << "\n";

            const std::vector<std::vector<DBdataView>>& dbPartitionsPerGpu = subPartitionsForGpus;

            // constexpr auto boundaries = getLengthPartitionBoundaries();
            // constexpr int numLengthPartitions = boundaries.size();
            const int numGpus = deviceIds.size();
            const bool useExtraThreadForBatchTransfer = numGpus > 1;
        
            size_t totalNumberOfSequencesToProcess = std::accumulate(numSequencesPerGpu.begin(), numSequencesPerGpu.end(), 0u);
            
            size_t totalNumberOfProcessedSequences = 0;
        
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                auto& ws = *workingSets[gpu];
        
                //cudaMemsetAsync(ws.d_total_overflow_number.data(), 0, sizeof(int), gpuStreams[gpu]);
                
                ws.resetMaxReduceArray(gpuStreams[gpu]);
                ws.resetTopNArrays(gpuStreams[gpu]);
        
                //create dependency on mainStream
                cudaEventRecord(ws.forkStreamEvent, gpuStreams[gpu]); CUERR;
                cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.forkStreamEvent, 0); CUERR;
                for(auto& stream : ws.workStreamsWithoutTemp){
                    cudaStreamWaitEvent(stream, ws.forkStreamEvent, 0); CUERR;
                }
                cudaStreamWaitEvent(ws.hostFuncStream, ws.forkStreamEvent, 0); CUERR;
            }       
        
            //variables per gpu to keep between loops
            struct Variables{
                int currentBuffer = 0;
                int previousBuffer = 0;
                cudaStream_t H2DcopyStream = cudaStreamLegacy;
                char* h_inputChars = nullptr;
                SequenceLengthT* h_inputLengths = nullptr;
                size_t* h_inputOffsets = nullptr;
                char* d_inputChars = nullptr;
                SequenceLengthT* d_inputLengths = nullptr;
                size_t* d_inputOffsets = nullptr;
                //int* d_overflow_number = nullptr;
                //ReferenceIdT* d_overflow_positions = nullptr;
                const std::vector<DeviceBatchCopyToPinnedPlan>* batchPlansPtr;
                const std::vector<DeviceBatchCopyToPinnedPlan>* batchPlansCachedDBPtr;
                const DeviceBatchCopyToPinnedPlan* currentPlanPtr;
                size_t processedSequences = 0;
                size_t processedBatches = 0;
            };
        
            std::vector<Variables> variables_vec(numGpus);
            //init variables
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                const auto& ws = *workingSets[gpu];
                auto& variables = variables_vec[gpu];
                variables.processedSequences = 0;
                variables.processedBatches = 0;
                variables.batchPlansPtr = &batchPlans[gpu];
                variables.batchPlansCachedDBPtr = &batchPlans_cachedDB[gpu];
            }
            
            while(totalNumberOfProcessedSequences < totalNumberOfSequencesToProcess){
                //set up gpu variables for current iteration
                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];
                    auto& variables = variables_vec[gpu];
                    if(variables.processedBatches < variables.batchPlansPtr->size()){

                        if(variables.processedBatches < ws.getNumBatchesInCachedDB()){
                            //will process a batch that could be cached in gpu memory
                            if(batchPlansDstInfoVec[gpu][variables.processedBatches].isUploaded == false){
                                //it is not cached, need upload
                                variables.currentBuffer = ws.copyBufferIndex;
                                if(variables.currentBuffer == 0){
                                    variables.previousBuffer = ws.numCopyBuffers - 1;
                                }else{
                                    variables.previousBuffer = (variables.currentBuffer - 1);
                                } 
                                variables.H2DcopyStream = ws.copyStreams[variables.currentBuffer];
                                if(ws.h_chardata_vec[variables.currentBuffer].size() > 0){
                                    variables.h_inputChars = ws.h_chardata_vec[variables.currentBuffer].data();
                                }else{
                                    variables.h_inputChars = nullptr;
                                }
                                if(ws.h_lengthdata_vec[variables.currentBuffer].size() > 0){
                                    variables.h_inputLengths = ws.h_lengthdata_vec[variables.currentBuffer].data();
                                }else{
                                    variables.h_inputLengths = nullptr;
                                }
                                if(ws.h_offsetdata_vec[variables.currentBuffer].size() > 0){
                                    variables.h_inputOffsets = ws.h_offsetdata_vec[variables.currentBuffer].data();
                                }else{
                                    variables.h_inputOffsets = nullptr;
                                }
                                variables.d_inputChars = batchPlansDstInfoVec[gpu][variables.processedBatches].charsPtr;
                                variables.d_inputLengths = batchPlansDstInfoVec[gpu][variables.processedBatches].lengthsPtr;
                                variables.d_inputOffsets = batchPlansDstInfoVec[gpu][variables.processedBatches].offsetsPtr;
                                //variables.d_overflow_number = ws.d_overflow_number.data() + variables.currentBuffer;
                                //variables.d_overflow_positions = ws.d_overflow_positions_vec[variables.currentBuffer].data();
                            }else{
                                //already uploaded. process all batches for cached db together
                                assert(variables.processedBatches == 0);
                                variables.currentBuffer = 0;
                                variables.previousBuffer = 0;
                                variables.H2DcopyStream = ws.copyStreams[0];
                                variables.h_inputChars = nullptr;
                                variables.h_inputLengths = nullptr;
                                variables.h_inputOffsets = nullptr;
                                variables.d_inputChars = ws.d_cacheddb->getCharData();
                                variables.d_inputLengths = ws.d_cacheddb->getLengthData();
                                variables.d_inputOffsets = ws.d_cacheddb->getOffsetData();
                                
                            }
                        }else{
                            //will process batch that cannot be cached
                            //upload to double buffer
                            variables.currentBuffer = ws.copyBufferIndex;
                            if(variables.currentBuffer == 0){
                                variables.previousBuffer = ws.numCopyBuffers - 1;
                            }else{
                                variables.previousBuffer = (variables.currentBuffer - 1);
                            } 
                            variables.H2DcopyStream = ws.copyStreams[variables.currentBuffer];
                            if(ws.h_chardata_vec[variables.currentBuffer].size() > 0){
                                variables.h_inputChars = ws.h_chardata_vec[variables.currentBuffer].data();
                            }else{
                                variables.h_inputChars = nullptr;
                            }
                            if(ws.h_lengthdata_vec[variables.currentBuffer].size() > 0){
                                variables.h_inputLengths = ws.h_lengthdata_vec[variables.currentBuffer].data();
                            }else{
                                variables.h_inputLengths = nullptr;
                            }
                            if(ws.h_offsetdata_vec[variables.currentBuffer].size() > 0){
                                variables.h_inputOffsets = ws.h_offsetdata_vec[variables.currentBuffer].data();
                            }else{
                                variables.h_inputOffsets = nullptr;
                            }
                            variables.d_inputChars = batchPlansDstInfoVec[gpu][variables.processedBatches].charsPtr;
                            variables.d_inputLengths = batchPlansDstInfoVec[gpu][variables.processedBatches].lengthsPtr;
                            variables.d_inputOffsets = batchPlansDstInfoVec[gpu][variables.processedBatches].offsetsPtr;
                            //variables.d_overflow_number = ws.d_overflow_number.data() + variables.currentBuffer;
                            //variables.d_overflow_positions = ws.d_overflow_positions_vec[variables.currentBuffer].data();
                        }
                    }
                }
                //upload batch
                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];
                    auto& variables = variables_vec[gpu];
                    if(variables.processedBatches < variables.batchPlansPtr->size()){
                        const bool needsUpload = !batchPlansDstInfoVec[gpu][variables.processedBatches].isUploaded;

                        variables.currentPlanPtr = [&](){
                            if(variables.processedBatches < ws.getNumBatchesInCachedDB()){
                                if(!needsUpload){
                                    return &(*variables.batchPlansCachedDBPtr)[0];
                                    //return &(*variables.batchPlansPtr)[variables.processedBatches];
                                }else{
                                    return &(*variables.batchPlansPtr)[variables.processedBatches];
                                }
                            }else{
                                return &(*variables.batchPlansPtr)[variables.processedBatches];
                            }
                        }();
                            
        
                        if(needsUpload){
                            //transfer data
                            //can only overwrite device buffer if it is no longer in use on workstream
                            cudaStreamWaitEvent(variables.H2DcopyStream, ws.deviceBufferEvents[variables.currentBuffer], 0); CUERR;
        
                            if(useExtraThreadForBatchTransfer){
                                assert(variables.h_inputChars != nullptr);
                                assert(variables.h_inputLengths != nullptr);
                                assert(variables.h_inputOffsets != nullptr);

                                cudaStreamWaitEvent(ws.hostFuncStream, ws.pinnedBufferEvents[variables.currentBuffer]); CUERR;
                                executePinnedCopyPlanWithHostCallback(
                                    *variables.currentPlanPtr, 
                                    variables.h_inputChars,
                                    variables.h_inputLengths,
                                    variables.h_inputOffsets,
                                    dbPartitionsPerGpu[gpu], 
                                    ws.hostFuncStream
                                );
                                cudaEventRecord(ws.forkStreamEvent, ws.hostFuncStream); CUERR;
                                cudaStreamWaitEvent(variables.H2DcopyStream, ws.forkStreamEvent, 0);
        
                                cudaMemcpyAsync(
                                    variables.d_inputChars,
                                    variables.h_inputChars,
                                    variables.currentPlanPtr->usedBytes,
                                    H2D,
                                    variables.H2DcopyStream
                                ); CUERR;
                                cudaMemcpyAsync(
                                    variables.d_inputLengths,
                                    variables.h_inputLengths,
                                    sizeof(SequenceLengthT) * variables.currentPlanPtr->usedSeq,
                                    H2D,
                                    variables.H2DcopyStream
                                ); CUERR;
                                cudaMemcpyAsync(
                                    variables.d_inputOffsets,
                                    variables.h_inputOffsets,
                                    sizeof(size_t) * (variables.currentPlanPtr->usedSeq+1),
                                    H2D,
                                    variables.H2DcopyStream
                                ); CUERR;
                            }else{
                                //synchronize to avoid overwriting pinned buffer of target before it has been fully transferred
                                cudaEventSynchronize(ws.pinnedBufferEvents[variables.currentBuffer]); CUERR;

                                executeCopyPlanH2DDirect(
                                    *variables.currentPlanPtr, 
                                    variables.d_inputChars,
                                    variables.d_inputLengths,
                                    variables.d_inputOffsets,
                                    dbPartitionsPerGpu[gpu], 
                                    variables.H2DcopyStream
                                );

                                // assert(variables.h_inputChars != nullptr);
                                // assert(variables.h_inputLengths != nullptr);
                                // assert(variables.h_inputOffsets != nullptr);
        
                                // executePinnedCopyPlanSerialAndTransferToGpu(
                                //     *variables.currentPlanPtr, 
                                //     variables.h_inputChars,
                                //     variables.h_inputLengths,
                                //     variables.h_inputOffsets,
                                //     variables.d_inputChars,
                                //     variables.d_inputLengths,
                                //     variables.d_inputOffsets,
                                //     dbPartitionsPerGpu[gpu], 
                                //     variables.H2DcopyStream
                                // );
                            }
                            
                            cudaEventRecord(ws.pinnedBufferEvents[variables.currentBuffer], variables.H2DcopyStream); CUERR;
                        }
                    }
                }

                //all data is ready for alignments. create dependencies for work streams
                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];
                    auto& variables = variables_vec[gpu];

                    if(variables.processedBatches < variables.batchPlansPtr->size()){
                        
                        cudaEventRecord(ws.forkStreamEvent, variables.H2DcopyStream); CUERR;
                        cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.forkStreamEvent, 0); CUERR;
                        for(auto& stream : ws.workStreamsWithoutTemp){
                            cudaStreamWaitEvent(stream, ws.forkStreamEvent, 0); CUERR;
                        }
                        //wait for previous batch to finish
                        cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.deviceBufferEvents[variables.previousBuffer], 0); CUERR;
                        for(auto& stream : ws.workStreamsWithoutTemp){
                            cudaStreamWaitEvent(stream, ws.deviceBufferEvents[variables.previousBuffer], 0); CUERR;
                        }

                    }
                }

                //determine maximum number of sequences to process over all gpus
                size_t maxNumSequencesInBatchForGpus = 0;
                for(int gpu = 0; gpu < numGpus; gpu++){
                    auto& variables = variables_vec[gpu];
                    if(variables.processedBatches < variables.batchPlansPtr->size()){
                        maxNumSequencesInBatchForGpus = std::max(maxNumSequencesInBatchForGpus, variables.currentPlanPtr->usedSeq);
                    }
                }
                const size_t seqsPerPass = maxReduceArraySize;

                for(size_t sequencePassOffset = 0; sequencePassOffset < maxNumSequencesInBatchForGpus; sequencePassOffset += seqsPerPass){
                    for(int gpu = 0; gpu < numGpus; gpu++){
                        cudaSetDevice(deviceIds[gpu]); CUERR;
                        auto& ws = *workingSets[gpu];
                        auto& variables = variables_vec[gpu];
    
                        if(variables.processedBatches < variables.batchPlansPtr->size()){
                            const size_t numSequencesInBatch = variables.currentPlanPtr->usedSeq;
                            
                            if(sequencePassOffset < numSequencesInBatch){
                                const char* const inputChars = variables.d_inputChars;
                                const SequenceLengthT* const inputLengths = variables.d_inputLengths;
                                const size_t* const inputOffsets = variables.d_inputOffsets;
                                auto d_selectedPositions = thrust::make_counting_iterator<ReferenceIdT>(sequencePassOffset);
                                const size_t numInPass = std::min(numSequencesInBatch - sequencePassOffset, seqsPerPass);
                                const cudaStream_t stream = ws.workStreamsWithoutTemp[0];

                                if(scanType == ScanType::Gapless){
                                    auto maxReduceArray = ws.getMaxReduceArray(variables.processedSequences + sequencePassOffset);

                                    runGaplessFilterKernels_PSSM(
                                        maxReduceArray,
                                        ws.gpuPermutedPSSMforGapless,
                                        inputChars,
                                        inputLengths,
                                        inputOffsets,
                                        d_selectedPositions,
                                        numInPass,
                                        ws.d_tempStorageHE.data(),
                                        ws.numTempBytes,
                                        stream
                                    );

                                    //db sequences are processed in ascending order. stable sort ensures that sequences with same score are sorted by ascending id without a custom comparator
                                    thrust::stable_sort_by_key(
                                        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                        ws.d_maxReduceArrayScores.data(),
                                        ws.d_maxReduceArrayScores.data() + numInPass,
                                        ws.d_maxReduceArrayIndices.data(),
                                        thrust::greater<float>()
                                    );

                                    if(sequencePassOffset > 0 || totalNumberOfProcessedSequences > 0){
                                        auto mergeInput1 = thrust::make_zip_iterator(
                                            ws.d_maxReduceArrayScores.data(),
                                            ws.d_maxReduceArrayIndices.data()
                                        );
                                        auto mergeInput2 = thrust::make_zip_iterator(
                                            ws.d_topN_scores.data(), 
                                            ws.d_topN_refIds.data()
                                        );
                                        auto mergeOutput = thrust::make_zip_iterator(
                                            ws.d_topN_scores_tmp.data(), 
                                            ws.d_topN_refIds_tmp.data()
                                        );
                                        thrust::merge(
                                            thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                            mergeInput1,
                                            mergeInput1 + std::min(numInPass, size_t(results_per_query)),
                                            mergeInput2,
                                            mergeInput2 + results_per_query,
                                            mergeOutput,
                                            CompareScoresDescendingRefIdsAscending{}
                                        );

                                        std::swap(ws.d_topN_scores, ws.d_topN_scores_tmp);
                                        std::swap(ws.d_topN_refIds, ws.d_topN_refIds_tmp);
                                    }else{
                                        cudaMemcpyAsync(
                                            ws.d_topN_scores.data(),
                                            ws.d_maxReduceArrayScores.data(), 
                                            sizeof(float) * results_per_query,
                                            cudaMemcpyDeviceToDevice,
                                            stream
                                        ); CUERR;
                                        cudaMemcpyAsync(
                                            ws.d_topN_refIds.data(),
                                            ws.d_maxReduceArrayIndices.data(), 
                                            sizeof(ReferenceIdT) * results_per_query,
                                            cudaMemcpyDeviceToDevice,
                                            stream
                                        ); CUERR;
                                    }
                                }else if(scanType == ScanType::SW_Endpos){
                                    constexpr bool subjectIsCaseSensitive = true;
                                    constexpr bool withEndPosition = true;

                                    auto maxReduceArray = ws.getMaxReduceArrayWithEndPositions(variables.processedSequences + sequencePassOffset);
        
                                    run_SW_endposition_kernels_PSSM<subjectIsCaseSensitive,withEndPosition>(
                                        maxReduceArray,
                                        ws.gpuPermutedPSSMforSW,
                                        inputChars,
                                        inputLengths,
                                        inputOffsets,
                                        d_selectedPositions,
                                        numInPass,
                                        ws.d_tempStorageHE.data(),
                                        ws.numTempBytes,
                                        stream
                                    );

                                    thrust::stable_sort_by_key(
                                        thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                        ws.d_maxReduceArrayScores.data(),
                                        ws.d_maxReduceArrayScores.data() + numInPass,
                                        thrust::make_zip_iterator(
                                            ws.d_maxReduceArrayIndices.data(),
                                            ws.d_maxReduceArrayExtras.data()
                                        ),
                                        thrust::greater<float>()
                                    );

                                    if(sequencePassOffset > 0 || totalNumberOfProcessedSequences > 0){
                                        auto mergeInput1 = thrust::make_zip_iterator(
                                            ws.d_maxReduceArrayScores.data(),
                                            ws.d_maxReduceArrayIndices.data(),
                                            ws.d_maxReduceArrayExtras.data()
                                        );
                                        auto mergeInput2 = thrust::make_zip_iterator(
                                            ws.d_topN_scores.data(), 
                                            ws.d_topN_refIds.data(),
                                            ws.d_topN_alignmentEndPositions.data()
                                        );
                                        auto mergeOutput = thrust::make_zip_iterator(
                                            ws.d_topN_scores_tmp.data(), 
                                            ws.d_topN_refIds_tmp.data(),
                                            ws.d_topN_alignmentEndPositions_tmp.data()
                                        );
                                        thrust::merge(
                                            thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                            mergeInput1,
                                            mergeInput1 + std::min(numInPass, size_t(results_per_query)),
                                            mergeInput2,
                                            mergeInput2 + results_per_query,
                                            mergeOutput,
                                            CompareScoresDescendingRefIdsAscending{}
                                        );

                                        // thrust::merge_by_key(thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                        //     ws.d_maxReduceArrayScores.data(), 
                                        //     ws.d_maxReduceArrayScores.data() + std::min(numInPass, size_t(results_per_query)),
                                        //     ws.d_topN_scores.data(), 
                                        //     ws.d_topN_scores.data() + results_per_query,
                                        //     thrust::make_zip_iterator(
                                        //         ws.d_maxReduceArrayIndices.data(),
                                        //         ws.d_maxReduceArrayExtras.data()
                                        //     ),
                                        //     thrust::make_zip_iterator(
                                        //         ws.d_topN_refIds.data(),
                                        //         ws.d_topN_alignmentEndPositions.data()
                                        //     ),
                                        //     ws.d_topN_scores_tmp.data(), 
                                        //     thrust::make_zip_iterator(
                                        //         ws.d_topN_refIds_tmp.data(),
                                        //         ws.d_topN_alignmentEndPositions_tmp.data()
                                        //     ),
                                        //     thrust::greater<float>()
                                        // );

                                        std::swap(ws.d_topN_scores, ws.d_topN_scores_tmp);
                                        std::swap(ws.d_topN_refIds, ws.d_topN_refIds_tmp);
                                        std::swap(ws.d_topN_alignmentEndPositions, ws.d_topN_alignmentEndPositions_tmp);
                                    }else{
                                        cudaMemcpyAsync(
                                            ws.d_topN_scores.data(),
                                            ws.d_maxReduceArrayScores.data(), 
                                            sizeof(float) * results_per_query,
                                            cudaMemcpyDeviceToDevice,
                                            stream
                                        ); CUERR;
                                        cudaMemcpyAsync(
                                            ws.d_topN_refIds.data(),
                                            ws.d_maxReduceArrayIndices.data(), 
                                            sizeof(ReferenceIdT) * results_per_query,
                                            cudaMemcpyDeviceToDevice,
                                            stream
                                        ); CUERR;
                                        cudaMemcpyAsync(
                                            ws.d_topN_alignmentEndPositions.data(),
                                            ws.d_maxReduceArrayExtras.data(), 
                                            sizeof(AlignmentEndPosition) * results_per_query,
                                            cudaMemcpyDeviceToDevice,
                                            stream
                                        ); CUERR;
                                    }
                                }
                            }
                        }
                    }
                }

                //alignments are done in workstreams. now, join all workstreams
                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];
                    auto& variables = variables_vec[gpu];

                    if(variables.processedBatches < variables.batchPlansPtr->size()){        
                        for(auto& stream : ws.workStreamsWithoutTemp){
                            cudaEventRecord(ws.forkStreamEvent, stream); CUERR;
                            cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.forkStreamEvent, 0); CUERR;    
                        }
                    }
                }
        
                //finish processing of batch
                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];
                    const auto& variables = variables_vec[gpu];
                    if(variables.processedBatches < variables.batchPlansPtr->size()){
                
                        //the batch is done and its data can be resused
                        cudaEventRecord(ws.deviceBufferEvents[variables.currentBuffer], ws.workStreamForTempUsage); CUERR;
        
                        //let other workstreams depend on temp usage stream
                        for(auto& stream : ws.workStreamsWithoutTemp){
                            cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.deviceBufferEvents[variables.currentBuffer], 0); CUERR;    
                        }
        
                        ws.copyBufferIndex = (ws.copyBufferIndex+1) % ws.numCopyBuffers;
                    }
                }
        
                //update running numbers
                for(int gpu = 0; gpu < numGpus; gpu++){
                    auto& variables = variables_vec[gpu];
                    if(variables.processedBatches < variables.batchPlansPtr->size()){

                        variables.processedSequences += variables.currentPlanPtr->usedSeq;
                        if(batchPlansDstInfoVec[gpu][variables.processedBatches].isUploaded){
                            variables.processedBatches += workingSets[gpu]->getNumBatchesInCachedDB();                            
                            //variables.processedBatches++;
                        }else{
                            variables.processedBatches++;
                        }
                        //std::cout << "variables.processedBatches: " << variables.processedBatches << "\n";
        
                        totalNumberOfProcessedSequences += variables.currentPlanPtr->usedSeq;
                    } 
                }
        
            } //while not done
        
        
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                auto& ws = *workingSets[gpu];

                if(batchPlansDstInfoVec[gpu].size() > 0){
                    if(!batchPlansDstInfoVec[gpu][0].isUploaded){
                        //all batches for cached db are now resident in gpu memory. update the flags
                        if(ws.getNumBatchesInCachedDB() > 0){
                            markCachedDBBatchesAsUploaded(gpu);

                            // current offsets in cached db store the offsets for each batch, i.e. for each batch the offsets will start again at 0
                            // compute prefix sum to obtain the single-batch offsets
                
                            cudaMemsetAsync(ws.d_cacheddb->getOffsetData(), 0, sizeof(size_t), ws.workStreamForTempUsage); CUERR;
                
                            auto d_paddedLengths = thrust::make_transform_iterator(
                                ws.d_cacheddb->getLengthData(),
                                RoundToNextMultiple<size_t, 4>{}
                            );
                
                            thrust::inclusive_scan(
                                thrust::cuda::par_nosync(thrust_async_allocator<char>(ws.workStreamForTempUsage)).on(ws.workStreamForTempUsage),
                                d_paddedLengths,
                                d_paddedLengths + ws.getNumSequencesInCachedDB(),
                                ws.d_cacheddb->getOffsetData() + 1
                            );
                        }
                    }
                }
            }
        
        
        
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                auto& ws = *workingSets[gpu];
                //create dependency for gpuStreams[gpu]
                cudaEventRecord(ws.forkStreamEvent, ws.workStreamForTempUsage); CUERR;
                cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
        
                for(auto& stream : ws.workStreamsWithoutTemp){
                    cudaEventRecord(ws.forkStreamEvent, stream); CUERR;
                    cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
                }
        
                // for(auto& stream : ws.copyStreams){
                //     cudaEventRecord(ws.forkStreamEvent, stream); CUERR;
                //     cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
                // }
        
                cudaEventRecord(ws.forkStreamEvent, ws.hostFuncStream); CUERR;
                cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
            }
            
            processingTheFirstQuery = false;
        }

        void processQueryOnGpusWithTargetSubjectIds(){


            assert(targetSubjectIds);
            const int numGpus = deviceIds.size();

            //todo: for proper multi-gpu we need to find the correct gpu for each target subject to re-use the cached data
            //for now, multi-gpu will always gather the data on the host and transfer it to gpu 0
            
            prefetchDBToGpus();
            cudaSetDevice(deviceIds[0]); CUERR;
            workingSets[0]->resetMaxReduceArray(gpuStreams[0]);
            workingSets[0]->resetTopNArrays(gpuStreams[0]);
            const size_t numCachedSubjects = workingSets[0]->d_cacheddb->getNumSubjects();
            auto cachedTargetSubjectIdsEnd = targetSubjectIds->begin();
            
            if(numGpus == 1){
                cachedTargetSubjectIdsEnd = std::lower_bound(targetSubjectIds->begin(), targetSubjectIds->end(), numCachedSubjects);
            }

            // [targetSubjectIds->begin() , cachedTargetSubjectIdsEnd) are in gpu mem and can be accessed directly via index, 
            // [cachedTargetSubjectIdsEnd, targetSubjectIds->end()) are in cpu mem and need to be processed in batches

            const size_t numCachedTargetSubjects = std::distance(targetSubjectIds->begin(), cachedTargetSubjectIdsEnd);
            const size_t numUncachedTargetSubjects = std::distance(cachedTargetSubjectIdsEnd, targetSubjectIds->end());
            if(verbose){
                std::cout << "numCachedTargetSubjects " << numCachedTargetSubjects << "\n";
                std::cout << "numUncachedTargetSubjects " << numUncachedTargetSubjects << "\n";
            }

            // size_t max300000 = 0;
            // for(size_t i = 0; i < 20000; i++){
            //     const auto& data = fullDB.getData();
            //     size_t index = fullDB.getData().numSequences() - 300000 + i;
            //     SequenceLengthT length = data.lengths()[index];
            //     max300000 += SDIV(length,4) * 4;
            // }
            // std::cout << "max300000 " << max300000 << "\n";


            //process subjects which reside in gpu memory
            if(numCachedTargetSubjects > 0){
                //cudaStream_t stream = ws.workStreamsWithoutTemp[0];
                cudaStream_t stream = gpuStreams[0];

                auto& ws = *workingSets[0];
                const char* const inputChars = ws.d_cacheddb->getCharData();
                const SequenceLengthT* const inputLengths = ws.d_cacheddb->getLengthData();
                const size_t* const inputOffsets = ws.d_cacheddb->getOffsetData();

                ReferenceIdT* d_selectedPositions;
                char* d_availableTempStorage = ws.d_tempStorageHE.data();
                size_t availableTempStorageBytes = ws.numTempBytes;

                size_t numBytesForSelectedPositions = SDIV(sizeof(ReferenceIdT) * numCachedTargetSubjects, 512) * 512;
                if(numBytesForSelectedPositions < availableTempStorageBytes * 0.4){
                    d_selectedPositions = (ReferenceIdT*)d_availableTempStorage;
                    d_availableTempStorage = ((char*)d_availableTempStorage) + numBytesForSelectedPositions;
                    availableTempStorageBytes -= numBytesForSelectedPositions;
                }else{
                    cudaMallocAsync(&d_selectedPositions, sizeof(ReferenceIdT) * numCachedTargetSubjects, stream); CUERR;
                }

                cudaMemcpyAsync(
                    d_selectedPositions, 
                    targetSubjectIds->subjectIds.data(), 
                    sizeof(ReferenceIdT) * numCachedTargetSubjects, 
                    cudaMemcpyHostToDevice, 
                    stream
                ); CUERR;

                if(scanType == ScanType::Gapless){
                    const size_t seqsPerPass = maxReduceArraySize;
                    for(size_t sequencePassOffset = 0; sequencePassOffset < numCachedTargetSubjects; sequencePassOffset += seqsPerPass){
                        const size_t numInPass = std::min(numCachedTargetSubjects - sequencePassOffset, seqsPerPass);

                        auto maxReduceArray = ws.getMaxReduceArray(sequencePassOffset);

                        runGaplessFilterKernels_PSSM(
                            maxReduceArray,
                            ws.gpuPermutedPSSMforGapless,
                            inputChars,
                            inputLengths,
                            inputOffsets,
                            d_selectedPositions + sequencePassOffset,
                            numInPass,
                            d_availableTempStorage,
                            availableTempStorageBytes,
                            stream
                        );

                        thrust::stable_sort_by_key(
                            thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                            ws.d_maxReduceArrayScores.data(),
                            ws.d_maxReduceArrayScores.data() + numInPass,
                            ws.d_maxReduceArrayIndices.data(),
                            thrust::greater<float>()
                        );

                        if(sequencePassOffset > 0){
                            auto mergeInput1 = thrust::make_zip_iterator(
                                ws.d_maxReduceArrayScores.data(),
                                ws.d_maxReduceArrayIndices.data()
                            );
                            auto mergeInput2 = thrust::make_zip_iterator(
                                ws.d_topN_scores.data(), 
                                ws.d_topN_refIds.data()
                            );
                            auto mergeOutput = thrust::make_zip_iterator(
                                ws.d_topN_scores_tmp.data(), 
                                ws.d_topN_refIds_tmp.data()
                            );
                            thrust::merge(
                                thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                mergeInput1,
                                mergeInput1 + std::min(numInPass, size_t(results_per_query)),
                                mergeInput2,
                                mergeInput2 + results_per_query,
                                mergeOutput,
                                CompareScoresDescendingRefIdsAscending{}
                            );
                            // thrust::merge_by_key(thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                            //     ws.d_maxReduceArrayScores.data(), 
                            //     ws.d_maxReduceArrayScores.data() + std::min(numInPass, size_t(results_per_query)),
                            //     ws.d_topN_scores.data(), 
                            //     ws.d_topN_scores.data() + results_per_query,
                            //     ws.d_maxReduceArrayIndices.data(), 
                            //     ws.d_topN_refIds.data(),
                            //     ws.d_topN_scores_tmp.data(), 
                            //     ws.d_topN_refIds_tmp.data(),
                            //     thrust::greater<float>()
                            // );

                            std::swap(ws.d_topN_scores, ws.d_topN_scores_tmp);
                            std::swap(ws.d_topN_refIds, ws.d_topN_refIds_tmp);
                        }else{
                            cudaMemcpyAsync(
                                ws.d_topN_scores.data(),
                                ws.d_maxReduceArrayScores.data(), 
                                sizeof(float) * results_per_query,
                                cudaMemcpyDeviceToDevice,
                                stream
                            ); CUERR;
                            cudaMemcpyAsync(
                                ws.d_topN_refIds.data(),
                                ws.d_maxReduceArrayIndices.data(), 
                                sizeof(ReferenceIdT) * results_per_query,
                                cudaMemcpyDeviceToDevice,
                                stream
                            ); CUERR;
                        }

                    }
                }else if(scanType == ScanType::SW_Endpos){
                    constexpr bool subjectIsCaseSensitive = true;
                    constexpr bool withEndPosition = true;                   

                    const size_t seqsPerPass = maxReduceArraySize;
                    for(size_t sequencePassOffset = 0; sequencePassOffset < numCachedTargetSubjects; sequencePassOffset += seqsPerPass){
                        const size_t numInPass = std::min(numCachedTargetSubjects - sequencePassOffset, seqsPerPass);

                        auto maxReduceArray = ws.getMaxReduceArrayWithEndPositions(sequencePassOffset);

                        run_SW_endposition_kernels_PSSM<subjectIsCaseSensitive,withEndPosition>(
                            maxReduceArray,
                            ws.gpuPermutedPSSMforSW,
                            inputChars,
                            inputLengths,
                            inputOffsets,
                            d_selectedPositions + sequencePassOffset,
                            numInPass,
                            d_availableTempStorage,
                            availableTempStorageBytes,
                            stream
                        );

                        thrust::stable_sort_by_key(
                            thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                            ws.d_maxReduceArrayScores.data(),
                            ws.d_maxReduceArrayScores.data() + numInPass,
                            thrust::make_zip_iterator(
                                ws.d_maxReduceArrayIndices.data(),
                                ws.d_maxReduceArrayExtras.data()
                            ),
                            thrust::greater<float>()
                        );

                        if(sequencePassOffset > 0){
                            auto mergeInput1 = thrust::make_zip_iterator(
                                ws.d_maxReduceArrayScores.data(),
                                ws.d_maxReduceArrayIndices.data(),
                                ws.d_maxReduceArrayExtras.data()
                            );
                            auto mergeInput2 = thrust::make_zip_iterator(
                                ws.d_topN_scores.data(), 
                                ws.d_topN_refIds.data(),
                                ws.d_topN_alignmentEndPositions.data()
                            );
                            auto mergeOutput = thrust::make_zip_iterator(
                                ws.d_topN_scores_tmp.data(), 
                                ws.d_topN_refIds_tmp.data(),
                                ws.d_topN_alignmentEndPositions_tmp.data()
                            );
                            thrust::merge(
                                thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                mergeInput1,
                                mergeInput1 + std::min(numInPass, size_t(results_per_query)),
                                mergeInput2,
                                mergeInput2 + results_per_query,
                                mergeOutput,
                                CompareScoresDescendingRefIdsAscending{}
                            );
                            // thrust::merge_by_key(thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                            //     ws.d_maxReduceArrayScores.data(), 
                            //     ws.d_maxReduceArrayScores.data() + std::min(numInPass, size_t(results_per_query)),
                            //     ws.d_topN_scores.data(), 
                            //     ws.d_topN_scores.data() + results_per_query,
                            //     thrust::make_zip_iterator(
                            //         ws.d_maxReduceArrayIndices.data(),
                            //         ws.d_maxReduceArrayExtras.data()
                            //     ),
                            //     thrust::make_zip_iterator(
                            //         ws.d_topN_refIds.data(),
                            //         ws.d_topN_alignmentEndPositions.data()
                            //     ),
                            //     ws.d_topN_scores_tmp.data(), 
                            //     thrust::make_zip_iterator(
                            //         ws.d_topN_refIds_tmp.data(),
                            //         ws.d_topN_alignmentEndPositions_tmp.data()
                            //     ),
                            //     thrust::greater<float>()
                            // );

                            std::swap(ws.d_topN_scores, ws.d_topN_scores_tmp);
                            std::swap(ws.d_topN_refIds, ws.d_topN_refIds_tmp);
                            std::swap(ws.d_topN_alignmentEndPositions, ws.d_topN_alignmentEndPositions_tmp);
                        }else{
                            cudaMemcpyAsync(
                                ws.d_topN_scores.data(),
                                ws.d_maxReduceArrayScores.data(), 
                                sizeof(float) * results_per_query,
                                cudaMemcpyDeviceToDevice,
                                stream
                            ); CUERR;
                            cudaMemcpyAsync(
                                ws.d_topN_refIds.data(),
                                ws.d_maxReduceArrayIndices.data(), 
                                sizeof(ReferenceIdT) * results_per_query,
                                cudaMemcpyDeviceToDevice,
                                stream
                            ); CUERR;
                            cudaMemcpyAsync(
                                ws.d_topN_alignmentEndPositions.data(),
                                ws.d_maxReduceArrayExtras.data(), 
                                sizeof(AlignmentEndPosition) * results_per_query,
                                cudaMemcpyDeviceToDevice,
                                stream
                            ); CUERR;
                        }

                    }
                }

                if(numBytesForSelectedPositions >= availableTempStorageBytes * 0.4){
                    cudaFreeAsync(d_selectedPositions, stream); CUERR;
                }
            }

            //process subjects which reside in host memory
            if(numUncachedTargetSubjects > 0){
                auto& ws = *workingSets[0];
                //cudaStream_t stream = ws.workStreamsWithoutTemp[0];
                cudaStream_t stream = gpuStreams[0];

                helpers::CpuTimer targetGatherTimer("targetGatherTimer");
                std::vector<char> targetChars;
                std::vector<size_t> targetOffsets(numUncachedTargetSubjects+1, 0);
                std::vector<SequenceLengthT> targetLengths(numUncachedTargetSubjects);
                for(size_t i = 0; i < numUncachedTargetSubjects; i++){
                    const auto& data = fullDB.getData();
                    const ReferenceIdT subjectId = *(cachedTargetSubjectIdsEnd + i);
                    const size_t offsetBegin = data.offsets()[subjectId];
                    const size_t offsetEnd = data.offsets()[subjectId+1];
                    SequenceLengthT length = data.lengths()[subjectId];
                    targetChars.insert(targetChars.end(), data.chars() + offsetBegin, data.chars() + offsetEnd);
                    targetOffsets[i+1] = targetChars.size();
                    targetLengths[i] = length;
                }
                if(verbose){
                    targetGatherTimer.print();
                }

                std::vector<DBdataView> targetDBPartition{DBdataView(
                    0,
                    numUncachedTargetSubjects,
                    0,
                    targetChars.data(),
                    targetLengths.data(),
                    targetOffsets.data(),
                    nullptr,
                    nullptr
                )};
                std::vector<DeviceBatchCopyToPinnedPlan> targetBatchPlans = computeDbCopyPlan(
                    targetDBPartition,
                    {0},
                    memoryConfig.maxBatchBytes,
                    memoryConfig.maxBatchSequences
                );


                char* d_targetBasePtr;
                char* d_targetChars;
                size_t* d_targetOffsets;
                SequenceLengthT* d_targetLengths;
                char* d_availableTempStorage = ws.d_tempStorageHE.data();
                size_t availableTempStorageBytes = ws.numTempBytes;
                
                size_t bytes[3]{
                    SDIV(std::min(memoryConfig.maxBatchBytes, targetChars.size()), 512) * 512, //chars
                    SDIV(sizeof(size_t) * std::min(memoryConfig.maxBatchSequences+1, targetLengths.size()+1), 512) * 512, //offsets
                    SDIV(sizeof(SequenceLengthT) * std::min(memoryConfig.maxBatchSequences, targetLengths.size()), 512) * 512, //lengths
                };
                size_t bytesSum = bytes[0] + bytes[1] + bytes[2];

                if(bytesSum < availableTempStorageBytes * 0.5){
                    d_targetChars = (char*)d_availableTempStorage;
                    d_availableTempStorage = ((char*)d_availableTempStorage) + bytes[0];                    
                    availableTempStorageBytes -= bytes[0];
                    d_targetOffsets = (size_t*)d_availableTempStorage;
                    d_availableTempStorage = ((char*)d_availableTempStorage) + bytes[1];                    
                    availableTempStorageBytes -= bytes[1];
                    d_targetLengths = (SequenceLengthT*)d_availableTempStorage;
                    d_availableTempStorage = ((char*)d_availableTempStorage) + bytes[2];                    
                    availableTempStorageBytes -= bytes[2];
                }else{
                    cudaMallocAsync(&d_targetBasePtr, bytesSum, stream); CUERR;
                    d_targetChars = (char*)d_targetBasePtr;
                    d_targetOffsets = (size_t*)(((char*)d_targetChars) + bytes[0]);
                    d_targetLengths = (SequenceLengthT*)(((char*)d_targetOffsets) + bytes[1]);
                }

                size_t numProcessed = 0;
                for(const auto& batchPlan : targetBatchPlans){

                    executeCopyPlanH2DDirect(
                        batchPlan, 
                        d_targetChars,
                        d_targetLengths,
                        d_targetOffsets,
                        targetDBPartition, 
                        stream
                    );
                
                    auto d_selectedPositions = thrust::make_counting_iterator<ReferenceIdT>(0);

                    if(scanType == ScanType::Gapless){

                        const size_t seqsPerPass = maxReduceArraySize;
                        for(size_t sequencePassOffset = 0; sequencePassOffset < numUncachedTargetSubjects; sequencePassOffset += seqsPerPass){
                            const size_t numInPass = std::min(numUncachedTargetSubjects - sequencePassOffset, seqsPerPass);

                            auto maxReduceArray = ws.getMaxReduceArray(numCachedTargetSubjects + numProcessed + sequencePassOffset);

                            runGaplessFilterKernels_PSSM(
                                maxReduceArray,
                                ws.gpuPermutedPSSMforGapless,
                                d_targetChars,
                                d_targetLengths,
                                d_targetOffsets,
                                d_selectedPositions + sequencePassOffset,
                                numInPass,
                                d_availableTempStorage,
                                availableTempStorageBytes,
                                stream
                            );

                            thrust::stable_sort_by_key(
                                thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                ws.d_maxReduceArrayScores.data(),
                                ws.d_maxReduceArrayScores.data() + numInPass,
                                ws.d_maxReduceArrayIndices.data(),
                                thrust::greater<float>()
                            );

                            //merge kernel results with previous results
                            if(sequencePassOffset > 0 || numCachedTargetSubjects > 0){
                                auto mergeInput1 = thrust::make_zip_iterator(
                                    ws.d_maxReduceArrayScores.data(),
                                    ws.d_maxReduceArrayIndices.data()
                                );
                                auto mergeInput2 = thrust::make_zip_iterator(
                                    ws.d_topN_scores.data(), 
                                    ws.d_topN_refIds.data()
                                );
                                auto mergeOutput = thrust::make_zip_iterator(
                                    ws.d_topN_scores_tmp.data(), 
                                    ws.d_topN_refIds_tmp.data()
                                );
                                thrust::merge(
                                    thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                    mergeInput1,
                                    mergeInput1 + std::min(numInPass, size_t(results_per_query)),
                                    mergeInput2,
                                    mergeInput2 + results_per_query,
                                    mergeOutput,
                                    CompareScoresDescendingRefIdsAscending{}
                                );
                                // thrust::merge_by_key(thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                //     ws.d_maxReduceArrayScores.data(), 
                                //     ws.d_maxReduceArrayScores.data() + std::min(numInPass, size_t(results_per_query)),
                                //     ws.d_topN_scores.data(), 
                                //     ws.d_topN_scores.data() + results_per_query,
                                //     ws.d_maxReduceArrayIndices.data(), 
                                //     ws.d_topN_refIds.data(),
                                //     ws.d_topN_scores_tmp.data(), 
                                //     ws.d_topN_refIds_tmp.data(),
                                //     thrust::greater<float>()
                                // );

                                std::swap(ws.d_topN_scores, ws.d_topN_scores_tmp);
                                std::swap(ws.d_topN_refIds, ws.d_topN_refIds_tmp);
                            }else{
                                cudaMemcpyAsync(
                                    ws.d_topN_scores.data(),
                                    ws.d_maxReduceArrayScores.data(), 
                                    sizeof(float) * results_per_query,
                                    cudaMemcpyDeviceToDevice,
                                    stream
                                ); CUERR;
                                cudaMemcpyAsync(
                                    ws.d_topN_refIds.data(),
                                    ws.d_maxReduceArrayIndices.data(), 
                                    sizeof(ReferenceIdT) * results_per_query,
                                    cudaMemcpyDeviceToDevice,
                                    stream
                                ); CUERR;
                            }

                        }
                    }else if(scanType == ScanType::SW_Endpos){
                        constexpr bool subjectIsCaseSensitive = true;
                        constexpr bool withEndPosition = true;

                        const size_t seqsPerPass = maxReduceArraySize;
                        for(size_t sequencePassOffset = 0; sequencePassOffset < numUncachedTargetSubjects; sequencePassOffset += seqsPerPass){
                            const size_t numInPass = std::min(numUncachedTargetSubjects - sequencePassOffset, seqsPerPass);

                            auto maxReduceArray = ws.getMaxReduceArrayWithEndPositions(numCachedTargetSubjects + numProcessed + sequencePassOffset);
    
                            run_SW_endposition_kernels_PSSM<subjectIsCaseSensitive,withEndPosition>(
                                maxReduceArray,
                                ws.gpuPermutedPSSMforSW,
                                d_targetChars,
                                d_targetLengths,
                                d_targetOffsets,
                                d_selectedPositions + sequencePassOffset,
                                numInPass,
                                d_availableTempStorage,
                                availableTempStorageBytes,
                                stream
                            );

                            thrust::stable_sort_by_key(
                                thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                ws.d_maxReduceArrayScores.data(),
                                ws.d_maxReduceArrayScores.data() + numInPass,
                                thrust::make_zip_iterator(
                                    ws.d_maxReduceArrayIndices.data(),
                                    ws.d_maxReduceArrayExtras.data()
                                ),
                                thrust::greater<float>()
                            );
                            
                            if(sequencePassOffset > 0 || numCachedTargetSubjects > 0){
                                auto mergeInput1 = thrust::make_zip_iterator(
                                    ws.d_maxReduceArrayScores.data(),
                                    ws.d_maxReduceArrayIndices.data(),
                                    ws.d_maxReduceArrayExtras.data()
                                );
                                auto mergeInput2 = thrust::make_zip_iterator(
                                    ws.d_topN_scores.data(), 
                                    ws.d_topN_refIds.data(),
                                    ws.d_topN_alignmentEndPositions.data()
                                );
                                auto mergeOutput = thrust::make_zip_iterator(
                                    ws.d_topN_scores_tmp.data(), 
                                    ws.d_topN_refIds_tmp.data(),
                                    ws.d_topN_alignmentEndPositions_tmp.data()
                                );
                                thrust::merge(
                                    thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                    mergeInput1,
                                    mergeInput1 + std::min(numInPass, size_t(results_per_query)),
                                    mergeInput2,
                                    mergeInput2 + results_per_query,
                                    mergeOutput,
                                    CompareScoresDescendingRefIdsAscending{}
                                );

                                // thrust::merge_by_key(thrust::cuda::par_nosync(thrust_async_allocator<char>(stream)).on(stream),
                                //     ws.d_maxReduceArrayScores.data(), 
                                //     ws.d_maxReduceArrayScores.data() + std::min(numInPass, size_t(results_per_query)),
                                //     ws.d_topN_scores.data(), 
                                //     ws.d_topN_scores.data() + results_per_query,
                                //     thrust::make_zip_iterator(
                                //         ws.d_maxReduceArrayIndices.data(),
                                //         ws.d_maxReduceArrayExtras.data()
                                //     ),
                                //     thrust::make_zip_iterator(
                                //         ws.d_topN_refIds.data(),
                                //         ws.d_topN_alignmentEndPositions.data()
                                //     ),
                                //     ws.d_topN_scores_tmp.data(), 
                                //     thrust::make_zip_iterator(
                                //         ws.d_topN_refIds_tmp.data(),
                                //         ws.d_topN_alignmentEndPositions_tmp.data()
                                //     ),
                                //     thrust::greater<float>()
                                // );

                                std::swap(ws.d_topN_scores, ws.d_topN_scores_tmp);
                                std::swap(ws.d_topN_refIds, ws.d_topN_refIds_tmp);
                                std::swap(ws.d_topN_alignmentEndPositions, ws.d_topN_alignmentEndPositions_tmp);
                            }else{
                                cudaMemcpyAsync(
                                    ws.d_topN_scores.data(),
                                    ws.d_maxReduceArrayScores.data(), 
                                    sizeof(float) * results_per_query,
                                    cudaMemcpyDeviceToDevice,
                                    stream
                                ); CUERR;
                                cudaMemcpyAsync(
                                    ws.d_topN_refIds.data(),
                                    ws.d_maxReduceArrayIndices.data(), 
                                    sizeof(ReferenceIdT) * results_per_query,
                                    cudaMemcpyDeviceToDevice,
                                    stream
                                ); CUERR;
                                cudaMemcpyAsync(
                                    ws.d_topN_alignmentEndPositions.data(),
                                    ws.d_maxReduceArrayExtras.data(), 
                                    sizeof(AlignmentEndPosition) * results_per_query,
                                    cudaMemcpyDeviceToDevice,
                                    stream
                                ); CUERR;
                            }

                        }
                    }

                    numProcessed += batchPlan.usedSeq;
                }

                if(bytesSum >= availableTempStorageBytes * 0.5){
                    cudaFreeAsync(d_targetBasePtr, stream); CUERR;
                }
            }

            // const int numGpus = deviceIds.size();
            // for(int gpu = 0; gpu < numGpus; gpu++){
            //     cudaSetDevice(deviceIds[gpu]); CUERR;
            //     auto& ws = *workingSets[gpu];
            //     //create dependency for gpuStreams[gpu]
            //     cudaEventRecord(ws.forkStreamEvent, ws.workStreamForTempUsage); CUERR;
            //     cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
        
            //     for(auto& stream : ws.workStreamsWithoutTemp){
            //         cudaEventRecord(ws.forkStreamEvent, stream); CUERR;
            //         cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
            //     }
        
            //     // for(auto& stream : ws.copyStreams){
            //     //     cudaEventRecord(ws.forkStreamEvent, stream); CUERR;
            //     //     cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
            //     // }
        
            //     cudaEventRecord(ws.forkStreamEvent, ws.hostFuncStream); CUERR;
            //     cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
            // }
            
            processingTheFirstQuery = false;
  
        }

        BenchmarkStats makeBenchmarkStats(double seconds, double cells, int overflows) const{
            BenchmarkStats stats;
            stats.seconds = seconds;
            stats.gcups = cells / 1000. / 1000. / 1000.;
            stats.gcups = stats.gcups / stats.seconds;
            stats.numOverflows = overflows;
            return stats;
        }

        void updateNumResultsPerQuery(){

            results_per_query = std::min(size_t(numTop), size_t(maxReduceArraySize));
            if(dbIsReady){
                results_per_query = std::min(size_t(results_per_query), fullDB.getData().numSequences());
            }
        }

        int getMaxSingleTileQueryLength_Gapless() const{
            auto largestconfig = *std::max_element(availableKernelConfigs_gapless_singletile.begin(),
                availableKernelConfigs_gapless_singletile.end(),
                [](const auto& l, const auto& r){
                    return l.tilesize < r.tilesize;
                }
            );  
            return largestconfig.tilesize;
            
            //must be a multiple of 64
            // return 2048;
            //return 0;
        }

        int getMaxSingleTileQueryLength_SW() const{
            auto largestconfig = *std::max_element(availableKernelConfigs_sw_singletile.begin(),
                availableKernelConfigs_sw_singletile.end(),
                [](const auto& l, const auto& r){
                    return l.tilesize < r.tilesize;
                }
            );  
            return largestconfig.tilesize;

            //must be a multiple of 32
            // return 1408;
            //return 0;
        }




        
    public:
        GaplessKernelConfig customKernelConfig_Gapless;
        bool useCustomKernelConfig_Gapless = false;

        SmithWatermanKernelConfig customKernelConfig_SW;
        bool useCustomKernelConfig_SW = false;

        void setCustomKernelConfig_Gapless(GaplessKernelConfig config){
            customKernelConfig_Gapless = config;
            useCustomKernelConfig_Gapless = true;
        }

        void setCustomKernelConfig_SW(SmithWatermanKernelConfig config){
            customKernelConfig_SW = config;
            useCustomKernelConfig_SW = true;
        }
    private:



        GaplessKernelConfig getSingleTileGroupRegConfigForPSSM_Gapless(int queryLength){
            if(useCustomKernelConfig_Gapless){
                return customKernelConfig_Gapless;
            }
            const auto& configs = availableKernelConfigs_gapless_singletile;
            auto it = std::lower_bound(configs.begin(), configs.end(), queryLength,
                [](const GaplessKernelConfig& l, int r){
                    return l.tilesize < r;
                }
            );
            if(it == configs.end()){
                throw std::runtime_error("kernel config does not exist");
            }else{
                return *it;
            }
        }

        GaplessKernelConfig getMultiTileGroupRegConfigForPSSM_Gapless(int queryLength){
            if(useCustomKernelConfig_Gapless){
                return customKernelConfig_Gapless;
            }
            const auto& configs = availableKernelConfigs_gapless_multitile;

            //find the config which best utilizes the last tile. larger tile sizes are preferred
            auto selectedConfig = configs[0];
            const int remainderInLastTile0 = queryLength % selectedConfig.tilesize;
            double utilization = remainderInLastTile0 == 0 ? 1.0 : double(remainderInLastTile0) / selectedConfig.tilesize;
            for(size_t i = 1; i < configs.size(); i++){
                const auto& newConfig = configs[i];
                const int remainderInLastTile = queryLength % newConfig.tilesize;
                const double newUtilization = remainderInLastTile == 0 ? 1.0 : double(remainderInLastTile) / newConfig.tilesize;
                if(newUtilization >= utilization){
                    utilization = newUtilization;
                    selectedConfig = newConfig;
                }
            }

            return selectedConfig;
        }

        SmithWatermanKernelConfig getSingleTileGroupRegConfigForPSSM_SW(int queryLength){
            if(useCustomKernelConfig_SW){
                return customKernelConfig_SW;
            }
            const auto& configs = availableKernelConfigs_sw_singletile;
            auto it = std::lower_bound(configs.begin(), configs.end(), queryLength,
                [](const SmithWatermanKernelConfig& l, int r){
                    return l.tilesize < r;
                }
            );
            if(it == configs.end()){
                throw std::runtime_error("kernel config does not exist");
            }else{
                return *it;
            }
        }

        SmithWatermanKernelConfig getMultiTileGroupRegConfigForPSSM_SW(int queryLength){
            if(useCustomKernelConfig_SW){
                return customKernelConfig_SW;
            }
            const auto& configs = availableKernelConfigs_sw_multitile;

            //find the config which best utilizes the last tile. larger tile sizes are preferred
            auto selectedConfig = configs[0];
            const int remainderInLastTile0 = queryLength % selectedConfig.tilesize;
            double utilization = remainderInLastTile0 == 0 ? 1.0 : double(remainderInLastTile0) / selectedConfig.tilesize;
            for(size_t i = 1; i < configs.size(); i++){
                const auto& newConfig = configs[i];
                const int remainderInLastTile = queryLength % newConfig.tilesize;
                const double newUtilization = remainderInLastTile == 0 ? 1.0 : double(remainderInLastTile) / newConfig.tilesize;
                if(newUtilization >= utilization){
                    utilization = newUtilization;
                    selectedConfig = newConfig;
                }
            }

            return selectedConfig;
        }

        std::vector<std::tuple<int,int>> getSupportedGroupRegConfigs_gapless_singletile() const{
            std::vector<std::tuple<int,int>> validRegConfigs;
            #define X(g,r)\
                validRegConfigs.push_back(std::make_tuple(g,r));
            
            PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

            #undef X
            return validRegConfigs;
        }

        std::vector<std::tuple<int,int>> getSupportedGroupRegConfigs_gapless_multitile() const{
            std::vector<std::tuple<int,int>> validRegConfigs;
            #define X(g,r)\
                validRegConfigs.push_back(std::make_tuple(g,r));
            
            PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

            #undef X
            return validRegConfigs;
        }

        std::vector<std::tuple<int,int>> getSupportedGroupRegConfigs_swendpos_singletile() const{
            std::vector<std::tuple<int,int>> validRegConfigs;
            #define X(g,r)\
                validRegConfigs.push_back(std::make_tuple(g,r));
            
            PSSM_SW_ENDPOS_SINGLETILE_FLOAT_OR_INT_FOR_EACH_VALID_CONFIG_DO_X

            #undef X
            return validRegConfigs;
        }

        std::vector<std::tuple<int,int>> getSupportedGroupRegConfigs_swendpos_multitile() const{
            std::vector<std::tuple<int,int>> validRegConfigs;
            #define X(g,r)\
                validRegConfigs.push_back(std::make_tuple(g,r));
            
            PSSM_SW_ENDPOS_MULTITILE_FLOAT_OR_INT_FOR_EACH_VALID_CONFIG_DO_X

            #undef X
            return validRegConfigs;
        }

        void initializeListOfAvailableKernelConfigs(const KernelConfigFilenames& kernelConfigFilenames){

            const auto configsGapless = [&](){
                if(kernelConfigFilenames.gapless){
                    return loadKernelConfigsFromFile_gapless(kernelConfigFilenames.gapless.value());
                }else{
                    return getOptimalKernelConfigs_gapless(deviceIds[0]);
                }
            }();

            {
                const auto supported = getSupportedGroupRegConfigs_gapless_singletile();

                for(const auto& config : configsGapless){
                    for(const auto& tup : supported){
                        if(config.groupsize == std::get<0>(tup) && config.numRegs == std::get<1>(tup)){
                            availableKernelConfigs_gapless_singletile.push_back(config);
                            break;
                        }
                    }
                }
                if(availableKernelConfigs_gapless_singletile.empty()){
                    throw std::runtime_error("availableKernelConfigs_gapless_singletile is empty");
                }
            }
            {
                const auto supported = getSupportedGroupRegConfigs_gapless_multitile();

                for(const auto& config : configsGapless){
                    for(const auto& tup : supported){
                        if(config.groupsize == std::get<0>(tup) && config.numRegs == std::get<1>(tup)){
                            availableKernelConfigs_gapless_multitile.push_back(config);
                            break;
                        }
                    }
                }
                if(availableKernelConfigs_gapless_multitile.empty()){
                    throw std::runtime_error("availableKernelConfigs_gapless_multitile is empty");
                }
            }

            const auto configsSW = [&](){
                if(kernelConfigFilenames.sw){
                    return loadKernelConfigsFromFile_sw(kernelConfigFilenames.sw.value());
                }else{
                    return getOptimalKernelConfigs_SW(deviceIds[0]);
                }
            }();

            {
                const auto supported = getSupportedGroupRegConfigs_swendpos_singletile();

                for(const auto& config : configsSW){
                    for(const auto& tup : supported){
                        if(config.groupsize == std::get<0>(tup) && config.numRegs == std::get<1>(tup)){
                            availableKernelConfigs_sw_singletile.push_back(config);
                            break;
                        }
                    }
                }
                if(availableKernelConfigs_sw_singletile.empty()){
                    throw std::runtime_error("availableKernelConfigs_sw_singletile is empty");
                }
            }
            {
                const auto supported = getSupportedGroupRegConfigs_swendpos_multitile();

                for(const auto& config : configsSW){
                    for(const auto& tup : supported){
                        if(config.groupsize == std::get<0>(tup) && config.numRegs == std::get<1>(tup)){
                            availableKernelConfigs_sw_multitile.push_back(config);
                            break;
                        }
                    }
                }
                if(availableKernelConfigs_sw_multitile.empty()){
                    throw std::runtime_error("availableKernelConfigs_sw_multitile is empty");
                }
            }
        }

        std::vector<GaplessKernelConfig> loadKernelConfigsFromFile_gapless(const std::string& filename){
            std::ifstream is(filename);
            if(!is){
                throw std::runtime_error("could not open file " + filename);
            }

            auto split = [](const std::string& str, char c){
                std::vector<std::string> result;

                std::stringstream ss(str);
                std::string s;

                while (std::getline(ss, s, c)) {
                    result.emplace_back(s);
                }

                return result;
            };

            std::vector<GaplessKernelConfig> result;

            std::string line;
            while(std::getline(is, line)){
                if(line.size() > 0){
                    if(line[0] == '#') continue;
                    auto tokens = split(line, ' ');
                    if(tokens.size() < 5) throw std::runtime_error("error parsing kernel configs file");

                    GaplessKernelConfig config;
                    config.tilesize = std::stoi(tokens[0]);
                    config.groupsize = std::stoi(tokens[1]);
                    config.numRegs = std::stoi(tokens[2]);
                    config.dpx = std::stoi(tokens[3]);
                    config.approach = GaplessKernelConfig::Approach(std::stoi(tokens[4]));
                    result.push_back(config);
                }
            }

            return result;
        }

        std::vector<SmithWatermanKernelConfig> loadKernelConfigsFromFile_sw(const std::string& filename){
            std::ifstream is(filename);
            if(!is){
                throw std::runtime_error("could not open file " + filename);
            }

            auto split = [](const std::string& str, char c){
                std::vector<std::string> result;

                std::stringstream ss(str);
                std::string s;

                while (std::getline(ss, s, c)) {
                    result.emplace_back(s);
                }

                return result;
            };

            std::vector<SmithWatermanKernelConfig> result;

            std::string line;
            while(std::getline(is, line)){
                if(line.size() > 0){
                    if(line[0] == '#') continue;
                    auto tokens = split(line, ' ');
                    if(tokens.size() < 5) throw std::runtime_error("error parsing kernel configs file");

                    SmithWatermanKernelConfig config;
                    config.tilesize = std::stoi(tokens[0]);
                    config.groupsize = std::stoi(tokens[1]);
                    config.numRegs = std::stoi(tokens[2]);
                    config.dpx = std::stoi(tokens[3]);
                    config.approach = SmithWatermanKernelConfig::Approach(std::stoi(tokens[4]));
                    result.push_back(config);
                }
            }

            return result;
        }

        int affine_local_DP_host_protein_blosum62(
            const char* seq1,
            const char* seq2,
            const int length1,
            const int length2,
            const int gap_open,
            const int gap_extend
        ) {
            const int NEGINFINITY = -10000;
            std::vector<int> penalty_H(2*(length2+1));
            std::vector<int> penalty_F(2*(length2+1));

            int E, F, maxi = 0, result;
            penalty_H[0] = 0;
            penalty_F[0] = NEGINFINITY;
            for (int index = 1; index <= length2; index++) {
                penalty_H[index] = 0;
                penalty_F[index] = NEGINFINITY;
            }

            auto convert_AA = cudasw4::ConvertAA_20{};

            auto BLOSUM = cudasw4::BLOSUM62_20::get2D();

            for (int row = 1; row <= length1; row++) {
                char seq1_char = seq1[row-1];
                char seq2_char;

                const int target_row = row & 1;
                const int source_row = !target_row;
                penalty_H[target_row*(length2+1)] = 0; //gap_open + (row-1)*gap_extend;
                penalty_F[target_row*(length2+1)] = gap_open + (row-1)*gap_extend;
                E = NEGINFINITY;
                for (int col = 1; col <= length2; col++) {
                    const int diag = penalty_H[source_row*(length2+1)+col-1];
                    const int abve = penalty_H[source_row*(length2+1)+col+0];
                    const int left = penalty_H[target_row*(length2+1)+col-1];
                    seq2_char = seq2[col-1];
                    const int residue = BLOSUM[convert_AA(seq1_char)][convert_AA(seq2_char)];
                    E = std::max(E+gap_extend, left+gap_open);
                    F = std::max(penalty_F[source_row*(length2+1)+col+0]+gap_extend, abve+gap_open);
                    result = std::max(0, std::max(diag + residue, std::max(E, F)));
                    penalty_H[target_row*(length2+1)+col] = result;
                    if (result > maxi) maxi = result;
                    penalty_F[target_row*(length2+1)+col] = F;
                }
            }
            return maxi;
        }

        //sequences must be in to ncbi converted format
        int affine_local_DP_host_protein_blosum62_converted(
            const char* seq1,
            const char* seq2,
            const int length1,
            const int length2,
            const int gap_open,
            const int gap_extend
        ) {
            const int NEGINFINITY = -10000;
            std::vector<int> penalty_H(2*(length2+1));
            std::vector<int> penalty_F(2*(length2+1));

            // std::cout << "length1 " << length1 << ", length2 " << length2 << "\n";

            // for(int i = 0; i < length1; i++){
            //     std::cout << int(seq1[i]) << " ";
            // }
            // std::cout << "\n";

            // for(int i = 0; i < length2; i++){
            //     std::cout << int(seq2[i]) << " ";
            // }
            // std::cout << "\n";

            int E, F, maxi = 0, result;
            penalty_H[0] = 0;
            penalty_F[0] = NEGINFINITY;
            for (int index = 1; index <= length2; index++) {
                penalty_H[index] = 0;
                penalty_F[index] = NEGINFINITY;
            }

            auto BLOSUM = cudasw4::BLOSUM62_20::get2D();

            for (int row = 1; row <= length1; row++) {
                int seq1_char = seq1[row-1];
                int seq2_char;

                const int target_row = row & 1;
                const int source_row = !target_row;
                penalty_H[target_row*(length2+1)] = 0; //gap_open + (row-1)*gap_extend;
                penalty_F[target_row*(length2+1)] = gap_open + (row-1)*gap_extend;
                E = NEGINFINITY;
                for (int col = 1; col <= length2; col++) {
                    const int diag = penalty_H[source_row*(length2+1)+col-1];
                    const int abve = penalty_H[source_row*(length2+1)+col+0];
                    const int left = penalty_H[target_row*(length2+1)+col-1];
                    seq2_char = seq2[col-1];
                    const int residue = BLOSUM[seq1_char][seq2_char];
                    E = std::max(E+gap_extend, left+gap_open);
                    F = std::max(penalty_F[source_row*(length2+1)+col+0]+gap_extend, abve+gap_open);
                    result = std::max(0, std::max(diag + residue, std::max(E, F)));
                    penalty_H[target_row*(length2+1)+col] = result;
                    if (result > maxi) maxi = result;
                    penalty_F[target_row*(length2+1)+col] = F;

                    //std::cout << maxi << " ";
                }
                //std::cout << "\n";
            }
            return maxi;
        }

        //sequences must be in to ncbi converted format
        int GaplessFilter_host_protein_converted_blosum62(
            const char* seq1,
            const char* seq2,
            const int length1,
            const int length2
        ) {

            //const int NEGINFINITY = -10000;
            std::vector<int> penalty_H(2*(length2+1));

            int maxi = 0, result;
            for (int index = 0; index <= length2; index++) {
                penalty_H[index] = 0;
            }

            auto BLOSUM = cudasw4::BLOSUM62_20::get2D();
            
            //std::cout << "CPU:\n";
            for (int row = 1; row <= length1; row++) {
                char seq1_char = seq1[row-1];
                char seq2_char;

                const int target_row = row & 1;
                const int source_row = !target_row;
                penalty_H[target_row*(length2+1)] = 0; //gap_open + (row-1)*gap_extend;
                for (int col = 1; col <= length2; col++) {
                    const int diag = penalty_H[source_row*(length2+1)+col-1];
                    seq2_char = seq2[col-1];

                    const int residue = BLOSUM[seq1_char][seq2_char];
                    result = std::max(0, diag + residue);
                    penalty_H[target_row*(length2+1)+col] = result;
                    if (result > maxi) maxi = result;
                }

                // for (int col = 1; col <= length2; col++) {
                //     printf("%2d ", penalty_H[target_row*(length2+1)+col]);
                // }
                // printf(", max %2d\n", maxi);
            }

            return maxi;
        }

        template<class OutputScores, class SelectedPositions>
        void runGaplessFilterKernels_PSSM(
            OutputScores& d_scores,
            GpuPermutedPSSMforGapless& permutedPSSM,
            const char* d_inputChars,
            const SequenceLengthT* d_inputLengths,
            const size_t* d_inputOffsets,
            SelectedPositions d_selectedPositions,
            size_t numSequences,
            char* d_tempStorage,
            size_t tempStorageBytes,
            cudaStream_t stream
        ){
            if(currentQueryLength <= getMaxSingleTileQueryLength_Gapless()){
                auto config = getSingleTileGroupRegConfigForPSSM_Gapless(currentQueryLength);

                if(!config.dpx){
                    if(config.approach == GaplessKernelConfig::Approach::hardcodedzero){
                        PSSM_2D_View<half2> strided_PSSM = permutedPSSM.makeHalf2View();
                        hardcodedzero::call_GaplessFilter_strided_PSSM_singletile_kernel<half2, 512>( 
                            config.groupsize, config.numRegs, d_inputChars, 
                            d_scores, d_inputOffsets, d_inputLengths, 
                            d_selectedPositions, numSequences, 
                            currentQueryLength, strided_PSSM, stream
                        );
                    }else{
                        PSSM_2D_View<half2> strided_PSSM = permutedPSSM.makeHalf2View();
                        kernelparamzero::call_GaplessFilter_strided_PSSM_singletile_kernel<half2, 512>(
                            config.groupsize, config.numRegs, d_inputChars, 
                            d_scores, d_inputOffsets, d_inputLengths, 
                            d_selectedPositions, numSequences, 
                            currentQueryLength, strided_PSSM, stream
                        );
                    }
                }else{
                    if(config.approach == GaplessKernelConfig::Approach::hardcodedzero){
                        PSSM_2D_View<short2> strided_PSSM = permutedPSSM.makeShort2View();
                        hardcodedzero::call_GaplessFilter_strided_PSSM_singletile_kernel<short2, 512>(
                            config.groupsize, config.numRegs, d_inputChars, 
                            d_scores, d_inputOffsets, d_inputLengths, 
                            d_selectedPositions, numSequences, 
                            currentQueryLength, strided_PSSM, stream
                        );
                    }else{
                        PSSM_2D_View<short2> strided_PSSM = permutedPSSM.makeShort2View();
                        kernelparamzero::call_GaplessFilter_strided_PSSM_singletile_kernel<short2, 512>(
                            config.groupsize, config.numRegs, d_inputChars, 
                            d_scores, d_inputOffsets, d_inputLengths, 
                            d_selectedPositions, numSequences, 
                            currentQueryLength, strided_PSSM, stream
                        );
                    }
                }
            }else{
                

                auto config = getMultiTileGroupRegConfigForPSSM_Gapless(currentQueryLength);

                int deviceId = 0;
                int numSMs = 0;
                cudaGetDevice(&deviceId);
                cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId);

                constexpr int threadBlockSize = 512;
                const int numGroupsPerBlock = threadBlockSize / config.groupsize;
                                            
                auto dbview = fullDB.getData();  //TODO only consider length of gpu partition, not full db
                const int maxSubjectLength = dbview.lengths()[dbview.numSequences()-1];

                //need to store 1 half value per subject position. kernel uses float2 for vectorized stores
                //4 halfs per float2
                const size_t tempStorageElementsPerGroup = SDIV(maxSubjectLength, 4);

                const size_t tempStorageElementsPerBlock = tempStorageElementsPerGroup * numGroupsPerBlock;
                const size_t tempStorageElementsAvailable = tempStorageBytes / sizeof(float2);

                const size_t maxNumBlocks = tempStorageElementsAvailable / tempStorageElementsPerBlock;
                if(maxNumBlocks == 0){
                    std::cout << "query with length " << currentQueryLength << " cannot be processed. ";
                    std::cout << "Not enough temp storage for a single threadblock. setting all scores to 0\n";
                    d_scores.setAllScoresToZero(stream);
                }else{
                    const int numThreadBlocks = numSMs; //std::min(maxNumBlocks, numSequences);

                    assert(sizeof(float2) * numThreadBlocks * tempStorageElementsPerBlock <= tempStorageBytes);

                    // std::cout << "maxSubjectLength: " << maxSubjectLength << ", numThreadBlocks: " << numThreadBlocks 
                    //     << ", tempStorageElementsPerBlock " << tempStorageElementsPerBlock
                    //     << ", tempStorage used: " << sizeof(float2) * size_t(numThreadBlocks) * tempStorageElementsPerBlock 
                    //     << " bytes" << "\n";

                    float2* const multiTileTempStorage = (float2*)d_tempStorage;

                    if(!config.dpx){
                        if(config.approach == GaplessKernelConfig::Approach::hardcodedzero){
                            PSSM_2D_View<half2> strided_PSSM = permutedPSSM.makeHalf2View();
                            hardcodedzero::call_GaplessFilter_strided_PSSM_multitile_kernel<half2, 512>(                              
                                numThreadBlocks, config.groupsize, config.numRegs, 
                                d_inputChars, 
                                d_scores, d_inputOffsets, d_inputLengths, 
                                d_selectedPositions, numSequences, 
                                currentQueryLength, strided_PSSM, 
                                multiTileTempStorage, tempStorageElementsPerGroup, stream
                            );
                        }else{
                            PSSM_2D_View<half2> strided_PSSM = permutedPSSM.makeHalf2View();
                            kernelparamzero::call_GaplessFilter_strided_PSSM_multitile_kernel<half2, 512>(
                                numThreadBlocks, config.groupsize, config.numRegs, 
                                d_inputChars, 
                                d_scores, d_inputOffsets, d_inputLengths, 
                                d_selectedPositions, numSequences, 
                                currentQueryLength, strided_PSSM, 
                                multiTileTempStorage, tempStorageElementsPerGroup, stream
                            );
                        }
                    }else{
                        if(config.approach == GaplessKernelConfig::Approach::hardcodedzero){
                            PSSM_2D_View<short2> strided_PSSM = permutedPSSM.makeShort2View();
                            hardcodedzero::call_GaplessFilter_strided_PSSM_multitile_kernel<short2, 512>(
                                numThreadBlocks, config.groupsize, config.numRegs, 
                                d_inputChars, 
                                d_scores, d_inputOffsets, d_inputLengths, 
                                d_selectedPositions, numSequences, 
                                currentQueryLength, strided_PSSM, 
                                multiTileTempStorage, tempStorageElementsPerGroup, stream
                            );
                        }else{
                            PSSM_2D_View<short2> strided_PSSM = permutedPSSM.makeShort2View();
                            kernelparamzero::call_GaplessFilter_strided_PSSM_multitile_kernel<short2, 512>(
                                numThreadBlocks, config.groupsize, config.numRegs, 
                                d_inputChars, 
                                d_scores, d_inputOffsets, d_inputLengths, 
                                d_selectedPositions, numSequences, 
                                currentQueryLength, strided_PSSM, 
                                multiTileTempStorage, tempStorageElementsPerGroup, stream
                            );
                        }
                    }
                }
            }
        }


        template<bool subjectIsCaseSensitive, bool withEndPosition, class OutputScores, class SelectedPositions>
        void run_SW_endposition_kernels_PSSM(
            OutputScores& d_scores,
            GpuPermutedPSSMforSW& permutedPSSM,
            const char* d_inputChars,
            const SequenceLengthT* d_inputLengths,
            const size_t* d_inputOffsets,
            SelectedPositions d_selectedPositions,
            size_t numSequences,
            char* d_tempStorage,
            size_t tempStorageBytes,
            cudaStream_t stream
        ){     

            if(currentQueryLength <= getMaxSingleTileQueryLength_SW()){
                auto config = getSingleTileGroupRegConfigForPSSM_SW(currentQueryLength);

                int deviceId = 0;
                int numSMs = 0;
                cudaGetDevice(&deviceId);
                cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId);

                const int numBlocks = std::min(size_t(numSMs), numSequences);

                if(!config.dpx){
                    PSSM_2D_View<float> strided_PSSM = permutedPSSM.makeView<float>();
                    call_amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_singletile<float, 512, withEndPosition, subjectIsCaseSensitive>(
                        numBlocks,
                        config.groupsize, config.numRegs, 
                        d_inputChars, 
                        d_scores, d_inputOffsets, d_inputLengths, 
                        d_selectedPositions, numSequences, 
                        currentQueryLength, strided_PSSM, 
                        gop,
                        gex,
                        stream
                    );
                }else{
                    PSSM_2D_View<int> strided_PSSM = permutedPSSM.makeView<int>();
                    call_amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_singletile<int, 512, withEndPosition, subjectIsCaseSensitive>(
                        numBlocks,
                        config.groupsize, config.numRegs, 
                        d_inputChars, 
                        d_scores, d_inputOffsets, d_inputLengths, 
                        d_selectedPositions, numSequences, 
                        currentQueryLength, strided_PSSM, 
                        gop,
                        gex,
                        stream
                    );
                }
            }else{
                auto config = getMultiTileGroupRegConfigForPSSM_SW(currentQueryLength);

                int deviceId = 0;
                int numSMs = 0;
                cudaGetDevice(&deviceId);
                cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId);

                constexpr int threadBlockSize = 512;
                const int numGroupsPerBlock = threadBlockSize / config.groupsize;
                                            
                auto dbview = fullDB.getData();  //TODO only consider length of gpu partition, not full db
                const int maxSubjectLength = dbview.lengths()[dbview.numSequences()-1];

                const int maxSubjectLengthPadded = maxSubjectLength + config.groupsize;
                const size_t tempBytesPerGroup = sizeof(float2) * maxSubjectLengthPadded;
                const size_t tempBytesPerBlock = tempBytesPerGroup * numGroupsPerBlock;

                //const size_t maxActiveGroups = std::min(numSequences, tempStorageBytes / tempBytesPerGroup);
                const size_t maxSimultaneousBlocks = tempStorageBytes / tempBytesPerBlock;
                if(maxSimultaneousBlocks == 0){
                    std::cout << "query with length " << currentQueryLength << " cannot be processed. ";
                    std::cout << "Not enough temp storage for a single threadblock. setting all scores to 0\n";
                    d_scores.setAllScoresToZero(stream);
                }else{
                    const int numBlocks = std::min(numSequences, std::min(size_t(numSMs), maxSimultaneousBlocks));
                    //const int numGroupsInGrid = numBlocks * numGroupsPerBlock;

                    // std::cout << "maxSubjectLengthPadded " << maxSubjectLengthPadded << "\n";
                    // std::cout << "tempBytesPerGroup " << tempBytesPerGroup << "\n";
                    // std::cout << "tempBytesPerBlock " << tempBytesPerBlock << "\n";
                    // std::cout << "maxSimultaneousBlocks " << maxSimultaneousBlocks << "\n";
                    // std::cout << "numBlocks " << numBlocks << "\n";
                    // std::cout << "numGroupsInGrid " << numGroupsInGrid << "\n";

                    if(!config.dpx){
                        PSSM_2D_View<float> strided_PSSM = permutedPSSM.makeView<float>();
                        call_amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_multitile<float, threadBlockSize, withEndPosition, subjectIsCaseSensitive>(
                            numBlocks,
                            config.groupsize, 
                            config.numRegs, 
                            d_inputChars, 
                            d_scores, 
                            d_inputOffsets, 
                            d_inputLengths, 
                            d_selectedPositions, 
                            numSequences, 
                            currentQueryLength, 
                            strided_PSSM, 
                            gop,
                            gex,
                            d_tempStorage, 
                            tempBytesPerGroup, 
                            stream
                        );
                    }else{
                        PSSM_2D_View<int> strided_PSSM = permutedPSSM.makeView<int>();
                        call_amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_multitile<int, threadBlockSize, withEndPosition, subjectIsCaseSensitive>(
                            numBlocks,
                            config.groupsize, 
                            config.numRegs, 
                            d_inputChars, 
                            d_scores, 
                            d_inputOffsets, 
                            d_inputLengths, 
                            d_selectedPositions, 
                            numSequences, 
                            currentQueryLength, 
                            strided_PSSM, 
                            gop,
                            gex,
                            d_tempStorage, 
                            tempBytesPerGroup, 
                            stream
                        );
                    }
                }
            }
        }

        void setNumTopNoCheck(int value){
            if(value >= 0){
                numTop = value;
                updateNumResultsPerQuery();

                cub::SwitchDevice sd(deviceIds[0]);
                const int numGpus = deviceIds.size();           

                h_finalAlignmentScores.resize(results_per_query);
                h_finalReferenceIds.resize(results_per_query);
                h_finalEndPositions.resize(results_per_query);
                d_finalAlignmentScores_allGpus.resize(results_per_query * numGpus);
                d_finalReferenceIds_allGpus.resize(results_per_query * numGpus);  
                d_finalEndPositions_allGpus.resize(results_per_query * numGpus);              
            }
        }


        std::vector<size_t> fullDB_numSequencesPerLengthPartition;
        std::vector<size_t> numSequencesPerGpu_total;
        std::vector<size_t> numSequencesPerGpuPrefixSum_total;

        //partition chars of whole DB amongst the gpus
        std::vector<size_t> numSequencesPerLengthPartitionPrefixSum;
        std::vector<DBdataView> dbPartitionsByLengthPartitioning;
        std::vector<std::vector<DBdataView>> subPartitionsForGpus;
        std::vector<std::vector<int>> lengthPartitionIdsForGpus;
        std::vector<size_t> numSequencesPerGpu;
        std::vector<size_t> numSequencesPerGpuPrefixSum;
        std::vector<CudaStream> gpuStreams;
        std::vector<CudaEvent> gpuEvents;
        std::vector<std::unique_ptr<GpuWorkingSet>> workingSets;  

        std::vector<std::vector<DeviceBatchCopyToPinnedPlan>> batchPlans;
        std::vector<std::vector<BatchDstInfo>> batchPlansDstInfoVec;

        std::vector<std::vector<DeviceBatchCopyToPinnedPlan>> batchPlans_cachedDB;
        std::vector<std::vector<BatchDstInfo>> batchPlansDstInfoVec_cachedDB;

        bool processingTheFirstQuery = true;
        int results_per_query;
        SequenceLengthT currentQueryLength;
        SequenceLengthT currentQueryLengthWithPadding;

        bool dbIsReady{};
        AnyDBWrapper fullDB;

        mutable std::unique_ptr<SequenceLengthStatistics> dbSequenceLengthStatistics;

        //final scan results. device data resides on gpu deviceIds[0]
        MyPinnedBuffer<float> h_finalAlignmentScores;
        MyPinnedBuffer<ReferenceIdT> h_finalReferenceIds;
        MyPinnedBuffer<AlignmentEndPosition> h_finalEndPositions;
        
        //MyPinnedBuffer<int> resultNumOverflows;
        MyDeviceBuffer<float> d_finalAlignmentScores_allGpus;
        MyDeviceBuffer<ReferenceIdT> d_finalReferenceIds_allGpus;
        MyDeviceBuffer<AlignmentEndPosition> d_finalEndPositions_allGpus;
        //MyDeviceBuffer<int> d_resultNumOverflows;
        std::unique_ptr<helpers::GpuTimer> scanTimer;

        size_t totalProcessedQueryLengths{};
        size_t totalNumOverflows{};
        std::unique_ptr<helpers::GpuTimer> totalTimer;

        HostGpuPartitionOffsets hostGpuPartitionOffsets;
        
        std::shared_ptr<TargetSubjectIds> targetSubjectIds;

        std::vector<GaplessKernelConfig> availableKernelConfigs_gapless_singletile;
        std::vector<GaplessKernelConfig> availableKernelConfigs_gapless_multitile;
        std::vector<SmithWatermanKernelConfig> availableKernelConfigs_sw_singletile;
        std::vector<SmithWatermanKernelConfig> availableKernelConfigs_sw_multitile;

        //--------------------------------------
        bool verbose = false;
        int gop = -11;
        int gex = -1;
        int numTop = 10;
        BlosumType blosumType = BlosumType::BLOSUM62_20;
        ScanType scanType = ScanType::Gapless;
        int maxReduceArraySize = MaxNumberOfResults::value();

        MemoryConfig memoryConfig;
        
        std::vector<int> deviceIds;

    };


} //namespace cudasw4

#endif


