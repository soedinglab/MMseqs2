#include "dbbatching.cuh"

#include "cuda_hip_rename.h"

#include "config.hpp"
#include "cuda_errorcheck.cuh"
#include "util.cuh"
#include "dbdata.hpp"

#include <vector>
#include <iostream>
#include <algorithm>

#include <thrust/for_each.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>

#include <omp.h>

#include "namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    struct ExecutePinnedCopyCallbackData{
        int numThreads;
        const DeviceBatchCopyToPinnedPlan* planPtr; 
        char* h_chardata;
        SequenceLengthT* h_lengthdata;
        size_t* h_offsetdata;
        const std::vector<DBdataView>* dbPartitionsPtr;
    };

    void executeCopyPlanH2DDirect(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* d_chardata,
        SequenceLengthT* d_lengthdata,
        size_t* d_offsetdata,
        const std::vector<DBdataView>& dbPartitions,
        cudaStream_t stream
    ){
        size_t usedBytes = 0;
        size_t usedSeq = 0;
        for(const auto& copyRange : plan.copyRanges){
            const auto& dbPartition = dbPartitions[copyRange.currentCopyPartition];
            const auto& firstSeq = copyRange.currentCopySeqInPartition;
            const auto& numToCopy = copyRange.numToCopy;
            size_t numBytesToCopy = dbPartition.offsets()[firstSeq + numToCopy] - dbPartition.offsets()[firstSeq];
    
            CUDACHECK(cudaMemcpyAsync(
                d_chardata + usedBytes,
                dbPartition.chars() + dbPartition.offsets()[firstSeq],
                numBytesToCopy,
                H2D,
                stream
            ));
            CUDACHECK(cudaMemcpyAsync(
                d_lengthdata + usedSeq,
                dbPartition.lengths() + firstSeq,
                sizeof(SequenceLengthT) * numToCopy,
                H2D,
                stream
            ));
            CUDACHECK(cudaMemcpyAsync(
                d_offsetdata + usedSeq,
                dbPartition.offsets() + firstSeq,
                sizeof(size_t) * (numToCopy+1),
                H2D,
                stream
            ));
            thrust::for_each(
                thrust::cuda::par_nosync.on(stream),
                d_offsetdata + usedSeq,
                d_offsetdata + usedSeq + (numToCopy+1),
                [
                    usedBytes,
                    firstOffset = dbPartition.offsets()[firstSeq]
                ] __device__ (size_t& off){
                    off = off - firstOffset + usedBytes;
                }
            );
    
            usedBytes += numBytesToCopy;
            usedSeq += numToCopy;
        }
    };

    void executePinnedCopyPlanParallel(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* h_chardata,
        SequenceLengthT* h_lengthdata,
        size_t* h_offsetdata,
        const std::vector<DBdataView>& dbPartitions,
        int numThreads
    ){
        size_t usedBytes = 0;
        size_t usedSeq = 0;
        for(const auto& copyRange : plan.copyRanges){
            const auto& dbPartition = dbPartitions[copyRange.currentCopyPartition];
            const auto& firstSeq = copyRange.currentCopySeqInPartition;
            const auto& numToCopy = copyRange.numToCopy;
            size_t numBytesToCopy = dbPartition.offsets()[firstSeq + numToCopy] - dbPartition.offsets()[firstSeq];

            #pragma omp parallel num_threads(numThreads)
            {
                const size_t workers = omp_get_num_threads();
                const size_t tid = omp_get_thread_num();
                const size_t chunksize = numBytesToCopy / workers;
                const size_t remainder = numBytesToCopy % workers;
                const size_t begin = tid * chunksize + std::min(tid, remainder);
                const size_t end = begin + chunksize + (tid < remainder);
                std::copy(
                    dbPartition.chars() + dbPartition.offsets()[firstSeq] + begin,
                    dbPartition.chars() + dbPartition.offsets()[firstSeq] + end,
                    h_chardata + usedBytes + begin
                );
            }

            std::copy(
                dbPartition.lengths() + firstSeq,
                dbPartition.lengths() + firstSeq+numToCopy,
                h_lengthdata + usedSeq
            );
            std::transform(
                dbPartition.offsets() + firstSeq,
                dbPartition.offsets() + firstSeq + (numToCopy+1),
                h_offsetdata + usedSeq,
                [&](size_t off){
                    return off - dbPartition.offsets()[firstSeq] + usedBytes;
                }
            );

            usedBytes += numBytesToCopy;
            usedSeq += numToCopy;
        }
    };
    
    void executePinnedCopyPlanSerial(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* h_chardata,
        SequenceLengthT* h_lengthdata,
        size_t* h_offsetdata,
        const std::vector<DBdataView>& dbPartitions
    ){
        size_t usedBytes = 0;
        size_t usedSeq = 0;
        for(const auto& copyRange : plan.copyRanges){
            const auto& dbPartition = dbPartitions[copyRange.currentCopyPartition];
            const auto& firstSeq = copyRange.currentCopySeqInPartition;
            const auto& numToCopy = copyRange.numToCopy;
            size_t numBytesToCopy = dbPartition.offsets()[firstSeq + numToCopy] - dbPartition.offsets()[firstSeq];
    
            auto end = std::copy(
                dbPartition.chars() + dbPartition.offsets()[firstSeq],
                dbPartition.chars() + dbPartition.offsets()[firstSeq + numToCopy],
                h_chardata + usedBytes
            );
            std::copy(
                dbPartition.lengths() + firstSeq,
                dbPartition.lengths() + firstSeq+numToCopy,
                h_lengthdata + usedSeq
            );
            std::transform(
                dbPartition.offsets() + firstSeq,
                dbPartition.offsets() + firstSeq + (numToCopy+1),
                h_offsetdata + usedSeq,
                [&](size_t off){
                    return off - dbPartition.offsets()[firstSeq] + usedBytes;
                }
            );
            usedBytes += std::distance(h_chardata + usedBytes, end);
            usedSeq += numToCopy;
        }
    };
    
    void executePinnedCopyPlanSerialAndTransferToGpu(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* h_chardata,
        SequenceLengthT* h_lengthdata,
        size_t* /*h_offsetdata*/,
        char* d_chardata,
        SequenceLengthT* d_lengthdata,
        size_t* d_offsetdata,
        const std::vector<DBdataView>& dbPartitions,
        cudaStream_t H2DcopyStream
    ){
    
        size_t usedBytes = 0;
        for(const auto& copyRange : plan.copyRanges){
            const auto& dbPartition = dbPartitions[copyRange.currentCopyPartition];
            const auto& firstSeq = copyRange.currentCopySeqInPartition;
            const auto& numToCopy = copyRange.numToCopy;
            size_t numBytesToCopy = dbPartition.offsets()[firstSeq + numToCopy] - dbPartition.offsets()[firstSeq];
            constexpr size_t maxTransferBatchSize = 8 * 1024 * 1024;
            for(size_t i = 0; i < numBytesToCopy; i += maxTransferBatchSize){
                const size_t x = std::min(numBytesToCopy - i, maxTransferBatchSize);
    
                std::copy_n(
                    dbPartition.chars() + dbPartition.offsets()[firstSeq] + i,
                    x,
                    h_chardata + usedBytes + i
                );
                CUDACHECK(cudaMemcpyAsync(
                    d_chardata + usedBytes + i,
                    h_chardata + usedBytes + i,
                    x,
                    H2D,
                    H2DcopyStream
                ));
            }
    
            // auto end = std::copy(
            //     dbPartition.chars() + dbPartition.offsets()[firstSeq],
            //     dbPartition.chars() + dbPartition.offsets()[firstSeq + numToCopy],
            //     h_chardata + usedBytes
            // );
            // CUDACHECK(cudaMemcpyAsync(
            //     d_chardata + usedBytes,
            //     h_chardata + usedBytes,
            //     numBytesToCopy,
            //     H2D,
            //     H2DcopyStream
            // ));
    
            usedBytes += numBytesToCopy;
        }
    
        size_t usedSeq = 0;
        for(const auto& copyRange : plan.copyRanges){
            const auto& dbPartition = dbPartitions[copyRange.currentCopyPartition];
            const auto& firstSeq = copyRange.currentCopySeqInPartition;
            const auto& numToCopy = copyRange.numToCopy;
    
            std::copy(
                dbPartition.lengths() + firstSeq,
                dbPartition.lengths() + firstSeq+numToCopy,
                h_lengthdata + usedSeq
            );
            // CUDACHECK(cudaMemcpyAsync(
            //     d_lengthdata + usedSeq,
            //     h_lengthdata + usedSeq,
            //     sizeof(size_t) * numToCopy,
            //     H2D,
            //     H2DcopyStream
            // ));
    
            usedSeq += numToCopy;
        }
        cudaMemcpyAsync(
            d_lengthdata,
            h_lengthdata,
            sizeof(SequenceLengthT) * plan.usedSeq,
            H2D,
            H2DcopyStream
        ); CUERR;
    
        CUDACHECK(cudaMemsetAsync(d_offsetdata, 0, sizeof(size_t), H2DcopyStream));
    
        auto d_paddedLengths = thrust::make_transform_iterator(
            d_lengthdata,
            [] __host__ __device__ (const SequenceLengthT& length){
                return size_t(SDIV(length, 4) * 4);
            }
        );
    
        thrust::inclusive_scan(
            thrust::cuda::par_nosync(thrust_async_allocator<char>(H2DcopyStream)).on(H2DcopyStream),
            d_paddedLengths,
            d_paddedLengths + plan.usedSeq,
            d_offsetdata + 1
        );
    
    };
        
    void executePinnedCopyPlanCallback(void* args){
        ExecutePinnedCopyCallbackData* callbackData = (ExecutePinnedCopyCallbackData*)args;
        const auto& plan = *callbackData->planPtr;
        auto& dbPartitions = *callbackData->dbPartitionsPtr;
        
        if(callbackData->numThreads > 1){
            executePinnedCopyPlanParallel(
                plan, 
                callbackData->h_chardata,
                callbackData->h_lengthdata,
                callbackData->h_offsetdata,
                dbPartitions,
                callbackData->numThreads
            );
        }else{
            executePinnedCopyPlanSerial(
                plan, 
                callbackData->h_chardata,
                callbackData->h_lengthdata,
                callbackData->h_offsetdata,
                dbPartitions
            );
        }
    
        delete callbackData;
    }
    
    void executePinnedCopyPlanWithHostCallback(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* h_chardata,
        SequenceLengthT* h_lengthdata,
        size_t* h_offsetdata,
        const std::vector<DBdataView>& dbPartitions, 
        cudaStream_t stream,
        int numThreads
    ){
        ExecutePinnedCopyCallbackData* data = new ExecutePinnedCopyCallbackData;
    
        data->planPtr = &plan;
        data->h_chardata = h_chardata,
        data->h_lengthdata = h_lengthdata,
        data->h_offsetdata = h_offsetdata,
        data->dbPartitionsPtr = &dbPartitions;
        data->numThreads = numThreads;
    
        CUDACHECK(cudaLaunchHostFunc(
            stream,
            executePinnedCopyPlanCallback,
            (void*)data
        ));
    }

LIBMARV_NAMESPACE_END
