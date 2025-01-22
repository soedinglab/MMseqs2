#ifndef DBBATCHING_CUH
#define DBBATCHING_CUH

#include "config.hpp"

#include <vector>
#include <iostream>
#include <algorithm>

#include <thrust/for_each.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/scan.h>

namespace cudasw4{

    struct DeviceBatchCopyToPinnedPlan{
        struct CopyRange{
            int lengthPartitionId;
            int currentCopyPartition;
            int currentCopySeqInPartition;
            int numToCopy;
        };
        size_t usedBytes = 0;
        size_t usedSeq = 0;
        std::vector<int> h_partitionIds;
        std::vector<size_t> h_numPerPartition;
        std::vector<CopyRange> copyRanges;

        friend std::ostream& operator<<(std::ostream& os, const DeviceBatchCopyToPinnedPlan& plan){
            os << "usedBytes " << plan.usedBytes << ", usedSeq " << plan.usedSeq << " ";
            for(int i = 0; i < int(plan.h_partitionIds.size()); i++){
                os << "(" << plan.h_partitionIds[i] << "," << plan.h_numPerPartition[i] << ") ";
            }
            
            return os;
        }
    };

    struct ExecutePinnedCopyCallbackData{
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
    
            cudaMemcpyAsync(
                d_chardata + usedBytes,
                dbPartition.chars() + dbPartition.offsets()[firstSeq],
                numBytesToCopy,
                H2D,
                stream
            ); CUERR;
            cudaMemcpyAsync(
                d_lengthdata + usedSeq,
                dbPartition.lengths() + firstSeq,
                sizeof(SequenceLengthT) * numToCopy,
                H2D,
                stream
            ); CUERR;
            cudaMemcpyAsync(
                d_offsetdata + usedSeq,
                dbPartition.offsets() + firstSeq,
                sizeof(size_t) * (numToCopy+1),
                H2D,
                stream
            ); CUERR;
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
                cudaMemcpyAsync(
                    d_chardata + usedBytes + i,
                    h_chardata + usedBytes + i,
                    x,
                    H2D,
                    H2DcopyStream
                ); CUERR;
            }
    
            // auto end = std::copy(
            //     dbPartition.chars() + dbPartition.offsets()[firstSeq],
            //     dbPartition.chars() + dbPartition.offsets()[firstSeq + numToCopy],
            //     h_chardata + usedBytes
            // );
            // cudaMemcpyAsync(
            //     d_chardata + usedBytes,
            //     h_chardata + usedBytes,
            //     numBytesToCopy,
            //     H2D,
            //     H2DcopyStream
            // ); CUERR;
    
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
            // cudaMemcpyAsync(
            //     d_lengthdata + usedSeq,
            //     h_lengthdata + usedSeq,
            //     sizeof(size_t) * numToCopy,
            //     H2D,
            //     H2DcopyStream
            // ); CUERR;
    
            usedSeq += numToCopy;
        }
        cudaMemcpyAsync(
            d_lengthdata,
            h_lengthdata,
            sizeof(SequenceLengthT) * plan.usedSeq,
            H2D,
            H2DcopyStream
        ); CUERR;
    
        cudaMemsetAsync(d_offsetdata, 0, sizeof(size_t), H2DcopyStream); CUERR;
    
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
        
    
        executePinnedCopyPlanSerial(
            plan, 
            callbackData->h_chardata,
            callbackData->h_lengthdata,
            callbackData->h_offsetdata,
            dbPartitions
        );
    
        delete callbackData;
    }
    
    void executePinnedCopyPlanWithHostCallback(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* h_chardata,
        SequenceLengthT* h_lengthdata,
        size_t* h_offsetdata,
        const std::vector<DBdataView>& dbPartitions, 
        cudaStream_t stream
    ){
        ExecutePinnedCopyCallbackData* data = new ExecutePinnedCopyCallbackData;
    
        data->planPtr = &plan;
        data->h_chardata = h_chardata,
        data->h_lengthdata = h_lengthdata,
        data->h_offsetdata = h_offsetdata,
        data->dbPartitionsPtr = &dbPartitions;
    
        cudaLaunchHostFunc(
            stream,
            executePinnedCopyPlanCallback,
            (void*)data
        ); CUERR;
    }

}

#endif