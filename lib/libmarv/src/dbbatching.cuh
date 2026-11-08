#ifndef DBBATCHING_CUH
#define DBBATCHING_CUH

#include "cuda_hip_rename.h"

#include "config.hpp"
#include "dbdata.hpp"

#include <vector>
#include <iostream>
#include <algorithm>

#include <thrust/for_each.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/scan.h>

#include <omp.h>

#include "namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

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


    void executeCopyPlanH2DDirect(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* d_chardata,
        SequenceLengthT* d_lengthdata,
        size_t* d_offsetdata,
        const std::vector<DBdataView>& dbPartitions,
        cudaStream_t stream
    );

    void executePinnedCopyPlanParallel(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* h_chardata,
        SequenceLengthT* h_lengthdata,
        size_t* h_offsetdata,
        const std::vector<DBdataView>& dbPartitions,
        int numThreads
    );
    
    void executePinnedCopyPlanSerial(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* h_chardata,
        SequenceLengthT* h_lengthdata,
        size_t* h_offsetdata,
        const std::vector<DBdataView>& dbPartitions
    );
    
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
    );
            
    void executePinnedCopyPlanWithHostCallback(
        const DeviceBatchCopyToPinnedPlan& plan, 
        char* h_chardata,
        SequenceLengthT* h_lengthdata,
        size_t* h_offsetdata,
        const std::vector<DBdataView>& dbPartitions, 
        cudaStream_t stream,
        int numThreads = 1
    );

LIBMARV_NAMESPACE_END

#endif