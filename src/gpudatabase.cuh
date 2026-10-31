#ifndef LIBMARV_GPU_DATABASE_CUH
#define LIBMARV_GPU_DATABASE_CUH

#include "cuda_hip_rename.h"

#if defined(__HIPCC__)
    #include "hip/hip_runtime.h"
#endif

#include <cuda/std/span>

#include "hpc_helpers/cuda_raiiwrappers.cuh"
#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/nvtx_markers.cuh"
#include "hpc_helpers/simple_allocation.cuh"

#include "config.hpp"
#include "backend_hooks.h"
#include "dbdata.hpp"
#include "length_partitions.hpp"
#include "util.cuh"
// #include "types.hpp"
#include "dbbatching.cuh"
#include "cuda_errorcheck.cuh"
#include "offset_iterator.cuh"
#include "gpudatabaseallocation.cuh"

#include <thrust/binary_search.h>
// #include <thrust/distance.h>
#include <thrust/scan.h>
#include <thrust/gather.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

#include <iostream>
// #include <string>
#include <vector>
#include <memory>
// #include <string_view>
// #include <optional>
// #include <sstream>

#include <cuda/cmath>

#include "namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN


    struct HostGpuPartitionOffsets; //forward declaration


    class GpuDatabase{
    public:
        template<class T>
        using MyPinnedBuffer = helpers::SimpleAllocationPinnedHost<T, 0>;
        template<class T>
        using MyDeviceBuffer = helpers::SimpleAllocationDevice<T, 0>;

        struct CallbackInputPerGpu{
            bool sequencesAreStoredConsecutively; //if true, all sequences reside in [d_inputChars, d_inputChars + numInputChars)
            int deviceId; //device to which pointers and stream belong
            int batchId; //it is the i-th batch of the current for-each call
            int sequenceLengthUpperBound; //elements of d_inputLengths are not greater than sequenceLengthUpperBound
            const char* d_inputChars; // size: numInputChars
            const SequenceLengthT* d_inputLengths; // size: numInputSequences. sequence i has length inputLengths[i]
            ReadOnlyCustomOffsetIterator<size_t*> d_inputOffsets; //size: (numInputSequences+1). sequence i begins at inputChars + (inputOffsets[i] - inputOffsets[0])
            size_t numInputSequences; // number of sequences in batch
            size_t numInputChars;
            cudaStream_t stream;
        };

        using ForEachCallback = std::function<void(cuda::std::span<CallbackInputPerGpu> callbackInput)>;

        struct GpuWorkingSet; //forward declaration

        struct SequenceLengthStatistics{
            SequenceLengthT max_length = 0;
            SequenceLengthT min_length = std::numeric_limits<SequenceLengthT>::max();
            size_t sumOfLengths = 0;
        };
    
        struct MemoryConfig{
            size_t maxBatchBytes = 128ull * 1024ull * 1024ull;
            size_t maxBatchSequences = 10'000'000;
            size_t maxGpuMem = std::numeric_limits<size_t>::max();
        };

    private:
        struct BatchDstInfo{
            bool isUploaded{};
            char* charsPtr{};
            SequenceLengthT* lengthsPtr{};
            size_t* offsetsPtr{};
        };

    public:

        GpuDatabase(
            std::vector<int> deviceIds_, 
            const MemoryConfig& memoryConfig_,
            bool verbose_
        );

        GpuDatabase() = delete;
        GpuDatabase(const GpuDatabase&) = delete;
        GpuDatabase(GpuDatabase&&) = default;
        GpuDatabase& operator=(const GpuDatabase&) = delete;
        GpuDatabase& operator=(GpuDatabase&&) = default;

        ~GpuDatabase();

        const std::vector<int>& getDeviceIds() const;


        void setDatabase(std::shared_ptr<DB> dbPtr);

        void setDatabase(std::shared_ptr<DBWithVectors> dbPtr);

        void setDatabase(std::shared_ptr<PseudoDB> dbPtr);

        void setDatabase(std::shared_ptr<MMseqsDB> dbPtr);

        void setDatabase(std::shared_ptr<ExternalDB> dbPtr);

        void setDatabase(std::shared_ptr<DB> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations);

        void setDatabase(std::shared_ptr<DBWithVectors> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations);

        void setDatabase(std::shared_ptr<PseudoDB> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations);

        void setDatabase(std::shared_ptr<MMseqsDB> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations);

        void setDatabase(std::shared_ptr<ExternalDB> dbPtr, const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations);

        std::vector<std::shared_ptr<GpuDatabaseAllocationBase>> getFullGpuDBAllocations();

        std::string_view getReferenceHeader(ReferenceIdT referenceId) const;

        int getReferenceLength(ReferenceIdT referenceId) const;

        std::string getReferenceSequence(ReferenceIdT referenceId) const;

        const char* getReferenceDatabaseSequenceHostPtr(ReferenceIdT referenceId);

        void prefetchDBToGpus();

        size_t getNumSequences() const;

        void printDBInfo() const;

        void printDBLengthPartitions() const;


    private:

        int fullscan_processCacheableBatches_notYetCached(
            ForEachCallback&& callback, 
            cuda::std::span<cudaStream_t> streams,
            int callbackBatchIdOffset
        );

        int fullscan_processCacheableBatches_alreadyCached(
            ForEachCallback&& callback, 
            cuda::std::span<cudaStream_t> streams,
            bool processCachedDBInSeparateBatches,
            int callbackBatchIdOffset
        );

        int fullscan_processCacheableBatches(
            ForEachCallback&& callback, 
            cuda::std::span<cudaStream_t> streams,
            bool processCachedDBInSeparateBatches,
            int callbackBatchIdOffset
        );

        //for each gpu, stream specified batches from host to device and process them
        int scan_processBatches_hoststreaming_pipelined(
            ForEachCallback&& callback, 
            cuda::std::span<cudaStream_t> streams,
            const std::vector<cuda::std::span<DeviceBatchCopyToPinnedPlan>>& batchPlansToProcessPerGpu,
            const std::vector<cuda::std::span<BatchDstInfo>>& destinationLocationsOfBatchesPerGpu,
            int callbackBatchIdOffset
        );

        int fullscan_processUncacheableBatches(
            ForEachCallback&& callback, 
            cuda::std::span<cudaStream_t> streams,
            int callbackBatchIdOffset
        );


    public:
        void forEachBatchOfSequencesInDatabase(
            ForEachCallback&& callback, 
            cuda::std::span<cudaStream_t> streams,
            bool processCachedDBInSeparateBatches
        );

        //streams[i] belongs to gpu deviceIds[i]
        void forEachBatchOfSequencesInDatabaseSubset(
            ForEachCallback&& callback, 
            cuda::std::span<ReferenceIdT> sortedSubsetRefIds,
            cuda::std::span<cudaStream_t> streams
        );

        SequenceLengthStatistics getSequenceLengthStatistics() const;

        //optional 256-entry table mapping stored bytes to alphabet codes (Metal and SIMD backends)
        void setLetterTranslation(const unsigned char* table){ letterTranslation.assign(table, table + 256); }
        const unsigned char* getLetterTranslation() const{ return letterTranslation.empty() ? nullptr : letterTranslation.data(); }

        //transform per gpu local sequence indices into global sequence indices
        void convertLocalIndicesToGlobalIndices(
            int gpuPartion, 
            ReferenceIdT* d_indices,
            int size,
            cudaStream_t stream
        );
        
    private:

        void initializeGpus();

        // check if we can use GPU unified memory
        bool canUseHostDbInPlace() const;

        // check if DB fits in host memory
        bool hostDbFitsInMemory() const;

        void resetDatabaseState();
        void makeReady();

        void makeReadyWithExistingFullGpuDB(const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations);

        void computeTotalNumSequencePerLengthPartition();

        void partitionDBAmongstGpus();

        void allocateGpuWorkingSets();

        void allocateGpuWorkingSetsWithExistingFullGpuDB(const std::vector<std::shared_ptr<GpuDatabaseAllocationBase>>& existingFullGpuDBAllocations);

        void createDBBatchesForGpus();

        void assignBatchesToGpuMem();

        void markCachedDBBatchesAsUploaded(int gpu);

        void printDBDataView(const DBdataView& view) const;

        void printDBDataViews(const std::vector<DBdataView>& views) const;

        std::vector<DeviceBatchCopyToPinnedPlan> computeDbCopyPlan(
            const std::vector<DBdataView>& dbPartitions,
            const std::vector<int>& lengthPartitionIds,
            size_t MAX_CHARDATA_BYTES,
            size_t MAX_SEQ
        ) const;

        bool shouldUseHostFuncForTransfer() const;



        int numCopyThreadsPerGpu = 1;


        std::vector<bool> allCacheableBatchesAreUploaded;

        std::vector<size_t> fullHostDB_numSequencesPerLengthPartition;
        std::vector<size_t> numSequencesPerGpu_total;
        std::vector<size_t> numSequencesPerGpuPrefixSum_total;

        //partition chars of whole DB amongst the gpus
        std::vector<size_t> numSequencesPerLengthPartitionPrefixSum;
        std::vector<DBdataView> dbPartitionsByLengthPartitioning;
        std::vector<std::vector<DBdataView>> subPartitionsForGpus;
        std::vector<std::vector<int>> lengthPartitionIdsForGpus;
        std::vector<size_t> numSequencesPerGpu;
        std::vector<size_t> numSequencesPerGpuPrefixSum;
        std::vector<std::unique_ptr<GpuWorkingSet>> workingSets;  

        std::vector<std::vector<DeviceBatchCopyToPinnedPlan>> batchPlans;
        std::vector<std::vector<BatchDstInfo>> batchPlansDstInfoVec;

        std::vector<std::vector<DeviceBatchCopyToPinnedPlan>> batchPlans_cachedDB;
        std::vector<std::vector<BatchDstInfo>> batchPlansDstInfoVec_cachedDB;

        bool dbIsReady{};
        AnyDBWrapper fullHostDB;
        // no-copy registrations of the host arrays while they are scanned in place
        std::shared_ptr<void> hostRegistration;

        mutable std::unique_ptr<SequenceLengthStatistics> dbSequenceLengthStatistics;

        std::unique_ptr<HostGpuPartitionOffsets> hostGpuPartitionOffsets;

        //--------------------------------------
        bool verbose = false;

        MemoryConfig memoryConfig;
        
        std::vector<int> deviceIds;

        std::vector<unsigned char> letterTranslation;

    };


LIBMARV_NAMESPACE_WITH_NESTING_END

#endif


