#ifndef LIBMARV_SCORE_ONLY_FOLDSEEK_SEARCH_CUH
#define LIBMARV_SCORE_ONLY_FOLDSEEK_SEARCH_CUH

#include "../cuda_hip_rename.h"

#if defined(__HIPCC__)
    #include "hip/hip_runtime.h"
#endif

#include <cuda/std/span>
#include <cuda/std/mdspan>
#if defined(__CUDACC__)
#include <cuda/std/tuple>
#else
#include <thrust/tuple.h>
#endif

#include <vector>
#include <memory>

#include "../aligners/score_only_alignment_interfaces.cuh"
#include "../aligners/score_with_endpos_alignment_interfaces.cuh"
#include "../config.hpp"
#include "../gpudatabase.cuh"
#include "search_common.cuh"

#include "../hpc_helpers/simple_allocation.cuh"
#include "../hpc_helpers/cuda_raiiwrappers.cuh"


#include "../namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN


struct ScoreOnlyFoldseekSearch{

    template<class T>
    using MyPinnedBuffer = helpers::SimpleAllocationPinnedHost<T, 0>;
    template<class T>
    using MyDeviceBuffer = helpers::SimpleAllocationDevice<T, 0>;

    //forward
    struct GpuWorkingSet;

    struct CompareScoresDescendingRefIdsAscending{
        template<class Tuple1, class Tuple2>
        __host__ __device__
        bool operator()(const Tuple1& lhs, const Tuple2& rhs) const{
            #if defined(__CUDACC__)
            const auto scoreL = cuda::std::get<0>(lhs);
            const auto refIdL = cuda::std::get<1>(lhs);
            const auto scoreR = cuda::std::get<0>(rhs);
            const auto refIdR = cuda::std::get<1>(rhs);
            #else
            const auto scoreL = thrust::get<0>(lhs);
            const auto refIdL = thrust::get<1>(lhs);
            const auto scoreR = thrust::get<0>(rhs);
            const auto refIdR = thrust::get<1>(rhs);
            #endif
            if(scoreL < scoreR) return false;
            if(scoreL > scoreR) return true;
            //scores are equal
            return refIdL < refIdR;
        }
    };

    struct ConvertCombinedTo3di{
        __host__ __device__
        char operator()(char combined) const{
            return combined / 12;
        }
    };

    struct ConvertCombinedTo3diDummy{
        __host__ __device__
        char operator()(char combined) const{
            return combined % 21;
        }
    };

    struct ConvertCombinedToAA12{
        __host__ __device__
        char operator()(char combined) const{
            return combined % 12;
        }
    };

    int gop = 0;
    int gex = 0;
    int kernelOutputArraySize = 1000000;
    int topNSize = 10000;
    int maxNumQueries = 1;
    std::vector<int> deviceIds;
    cuda::std::span<MultiConfigScoreOnlyAlignerInterface*> standardAlignerPerGpu;

    cuda::std::span<MultiConfigScoreWithEndPosAlignerInterface*> standardEndposAlignerPerGpu;

    OverflowSettings scoreOnlyOverflowSettings;
    OverflowSettings endposOverflowSettings;
    MaskingOptions maskingOptions;

    std::vector<std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>>> queryPssmPerGpuPerQuery_forStandardScoreOnly;
    std::vector<std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>>> queryPssmPerGpuPerQuery_forOverflowScoreOnly;
    std::vector<std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>>> queryPssmPerGpuPerQuery_forStandardEndPos;
    std::vector<std::vector<std::unique_ptr<libmarv::GpuConvertedPSSM>>> queryPssmPerGpuPerQuery_forOverflowEndPos;

    std::vector<std::vector<std::unique_ptr<libmarv::GpuPlainPSSM>>> queryPssmPerGpuPerQuery_aa12;

    cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> hostQueryPssmViews_3di;
    cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> hostQueryPssmViews_aa12;
    

    std::vector<std::unique_ptr<GpuWorkingSet>> workingSets;
    std::vector<MyPinnedBuffer<int>> h_finalAlignmentScoresPerQuery;
    std::vector<MyPinnedBuffer<ReferenceIdT>> h_finalReferenceIdsPerQuery;
    std::unique_ptr<helpers::GpuTimer> scanTimer;

    MyPinnedBuffer<int> h_numOverflowsPerGpuPerQuery;
    std::vector<CudaStream> streamsPerGpuPerQuery;

    std::vector<std::vector<int>> resultSizePerQueryPerGpu;
    std::vector<int> processedSequencesPerGpu;

    std::vector<int> currentQueryLengths;
    std::vector<std::pair<int,int>> currentMinMaxValuesInPssms;
    bool int8IsAllowed;

    GpuDatabase* gpuDatabase;

    bool isOrdinaryAA20db = true; //for testing with a standard protein db.

    ScoreOnlyFoldseekSearch(bool isAA20db, //for testing with a standard protein db
        GpuDatabase& gpuDatabase_, int kernelOutputArraySize_, int topNSize_, int maxNumQueries_ = 1);

    ~ScoreOnlyFoldseekSearch();


    void reset();

    void destroyWorkingSets();

    void destroyPssms_3di();

    void destroyPssms_aa12();

    std::vector<SearchResult> execute(
        cuda::std::span<MultiConfigScoreOnlyAlignerInterface*> alignersScoreOnly,
        const OverflowSettings& scoreOnlyOverflowSettings_,
        cuda::std::span<MultiConfigScoreWithEndPosAlignerInterface*> alignersEndPos,
        const OverflowSettings& endposOverflowSettings_,
        cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews_3di,
        cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> queryPssmViews_aa12,
        MaskingOptions maskingOptions_,
        int gapopenscore = 0,
        int gapextendscore = 0
    );

public: //public to enable extended device lambda , should be private

    void forEachCallback_fullscan(cuda::std::span<GpuDatabase::CallbackInputPerGpu> callbackInput);

private:

    void consumeKernelOutputArraysFullScan(int queryIndex, int gpu, cudaStream_t stream);

    void finalizeFullScanResults(cuda::std::span<cudaStream_t> streams);

    void consumeKernelOutputArraysSubScan(int queryIndex, int gpu, cudaStream_t stream);

    void finalizeSubScanResults(int queryIndex, cuda::std::span<cudaStream_t> streams);

public: //public to enable extended device lambda, should be private

    void forEachCallback_subset(cuda::std::span<GpuDatabase::CallbackInputPerGpu> callbackInput, int queryIndex);

private:

    void initWorkingSets();

    void setAligner(
        cuda::std::span<MultiConfigScoreOnlyAlignerInterface*> alignersScoreOnly,
        const OverflowSettings& scoreOnlyOverflowSettings_,
        cuda::std::span<MultiConfigScoreWithEndPosAlignerInterface*> alignersEndPos,
        const OverflowSettings& endposOverflowSettings_
    );


    int getOverflowThreshold(int queryIndex, Datatype datatype);

    void setPssms_3di(cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> hostQueryPssmViews_3di_);

    
    void setPssms_aa12(cuda::std::span<cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>>> hostQueryPssmViews_aa12_);

    void clearProcessedSequencesAndCachedResults();

    void clearProcessedSequencesAndCachedResults(int queryIndex);

};


LIBMARV_NAMESPACE_WITH_NESTING_END

#endif