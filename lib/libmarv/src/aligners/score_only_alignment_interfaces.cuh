#ifndef LIBMARV_SCOREONLY_ALIGNMENT_INTERFACES_CUH
#define LIBMARV_SCOREONLY_ALIGNMENT_INTERFACES_CUH

#include "../cuda_hip_rename.h"

#if defined(__HIPCC__)
    #include "hip/hip_runtime.h"
#endif


#include "../pssm.cuh"
#include "../util.cuh"
#include "../cuda_errorcheck.cuh"
#include "../alignment_algorithms/alignment_interface_data.cuh"
#include "../config.hpp"
#include "../offset_iterator.cuh"

#include <map>
#include <optional>
#include <cuda/std/mdspan>
#include <cuda/std/span>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

    struct ScoreOnlyAlignmentInterface{
        virtual ~ScoreOnlyAlignmentInterface() = default;

        virtual bool isSingleTile(const OneToAllInputDataPSSM& inputData) const = 0;

        virtual int getTileSize() const = 0;

        virtual AlignmentAlgorithmE getAlignmentAlgorithm() const = 0;

        virtual Datatype getScoreDatatype() const = 0;

        virtual void makeGpuPssm(
            GpuConvertedPSSM& result,
            const cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>> hostPssmView,
            cudaStream_t stream
        ) = 0;

        virtual size_t getMinimumSuggestedTempBytes_multiTile(const OneToAllInputDataPSSM& inputData) const = 0;

        virtual void scoreOnly_singleTile(
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs gapArgs,
            cudaStream_t stream
        ) = 0;

        virtual void scoreOnly_multiTile(
            char* d_temp,
            size_t tempBytes,
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs gapArgs,
            cudaStream_t stream
        ) = 0;
    };

    struct ScoreOnlyBatchQuery{
        int* scores;
        OneToAllInputDataPSSM input;
        const GpuConvertedPSSM* pssm;
    };

    struct MultiConfigScoreOnlyAlignerInterface{
        virtual ~MultiConfigScoreOnlyAlignerInterface() = default;

        // number of subjects scored per pass (size of the per-query score/id buffers)
        virtual int preferredKernelOutputArraySize() const { return 10'000'000; }

        // scores all queries of a pass together; returns false if the backend cannot, the caller then scores per query
        virtual bool scoreOnly_batch(cuda::std::span<const ScoreOnlyBatchQuery>, GapScoreArgs, cudaStream_t) { return false; }

        virtual bool isSingleTile(const OneToAllInputDataPSSM& inputData) const = 0;

        virtual int getTileSize(const OneToAllInputDataPSSM& inputData) const = 0;

        virtual AlignmentAlgorithmE getAlignmentAlgorithm() const = 0;

        virtual Datatype getScoreDatatype(int queryLength) const = 0;

        virtual void makeGpuPssm(
            GpuConvertedPSSM& result,
            const cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>> hostPssmView,
            cudaStream_t stream
        ) = 0;

        virtual size_t getMinimumSuggestedTempBytes_multiTile(const OneToAllInputDataPSSM& inputData) const = 0;

       
        virtual void scoreOnly_singleTile(
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs gapArgs,
            cudaStream_t stream
        ) = 0;

        virtual void scoreOnly_multiTile(
            char* d_temp,
            size_t tempBytes,
            int* const devAlignmentScores,
            const OneToAllInputDataPSSM& inputData,
            const GpuConvertedPSSM& strided_PSSM,
            GapScoreArgs gapArgs,
            cudaStream_t stream
        ) = 0;
    };
    


LIBMARV_NAMESPACE_END


#endif //LIBMARV_SCOREONLY_ALIGNMENT_INTERFACES_CUH