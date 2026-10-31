#ifndef LIBMARV_SEARCH_COMMON_CUH
#define LIBMARV_SEARCH_COMMON_CUH

#include "../cuda_hip_rename.h"

#if defined(__HIPCC__)
    #include "hip/hip_runtime.h"
#endif

#include <memory>
#include <vector>

#include "../aligners/score_only_alignment_interfaces.cuh"
#include "../aligners/score_with_endpos_alignment_interfaces.cuh"

#include "../config.hpp"

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

struct BenchmarkStats{
    double seconds{};
    double gcups{};
};

struct SearchResult{
    std::vector<int> scores{};
    std::vector<ReferenceIdT> referenceIds{};
    std::vector<AlignmentEndPosition> endPositions{};
    std::unique_ptr<BenchmarkStats> stats{};
};

struct IsOverflowScore{
    int threshold;

    __host__ __device__
    IsOverflowScore(int t) : threshold(t){}

    __host__ __device__
    bool operator()(int score){
        return score >= threshold;
    }
};

struct IsIndexToOverflowScore{
    int threshold;
    const int* scores;

    __host__ __device__
    IsIndexToOverflowScore(int t, const int* s) : threshold(t), scores(s){}

    __host__ __device__
    bool operator()(int index){
        return scores[index] >= threshold;
    }
};


struct OverflowSettings{
    enum class AlignerType{
        ScoreOnly,
        ScoreEndpos
    };

    bool checkForOverflow{};
    AlignerType alignerType{};
    cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*> overflowAlignersScoreOnly{};
    cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> overflowAlignersScoreWithEndpos{};


    static OverflowSettings withoutOverflowCorrection(){
        OverflowSettings result;
        result.checkForOverflow = false;
        return result;
    }

    //re-compute all alignments with a score >= threshold
    static OverflowSettings withOverflowCorrection(cuda::std::span<libmarv::MultiConfigScoreOnlyAlignerInterface*> aligners){
        OverflowSettings result;
        result.checkForOverflow = true;
        result.overflowAlignersScoreOnly = aligners;
        result.alignerType = AlignerType::ScoreOnly;
        return result;
    }

    //re-compute all alignments with a score >= threshold
    static OverflowSettings withOverflowCorrection(cuda::std::span<libmarv::MultiConfigScoreWithEndPosAlignerInterface*> aligners){
        OverflowSettings result;
        result.checkForOverflow = true;
        result.overflowAlignersScoreWithEndpos = aligners;
        result.alignerType = AlignerType::ScoreEndpos;
        return result;
    }
};

struct MaskingOptions{
    MaskingOptions() = default;
    MaskingOptions(int threshold) : MaskingOptions(threshold, 20){}
    MaskingOptions(int threshold, char letter) : maskingLetter(letter), maskingThreshold(threshold){}

    char maskingLetter = 20;
    int maskingThreshold = 6; //mask runs of length >= maskingThreshold
};



LIBMARV_NAMESPACE_END

#endif