#include "pssmkernels_gapless.cuh"

namespace cudasw4{

namespace hardcodedzero{

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        template void call_GaplessFilter_strided_PSSM_singletile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        template void call_GaplessFilter_strided_PSSM_singletile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator


    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        template void call_GaplessFilter_strided_PSSM_multitile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator


    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        template void call_GaplessFilter_strided_PSSM_multitile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator


} //namespace hardcodedzero



} //namespace cudasw4