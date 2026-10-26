#include "pssmkernels_gapless_int8.cuh"

namespace cudasw4{

    #define ScoreOutputIterator TopNMaximaArray
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        template void uint8x4::call_GaplessFilter_strided_PSSM_singletile_uint8x4_kernel<ScoreType_u8x4, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<ScoreType_u8x4>&, \
            const int, \
            cudaStream_t \
        );

        PSSM_GAPLESS_INT8_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef ScoreOutputIterator


    #define ScoreOutputIterator TopNMaximaArray
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        template void uint8x4::call_GaplessFilter_strided_PSSM_multitile_uint8x4_kernel<ScoreType_u8x4, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<ScoreType_u8x4>&, \
            const int, \
            std::uint32_t*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_INT8_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef ScoreOutputIterator





} //namespace cudasw4