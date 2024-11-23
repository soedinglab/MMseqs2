//#include "validtileconfigs.hpp"
#include "util.cuh"
#include "config.hpp"
#include "pssmkernels_smithwaterman.cuh"

#include <thrust/iterator/counting_iterator.h>

namespace cudasw4{


#define ScoreOutputIterator TopNMaximaArrayWithExtra<AlignmentEndPosition>
#define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
#define withEndPosition true
#define subjectIsCaseSensitive true
#define X(g,r) \
    template void call_amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_singletile<float, 512, g, r, withEndPosition, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
        int, \
        const char * const, \
        ScoreOutputIterator const, \
        const size_t* const, \
        const SequenceLengthT* const, \
        PositionsIterator const, \
        const int, \
        const SequenceLengthT, \
        const PSSM_2D_View<float>&, \
        const float,  \
        const float, \
        cudaStream_t \
    );

PSSM_SW_ENDPOS_SINGLETILE_FLOAT_OR_INT_FOR_EACH_VALID_CONFIG_DO_X

#undef X
#undef subjectIsCaseSensitive
#undef withEndPosition
#undef SequenceLengthT
#undef PositionsIterator
#undef ScoreOutputIterator


#define ScoreOutputIterator TopNMaximaArrayWithExtra<AlignmentEndPosition>
#define PositionsIterator ReferenceIdT*
#define withEndPosition true
#define subjectIsCaseSensitive true
#define X(g,r) \
    template void call_amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_singletile<float, 512, g, r, withEndPosition, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
        int, \
        const char * const, \
        ScoreOutputIterator const, \
        const size_t* const, \
        const SequenceLengthT* const, \
        PositionsIterator const, \
        const int, \
        const SequenceLengthT, \
        const PSSM_2D_View<float>&, \
        const float,  \
        const float, \
        cudaStream_t \
    );

PSSM_SW_ENDPOS_SINGLETILE_FLOAT_OR_INT_FOR_EACH_VALID_CONFIG_DO_X

#undef X
#undef subjectIsCaseSensitive
#undef withEndPosition
#undef SequenceLengthT
#undef PositionsIterator
#undef ScoreOutputIterator



#define ScoreOutputIterator TopNMaximaArrayWithExtra<AlignmentEndPosition>
#define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
#define withEndPosition true
#define subjectIsCaseSensitive true
#define X(g,r) \
    template void call_amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_multitile<float, 512, g, r, withEndPosition, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
        int, \
        const char * const, \
        ScoreOutputIterator const, \
        const size_t* const, \
        const SequenceLengthT* const, \
        PositionsIterator const, \
        const int, \
        const SequenceLengthT, \
        const PSSM_2D_View<float>&, \
        const float,  \
        const float, \
        char* const, \
        const size_t, \
        cudaStream_t \
    );

    PSSM_SW_ENDPOS_MULTITILE_FLOAT_OR_INT_FOR_EACH_VALID_CONFIG_DO_X

#undef X
#undef subjectIsCaseSensitive
#undef withEndPosition
#undef PositionsIterator
#undef ScoreOutputIterator



#define ScoreOutputIterator TopNMaximaArrayWithExtra<AlignmentEndPosition>
#define PositionsIterator ReferenceIdT*
#define withEndPosition true
#define subjectIsCaseSensitive true
#define X(g,r) \
    template void call_amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_multitile<float, 512, g, r, withEndPosition, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
        int, \
        const char * const, \
        ScoreOutputIterator const, \
        const size_t* const, \
        const SequenceLengthT* const, \
        PositionsIterator const, \
        const int, \
        const SequenceLengthT, \
        const PSSM_2D_View<float>&, \
        const float,  \
        const float, \
        char* const, \
        const size_t, \
        cudaStream_t \
    );

    PSSM_SW_ENDPOS_MULTITILE_FLOAT_OR_INT_FOR_EACH_VALID_CONFIG_DO_X

#undef X
#undef subjectIsCaseSensitive
#undef withEndPosition
#undef PositionsIterator
#undef ScoreOutputIterator



} //namespace cudasw4