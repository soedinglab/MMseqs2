#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

    // AMD/HIP gapless-endpos aligner. Uses 32-bit float configs because RDNA has
    // no DPX packed-int16 instructions (Short2) and fp16 (Half2) loses precision
    // for large scores (e.g. foldseek 3di self-alignments).
    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_amd(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21,Approach::unused, Datatype::Float,128,4,16>,
                GaplessConfig<21,Approach::unused, Datatype::Float,128,4,32>,
                GaplessConfig<21,Approach::unused, Datatype::Float,128,8,32>,
                GaplessConfig<21,Approach::unused, Datatype::Float,128,16,32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        using MultiTileConfigs = SingleTileConfigs;


        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

LIBMARV_NAMESPACE_END
