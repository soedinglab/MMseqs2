#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN



    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_default_ptx_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                SWConfig<21,Approach::unused, Datatype::Float,128,4,16>,
                SWConfig<21,Approach::unused, Datatype::Float,128,4,32>,
                SWConfig<21,Approach::unused, Datatype::Float,128,8,32>,
                SWConfig<21,Approach::unused, Datatype::Float,128,16,32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        using MultiTileConfigs = SingleTileConfigs;
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::SW);
    }

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_default_ptx_32bit(){
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

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_default_ptx_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,4,16>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,4,32>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,8,32>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,16,32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        using MultiTileConfigs = SingleTileConfigs;
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }



LIBMARV_NAMESPACE_END