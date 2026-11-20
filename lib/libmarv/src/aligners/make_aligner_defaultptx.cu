#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN



    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_default_ptx_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,4,16>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,4,32>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,4,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,8,32>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,8,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2,128,16,32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        using MultiTileConfigs = SingleTileConfigs;
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_default_ptx_32bit(){
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
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }





    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_sw_aligner_default_ptx_32bit(){
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
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::SW);
    }


LIBMARV_NAMESPACE_END