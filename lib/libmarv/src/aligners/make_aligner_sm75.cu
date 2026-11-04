#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm75_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,4>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,8>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,12>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,16>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,20>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,24>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,28>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,32>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,64>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,64>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,48>
                //larger tiles are not supported because shared memory size is too small
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,4>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,8>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,12>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,16>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,20>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,24>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,28>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,32>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,64>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,64>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,48>
                //larger tiles are not supported because shared memory size is too small
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }



    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm75_32bit(){
        //not tuned
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,4>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,8>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,12>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,16>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,20>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,24>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,28>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,32>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,36>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,40>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,44>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,48>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,52>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,56>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,60>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,64>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,36>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,40>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,44>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,48>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,52>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,56>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,60>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,64>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 16,36>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 16,40>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 16,44>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 16,48>
                //larger tiles are not supported because shared memory size is too small
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        //not tuned
        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,4>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,8>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,12>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,16>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,20>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,24>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,28>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,32>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,36>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,40>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,44>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,48>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,52>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,56>,
                // GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,60>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 4,64>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,36>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,40>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,44>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,48>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,52>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,56>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,60>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 8,64>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 16,36>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 16,40>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 16,44>,
                GaplessConfig<21,Approach::unused, Datatype::Float, 512, 16,48>
                //larger tiles are not supported because shared memory size is too small
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


LIBMARV_NAMESPACE_END