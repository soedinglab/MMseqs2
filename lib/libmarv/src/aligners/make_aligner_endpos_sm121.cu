#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_sm121_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 4>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 8>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 12>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 32>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 32>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 40>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 32>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 32>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 40>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 32>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::SW);
    }

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm121_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 58>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 46>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 54>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 62>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 42>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 50>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 34>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 38>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 58>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 36>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 4>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 12>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 42>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 54>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 34>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 38>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 46>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 50>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 58>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 62>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 42>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 36>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm121_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 4, 18>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 58>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 62>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 42>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 46>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 50>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 34>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 38>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 54>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 58>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>
                >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 4>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 20>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 28>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 46>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 48>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 50>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 52>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 54>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 56>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 58>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 32>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 38>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 42>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 34>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 58>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 62>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 64>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

LIBMARV_NAMESPACE_END