#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm110_16bit(){
        //NOTE: copied from a different arch, not tuned.
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,4>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,8>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,12>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,16>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,20>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,24>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,28>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,32>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,36>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,48>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,56>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,64>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,36>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,40>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,44>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,48>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,52>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 8,60>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,64>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,36>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 16,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 16,48>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,52>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,56>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,60>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,64>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        //NOTE: copied from a different arch, not tuned.
        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,4>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,8>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,12>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,16>,
                // GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,20>,
                // GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,24>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,28>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,32>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,36>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,48>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,56>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 4,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 4,64>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,36>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,40>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,44>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,48>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,52>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 8,60>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 8,64>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,36>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 16,44>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Short2, 512, 16,48>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,52>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,56>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,60>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Short2, 512, 16,64>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm110_32bit(){
        //NOTE: copied from a different arch, not tuned.
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{723,1539,2360,3144,3927,4708,5503,6274,6976,7745,
            //             8388,8690,8879,10021,10092,10393,10929,10932,11628,11323,11276,11626,11705,12036,
            //             11994,11370,12412,11952,12087,12197,12362,12445,11516,11853,12183,11997,12291,12627,
            //             12721,12639,11686};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 36>
            >;


            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


LIBMARV_NAMESPACE_END