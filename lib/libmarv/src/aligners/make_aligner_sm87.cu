#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm87_16bit(){
        //NOTE: copied from a different arch, not tuned.
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,4>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,8>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,12>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,16>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,20>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,24>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,28>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,32>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,36>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,40>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,44>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,48>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,64>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 8,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,44>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 8,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,52>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 8,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,64>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,40>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,44>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,56>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,60>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,64>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        //NOTE: copied from a different arch, not tuned.
        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,4>,
                // GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,8>,
                // GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,12>,
                // GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,16>,
                // GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,20>,
                // GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,24>,
                // GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,28>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,32>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,36>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,40>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,44>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,48>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 4,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 4,64>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 8,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,40>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,44>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 8,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,52>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 8,56>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,60>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 8,64>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,36>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,40>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,44>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,48>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,52>,
                GaplessConfig<21,Approach::hardcodedzero, Datatype::Half2, 512, 16,56>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,60>,
                GaplessConfig<21,Approach::kernelparamzero, Datatype::Half2, 512, 16,64>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm87_32bit(){
        //NOTE: copied from a different arch, not tuned.
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{252,520,788,1051,1307,1561,1810,2052,
            //             2291,2472,2667,2873,3057,3272,3329,3427,3342,3380,3430,3511,3553,3630,3644,
            //             3668,3483,3513,3535,3552,3556,3580,3585,3565,3269,3349,3305,
            //             3396,3363,3401,3408,3429,3261};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>,
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