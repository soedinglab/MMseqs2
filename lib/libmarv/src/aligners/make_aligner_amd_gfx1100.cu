#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"


#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_amd_gfx1100_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 4>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 8>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 12>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 16>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>
            >;


            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{1412,2078,2535,2960,3266,3564,3680,3870,4012,4163,
                    4384,4481,4646,4704,4699,4712,4830,4913,5070,5093,5140,5122,5155,5210,5181,5286,
                    5385,5382,5411,5390,5420,5454,5171,5291,5365,5407};



                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 4>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 8>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 12>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 16>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 20>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 24>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 28>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 32>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 36>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 44>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 48>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 52>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 56>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 60>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 64>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 36>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 44>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 48>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 52>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 56>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>
            >;


            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{1014,2116,2782,3114,3432,3639,3727,3809,3873,
            //         3934,3997,4055,4105,4175,4205,4282,3877,3959,4047,4117,4191,4234,4264,4295,
            //         4037,4118,4208,4276,4336,4391,4418,4435,4093,4167,4250,4308};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_amd_gfx1100_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 48>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{577,1308,1496,1597,1801,1994,2182,2279,
                    2374,2466,2538,2593,2639,2677,2743,2894,2989,3079,3116,3169,3192,3203,
                    3233,3312,3331,3376,3407,3431,3444,3433,3482,3531,3429,3474,3507,3537};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 4>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 8>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 12>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 16>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 20>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 24>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 28>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 32>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 40>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 48>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 28>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 32>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 24>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{453,1140,1855,2368,2478,2871,2872,
            //         2995,3116,3125,3222,3185,3123,3182,3259,3335,3457,3143,3306,3407,
            //         3551,3622,3663,3310,3470,3523,3684,3730,3776,3388};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };


        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_sw_aligner_amd_gfx1100_32bit(){
        //NOT TUNED
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