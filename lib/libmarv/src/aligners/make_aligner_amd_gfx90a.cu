#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"


#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_amd_gfx90a_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 8>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 20>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{2850,4418,5125,5443,5663,5810,5881,5989,5625,5803,
                    5751,5932,5998,6088,5618,5858,5815,5975,6041,6129,
                    5662,5865,5849,5921,4225,4278,4301,4406,3963,4240,
                    4132,4296};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        struct MultiTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 48>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{3711,3770,3816,3881,3891,3899,3952,3722,3791,3850,
                    3882,3956,3954,3974,4028,3740,3822,3889,3922};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };


        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_amd_gfx90a_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 48>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{1258,2313,2890,3211,3134,3156,3077,3330,3374,3249,
                    3373,3214,3426,3369,3460,3434,3443,3395,3133,3416,
                    3498,3003,3124,2700,2810,2793,2879};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        struct MultiTileConfigs{
            using type = std::tuple<
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

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{2537,2567,2554,2494,2439,2555,2591,2646,2618,2548,
                    2464};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };


        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_sw_aligner_amd_gfx90a_32bit(){
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
