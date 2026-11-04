#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm80_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 8>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 12>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 16>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 20>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 4, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 4, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 8, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 60>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{1860,3029,4363,4940,5278,5440,5541,5537,
                        5614,5616,5636,5659,5718,5728,5705,5754,5692,5711,5734,5748,5830,6120,6566,
                        6952,6966,6993,7026,7065,7151,7111,7073,7112,6991,7004,7029,7010,7159,7120,
                        7075,7113,6984,7049,7058,7018,7102,7171,7130};

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
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 8>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 12>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 16>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 20>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 24>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 28>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 32>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 44>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 48>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 52>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 56>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 60>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{1685,2789,4732,5420,5888,6142,6314,6474,6682,
            //             6721,6743,6730,6735,6864,6907,6925,6797,6838,6823,6743,6849,6866,6969,7093,6837,
            //             6876,6852,6812,6865,6863,6913,6976,6862,6882,6856,6822,6873,6876,6924,6984,6883,
            //             6901,6950,6967,6983,7004,6977};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm80_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 8, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 60>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{262,538,662,875,1072,1269,1441,1630,1819,
                        1988,2147,2317,2433,2563,2684,2886,2971,2985,3019,3092,3104,3165,3319,3535,
                        3748,3771,3828,3872,3892,3941,3970,3963,3778,3801,3830,3898,3916,3959,3993,
                        3956,3703,3802,3908,3833,3936,3964,3999};

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
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 52>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 56>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 32>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 48>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 52>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 56>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 60>,
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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 60>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{256,531,798,1057,1302,1535,1754,1988,2175,
            //             2382,2573,2770,2964,3133,3251,3277,3485,3551,3615,3692,3630,3647,3621,
            //             3675,3593,3646,3658,3747,3788,3717,3676,3726,3643,3691,3718,3810,
            //             3846,3769,3709,3760,3749,3824,3791,3869,3811,3835,3865};

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