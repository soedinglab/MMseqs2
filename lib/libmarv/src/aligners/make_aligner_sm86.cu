#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm86_16bit(){

        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 8>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 20>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 36>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{1790,3714,5069,5359,5681,5981,6116,6229,6469,6646,
                        6826,6800,6879,6959,7032,7077,6773,6792,6981,6983,7089,7125,7192,7252,6809,6819,6919,
                        6885,6890,6856,6887,6855,6398,6406,6508,6512,6527,6543,6537,6602,6286};

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
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 12>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 16>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 20>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 24>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 28>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 32>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 44>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 48>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 52>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 56>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{1594,3238,4164,4697,5062,5134,5318,5485,5557,5727,
            //             5802,5715,5807,5868,5912,5946,5880,5998,6094,6047,6119,6189,6176,6243,5937,6064,6130,
            //             6082,6173,6203,6233,6290,5979,6096,6154,6112,6184,6221,6239,6300,6040};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm86_32bit(){
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

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{252,520,788,1051,1307,1561,1810,2052,
                        2291,2472,2667,2873,3057,3272,3329,3427,3342,3380,3430,3511,3553,3630,3644,
                        3668,3483,3513,3535,3552,3556,3580,3585,3565,3269,3349,3305,
                        3396,3363,3401,3408,3429,3261};

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
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 28>,
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
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{245,508,768,1021,1266,1507,1712,1939,2089,
            //             2300,2463,2543,2601,2734,2595,2811,2895,2937,3003,3046,3057,3106,3069,3085,
            //             3017,3050,3098,3144,3185,3182,3108,3118,3062,3109,3158,3196,3231,3214,
            //             3146,3151,3211};

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