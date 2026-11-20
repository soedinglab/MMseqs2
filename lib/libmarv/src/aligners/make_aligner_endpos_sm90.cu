#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm90_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 8>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 12>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 8, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 32, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 32, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 64>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{2830,3709,3975,4165,4269,4310,4335,4394,4400,4449,4443,
                        4450,4476,4433,4432,4428,4425,4470,4464,4461,4493,4467,4446,4443,4453,4487,4478,4474,4505,
                        4475,4464,4461,4459,4490,4482,4475,4507,4478,4430,4419,4464,4496,4477,
                        4506,4487,4505,4429,4423};

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
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 16>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 20>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 24>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 28>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 32>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 20>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 44>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 24>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 20>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 32, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 32, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 32, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 64>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{2146,2944,3265,3186,3342,3374,3516,3510,3523,3420,3090,3515,
            //             3520,3586,3402,3473,3021,3515,3564,3619,3546,3486,3111,3517,3571,
            //             3617,3549,3476,3551,3590,3585,3600,3316,3143,3023,2892,2807,2545,2431};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm90_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 64>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{427,904,1381,1857,2326,2740,2922,2995,2986,3042,
                        3040,3094,3098,3113,3108,3121,3059,3102,3120,3148,3155,3160,3151,3152,3107,3134,
                        3152,3168,3178,3178,3173,3170,3111,3141,3156,3171,3185,3182,3176,3174,3189,3195,
                        3211,3212,3230,3242,3222,3207};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            
            using type = std::tuple<
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 4>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 8>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 12>,
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
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 64>,
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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 64>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{431,908,1366,1821,2222,2611,2888,3054,3168,3218,
            //             3253,3299,3364,3356,3370,3357,3319,3347,3387,3422,3458,3459,3452,3400,3393,3446,
            //             3467,3483,3511,3518,3510,3482,3413,3456,3479,3482,3520,3529,3517,3490,3522,
            //             3566,3526,3573,3584,3595,3603,3550};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_sm90_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 4>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 8>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 12>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 16>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 20>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 24>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 28>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 32>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 36>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 40>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 44>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 24>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 28>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 32>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 36>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 40>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 44>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 24>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 28>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 32>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 36>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 32, 20>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 44>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 32, 24>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 32, 28>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 32, 32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 32>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 36>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 40>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 44>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 24>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 28>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 32>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 36>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 32, 20>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 44>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 32, 24>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 32, 28>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 32, 32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::SW);
    }

LIBMARV_NAMESPACE_END