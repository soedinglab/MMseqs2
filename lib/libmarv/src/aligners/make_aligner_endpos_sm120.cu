#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_sm120_32bit(){
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
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 20>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 4, 44>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 24>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 28>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 32>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 8, 36>,
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 20>,
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
                SWConfig<21, Approach::unused, Datatype::Int, 512, 16, 20>,
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


#if 0    
    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm120_8bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 4>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 8>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 12>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 16>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 512, 4, 20>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 24>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 28>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 128, 4, 32>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 36>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 40>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 44>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 48>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 52>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 56>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 60>,
                GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 4, 64>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 8, 36>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 8, 40>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 8, 44>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 8, 48>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 8, 52>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 8, 56>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 8, 60>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 8, 64>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 256, 16, 36>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 512, 16, 40>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 512, 16, 44>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 512, 16, 48>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 512, 16, 52>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 512, 16, 56>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 512, 16, 60>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 512, 16, 64>,
                // GaplessConfig<21,Approach::unused, Datatype::UInt8x4, 512, 32, 36>,
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
               GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 4, 16>,
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }
#endif

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm120_16bit(){

        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 8>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 12>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 20>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 36>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{6452,9090,9634,10429,10242,10704,10207,10761,
                        10633,10660,10660,10952,11486,10886,10569,10785,11064,10876,10842,11211,11642,
                        10897,10619,10867,11102,10898,10934,11053,11701,10906,10650,10868,10796,10759,
                        10749,10794,10913,10815,10873,10871,10966};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 4>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 8>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 12>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 16>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 20>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 24>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 28>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 32>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 36>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 24>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 28>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 20>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{4714,6262,6820,7004,7244,7317,7351,7550,7444,
            //             7215,7408,7363,6910,7371,7272,7281,6923,7268,6933,7412,7398,
            //             7277,7296,7003,7278,7131,7414,7385,7186,7239,6911,7197,
            //             6855,7322,6976,7287,7210};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm120_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{753,1568,2376,3165,3926,4692,5397,6106,5985,7008,
                        6688,7551,7079,7460,7305,7769,7285,7546,7568,7767,7754,7861,7851,7939,7690,7814,
                        7860,7987,8018,8061,8010,8050,7803,7852,7925,8007,8032,8095,8110,8087,8115};

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
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 8>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 12>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 16>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 20>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 24>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 28>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 32>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 24>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 52>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 28>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 60>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 32>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 24>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 52>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 28>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{754,1556,2321,2982,3538,4087,4296,4510,4720,4580,
            //             4565,4830,4494,4988,4161,5016,5168,4934,4961,5022,4812,5103,4307,5163,5283,5153,
            //             5128,5124,4979,5221,4363,5354,5322,5297,5226,5402,4809,5400,4062,5379,5329};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };   
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }




LIBMARV_NAMESPACE_END