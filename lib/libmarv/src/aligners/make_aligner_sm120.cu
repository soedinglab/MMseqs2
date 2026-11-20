#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"


#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

#if defined(__CUDACC__) && !(__CUDACC_VER_MAJOR__ < 13 || (__CUDACC_VER_MAJOR__ == 13 && __CUDACC_VER_MINOR__ < 2))
    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm120_8bit(){

        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 128, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 56>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 128, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 128, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 32, 36>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{22408,32484,37096,37791,41848,42785,42098,
                        42242,43157,43745,44487,44468,45521,46303,46238,46591,45366,45497,46209,46425,
                        47158,47578,47707,48251,46059,46395,46830,47420,50063,48106,50665,51225,48251,
                        48439,49266,48599,49996,50561,49619,51125,47545};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }

            
        };
        
    
        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 4>,
                // GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 8>,
                // GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 12>,
                // GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 16>,
                // GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 20>,
                // GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 24>,
                // GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 2, 56>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 256, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::UInt8x4, 512, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{17050,27868,32386,33633,36155,36783,38412,
            //             41279,41795,42799,43866,43564,45344,45213,44840,45088,44466,45559,45914,45779,
            //             47407,47802,45505,47120,45306,46466,46749,46343,47762,48035,46529,47040,46010,
            //             46954,46840,46372,47375,48178,46833,47808,45754};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }

            
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }
#endif    

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm120_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 8>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 256, 2, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 8, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 32, 36>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{8983,15629,18041,19502,20843,21705,21976,22128,21717,
                        22415,22812,22912,23202,23591,23839,23879,22081,22045,22680,22567,22832,24077,
                        24116,24245,23677,23865,23986,24233,24467,24650,24937,24860,23844,23824,24047,
                        24142,24414,24669,24737,24867,23787};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 4>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 8>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 12>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 16>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 20>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 24>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 28>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 32>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 44>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 48>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 52>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 56>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 32, 36>
            >;


            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{7354,13357,15633,17035,18792,19762,20383,20775,
            //             21903,22110,22513,22363,22957,23318,22915,22721,22609,23000,23349,23423,23901,24015,
            //             23791,23489,22915,23401,23492,23812,23951,24242,24034,23476,23153,23393,23436,23422,
            //             23597,23466,23547,23184,22657};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm120_32bit(){
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

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{723,1539,2360,3144,3927,4708,5503,6274,6976,7745,
                        8388,8690,8879,10021,10092,10393,10929,10932,11628,11323,11276,11626,11705,12036,
                        11994,11370,12412,11952,12087,12197,12362,12445,11516,11853,12183,11997,12291,12627,
                        12721,12639,11686};

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
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 12>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 16>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 20>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 24>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 28>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 32>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 40>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 48>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 52>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 56>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 64>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 48>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 52>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 56>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 60>,
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
            //     constexpr auto result = std::array{762,1567,2381,3181,3825,4562,5236,5900,6665,
            //             7219,7947,8397,9169,9581,9946,10181,10115,10603,11096,11133,11380,11549,11645,
            //             11677,11137,11402,11785,11739,12054,11774,11782,12271,11240,11556,11851,11729,
            //             11931,11946,11856,12170,11523};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm120_8bit_if_available_else_16bit(){


        #if defined(__CUDACC__) && !(__CUDACC_VER_MAJOR__ < 13 || (__CUDACC_VER_MAJOR__ == 13 && __CUDACC_VER_MINOR__ < 2))

        return make_multiconfig_gapless_aligner_sm120_8bit();

        #else

        return make_multiconfig_gapless_aligner_sm120_16bit();

        #endif
    }




#if 0
    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_sw_aligner_sm120_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                SWConfig<21, Approach::unused, Datatype::Float, 256, 8, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 256, 16, 16>,
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                SWConfig<21, Approach::unused, Datatype::Float, 256, 16, 16>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::SW);
    }
#endif


LIBMARV_NAMESPACE_END