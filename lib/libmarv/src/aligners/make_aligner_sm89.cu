#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm89_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 4>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 8>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 12>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 36>
            >;


            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{4332,8190,11363,11740,12737,13318,13598,13839,14413,
                        14673,15013,15164,15269,15590,15829,15683,15031,15155,15631,15770,15899,15992,16196,
                        16250,15253,15430,15749,16048,16115,16211,16354,16439,15580,16140,16701,17169,16662,
                        16280,17163,17454,15745};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 4>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 4>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 12>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 16>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 20>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 24>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 28>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 32>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 44>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 48>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 52>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 56>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 36>
            >;


            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{2725,5860,7322,9343,12106,12956,13835,13807,14751,
            //             14753,14712,14814,15449,14019,14313,14087,14363,14336,14799,14861,15149,15291,15085,
            //             15947,14830,15360,15615,15565,15903,15137,16069,16356,15618,15845,16041,15940,16359,
            //             16202,16395,16617,15919};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }



    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm89_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 20>,
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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>
            >;


            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{612,1266,1744,2325,2901,3459,4009,4547,5071,
                        5472,5884,6393,6770,7181,7408,7671,7480,7429,7618,7732,7884,8098,8116,8177,7840,
                        7780,7928,7958,8150,8297,8355,8403,8255,8604,8777,8353,8556,8333,8680,8785,8692};

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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>
            >;


            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{595,1229,1689,2313,3084,3675,4284,4864,5400,
            //             5916,6424,6965,7244,7089,5938,7509,6974,7216,7163,7506,7461,7483,7618,7904,7929,
            //             7783,7909,8094,8071,8175,8196,8098,7876,8103,8229,8400,8396,8391,8253,8359,8395};

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