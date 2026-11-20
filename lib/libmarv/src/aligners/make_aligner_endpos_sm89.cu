#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_sm89_32bit(){

        struct SingleTileConfigs{
            using type = std::tuple<
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 4>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 8>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 12>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 32>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
        struct MultiTileConfigs{
            using type = std::tuple<
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 32>
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::SW);
    }


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm89_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 4>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 8>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 12>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 8, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 36>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{3714,4597,4984,5202,5336,5393,5432,5500,5509,5555,
                        5518,5525,5605,5607,5563,5529,5561,5590,5576,5588,5635,6049,6157,6165,6136,6122,6137,
                        6143,6189,6159,6155,6155,6132,6085,6106,6108,6178,6148,6142,6149,6173};

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
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 12>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 16>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 20>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 24>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 28>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 32>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 44>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 48>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 52>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 56>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 32, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{2762,3693,4516,4727,4815,4658,4718,4824,4692,4724,
            //             4913,4955,4999,4845,4675,4866,4717,4871,4976,5043,4944,4722,4629,4833,4788,4889,5003,
            //             5043,5027,4853,4705,4851,4821,4913,5012,5055,5059,4907,4702,4994,4842};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };   
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm89_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 56>,
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

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{577,1227,1700,2277,2852,3322,3570,3599,3634,
                        3712,3740,3806,3825,3853,3830,3832,3748,3800,3823,3854,3876,3898,3868,3874,3875,
                        4196,4236,4267,4289,4306,4267,4275,4182,4199,4229,4260,4285,4289,4262,4266,4266};

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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{581,1215,1668,2462,2984,3499,3865,4087,4142,4277,
            //             4349,4388,4443,4497,4536,4483,4401,4535,4505,4557,4664,4636,4679,4618,4561,4625,4645,
            //             4691,4748,4758,4748,4749,4585,4648,4669,4714,4767,4781,4760,4766,4766};

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