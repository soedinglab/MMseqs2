#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_sm80_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 4>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 8>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 12>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 32>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 32>,
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
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 32>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 20>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 24>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 28>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 16, 32>,
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


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm80_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 20>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 60>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{1640,2526,2584,2717,2627,2655,2655,2804,2682,
                        2712,2693,2679,2804,3079,3270,3267,3294,3329,3306,3276,3295,3290,3290,3279,3306,
                        3311,3313,3299,3299,3299,3300,3322,3308,3315,3315,3338,3301,3351,3303,3329,3319,
                        3301,3308,3294,3279,3329,3340};

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
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 24>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 28>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 16>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 20>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 44>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 12>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 52>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 56>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 16>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 20>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 32, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 32, 60>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{1413,2174,2341,2483,2555,2545,2410,2531,2300,
            //             2615,2347,2490,2395,2413,2378,2550,2232,2598,2474,2570,2431,2450,2320,2561,2357,
            //             2612,2490,2580,2453,2458,2359,2580,2368,2618,2518,2601,2466,2466,2367,2528,2372,
            //             2422,2525,2504,2536,2450,2420};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };   
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm80_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 52>,
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
                constexpr auto result = std::array{246,521,644,854,1064,1256,1370,1439,1424,
                        1445,1458,1475,1482,1513,1600,1792,1779,1797,1807,1818,1829,1836,1849,1828,
                        1803,1812,1824,1829,1842,1848,1860,1842,1808,1817,1829,1830,1844,1849,1863,
                        1839,1846,1859,1865,1866,1872,1877,1886};

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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 60>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{251,526,783,1043,1270,1488,1657,1720,1753,1799,
            //             1836,1842,1863,1874,1875,1875,1838,1856,1865,1889,1907,1908,1915,1897,1893,1911,
            //             1914,1925,1939,1945,1956,1936,1902,1916,1922,1933,1943,1951,1961,1938,1963,1981,
            //             1982,1985,1996,1994,1993};

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