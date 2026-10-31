#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"


#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_amd_gfx1150_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 4>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 8>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 16>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{905,1694,2214,2433,2601,2706,2705,2729,2752,2793,
                        2825,2872,2832,2808,2749,2849,2864,2918,2874,2919,2876,2828,2819,2856,2864,2915,2873,
                        2889,2883,2800,2802,2819,2764,2804,2860,2846};

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
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 8>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 12>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 16>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 20>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 24>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 28>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 32>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 36>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 40>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 44>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 48>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 52>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 56>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 60>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 64>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 44>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 48>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 52>,
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
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{564,1022,1341,1541,1682,1758,1839,1905,1938,1967,
            //             1998,2020,2042,2066,2084,2091,1824,1867,1909,1945,1971,1988,2000,2008,1844,1886,1930,
            //             1960,1986,2001,2013,2018,1856,1899,1941,1971};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };


        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_amd_gfx1150_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 128, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 28>,
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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 48>,
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
                constexpr auto result = std::array{263,843,1067,1288,1443,1565,1621,1658,1666,
                        1748,1696,1702,1698,1691,1692,1718,1724,1769,1746,1765,1764,1734,
                        1740,1787,1747,1765,1759,1765,1761,1736,1739,1765,1726,1735,1753,1755};

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
            //     constexpr auto result = std::array{296,601,911,1166,1278,1381,1464,1437,1503,
            //             1475,1477,1568,1450,1496,1567,1590,1627,1462,1526,1559,1615,
            //             1631,1653,1497,1552,1582,1628,1645,1664,1506};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };


        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_sw_aligner_amd_gfx1150_32bit(){
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