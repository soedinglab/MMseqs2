#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm86_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 8>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 128, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 24>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 128, 2, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 32, 36>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{1554,2103,2261,2323,2391,2424,2434,2461,2475,2487,
                        2486,2493,2495,2527,2522,2487,2490,2502,2498,2501,2503,2532,2526,2496,2494,2505,2504,
                        2508,2509,2538,2522,2501,2496,2497,2505,2498,2502,2516,2507,2527,2516};

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
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 12>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 16>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 20>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 24>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 28>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 32>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 36>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 40>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 44>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 48>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 52>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 56>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 16>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 20>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 32, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 16, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 32, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{1164,1655,1829,1918,1948,1879,1845,1963,1883,1949,
            //             1925,1994,1989,1986,1908,1952,1900,1978,2007,2022,2022,1949,1867,1971,1914,1981,2008,
            //             2032,2056,2000,1909,1982,1923,1997,2022,2041,2060,2011,1916,2019,1959};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm86_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 60>,
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
                constexpr auto result = std::array{236,503,768,1028,1287,1518,1626,1668,1658,1693,
                        1709,1724,1734,1749,1756,1742,1711,1730,1743,1749,1760,1767,1777,1759,1735,1744,
                        1759,1761,1772,1779,1787,1768,1737,1746,1761,1763,1771,1781,1785,1770,1780};

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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 36>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{239,502,754,1008,1230,1444,1576,1652,1676,
            //             1734,1773,1768,1777,1811,1794,1800,1841,1868,1875,1891,
            //             1914,1933,1935,1917,1912,1930,1942,1952,1968,1986,1991,
            //             1978,1924,1936,1954,1969,1982,2000,2001,1980,1984};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_sm86_32bit(){
        //NOTE: copied from a different arch,
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

LIBMARV_NAMESPACE_END