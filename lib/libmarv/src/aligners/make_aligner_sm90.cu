#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm90_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 2, 8>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 2, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 2, 32>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 256, 2, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 256, 2, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 2, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 2, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 4, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 4, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 256, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 4, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 4, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 8, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 128, 8, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 16, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 16, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 256, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 16, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 128, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 32, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 32, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 32, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 32, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 32, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 32, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 32, 60>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 256, 32, 64>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{4202,9227,11513,12064,13195,13644,13797,13122,14340,14495,
                        14642,14630,14949,15047,15065,14510,14822,14856,15027,14903,15258,15292,15410,15150,14971,14990,
                        15168,15259,15375,15392,15506,15424,14955,15105,15208,15315,15399,15482,15466,15094,15023,15212,
                        15009,15133,15146,15471,15248,14991};

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
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 16>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 20>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 24>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 28>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 32>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 36>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 40>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 44>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 48>,
                // GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 2, 52>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 56>,
                // GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 2, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 2, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 8, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 8, 64>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 16, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 16, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 16, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 32, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 32, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 32, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 512, 32, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 32, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Short2, 512, 32, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 32, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Short2, 256, 32, 64>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{4141,7830,10194,11021,12217,12739,13247,12587,13725,13827,
            //             14011,14033,14259,14241,14328,13893,14249,14202,14417,14373,14644,14628,14576,14447,14336,
            //             14382,14550,14628,14773,14864,14765,14609,14436,14565,14657,14805,14900,15063,14866,14737,
            //             14621,14731,14931,15055,15050,15122,14880,14927};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_sm90_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 2, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 128, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 256, 32, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 32, 64>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{409,884,1356,1841,2319,2784,3254,3710,4168,4624,
                        5032,5509,5883,6273,6532,6836,7001,7308,7211,7446,7512,7626,7560,7402,7391,7560,7532,
                        7645,7714,7769,7729,7607,7311,7524,7598,7727,7779,7809,7754,7700,7627,7746,7778,7806,
                        7847,7883,7892,7874};

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
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 8>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 12>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 16>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 20>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 24>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 28>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 32>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 36>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 40>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 44>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 48>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 52>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 2, 56>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 32>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 36>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 40>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 44>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 48>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 52>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 56>,
                // GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Int, 512, 16, 60>,
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

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{435,918,1392,1860,2305,2778,3237,3706,4139,
            //             4586,5019,5451,5864,6264,5959,6357,6724,6990,7120,7145,7290,7223,7126,6767,6993,
            //             7147,7280,7315,7331,7460,7326,7250,7141,7293,7433,7496,7461,7606,7456,7390,7456,
            //             7558,7614,7662,7696,7692,7581,7601};

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