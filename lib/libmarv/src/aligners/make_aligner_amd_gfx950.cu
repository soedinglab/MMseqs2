#include "score_only_alignment.cuh"
#include "score_only_alignment_interfaces.cuh"


#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_amd_gfx950_16bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 4>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 8>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 12>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 16>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 20>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 24>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 28>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 2, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 40>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 2, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 2, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 4, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 256, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 256, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 28>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 32>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{10888,16639,18550,19710,20081,20422,20558,20920,20848,20879,
                    20615,20854,20982,21034,21485,21249,21680,21064,21408,21452,
                    21618,21966,21822,22215,21680,22029,22108,22265,22861,22340,
                    22595,21916,22216};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        struct MultiTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 2, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 48>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 52>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 56>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 4, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 4, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 36>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 48>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 52>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 8, 56>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 60>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 8, 64>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 36>,
                GaplessConfig<21, Approach::kernelparamzero, Datatype::Half2, 512, 16, 40>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 44>,
                GaplessConfig<21, Approach::hardcodedzero, Datatype::Half2, 512, 16, 48>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{14382,13301,13303,13638,13811,13916,14113,14090,14502,13714,
                    13685,13945,14101,14209,14397,14359,14761,13926,13936,14171,
                    14334};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };


        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }


    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_gapless_aligner_amd_gfx950_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 64>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 48>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{5718,9609,11025,11426,12238,12184,13140,12575,12161,12952,
                    12777,13569,14187,14721,12970,13482,13420,14178,13800,14900,
                    15144,13332,14230,13836,14741,14377,15494,15581,14321,15110,
                    15660};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        struct MultiTileConfigs{
            using type = std::tuple<
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

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{10928,10932,10539,10663,10627,11169,11472,11512,11456,11315,
                    10774};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };


        return std::make_unique<MultiConfigScoreOnlyAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreOnlyAlignerInterface> make_multiconfig_sw_aligner_amd_gfx950_32bit(){
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
