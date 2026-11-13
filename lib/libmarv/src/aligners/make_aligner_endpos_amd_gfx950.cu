#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"


#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_amd_gfx950_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 16, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{3497,6184,6921,7325,7484,7483,7664,7635,7133,7843,
                    7378,7800,7888,7931,8106,8077,8114,8179,8444,8443,
                    8517,8464,8616,8448};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };

        struct MultiTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 20>
            >;

            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{5367,5646,5411,5683};

                static_assert(result.size() == std::tuple_size_v<type>);

                return result;
            };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };


        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_amd_gfx950_32bit(){
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


        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::SW);
    }

LIBMARV_NAMESPACE_END
