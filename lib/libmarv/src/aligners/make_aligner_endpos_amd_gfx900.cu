#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"


#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_amd_gfx900_32bit(){
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 48>
            >;
            
            static constexpr auto getPerformanceHints(){
                constexpr auto result = std::array{3968,7370,8766,9493,9827,9397,
                    8806,8604,7719,8274,7883,8416,8415,8521,7817,8401,8136,8645,8650,
                    8786,8172,8292,8486,8343,6109,5971,6626,5901,6393,6378,5851};

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
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 12>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 16>,
                // GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 20>
            >;
            
            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{2210,3393,4206,4764,5157,4605,5213,5587,
            //         4873,5440,5762,4922,5522,5864,4987,5560,5918};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };


        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::Gapless);
    }

    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_amd_gfx900_32bit(){
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