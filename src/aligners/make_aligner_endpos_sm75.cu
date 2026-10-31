#include "score_with_endpos_alignment.cuh"
#include "score_with_endpos_alignment_interfaces.cuh"

#include <memory>

#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN


    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_sw_endpos_aligner_sm75_32bit(){
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
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 16>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                SWConfig<21, Approach::unused, Datatype::Float, 512, 8, 20>,
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
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>
                //larger tiles are not supported because shared memory size is too small
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
                SWConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>
                //larger tiles are not supported because shared memory size is too small
            >;

            static constexpr bool hasPerformanceHints(){
                return false;
            }
        };
    
                    
        return std::make_unique<MultiConfigScoreWithEndPosAligner<SingleTileConfigs, MultiTileConfigs>>(AlignmentAlgorithmE::SW);
    }



    std::unique_ptr<MultiConfigScoreWithEndPosAlignerInterface> make_multiconfig_gapless_endpos_aligner_sm75_32bit(){
        //NOTE: copied from a different arch, not tuned.
        struct SingleTileConfigs{
            using type = std::tuple<
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 4>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 8>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 12>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 16>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 2, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 20>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 24>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 2, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 4, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 56>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 4, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 256, 8, 32>,
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
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{};

            //     static_assert(result.size() == std::tuple_size_v<type>);

            //     return result;
            // };

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
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 48>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 52>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 28>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 8, 60>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 32>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 36>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 40>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 16, 44>,
                GaplessConfig<21, Approach::unused, Datatype::Float, 512, 32, 24>
            >;

            // static constexpr auto getPerformanceHints(){
            //     constexpr auto result = std::array{};

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