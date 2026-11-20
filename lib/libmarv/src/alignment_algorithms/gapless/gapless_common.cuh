#ifndef LIBMARV_GAPLESS_COMMON_CUH
#define LIBMARV_GAPLESS_COMMON_CUH


#include "../../namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN

    template<class ScoreType>
    struct UpdateMax_gapless_scoreOnly{
        using DefaultMathOps = MathOps<ScoreType>;

        ScoreType maximum{};

        template<int numRegs>
        __device__
        void operator()(const ScoreType (&values)[numRegs], int /*tileNr*/, int /*row*/){
            constexpr int numPairs = numRegs / 2;
            #pragma unroll
            for(int i = 0; i < numPairs; i++){
                maximum = DefaultMathOps::max3(maximum, values[2*i + 0], values[2*i + 1]);
            }
            if constexpr(numRegs % 2 == 1){
                maximum = DefaultMathOps::max(maximum, values[numRegs-1]);
            }
        }

    };

    template<class ScoreType>
    struct UpdateMax_gapless_scalar_endPos_singleTile{
        using DefaultMathOps = MathOps<ScoreType>;

        ScoreType maximum{};
        int positionOfMaxObserved_reg{-1};
        int positionOfMaxObserved_row{0};

        template<int numRegs>
        __device__
        void operator()(const ScoreType (&values)[numRegs], int /*tileNr*/, int row){
            #pragma unroll
            for(int i = 0; i < numRegs; i++){
                bool oldMaxIsGEQ;
                maximum = DefaultMathOps::max(maximum, values[i], &oldMaxIsGEQ);

                if(!oldMaxIsGEQ){
                    positionOfMaxObserved_reg = i;
                    positionOfMaxObserved_row = row;
                }
            }
        }

    };

    template<class ScoreType>
    struct UpdateMax_gapless_scalar_endPos_multiTile{
        using DefaultMathOps = MathOps<ScoreType>;

        ScoreType maximum{};
        int positionOfMaxObserved_reg{-1};
        int positionOfMaxObserved_row{0};
        int positionOfMaxObserved_tileNr{0};

        template<int numRegs>
        __device__
        void operator()(const ScoreType (&values)[numRegs], int tileNr, int row){
            #pragma unroll
            for(int i = 0; i < numRegs; i++){
                bool oldMaxIsGEQ;
                maximum = DefaultMathOps::max(maximum, values[i], &oldMaxIsGEQ);
                if(!oldMaxIsGEQ){
                    positionOfMaxObserved_reg = i;
                    positionOfMaxObserved_row = row;
                    positionOfMaxObserved_tileNr = tileNr;
                }
            }
        }

    };


    template<class ScoreType>
    struct UpdateMax_gapless_vec2_endPos_singleTile{
        using DefaultMathOps = MathOps<ScoreType>;

        ScoreType maximum{};
        int2 positionOfMaxObserved_regAndLane{-1,-1};
        int2 positionOfMaxObserved_row{0,0};

        template<int numRegs>
        __device__
        void operator()(const ScoreType (&values)[numRegs], int /*tileNr*/, int row){
            #pragma unroll
            for(int i = 0; i < numRegs; i++){
                bool oldMaxIsGEQ_x;
                bool oldMaxIsGEQ_y;
                maximum = DefaultMathOps::max(maximum, values[i], &oldMaxIsGEQ_x, &oldMaxIsGEQ_y);

                if(!oldMaxIsGEQ_x){
                    positionOfMaxObserved_regAndLane.x = 2*i+0;
                    positionOfMaxObserved_row.x = row;
                }
                if(!oldMaxIsGEQ_y){
                    positionOfMaxObserved_regAndLane.y = 2*i+1;
                    positionOfMaxObserved_row.y = row;
                }
            }
        }

    };

    template<class ScoreType>
    struct UpdateMax_gapless_vec2_endPos_multiTile{
        using DefaultMathOps = MathOps<ScoreType>;

        ScoreType maximum{};
        int2 positionOfMaxObserved_regAndLane{-1,-1};
        int2 positionOfMaxObserved_row{0,0};
        int2 positionOfMaxObserved_tileNr{0,0};

        template<int numRegs>
        __device__
        void operator()(const ScoreType (&values)[numRegs], int tileNr, int row){
            #pragma unroll
            for(int i = 0; i < numRegs; i++){
                bool oldMaxIsGEQ_x;
                bool oldMaxIsGEQ_y;
                maximum = DefaultMathOps::max(maximum, values[i], &oldMaxIsGEQ_x, &oldMaxIsGEQ_y);
                if(!oldMaxIsGEQ_x){
                    positionOfMaxObserved_regAndLane.x = 2*i+0;
                    positionOfMaxObserved_row.x = row;
                    positionOfMaxObserved_tileNr.x = tileNr;
                }
                if(!oldMaxIsGEQ_y){
                    positionOfMaxObserved_regAndLane.y = 2*i+1;
                    positionOfMaxObserved_row.y = row;
                    positionOfMaxObserved_tileNr.y = tileNr;
                }
            }
        }

    };


    template<class ScoreType>
    struct UpdateMax_gapless_vec4_endPos_singleTile{
        using DefaultMathOps = MathOps<ScoreType>;

        ScoreType maximum{};
        int4 positionOfMaxObserved_regAndLane{-1,-1,-1,-1};
        int4 positionOfMaxObserved_row{0,0,0,0};

        template<int numRegs>
        __device__
        void operator()(const ScoreType (&values)[numRegs], int /*tileNr*/, int row){
            #pragma unroll
            for(int i = 0; i < numRegs; i++){
                bool oldMaxIsGEQ_x;
                bool oldMaxIsGEQ_y;
                bool oldMaxIsGEQ_z;
                bool oldMaxIsGEQ_w;
                maximum = DefaultMathOps::max(maximum, values[i], &oldMaxIsGEQ_x, &oldMaxIsGEQ_y, &oldMaxIsGEQ_z, &oldMaxIsGEQ_w);
                if(!oldMaxIsGEQ_x){
                    positionOfMaxObserved_regAndLane.x = 4*i+0;
                    positionOfMaxObserved_row.x = row;
                }
                if(!oldMaxIsGEQ_y){
                    positionOfMaxObserved_regAndLane.y = 4*i+1;
                    positionOfMaxObserved_row.y = row;
                }
                if(!oldMaxIsGEQ_z){
                    positionOfMaxObserved_regAndLane.z = 4*i+2;
                    positionOfMaxObserved_row.z = row;
                }
                if(!oldMaxIsGEQ_w){
                    positionOfMaxObserved_regAndLane.w = 4*i+3;
                    positionOfMaxObserved_row.w = row;
                }
            }
        }

    };

    template<class ScoreType>
    struct UpdateMax_gapless_vec4_endPos_multiTile{
        using DefaultMathOps = MathOps<ScoreType>;

        ScoreType maximum{};
        int4 positionOfMaxObserved_regAndLane{-1,-1,-1,-1};
        int4 positionOfMaxObserved_row{0,0,0,0};
        int4 positionOfMaxObserved_tileNr{0,0,0,0};

        template<int numRegs>
        __device__
        void operator()(const ScoreType (&values)[numRegs], int tileNr, int row){
            #pragma unroll
            for(int i = 0; i < numRegs; i++){
                bool oldMaxIsGEQ_x;
                bool oldMaxIsGEQ_y;
                bool oldMaxIsGEQ_z;
                bool oldMaxIsGEQ_w;
                maximum = DefaultMathOps::max(maximum, values[i], &oldMaxIsGEQ_x, &oldMaxIsGEQ_y, &oldMaxIsGEQ_z, &oldMaxIsGEQ_w);
                if(!oldMaxIsGEQ_x){
                    positionOfMaxObserved_regAndLane.x = 4*i+0;
                    positionOfMaxObserved_row.x = row;
                    positionOfMaxObserved_tileNr.x = tileNr;
                }
                if(!oldMaxIsGEQ_y){
                    positionOfMaxObserved_regAndLane.y = 4*i+1;
                    positionOfMaxObserved_row.y = row;
                    positionOfMaxObserved_tileNr.y = tileNr;
                }
                if(!oldMaxIsGEQ_z){
                    positionOfMaxObserved_regAndLane.z = 4*i+2;
                    positionOfMaxObserved_row.z = row;
                    positionOfMaxObserved_tileNr.z = tileNr;
                }
                if(!oldMaxIsGEQ_w){
                    positionOfMaxObserved_regAndLane.w = 4*i+3;
                    positionOfMaxObserved_row.w = row;
                    positionOfMaxObserved_tileNr.w = tileNr;
                }
            }
        }

    };






LIBMARV_NAMESPACE_WITH_NESTING_END



#endif// LIBMARV_GAPLESS_COMMON_CUH