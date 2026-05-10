#ifndef PSSM_KERNELS_INT8_CUH
#define PSSM_KERNELS_INT8_CUH

#include <map>

#include "pssm.cuh"
#include "convert.cuh"
#include "mathops.cuh"
#include "util.cuh"

#include "custom_score_types.cuh"
#include "ptx_wrappers.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;


#define USE_IMPROVED_SMEM_INT8


#if 0
#define PSSM_GAPLESS_INT8_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X \
    X(4,4)

#define PSSM_GAPLESS_INT8_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X \
    X(4,4)
#endif


#if 1
#define PSSM_GAPLESS_INT8_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X \
    X(4,4) X(4,8) X(4,12) X(4,16) X(4,20) X(4,24) X(4,28) X(4,32) \
    X(4,36) X(4,40) X(4,44) X(4,48) X(4,52) X(4,56) X(4,60) X(4,64) \
    X(8,36) X(8,40) X(8,44) X(8,48) X(8,52) X(8,56) X(8,60) X(8,64) \
    X(16,36) X(16,40) X(16,44) X(16,48) X(16,52) X(16,56) X(16,60) X(16,64)

#define PSSM_GAPLESS_INT8_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X \
    X(4,4) X(4,8) X(4,12) X(4,16) X(4,20) X(4,24) X(4,28) X(4,32) \
    X(4,36) X(4,40) X(4,44) X(4,48) X(4,52) X(4,56) X(4,60) X(4,64) \
    X(8,36) X(8,40) X(8,44) X(8,48) X(8,52) X(8,56) X(8,60) X(8,64) \
    X(16,36) X(16,40) X(16,44) X(16,48) X(16,52) X(16,56) X(16,60) X(16,64)
#endif

namespace cudasw4{




namespace uint8x4{

    __device__ __forceinline__
    unsigned int max3_u8x4(unsigned int a, unsigned int b, unsigned int c){
        #ifdef HAS_BLACKWELL_INT8_PTX
        return ptx_max_u8x4(ptx_max_u8x4(a,b), c);
        #else
        return 0;
        #endif
    }

    __device__ __forceinline__
    unsigned int max3_s8x4(unsigned int a, unsigned int b, unsigned int c){
        #ifdef HAS_BLACKWELL_INT8_PTX
        return ptx_max_s8x4(ptx_max_s8x4(a,b), c);
        #else
        return 0;
        #endif
    }


    

    template<class ScoreType> struct ScalarScoreType{};
    // template<> struct ScalarScoreType<ScoreType_s8x4>{ using type = cuda::std::int8_t; };
    template<> struct ScalarScoreType<ScoreType_u8x4>{ using type = cuda::std::uint8_t; };


 
    template<class ScoreType, int numRegs, class Group, class SharedPSSM, class SmemIndexCalculator>
    struct GaplessPSSMState{
        static_assert(
            std::is_same_v<ScoreType, ScoreType_u8x4> 
            // || std::is_same_v<ScoreType, ScoreType_s8x4>
        );

        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOps = MathOps<ScoreType>;

        ScoreType substitutionScoreBias;
        ScoreType penalty_here_array[numRegs];
        ScoreType maximum{}; //0
        ScoreType penalty_diag{}; //0
        SharedPSSM& shared_strided_PSSM;
        Group& group;

        __device__
        GaplessPSSMState(SharedPSSM& s, Group& g, ScoreType b) : substitutionScoreBias(b), shared_strided_PSSM(s), group(g) {

        }

        __device__
        void resetScores(){
            #pragma unroll
            for(int i = 0; i < numRegs; i++){
                penalty_here_array[i] = ScoreType{};
            }
            
            penalty_diag = ScoreType{};
        }

        __device__
        void resetMaximum(){
            maximum = ScoreType{};
        }

        #if 1
        //mathops
        __device__
        void relax(int subject_letter, int subjectPos = 0){
            SmemIndexCalculator smemIndexCalculator;

            ScoreType score2;
            ScoreType penalty_temp0;
            ScoreType penalty_temp1;

            const auto* row = &shared_strided_PSSM.data[subject_letter][0];

            float4 foo = *((float4*)&row[smemIndexCalculator.getIndex(0)]);
            memcpy(&score2, &foo.x, sizeof(ScoreType));
            //if(threadIdx.x == 0)
            // {
            //     int8_t tmp[16];
            //     memcpy(&tmp[0], &foo, 16);
            //     //printf("%d %d %d %d ", s.x(), s.y(), s.z(), s.w());

            //     //printf("letter %d, loaded %d %d %d %d\n", subject_letter, score2.x(), score2.y(), score2.z(), score2.w());
            //     printf("thread %d, letter %d, loaded (%d %d %d %d), (%d %d %d %d), (%d %d %d %d), (%d %d %d %d)\n", 
            //         threadIdx.x, subject_letter, 
            //         tmp[0], tmp[1], tmp[2], tmp[3], 
            //         tmp[4], tmp[5], tmp[6], tmp[7], 
            //         tmp[8], tmp[9], tmp[10], tmp[11], 
            //         tmp[12], tmp[13], tmp[14], tmp[15]);
            // }
            penalty_temp0 = penalty_here_array[0];
            penalty_here_array[0] = MathOps::sub_sat(MathOps::add_sat(penalty_diag, score2), substitutionScoreBias);

            memcpy(&score2, &foo.y, sizeof(ScoreType));
            penalty_temp1 = penalty_here_array[1];
            penalty_here_array[1] = MathOps::sub_sat(MathOps::add_sat(penalty_temp0, score2), substitutionScoreBias);
            maximum = MathOps::max3(maximum, penalty_here_array[1], penalty_here_array[0]);

            memcpy(&score2, &foo.z, sizeof(ScoreType));
            penalty_temp0 = penalty_here_array[2];
            penalty_here_array[2] = MathOps::sub_sat(MathOps::add_sat(penalty_temp1, score2), substitutionScoreBias);

            memcpy(&score2, &foo.w, sizeof(ScoreType));
            penalty_temp1 = penalty_here_array[3];
            penalty_here_array[3] = MathOps::sub_sat(MathOps::add_sat(penalty_temp0, score2), substitutionScoreBias);
            maximum = MathOps::max3(maximum, penalty_here_array[3], penalty_here_array[2]);


            #pragma unroll
            for (int i=1; i<numRegs/4; i++) {
                foo = *((float4*)&row[smemIndexCalculator.getIndex(i)]);
                memcpy(&score2, &foo.x, sizeof(ScoreType));
                penalty_temp0 = penalty_here_array[4*i];
                penalty_here_array[4*i] = MathOps::sub_sat(MathOps::add_sat(penalty_temp1, score2), substitutionScoreBias);

                memcpy(&score2, &foo.y, sizeof(ScoreType));
                penalty_temp1 = penalty_here_array[4*i+1];
                penalty_here_array[4*i+1] = MathOps::sub_sat(MathOps::add_sat(penalty_temp0, score2), substitutionScoreBias);
                maximum = MathOps::max3(maximum, penalty_here_array[4*i+1], penalty_here_array[4*i]);

                memcpy(&score2, &foo.z, sizeof(ScoreType));
                penalty_temp0 = penalty_here_array[4*i+2];
                penalty_here_array[4*i+2] = MathOps::sub_sat(MathOps::add_sat(penalty_temp1, score2), substitutionScoreBias);

                memcpy(&score2, &foo.w, sizeof(ScoreType));
                penalty_temp1 = penalty_here_array[4*i+3];
                penalty_here_array[4*i+3] = MathOps::sub_sat(MathOps::add_sat(penalty_temp0, score2), substitutionScoreBias);
                maximum = MathOps::max3(maximum, penalty_here_array[4*i+3], penalty_here_array[4*i+2]);
            }

            // {
            //     int8_t tmp[16];
            //     memcpy(&tmp[0], &penalty_here_array, 16);
            //     //printf("%d %d %d %d ", s.x(), s.y(), s.z(), s.w());

            //     //printf("letter %d, loaded %d %d %d %d\n", subject_letter, score2.x(), score2.y(), score2.z(), score2.w());
            //     printf("thread %d, letter %d, computed (%d %d %d %d), (%d %d %d %d), (%d %d %d %d), (%d %d %d %d)\n", 
            //         threadIdx.x, subject_letter, 
            //         tmp[0], tmp[1], tmp[2], tmp[3], 
            //         tmp[4], tmp[5], tmp[6], tmp[7], 
            //         tmp[8], tmp[9], tmp[10], tmp[11], 
            //         tmp[12], tmp[13], tmp[14], tmp[15]);
            // }
        };
        #endif

        #if 0
        //no mathops
        __device__
        void relax(int subject_letter, int subjectPos = 0){
            static_assert(std::is_same_v<ScoreType, ScoreType_u8x4>);
            SmemIndexCalculator smemIndexCalculator;

            ScoreType score2;
            ScoreType penalty_temp0;
            ScoreType penalty_temp1;

            const auto* row = &shared_strided_PSSM.data[subject_letter][0];

            float4 foo = *((float4*)&row[smemIndexCalculator.getIndex(0)]);
            memcpy(&score2, &foo.x, sizeof(ScoreType));
            penalty_temp0 = penalty_here_array[0];
            penalty_here_array[0] = ptx_add_sat_u8x4(penalty_diag, score2);
            penalty_here_array[0] = ptx_sub_sat_u8x4(penalty_here_array[0], substitutionScoreBias);

            memcpy(&score2, &foo.y, sizeof(ScoreType));
            penalty_temp1 = penalty_here_array[1];
            penalty_here_array[1] = ptx_add_sat_u8x4(penalty_temp0, score2);
            penalty_here_array[1] = ptx_sub_sat_u8x4(penalty_here_array[1], substitutionScoreBias);
            maximum = max3_u8x4(maximum, penalty_here_array[1], penalty_here_array[0]);

            memcpy(&score2, &foo.z, sizeof(ScoreType));
            penalty_temp0 = penalty_here_array[2];
            penalty_here_array[2] = ptx_add_sat_u8x4(penalty_temp1, score2);
            penalty_here_array[2] = ptx_sub_sat_u8x4(penalty_here_array[2], substitutionScoreBias);

            memcpy(&score2, &foo.w, sizeof(ScoreType));
            penalty_temp1 = penalty_here_array[3];
            penalty_here_array[3] = ptx_add_sat_u8x4(penalty_temp0, score2);
            penalty_here_array[3] = ptx_sub_sat_u8x4(penalty_here_array[3], substitutionScoreBias);
            maximum = max3_u8x4(maximum, penalty_here_array[3], penalty_here_array[2]);


            #pragma unroll
            for (int i=1; i<numRegs/4; i++) {
                foo = *((float4*)&row[smemIndexCalculator.getIndex(i)]);
                memcpy(&score2, &foo.x, sizeof(ScoreType));
                penalty_temp0 = penalty_here_array[4*i];
                penalty_here_array[4*i] = ptx_add_sat_u8x4(penalty_temp1, score2);
                penalty_here_array[4*i] = ptx_sub_sat_u8x4(penalty_here_array[4*i], substitutionScoreBias);

                memcpy(&score2, &foo.y, sizeof(ScoreType));
                penalty_temp1 = penalty_here_array[4*i+1];
                penalty_here_array[4*i+1] = ptx_add_sat_u8x4(penalty_temp0, score2);
                penalty_here_array[4*i+1] = ptx_sub_sat_u8x4(penalty_here_array[4*i+1], substitutionScoreBias);
                maximum = max3_u8x4(maximum, penalty_here_array[4*i+1], penalty_here_array[4*i]);

                memcpy(&score2, &foo.z, sizeof(ScoreType));
                penalty_temp0 = penalty_here_array[4*i+2];
                penalty_here_array[4*i+2] = ptx_add_sat_u8x4(penalty_temp1, score2);
                penalty_here_array[4*i+2] = ptx_sub_sat_u8x4(penalty_here_array[4*i+2], substitutionScoreBias);

                memcpy(&score2, &foo.w, sizeof(ScoreType));
                penalty_temp1 = penalty_here_array[4*i+3];
                penalty_here_array[4*i+3] = ptx_add_sat_u8x4(penalty_temp0, score2);
                penalty_here_array[4*i+3] = ptx_sub_sat_u8x4(penalty_here_array[4*i+3], substitutionScoreBias);
                maximum = max3_u8x4(maximum, penalty_here_array[4*i+3], penalty_here_array[4*i+2]);
            }
        };
        #endif



        __device__
        void shuffleScores(const Scalar& border_in){
            penalty_diag = group.shfl_up(penalty_here_array[numRegs-1], 1);
            const ScoreType penalty_temp0 = group.shfl_down(penalty_here_array[numRegs-1], group.size()-1);

            if (group.thread_rank() == 0) {

                // int8_t tmp[4];

                // memcpy(&tmp[0], &penalty_temp0, 4);
                // printf("border in, %d, penalty_temp0 %d %d %d %d\n",
                //     border_in,
                //     tmp[0], tmp[1], tmp[2], tmp[3]
                // );
                
                // memcpy(&tmp[0], &penalty_diag, 4);
                // printf("penalty_diag before %d %d %d %d\n",
                //     tmp[0], tmp[1], tmp[2], tmp[3]
                // );

                penalty_diag = ScoreType(border_in, penalty_temp0.x(), penalty_temp0.y(), penalty_temp0.z());
                
                // memcpy(&tmp[0], &penalty_diag, 4);
                // printf("penalty_diag after %d %d %d %d\n",
                //     tmp[0], tmp[1], tmp[2], tmp[3]
                // );
            }
        }

        __device__
        void stepSingleTile(int subject_letter, int subjectPos){
            relax(subject_letter, subjectPos);
            shuffleScores(Scalar{});
        }

        __device__
        void stepFirstTile(int subject_letter, Scalar& border_out){
            relax(subject_letter);
            shuffleScores(Scalar{});
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1].w();
            }

            // if(group.thread_rank() == group.size() - 1){
            //     printf("thread %d, penalty_out = %d\n", threadIdx.x, border_out);
            // }
        }

        __device__
        void stepIntermediateTile(int subject_letter, const Scalar& border_in, Scalar& border_out){
            relax(subject_letter);
            shuffleScores(border_in);
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1].w();
            }
            // if(group.thread_rank() == group.size() - 1){
            //     printf("thread %d, penalty_out = %d\n", threadIdx.x, border_out);
            // }
        }

        __device__
        void stepLastTile(int subject_letter, const Scalar& border_in){
            relax(subject_letter);
            shuffleScores(border_in);
        }

        __device__
        void reduceMaximumScore(){
            maximum = MathOps::reduce_max(group, maximum);
            // maximum = cooperative_groups::reduce(group, maximum.raw, [](const auto& l, const auto& r){return ptx_max_u8x4(l,r);});
        }
    };


    /*
    PSSM kernel for a query of max length (4 * group_size * numRegs)
    */
    template<
        class ScoreType, 
        int blocksize, 
        int group_size, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator
    >
    __global__
    __launch_bounds__(512,1)
    void GaplessFilter_strided_PSSM_singletile_uint8x4_kernel(
        __grid_constant__ const char * const devChars,
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const size_t* const devOffsets,
        __grid_constant__ const SequenceLengthT* const devLengths,
        __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
        __grid_constant__ const int numSelected,
        __grid_constant__ const SequenceLengthT queryLength,
        __grid_constant__ const PSSM_2D_View<ScoreType> strided_PSSM,
        __grid_constant__ const ScoreType substitutionScoreBias // values in strided_PSSM are computed from original PSSM + bias

    ) {
        #if !defined(HAS_BLACKWELL_INT8_PTX)
            return;
        #endif

        static_assert(numRegs % 4 == 0);
        static_assert(blocksize % group_size == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % group_size == 0);

        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM_INT8
        constexpr int numColumnsPSSM = std::max(group_size,8) * numRegs;
        #else
        constexpr int numColumnsPSSM = group_size * numRegs;
        #endif

        using SharedPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;
        using MathOps = MathOps<ScoreType>;

        extern  __shared__ char externalSmem[];

        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);


        auto group = cg::tiled_partition<group_size>(cg::this_thread_block());
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / group_size;
        //const int numGroupsInGrid = (blockDim.x * gridDim.x) / group_size;

        #ifdef USE_IMPROVED_SMEM_INT8
        using SmemIndexCalculator = typename std::conditional<
            group_size == 4, 
            SmemIndexCalculator<group_size, 2>,
            SmemIndexCalculator<group_size, 1>
        >::type;
        #else
        using SmemIndexCalculator = SmemIndexCalculator<group_size, 1>;
        #endif

        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator>;
        State state(shared_strided_PSSM, group, substitutionScoreBias);

        auto load_PSSM_single = [&]() {
            for (int i=threadIdx.x; i<21*group_size*numRegs; i+=blockDim.x) {
                const int letter = i/(group_size*numRegs);
                const int col = i%(group_size*numRegs);
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][col];
            }
            __syncthreads();
        };

        auto load_PSSM_double = [&]() {
            for (int i=threadIdx.x; i<21*group_size*numRegs; i+=blockDim.x) {
                const int letter = i/(group_size*numRegs);
                const int col = i%(group_size*numRegs);
                auto value = strided_PSSM[letter][col];

                const int float4Index = col / 4;
                const int offsetWithinFloat4 = col % 4;

                const int ithChunkOfFour = float4Index / group_size;
                const int float4PositionInChunkOfFour = float4Index % group_size;

                const int outputFloat4Index0 = (ithChunkOfFour*2*group_size + 0*group_size) + float4PositionInChunkOfFour;
                const int outputFloat4Index1 = (ithChunkOfFour*2*group_size + 1*group_size) + float4PositionInChunkOfFour;

                shared_strided_PSSM.data[letter][4*outputFloat4Index0 + offsetWithinFloat4] = value;
                shared_strided_PSSM.data[letter][4*outputFloat4Index1 + offsetWithinFloat4] = value;
            }
            __syncthreads();

            // for (int i=threadIdx.x; i<21*group_size*numRegs*2; i+=blockDim.x) {
            //     const int letter = i/(group_size*numRegs*2);
            //     const int outputCol = i%(group_size*numRegs*2);

            //     const int outputFloat4Index = outputCol / 4;
            //     const int offsetWithinFloat4 = outputCol % 4;
            //     const int inputIthChunkOfFour = outputFloat4Index / (2*group_size);
            //     const int remainder1 = outputFloat4Index % (2*group_size); // 0*group_size + float4PositionInChunkOfFour or 1*group_size + float4PositionInChunkOfFour
            //     const int float4PositionInChunkOfFour = remainder1 % group_size;
            //     const int float4Index = 4*inputIthChunkOfFour+float4PositionInChunkOfFour;

            //     // const int inputIthChunkOfFour = outputCol / 2*group_size;
            //     // const int remainder1 = outputCol % 2*group_size; // 0*group_size + float4PositionInChunkOfFour or 1*group_size + float4PositionInChunkOfFour
            //     // const int float4PositionInChunkOfFour = remainder1 % group_size;
            //     // const int float4Index = 4*inputIthChunkOfFour+float4PositionInChunkOfFour;
            //     // const int offsetWithinFloat4 = col % 4;

            //     const int inputCol = 4*float4Index + offsetWithinFloat4;
            //     shared_strided_PSSM.data[letter][outputCol] = strided_PSSM[letter][inputCol];
            // }
            // __syncthreads();
        };

        auto load_PSSM = [&](){
            if constexpr(SmemIndexCalculator::factor == 2){
                load_PSSM_double();
            }else{
                load_PSSM_single();
            }

            // __syncthreads();
            // if(threadIdx.x == 0){
            //     printf("smem pssm\n");
            //     for(int r = 0; r < 21; r++){
            //         for(int c = 0; c < numColumnsPSSM; c++){
            //             auto s = shared_strided_PSSM.data[r][c];
            //             int8_t tmp[4];
            //             memcpy(&tmp[0], &s, 4);
            //             //printf("%d %d %d %d ", s.x(), s.y(), s.z(), s.w());
            //             printf("%d %d %d %d ", tmp[0], tmp[1], tmp[2], tmp[3]);
            //         }
            //         printf("\n");
            //     }
            // }
            // __syncthreads();
        };

        const char4* subjectAsChar4;
        char4 new_subject_letter4;

        auto makeCaseInsensitive4 = [](char4 encoded4){
            unsigned int asUint;
            memcpy(&asUint, &encoded4, sizeof(unsigned int));

            if constexpr(subjectIsCaseSensitive){
                // asUint = CaseSensitive_to_CaseInsensitive{}(asUint);
                asUint = ClampToInvalid{}(asUint);
            }

            memcpy(&encoded4, &asUint, sizeof(unsigned int));
            return encoded4;
        };

        load_PSSM();


        //for(int alignmentId = idOfGroupInGrid; alignmentId < numSelected; alignmentId += numGroupsInGrid){
        const int alignmentId = idOfGroupInGrid;
        if(alignmentId < numSelected){
            const auto subjectId = d_positions_of_selected_lengths[alignmentId];
            const SequenceLengthT subjectLength = devLengths[subjectId];
            const size_t base_S = devOffsets[subjectId]-devOffsets[0];
            // const size_t base_S = devOffsets[subjectId]; //for debugging

            state.resetScores();
            state.resetMaximum();

            subjectAsChar4 = reinterpret_cast<const char4*>(&devChars[base_S]);


            int k;
            for (k=0; k<subjectLength-3; k+=4) {
                new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                state.stepSingleTile(new_subject_letter4.x, k);
                state.stepSingleTile(new_subject_letter4.y, k+1);
                state.stepSingleTile(new_subject_letter4.z, k+2);
                state.stepSingleTile(new_subject_letter4.w, k+3);
            }

            if (subjectLength%4 >= 1) {
                new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                state.stepSingleTile(new_subject_letter4.x, k);
            }

            if (subjectLength%4 >= 2) {
                state.stepSingleTile(new_subject_letter4.y, k+1);
            }

            if (subjectLength%4 >= 3) {
                state.stepSingleTile(new_subject_letter4.z, k+2);
            }

            state.reduceMaximumScore();

            const unsigned int overall_max = max(state.maximum.x(), max(state.maximum.y(), max(state.maximum.z(), state.maximum.w())));


            // const unsigned int overall_max = max(state.maximum.x(), max(state.maximum.y(), max(state.maximum.z(), state.maximum.w())))
            //     + state.bestSubjectPos.x
            //     + state.bestSubjectPos.y
            //     + state.bestSubjectPos.z
            //     + state.bestSubjectPos.w
            //     + state.bestPssmColumn.x
            //     + state.bestPssmColumn.y
            //     + state.bestPssmColumn.z
            //     + state.bestPssmColumn.w;

            if(group.thread_rank() == 0){
                // printf("overall_max = %u\n", overall_max);
                devAlignmentScores[alignmentId] = overall_max;
            }
        }
    }



    /*
    PSSM kernel for a query of max length (4 * group_size * numRegs)
    */
    template<
        class ScoreType,
        int blocksize, 
        int group_size, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator
    >
    void call_GaplessFilter_strided_PSSM_singletile_uint8x4_kernel(
        const char * const devChars,
        ScoreOutputIterator const devAlignmentScores,
        const size_t* const devOffsets,
        const SequenceLengthT* const devLengths,
        PositionsIterator const d_positions_of_selected_lengths,
        const int numSelected,
        const SequenceLengthT queryLength,
        const PSSM_2D_View<ScoreType>& strided_PSSM,
        const int substitutionScoreBias, // values in strided_PSSM are computed from original PSSM + bias
        cudaStream_t stream
    ){
        if(substitutionScoreBias > 255){
            throw std::runtime_error("substitutionScoreBias cannot be represented with 8 bits.");
        }
        const ScoreType substitutionScoreBias_vec4(
            substitutionScoreBias,
            substitutionScoreBias,
            substitutionScoreBias,
            substitutionScoreBias
        );

        constexpr int groupsPerBlock = blocksize / group_size;
        constexpr int alignmentsPerGroup = 1;
        constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
        // std::cout << "blocksize " << blocksize << ", group_size " << group_size 
        //     << ", alignmentsPerBlock " << alignmentsPerBlock << ", numSelected " << numSelected << "\n";

        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM_INT8
        constexpr int numColumnsPSSM = std::max(group_size,8) * numRegs;
        #else
        constexpr int numColumnsPSSM = group_size * numRegs;
        #endif
        using SharedPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;
        
        int smem = sizeof(SharedPSSM);
        auto kernel = GaplessFilter_strided_PSSM_singletile_uint8x4_kernel<
            ScoreType,
            blocksize, 
            group_size, 
            numRegs, 
            subjectIsCaseSensitive,
            ScoreOutputIterator>;

        auto setSmemKernelAttribute = [&](){
            static std::map<int, bool> isSet;
            if(smem > 48*1024){
                int deviceId;
                cudaGetDevice(&deviceId); CUERR;
                if(!isSet[deviceId]){
                    cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, smem); CUERR;
                    isSet[deviceId] = true;
                }
            }
        };
        setSmemKernelAttribute();

        dim3 grid = (numSelected + alignmentsPerBlock - 1) / alignmentsPerBlock;

        kernel<<<grid, blocksize, smem, stream>>>(
            devChars,
            devAlignmentScores,
            devOffsets,
            devLengths,
            d_positions_of_selected_lengths,
            numSelected,       
            queryLength,
            strided_PSSM,
            substitutionScoreBias_vec4
        ); CUERR;
    }

    template<class ScoreType, int blocksize, class ScoreOutputIterator>
    void call_GaplessFilter_strided_PSSM_singletile_uint8x4_kernel(
        int group_size,
        int numRegs,
        const char * const devChars,
        ScoreOutputIterator const devAlignmentScores,
        const size_t* const devOffsets,
        const SequenceLengthT* const devLengths,
        PositionsIterator const d_positions_of_selected_lengths,
        const int numSelected,
        const SequenceLengthT queryLength,
        const PSSM_2D_View<ScoreType>& strided_PSSM,
        const int substitutionScoreBias, // values in strided_PSSM are computed from original PSSM + bias
        cudaStream_t stream
    ){
        constexpr bool subjectIsCaseSensitive = true;

        #define X(g,r) \
            if(group_size == g && numRegs == r){ \
                call_GaplessFilter_strided_PSSM_singletile_uint8x4_kernel<ScoreType, blocksize,g,r,subjectIsCaseSensitive>( \
                    devChars, devAlignmentScores, devOffsets, devLengths, d_positions_of_selected_lengths, \
                    numSelected, queryLength, strided_PSSM, substitutionScoreBias, stream \
                ); \
            } else 

            PSSM_GAPLESS_INT8_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X
        { throw std::runtime_error("invalid groupsize/numregs config");}

        #undef X
    }


    #define ScoreOutputIterator TopNMaximaArray
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void uint8x4::call_GaplessFilter_strided_PSSM_singletile_uint8x4_kernel<ScoreType_u8x4, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<ScoreType_u8x4>&, \
            const int, \
            cudaStream_t \
        );

        PSSM_GAPLESS_INT8_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef ScoreOutputIterator


    /*
    PSSM kernel for arbitrary query length
    */
    template<
        class ScoreType, 
        int blocksize, 
        int group_size, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator
    >
    __global__
    __launch_bounds__(512,1)
    void GaplessFilter_strided_PSSM_multitile_uint8x4_kernel(
        __grid_constant__ const char * const devChars,
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const size_t* const devOffsets,
        __grid_constant__ const SequenceLengthT* const devLengths,
        __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
        __grid_constant__ const int numSelected,
        __grid_constant__ const SequenceLengthT queryLength,
        __grid_constant__ const PSSM_2D_View<ScoreType> strided_PSSM,
        __grid_constant__ const ScoreType substitutionScoreBias, // values in strided_PSSM are computed from original PSSM + bias
        __grid_constant__ std::uint32_t* const multiTileTempStorage,
        __grid_constant__ const size_t tempStorageElementsPerGroup
    ) {
        #if !defined(HAS_BLACKWELL_INT8_PTX)
            return;
        #endif

        static_assert(numRegs % 4 == 0);
        static_assert(blocksize % group_size == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % group_size == 0);

        extern  __shared__ char externalSmem[];
        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM_INT8
        constexpr int numColumnsPSSM = std::max(group_size,8) * numRegs;
        #else
        constexpr int numColumnsPSSM = group_size * numRegs;
        #endif
        using SharedPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;

        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);

        using MathOps = MathOps<ScoreType>;
        using Scalar = typename ScalarScoreType<ScoreType>::type;

        auto group = cg::tiled_partition<group_size>(cg::this_thread_block());
        const int numGroupsInBlock = blockDim.x / group_size;
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / group_size;
        const int numGroupsInGrid = (blockDim.x * gridDim.x) / group_size;
        
        const size_t groupTempStorageOffset = idOfGroupInGrid * tempStorageElementsPerGroup;
        std::uint32_t* const groupTempStorage = multiTileTempStorage + groupTempStorageOffset;
        
        const int numTiles = SDIV(queryLength, 4 * group_size * numRegs);

        #ifdef USE_IMPROVED_SMEM_INT8
        using SmemIndexCalculator = typename std::conditional<
            group_size == 4, 
            SmemIndexCalculator<group_size, 2>,
            SmemIndexCalculator<group_size, 1>
        >::type;
        #else
        using SmemIndexCalculator = SmemIndexCalculator<group_size, 1>;
        #endif

        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator>;
        State state(shared_strided_PSSM, group, substitutionScoreBias);


        alignas(8) Scalar penalty_in[4]{};
        alignas(8) Scalar penalty_out[4]{};

        auto load_PSSM_single = [&](int tileNr) {
            const int columnOffset = tileNr * group_size * numRegs;
            __syncthreads(); //wait for all groups before overwriting pssm

            for (int i=threadIdx.x; i<21*group_size*numRegs; i+=blockDim.x) {
                int letter = i/(group_size*numRegs);
                int col = i%(group_size*numRegs);
                //shared_strided_PSSM.data[letter][col] = strided_PSSM_1d[i];
                //shared_strided_PSSM.data[letter][col] = strided_PSSM.data[i];
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][columnOffset + col];
            }
            __syncthreads();
        };

        auto load_PSSM_double = [&](int tileNr) {
            const int columnOffset = tileNr * group_size * numRegs;
            __syncthreads(); //wait for all groups before overwriting pssm

            for (int i=threadIdx.x; i<21*group_size*numRegs; i+=blockDim.x) {
                const int letter = i/(group_size*numRegs);
                const int col = i%(group_size*numRegs);
                auto value = strided_PSSM[letter][columnOffset + col];

                const int float4Index = col / 4;
                const int offsetWithinFloat4 = col % 4;

                const int ithChunkOfFour = float4Index / group_size;
                const int float4PositionInChunkOfFour = float4Index % group_size;

                const int outputFloat4Index0 = (ithChunkOfFour*2*group_size + 0*group_size) + float4PositionInChunkOfFour;
                const int outputFloat4Index1 = (ithChunkOfFour*2*group_size + 1*group_size) + float4PositionInChunkOfFour;

                shared_strided_PSSM.data[letter][4*outputFloat4Index0 + offsetWithinFloat4] = value;
                shared_strided_PSSM.data[letter][4*outputFloat4Index1 + offsetWithinFloat4] = value;
            }
            __syncthreads();
        };

        auto load_PSSM = [&](int tileNr){
            if constexpr(SmemIndexCalculator::factor == 2){
                load_PSSM_double(tileNr);
            }else{
                load_PSSM_single(tileNr);
            }

            // __syncthreads();
            // if(threadIdx.x == 0){
            //     printf("smem pssm\n");
            //     for(int r = 0; r < 21; r++){
            //         for(int c = 0; c < numColumnsPSSM; c++){
            //             auto s = shared_strided_PSSM.data[r][c];
            //             int8_t tmp[4];
            //             memcpy(&tmp[0], &s, 4);
            //             //printf("%d %d %d %d ", s.x(), s.y(), s.z(), s.w());
            //             printf("%d %d %d %d ", tmp[0], tmp[1], tmp[2], tmp[3]);
            //         }
            //         printf("\n");
            //     }
            // }
            // __syncthreads();
        };

        char4 new_subject_letter4;

        auto makeCaseInsensitive4 = [](char4 encoded4){
            unsigned int asUint;
            memcpy(&asUint, &encoded4, sizeof(unsigned int));

            if constexpr(subjectIsCaseSensitive){
                //asUint = CaseSensitive_to_CaseInsensitive{}(asUint);
                asUint = ClampToInvalid{}(asUint);
            }

            memcpy(&encoded4, &asUint, sizeof(unsigned int));
            return encoded4;
        };


        //need to round up to blocks because loading pssm is a block-wide operation
        const int numSelectedRoundedUp = SDIV(numSelected, numGroupsInBlock) * numGroupsInBlock;

        for(int alignmentId = idOfGroupInGrid; alignmentId < numSelectedRoundedUp; alignmentId += numGroupsInGrid){

            size_t subjectId;
            SequenceLengthT subjectLength;
            size_t base_S;
            const char4* subjectAsChar4;

            //first tile
            {
                /* 
                    -----------------------
                    Process tile 0
                    ----------------------- 
                */

                //load pssm for tile 0. blockwide operation
                load_PSSM(0);

                if(alignmentId < numSelected){
                    subjectId = d_positions_of_selected_lengths[alignmentId];
                    subjectLength = devLengths[subjectId];
                    base_S = devOffsets[subjectId]-devOffsets[0];
                    // base_S = devOffsets[subjectId]; //for debugging

                    // if(threadIdx.x == 0){
                    //     printf("in kernel, subject:\n");
                    //     for(int x = 0; x < subjectLength; x++){
                    //         printf("%d", (devChars[base_S+x]));
                    //     }
                    //     printf("\n");
                    // }

                    state.resetScores();
                    state.resetMaximum();
                    subjectAsChar4 = reinterpret_cast<const char4*>(&devChars[base_S]);

                    int k;

                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {

                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);

                        state.stepFirstTile(new_subject_letter4.x, penalty_out[0]);
                        state.stepFirstTile(new_subject_letter4.y, penalty_out[1]);
                        state.stepFirstTile(new_subject_letter4.z, penalty_out[2]);
                        state.stepFirstTile(new_subject_letter4.w, penalty_out[3]);
                        
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((std::uint32_t*)&penalty_out[0]);
                        }
                    }

                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        state.stepFirstTile(new_subject_letter4.x, penalty_out[0]);
                        // if(group.thread_rank() == group.size() - 1){
                        //     printf("thread %d, penalty_out[0] = %d\n", threadIdx.x, penalty_out[0]);
                        // }
                    }

                    if (subjectLength%4 >= 2) {
                        state.stepFirstTile(new_subject_letter4.y, penalty_out[1]);
                    }

                    if (subjectLength%4 >= 3) {
                        state.stepFirstTile(new_subject_letter4.z, penalty_out[2]);
                    }

                    //if there were remaining rows, update temp storage
                    if(subjectLength % 4 > 0){
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((std::uint32_t*)&penalty_out[0]);
                        }
                    }
                }
            }

            //intermediate tiles
            for(int tileNr = 1; tileNr < numTiles - 1; tileNr++){
                /* 
                    -----------------------
                    Process tile tileNr
                    ----------------------- 
                */

                //load pssm for tile tileNr. blockwide operation
                load_PSSM(tileNr);

                if(alignmentId < numSelected){    
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((std::uint32_t*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepIntermediateTile(new_subject_letter4.x, penalty_in[0], penalty_out[0]);
                        state.stepIntermediateTile(new_subject_letter4.y, penalty_in[1], penalty_out[1]);
                        state.stepIntermediateTile(new_subject_letter4.z, penalty_in[2], penalty_out[2]);
                        state.stepIntermediateTile(new_subject_letter4.w, penalty_in[3], penalty_out[3]);
            
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((std::uint32_t*)&penalty_out[0]);
                        }
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((std::uint32_t*)&penalty_in[0]) = groupTempStorage[k/4];
                        }
                        state.stepIntermediateTile(new_subject_letter4.x, penalty_in[0], penalty_out[0]);
                    }
        
                    if (subjectLength%4 >= 2) {
                        state.stepIntermediateTile(new_subject_letter4.y, penalty_in[1], penalty_out[1]);
                    }
        
                    if (subjectLength%4 >= 3) {
                        state.stepIntermediateTile(new_subject_letter4.z, penalty_in[2], penalty_out[2]);
                    }
        
                    //if there were remaining rows, update temp storage
                    if(subjectLength % 4 > 0){
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((std::uint32_t*)&penalty_out[0]);
                        }
                    }
                }
            }

            //last tile
            if(numTiles > 1){
                /* 
                    -----------------------
                    Process last tile (numTiles-1)
                    ----------------------- 
                */

                //load pssm for tile (numTiles-1). blockwide operation
                load_PSSM(numTiles-1);

                if(alignmentId < numSelected){
        
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((std::uint32_t*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepLastTile(new_subject_letter4.x, penalty_in[0]);
                        state.stepLastTile(new_subject_letter4.y, penalty_in[1]);
                        state.stepLastTile(new_subject_letter4.z, penalty_in[2]);
                        state.stepLastTile(new_subject_letter4.w, penalty_in[3]);
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((std::uint32_t*)&penalty_in[0]) = groupTempStorage[k/4];
                        }
                        state.stepLastTile(new_subject_letter4.x, penalty_in[0]);
                    }
        
                    if (subjectLength%4 >= 2) {
                        state.stepLastTile(new_subject_letter4.y, penalty_in[1]);
                    }
        
                    if (subjectLength%4 >= 3) {
                        state.stepLastTile(new_subject_letter4.z, penalty_in[2]);
                    }
                }
            }

            if(alignmentId < numSelected){
                state.reduceMaximumScore();
                const unsigned int overall_max = max(state.maximum.x(), max(state.maximum.y(), max(state.maximum.z(), state.maximum.w())));

                if(group.thread_rank() == 0){
                    devAlignmentScores[alignmentId] = overall_max;
                }
            }
        }

    }



    /*
    PSSM kernel for a query of max length (2 * group_size * numRegs)
    */
    template<
        class ScoreType, 
        int blocksize, 
        int group_size, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator
    >
    void call_GaplessFilter_strided_PSSM_multitile_uint8x4_kernel(
        int numThreadBlocks,
        const char * const devChars,
        ScoreOutputIterator const devAlignmentScores,
        const size_t* const devOffsets,
        const SequenceLengthT* const devLengths,
        PositionsIterator const d_positions_of_selected_lengths,
        const int numSelected,
        const SequenceLengthT queryLength,
        const PSSM_2D_View<ScoreType>& strided_PSSM,
        const int substitutionScoreBias, // values in strided_PSSM are computed from original PSSM + bias
        std::uint32_t* const multiTileTempStorage,
        size_t tempStorageElementsPerGroup, //number of std::uint32_t per group
        cudaStream_t stream
    ){
        if(substitutionScoreBias > 255){
            throw std::runtime_error("substitutionScoreBias cannot be represented with 8 bits.");
        }
        const ScoreType substitutionScoreBias_vec4(
            substitutionScoreBias,
            substitutionScoreBias,
            substitutionScoreBias,
            substitutionScoreBias
        );

        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM_INT8
        constexpr int numColumnsPSSM = std::max(group_size,8) * numRegs;
        #else
        constexpr int numColumnsPSSM = group_size * numRegs;
        #endif
        using SharedPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;

        int smem = sizeof(SharedPSSM);
        auto kernel = GaplessFilter_strided_PSSM_multitile_uint8x4_kernel<ScoreType, blocksize, group_size, numRegs, subjectIsCaseSensitive,
            ScoreOutputIterator>;

        auto setSmemKernelAttribute = [&](){
            static std::map<int, bool> isSet;
            if(smem > 48*1024){
                int deviceId;
                cudaGetDevice(&deviceId); CUERR;
                if(!isSet[deviceId]){
                    cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, smem); CUERR;
                    isSet[deviceId] = true;
                }
            }
        };
        setSmemKernelAttribute();

        dim3 grid = std::min(numSelected, numThreadBlocks);

        kernel<<<grid, blocksize, smem, stream>>>(
            devChars,
            devAlignmentScores,
            devOffsets,
            devLengths,
            d_positions_of_selected_lengths,
            numSelected,       
            queryLength,
            strided_PSSM,
            substitutionScoreBias_vec4,
            multiTileTempStorage,
            tempStorageElementsPerGroup
        ); CUERR;
    }



    template<class ScoreType, int blocksize, class ScoreOutputIterator>
    void call_GaplessFilter_strided_PSSM_multitile_uint8x4_kernel(
        int numThreadBlocks,
        int group_size,
        int numRegs,
        const char * const devChars,
        ScoreOutputIterator const devAlignmentScores,
        const size_t* const devOffsets,
        const SequenceLengthT* const devLengths,
        PositionsIterator const d_positions_of_selected_lengths,
        const int numSelected,
        const SequenceLengthT queryLength,
        const PSSM_2D_View<ScoreType>& strided_PSSM,
        const int substitutionScoreBias, // values in strided_PSSM are computed from original PSSM + bias
        std::uint32_t* const multiTileTempStorage,
        size_t tempStorageElementsPerGroup, //number of std::uint32_t per group
        cudaStream_t stream
    ){
        constexpr bool subjectIsCaseSensitive = true;

        #define X(g,r) \
            if(group_size == g && numRegs == r){ \
                call_GaplessFilter_strided_PSSM_multitile_uint8x4_kernel<ScoreType, blocksize,g,r,subjectIsCaseSensitive>( \
                    numThreadBlocks, devChars, devAlignmentScores, devOffsets, devLengths, d_positions_of_selected_lengths, \
                    numSelected, queryLength, strided_PSSM, substitutionScoreBias, multiTileTempStorage, tempStorageElementsPerGroup, stream \
                ); \
            } else 

            PSSM_GAPLESS_INT8_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X
        { throw std::runtime_error("invalid groupsize/numregs config");}

        #undef X
    }









    #define ScoreOutputIterator TopNMaximaArray
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void uint8x4::call_GaplessFilter_strided_PSSM_multitile_uint8x4_kernel<ScoreType_u8x4, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<ScoreType_u8x4>&, \
            const int, \
            std::uint32_t*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_INT8_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef ScoreOutputIterator




} 



} //namespace cudasw4

#endif