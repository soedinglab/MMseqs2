#ifndef PSSM_KERNELS_CUH
#define PSSM_KERNELS_CUH

#include "../../cuda_hip_rename.h"

#if defined(__CUDACC__)
#include <cuda_fp16.h>
#endif

#include <map>

#include "../../pssm.cuh"
#include "../../convert.cuh"
#include "../../mathops.cuh"
#include "../../util.cuh"
#include "gapless_common.cuh"
#include "../../offset_iterator.cuh"
#include "../alignment_interface_data.cuh"

#if defined(__CUDACC__)
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#endif
#if defined(__HIPCC__)
    #include "hip/hip_runtime.h"
    #include <hip/hip_cooperative_groups.h>
#endif
namespace cg = cooperative_groups;
#include "../../cuda_hip_compatibility.cuh"



#include "../../namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN



namespace hardcodedzero{


    template<
        class ScoreType, 
        int numRegs, 
        class Group, 
        class SharedPSSM, 
        class SmemIndexCalculator,
        int relaxChunkSize_,
        class UpdateMaxOp
    >
    struct GaplessPSSMState{
        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOps = MathOps<ScoreType>;

        static constexpr int relaxChunkSize = 4; // 16-byte loads load 4 half2/short2 quantities
        static_assert(relaxChunkSize == 4);
        static_assert(relaxChunkSize_ == 4);

        static constexpr int numRelaxChunks = SDIV(numRegs, relaxChunkSize);
        static constexpr int numRegsInLastRelaxChunk = numRegs - (numRelaxChunks-1)*relaxChunkSize;

        ScoreType penalty_here_array[numRegs];
        ScoreType penalty_diag{}; //0
        SharedPSSM& shared_strided_PSSM;
        Group& group;
        UpdateMaxOp updateMaxOp;

        __device__
        GaplessPSSMState(SharedPSSM& s, Group& g, UpdateMaxOp updateMaxOp_)
            : shared_strided_PSSM(s), group(g), updateMaxOp(updateMaxOp_){}

        __device__
        void resetScores(){
            #pragma unroll
            for(int i = 0; i < numRegs; i++){
                penalty_here_array[i] = ScoreType{};
            }
            
            penalty_diag = ScoreType{};
        }

        __device__
        void relax(int subject_letter, int row, int tileNr){
            SmemIndexCalculator smemIndexCalculator;

            ScoreType penalty_temp[2];

            const auto* pssmptr = &shared_strided_PSSM.data[subject_letter][0];
            float4 foo = *((float4*)&pssmptr[smemIndexCalculator.getIndexOfFirstElementInVectorType(0)]);

            ScoreType foo2[relaxChunkSize];
            memcpy(&foo2[0], &foo, sizeof(float4));

            // printf("thread %d, subject_letter %d, loaded (%d %d) (%d %d) (%d %d) (%d %d)\n", 
            //     group.thread_rank(), subject_letter ,
            //     int(foo2[0].x), int(foo2[0].y),
            //     int(foo2[1].x), int(foo2[1].y),
            //     int(foo2[2].x), int(foo2[2].y),
            //     int(foo2[3].x), int(foo2[3].y)
            // );

            penalty_temp[0] = penalty_here_array[0];
            penalty_here_array[0] = MathOps::add_relu(penalty_diag, foo2[0]);

            constexpr int loopLimit1 = (numRelaxChunks == 1) ? numRegsInLastRelaxChunk : relaxChunkSize;
            #pragma unroll
            for(int k = 1; k < loopLimit1; k++){
                penalty_temp[k % 2] = penalty_here_array[k];
                penalty_here_array[k] = MathOps::add_relu(penalty_temp[(k+1)%2], foo2[k]);
            }

            #pragma unroll
            for(int i = 1; i < numRelaxChunks; i++){
                foo = *((float4*)&pssmptr[smemIndexCalculator.getIndexOfFirstElementInVectorType(i)]);
                memcpy(&foo2[0], &foo, sizeof(float4));

                const int loopLimit2 = (i == numRelaxChunks-1) ? numRegsInLastRelaxChunk : relaxChunkSize;
                #pragma unroll
                for(int k = 0; k < relaxChunkSize; k++){
                    if(k < loopLimit2){
                        penalty_temp[k % 2] = penalty_here_array[relaxChunkSize*i + k];
                        penalty_here_array[relaxChunkSize*i + k] = MathOps::add_relu(penalty_temp[(k+1)%2], foo2[k]);
                    }
                }
            }

            // printf("thread %d, subject_letter %d, computed (%d %d) (%d %d) (%d %d) (%d %d)\n", 
            //     group.thread_rank(), subject_letter ,
            //     int(penalty_here_array[0].x), int(penalty_here_array[0].y),
            //     int(penalty_here_array[1].x), int(penalty_here_array[1].y),
            //     int(penalty_here_array[2].x), int(penalty_here_array[2].y),
            //     int(penalty_here_array[3].x), int(penalty_here_array[3].y)
            // );

            updateMaxOp(penalty_here_array, tileNr, row);
        }

        __device__
        void shuffleScores(const Scalar& border_in){
            penalty_diag = group.shfl_up(penalty_here_array[numRegs-1], 1);
            const ScoreType penalty_temp0 = compatibility::group_shfl_down(group, penalty_here_array[numRegs-1], group.size()-1);

            if (group.thread_rank() == 0) {
                penalty_diag.x = border_in;
                penalty_diag.y = penalty_temp0.x;
            }
        }

        __device__
        void stepSingleTile(int subject_letter, int row, bool print = false){
            if(print && group.thread_rank() == 0){
                printf("pos %d, letter %d\n", -1, subject_letter);
            }
            relax(subject_letter, row, 0);
            shuffleScores(Scalar{});
        }

        __device__
        void stepFirstTile(int subject_letter, int row, Scalar& border_out){
            relax(subject_letter, row, 0);
            shuffleScores(Scalar{});
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1].y;
            }
        }

        __device__
        void stepIntermediateTile(int subject_letter, int row, int tileNr, const Scalar& border_in, Scalar& border_out){
            relax(subject_letter, row, tileNr);
            shuffleScores(border_in);
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1].y;
            }
        }

        __device__
        void stepLastTile(int subject_letter, int row, int tileNr, const Scalar& border_in){
            relax(subject_letter, row, tileNr);
            shuffleScores(border_in);
        }
    };


    /*
    PSSM kernel for a query of max length (2 * groupsize * numRegs)
    */
    template<
        class ScoreType, 
        int blocksize, 
        int groupsize, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PssmScoreType,
        class InputData
    >
    __global__
    __launch_bounds__(blocksize,1)
    void GaplessFilter_strided_PSSM_singletile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(std::is_same_v<PssmScoreType, ScoreType>);

        static_assert(blocksize % groupsize == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % groupsize == 0);

        constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
        constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();

        constexpr int numRowsPSSM = 21;
        constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
        constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;

        using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;

        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOpsScalar = MathOps<Scalar>;
        using MathOps = MathOps<ScoreType>;

        extern  __shared__ char externalSmem[];

        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);


        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        //const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;
        using SmemIndexCalculator = PssmSmemIndexCalculator<groupsize, pssmReplicationFactor, relaxChunkSize>;

        using UpdateMaxOp = UpdateMax_gapless_scoreOnly<ScoreType>;
        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize, UpdateMaxOp>;


        auto load_PSSM_single = [&]() {
            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][col];
            }
            __syncthreads();
        };

        auto load_PSSM_replicated = [&]() {
            SmemIndexCalculator pssmSmemIndexCalculator;
            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                auto value = strided_PSSM[letter][col];

                auto repIndices = pssmSmemIndexCalculator.computeReplicatedColumnIndices(col);

                for(auto repIndex : repIndices){
                    shared_strided_PSSM.data[letter][repIndex] = value;
                }
            }
            __syncthreads();
        };

        auto load_PSSM = [&](){
            if constexpr(SmemIndexCalculator::factor > 1){
                load_PSSM_replicated();
            }else{
                load_PSSM_single();
            }

            // __syncthreads();
            // if(threadIdx.x == 0){
            //     printf("smem pssm\n");
            //     for(int r = 0; r < 21; r++){
            //         for(int c = 0; c < numColumnsPSSM; c++){
            //             auto s = shared_strided_PSSM.data[r][c];
            //             half2 hs = static_cast<half2>(s);
            //             printf("%d %d ", int(hs.x), int(hs.y));
            //         }
            //         printf("\n");
            //     }
            // }
            // __syncthreads();
        };



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

        const SequenceLengthT queryLength = inputData.getQueryLength();
        const int numAlignments = inputData.getNumAlignments();

        const int alignmentId = idOfGroupInGrid;
        if(alignmentId < numAlignments){
            char4 new_subject_letter4;
            const int8_t* subjectData = inputData.getSubject(alignmentId);
            const SequenceLengthT subjectLength = inputData.getSubjectLength(alignmentId);
            const char4* subjectAsChar4 = reinterpret_cast<const char4*>(subjectData);
           

            State state(shared_strided_PSSM, group, UpdateMaxOp{});
            state.resetScores();

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

            auto maximumTracker = state.updateMaxOp;
            int myMaximum = MathOpsScalar::max(maximumTracker.maximum.x, maximumTracker.maximum.y);
               
            const int alignmentScore = compatibility::group_reduce_max(group, myMaximum);

            if(group.thread_rank() == 0){
                devAlignmentScores[alignmentId] = alignmentScore;
            }
        }
    }





    /*
    PSSM kernel for arbitrary query length
    */
    template<
        class ScoreType, 
        int blocksize, 
        int groupsize, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PssmScoreType,
        class InputData
    >
    __global__
    __launch_bounds__(blocksize,1)
    void GaplessFilter_strided_PSSM_multitile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM,
        __grid_constant__ float2* const multiTileTempStorage,
        __grid_constant__ const size_t tempStorageElementsPerGroup
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(std::is_same_v<PssmScoreType, ScoreType>);

        static_assert(blocksize % groupsize == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % groupsize == 0);

        constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
        constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();
        
        constexpr int numRowsPSSM = 21;
        constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
        constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;
        using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
        
        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOpsScalar = MathOps<Scalar>;
        using MathOps = MathOps<ScoreType>;
        
        extern  __shared__ char externalSmem[];
        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);

        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int numGroupsInBlock = blockDim.x / groupsize;
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;
        
        const size_t groupTempStorageOffset = idOfGroupInGrid * tempStorageElementsPerGroup;
        float2* const groupTempStorage = multiTileTempStorage + groupTempStorageOffset;
        
        const SequenceLengthT queryLength = inputData.getQueryLength();
        const int numTiles = SDIV(queryLength, 2 * groupsize * numRegs);

        using SmemIndexCalculator = PssmSmemIndexCalculator<groupsize, pssmReplicationFactor, relaxChunkSize>;

        using UpdateMaxOp = UpdateMax_gapless_scoreOnly<ScoreType>;
        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize, UpdateMaxOp>;

        alignas(8) Scalar penalty_in[4];
        alignas(8) Scalar penalty_out[4];

        auto load_PSSM_single = [&](int tileNr) {
            const int columnOffset = tileNr * pssmTileSize;
            __syncthreads(); //wait for all groups before overwriting pssm

            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                int letter = i/(pssmTileSize);
                int col = i%(pssmTileSize);
                //shared_strided_PSSM.data[letter][col] = strided_PSSM_1d[i];
                //shared_strided_PSSM.data[letter][col] = strided_PSSM.data[i];
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][columnOffset + col];
            }
            __syncthreads();
        };

        auto load_PSSM_replicated = [&](int tileNr) {
            const int columnOffset = tileNr * pssmTileSize;
            __syncthreads(); //wait for all groups before overwriting pssm

            SmemIndexCalculator pssmSmemIndexCalculator;

            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                auto value = strided_PSSM[letter][columnOffset + col];

                auto repIndices = pssmSmemIndexCalculator.computeReplicatedColumnIndices(col);

                for(auto repIndex : repIndices){
                    shared_strided_PSSM.data[letter][repIndex] = value;
                }
            }
            __syncthreads();
        };

        auto load_PSSM = [&](int tileNr){
            if constexpr(SmemIndexCalculator::factor > 1){
                load_PSSM_replicated(tileNr);
            }else{
                load_PSSM_single(tileNr);
            }
        };



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

        const int numAlignments = inputData.getNumAlignments();
        //need to round up to blocks because loading pssm is a block-wide operation
        const int numAlignmentsRoundedUp = SDIV(numAlignments, numGroupsInBlock) * numGroupsInBlock;

        for(int alignmentId = idOfGroupInGrid; alignmentId < numAlignmentsRoundedUp; alignmentId += numGroupsInGrid){
            char4 new_subject_letter4;
            SequenceLengthT subjectLength;
            const char4* subjectAsChar4;

            State state(shared_strided_PSSM, group, UpdateMaxOp{});

            //first tile
            {
                /* 
                    -----------------------
                    Process tile 0
                    ----------------------- 
                */

                //load pssm for tile 0. blockwide operation
                load_PSSM(0);

                if(alignmentId < numAlignments){
                    const int8_t* subjectData = inputData.getSubject(alignmentId);
                    subjectLength = inputData.getSubjectLength(alignmentId);
                    subjectAsChar4 = reinterpret_cast<const char4*>(subjectData);
                    

                    state.resetScores();

                    int k;

                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {

                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);

                        state.stepFirstTile(new_subject_letter4.x, k, penalty_out[0]);
                        state.stepFirstTile(new_subject_letter4.y, k+1, penalty_out[1]);
                        state.stepFirstTile(new_subject_letter4.z, k+2, penalty_out[2]);
                        state.stepFirstTile(new_subject_letter4.w, k+3, penalty_out[3]);
                        
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
                        }
                    }

                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        state.stepFirstTile(new_subject_letter4.x, k, penalty_out[0]);
                    }

                    if (subjectLength%4 >= 2) {
                        state.stepFirstTile(new_subject_letter4.y, k+1, penalty_out[1]);
                    }

                    if (subjectLength%4 >= 3) {
                        state.stepFirstTile(new_subject_letter4.z, k+2, penalty_out[2]);
                    }

                    //if there were remaining rows, update temp storage
                    if(subjectLength % 4 > 0){
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
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

                if(alignmentId < numAlignments){    
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepIntermediateTile(new_subject_letter4.x, k, tileNr, penalty_in[0], penalty_out[0]);
                        state.stepIntermediateTile(new_subject_letter4.y, k+1, tileNr, penalty_in[1], penalty_out[1]);
                        state.stepIntermediateTile(new_subject_letter4.z, k+2, tileNr, penalty_in[2], penalty_out[2]);
                        state.stepIntermediateTile(new_subject_letter4.w, k+3, tileNr, penalty_in[3], penalty_out[3]);
            
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
                        }
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }
                        state.stepIntermediateTile(new_subject_letter4.x, k, tileNr, penalty_in[0], penalty_out[0]);
                    }
        
                    if (subjectLength%4 >= 2) {
                        state.stepIntermediateTile(new_subject_letter4.y, k+1, tileNr, penalty_in[1], penalty_out[1]);
                    }
        
                    if (subjectLength%4 >= 3) {
                        state.stepIntermediateTile(new_subject_letter4.z, k+2, tileNr, penalty_in[2], penalty_out[2]);
                    }
        
                    //if there were remaining rows, update temp storage
                    if(subjectLength % 4 > 0){
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
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
                const int tileNr = numTiles-1;

                //load pssm for tile (numTiles-1). blockwide operation
                load_PSSM(numTiles-1);

                if(alignmentId < numAlignments){
        
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepLastTile(new_subject_letter4.x, k, tileNr, penalty_in[0]);
                        state.stepLastTile(new_subject_letter4.y, k+1, tileNr, penalty_in[1]);
                        state.stepLastTile(new_subject_letter4.z, k+2, tileNr, penalty_in[2]);
                        state.stepLastTile(new_subject_letter4.w, k+3, tileNr, penalty_in[3]);
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }
                        state.stepLastTile(new_subject_letter4.x, k, tileNr, penalty_in[0]);
                    }
        
                    if (subjectLength%4 >= 2) {
                        state.stepLastTile(new_subject_letter4.y, k+1, tileNr, penalty_in[1]);
                    }
        
                    if (subjectLength%4 >= 3) {
                        state.stepLastTile(new_subject_letter4.z, k+2, tileNr, penalty_in[2]);
                    }
                }
            }

            if(alignmentId < numAlignments){
                auto maximumTracker = state.updateMaxOp;
                int myMaximum = MathOpsScalar::max(maximumTracker.maximum.x, maximumTracker.maximum.y);
                
                const int alignmentScore = compatibility::group_reduce_max(group, myMaximum);

                if(group.thread_rank() == 0){
                    devAlignmentScores[alignmentId] = alignmentScore;
                }
            }
        }

    }










    /*
    PSSM kernel for a query of max length (2 * groupsize * numRegs)
    */
    template<
        class ScoreType, 
        int blocksize, 
        int groupsize, 
        int numRegs, 
        bool withEndPosition,
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PssmScoreType,
        class InputData
    >
    __global__
    __launch_bounds__(blocksize,1)
    void GaplessFilter_strided_PSSM_endPos_singletile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ int* const subjectEndPositions_exclusive, //not accessed if withEndPosition == false
        __grid_constant__ int* const queryEndPositions_exclusive, //not accessed if withEndPosition == false
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(std::is_same_v<PssmScoreType, ScoreType>);

        static_assert(blocksize % groupsize == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % groupsize == 0);

        constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
        constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();

        constexpr int numRowsPSSM = 21;
        constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
        constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;

        using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;

        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOpsScalar = MathOps<Scalar>;
        using MathOps = MathOps<ScoreType>;

        extern  __shared__ char externalSmem[];

        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);


        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        //const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;
        using SmemIndexCalculator = PssmSmemIndexCalculator<groupsize, pssmReplicationFactor, relaxChunkSize>;

        using UpdateMaxOp = typename std::conditional<
            withEndPosition,
            UpdateMax_gapless_vec2_endPos_singleTile<ScoreType>,
            UpdateMax_gapless_scoreOnly<ScoreType>
        >::type;
        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize, UpdateMaxOp>;


        auto load_PSSM_single = [&]() {
            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][col];
            }
            __syncthreads();
        };

        auto load_PSSM_replicated = [&]() {
            SmemIndexCalculator pssmSmemIndexCalculator;
            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                auto value = strided_PSSM[letter][col];

                auto repIndices = pssmSmemIndexCalculator.computeReplicatedColumnIndices(col);

                for(auto repIndex : repIndices){
                    shared_strided_PSSM.data[letter][repIndex] = value;
                }
            }
            __syncthreads();
        };

        auto load_PSSM = [&](){
            if constexpr(SmemIndexCalculator::factor > 1){
                load_PSSM_replicated();
            }else{
                load_PSSM_single();
            }

            // __syncthreads();
            // if(threadIdx.x == 0){
            //     printf("smem pssm\n");
            //     for(int r = 0; r < 21; r++){
            //         printf("r %d\n", r);
            //         for(int c = 0; c < numColumnsPSSM; c++){
            //             auto s = shared_strided_PSSM.data[r][c];
            //             auto hs = static_cast<ScoreType>(s);
            //             printf("%d %d ", int(hs.x), int(hs.y));
            //         }
            //         printf("\n");
            //     }
            // }
            // __syncthreads();
        };



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

        const SequenceLengthT queryLength = inputData.getQueryLength();
        const int numAlignments = inputData.getNumAlignments();

        const int alignmentId = idOfGroupInGrid;
        if(alignmentId < numAlignments){
            char4 new_subject_letter4;
            const int8_t* subjectData = inputData.getSubject(alignmentId);
            const SequenceLengthT subjectLength = inputData.getSubjectLength(alignmentId);
            const char4* subjectAsChar4 = reinterpret_cast<const char4*>(subjectData);
            

            State state(shared_strided_PSSM, group, UpdateMaxOp{});
            state.resetScores();

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

            //find thread-local optimum
            auto maximumTracker = state.updateMaxOp;
            if constexpr(withEndPosition){
                int myMaximum = -1;
                int myQueryEndInclusive = 0;
                int mySubjectEndInclusive = 0;

                // printf("begin. thread %d, maximum %d %d, regAndLane %d %d, row %d %d\n", 
                //     group.thread_rank(), int(maximumTracker.maximum.x), int(maximumTracker.maximum.y),
                //    int(maximumTracker.positionOfMaxObserved_regAndLane.x), int(maximumTracker.positionOfMaxObserved_regAndLane.y),
                //    int(maximumTracker.positionOfMaxObserved_row.x), int(maximumTracker.positionOfMaxObserved_row.y));

                

                int reg = maximumTracker.positionOfMaxObserved_regAndLane.x / 2;
                int lane = maximumTracker.positionOfMaxObserved_regAndLane.x % 2;
                int regInGroup = group.thread_rank() * numRegs + reg;
                int queryPos = regInGroup + lane * (group.num_threads() * numRegs);
                int subjectPos = maximumTracker.positionOfMaxObserved_row.x;
                int maximum = maximumTracker.maximum.x;
                
                if(queryPos < queryLength && subjectPos < subjectLength){
                    if(maximum > myMaximum){
                        myMaximum = maximum;
                        myQueryEndInclusive = queryPos;
                        mySubjectEndInclusive = subjectPos;
                    }
                }

                reg = maximumTracker.positionOfMaxObserved_regAndLane.y / 2;
                lane = maximumTracker.positionOfMaxObserved_regAndLane.y % 2;
                regInGroup = group.thread_rank() * numRegs + reg;
                queryPos = regInGroup + lane * (group.num_threads() * numRegs);
                subjectPos = maximumTracker.positionOfMaxObserved_row.y;
                maximum = maximumTracker.maximum.y;
                if(queryPos < queryLength && subjectPos < subjectLength){
                    if(maximum > myMaximum){
                        myMaximum = maximum;
                        myQueryEndInclusive = queryPos;
                        mySubjectEndInclusive = subjectPos;
                    }
                }

                //find group-wide optimum
                const int3 packed = make_int3(myMaximum, 
                    myQueryEndInclusive,
                    mySubjectEndInclusive);
                const int3 maxPacked = compatibility::group_reduce(group, packed, [](int3 l, int3 r){
                    if(l.x > r.x){
                        return l;
                    }else{
                        return r;
                    }
                });

                const int alignmentScore = maxPacked.x;
                const int queryEndExclusive = maxPacked.y+1;
                const int subjectEndExclusive = maxPacked.z+1;

                if(group.thread_rank() == 0){
                    devAlignmentScores[alignmentId] = alignmentScore;
                    subjectEndPositions_exclusive[alignmentId] = subjectEndExclusive;
                    queryEndPositions_exclusive[alignmentId] = queryEndExclusive;
                }            
            }else{
                int myMaximum = MathOpsScalar::max(maximumTracker.maximum.x, maximumTracker.maximum.y);
                
                const int alignmentScore = compatibility::group_reduce_max(group, myMaximum);

                if(group.thread_rank() == 0){
                    devAlignmentScores[alignmentId] = alignmentScore;
                }
            }
        }
    }





    /*
    PSSM kernel for arbitrary query length
    */
    template<
        class ScoreType, 
        int blocksize, 
        int groupsize, 
        int numRegs, 
        bool withEndPosition,
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PssmScoreType,
        class InputData
    >
    __global__
    __launch_bounds__(blocksize,1)
    void GaplessFilter_strided_PSSM_endPos_multitile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ int* const subjectEndPositions_exclusive, //not accessed if withEndPosition == false
        __grid_constant__ int* const queryEndPositions_exclusive, //not accessed if withEndPosition == false
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM,
        __grid_constant__ float2* const multiTileTempStorage,
        __grid_constant__ const size_t tempStorageElementsPerGroup
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(std::is_same_v<PssmScoreType, ScoreType>);

        static_assert(blocksize % groupsize == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % groupsize == 0);

        constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
        constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();
        
        constexpr int numRowsPSSM = 21;
        constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
        constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;
        using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
        
        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOpsScalar = MathOps<Scalar>;
        using MathOps = MathOps<ScoreType>;
        
        extern  __shared__ char externalSmem[];
        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);

        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int numGroupsInBlock = blockDim.x / groupsize;
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;
        
        const size_t groupTempStorageOffset = idOfGroupInGrid * tempStorageElementsPerGroup;
        float2* const groupTempStorage = multiTileTempStorage + groupTempStorageOffset;
        
        const SequenceLengthT queryLength = inputData.getQueryLength();
        const int numTiles = SDIV(queryLength, 2 * groupsize * numRegs);

        using SmemIndexCalculator = PssmSmemIndexCalculator<groupsize, pssmReplicationFactor, relaxChunkSize>;

        using UpdateMaxOp = typename std::conditional<
            withEndPosition,
            UpdateMax_gapless_vec2_endPos_multiTile<ScoreType>,
            UpdateMax_gapless_scoreOnly<ScoreType>
        >::type;
        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize, UpdateMaxOp>;

        alignas(8) Scalar penalty_in[4];
        alignas(8) Scalar penalty_out[4];

        auto load_PSSM_single = [&](int tileNr) {
            const int columnOffset = tileNr * pssmTileSize;
            __syncthreads(); //wait for all groups before overwriting pssm

            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                int letter = i/(pssmTileSize);
                int col = i%(pssmTileSize);
                //shared_strided_PSSM.data[letter][col] = strided_PSSM_1d[i];
                //shared_strided_PSSM.data[letter][col] = strided_PSSM.data[i];
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][columnOffset + col];
            }
            __syncthreads();
        };

        auto load_PSSM_replicated = [&](int tileNr) {
            const int columnOffset = tileNr * pssmTileSize;
            __syncthreads(); //wait for all groups before overwriting pssm

            SmemIndexCalculator pssmSmemIndexCalculator;

            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                auto value = strided_PSSM[letter][columnOffset + col];

                auto repIndices = pssmSmemIndexCalculator.computeReplicatedColumnIndices(col);

                for(auto repIndex : repIndices){
                    shared_strided_PSSM.data[letter][repIndex] = value;
                }
            }
            __syncthreads();
        };

        auto load_PSSM = [&](int tileNr){
            if constexpr(SmemIndexCalculator::factor > 1){
                load_PSSM_replicated(tileNr);
            }else{
                load_PSSM_single(tileNr);
            }
        };



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


        const int numAlignments = inputData.getNumAlignments();
        //need to round up to blocks because loading pssm is a block-wide operation
        const int numAlignmentsRoundedUp = SDIV(numAlignments, numGroupsInBlock) * numGroupsInBlock;

        for(int alignmentId = idOfGroupInGrid; alignmentId < numAlignmentsRoundedUp; alignmentId += numGroupsInGrid){
            char4 new_subject_letter4;
            SequenceLengthT subjectLength;
            const char4* subjectAsChar4;

            State state(shared_strided_PSSM, group, UpdateMaxOp{});

            //first tile
            {
                /* 
                    -----------------------
                    Process tile 0
                    ----------------------- 
                */

                //load pssm for tile 0. blockwide operation
                load_PSSM(0);

                if(alignmentId < numAlignments){
                    const int8_t* subjectData = inputData.getSubject(alignmentId);
                    subjectLength = inputData.getSubjectLength(alignmentId);
                    subjectAsChar4 = reinterpret_cast<const char4*>(subjectData);
                    

                    state.resetScores();

                    int k;

                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {

                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);

                        state.stepFirstTile(new_subject_letter4.x, k, penalty_out[0]);
                        state.stepFirstTile(new_subject_letter4.y, k+1, penalty_out[1]);
                        state.stepFirstTile(new_subject_letter4.z, k+2, penalty_out[2]);
                        state.stepFirstTile(new_subject_letter4.w, k+3, penalty_out[3]);
                        
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
                        }
                    }

                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        state.stepFirstTile(new_subject_letter4.x, k, penalty_out[0]);
                    }

                    if (subjectLength%4 >= 2) {
                        state.stepFirstTile(new_subject_letter4.y, k+1, penalty_out[1]);
                    }

                    if (subjectLength%4 >= 3) {
                        state.stepFirstTile(new_subject_letter4.z, k+2, penalty_out[2]);
                    }

                    //if there were remaining rows, update temp storage
                    if(subjectLength % 4 > 0){
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
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

                if(alignmentId < numAlignments){    
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepIntermediateTile(new_subject_letter4.x, k, tileNr, penalty_in[0], penalty_out[0]);
                        state.stepIntermediateTile(new_subject_letter4.y, k+1, tileNr, penalty_in[1], penalty_out[1]);
                        state.stepIntermediateTile(new_subject_letter4.z, k+2, tileNr, penalty_in[2], penalty_out[2]);
                        state.stepIntermediateTile(new_subject_letter4.w, k+3, tileNr, penalty_in[3], penalty_out[3]);
            
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
                        }
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }
                        state.stepIntermediateTile(new_subject_letter4.x, k, tileNr, penalty_in[0], penalty_out[0]);
                    }
        
                    if (subjectLength%4 >= 2) {
                        state.stepIntermediateTile(new_subject_letter4.y, k+1, tileNr, penalty_in[1], penalty_out[1]);
                    }
        
                    if (subjectLength%4 >= 3) {
                        state.stepIntermediateTile(new_subject_letter4.z, k+2, tileNr, penalty_in[2], penalty_out[2]);
                    }
        
                    //if there were remaining rows, update temp storage
                    if(subjectLength % 4 > 0){
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
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
                const int tileNr = numTiles-1;

                //load pssm for tile (numTiles-1). blockwide operation
                load_PSSM(numTiles-1);

                if(alignmentId < numAlignments){
        
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepLastTile(new_subject_letter4.x, k, tileNr, penalty_in[0]);
                        state.stepLastTile(new_subject_letter4.y, k+1, tileNr, penalty_in[1]);
                        state.stepLastTile(new_subject_letter4.z, k+2, tileNr, penalty_in[2]);
                        state.stepLastTile(new_subject_letter4.w, k+3, tileNr, penalty_in[3]);
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }
                        state.stepLastTile(new_subject_letter4.x, k, tileNr, penalty_in[0]);
                    }
        
                    if (subjectLength%4 >= 2) {
                        state.stepLastTile(new_subject_letter4.y, k+1, tileNr, penalty_in[1]);
                    }
        
                    if (subjectLength%4 >= 3) {
                        state.stepLastTile(new_subject_letter4.z, k+2, tileNr, penalty_in[2]);
                    }
                }
            }

            if(alignmentId < numAlignments){
                //find thread-local optimum
                auto maximumTracker = state.updateMaxOp;
                if constexpr(withEndPosition){
                    int myMaximum = -1;
                    int myQueryEndInclusive = 0;
                    int mySubjectEndInclusive = 0;


                    int reg = maximumTracker.positionOfMaxObserved_regAndLane.x / 2;
                    int lane = maximumTracker.positionOfMaxObserved_regAndLane.x % 2;
                    int regInGroup = group.thread_rank() * numRegs + reg;
                    int queryPos = maximumTracker.positionOfMaxObserved_tileNr.x * 2 * (group.num_threads() * numRegs)
                        + regInGroup + lane * (group.num_threads() * numRegs);
                    int subjectPos = maximumTracker.positionOfMaxObserved_row.x;
                    int maximum = maximumTracker.maximum.x;
                    if(queryPos < queryLength && subjectPos < subjectLength){
                        if(maximum > myMaximum){
                            myMaximum = maximum;
                            myQueryEndInclusive = queryPos;
                            mySubjectEndInclusive = subjectPos;
                        }
                    }

                    reg = maximumTracker.positionOfMaxObserved_regAndLane.y / 2;
                    lane = maximumTracker.positionOfMaxObserved_regAndLane.y % 2;
                    regInGroup = group.thread_rank() * numRegs + reg;
                    queryPos = maximumTracker.positionOfMaxObserved_tileNr.y * 2 * (group.num_threads() * numRegs)
                        + regInGroup + lane * (group.num_threads() * numRegs);
                    subjectPos = maximumTracker.positionOfMaxObserved_row.y;
                    maximum = maximumTracker.maximum.y;
                    if(queryPos < queryLength && subjectPos < subjectLength){
                        if(maximum > myMaximum){
                            myMaximum = maximum;
                            myQueryEndInclusive = queryPos;
                            mySubjectEndInclusive = subjectPos;
                        }
                    }

                    //find group-wide optimum
                    const int3 packed = make_int3(myMaximum, 
                        myQueryEndInclusive,
                        mySubjectEndInclusive);
                    const int3 maxPacked = compatibility::group_reduce(group, packed, [](int3 l, int3 r){
                        if(l.x > r.x){
                            return l;
                        }else{
                            return r;
                        }
                    });

                    const int alignmentScore = maxPacked.x;
                    const int queryEndExclusive = maxPacked.y+1;
                    const int subjectEndExclusive = maxPacked.z+1;

                    if(group.thread_rank() == 0){
                        devAlignmentScores[alignmentId] = alignmentScore;
                        subjectEndPositions_exclusive[alignmentId] = subjectEndExclusive;
                        queryEndPositions_exclusive[alignmentId] = queryEndExclusive;
                    }            
                }else{
                    int myMaximum = MathOpsScalar::max(maximumTracker.maximum.x, maximumTracker.maximum.y);
                    
                    const int alignmentScore = compatibility::group_reduce_max(group, myMaximum);

                    if(group.thread_rank() == 0){
                        devAlignmentScores[alignmentId] = alignmentScore;
                    }
                }
            }
        }

    }


} //namespace hardcodedzero


namespace kernelparamzero{    


    template<
        class ScoreType, 
        int numRegs, 
        class Group, 
        class SharedPSSM, 
        class SmemIndexCalculator,
        int relaxChunkSize_
    >
    struct GaplessPSSMState{
        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOps = MathOps<ScoreType>;

        static constexpr int relaxChunkSize = 4; // 16-byte loads load 4 half2/short2 quantities
        static_assert(relaxChunkSize == 4);
        static_assert(relaxChunkSize_ == 4);

        static constexpr int numRelaxChunks = SDIV(numRegs, relaxChunkSize);
        static constexpr int numRegsInLastRelaxChunk = numRegs - (numRelaxChunks-1)*relaxChunkSize;

        ScoreType penalty_here_array[numRegs];
        ScoreType maximum{}; //0
        ScoreType penalty_diag{}; //0
        SharedPSSM& shared_strided_PSSM;
        Group& group;

        __device__
        GaplessPSSMState(SharedPSSM& s, Group& g) : shared_strided_PSSM(s), group(g) {}

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

        __device__
        void relax(int subject_letter, ScoreType zero){
            SmemIndexCalculator smemIndexCalculator;

            ScoreType penalty_temp[2];

            const auto* row = &shared_strided_PSSM.data[subject_letter][0];
            float4 foo = *((float4*)&row[smemIndexCalculator.getIndexOfFirstElementInVectorType(0)]);

            ScoreType foo2[relaxChunkSize];
            memcpy(&foo2[0], &foo, sizeof(float4));

            penalty_temp[0] = penalty_here_array[0];
            penalty_here_array[0] = MathOps::add_relu(penalty_diag, foo2[0], zero);

            constexpr int loopLimit1 = (numRelaxChunks == 1) ? numRegsInLastRelaxChunk : relaxChunkSize;
            #pragma unroll
            for(int k = 1; k < loopLimit1; k++){
                penalty_temp[k % 2] = penalty_here_array[k];
                penalty_here_array[k] = MathOps::add_relu(penalty_temp[(k+1)%2], foo2[k], zero);
            }

            #pragma unroll
            for(int i = 1; i < numRelaxChunks; i++){
                foo = *((float4*)&row[smemIndexCalculator.getIndexOfFirstElementInVectorType(i)]);
                memcpy(&foo2[0], &foo, sizeof(float4));

                const int loopLimit2 = (i == numRelaxChunks-1) ? numRegsInLastRelaxChunk : relaxChunkSize;
                #pragma unroll
                for(int k = 0; k < relaxChunkSize; k++){
                    if(k < loopLimit2){
                        penalty_temp[k % 2] = penalty_here_array[relaxChunkSize*i + k];
                        penalty_here_array[relaxChunkSize*i + k] = MathOps::add_relu(penalty_temp[(k+1)%2], foo2[k], zero);
                    }
                }
            }

            #pragma unroll
            for(int k = 0; k < numRegs/2; k++){
                maximum = MathOps::max3(maximum, penalty_here_array[2*k+0], penalty_here_array[2*k+1]);
            }
            if constexpr(numRegs % 2 == 1){
                maximum = MathOps::max(maximum, penalty_here_array[numRegs-1]);
            }
        }

        __device__
        void shuffleScores(const Scalar& border_in){
            penalty_diag = group.shfl_up(penalty_here_array[numRegs-1], 1);
            const ScoreType penalty_temp0 = compatibility::group_shfl_down(group, penalty_here_array[numRegs-1], group.size()-1);

            if (group.thread_rank() == 0) {
                penalty_diag.x = border_in;
                penalty_diag.y = penalty_temp0.x;
            }
        }

        __device__
        void stepSingleTile(int subject_letter, ScoreType zero){
            relax(subject_letter, zero);
            shuffleScores(Scalar{});
        }

        __device__
        void stepFirstTile(int subject_letter, Scalar& border_out, ScoreType zero){
            relax(subject_letter, zero);
            shuffleScores(Scalar{});
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1].y;
            }
        }

        __device__
        void stepIntermediateTile(int subject_letter, const Scalar& border_in, Scalar& border_out, ScoreType zero){
            relax(subject_letter, zero);
            shuffleScores(border_in);
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1].y;
            }
        }

        __device__
        void stepLastTile(int subject_letter, const Scalar& border_in, ScoreType zero){
            relax(subject_letter, zero);
            shuffleScores(border_in);
        }

        __device__
        void reduceMaximumScore(){
            maximum = MathOps::reduce_max(group, maximum);
        }
    };


    /*
    PSSM kernel for a query of max length (2 * groupsize * numRegs)
    */
    template<
        class ScoreType, 
        int blocksize, 
        int groupsize, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PssmScoreType,
        class InputData
    >
    __global__
    __launch_bounds__(blocksize,1)
    void GaplessFilter_strided_PSSM_singletile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM,
        __grid_constant__ const ScoreType zero
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(std::is_same_v<PssmScoreType, ScoreType>);
        static_assert(blocksize % groupsize == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % groupsize == 0);

        constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
        constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();

        constexpr int numRowsPSSM = 21;
        constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
        constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;

        using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
        using MathOps = MathOps<ScoreType>;

        extern  __shared__ char externalSmem[];

        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);


        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        //const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;

        using SmemIndexCalculator = PssmSmemIndexCalculator<groupsize, pssmReplicationFactor, relaxChunkSize>;

        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize>;
        
        auto load_PSSM_single = [&]() {
            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][col];
            }
            __syncthreads();
        };

        auto load_PSSM_replicated = [&]() {
            SmemIndexCalculator pssmSmemIndexCalculator;
            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                auto value = strided_PSSM[letter][col];

                auto repIndices = pssmSmemIndexCalculator.computeReplicatedColumnIndices(col);

                for(auto repIndex : repIndices){
                    shared_strided_PSSM.data[letter][repIndex] = value;
                }
            }
            __syncthreads();
        };

        auto load_PSSM = [&](){
            if constexpr(SmemIndexCalculator::factor > 1){
                load_PSSM_replicated();
            }else{
                load_PSSM_single();
            }
        };



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

        const SequenceLengthT queryLength = inputData.getQueryLength();
        const int numAlignments = inputData.getNumAlignments();

        const int alignmentId = idOfGroupInGrid;
        if(alignmentId < numAlignments){
            char4 new_subject_letter4;
            const int8_t* subjectData = inputData.getSubject(alignmentId);
            const SequenceLengthT subjectLength = inputData.getSubjectLength(alignmentId);
            const char4* subjectAsChar4 = reinterpret_cast<const char4*>(subjectData);
            

            State state(shared_strided_PSSM, group);
            state.resetScores();
            state.resetMaximum();

            int k;
            for (k=0; k<subjectLength-3; k+=4) {
                new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                state.stepSingleTile(new_subject_letter4.x, zero);
                state.stepSingleTile(new_subject_letter4.y, zero);
                state.stepSingleTile(new_subject_letter4.z, zero);
                state.stepSingleTile(new_subject_letter4.w, zero);
            }

            if (subjectLength%4 >= 1) {
                new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                state.stepSingleTile(new_subject_letter4.x, zero);
            }

            if (subjectLength%4 >= 2) {
                state.stepSingleTile(new_subject_letter4.y, zero);
            }

            if (subjectLength%4 >= 3) {
                state.stepSingleTile(new_subject_letter4.z, zero);
            }

            state.reduceMaximumScore();
            const float overall_max = MathOps::max(state.maximum.x, state.maximum.y);

            if(group.thread_rank() == 0){
                devAlignmentScores[alignmentId] = overall_max;
            }
        }
    }




    /*
    PSSM kernel for arbitrary query length
    */
    template<
        class ScoreType, 
        int blocksize, 
        int groupsize, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PssmScoreType,
        class InputData
    >
    __global__
    __launch_bounds__(blocksize,1)
    void GaplessFilter_strided_PSSM_multitile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM,
        __grid_constant__ float2* const multiTileTempStorage,
        __grid_constant__ const size_t tempStorageElementsPerGroup,
        __grid_constant__ const ScoreType zero
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(std::is_same_v<PssmScoreType, ScoreType>);
        static_assert(blocksize % groupsize == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % groupsize == 0);

        extern  __shared__ char externalSmem[];

        constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
        constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();
        
        constexpr int numRowsPSSM = 21;
        constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
        constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;
        using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
        
        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);

        using MathOps = MathOps<ScoreType>;
        using Scalar = typename ScalarScoreType<ScoreType>::type;

        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int numGroupsInBlock = blockDim.x / groupsize;
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;
        
        const size_t groupTempStorageOffset = idOfGroupInGrid * tempStorageElementsPerGroup;
        float2* const groupTempStorage = multiTileTempStorage + groupTempStorageOffset;
        
        const SequenceLengthT queryLength = inputData.getQueryLength();
        const int numTiles = SDIV(queryLength, 2 * groupsize * numRegs);

        using SmemIndexCalculator = PssmSmemIndexCalculator<groupsize, pssmReplicationFactor, relaxChunkSize>;

        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize>;

        alignas(8) Scalar penalty_in[4];
        alignas(8) Scalar penalty_out[4];

        auto load_PSSM_single = [&](int tileNr) {
            const int columnOffset = tileNr * pssmTileSize;
            __syncthreads(); //wait for all groups before overwriting pssm

            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                int letter = i/(pssmTileSize);
                int col = i%(pssmTileSize);
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][columnOffset + col];
            }
            __syncthreads();
        };

        auto load_PSSM_replicated = [&](int tileNr) {
            const int columnOffset = tileNr * pssmTileSize;
            __syncthreads(); //wait for all groups before overwriting pssm

            SmemIndexCalculator pssmSmemIndexCalculator;

            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                auto value = strided_PSSM[letter][columnOffset + col];

                auto repIndices = pssmSmemIndexCalculator.computeReplicatedColumnIndices(col);

                for(auto repIndex : repIndices){
                    shared_strided_PSSM.data[letter][repIndex] = value;
                }
            }
            __syncthreads();
        };

        auto load_PSSM = [&](int tileNr){
            if constexpr(SmemIndexCalculator::factor > 1){
                load_PSSM_replicated(tileNr);
            }else{
                load_PSSM_single(tileNr);
            }
        };


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

        const int numAlignments = inputData.getNumAlignments();
        //need to round up to blocks because loading pssm is a block-wide operation
        const int numAlignmentsRoundedUp = SDIV(numAlignments, numGroupsInBlock) * numGroupsInBlock;

        for(int alignmentId = idOfGroupInGrid; alignmentId < numAlignmentsRoundedUp; alignmentId += numGroupsInGrid){
            char4 new_subject_letter4;
            SequenceLengthT subjectLength;
            const char4* subjectAsChar4;

            State state(shared_strided_PSSM, group);

            //first tile
            {
                /* 
                    -----------------------
                    Process tile 0
                    ----------------------- 
                */

                //load pssm for tile 0. blockwide operation
                load_PSSM(0);

                if(alignmentId < numAlignments){
                    const int8_t* subjectData = inputData.getSubject(alignmentId);
                    subjectLength = inputData.getSubjectLength(alignmentId);
                    subjectAsChar4 = reinterpret_cast<const char4*>(subjectData);
            

                    state.resetScores();
                    state.resetMaximum();
                    
                    int k;

                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {

                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);

                        state.stepFirstTile(new_subject_letter4.x, penalty_out[0], zero);
                        state.stepFirstTile(new_subject_letter4.y, penalty_out[1], zero);
                        state.stepFirstTile(new_subject_letter4.z, penalty_out[2], zero);
                        state.stepFirstTile(new_subject_letter4.w, penalty_out[3], zero);
                        
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
                        }
                    }

                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        state.stepFirstTile(new_subject_letter4.x, penalty_out[0], zero);
                    }

                    if (subjectLength%4 >= 2) {
                        state.stepFirstTile(new_subject_letter4.y, penalty_out[1], zero);
                    }

                    if (subjectLength%4 >= 3) {
                        state.stepFirstTile(new_subject_letter4.z, penalty_out[2], zero);
                    }

                    //if there were remaining rows, update temp storage
                    if(subjectLength % 4 > 0){
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
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

                if(alignmentId < numAlignments){    
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepIntermediateTile(new_subject_letter4.x, penalty_in[0], penalty_out[0], zero);
                        state.stepIntermediateTile(new_subject_letter4.y, penalty_in[1], penalty_out[1], zero);
                        state.stepIntermediateTile(new_subject_letter4.z, penalty_in[2], penalty_out[2], zero);
                        state.stepIntermediateTile(new_subject_letter4.w, penalty_in[3], penalty_out[3], zero);
            
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
                        }
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }
                        state.stepIntermediateTile(new_subject_letter4.x, penalty_in[0], penalty_out[0], zero);
                    }
        
                    if (subjectLength%4 >= 2) {
                        state.stepIntermediateTile(new_subject_letter4.y, penalty_in[1], penalty_out[1], zero);
                    }
        
                    if (subjectLength%4 >= 3) {
                        state.stepIntermediateTile(new_subject_letter4.z, penalty_in[2], penalty_out[2], zero);
                    }
        
                    //if there were remaining rows, update temp storage
                    if(subjectLength % 4 > 0){
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
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

                if(alignmentId < numAlignments){
        
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepLastTile(new_subject_letter4.x, penalty_in[0], zero);
                        state.stepLastTile(new_subject_letter4.y, penalty_in[1], zero);
                        state.stepLastTile(new_subject_letter4.z, penalty_in[2], zero);
                        state.stepLastTile(new_subject_letter4.w, penalty_in[3], zero);
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }
                        state.stepLastTile(new_subject_letter4.x, penalty_in[0], zero);
                    }
        
                    if (subjectLength%4 >= 2) {
                        state.stepLastTile(new_subject_letter4.y, penalty_in[1], zero);
                    }
        
                    if (subjectLength%4 >= 3) {
                        state.stepLastTile(new_subject_letter4.z, penalty_in[2], zero);
                    }
                }
            }

            if(alignmentId < numAlignments){
                state.reduceMaximumScore();
                const float overall_max = MathOps::max(state.maximum.x, state.maximum.y);

                if(group.thread_rank() == 0){
                    devAlignmentScores[alignmentId] = overall_max;
                }
            }
        }

    }



} //namespace kernelparamzero










namespace half2fmarelu{

#if defined(__CUDACC__)
    template<class ScoreType, int numRegs, class Group, class SharedPSSM, class SmemIndexCalculator, int relaxChunkSize_>
    struct GaplessPSSMState{
        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOps = MathOps<ScoreType>;

        static constexpr int relaxChunkSize = 8; // 16-byte loads load 8 fp8x2 quantities
        static_assert(relaxChunkSize == 8);
        static_assert(relaxChunkSize_ == 8);

        static constexpr int numRelaxChunks = SDIV(numRegs, relaxChunkSize);
        static constexpr int numRegsInLastRelaxChunk = numRegs - (numRelaxChunks-1)*relaxChunkSize;

        ScoreType penalty_here_array[numRegs];
        ScoreType maximum{}; //0
        ScoreType penalty_diag{}; //0
        SharedPSSM& shared_strided_PSSM;
        Group& group;

        __device__
        GaplessPSSMState(SharedPSSM& s, Group& g) : shared_strided_PSSM(s), group(g) {}

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

        __device__
        half2 add_relu(half2 a, half2 b){
#if defined(__CUDACC__)
            return __hfma2_relu(a,make_half2(1,1),b);
#else
            const half2 one = __floats2half2_rn(1.f, 1.f);
            const half2 zero = __floats2half2_rn(0.f, 0.f);
            return compatibility::hmax2(__hfma2(a, one, b), zero);
#endif
        }

        __device__
        void relax(int subject_letter){
            SmemIndexCalculator smemIndexCalculator;

            ScoreType penalty_temp[2];

            const auto* row = &shared_strided_PSSM.data[subject_letter][0];

            float4 foo = *((float4*)&row[smemIndexCalculator.getIndexOfFirstElementInVectorType(0)]);
            __nv_fp8x2_e4m3 foofp8x2[relaxChunkSize];
            memcpy(&foofp8x2[0], &foo, sizeof(float4));



            penalty_temp[0] = penalty_here_array[0];
            penalty_here_array[0] = add_relu(penalty_diag, static_cast<half2>(foofp8x2[0]));

            auto xxx = static_cast<half2>(foofp8x2[0]);
            // printf("thread %d, k %d, load %d %d", group.thread_rank(), 0, int(xxx.x), int(xxx.y));

            constexpr int loopLimit1 = (numRelaxChunks == 1) ? numRegsInLastRelaxChunk : relaxChunkSize;
            #pragma unroll
            for(int k = 1; k < loopLimit1; k++){
                penalty_temp[k % 2] = penalty_here_array[k];
                penalty_here_array[k] = add_relu(penalty_temp[(k+1)%2], static_cast<half2>(foofp8x2[k]));
                xxx = static_cast<half2>(foofp8x2[k]);
                // printf("thread %d, k %d, load %f %f\n", group.thread_rank(), k, float(xxx.x), float(xxx.y));
            }

            #pragma unroll
            for(int i = 1; i < numRelaxChunks; i++){
                foo = *((float4*)&row[smemIndexCalculator.getIndexOfFirstElementInVectorType(i)]);
                memcpy(&foofp8x2[0], &foo, sizeof(float4));

                const int loopLimit2 = (i == numRelaxChunks-1) ? numRegsInLastRelaxChunk : relaxChunkSize;
                #pragma unroll
                for(int k = 0; k < relaxChunkSize; k++){
                    if(k < loopLimit2){
                        penalty_temp[k % 2] = penalty_here_array[relaxChunkSize*i + k];
                        penalty_here_array[relaxChunkSize*i + k] = add_relu(penalty_temp[(k+1)%2], static_cast<half2>(foofp8x2[k]));
                        xxx = static_cast<half2>(foofp8x2[k]);
                        // printf("thread %d, k %d, load %f %f\n", group.thread_rank(), relaxChunkSize*i + k, float(xxx.x), float(xxx.y));
                    }
                }
            }

            // printf("thread %d, subject_letter %d, computed (%d %d) (%d %d) \n", 
            //     group.thread_rank(), subject_letter ,
            //     int(penalty_here_array[0].x), int(penalty_here_array[0].y),
            //     int(penalty_here_array[1].x), int(penalty_here_array[1].y)
            // );

            #pragma unroll
            for(int k = 0; k < numRegs/2; k++){
                maximum = MathOps::max3(maximum, penalty_here_array[2*k+0], penalty_here_array[2*k+1]);
                // printf("thread %d, k %d, maximum %f %f\n", group.thread_rank(), k, float(maximum.x), float(maximum.y));
            }
            if constexpr(numRegs % 2 == 1){
                maximum = MathOps::max(maximum, penalty_here_array[numRegs-1]);
                // printf("thread %d, k %d, maximum %f %f\n", group.thread_rank(), numRegs-1, float(maximum.x), float(maximum.y));
            }
        };

        __device__
        void shuffleScores(const Scalar& border_in){
            penalty_diag = group.shfl_up(penalty_here_array[numRegs-1], 1);
            const ScoreType penalty_temp0 = compatibility::group_shfl_down(group, penalty_here_array[numRegs-1], group.size()-1);

            if (group.thread_rank() == 0) {
                penalty_diag.x = border_in;
                penalty_diag.y = penalty_temp0.x;
            }
        }

        __device__
        void stepSingleTile(int subject_letter, bool print = false){
            if(print && group.thread_rank() == 0){
                printf("pos %d, letter %d\n", -1, subject_letter);
            }
            relax(subject_letter);
            shuffleScores(Scalar{});
        }

        __device__
        void stepFirstTile(int subject_letter, Scalar& border_out){
            relax(subject_letter);
            shuffleScores(Scalar{});
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1].y;
            }
        }

        __device__
        void stepIntermediateTile(int subject_letter, const Scalar& border_in, Scalar& border_out){
            relax(subject_letter);
            shuffleScores(border_in);
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1].y;
            }
        }

        __device__
        void stepLastTile(int subject_letter, const Scalar& border_in){
            relax(subject_letter);
            shuffleScores(border_in);
        }

        __device__
        void reduceMaximumScore(){
            maximum = MathOps::reduce_max(group, maximum);
        }
    };
#endif

    /*
    PSSM kernel for a query of max length (2 * groupsize * numRegs)
    */
    template<
        class ScoreType, 
        int blocksize, 
        int groupsize, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PssmScoreType,
        class InputData
    >
    __global__
    __launch_bounds__(blocksize,1)
    void GaplessFilter_strided_PSSM_singletile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM
    ) {
#if defined(__CUDACC__)        
        static_assert(std::is_same_v<PssmScoreType, __nv_fp8x2_e4m3>);
        static_assert(std::is_same_v<ScoreType, half2>);

        static_assert(blocksize % groupsize == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % groupsize == 0);

        constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
        constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();

        constexpr int numRowsPSSM = 21;
        constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
        constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;

        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("device: pssmTileSize = %d, numColumnsPSSM = %d, groupsize %d, numRegs %d\n", 
        //         pssmTileSize, numColumnsPSSM, groupsize, numRegs);
        // }


        using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
        // if(threadIdx.x + blockIdx.x == 0){
        //     printf("device: smem pssm takes %lu bytes\n", sizeof(SharedPSSM));
        // }
        using MathOps = MathOps<ScoreType>;

        extern  __shared__ char externalSmem[];

        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);


        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        //const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;

        using SmemIndexCalculator = PssmSmemIndexCalculator<groupsize, pssmReplicationFactor, relaxChunkSize>;

        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize>;

        auto load_PSSM_single = [&]() {
            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][col];
            }
            __syncthreads();
        };

        auto load_PSSM_replicated = [&]() {
            SmemIndexCalculator pssmSmemIndexCalculator;
            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/(pssmTileSize);
                const int col = i%(pssmTileSize);
                auto value = strided_PSSM[letter][col];

                auto repIndices = pssmSmemIndexCalculator.computeReplicatedColumnIndices(col);

                for(auto repIndex : repIndices){
                    shared_strided_PSSM.data[letter][repIndex] = value;
                }
            }
            __syncthreads();

        };

        auto load_PSSM = [&](){
            // __syncthreads();
            // if(threadIdx.x == 0){
            //     printf("gmem pssm\n");
            //     for(int r = 0; r < 21; r++){
            //         for(int c = 0; c < numColumnsPSSM; c++){
            //             auto s = strided_PSSM[r][c];
            //             half2 hs = static_cast<half2>(s);
            //             printf("%d %d ", int(hs.x), int(hs.y));
            //         }
            //         printf("\n");
            //     }
            // }
            // __syncthreads();

            if constexpr(SmemIndexCalculator::factor > 1){
                load_PSSM_replicated();
            }else{
                load_PSSM_single();
            }

            // __syncthreads();
            // if(threadIdx.x == 0){
            //     printf("smem pssm\n");
            //     for(int r = 0; r < 21; r++){
            //         for(int c = 0; c < numColumnsPSSM; c++){
            //             auto s = shared_strided_PSSM.data[r][c];
            //             half2 hs = static_cast<half2>(s);
            //             printf("%d %d ", int(hs.x), int(hs.y));
            //         }
            //         printf("\n");
            //     }
            // }
            // __syncthreads();
        };


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

        const SequenceLengthT queryLength = inputData.getQueryLength();
        const int numAlignments = inputData.getNumAlignments();

        const int alignmentId = idOfGroupInGrid;
        if(alignmentId < numAlignments){
            char4 new_subject_letter4;
            const int8_t* subjectData = inputData.getSubject(alignmentId);
            const SequenceLengthT subjectLength = inputData.getSubjectLength(alignmentId);
            const char4* subjectAsChar4 = reinterpret_cast<const char4*>(subjectData);
            

            State state(shared_strided_PSSM, group);

            state.resetScores();
            state.resetMaximum();

            int k;
            for (k=0; k<subjectLength-3; k+=4) {
                new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                state.stepSingleTile(new_subject_letter4.x);
                state.stepSingleTile(new_subject_letter4.y);
                state.stepSingleTile(new_subject_letter4.z);
                state.stepSingleTile(new_subject_letter4.w);
            }

            if (subjectLength%4 >= 1) {
                new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                state.stepSingleTile(new_subject_letter4.x);
            }

            if (subjectLength%4 >= 2) {
                state.stepSingleTile(new_subject_letter4.y);
            }

            if (subjectLength%4 >= 3) {
                state.stepSingleTile(new_subject_letter4.z);
            }

            state.reduceMaximumScore();
            const float overall_max = MathOps::max(state.maximum.x, state.maximum.y);

            if(group.thread_rank() == 0){
                // printf("%f %f %f\n", overall_max, float(state.maximum.x), float(state.maximum.y));
                devAlignmentScores[alignmentId] = overall_max;
            }
        }
#endif        
    }




    /*
    PSSM kernel for arbitrary query length
    */
    template<
        class ScoreType, 
        int blocksize, 
        int groupsize, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PssmScoreType,
        class InputData
    >
    __global__
    __launch_bounds__(blocksize,1)
    void GaplessFilter_strided_PSSM_multitile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM,
        __grid_constant__ float2* const multiTileTempStorage,
        __grid_constant__ const size_t tempStorageElementsPerGroup
    ) {
#if defined(__CUDACC__)
        static_assert(std::is_same_v<PssmScoreType, __nv_fp8x2_e4m3>);
        static_assert(std::is_same_v<ScoreType, half2>);

        static_assert(blocksize % groupsize == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % groupsize == 0);



        
        constexpr int relaxChunkSize = 16 / sizeof(PssmScoreType);
        constexpr int pssmReplicationFactor = PssmReplicationFactor<groupsize, PssmScoreType, relaxChunkSize>::value();
        
        constexpr int numRowsPSSM = 21;
        constexpr int pssmTileSize = groupsize * relaxChunkSize * SDIV(numRegs, relaxChunkSize);
        constexpr int numColumnsPSSM = pssmTileSize * pssmReplicationFactor;
        using SharedPSSM = SharedPSSM_singletile<PssmScoreType, numRowsPSSM, numColumnsPSSM>;
        
        extern  __shared__ char externalSmem[];
        SharedPSSM& shared_strided_PSSM = *((SharedPSSM*)externalSmem);

        using MathOps = MathOps<ScoreType>;
        using Scalar = typename ScalarScoreType<ScoreType>::type;

        auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
        const int numGroupsInBlock = blockDim.x / groupsize;
        const int idOfGroupInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
        const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;
        
        const size_t groupTempStorageOffset = idOfGroupInGrid * tempStorageElementsPerGroup;
        float2* const groupTempStorage = multiTileTempStorage + groupTempStorageOffset;
        
        const SequenceLengthT queryLength = inputData.getQueryLength();
        const int numTiles = SDIV(queryLength, 2 * groupsize * numRegs);

        using SmemIndexCalculator = PssmSmemIndexCalculator<groupsize, pssmReplicationFactor, relaxChunkSize>;

        using State = GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize>;

        alignas(8) Scalar penalty_in[4];
        alignas(8) Scalar penalty_out[4];

        auto load_PSSM_single = [&](int tileNr) {
            const int columnOffset = tileNr * pssmTileSize;
            __syncthreads(); //wait for all groups before overwriting pssm

            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                int letter = i/pssmTileSize;
                int col = i%pssmTileSize;
                shared_strided_PSSM.data[letter][col] = strided_PSSM[letter][columnOffset + col];
            }
            __syncthreads();
        };

        auto load_PSSM_replicated = [&](int tileNr) {
            SmemIndexCalculator pssmSmemIndexCalculator;
            const int columnOffset = tileNr * pssmTileSize;
            __syncthreads(); //wait for all groups before overwriting pssm

            for (int i=threadIdx.x; i<21*pssmTileSize; i+=blockDim.x) {
                const int letter = i/pssmTileSize;
                const int col = i%pssmTileSize;
                auto value = strided_PSSM[letter][columnOffset + col];

                auto repIndices = pssmSmemIndexCalculator.computeReplicatedColumnIndices(col);

                for(auto repIndex : repIndices){
                    shared_strided_PSSM.data[letter][repIndex] = value;
                }
            }
            __syncthreads();
        };

        auto load_PSSM = [&](int tileNr){
            if constexpr(SmemIndexCalculator::factor > 1){
                load_PSSM_replicated(tileNr);
            }else{
                load_PSSM_single(tileNr);
            }
        };



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

        const int numAlignments = inputData.getNumAlignments();
        //need to round up to blocks because loading pssm is a block-wide operation
        const int numAlignmentsRoundedUp = SDIV(numAlignments, numGroupsInBlock) * numGroupsInBlock;

        for(int alignmentId = idOfGroupInGrid; alignmentId < numAlignmentsRoundedUp; alignmentId += numGroupsInGrid){
            char4 new_subject_letter4;
            SequenceLengthT subjectLength;
            const char4* subjectAsChar4;

            State state(shared_strided_PSSM, group);

            //first tile
            {
                /* 
                    -----------------------
                    Process tile 0
                    ----------------------- 
                */

                //load pssm for tile 0. blockwide operation
                load_PSSM(0);

                if(alignmentId < numAlignments){
                    const int8_t* subjectData = inputData.getSubject(alignmentId);
                    subjectLength = inputData.getSubjectLength(alignmentId);
                    subjectAsChar4 = reinterpret_cast<const char4*>(subjectData);
                    

                    state.resetScores();
                    state.resetMaximum();
                    
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
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
                        }
                    }

                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        state.stepFirstTile(new_subject_letter4.x, penalty_out[0]);
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
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
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

                if(alignmentId < numAlignments){    
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepIntermediateTile(new_subject_letter4.x, penalty_in[0], penalty_out[0]);
                        state.stepIntermediateTile(new_subject_letter4.y, penalty_in[1], penalty_out[1]);
                        state.stepIntermediateTile(new_subject_letter4.z, penalty_in[2], penalty_out[2]);
                        state.stepIntermediateTile(new_subject_letter4.w, penalty_in[3], penalty_out[3]);
            
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
                        }
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
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
                            groupTempStorage[k/4] = *((float2*)&penalty_out[0]);
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

                if(alignmentId < numAlignments){
        
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
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
                            *((float2*)&penalty_in[0]) = groupTempStorage[k/4];
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

            if(alignmentId < numAlignments){
                state.reduceMaximumScore();
                const float overall_max = MathOps::max(state.maximum.x, state.maximum.y);

                if(group.thread_rank() == 0){
                    devAlignmentScores[alignmentId] = overall_max;
                }
            }
        }
#endif
    }


} //namespace half2fmarelu







LIBMARV_NAMESPACE_WITH_NESTING_END

#endif