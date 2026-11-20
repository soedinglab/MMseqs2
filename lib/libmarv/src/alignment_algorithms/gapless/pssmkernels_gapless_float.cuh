#ifndef PSSM_KERNELS_GAPLESS_FLOAT_CUH
#define PSSM_KERNELS_GAPLESS_FLOAT_CUH

#include "../../cuda_hip_rename.h"

#include <map>

#include "../../pssm.cuh"
#include "../../convert.cuh"
#include "../../mathops.cuh"
#include "../../util.cuh"
#include "../../offset_iterator.cuh"
#include "../../custom_score_types.cuh"
#include "../../ptx_wrappers.cuh"
#include "gapless_common.cuh"
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







    template<
        class ScoreType, 
        int numRegs, 
        class Group, 
        class SharedPSSM, 
        class SmemIndexCalculator,
        int relaxChunkSize_,
        class UpdateMaxOp
    >
    struct GaplessAlignmentState_scalar{
        using Scalar = ScoreType;
        using MathOps = MathOps<ScoreType>;
        using PssmScoreType = typename SharedPSSM::value_type;

        static constexpr int relaxChunkSize = relaxChunkSize_;

        static constexpr int numRelaxChunks = SDIV(numRegs, relaxChunkSize);
        static constexpr int numRegsInLastRelaxChunk = numRegs - (numRelaxChunks-1)*relaxChunkSize;

        ScoreType penalty_here_array[numRegs];
        ScoreType penalty_diag{}; //0
        SharedPSSM& shared_strided_PSSM;
        Group& group;
        UpdateMaxOp updateMaxOp;

        __device__
        GaplessAlignmentState_scalar(SharedPSSM& s, Group& g, UpdateMaxOp updateMaxOp_) 
            : shared_strided_PSSM(s), group(g), updateMaxOp(updateMaxOp_) {

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
        void relax(int subject_letter, int row, int tileNr){
            SmemIndexCalculator smemIndexCalculator;

            ScoreType penalty_temp[2];

            const auto* pssmptr = &shared_strided_PSSM.data[subject_letter][0];            

            float4 foo = *((float4*)&pssmptr[smemIndexCalculator.getIndexOfFirstElementInVectorType(0)]);

            PssmScoreType foo2[relaxChunkSize];
            memcpy(&foo2[0], &foo, sizeof(float4));

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

            // printf("row %d, penalty_here_array: %d %d %d %d\n",
            //     row,
            //     penalty_here_array[0], penalty_here_array[1], penalty_here_array[2], penalty_here_array[3]);

            updateMaxOp(penalty_here_array, tileNr, row);

        };



        __device__
        void shuffleScores(const Scalar& border_in){
            penalty_diag = group.shfl_up(penalty_here_array[numRegs-1], 1);
            const ScoreType penalty_temp0 = compatibility::group_shfl_down(group, penalty_here_array[numRegs-1], group.size()-1);

            if (group.thread_rank() == 0) {
                penalty_diag = border_in;
            }
        }

        __device__
        void stepSingleTile(int subject_letter, int row){
            relax(subject_letter, row, 0);
            shuffleScores(Scalar{});
        }

        __device__
        void stepFirstTile(int subject_letter, int row, Scalar& border_out){
            relax(subject_letter, row, 0);
            shuffleScores(Scalar{});
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1];
            }

            // if(group.thread_rank() == group.size() - 1){
            //     printf("thread %d, penalty_out = %d\n", threadIdx.x, border_out);
            // }
        }

        __device__
        void stepIntermediateTile(int subject_letter, int row, int tileNr, const Scalar& border_in, Scalar& border_out){
            relax(subject_letter, row, tileNr);
            shuffleScores(border_in);
            if(group.thread_rank() == group.size() - 1){
                border_out = penalty_here_array[numRegs-1];
            }
            // if(group.thread_rank() == group.size() - 1){
            //     printf("thread %d, penalty_out = %d\n", threadIdx.x, border_out);
            // }
        }

        __device__
        void stepLastTile(int subject_letter, int row, int tileNr, const Scalar& border_in){
            relax(subject_letter, row, tileNr);
            shuffleScores(border_in);
        }
    };


    /*
    PSSM kernel for a query of max length (groupsize * numRegs)
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
    void GaplessAlignment_PSSM_scalar_singletile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ int* const subjectEndPositions_exclusive, //not accessed if withEndPosition == false
        __grid_constant__ int* const queryEndPositions_exclusive, //not accessed if withEndPosition == false
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM

    ) {

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

        using UpdateMaxOp = typename std::conditional<
            withEndPosition,
            UpdateMax_gapless_scalar_endPos_singleTile<ScoreType>,
            UpdateMax_gapless_scoreOnly<ScoreType>
        >::type;
        using State = GaplessAlignmentState_scalar<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize, UpdateMaxOp>;

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
            //             printf("%d ", int(s));
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

                int reg = maximumTracker.positionOfMaxObserved_reg;
                int regInGroup = group.thread_rank() * numRegs + reg;
                int queryPos = regInGroup;
                int subjectPos = maximumTracker.positionOfMaxObserved_row;
                int maximum = maximumTracker.maximum;
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
                int myMaximum = maximumTracker.maximum;
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
    void GaplessAlignment_PSSM_scalar_multitile_kernel(
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ int* const subjectEndPositions_exclusive, //not accessed if withEndPosition == false
        __grid_constant__ int* const queryEndPositions_exclusive, //not accessed if withEndPosition == false
        __grid_constant__ const InputData inputData,
        __grid_constant__ const PSSM_2D_View<PssmScoreType> strided_PSSM,
        __grid_constant__ float4* const multiTileTempStorage,
        __grid_constant__ const size_t tempStorageElementsPerGroup
    ) {

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
        float4* const groupTempStorage = multiTileTempStorage + groupTempStorageOffset;
        
        const SequenceLengthT queryLength = inputData.getQueryLength();
        const int numTiles = SDIV(queryLength, groupsize * numRegs);

        using SmemIndexCalculator = PssmSmemIndexCalculator<groupsize, pssmReplicationFactor, relaxChunkSize>;

        using UpdateMaxOp = typename std::conditional<
            withEndPosition,
            UpdateMax_gapless_scalar_endPos_multiTile<ScoreType>,
            UpdateMax_gapless_scoreOnly<ScoreType>
        >::type;
        using State = GaplessAlignmentState_scalar<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator, relaxChunkSize, UpdateMaxOp>;


        alignas(16) Scalar penalty_in[4]{};
        alignas(16) Scalar penalty_out[4]{};

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

            SequenceLengthT subjectLength;
            const char4* subjectAsChar4;
            char4 new_subject_letter4;

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
                            groupTempStorage[k/4] = *((float4*)&penalty_out[0]);
                        }
                    }

                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        state.stepFirstTile(new_subject_letter4.x, k, penalty_out[0]);
                        // if(group.thread_rank() == group.size() - 1){
                        //     printf("thread %d, penalty_out[0] = %d\n", threadIdx.x, penalty_out[0]);
                        // }
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
                            groupTempStorage[k/4] = *((float4*)&penalty_out[0]);
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
                            *((float4*)&penalty_in[0]) = groupTempStorage[k/4];
                        }

                        state.stepIntermediateTile(new_subject_letter4.x, k, tileNr, penalty_in[0], penalty_out[0]);
                        state.stepIntermediateTile(new_subject_letter4.y, k+1, tileNr, penalty_in[1], penalty_out[1]);
                        state.stepIntermediateTile(new_subject_letter4.z, k+2, tileNr, penalty_in[2], penalty_out[2]);
                        state.stepIntermediateTile(new_subject_letter4.w, k+3, tileNr, penalty_in[3], penalty_out[3]);
            
                        //update temp storage for next tile
                        if(group.thread_rank() == group.size() - 1){
                            groupTempStorage[k/4] = *((float4*)&penalty_out[0]);
                        }
                    }
        
                    //process at most 3 remaining rows
                    if (subjectLength%4 >= 1) {
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        //load input penalty for remaining rows
                        if (group.thread_rank() == 0){
                            *((float4*)&penalty_in[0]) = groupTempStorage[k/4];
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
                            groupTempStorage[k/4] = *((float4*)&penalty_out[0]);
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
                load_PSSM(tileNr);

                if(alignmentId < numAlignments){
        
                    state.resetScores();
        
                    int k;
        
                    //process rows in chunks of 4 rows
                    for (k=0; k<subjectLength-3; k+=4) {
        
                        new_subject_letter4 = makeCaseInsensitive4(subjectAsChar4[k/4]);
                        if (group.thread_rank() == 0){
                            *((float4*)&penalty_in[0]) = groupTempStorage[k/4];
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
                            *((float4*)&penalty_in[0]) = groupTempStorage[k/4];
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

                    int reg = maximumTracker.positionOfMaxObserved_reg;
                    int regInGroup = group.thread_rank() * numRegs + reg;
                    int queryPos = maximumTracker.positionOfMaxObserved_tileNr * (group.num_threads() * numRegs)
                        + regInGroup;
                    int subjectPos = maximumTracker.positionOfMaxObserved_row;
                    int maximum = maximumTracker.maximum;
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
                    int myMaximum = maximumTracker.maximum;
               
                    const int alignmentScore = compatibility::group_reduce_max(group, myMaximum);

                    if(group.thread_rank() == 0){
                        devAlignmentScores[alignmentId] = alignmentScore;
                    }
                }
            }
        }

    }






LIBMARV_NAMESPACE_WITH_NESTING_END

#endif