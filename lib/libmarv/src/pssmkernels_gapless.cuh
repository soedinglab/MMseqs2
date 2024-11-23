#ifndef PSSM_KERNELS_CUH
#define PSSM_KERNELS_CUH

#include <cuda_fp16.h>

#include <map>

#include "pssm.cuh"
#include "convert.cuh"
#include "mathops.cuh"
#include "util.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;


#define USE_IMPROVED_SMEM


#define PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X \
    X(4,4) X(4,8) X(4,12) X(4,16) X(4,20) X(4,24) X(4,28) X(4,32) \
    X(4,36) X(4,40) X(4,44) X(4,48) X(4,52) X(4,56) X(4,60) X(4,64) \
    X(8,4) X(8,8) X(8,12) X(8,16) X(8,20) X(8,24) X(8,28) X(8,32) \
    X(8,36) X(8,40) X(8,44) X(8,48) X(8,52) X(8,56) X(8,60) X(8,64) \
    X(16,4) X(16,8) X(16,12) X(16,16) X(16,20) X(16,24) X(16,28) X(16,32) \
    X(16,36) X(16,40) X(16,44) X(16,48) X(16,52) X(16,56) X(16,60) X(16,64)

#define PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X \
    X(8,4) X(8,8) X(8,12) X(8,16) X(8,20) X(8,24) X(8,28) X(8,32) \
    X(8,36) X(8,40) X(8,44) X(8,48) X(8,52) X(8,56) X(8,60) X(8,64) \
    X(16,4) X(16,8) X(16,12) X(16,16) X(16,20) X(16,24) X(16,28) X(16,32) \
    X(16,36) X(16,40) X(16,44) X(16,48) X(16,52) X(16,56) X(16,60) X(16,64)


namespace cudasw4{




namespace hardcodedzero{

    template<class ScoreType> struct ScalarScoreType{};
    template<> struct ScalarScoreType<half2>{ using type = half; };
    template<> struct ScalarScoreType<short2>{ using type = short; };
    template<> struct ScalarScoreType<int>{ using type = int; };
    template<> struct ScalarScoreType<float>{ using type = float; };

    template<class ScoreType, int numRegs, class Group, class SharedPSSM, class SmemIndexCalculator>
    struct GaplessPSSMState{
        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOps = MathOps<ScoreType>;

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
        void relax(int subject_letter){
            SmemIndexCalculator smemIndexCalculator;

            ScoreType score2;
            ScoreType penalty_temp0;
            ScoreType penalty_temp1;

            const auto* row = &shared_strided_PSSM.data[subject_letter][0];

            float4 foo = *((float4*)&row[smemIndexCalculator.getIndex(0)]);
            memcpy(&score2, &foo.x, sizeof(ScoreType));
            penalty_temp0 = penalty_here_array[0];
            penalty_here_array[0] = MathOps::add_relu(penalty_diag, score2);

            memcpy(&score2, &foo.y, sizeof(ScoreType));
            penalty_temp1 = penalty_here_array[1];
            penalty_here_array[1] = MathOps::add_relu(penalty_temp0, score2);
            maximum = MathOps::max3(maximum, penalty_here_array[1], penalty_here_array[0]);

            memcpy(&score2, &foo.z, sizeof(ScoreType));
            penalty_temp0 = penalty_here_array[2];
            penalty_here_array[2] = MathOps::add_relu(penalty_temp1, score2);

            memcpy(&score2, &foo.w, sizeof(ScoreType));
            penalty_temp1 = penalty_here_array[3];
            penalty_here_array[3] = MathOps::add_relu(penalty_temp0, score2);
            maximum = MathOps::max3(maximum, penalty_here_array[3], penalty_here_array[2]);


            #pragma unroll
            for (int i=1; i<numRegs/4; i++) {
                foo = *((float4*)&row[smemIndexCalculator.getIndex(i)]);
                memcpy(&score2, &foo.x, sizeof(ScoreType));
                penalty_temp0 = penalty_here_array[4*i];
                penalty_here_array[4*i] = MathOps::add_relu(penalty_temp1, score2);

                memcpy(&score2, &foo.y, sizeof(ScoreType));
                penalty_temp1 = penalty_here_array[4*i+1];
                penalty_here_array[4*i+1] = MathOps::add_relu(penalty_temp0, score2);
                maximum = MathOps::max3(maximum, penalty_here_array[4*i+1], penalty_here_array[4*i]);

                memcpy(&score2, &foo.z, sizeof(ScoreType));
                penalty_temp0 = penalty_here_array[4*i+2];
                penalty_here_array[4*i+2] = MathOps::add_relu(penalty_temp1, score2);

                memcpy(&score2, &foo.w, sizeof(ScoreType));
                penalty_temp1 = penalty_here_array[4*i+3];
                penalty_here_array[4*i+3] = MathOps::add_relu(penalty_temp0, score2);
                maximum = MathOps::max3(maximum, penalty_here_array[4*i+3], penalty_here_array[4*i+2]);
            }
        };

        __device__
        void shuffleScores(const Scalar& border_in){
            penalty_diag = group.shfl_up(penalty_here_array[numRegs-1], 1);
            const ScoreType penalty_temp0 = group.shfl_down(penalty_here_array[numRegs-1], group.size()-1);

            if (group.thread_rank() == 0) {
                penalty_diag.x = border_in;
                penalty_diag.y = penalty_temp0.x;
            }
        }

        __device__
        void stepSingleTile(int subject_letter){
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


    /*
    PSSM kernel for a query of max length (2 * group_size * numRegs)
    */
    template<
        class ScoreType, 
        int blocksize, 
        int group_size, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PositionsIterator
    >
    __global__
    __launch_bounds__(512,1)
    void GaplessFilter_strided_PSSM_singletile_kernel(
        __grid_constant__ const char * const devChars,
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const size_t* const devOffsets,
        __grid_constant__ const SequenceLengthT* const devLengths,
        __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
        __grid_constant__ const int numSelected,
        __grid_constant__ const SequenceLengthT queryLength,
        __grid_constant__ const PSSM_2D_View<ScoreType> strided_PSSM
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(numRegs % 4 == 0);
        static_assert(blocksize % group_size == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % group_size == 0);

        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM
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

        #ifdef USE_IMPROVED_SMEM
        using SmemIndexCalculator = typename std::conditional<
            group_size == 4, 
            SmemIndexCalculator<group_size, 2>,
            SmemIndexCalculator<group_size, 1>
        >::type;
        #else
        using SmemIndexCalculator = SmemIndexCalculator<group_size, 1>;
        #endif
        GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator> state(shared_strided_PSSM, group);

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

            state.resetScores();
            state.resetMaximum();

            subjectAsChar4 = reinterpret_cast<const char4*>(&devChars[base_S]);

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
                devAlignmentScores[alignmentId] = overall_max;
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
        class ScoreOutputIterator, 
        class PositionsIterator
    >
    void call_GaplessFilter_strided_PSSM_singletile_kernel(
        const char * const devChars,
        ScoreOutputIterator const devAlignmentScores,
        const size_t* const devOffsets,
        const SequenceLengthT* const devLengths,
        PositionsIterator const d_positions_of_selected_lengths,
        const int numSelected,
        const SequenceLengthT queryLength,
        const PSSM_2D_View<ScoreType>& strided_PSSM,
        cudaStream_t stream
    ){
        constexpr int groupsPerBlock = blocksize / group_size;
        constexpr int alignmentsPerGroup = 1;
        constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
        // std::cout << "blocksize " << blocksize << ", group_size " << group_size 
        //     << ", alignmentsPerBlock " << alignmentsPerBlock << ", numSelected " << numSelected << "\n";

        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM
        constexpr int numColumnsPSSM = std::max(group_size,8) * numRegs;
        #else
        constexpr int numColumnsPSSM = group_size * numRegs;
        #endif
        using SharedPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;
        
        int smem = sizeof(SharedPSSM);
        auto kernel = GaplessFilter_strided_PSSM_singletile_kernel<
            ScoreType,
            blocksize, 
            group_size, 
            numRegs, 
            subjectIsCaseSensitive,
            ScoreOutputIterator, 
            PositionsIterator>;

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
            strided_PSSM
        ); CUERR;
    }






    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_singletile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_singletile_kernel<short2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<short2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator




    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_singletile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_singletile_kernel<short2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<short2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator



    template<class ScoreType, int blocksize, class ScoreOutputIterator, class PositionsIterator>
    void call_GaplessFilter_strided_PSSM_singletile_kernel(
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
        cudaStream_t stream
    ){
        constexpr bool subjectIsCaseSensitive = true;

        #define X(g,r) \
            if(group_size == g && numRegs == r){ \
                call_GaplessFilter_strided_PSSM_singletile_kernel<ScoreType, blocksize,g,r,subjectIsCaseSensitive>( \
                    devChars, devAlignmentScores, devOffsets, devLengths, d_positions_of_selected_lengths, \
                    numSelected, queryLength, strided_PSSM, stream \
                ); \
            } else 

            PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X
        { throw std::runtime_error("invalid groupsize/numregs config");}

        #undef X
    }















    /*
    PSSM kernel for arbitrary query length
    */
    template<
        class ScoreType, 
        int blocksize, 
        int group_size, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PositionsIterator
    >
    __global__
    __launch_bounds__(512,1)
    void GaplessFilter_strided_PSSM_multitile_kernel(
        __grid_constant__ const char * const devChars,
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const size_t* const devOffsets,
        __grid_constant__ const SequenceLengthT* const devLengths,
        __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
        __grid_constant__ const int numSelected,
        __grid_constant__ const SequenceLengthT queryLength,
        __grid_constant__ const PSSM_2D_View<ScoreType> strided_PSSM,
        __grid_constant__ float2* const multiTileTempStorage,
        __grid_constant__ const size_t tempStorageElementsPerGroup
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(numRegs % 4 == 0);
        static_assert(blocksize % group_size == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % group_size == 0);

        extern  __shared__ char externalSmem[];
        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM
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
        float2* const groupTempStorage = multiTileTempStorage + groupTempStorageOffset;
        
        const int numTiles = SDIV(queryLength, 2 * group_size * numRegs);

        #ifdef USE_IMPROVED_SMEM
        using SmemIndexCalculator = typename std::conditional<
            group_size == 4, 
            SmemIndexCalculator<group_size, 2>,
            SmemIndexCalculator<group_size, 1>
        >::type;
        #else
        using SmemIndexCalculator = SmemIndexCalculator<group_size, 1>;
        #endif
        GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator> state(shared_strided_PSSM, group);

        alignas(8) Scalar penalty_in[4];
        alignas(8) Scalar penalty_out[4];

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

                if(alignmentId < numSelected){    
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

                if(alignmentId < numSelected){
        
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

            if(alignmentId < numSelected){
                state.reduceMaximumScore();
                const float overall_max = MathOps::max(state.maximum.x, state.maximum.y);

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
        class ScoreOutputIterator,
        class PositionsIterator
    >
    void call_GaplessFilter_strided_PSSM_multitile_kernel(
        int numThreadBlocks,
        const char * const devChars,
        ScoreOutputIterator const devAlignmentScores,
        const size_t* const devOffsets,
        const SequenceLengthT* const devLengths,
        PositionsIterator const d_positions_of_selected_lengths,
        const int numSelected,
        const SequenceLengthT queryLength,
        const PSSM_2D_View<ScoreType>& strided_PSSM,
        float2* const multiTileTempStorage,
        size_t tempStorageElementsPerGroup, //number of float2s per group
        cudaStream_t stream
    ){
        //constexpr int groupsPerBlock = blocksize / group_size;
        //constexpr int alignmentsPerGroup = 1;
        //constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
        // std::cout << "blocksize " << blocksize << ", group_size " << group_size 
        //     << ", alignmentsPerBlock " << alignmentsPerBlock << ", numSelected " << numSelected << "\n";

        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM
        constexpr int numColumnsPSSM = std::max(group_size,8) * numRegs;
        #else
        constexpr int numColumnsPSSM = group_size * numRegs;
        #endif
        using SharedPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;

        int smem = sizeof(SharedPSSM);
        auto kernel = GaplessFilter_strided_PSSM_multitile_kernel<ScoreType, blocksize, group_size, numRegs, subjectIsCaseSensitive,
            ScoreOutputIterator, PositionsIterator>;

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
            multiTileTempStorage,
            tempStorageElementsPerGroup
        ); CUERR;
    }




    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_multitile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_multitile_kernel<short2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<short2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator


    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_multitile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_multitile_kernel<short2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<short2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator


    template<
        class ScoreType,
        int blocksize, 
        class ScoreOutputIterator, 
        class PositionsIterator
    >
    void call_GaplessFilter_strided_PSSM_multitile_kernel(
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
        float2* const multiTileTempStorage,
        size_t tempStorageElementsPerGroup,
        cudaStream_t stream
    ){
        constexpr bool subjectIsCaseSensitive = true;

        #define X(g,r) \
            if(group_size == g && numRegs == r){ \
                call_GaplessFilter_strided_PSSM_multitile_kernel<ScoreType, blocksize,g,r,subjectIsCaseSensitive>( \
                    numThreadBlocks, devChars, devAlignmentScores, devOffsets, devLengths, \
                    d_positions_of_selected_lengths, numSelected, queryLength, strided_PSSM, \
                    multiTileTempStorage, tempStorageElementsPerGroup, stream \
                ); \
            } else 

            PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X
        { throw std::runtime_error("invalid groupsize/numregs config");}

        #undef X
    }

} //namespace hardcodedzero


namespace kernelparamzero{    

    template<class ScoreType> struct ScalarScoreType{};
    template<> struct ScalarScoreType<half2>{ using type = half; };
    template<> struct ScalarScoreType<short2>{ using type = short; };
    template<> struct ScalarScoreType<int>{ using type = int; };
    template<> struct ScalarScoreType<float>{ using type = float; };

    template<class ScoreType, int numRegs, class Group, class SharedPSSM, class SmemIndexCalculator>
    struct GaplessPSSMState{
        using Scalar = typename ScalarScoreType<ScoreType>::type;
        using MathOps = MathOps<ScoreType>;

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

            ScoreType score2;
            ScoreType penalty_temp0;
            ScoreType penalty_temp1;

            const auto* row = &shared_strided_PSSM.data[subject_letter][0];

            float4 foo = *((float4*)&row[smemIndexCalculator.getIndex(0)]);
            memcpy(&score2, &foo.x, sizeof(ScoreType));
            penalty_temp0 = penalty_here_array[0];
            penalty_here_array[0] = MathOps::add_relu(penalty_diag, score2, zero);

            memcpy(&score2, &foo.y, sizeof(ScoreType));
            penalty_temp1 = penalty_here_array[1];
            penalty_here_array[1] = MathOps::add_relu(penalty_temp0, score2, zero);
            maximum = MathOps::max3(maximum, penalty_here_array[1], penalty_here_array[0]);

            memcpy(&score2, &foo.z, sizeof(ScoreType));
            penalty_temp0 = penalty_here_array[2];
            penalty_here_array[2] = MathOps::add_relu(penalty_temp1, score2, zero);

            memcpy(&score2, &foo.w, sizeof(ScoreType));
            penalty_temp1 = penalty_here_array[3];
            penalty_here_array[3] = MathOps::add_relu(penalty_temp0, score2, zero);
            maximum = MathOps::max3(maximum, penalty_here_array[3], penalty_here_array[2]);


            #pragma unroll
            for (int i=1; i<numRegs/4; i++) {
                foo = *((float4*)&row[smemIndexCalculator.getIndex(i)]);
                memcpy(&score2, &foo.x, sizeof(ScoreType));
                penalty_temp0 = penalty_here_array[4*i];
                penalty_here_array[4*i] = MathOps::add_relu(penalty_temp1, score2, zero);

                memcpy(&score2, &foo.y, sizeof(ScoreType));
                penalty_temp1 = penalty_here_array[4*i+1];
                penalty_here_array[4*i+1] = MathOps::add_relu(penalty_temp0, score2, zero);
                maximum = MathOps::max3(maximum, penalty_here_array[4*i+1], penalty_here_array[4*i]);

                memcpy(&score2, &foo.z, sizeof(ScoreType));
                penalty_temp0 = penalty_here_array[4*i+2];
                penalty_here_array[4*i+2] = MathOps::add_relu(penalty_temp1, score2, zero);

                memcpy(&score2, &foo.w, sizeof(ScoreType));
                penalty_temp1 = penalty_here_array[4*i+3];
                penalty_here_array[4*i+3] = MathOps::add_relu(penalty_temp0, score2, zero);
                maximum = MathOps::max3(maximum, penalty_here_array[4*i+3], penalty_here_array[4*i+2]);
            }
        };

        __device__
        void shuffleScores(const Scalar& border_in){
            penalty_diag = group.shfl_up(penalty_here_array[numRegs-1], 1);
            const ScoreType penalty_temp0 = group.shfl_down(penalty_here_array[numRegs-1], group.size()-1);

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
    PSSM kernel for a query of max length (2 * group_size * numRegs)
    */
    template<
        class ScoreType, 
        int blocksize, 
        int group_size, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PositionsIterator
    >
    __global__
    __launch_bounds__(512,1)
    void GaplessFilter_strided_PSSM_singletile_kernel(
        __grid_constant__ const char * const devChars,
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const size_t* const devOffsets,
        __grid_constant__ const SequenceLengthT* const devLengths,
        __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
        __grid_constant__ const int numSelected,
        __grid_constant__ const SequenceLengthT queryLength,
        __grid_constant__ const PSSM_2D_View<ScoreType> strided_PSSM,
        __grid_constant__ const ScoreType zero
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(numRegs % 4 == 0);
        static_assert(blocksize % group_size == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % group_size == 0);

        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM
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

        #ifdef USE_IMPROVED_SMEM
        using SmemIndexCalculator = typename std::conditional<
            group_size == 4, 
            SmemIndexCalculator<group_size, 2>,
            SmemIndexCalculator<group_size, 1>
        >::type;
        #else
        using SmemIndexCalculator = SmemIndexCalculator<group_size, 1>;
        #endif
        GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator> state(shared_strided_PSSM, group);

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
        };

        auto load_PSSM = [&](){
            if constexpr(SmemIndexCalculator::factor == 2){
                load_PSSM_double();
            }else{
                load_PSSM_single();
            }
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

            state.resetScores();
            state.resetMaximum();

            subjectAsChar4 = reinterpret_cast<const char4*>(&devChars[base_S]);

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
    PSSM kernel for a query of max length (2 * group_size * numRegs)
    */
    template<
        class ScoreType,
        int blocksize, 
        int group_size, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator, 
        class PositionsIterator
    >
    void call_GaplessFilter_strided_PSSM_singletile_kernel(
        const char * const devChars,
        ScoreOutputIterator const devAlignmentScores,
        const size_t* const devOffsets,
        const SequenceLengthT* const devLengths,
        PositionsIterator const d_positions_of_selected_lengths,
        const int numSelected,
        const SequenceLengthT queryLength,
        const PSSM_2D_View<ScoreType>& strided_PSSM,
        cudaStream_t stream
    ){
        constexpr int groupsPerBlock = blocksize / group_size;
        constexpr int alignmentsPerGroup = 1;
        constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
        // std::cout << "blocksize " << blocksize << ", group_size " << group_size 
        //     << ", alignmentsPerBlock " << alignmentsPerBlock << ", numSelected " << numSelected << "\n";

        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM
        constexpr int numColumnsPSSM = std::max(group_size,8) * numRegs;
        #else
        constexpr int numColumnsPSSM = group_size * numRegs;
        #endif
        using SharedPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;
        
        int smem = sizeof(SharedPSSM);
        auto kernel = GaplessFilter_strided_PSSM_singletile_kernel<
            ScoreType,
            blocksize, 
            group_size, 
            numRegs, 
            subjectIsCaseSensitive,
            ScoreOutputIterator, 
            PositionsIterator>;

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
            MathOps<ScoreType>::zero_score()
        ); CUERR;
    }






    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_singletile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_singletile_kernel<short2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<short2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator




    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_singletile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_singletile_kernel<short2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<short2>&, \
            cudaStream_t \
        );

        PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator



    template<class ScoreType, int blocksize, class ScoreOutputIterator, class PositionsIterator>
    void call_GaplessFilter_strided_PSSM_singletile_kernel(
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
        cudaStream_t stream
    ){
        constexpr bool subjectIsCaseSensitive = true;

        #define X(g,r) \
            if(group_size == g && numRegs == r){ \
                call_GaplessFilter_strided_PSSM_singletile_kernel<ScoreType, blocksize,g,r,subjectIsCaseSensitive>( \
                    devChars, devAlignmentScores, devOffsets, devLengths, d_positions_of_selected_lengths, \
                    numSelected, queryLength, strided_PSSM, stream \
                ); \
            } else 

            PSSM_GAPLESS_SINGLETILE_FOR_EACH_VALID_CONFIG_DO_X
        { throw std::runtime_error("invalid groupsize/numregs config");}

        #undef X
    }















    /*
    PSSM kernel for arbitrary query length
    */
    template<
        class ScoreType, 
        int blocksize, 
        int group_size, 
        int numRegs, 
        bool subjectIsCaseSensitive, 
        class ScoreOutputIterator,
        class PositionsIterator
    >
    __global__
    __launch_bounds__(512,1)
    void GaplessFilter_strided_PSSM_multitile_kernel(
        __grid_constant__ const char * const devChars,
        __grid_constant__ ScoreOutputIterator const devAlignmentScores,
        __grid_constant__ const size_t* const devOffsets,
        __grid_constant__ const SequenceLengthT* const devLengths,
        __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
        __grid_constant__ const int numSelected,
        __grid_constant__ const SequenceLengthT queryLength,
        __grid_constant__ const PSSM_2D_View<ScoreType> strided_PSSM,
        __grid_constant__ float2* const multiTileTempStorage,
        __grid_constant__ const size_t tempStorageElementsPerGroup,
        __grid_constant__ const ScoreType zero
    ) {
        if constexpr (std::is_same_v<ScoreType, short2>) {
            #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
            return;
            #endif
        }
        static_assert(numRegs % 4 == 0);
        static_assert(blocksize % group_size == 0);
        __builtin_assume(blockDim.x == blocksize);
        __builtin_assume(blockDim.x % group_size == 0);

        extern  __shared__ char externalSmem[];
        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM
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
        float2* const groupTempStorage = multiTileTempStorage + groupTempStorageOffset;
        
        const int numTiles = SDIV(queryLength, 2 * group_size * numRegs);

        #ifdef USE_IMPROVED_SMEM
        using SmemIndexCalculator = typename std::conditional<
            group_size == 4, 
            SmemIndexCalculator<group_size, 2>,
            SmemIndexCalculator<group_size, 1>
        >::type;
        #else
        using SmemIndexCalculator = SmemIndexCalculator<group_size, 1>;
        #endif
        GaplessPSSMState<ScoreType, numRegs, decltype(group), SharedPSSM, SmemIndexCalculator> state(shared_strided_PSSM, group);

        alignas(8) Scalar penalty_in[4];
        alignas(8) Scalar penalty_out[4];

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

                    state.resetScores();
                    state.resetMaximum();
                    subjectAsChar4 = reinterpret_cast<const char4*>(&devChars[base_S]);

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

                if(alignmentId < numSelected){    
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

                if(alignmentId < numSelected){
        
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

            if(alignmentId < numSelected){
                state.reduceMaximumScore();
                const float overall_max = MathOps::max(state.maximum.x, state.maximum.y);

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
        class ScoreOutputIterator,
        class PositionsIterator
    >
    void call_GaplessFilter_strided_PSSM_multitile_kernel(
        int numThreadBlocks,
        const char * const devChars,
        ScoreOutputIterator const devAlignmentScores,
        const size_t* const devOffsets,
        const SequenceLengthT* const devLengths,
        PositionsIterator const d_positions_of_selected_lengths,
        const int numSelected,
        const SequenceLengthT queryLength,
        const PSSM_2D_View<ScoreType>& strided_PSSM,
        float2* const multiTileTempStorage,
        size_t tempStorageElementsPerGroup, //number of float2s per group
        cudaStream_t stream
    ){
        //constexpr int groupsPerBlock = blocksize / group_size;
        //constexpr int alignmentsPerGroup = 1;
        //constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
        // std::cout << "blocksize " << blocksize << ", group_size " << group_size 
        //     << ", alignmentsPerBlock " << alignmentsPerBlock << ", numSelected " << numSelected << "\n";

        constexpr int numRowsPSSM = 21;
        #ifdef USE_IMPROVED_SMEM
        constexpr int numColumnsPSSM = std::max(group_size,8) * numRegs;
        #else
        constexpr int numColumnsPSSM = group_size * numRegs;
        #endif
        using SharedPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;

        int smem = sizeof(SharedPSSM);
        auto kernel = GaplessFilter_strided_PSSM_multitile_kernel<ScoreType, blocksize, group_size, numRegs, subjectIsCaseSensitive,
            ScoreOutputIterator, PositionsIterator>;

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
            multiTileTempStorage,
            tempStorageElementsPerGroup,
            MathOps<ScoreType>::zero_score()
        ); CUERR;
    }



    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_multitile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator decltype(thrust::make_counting_iterator<ReferenceIdT>(0))
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_multitile_kernel<short2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<short2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator


    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_multitile_kernel<half2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<half2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator

    #define ScoreOutputIterator TopNMaximaArray
    #define PositionsIterator ReferenceIdT*
    #define subjectIsCaseSensitive true
    #define X(g,r) \
        extern template void call_GaplessFilter_strided_PSSM_multitile_kernel<short2, 512, g, r, subjectIsCaseSensitive, ScoreOutputIterator, PositionsIterator>( \
            int, \
            const char * const, \
            ScoreOutputIterator const, \
            const size_t* const, \
            const SequenceLengthT* const, \
            PositionsIterator const, \
            const int, \
            const SequenceLengthT, \
            const PSSM_2D_View<short2>&, \
            float2*, \
            size_t, \
            cudaStream_t \
        );

        PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X

    #undef X
    #undef subjectIsCaseSensitive
    #undef PositionsIterator
    #undef ScoreOutputIterator


    template<
        class ScoreType,
        int blocksize, 
        class ScoreOutputIterator, 
        class PositionsIterator
    >
    void call_GaplessFilter_strided_PSSM_multitile_kernel(
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
        float2* const multiTileTempStorage,
        size_t tempStorageElementsPerGroup,
        cudaStream_t stream
    ){
        constexpr bool subjectIsCaseSensitive = true;

        #define X(g,r) \
            if(group_size == g && numRegs == r){ \
                call_GaplessFilter_strided_PSSM_multitile_kernel<ScoreType, blocksize,g,r,subjectIsCaseSensitive>( \
                    numThreadBlocks, devChars, devAlignmentScores, devOffsets, devLengths, \
                    d_positions_of_selected_lengths, numSelected, queryLength, strided_PSSM, \
                    multiTileTempStorage, tempStorageElementsPerGroup, stream \
                ); \
            } else 

            PSSM_GAPLESS_MULTITILE_FOR_EACH_VALID_CONFIG_DO_X
        { throw std::runtime_error("invalid groupsize/numregs config");}

        #undef X
    }


} //namespace kernelparamzero


} //namespace cudasw4

#endif