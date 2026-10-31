#ifndef PSSM_KERNELS_SMITH_WATERMAN_CUH
#define PSSM_KERNELS_SMITH_WATERMAN_CUH

#include "../../cuda_hip_rename.h"

#include "../../config.hpp"
#include "../../pssm.cuh"
#include "../../convert.cuh"
#include "../../mathops.cuh"
#include "../../util.cuh"
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



template<
    class ScoreType,
    int blocksize, 
    int groupsize, 
    int numItems, 
    bool withEndPosition,
    bool subjectIsCaseSensitive,
    class InputData
>
__global__
__launch_bounds__(blocksize,1)
void amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_singletile(
    __grid_constant__ int* const devAlignmentScores,
    __grid_constant__ int* const subjectEndPositions_exclusive, //not accessed if withEndPositions == false
    __grid_constant__ int* const queryEndPositions_exclusive, //not accessed if withEndPositions == false
    __grid_constant__ const InputData inputData,
    __grid_constant__ const PSSM_2D_View<ScoreType> strided_PSSM,
    __grid_constant__ const ScoreType gapopenscore, 
    __grid_constant__ const ScoreType gapextendscore
){
    if constexpr (std::is_same_v<ScoreType, int>) {
        #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
        return;
        #endif
    }
    static_assert(std::is_same_v<ScoreType, float> || std::is_same_v<ScoreType, int>);

    static_assert(groupsize >= 4);
    static_assert(groupsize <= 32);
    static_assert(blocksize % groupsize == 0);

    __builtin_assume(blockDim.x == blocksize);
    __builtin_assume(blockDim.x % groupsize == 0);
    __builtin_assume(groupsize >= 4);
    __builtin_assume(groupsize <= 32);
    
    auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
    
    const int groupIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
    const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;

    constexpr ScoreType oobscore = -9999999;
    constexpr int badLetter = 20;

    constexpr int numRowsPSSM = 21;
    constexpr int numColumnsPSSM = groupsize * numItems;

    using MathOps = MathOps<ScoreType>;

    using SPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;
    extern __shared__ float4 externalSharedMem[];
    SPSSM& shared_pssm = *((SPSSM*)((char*)&externalSharedMem[0]));

    auto load_PSSM = [&](){
        for(int i = threadIdx.x; i < 21 * groupsize * numItems; i += blockDim.x){
            const int row = i / (groupsize * numItems);
            const int col = i % (groupsize * numItems);
            shared_pssm.data[row][col] = strided_PSSM[row][col];
        }
        __syncthreads();
    };

    auto makeCaseInsensitive4 = [](char4 encoded4){
        unsigned int asUint;
        memcpy(&asUint, &encoded4, sizeof(unsigned int));

        if constexpr(subjectIsCaseSensitive){
            asUint = CaseSensitive_to_CaseInsensitive{}(asUint);
            // asUint = ClampToInvalid{}(asUint);
        }

        memcpy(&encoded4, &asUint, sizeof(unsigned int));
        return encoded4;
    };

    //load PSSM to smem
    for(int i = threadIdx.x; i < 21 * groupsize * numItems; i += blockDim.x){
        const int row = i / (groupsize * numItems);
        const int col = i % (groupsize * numItems);
        shared_pssm.data[row][col] = strided_PSSM[row][col];
    }
    __syncthreads();

    const SequenceLengthT queryLength = inputData.getQueryLength();
    const int numAlignments = inputData.getNumAlignments();
    for(int alignmentId = groupIdInGrid; alignmentId < numAlignments; alignmentId += numGroupsInGrid){

        ScoreType scoresF[numItems]{};
        ScoreType scoresM[numItems]{};
        ScoreType scoreLeft;
        ScoreType scoreDiag;
        ScoreType E;
        ScoreType maxObserved = oobscore;
        int positionOfMaxObserved_y = 0;
        int positionOfMaxObserved_itemIndex = 0;


        auto printState = [&](int row){
            // for(int g = 0; g < 1; g++){
            //     if(groupIdInGrid == g){
            //         printf("printstate row %d, groupIdInGrid %d\n", row, groupIdInGrid);
            //         if(group.thread_rank() == 0){
            //             printf("M\n");
            //         }
            //         for(int t = 0; t < groupsize; t++){
            //             if(t == group.thread_rank()){
            //                 for(int i = 0; i < numItems; i++){
            //                     printf("%3f ", scoresM[i]);
            //                 }
            //                 printf("\n");
            //             }
            //             group.sync();
            //         }
            //         if(group.thread_rank() == 0){
            //             printf("\n");
            //         }

            //         if(group.thread_rank() == 0){
            //             printf("F\n");
            //         }
            //         for(int t = 0; t < groupsize; t++){
            //             if(t == group.thread_rank()){
            //                 for(int i = 0; i < numItems; i++){
            //                     printf("%3f ", scoresF[i]);
            //                 }
            //                 printf("\n");
            //             }
            //             group.sync();
            //         }
            //         if(group.thread_rank() == 0){
            //             printf("\n");
            //         }
            //     }
            // }
        };

        SequenceLengthT subjectLength = inputData.getSubjectLength(alignmentId);
        const int8_t* subjectData = inputData.getSubject(alignmentId);
        const char4* groupSubjectChar4 = reinterpret_cast<const char4*>(subjectData);

        int loadOffsetLimit = SDIV(subjectLength, 4);
        int subjectLoadOffset = group.thread_rank();                
        char4 current4Letters;
        int currentLetter = badLetter;

        auto loadNext4Letters = [&](){
            if(subjectLoadOffset < loadOffsetLimit){
                current4Letters = makeCaseInsensitive4(groupSubjectChar4[subjectLoadOffset]);
                subjectLoadOffset += group.size();
            }else{
                current4Letters = makeCaseInsensitive4(make_char4(badLetter, badLetter, badLetter, badLetter));
            }
        };

        auto shuffleCurrentLetter = [&](){
            currentLetter = group.shfl_up(currentLetter, 1);
        };

        auto shuffle4Letters = [&](){
            static_assert(sizeof(char4) == sizeof(int));
            int temp;
            memcpy(&temp, &current4Letters, sizeof(char4));
            temp = compatibility::group_shfl_down(group, temp, 1);
            memcpy(&current4Letters, &temp, sizeof(int));
        };

        auto relaxFirstDiagonal = [&](int row, bool isFirstTile){
            static_assert(numItems % 4 == 0);

            using Vec4T = typename Vectorized4<ScoreType>::type;
            const Vec4T* const pssmRow4 = reinterpret_cast<const Vec4T*>(&shared_pssm.data[currentLetter][0]);
            //const Vec4T* const pssmRow4 = reinterpret_cast<const Vec4T*>(&strided_PSSM[currentLetter][0]);
            Vec4T foo = pssmRow4[0 * groupsize + group.thread_rank()];
            ScoreType fooArray[4];
            memcpy(&fooArray[0], &foo, sizeof(Vec4T));

            //in the first tile E is always computed. In succeeding tiles, E is already computed for the first thread (loaded from temp storage)
            if(isFirstTile){
                E = MathOps::add_max(scoreLeft, gapopenscore, MathOps::add(E, gapextendscore));
            }else{
                if(group.thread_rank() > 0){
                    E = MathOps::add_max(scoreLeft, gapopenscore, MathOps::add(E, gapextendscore));
                }
            }

            scoresF[0] = MathOps::add_max(scoresM[0], gapopenscore, MathOps::add(scoresF[0], gapextendscore));
            ScoreType upTempScore = scoresM[0];
            scoresM[0] = MathOps::add_max_relu(scoreDiag, fooArray[0], MathOps::max(E, scoresF[0]));
            scoreDiag = upTempScore;
            if constexpr (withEndPosition){
                if(maxObserved < scoresM[0]){
                    maxObserved = scoresM[0];
                    positionOfMaxObserved_itemIndex = 0;
                    positionOfMaxObserved_y = row;
                }
            }else{
                maxObserved = MathOps::max(maxObserved, scoresM[0]);
            }


            #pragma unroll
            for(int k = 1; k < 4; k++){
                E = MathOps::add_max(scoresM[k-1], gapopenscore, MathOps::add(E, gapextendscore));
                scoresF[k] = MathOps::add_max(scoresM[k], gapopenscore, MathOps::add(scoresF[k], gapextendscore));
                ScoreType upTempScore = scoresM[k];
                scoresM[k] = MathOps::add_max_relu(scoreDiag, fooArray[k], MathOps::max(E, scoresF[k]));
                scoreDiag = upTempScore;
                if constexpr (withEndPosition){
                    if(maxObserved < scoresM[k]){
                        maxObserved = scoresM[k];
                        positionOfMaxObserved_itemIndex = k;
                        positionOfMaxObserved_y = row;
                    }
                }else{
                    maxObserved = MathOps::max(maxObserved, scoresM[k]);
                }
            }

            #pragma unroll
            for(int i = 1; i < numItems/4; i++){
                foo = pssmRow4[i * group.size() + group.thread_rank()];
                memcpy(&fooArray[0], &foo, sizeof(ScoreType) * 4);

                #pragma unroll
                for(int k = 0; k < 4; k++){
                    E = MathOps::add_max(scoresM[4*i + k-1], gapopenscore, MathOps::add(E, gapextendscore));
                    scoresF[4*i + k] = MathOps::add_max(scoresM[4*i + k], gapopenscore, MathOps::add(scoresF[4*i + k], gapextendscore));
                    ScoreType upTempScore = scoresM[4*i + k];
                    scoresM[4*i + k] = MathOps::add_max_relu(scoreDiag, fooArray[k], MathOps::max(E, scoresF[4*i + k]));
                    scoreDiag = upTempScore;
                    if constexpr (withEndPosition){
                        if(maxObserved < scoresM[4*i + k]){
                            maxObserved = scoresM[4*i + k];
                            positionOfMaxObserved_itemIndex = 4*i + k;
                            positionOfMaxObserved_y = row;
                        }
                    }else{
                        maxObserved = MathOps::max(maxObserved, scoresM[4*i + k]);
                    }
                }
            }

            //advance E by 1 column and F by 1 row to allow for optimized computations of remaining diagonals
            E = MathOps::add_max(scoresM[numItems-1], gapopenscore, MathOps::add(E, gapextendscore));
            for(int k = 0; k < numItems; k++){
                scoresF[k] = MathOps::add_max(scoresM[k], gapopenscore, MathOps::add(scoresF[k], gapextendscore));
            }

            //printState(row);
        };

        auto relax = [&](int row){
            static_assert(numItems % 4 == 0);

            using Vec4T = typename Vectorized4<ScoreType>::type;
            const Vec4T* const pssmRow4 = reinterpret_cast<const Vec4T*>(&shared_pssm.data[currentLetter][0]);
            //const Vec4T* const pssmRow4 = reinterpret_cast<const Vec4T*>(&strided_PSSM[currentLetter][0]);
            Vec4T foo = pssmRow4[0 * groupsize + group.thread_rank()];
            ScoreType fooArray[4];
            memcpy(&fooArray[0], &foo, sizeof(Vec4T));

            // E of current column and scoresF of current row are already computed

            ScoreType tempM = scoresM[0];
            scoresM[0] = MathOps::add_max_relu(scoreDiag, fooArray[0], MathOps::max(E, scoresF[0]));     
            if constexpr (withEndPosition){
                if(maxObserved < scoresM[0]){
                    maxObserved = scoresM[0];
                    positionOfMaxObserved_itemIndex = 0;
                    positionOfMaxObserved_y = row;
                }
            }else{
                maxObserved = MathOps::max(maxObserved, scoresM[0]);
            }
            E = MathOps::add_max(scoresM[0], gapopenscore, MathOps::add(E, gapextendscore));
            scoresF[0] = MathOps::add_max(scoresM[0], gapopenscore, MathOps::add(scoresF[0],gapextendscore)); //this computes F of the next row !
            scoreDiag = tempM;

            #pragma unroll
            for(int i = 1; i < 4; i++){
                tempM = scoresM[i];
                scoresM[i] = MathOps::add_max_relu(scoreDiag, fooArray[i], MathOps::max(E, scoresF[i]));
                if constexpr (withEndPosition){
                    if(maxObserved < scoresM[i]){
                        maxObserved = scoresM[i];
                        positionOfMaxObserved_itemIndex = i;
                        positionOfMaxObserved_y = row;
                    }
                }else{
                    maxObserved = MathOps::max(maxObserved, scoresM[i]);
                }
                E = MathOps::add_max(scoresM[i], gapopenscore, MathOps::add(E, gapextendscore));
                scoresF[i] = MathOps::add_max(scoresM[i], gapopenscore, MathOps::add(scoresF[i], gapextendscore)); //this computes F of the next row !
                scoreDiag = tempM;
            }

            #pragma unroll
            for(int k = 1; k < numItems/4; k++){
                foo = pssmRow4[k * groupsize + group.thread_rank()];
                memcpy(&fooArray[0], &foo, sizeof(Vec4T));

                #pragma unroll
                for(int i = 0; i < 4; i++){
                    const int index = k*4+i;
                    tempM = scoresM[index];
                    scoresM[index] = MathOps::add_max_relu(scoreDiag, fooArray[i], MathOps::max(E, scoresF[index]));
                    if constexpr (withEndPosition){
                        if(maxObserved < scoresM[index]){
                            maxObserved = scoresM[index];
                            positionOfMaxObserved_itemIndex = index;
                            positionOfMaxObserved_y = row;
                        }
                    }else{
                        maxObserved = MathOps::max(maxObserved, scoresM[index]);
                    }
                    E = MathOps::add_max(scoresM[index], gapopenscore, MathOps::add(E, gapextendscore));
                    scoresF[index] = MathOps::add_max(scoresM[index], gapopenscore, MathOps::add(scoresF[index], gapextendscore)); //this computes F of the next row !
                    scoreDiag = tempM;
                }

            }

            //printState(row);
        };

        auto initScores = [&](){
            if(group.thread_rank() == 0){
                #pragma unroll
                for (int i=0; i < numItems; i++) {
                    scoresM[i] = 0;
                    scoresF[i] = oobscore;
                }
                scoreDiag = 0;
                scoreLeft = 0;
                E = oobscore;
            }else{
                #pragma unroll
                for (int i=0; i < numItems; i++) {
                    scoresM[i] = oobscore;
                    scoresF[i] = oobscore;
                }
                scoreDiag = oobscore;
                scoreLeft = group.thread_rank() == 1 ? 0 : oobscore;
                E = oobscore;
            }
        };

        auto shuffleScores = [&](){
            scoreDiag = scoreLeft;
            const ScoreType newscoreLeft = group.shfl_up(scoresM[numItems-1], 1);
            const ScoreType newE = group.shfl_up(E, 1);
            if(group.thread_rank() == 0){
                //scoreLeft is only modified in this function and is initialized with 0 for thread 0
                // assert(scoreLeft == 0);
                //scoreLeft = 0;

                // E = oobscore;
                E = gapopenscore; // After first diagonal was processed, thread 0 needs E of matrix column 1, not -infty
            }else{
                scoreLeft = newscoreLeft;
                E = newE;
            }
        };

        loadNext4Letters();
        initScores();

        const int outputThreadRank = (queryLength-1) / numItems;
        const int numRows = subjectLength + outputThreadRank + 1;


        //printState(0);

        //process 4 letters per iteration
        int r = 1;
        if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
        constexpr bool isFirstTile = true;
        relaxFirstDiagonal(r, isFirstTile); //x
        shuffleScores();

        if(r+1 < numRows){
            shuffleCurrentLetter();
            if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
            relax(r+1); //y
            shuffleScores();
        }
        if(r+2 < numRows){
            shuffleCurrentLetter();
            if(group.thread_rank() == 0){ currentLetter = current4Letters.z; }
            relax(r+2); //z
            shuffleScores();        
        }
        if(r+3 < numRows){
            shuffleCurrentLetter();

            if(group.thread_rank() == 0){ currentLetter = current4Letters.w; }
            relax(r+3); //w
            shuffleScores();

            shuffleCurrentLetter(); 
            if((r + 3) % (4*group.size()) == 0){
                //used up all query letters stored across the group. reload
                loadNext4Letters();
            }else{
                //get next 4 letters from neighbor
                shuffle4Letters();
            }
        }
        r = 5;
        for(; r < numRows - 3; r += 4){

            if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
            relax(r); //x
            shuffleScores();

            shuffleCurrentLetter();
            if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
            relax(r+1); //y
            shuffleScores();

            shuffleCurrentLetter();
            if(group.thread_rank() == 0){ currentLetter = current4Letters.z; }
            relax(r+2); //z
            shuffleScores();

            shuffleCurrentLetter();
            if(group.thread_rank() == 0){ currentLetter = current4Letters.w; }
            relax(r+3); //w
            shuffleScores();

            shuffleCurrentLetter(); 
            if((r + 3) % (4*group.size()) == 0){
                //used up all query letters stored across the group. reload
                loadNext4Letters();
            }else{
                //get next 4 letters from neighbor
                shuffle4Letters();
            }     
        }

        //can have at most 3 remaining rows
        if(r < numRows){
            if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }   
            relax(r); //x
            shuffleScores();
            shuffleCurrentLetter();
            
        }
        if(r+1 < numRows){
            if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
            relax(r+1); //y
            shuffleScores();
            shuffleCurrentLetter();
            
        }
        if(r+2 < numRows){
            if(group.thread_rank() == 0){ currentLetter = current4Letters.z; }
            relax(r+2); //z
        }

        if constexpr (withEndPosition){
            if(alignmentId < numAlignments){
                const int myQueryEndInclusive = group.thread_rank() * numItems + positionOfMaxObserved_itemIndex;
                const int mySubjectEndInclusive = positionOfMaxObserved_y - group.thread_rank() - 1;

                //disable any maxima that come from out-of-bounds cells
                if(myQueryEndInclusive >= queryLength){
                    maxObserved = -1;
                }
                if(mySubjectEndInclusive < 0 || mySubjectEndInclusive >= subjectLength){
                    maxObserved = -1;
                }
                const int3 packed = make_int3(maxObserved, 
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
            }
        }else{
            if(alignmentId < numAlignments){
                const int alignmentScore = compatibility::group_reduce_max(group, maxObserved);

                if(group.thread_rank() == 0){
                    devAlignmentScores[alignmentId] = alignmentScore;
                }
            }
        }
    }

}





template<
    class ScoreType,
    int blocksize, 
    int groupsize, 
    int numItems, 
    bool withEndPosition,
    bool subjectIsCaseSensitive,
    class InputData
>
__global__
__launch_bounds__(blocksize,1)
void amino_gpu_localAlignmentKernel_affinegap_floatOrInt_pssm_multitile(
    __grid_constant__ int* const devAlignmentScores,
    __grid_constant__ int* const subjectEndPositions_exclusive, //not accessed if withEndPositions == false
    __grid_constant__ int* const queryEndPositions_exclusive, //not accessed if withEndPositions == false
    __grid_constant__ const InputData inputData,
    __grid_constant__ const PSSM_2D_View<ScoreType> strided_PSSM,
    __grid_constant__ const ScoreType gapopenscore, 
    __grid_constant__ const ScoreType gapextendscore,
    __grid_constant__ char* const tempStorage,
    __grid_constant__ const size_t tempBytesPerGroup
){
    if constexpr (std::is_same_v<ScoreType, int>) {
        #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 900
        return;
        #endif
    }
    static_assert(std::is_same_v<ScoreType, float> || std::is_same_v<ScoreType, int>);

    static_assert(groupsize >= 4);
    static_assert(groupsize <= 32);
    static_assert(blocksize % groupsize == 0);

    __builtin_assume(blockDim.x == blocksize);
    __builtin_assume(blockDim.x % groupsize == 0);
    __builtin_assume(groupsize >= 4);
    __builtin_assume(groupsize <= 32);
    
    auto group = cg::tiled_partition<groupsize>(cg::this_thread_block());
    
    const int groupIdInGrid = (threadIdx.x + blockIdx.x * blockDim.x) / groupsize;
    const int numGroupsInGrid = (blockDim.x * gridDim.x) / groupsize;
    constexpr int numGroupsPerBlock = blocksize / groupsize;

    constexpr ScoreType oobscore = -9999999;
    constexpr int badLetter = 20;

    constexpr int numRowsPSSM = 21;
    constexpr int numColumnsPSSM = groupsize * numItems;

    using MathOps = MathOps<ScoreType>;
    using SPSSM = SharedPSSM_singletile<ScoreType, numRowsPSSM, numColumnsPSSM>;
    extern __shared__ float4 externalSharedMem[];
    SPSSM& shared_pssm = *((SPSSM*)((char*)&externalSharedMem[0]));

    auto load_PSSM = [&](int tileNr){
        __syncthreads();
        // if(threadIdx.x == 0){
        //     printf("load_PSSM tileNr %d\n", tileNr);
        // }
        const int columnOffset = tileNr * groupsize * numItems;
        for(int i = threadIdx.x; i < 21 * groupsize * numItems; i += blockDim.x){
            const int row = i / (groupsize * numItems);
            const int col = i % (groupsize * numItems);
            shared_pssm.data[row][col] = strided_PSSM[row][columnOffset + col];
        }
        __syncthreads();
        // if(blockIdx.x == 0 && threadIdx.x == 0){
        //     printf("in gmem:\n");
        //     for(int i = 0; i < groupsize * numItems; i++){
        //         printf("%3f ", strided_PSSM[0][columnOffset + i]);
        //     }
        //     printf("\n");

        //     printf("in smem:\n");
        //     for(int i = 0; i < groupsize * numItems; i++){
        //         printf("%3f ", shared_pssm.data[0][i]);
        //     }
        //     printf("\n");
        // }
        // __syncthreads();
    };

    auto makeCaseInsensitive4 = [](char4 encoded4){
        unsigned int asUint;
        memcpy(&asUint, &encoded4, sizeof(unsigned int));

        if constexpr(subjectIsCaseSensitive){
            asUint = CaseSensitive_to_CaseInsensitive{}(asUint);
            // asUint = ClampToInvalid{}(asUint);
        }

        memcpy(&encoded4, &asUint, sizeof(unsigned int));
        return encoded4;
    };

    using Vec2T = typename Vectorized2<ScoreType>::type;

    Vec2T* groupTempStorage = (Vec2T*)(((char*)tempStorage) + tempBytesPerGroup * groupIdInGrid);

    auto clearOutOfTileTempStorage = [&](int subjectLength){
        if(group.thread_rank() < group.size() - 1){
            groupTempStorage[subjectLength + group.thread_rank()] = Vec2T{};
        }
    };

    const SequenceLengthT queryLength = inputData.getQueryLength();
    const int numAlignments = inputData.getNumAlignments();
    const int numAlignmentsRounded = SDIV(numAlignments, numGroupsPerBlock) * numGroupsPerBlock;
    const int numTiles = SDIV(queryLength, groupsize * numItems);

    for(int alignmentId = groupIdInGrid; alignmentId < numAlignmentsRounded; alignmentId += numGroupsInGrid){

        ScoreType scoresF[numItems]{};
        ScoreType scoresM[numItems]{};
        ScoreType scoreLeft;
        ScoreType scoreDiag;
        ScoreType E;
        ScoreType maxObserved = oobscore;
        int positionOfMaxObserved_y = 0;
        int positionOfMaxObserved_tileNr = 0;
        int positionOfMaxObserved_itemIndex = 0;

        Vec2T tileLastColumnM_E;
        Vec2T leftBorderM_E;

        // #define PRINT_WRITE
        // #define PRINT_LOAD

        auto printState = [&](int row){
            // for(int g = 0; g < 1; g++){
            //     if(groupIdInGrid == g){
            //         printf("printstate row %d, groupIdInGrid %d\n", row, groupIdInGrid);
            //         if(group.thread_rank() == 0){
            //             printf("M\n");
            //         }
            //         for(int t = 0; t < groupsize; t++){
            //             if(t == group.thread_rank()){
            //                 for(int i = 0; i < numItems; i++){
            //                     printf("%3f ", scoresM[i]);
            //                 }
            //                 printf("\n");
            //             }
            //             group.sync();
            //         }
            //         if(group.thread_rank() == 0){
            //             printf("\n");
            //         }

            //         if(group.thread_rank() == 0){
            //             printf("F\n");
            //         }
            //         for(int t = 0; t < groupsize; t++){
            //             if(t == group.thread_rank()){
            //                 for(int i = 0; i < numItems; i++){
            //                     printf("%3f ", scoresF[i]);
            //                 }
            //                 printf("\n");
            //             }
            //             group.sync();
            //         }
            //         if(group.thread_rank() == 0){
            //             printf("\n");
            //         }
            //     }
            // }
        };

        SequenceLengthT subjectLength = 0;
        const char4* groupSubjectChar4 = nullptr;
        int loadOffsetLimit = 0;
        int subjectLoadOffset = 0;
        char4 current4Letters;
        int currentLetter;
        int tempLoadOffset = 0;
        int tempWriteOffset = 0;

        auto loadNext4Letters = [&](){
            if(subjectLoadOffset < loadOffsetLimit){
                current4Letters = makeCaseInsensitive4(groupSubjectChar4[subjectLoadOffset]);
                subjectLoadOffset += group.size();
            }else{
                current4Letters = makeCaseInsensitive4(make_char4(badLetter, badLetter, badLetter, badLetter));
            }
        };

        auto shuffleCurrentLetter = [&](){
            currentLetter = group.shfl_up(currentLetter, 1);
        };

        auto shuffle4Letters = [&](){
            static_assert(sizeof(char4) == sizeof(int));
            int temp;
            memcpy(&temp, &current4Letters, sizeof(char4));
            temp = compatibility::group_shfl_down(group, temp, 1);
            memcpy(&current4Letters, &temp, sizeof(int));
        };

        auto setTileLastColumn = [&](){
            if(group.thread_rank() == group.size() - 1){
                tileLastColumnM_E.x = scoresM[numItems-1];
                tileLastColumnM_E.y = E;
            }
        };

        auto shuffleTileLastColumn = [&](){
            tileLastColumnM_E = compatibility::group_shfl_down(group, tileLastColumnM_E, 1);
        };
        auto shuffleLeftBorder = [&](){
            leftBorderM_E = compatibility::group_shfl_down(group, leftBorderM_E, 1);
        };

        auto relaxFirstDiagonal = [&](int row, int tileNr){
            static_assert(numItems % 4 == 0);

            using Vec4T = typename Vectorized4<ScoreType>::type;
            const Vec4T* const pssmRow4 = reinterpret_cast<const Vec4T*>(&shared_pssm.data[currentLetter][0]);
            //const Vec4T* const pssmRow4 = reinterpret_cast<const Vec4T*>(&strided_PSSM[currentLetter][0]);
            Vec4T foo = pssmRow4[0 * groupsize + group.thread_rank()];
            ScoreType fooArray[4];
            memcpy(&fooArray[0], &foo, sizeof(Vec4T));

            //in the first tile E is always computed. In succeeding tiles, E is already computed for the first thread (loaded from temp storage)
            if(tileNr == 0){
                E = MathOps::add_max(scoreLeft, gapopenscore, MathOps::add(E, gapextendscore));
            }else{
                if(group.thread_rank() > 0){
                    E = MathOps::add_max(scoreLeft, gapopenscore, MathOps::add(E, gapextendscore));
                }
            }

            scoresF[0] = MathOps::add_max(scoresM[0], gapopenscore, MathOps::add(scoresF[0], gapextendscore));
            ScoreType upTempScore = scoresM[0];
            scoresM[0] = MathOps::add_max_relu(scoreDiag, fooArray[0], MathOps::max(E, scoresF[0]));
            scoreDiag = upTempScore;
            if constexpr (withEndPosition){
                if(maxObserved < scoresM[0]){
                    maxObserved = scoresM[0];
                    positionOfMaxObserved_tileNr = tileNr;
                    positionOfMaxObserved_itemIndex = 0;
                    positionOfMaxObserved_y = row;
                }
            }else{
                maxObserved = MathOps::max(maxObserved, scoresM[0]);
            }


            #pragma unroll
            for(int k = 1; k < 4; k++){
                E = MathOps::add_max(scoresM[k-1], gapopenscore, MathOps::add(E, gapextendscore));
                scoresF[k] = MathOps::add_max(scoresM[k], gapopenscore, MathOps::add(scoresF[k], gapextendscore));
                ScoreType upTempScore = scoresM[k];
                scoresM[k] = MathOps::add_max_relu(scoreDiag, fooArray[k], MathOps::max(E, scoresF[k]));
                scoreDiag = upTempScore;
                if constexpr (withEndPosition){
                    if(maxObserved < scoresM[k]){
                        maxObserved = scoresM[k];
                        positionOfMaxObserved_tileNr = tileNr;
                        positionOfMaxObserved_itemIndex = k;
                        positionOfMaxObserved_y = row;
                    }
                }else{
                    maxObserved = MathOps::max(maxObserved, scoresM[k]);
                }
            }

            #pragma unroll
            for(int i = 1; i < numItems/4; i++){
                foo = pssmRow4[i * group.size() + group.thread_rank()];
                memcpy(&fooArray[0], &foo, sizeof(ScoreType) * 4);

                #pragma unroll
                for(int k = 0; k < 4; k++){
                    E = MathOps::add_max(scoresM[4*i + k-1], gapopenscore, MathOps::add(E, gapextendscore));
                    scoresF[4*i + k] = MathOps::add_max(scoresM[4*i + k], gapopenscore, MathOps::add(scoresF[4*i + k], gapextendscore));
                    ScoreType upTempScore = scoresM[4*i + k];
                    scoresM[4*i + k] = MathOps::add_max_relu(scoreDiag, fooArray[k], MathOps::max(E, scoresF[4*i + k]));
                    scoreDiag = upTempScore;
                    if constexpr (withEndPosition){
                        if(maxObserved < scoresM[4*i + k]){
                            maxObserved = scoresM[4*i + k];
                            positionOfMaxObserved_tileNr = tileNr;
                            positionOfMaxObserved_itemIndex = 4*i + k;
                            positionOfMaxObserved_y = row;
                        }
                    }else{
                        maxObserved = MathOps::max(maxObserved, scoresM[4*i + k]);
                    }
                }
            }

            //advance E by 1 column and F by 1 row to allow for optimized computations of remaining diagonals
            E = MathOps::add_max(scoresM[numItems-1], gapopenscore, MathOps::add(E, gapextendscore));
            for(int k = 0; k < numItems; k++){
                scoresF[k] = MathOps::add_max(scoresM[k], gapopenscore, MathOps::add(scoresF[k], gapextendscore));
            }

            //printState(row);
        };

        auto relax = [&](int row, int tileNr){
            static_assert(numItems % 4 == 0);

            using Vec4T = typename Vectorized4<ScoreType>::type;
            const Vec4T* const pssmRow4 = reinterpret_cast<const Vec4T*>(&shared_pssm.data[currentLetter][0]);
            //const Vec4T* const pssmRow4 = reinterpret_cast<const Vec4T*>(&strided_PSSM[currentLetter][0]);
            Vec4T foo = pssmRow4[0 * groupsize + group.thread_rank()];
            ScoreType fooArray[4];
            memcpy(&fooArray[0], &foo, sizeof(Vec4T));

            // E of current column and scoresF of current row are already computed

            ScoreType tempM = scoresM[0];
            scoresM[0] = MathOps::add_max_relu(scoreDiag, fooArray[0], MathOps::max(E, scoresF[0]));
            if constexpr (withEndPosition){
                if(maxObserved < scoresM[0]){
                    maxObserved = scoresM[0];
                    positionOfMaxObserved_tileNr = tileNr;
                    positionOfMaxObserved_itemIndex = 0;
                    positionOfMaxObserved_y = row;
                }
            }else{
                maxObserved = MathOps::max(maxObserved, scoresM[0]);
            }
            E = MathOps::add_max(scoresM[0],gapopenscore, MathOps::add(E, gapextendscore));
            scoresF[0] = MathOps::add_max(scoresM[0], gapopenscore, MathOps::add(scoresF[0], gapextendscore)); //this computes F of the next row !
            scoreDiag = tempM;

            #pragma unroll
            for(int i = 1; i < 4; i++){
                tempM = scoresM[i];
                scoresM[i] = MathOps::add_max_relu(scoreDiag, fooArray[i], MathOps::max(E, scoresF[i]));
                if constexpr (withEndPosition){
                    if(maxObserved < scoresM[i]){
                        maxObserved = scoresM[i];
                        positionOfMaxObserved_tileNr = tileNr;
                        positionOfMaxObserved_itemIndex = i;
                        positionOfMaxObserved_y = row;
                    }
                }else{
                    maxObserved = MathOps::max(maxObserved, scoresM[i]);
                }
                E = MathOps::add_max(scoresM[i], gapopenscore, MathOps::add(E, gapextendscore));
                scoresF[i] = MathOps::add_max(scoresM[i], gapopenscore, MathOps::add(scoresF[i], gapextendscore)); //this computes F of the next row !
                scoreDiag = tempM;
            }

            #pragma unroll
            for(int k = 1; k < numItems/4; k++){
                foo = pssmRow4[k * groupsize + group.thread_rank()];
                memcpy(&fooArray[0], &foo, sizeof(Vec4T));

                #pragma unroll
                for(int i = 0; i < 4; i++){
                    const int index = k*4+i;
                    tempM = scoresM[index];
                    scoresM[index] = MathOps::add_max_relu(scoreDiag, fooArray[i], MathOps::max(E, scoresF[index]));
                    if constexpr (withEndPosition){
                        if(maxObserved < scoresM[index]){
                            maxObserved = scoresM[index];
                            positionOfMaxObserved_tileNr = tileNr;
                            positionOfMaxObserved_itemIndex = index;
                            positionOfMaxObserved_y = row;
                        }
                    }else{
                        maxObserved = MathOps::max(maxObserved, scoresM[index]);
                    }
                    E = MathOps::add_max(scoresM[index], gapopenscore, MathOps::add(E, gapextendscore));
                    scoresF[index] = MathOps::add_max(scoresM[index], gapopenscore, MathOps::add(scoresF[index], gapextendscore)); //this computes F of the next row !
                    scoreDiag = tempM;
                }

            }

            //printState(row);
        };

        auto initScoresFirstTile = [&](){
            if(group.thread_rank() == 0){
                #pragma unroll
                for (int i=0; i < numItems; i++) {
                    scoresM[i] = 0;
                    scoresF[i] = oobscore;
                }
                scoreDiag = 0;
                scoreLeft = 0;
                E = oobscore;
            }else{
                #pragma unroll
                for (int i=0; i < numItems; i++) {
                    scoresM[i] = oobscore;
                    scoresF[i] = oobscore;
                }
                scoreDiag = oobscore;
                scoreLeft = group.thread_rank() == 1 ? 0 : oobscore;
                E = oobscore;
            }
        };

        auto shuffleScoresFirstTile = [&](){
            scoreDiag = scoreLeft;
            const ScoreType newscoreLeft = group.shfl_up(scoresM[numItems-1], 1);
            const ScoreType newE = group.shfl_up(E, 1);
            if(group.thread_rank() == 0){
                //scoreLeft is only modified in this function and is initialized with 0 for thread 0
                // assert(scoreLeft == 0);
                //scoreLeft = 0;

                // E = oobscore;
                E = gapopenscore; // After first diagonal was processed, thread 0 needs E of matrix column 1, not -infty
            }else{
                scoreLeft = newscoreLeft;
                E = newE;
            }
        };

        auto initScoresNotFirstTile = [&](int tileNr){
            if(group.thread_rank() == 0){
                #pragma unroll
                for (int i=0; i < numItems; i++) {
                    scoresM[i] = 0;
                    scoresF[i] = oobscore;
                }
                scoreDiag = 0;
                scoreLeft = leftBorderM_E.x;
                E = leftBorderM_E.y;
            }else{
                #pragma unroll
                for (int i=0; i < numItems; i++) {
                    scoresM[i] = oobscore;
                    scoresF[i] = oobscore;
                }
                scoreDiag = oobscore;
                scoreLeft = group.thread_rank() == 1 ? 0 : oobscore;
                E = oobscore;
            }
        };

        auto shuffleScoresNotFirstTile = [&](){
            scoreDiag = scoreLeft;
            const ScoreType newscoreLeft = group.shfl_up(scoresM[numItems-1], 1);
            const ScoreType newE = group.shfl_up(E, 1);
            if(group.thread_rank() == 0){
                scoreLeft = leftBorderM_E.x;
                E = leftBorderM_E.y;
            }else{
                scoreLeft = newscoreLeft;
                E = newE;
            }
        };

        //first tile
        {
            /* 
                -----------------------
                Process tile 0
                ----------------------- 
            */
            constexpr int tileNr = 0;

            //load pssm for tile 0. blockwide operation
            load_PSSM(0);

            if(alignmentId < numAlignments){

                const int8_t* subjectData = inputData.getSubject(alignmentId);
                subjectLength = inputData.getSubjectLength(alignmentId);
                groupSubjectChar4 = reinterpret_cast<const char4*>(subjectData);
                clearOutOfTileTempStorage(subjectLength);

                // if(threadIdx.x == 0){
                //     printf("subjectLength %d, queryLength %d, globalIndex %d, offset %lu\n", subjectLength, queryLength, globalIndex, charOffset);
                // }
                
                loadOffsetLimit = SDIV(subjectLength, 4);
                subjectLoadOffset = group.thread_rank();
                loadNext4Letters();
                currentLetter = badLetter;

                tempWriteOffset = group.thread_rank();

                initScoresFirstTile();

                const int numRows = (subjectLength + 1) + (groupsize-1);
                int r = 1;

                //process first groupsize - 1 diagonals which contain out-of-bound threads
                {
                    if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
                    relaxFirstDiagonal(r, tileNr); //x
                    shuffleScoresFirstTile();
                    shuffleCurrentLetter();

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
                    relax(r+1, tileNr); //y
                    shuffleScoresFirstTile();
                    shuffleCurrentLetter();

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.z; }
                    relax(r+2, tileNr); //z
                    shuffleScoresFirstTile();
                    shuffleCurrentLetter();

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.w; }
                    shuffle4Letters();

                    r = 4;
                    for(; r < groupsize - 1; r += 4){                    
                        relax(r, tileNr); //w
                        shuffleScoresFirstTile();
                        shuffleCurrentLetter(); 

                        if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
                        relax(r+1, tileNr); //x
                        shuffleScoresFirstTile();
                        shuffleCurrentLetter(); 

                        if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
                        relax(r+2, tileNr); //y
                        shuffleScoresFirstTile();
                        shuffleCurrentLetter(); 

                        if(group.thread_rank() == 0){ currentLetter = current4Letters.z; }
                        relax(r+3, tileNr); //z
                        shuffleScoresFirstTile();
                        shuffleCurrentLetter(); 

                        if(group.thread_rank() == 0){ currentLetter = current4Letters.w; }
                        shuffle4Letters();
                    }
                }
 
                //process remaining diagonals. process in chunks of 4 diagonals.
                //for those diagonals we need to store the last column of the tile to temp memory
                //last column is stored in "rightBorder"

                //r starts with r=max(4, groupsize)
                for(; r < numRows - 3; r += 4){

                    relax(r, tileNr); //w
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    shuffleScoresFirstTile();
                    shuffleCurrentLetter(); 


                    if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
                    relax(r+1, tileNr); //x
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    shuffleScoresFirstTile();
                    shuffleCurrentLetter(); 


                    if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
                    relax(r+2, tileNr); //y
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    shuffleScoresFirstTile();
                    shuffleCurrentLetter();


                    if(group.thread_rank() == 0){ currentLetter = current4Letters.z; }
                    relax(r+3, tileNr); //z 
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    shuffleScoresFirstTile();
                    shuffleCurrentLetter(); 

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.w; }

                    if((r + 4) % (4*group.size()) == 0){
                        //used up all query letters stored across the group. reload
                        loadNext4Letters();
                    }else{
                        //get next 4 letters from neighbor
                        shuffle4Letters();
                    }

                    if((r + 4) % (group.size()) == 0){
                        #ifdef PRINT_WRITE
                        printf("tid %d, write %f %f to %d\n", group.thread_rank(), tileLastColumnM_E.x, tileLastColumnM_E.y, tempWriteOffset);
                        #endif
                        groupTempStorage[tempWriteOffset] = tileLastColumnM_E;
                        tempWriteOffset += group.size();
                    }                    
                }

                //can have at most 3 remaining rows
                if(r < numRows){
                    relax(r, tileNr); //w
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    shuffleScoresFirstTile();
                    shuffleCurrentLetter();

                }
                if(r+1 < numRows){
                    if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
                    relax(r+1, tileNr); //x
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores                    
                    shuffleScoresFirstTile();
                    shuffleCurrentLetter();
                }
                if(r+2 < numRows){
                    if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
                    relax(r+2, tileNr); //y
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                }

                const int totalChunksOfFour = subjectLength / 4;
                const int unsavedChunksOfFour = totalChunksOfFour % (group.size() / 4);
                const int numThreadsWithValidTileLastColumn = unsavedChunksOfFour * 4 + subjectLength % 4;
                if(numThreadsWithValidTileLastColumn > 0){
                    const int firstValidThread = group.size() - numThreadsWithValidTileLastColumn;
                    if(group.thread_rank() >= firstValidThread){
                        #ifdef PRINT_WRITE
                        printf("last write. tid %d, write %f %f to %d\n", group.thread_rank(), tileLastColumnM_E.x, tileLastColumnM_E.y, tempWriteOffset - firstValidThread);
                        #endif
                        groupTempStorage[tempWriteOffset - firstValidThread] = tileLastColumnM_E;
                    }
                }
            }
        }



        for(int tileNr = 1; tileNr < numTiles; tileNr++){
            load_PSSM(tileNr);

            /* 
                -----------------------
                Process tile tileNr
                ----------------------- 
            */

            if(alignmentId < numAlignments){

                subjectLoadOffset = group.thread_rank();
                loadNext4Letters();
                currentLetter = badLetter;

                tempWriteOffset = group.thread_rank();

                #ifdef PRINT_LOAD
                printf("tid %d, load %f %f from %d\n", group.thread_rank(), groupTempStorage[group.thread_rank()].x, groupTempStorage[group.thread_rank()].y, group.thread_rank());
                #endif
                leftBorderM_E = groupTempStorage[group.thread_rank()];
                tempLoadOffset = group.size() + group.thread_rank();


                initScoresNotFirstTile(tileNr);

                const int numRows = (subjectLength + 1) + (groupsize-1);
                int r = 1;

                //process first groupsize - 1 diagonals which contain out-of-bound threads
                {
                    if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
                    relaxFirstDiagonal(r, tileNr); //x
                    shuffleLeftBorder(); //must be called before shuffleScores
                    shuffleScoresNotFirstTile();
                    shuffleCurrentLetter();

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
                    relax(r+1, tileNr); //y
                    shuffleLeftBorder(); //must be called before shuffleScores
                    shuffleScoresNotFirstTile();
                    shuffleCurrentLetter();

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.z; }
                    relax(r+2, tileNr); //z
                    shuffleLeftBorder(); //must be called before shuffleScores
                    shuffleScoresNotFirstTile();
                    shuffleCurrentLetter();

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.w; }
                    shuffle4Letters();

                    r = 4;
                    for(; r < groupsize - 1; r += 4){                    
                        relax(r, tileNr); //w
                        shuffleLeftBorder(); //must be called before shuffleScores
                        shuffleScoresNotFirstTile();
                        shuffleCurrentLetter(); 

                        if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
                        relax(r+1, tileNr); //x
                        shuffleLeftBorder(); //must be called before shuffleScores
                        shuffleScoresNotFirstTile();
                        shuffleCurrentLetter(); 

                        if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
                        relax(r+2, tileNr); //y
                        shuffleLeftBorder(); //must be called before shuffleScores
                        shuffleScoresNotFirstTile();
                        shuffleCurrentLetter(); 

                        if(group.thread_rank() == 0){ currentLetter = current4Letters.z; }
                        relax(r+3, tileNr); //z
                        shuffleLeftBorder(); //must be called before shuffleScores
                        shuffleScoresNotFirstTile();
                        shuffleCurrentLetter(); 

                        if(group.thread_rank() == 0){ currentLetter = current4Letters.w; }
                        shuffle4Letters();
                    }
                }

                //process remaining diagonals. process in chunks of 4 diagonals.
                //for those diagonals we need to store the last column of the tile to temp memory
                //last column is stored in "rightBorder"

                //r starts with r=max(4, groupsize)
                for(; r < numRows - 3; r += 4){

                    relax(r, tileNr); //w
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    if(r % group.size() == 0 && r < subjectLength){
                        #ifdef PRINT_LOAD
                        printf("tid %d, load %f %f from %d\n", group.thread_rank(), groupTempStorage[tempLoadOffset].x, groupTempStorage[tempLoadOffset].y, tempLoadOffset);
                        #endif
                        leftBorderM_E = groupTempStorage[tempLoadOffset];
                        tempLoadOffset += group.size();
                    }else{
                        shuffleLeftBorder(); //must be called before shuffleScores
                    }
                    shuffleScoresNotFirstTile();
                    shuffleCurrentLetter(); 
                    



                    if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
                    relax(r+1, tileNr); //x
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    shuffleLeftBorder(); //must be called before shuffleScores
                    shuffleScoresNotFirstTile();
                    shuffleCurrentLetter(); 
                    

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
                    relax(r+2, tileNr); //y
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    shuffleLeftBorder(); //must be called before shuffleScores
                    shuffleScoresNotFirstTile();
                    shuffleCurrentLetter(); 
                    

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.z; }
                    relax(r+3, tileNr); //z
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    shuffleLeftBorder(); //must be called before shuffleScores
                    shuffleScoresNotFirstTile();
                    shuffleCurrentLetter(); 

                    if(group.thread_rank() == 0){ currentLetter = current4Letters.w; }

                    if((r + 4) % (4*group.size()) == 0){
                        //used up all query letters stored across the group. reload
                        loadNext4Letters();
                    }else{
                        //get next 4 letters from neighbor
                        shuffle4Letters();
                    }

                    if((r + 4) % (group.size()) == 0){
                        #ifdef PRINT_WRITE
                        printf("tid %d, write %f %f to %d\n", group.thread_rank(), tileLastColumnM_E.x, tileLastColumnM_E.y, tempWriteOffset);
                        #endif
                        groupTempStorage[tempWriteOffset] = tileLastColumnM_E;
                        tempWriteOffset += group.size();
                    }
                }

                //can have at most 3 remaining rows
                if(r < numRows){
                    relax(r, tileNr); //w
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    if(r % group.size() == 0 && r < subjectLength){
                        #ifdef PRINT_LOAD
                        printf("last load. tid %d, load %f %f from %d\n", group.thread_rank(), groupTempStorage[tempLoadOffset].x, groupTempStorage[tempLoadOffset].y, tempLoadOffset);
                        #endif
                        leftBorderM_E = groupTempStorage[tempLoadOffset];
                        tempLoadOffset += group.size();
                    }else{
                        shuffleLeftBorder(); //must be called before shuffleScores
                    }
                    shuffleScoresNotFirstTile();
                    shuffleCurrentLetter();
                    

                }
                if(r+1 < numRows){
                    if(group.thread_rank() == 0){ currentLetter = current4Letters.x; }
                    relax(r+1, tileNr); //x
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                    shuffleLeftBorder(); //must be called before shuffleScores
                    shuffleScoresNotFirstTile();
                    shuffleCurrentLetter();
                    
                }
                if(r+2 < numRows){
                    if(group.thread_rank() == 0){ currentLetter = current4Letters.y; }
                    relax(r+2, tileNr); //y
                    shuffleTileLastColumn(); //must be called before setTileLastColumn
                    setTileLastColumn(); //must be called before shuffleScores
                }

                const int totalChunksOfFour = subjectLength / 4;
                const int unsavedChunksOfFour = totalChunksOfFour % (group.size() / 4);
                const int numThreadsWithValidTileLastColumn = unsavedChunksOfFour * 4 + subjectLength % 4;
                if(numThreadsWithValidTileLastColumn > 0){
                    const int firstValidThread = group.size() - numThreadsWithValidTileLastColumn;
                    if(group.thread_rank() >= firstValidThread){
                        #ifdef PRINT_WRITE
                        printf("last write. tid %d, write %f %f\n", group.thread_rank(), tileLastColumnM_E.x, tileLastColumnM_E.y);
                        #endif
                        groupTempStorage[tempWriteOffset - firstValidThread] = tileLastColumnM_E;
                    }
                }
            }
        }
        //printState(r+3);

        if constexpr (withEndPosition){
 
            if(alignmentId < numAlignments){
                const int myQueryEndInclusive = positionOfMaxObserved_tileNr * groupsize * numItems 
                    + group.thread_rank() * numItems + positionOfMaxObserved_itemIndex;
                const int mySubjectEndInclusive = positionOfMaxObserved_y - group.thread_rank() - 1;

                //disable any maxima that come from out-of-bounds cells
                if(myQueryEndInclusive >= queryLength){
                    maxObserved = -1;
                }
                if(mySubjectEndInclusive < 0 || mySubjectEndInclusive >= subjectLength){
                    maxObserved = -1;
                }
                const int3 packed = make_int3(maxObserved, 
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
            }
        }else{
            if(alignmentId < numAlignments){
                const int alignmentScore = compatibility::group_reduce_max(group, maxObserved);

                if(group.thread_rank() == 0){
                    devAlignmentScores[alignmentId] = alignmentScore;
                }
            }
        }
    }

}






LIBMARV_NAMESPACE_WITH_NESTING_END



#endif