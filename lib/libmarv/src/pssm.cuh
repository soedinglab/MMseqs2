#ifndef PSSM_CUH
#define PSSM_CUH

#include "cuda_hip_rename.h"
#include "stream_ordered_allocation_helper.cuh"

#include "config.hpp"
#include "types.hpp"
#include "convert.cuh"
#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/simple_allocation.cuh"
#include "custom_score_types.cuh"
#include "util.cuh"
#include "cuda_errorcheck.cuh"

#if !defined(__HIPCC__)
#include <cuda_fp16.h>
#include <cuda_fp8.h>
#else
#include <hip/hip_fp16.h>
#include <hip/hip_fp8.h>
#endif

#include <optional>
#include <vector>
#include <cassert>
#include <cuda/std/mdspan>

#include "namespace.hpp"
LIBMARV_NAMESPACE_WITH_NESTING_BEGIN

template<class T>
struct PSSM_2D_View{
    int numRows;
    int numColumns;
    int stride;
    const T* data;

    __host__ __device__
    const T* operator[](int encodedSubjectLetter) const{
        return data + encodedSubjectLetter * stride;
    }
};

template<class T>
struct PSSM_2D_ModifiableView{
    int numRows;
    int numColumns;
    int stride;
    T* data;

    __host__ __device__
    const T* operator[](int encodedSubjectLetter) const{
        return data + encodedSubjectLetter * stride;
    }

    __host__ __device__
    T* operator[](int encodedSubjectLetter){
        return data + encodedSubjectLetter * stride;
    }
};

struct PSSM{
    int alphabetSize;
    SequenceLengthT queryLength;
    std::vector<int> data;

    PSSM(int queryLength_, int alphabetSize_) : 
        alphabetSize(alphabetSize_),
        queryLength(queryLength_),
        data(alphabetSize * queryLength){

    }

    int* operator[](int encodedSubjectLetter){
        return data.data() + encodedSubjectLetter * queryLength;
    }
    const int* operator[](int encodedSubjectLetter) const{
        return data.data() + encodedSubjectLetter * queryLength;
    }

    int getMinElement() const{
        if(data.size() == 0){
            return 0;
        }else{
            return *std::min_element(data.begin(), data.end());
        }
    }

    PSSM_2D_View<int> makeView(int startQueryPos = 0) const{
        assert(startQueryPos < queryLength);
        PSSM_2D_View<int> result;
        result.numRows = alphabetSize;
        result.numColumns = queryLength - startQueryPos;
        result.stride = queryLength;
        result.data = data.data() + startQueryPos;
        return result;
    }

    //Generator::operator()(int queryPosition, int alphabetIndex) should return the pssm score
    //for specific query position and alphabet letter
    template<class Generator>
    static PSSM fromGenerator(int alphabetSize, int queryLength, Generator generator){
        PSSM pssm(queryLength, alphabetSize);

        for (int subjectLetter = 0; subjectLetter < alphabetSize; subjectLetter++) {
            for (int col = 0; col < queryLength; col++){
                pssm[subjectLetter][col] = generator(col, subjectLetter);
            }
        }

        return pssm;
    }

    static PSSM fromPSSM(const char* /*encodedQuery*/, int queryLength, const int8_t * pssm, int alphabetSize){
        PSSM retPssm(queryLength, alphabetSize);
        for (int subjectLetter = 0; subjectLetter < alphabetSize; subjectLetter++) {
            for (int col = 0; col < queryLength; col++){
                retPssm[subjectLetter][col] = static_cast<int>(pssm[subjectLetter * queryLength + col]);
            }
        }
        return retPssm;
    }

    // query must have been encoded with ConvertAA_20
    static PSSM fromBlosum(BlosumType blosumType, const char* encodedQuery, int queryLength){
        auto make = [](const auto& blosum2D, const char* encodedQuery, int queryLength){
            const int alphabetSize = blosum2D.size();


            // auto generator = [&](int queryPosition, int alphabetIndex){
            //     const int queryLetter = encodedQuery[queryPosition];
            //     return blosum2D[queryLetter][alphabetIndex];
            // };

            auto generator_mmseqs_conversion = [&](int queryPosition, int alphabetIndex){
                //blosum layout is for ncbi encoded letters, but input is mmseqs encoded.
                //convert both queryletter and alphabetIndex to ncbi format
                const int queryLetter_ncbi = ConvertAA_20_mmseqs_to_ncbi{}(encodedQuery[queryPosition]);
                const int alphabetIndex_ncbi = ConvertAA_20_mmseqs_to_ncbi{}(alphabetIndex);
                return blosum2D[queryLetter_ncbi][alphabetIndex_ncbi];
            };
            return fromGenerator(alphabetSize, queryLength, generator_mmseqs_conversion);
        };

        switch(blosumType){
            case BlosumType::BLOSUM45_20: {
                BLOSUM45_20 blosum;
                return make(blosum.get2D(), encodedQuery, queryLength);
            }
            case BlosumType::BLOSUM50_20: {
                BLOSUM50_20 blosum;
                return make(blosum.get2D(), encodedQuery, queryLength);
            }
            case BlosumType::BLOSUM62_20: {
                BLOSUM62_20 blosum;
                return make(blosum.get2D(), encodedQuery, queryLength);
            }
            case BlosumType::BLOSUM80_20: {
                BLOSUM80_20 blosum;
                return make(blosum.get2D(), encodedQuery, queryLength);
            }
            default:
                throw std::runtime_error("PSSM::fromBlosum invalid blosum type");
        }

    }

};


struct GpuPSSM{
    int alphabetSize;
    SequenceLengthT queryLength;
    helpers::SimpleAllocationDevice<int, 0> data;

    GpuPSSM() = default;

    GpuPSSM(int queryLength_, int alphabetSize_) : 
        alphabetSize(alphabetSize_),
        queryLength(queryLength_),
        data(alphabetSize * queryLength){

    }

    GpuPSSM(const PSSM& rhs, cudaStream_t stream){
        upload(rhs, stream);
    }

    void resize(int queryLength_, int alphabetSize_){
        alphabetSize = alphabetSize_;
        queryLength = queryLength_;
        data.resize(alphabetSize * queryLength);
    }

    void upload(const PSSM& rhs, cudaStream_t stream){
        alphabetSize = rhs.alphabetSize;
        queryLength = rhs.queryLength;
        data.resize(rhs.data.size());
        cudaMemcpyAsync(
            data.data(), 
            rhs.data.data(), 
            sizeof(int) * rhs.data.size(), 
            cudaMemcpyHostToDevice, 
            stream
        ); CUERR
    }

    PSSM_2D_View<int> makeView(int startQueryPos = 0) const{
        assert(startQueryPos < queryLength);
        PSSM_2D_View<int> result;
        result.numRows = alphabetSize;
        result.numColumns = queryLength - startQueryPos;
        result.stride = queryLength;
        result.data = data.data() + startQueryPos;
        return result;
    }
};


//PSSM alignment kernel will use vectorized loads to load multiple half elements
// accessSizeBytes specifies the vector size in bytes, e.g. 16 for float4
template<int accessSizeBytes, class InputT, class OutputT>
__global__
void permute_PSSM_for_gapless_kernel(
    PSSM_2D_ModifiableView<OutputT> resultView,
    PSSM_2D_View<InputT> inputView,
    const int numRegs,
    const int groupsize
) {
    static_assert(sizeof(OutputT) == 4 || sizeof(OutputT) == 2 || sizeof(OutputT) == 1);
    static_assert(accessSizeBytes == 4 || accessSizeBytes == 8 || accessSizeBytes == 16); //float, float2, float4

    //2 for half / short, 4 for int8
    constexpr int numOutputElemsPer4Bytes = 4 / sizeof(OutputT);

    constexpr int numFloatsPerAccess = accessSizeBytes / 4;
    constexpr int numOutputElemsPerAccess = numFloatsPerAccess * numOutputElemsPer4Bytes;

    const int thid = threadIdx.x + blockIdx.x*blockDim.x;

    const int numColumns = inputView.numColumns;
    const int numRows = inputView.numRows;

    const int tileSize = (numRegs * groupsize * numOutputElemsPer4Bytes);

    for(int inputRow = blockIdx.y; inputRow < numRows; inputRow += gridDim.y){

        for (int inputCol = thid; inputCol < numColumns; inputCol += blockDim.x*gridDim.x) {
            const int tileId = inputCol / tileSize;
            const int tileColumnOffset = tileId * tileSize;
            const int columnInTile = inputCol - tileColumnOffset;

            const int stride = tileSize/numOutputElemsPer4Bytes;
            const int offset = (numOutputElemsPer4Bytes*(columnInTile%stride) + columnInTile/stride)%numOutputElemsPerAccess;
            const int thread = (columnInTile%stride)/numRegs;
            const int part = (columnInTile%numRegs)/numFloatsPerAccess;
            const int resultCol = numOutputElemsPerAccess*thread + offset + part*numOutputElemsPerAccess*groupsize;

            resultView[inputRow][tileColumnOffset + resultCol] = inputView[inputRow][inputCol];
        }
    }
}



template<class InputT, class OutputT>
__global__
void permute_PSSM_for_gapless_kernel_withRelaxChunkSize(
    PSSM_2D_ModifiableView<OutputT> resultView,
    PSSM_2D_View<InputT> inputView,
    const int numRegs,
    const int groupsize,
    int relaxChunkSize, // number of regs per thread for which data is loaded per load instruction
    int workVectorLanesPerReg, //number of packed elements per reg for gapless computation
    int loadVectorLanesPerReg //number of packed elements per reg for data loading during computation
) {
    // if(threadIdx.x + blockIdx.x + threadIdx.y + blockIdx.y == 0){
    //     printf("numRegs %d, groupsize %d, relaxChunkSize %d, workVectorLanesPerReg %d, loadVectorLanesPerReg %d\n", 
    //         numRegs,
    //         groupsize,
    //         relaxChunkSize,
    //         workVectorLanesPerReg,
    //         loadVectorLanesPerReg
    //     );
    // }
    const int workTileSize = (groupsize * numRegs * workVectorLanesPerReg);
    const int numRelaxChunks = SDIV(numRegs, relaxChunkSize);
    const int loadTileSize = (groupsize * numRelaxChunks * relaxChunkSize * loadVectorLanesPerReg);

    const int thid = threadIdx.x + blockIdx.x*blockDim.x;

    const int numColumns = inputView.numColumns;
    const int numRows = inputView.numRows;

    const int strideForOutputReg = numRegs * groupsize;

    for(int inputRow = blockIdx.y; inputRow < numRows; inputRow += gridDim.y){

        for (int inputCol = thid; inputCol < numColumns; inputCol += blockDim.x*gridDim.x) {
            //current tile in input pssm
            const int tileId = inputCol / workTileSize;
            //global column number of first column in current input tile
            const int tileColumnOffset = tileId * workTileSize;
            //column number local to the current input tile
            const int columnInInputTile = inputCol - tileColumnOffset;
            //global column number of first column of current tile in output tile
            const int outputTileColumnOffset = tileId * loadTileSize;

            const int outputReg = columnInInputTile % strideForOutputReg;
            const int indexInOutputReg = columnInInputTile / strideForOutputReg;
            const int columnInUnpermutedOutputTile = workVectorLanesPerReg*outputReg + indexInOutputReg;


            //permute output registers for shared memory access

            const int threadId = columnInUnpermutedOutputTile / (numRegs * workVectorLanesPerReg);
            const int columnInThread = columnInUnpermutedOutputTile % (numRegs * workVectorLanesPerReg);
            const int relaxChunk = columnInThread / (relaxChunkSize * loadVectorLanesPerReg);
            const int columnInRelaxChunk = columnInThread % (relaxChunkSize * loadVectorLanesPerReg);

            const int outputAccessChunk = relaxChunk * groupsize + threadId;
            const int outputCol = outputAccessChunk * relaxChunkSize * loadVectorLanesPerReg + columnInRelaxChunk;
            resultView[inputRow][outputTileColumnOffset + outputCol] = convert<OutputT>(inputView[inputRow][inputCol]);
            // if(inputRow == 0){
            //     printf("%d <- %d\n", outputCol, inputCol);
            // }
        }
    }
}





template<int accessSizeBytes, class InputT, class OutputT>
__global__
void permute_PSSM_for_SW_kernel(
    PSSM_2D_ModifiableView<OutputT> resultView,
    PSSM_2D_View<InputT> inputView,
    int elementsPerThread,
    int groupsize
) {
    static_assert(accessSizeBytes == 4 || accessSizeBytes == 8 || accessSizeBytes == 16); //float, float2, float4
    static_assert(accessSizeBytes % sizeof(OutputT) == 0);
    constexpr int numElementsPerAccess = accessSizeBytes / sizeof(OutputT);
    assert(elementsPerThread % numElementsPerAccess == 0);

    const int tileSize = (elementsPerThread * groupsize);
    const int numAccesses = elementsPerThread / numElementsPerAccess;

    const int numColumns = inputView.numColumns;
    const int numRows = inputView.numRows;

    for(int inputRow = blockIdx.y; inputRow < numRows; inputRow += gridDim.y){

        for (int inputCol = threadIdx.x + blockIdx.x * blockDim.x; inputCol < numColumns; inputCol += blockDim.x * gridDim.x) {
            const int tileId = inputCol / tileSize;
            const int tileColumnOffset = tileId * tileSize;
            const int columnInTile = inputCol - tileColumnOffset;

            const int accessChunk = columnInTile / numElementsPerAccess;
            const int elementIdInAccessChunk = columnInTile % numElementsPerAccess;
            const int accessChunkIdInThread = accessChunk % numAccesses;
            const int threadId = accessChunk / numAccesses;

            const int outputAccessChunk = accessChunkIdInThread * groupsize + threadId;
            const int outputCol = outputAccessChunk * numElementsPerAccess + elementIdInAccessChunk;
            resultView[inputRow][tileColumnOffset + outputCol] = inputView[inputRow][inputCol];
        }
    }
}

struct GpuPermutedPSSMforGapless{
    int alphabetSize;
    int numRegs;
    int groupsize;
    int columnstride;
    int allocatedNumRegs;
    SequenceLengthT queryLength;
    helpers::SimpleAllocationDevice<char, 0> data;

    std::optional<int> substitutionScoreBias_optional;

    template<class ScoreType, class PssmScoreType>
    void resize(int groupsize_, int numRegs_, int alphabetSize_, int queryLength_, int allocatedNumRegs_, cudaStream_t stream){
        assert(512 % sizeof(PssmScoreType) == 0);

        if((groupsize_ * numRegs_) % 2 == 1){
            throw std::runtime_error("GpuPermutedPSSMforGapless resize error. elements per row must be even");
        }
        groupsize = groupsize_;
        numRegs = numRegs_;
        alphabetSize = alphabetSize_;
        queryLength = queryLength_;
        allocatedNumRegs = allocatedNumRegs_;

        int lanes = 1;
        if(IsVec2Type<ScoreType>::value){
            lanes = 2;
        }else if(IsVec4Type<ScoreType>::value){
            lanes = 4;
        }

        const int workTileSize = (numRegs * groupsize * lanes);
        const int loadTileSize = (allocatedNumRegs * groupsize * lanes);
        const int numTiles = SDIV(queryLength, workTileSize);

        const int numAllocatedColumns = numTiles * loadTileSize;

        columnstride = numAllocatedColumns;
        // std::cout << "allocatedNumRegs " << allocatedNumRegs << ", workTileSize " << workTileSize << ", loadTileSize " << loadTileSize 
        //     << ", numTiles " << numTiles 
        //     << ", columnstride " << columnstride 
        //     << ", queryLength " << queryLength
        //     << "\n";
        //numPaddedColumns = (SDIV(groupsize * numItems * sizeof(PssmScoreType), 512) * 512) / sizeof(PssmScoreType);

        data.resize(sizeof(PssmScoreType) * alphabetSize * numAllocatedColumns);
        //init with 0 so oob elements won't contribute to the score
        cudaMemsetAsync(data.data(), 0, sizeof(char) * data.size(), stream);
    }

    // template<class T>
    // void resize(int groupsize_, int numRegs_, int alphabetSize_, int queryLength_, cudaStream_t stream){
    //     resize<T>(groupsize_, numRegs_, alphabetSize_, queryLength_, numRegs_, stream);
    // }

    template<class T>
    PSSM_2D_View<T> makeScalarView() const{
        PSSM_2D_View<T> view;
        view.numRows = alphabetSize;
        view.numColumns = columnstride;
        view.stride = columnstride;
        view.data = reinterpret_cast<const T*>(data.data());

        return view;
    }

    template<class T>
    PSSM_2D_View<T> makeVec2View() const{
        assert(columnstride % 2 == 0);

        PSSM_2D_View<T> view;
        view.numRows = alphabetSize;
        view.numColumns = columnstride / 2;
        view.stride = columnstride / 2;
        view.data = reinterpret_cast<const T*>(data.data());

        return view;
    }

    template<class T>
    PSSM_2D_View<T> makeVec4View() const{
        assert(columnstride % 4 == 0);

        PSSM_2D_View<T> view;
        view.numRows = alphabetSize;
        view.numColumns = columnstride / 4;
        view.stride = columnstride / 4;
        view.data = reinterpret_cast<const T*>(data.data());

        return view;
    }


    PSSM_2D_View<half2> makeHalf2View() const{
        return makeVec2View<half2>();
    }

    PSSM_2D_View<short2> makeShort2View() const{
        return makeVec2View<short2>();
    }

    PSSM_2D_View<ScoreType_u8x4> makeUint8x4View() const{
        return makeVec4View<ScoreType_u8x4>();
    }

    // template<class OutputT, int accessSizeBytes, class InputT>
    // void fromGpuPSSMView(PSSM_2D_View<InputT> inputView, int groupsize_, int numRegs_, std::optional<int> substitutionScoreBias, cudaStream_t stream){
    //     resize<OutputT>(groupsize_, numRegs_, inputView.numRows, inputView.numColumns, stream);

    //     substitutionScoreBias_optional = substitutionScoreBias;

    //     PSSM_2D_ModifiableView<OutputT> resultView;
    //     resultView.numRows = alphabetSize;
    //     resultView.numColumns = queryLength;
    //     resultView.stride = columnstride;
    //     resultView.data = reinterpret_cast<OutputT*>(data.data());

    //     dim3 block(128,1,1);
    //     dim3 grid(SDIV(inputView.numColumns, block.x), inputView.numRows, 1);

    //     permute_PSSM_for_gapless_kernel<accessSizeBytes><<<grid, block, 0, stream>>>(
    //         resultView,
    //         inputView,
    //         numRegs,
    //         groupsize
    //     ); CUERR;
    // }

    template<class ScoreType, class PssmScoreType, class InputT>
    void fromGpuPSSMView_withRelaxChunkSize(PSSM_2D_View<InputT> inputView, int groupsize_, int numRegs_, std::optional<int> substitutionScoreBias, cudaStream_t stream){
        
        int workVectorLanesPerReg = 1;
        if(IsVec2Type<ScoreType>::value){
            workVectorLanesPerReg = 2;
        }else if(IsVec4Type<ScoreType>::value){
            workVectorLanesPerReg = 4;
        }

        //score half2, pssm half2 -> 2
        //score half2, pssm fp8x2 -> 2
        //score int8x4, pssm int8x4 -> 4
        int loadVectorLanesPerReg = 1;
        if(IsVec2Type<ScoreType>::value){
            loadVectorLanesPerReg = 2;
        }else if(IsVec4Type<ScoreType>::value){
            loadVectorLanesPerReg = 4;
        }

        using ScalarPssmScoreType = typename ScalarScoreType<PssmScoreType>::type;

        constexpr int accessSizeBytes = 16;
        int relaxChunkSize = accessSizeBytes / sizeof(PssmScoreType);
        const int numRelaxChunks = SDIV(numRegs_, relaxChunkSize);

        //resize<OutputT>(groupsize_, numRegs_, inputView.numRows, inputView.numColumns, stream);
        resize<ScoreType, ScalarPssmScoreType>(groupsize_, numRegs_, inputView.numRows, inputView.numColumns, relaxChunkSize * numRelaxChunks, stream);

        substitutionScoreBias_optional = substitutionScoreBias;

        PSSM_2D_ModifiableView<ScalarPssmScoreType> resultView;
        resultView.numRows = alphabetSize;
        resultView.numColumns = queryLength;
        resultView.stride = columnstride;
        resultView.data = reinterpret_cast<ScalarPssmScoreType*>(data.data());

        dim3 block(128,1,1);
        dim3 grid(SDIV(inputView.numColumns, block.x), inputView.numRows, 1);

#if !defined(__METAL_BACKEND__) && !defined(__SIMD_BACKEND__)
        permute_PSSM_for_gapless_kernel_withRelaxChunkSize<<<grid, block, 0, stream>>>(
            resultView,
            inputView,
            numRegs,
            groupsize,
            relaxChunkSize, // number of regs per thread for which data is loaded per load instruction
            workVectorLanesPerReg, //number of packed elements per reg for gapless computation
            loadVectorLanesPerReg //number of packed elements per reg for data loading during computation
        ); CUERR;
#endif
    }
    

};

struct GpuPermutedPSSMforSW{
    int alphabetSize;
    int numRegs;
    int groupsize;
    int columnstride;
    SequenceLengthT queryLength;
    helpers::SimpleAllocationDevice<char, 0> data;

    std::optional<int> substitutionScoreBias_optional;

    template<class T>
    void resize(int groupsize_, int numRegs_, int alphabetSize_, int queryLength_, cudaStream_t stream){
        assert(512 % sizeof(T) == 0);

        groupsize = groupsize_;
        numRegs = numRegs_;
        alphabetSize = alphabetSize_;
        queryLength = queryLength_;

        const int tileSize = (numRegs * groupsize);
        const int numTiles = SDIV(queryLength, tileSize);
        columnstride = numTiles * tileSize;
        //numPaddedColumns = (SDIV(groupsize * numItems * sizeof(PssmScoreType), 512) * 512) / sizeof(PssmScoreType);

        data.resize(sizeof(T) * alphabetSize * tileSize * numTiles);
        //init with 0 so oob elements won't contribute to the score
        cudaMemsetAsync(data.data(), 0, sizeof(char) * data.size(), stream);
    }

    template<class T>
    PSSM_2D_View<T> makeView() const{
        PSSM_2D_View<T> view;
        view.numRows = alphabetSize;
        view.numColumns = columnstride;
        view.stride = columnstride;
        view.data = reinterpret_cast<const T*>(data.data());

        return view;
    }

    template<int accessSizeBytes, class OutputT, class InputT>
    void fromGpuPSSMView(PSSM_2D_View<InputT> inputView, int groupsize_, int numRegs_, cudaStream_t stream){
        resize<OutputT>(groupsize_, numRegs_, inputView.numRows, inputView.numColumns, stream);

        PSSM_2D_ModifiableView<OutputT> resultView;
        resultView.numRows = alphabetSize;
        resultView.numColumns = queryLength;
        resultView.stride = columnstride;
        resultView.data = reinterpret_cast<OutputT*>(data.data());

        dim3 block(128,1,1);
        dim3 grid(SDIV(inputView.numColumns, block.x), inputView.numRows, 1);

#if !defined(__METAL_BACKEND__) && !defined(__SIMD_BACKEND__)
        permute_PSSM_for_SW_kernel<accessSizeBytes><<<grid, block, 0, stream>>>(
            resultView,
            inputView,
            numRegs,
            groupsize
        ); CUERR;
#endif
    }
};


struct GpuPermutedPSSM{

    int queryLength;

    GpuPermutedPSSM() = default;
    GpuPermutedPSSM(int queryLength_, GpuPermutedPSSMforGapless* gapless_, GpuPermutedPSSMforSW* sw_) : queryLength(queryLength_), gapless_ptr(gapless_), sw_ptr(sw_){}

    GpuPermutedPSSMforGapless& gapless(){
        return *gapless_ptr;
    } 
    const GpuPermutedPSSMforGapless& gapless() const{
        return *gapless_ptr;
    }

    GpuPermutedPSSMforSW& sw(){
        return *sw_ptr;
    } 
    const GpuPermutedPSSMforSW& sw() const{
        return *sw_ptr;
    } 

private:    
    GpuPermutedPSSMforGapless* gapless_ptr;
    GpuPermutedPSSMforSW* sw_ptr;
};



template<class ScoreType, class PssmScoreType, class TT, int groupsize, int numRegs>
void makeConvertedGpuPssm_gapless(
    GpuPermutedPSSMforGapless& result,
    const cuda::std::mdspan<TT, cuda::std::dextents<int,2>> hostPssmView, //numRows * queryLength elements, score to compare encoded subject letter s and query position q < queryLength must be at pssm[s][q]
    cudaStream_t stream
){
    using T = typename std::remove_cv<TT>::type;
    const int numRows = hostPssmView.extent(0);
    const int queryLength = hostPssmView.extent(1);
 
    int workVectorLanesPerReg = 1;
    if(IsVec2Type<ScoreType>::value){
        workVectorLanesPerReg = 2;
    }else if(IsVec4Type<ScoreType>::value){
        workVectorLanesPerReg = 4;
    }

    //score half2, pssm half2 -> 2
    //score half2, pssm fp8x2 -> 2
    //score int8x4, pssm int8x4 -> 4
    int loadVectorLanesPerReg = 1;
    if(IsVec2Type<ScoreType>::value){
        loadVectorLanesPerReg = 2;
    }else if(IsVec4Type<ScoreType>::value){
        loadVectorLanesPerReg = 4;
    }

    using ScalarPssmScoreType = typename ScalarScoreType<PssmScoreType>::type;


    std::vector<T> h_pssm(numRows * queryLength);
    int substitutionScoreBias = 0;
    if constexpr(IsVec4Type<ScoreType>::value){
        for(int s = 0; s < numRows; s++){
            for(int p = 0; p < queryLength; p++){
                int val = hostPssmView(s,p);
                substitutionScoreBias = std::min(substitutionScoreBias, val);
            }
        }
        substitutionScoreBias = std::abs(substitutionScoreBias);
        result.substitutionScoreBias_optional = substitutionScoreBias;
    }
    for(int s = 0; s < numRows; s++){
        for(int p = 0; p < queryLength; p++){
            h_pssm[s * queryLength + p] = hostPssmView(s,p) + substitutionScoreBias;
        }
    }
    T* d_pssm;
    CUDACHECK(helper_gpuMallocAsync((void**)&d_pssm, sizeof(T) * numRows * queryLength, stream));
    CUDACHECK(cudaMemcpyAsync(d_pssm, h_pssm.data(), sizeof(T) * numRows * queryLength, cudaMemcpyHostToDevice, stream));


    constexpr int accessSizeBytes = 16;
    int relaxChunkSize = accessSizeBytes / sizeof(PssmScoreType);
    const int numRelaxChunks = SDIV(numRegs, relaxChunkSize);

    result.resize<ScoreType, ScalarPssmScoreType>(groupsize, numRegs, numRows, queryLength, relaxChunkSize * numRelaxChunks, stream);


    PSSM_2D_View<T> unpermutedView;
    unpermutedView.numRows = numRows;
    unpermutedView.numColumns = queryLength;
    unpermutedView.stride = queryLength;
    unpermutedView.data = reinterpret_cast<T*>(d_pssm);

    PSSM_2D_ModifiableView<ScalarPssmScoreType> resultView;
    resultView.numRows = numRows;
    resultView.numColumns = queryLength;
    resultView.stride = result.columnstride;
    resultView.data = reinterpret_cast<ScalarPssmScoreType*>(result.data.data());

    dim3 block(128,1,1);
    dim3 grid(SDIV(unpermutedView.numColumns, block.x), unpermutedView.numRows, 1);

#if !defined(__METAL_BACKEND__) && !defined(__SIMD_BACKEND__)
    permute_PSSM_for_gapless_kernel_withRelaxChunkSize<<<grid, block, 0, stream>>>(
        resultView,
        unpermutedView,
        numRegs,
        groupsize,
        relaxChunkSize, // number of regs per thread for which data is loaded per load instruction
        workVectorLanesPerReg, //number of packed elements per reg for gapless computation
        loadVectorLanesPerReg //number of packed elements per reg for data loading during computation
    );
#endif
    CUDACHECKASYNC;

    CUDACHECK(helper_gpuFreeAsync((void*)d_pssm, stream));
}

template<class ScoreType, class PssmScoreType, class TT, int groupsize, int numRegs>
void makeConvertedGpuPssm(
    GpuPermutedPSSMforGapless& result,
    const cuda::std::mdspan<TT, cuda::std::dextents<int,2>> hostPssmView, //numRows * queryLength elements, score to compare encoded subject letter s and query position q < queryLength must be at pssm[s][q]
    cudaStream_t stream
){
    makeConvertedGpuPssm_gapless<ScoreType, PssmScoreType, TT, groupsize, numRegs>(result, hostPssmView, stream);
}



template<class ScoreType, class PssmScoreType, class TT, int groupsize, int numRegs>
void makeConvertedGpuPssm_sw(
    GpuPermutedPSSMforSW& result,
    const cuda::std::mdspan<TT, cuda::std::dextents<int,2>> hostPssmView, //numRows * queryLength elements, score to compare encoded subject letter s and query position q < queryLength must be at pssm[s][q]
    cudaStream_t stream
){
    using T = typename std::remove_cv<TT>::type;
    const int numRows = hostPssmView.extent(0);
    const int queryLength = hostPssmView.extent(1);

    std::vector<T> h_pssm(numRows * queryLength);
    int substitutionScoreBias = 0;
    if constexpr(IsVec4Type<ScoreType>::value){
        for(int s = 0; s < numRows; s++){
            for(int p = 0; p < queryLength; p++){
                int val = hostPssmView(s,p);
                substitutionScoreBias = std::min(substitutionScoreBias, val);
            }
        }
        substitutionScoreBias = std::abs(substitutionScoreBias);
        result.substitutionScoreBias_optional = substitutionScoreBias;
    }
    for(int s = 0; s < numRows; s++){
        for(int p = 0; p < queryLength; p++){
            h_pssm[s * queryLength + p] = hostPssmView(s,p) + substitutionScoreBias;
        }
    }
    T* d_pssm;
    CUDACHECK(helper_gpuMallocAsync((void**)&d_pssm, sizeof(T) * numRows * queryLength, stream));
    CUDACHECK(cudaMemcpyAsync(d_pssm, h_pssm.data(), sizeof(T) * numRows * queryLength, cudaMemcpyHostToDevice, stream));

    result.resize<PssmScoreType>(groupsize, numRegs, numRows, queryLength, stream);

    PSSM_2D_View<T> unpermutedView;
    unpermutedView.numRows = numRows;
    unpermutedView.numColumns = queryLength;
    unpermutedView.stride = queryLength;
    unpermutedView.data = reinterpret_cast<T*>(d_pssm);

    PSSM_2D_ModifiableView<PssmScoreType> resultView;
    resultView.numRows = numRows;
    resultView.numColumns = queryLength;
    resultView.stride = result.columnstride;
    resultView.data = reinterpret_cast<PssmScoreType*>(result.data.data());

    dim3 block(128,1,1);
    dim3 grid(SDIV(unpermutedView.numColumns, block.x), unpermutedView.numRows, 1);

    constexpr int accessSizeBytes = 16;
#if !defined(__METAL_BACKEND__) && !defined(__SIMD_BACKEND__)
    permute_PSSM_for_SW_kernel<accessSizeBytes><<<grid, block, 0, stream>>>(
        resultView,
        unpermutedView,
        numRegs,
        groupsize
    );
#endif
    CUDACHECKASYNC;

    CUDACHECK(helper_gpuFreeAsync((void*)d_pssm, stream));
}

template<class ScoreType, class PssmScoreType, class TT, int groupsize, int numRegs>
void makeConvertedGpuPssm(
    GpuPermutedPSSMforSW& result,
    const cuda::std::mdspan<TT, cuda::std::dextents<int,2>> hostPssmView, //numRows * queryLength elements, score to compare encoded subject letter s and query position q < queryLength must be at pssm[s][q]
    cudaStream_t stream
){
    makeConvertedGpuPssm_sw<ScoreType, PssmScoreType, TT, groupsize, numRegs>(result, hostPssmView, stream);
}


struct GpuConvertedPSSM{
    enum class Kind{Gapless, SW};

    Kind kindOfPssm;
    int alphabetSize;
    int numRegs;
    int groupsize;
    int columnstride;
    int allocatedNumRegs;
    SequenceLengthT queryLength;
    helpers::SimpleAllocationDevice<char, 0> data;

    std::optional<int> substitutionScoreBias_optional;
    int gaplessScoreBound = 0; // sum of positive column maxima, an upper bound on every gapless score

    template<class ScoreType, class ScalarPssmScoreType>
    void resize_gapless(int groupsize_, int numRegs_, int alphabetSize_, int queryLength_, int allocatedNumRegs_, cudaStream_t stream){
        kindOfPssm = Kind::Gapless;
        assert(512 % sizeof(ScalarPssmScoreType) == 0);
        static_assert(IsScalarType<ScalarPssmScoreType>::value);

        if((groupsize_ * numRegs_) % 2 == 1){
            throw std::runtime_error("GpuPermutedPSSMforGapless resize error. elements per row must be even");
        }
        groupsize = groupsize_;
        numRegs = numRegs_;
        alphabetSize = alphabetSize_;
        queryLength = queryLength_;
        allocatedNumRegs = allocatedNumRegs_;

        int lanes = 1;
        if(IsVec2Type<ScoreType>::value){
            lanes = 2;
        }else if(IsVec4Type<ScoreType>::value){
            lanes = 4;
        }

        const int workTileSize = (numRegs * groupsize * lanes);
        const int loadTileSize = (allocatedNumRegs * groupsize * lanes);
        const int numTiles = SDIV(queryLength, workTileSize);

        const int numAllocatedColumns = numTiles * loadTileSize;

        columnstride = numAllocatedColumns;
        // std::cout << "allocatedNumRegs " << allocatedNumRegs << ", workTileSize " << workTileSize << ", loadTileSize " << loadTileSize 
        //     << ", numTiles " << numTiles 
        //     << ", columnstride " << columnstride 
        //     << ", queryLength " << queryLength
        //     << "\n";
        //numPaddedColumns = (SDIV(groupsize * numItems * sizeof(PssmScoreType), 512) * 512) / sizeof(PssmScoreType);

        data.resize(sizeof(ScalarPssmScoreType) * alphabetSize * numAllocatedColumns);
        //init with 0 so oob elements won't contribute to the score
        cudaMemsetAsync(data.data(), 0, sizeof(char) * data.size(), stream);
    }

    template<class PssmScoreType>
    void resize_sw(int groupsize_, int numRegs_, int alphabetSize_, int queryLength_, cudaStream_t stream){
        kindOfPssm = Kind::SW;
        assert(512 % sizeof(PssmScoreType) == 0);

        groupsize = groupsize_;
        numRegs = numRegs_;
        alphabetSize = alphabetSize_;
        queryLength = queryLength_;

        const int tileSize = (numRegs * groupsize);
        const int numTiles = SDIV(queryLength, tileSize);
        columnstride = numTiles * tileSize;
        allocatedNumRegs = columnstride;
        //numPaddedColumns = (SDIV(groupsize * numItems * sizeof(PssmScoreType), 512) * 512) / sizeof(PssmScoreType);

        data.resize(sizeof(PssmScoreType) * alphabetSize * tileSize * numTiles);
        //init with 0 so oob elements won't contribute to the score
        cudaMemsetAsync(data.data(), 0, sizeof(char) * data.size(), stream);
    }

    template<class T>
    PSSM_2D_View<T> makeScalarView() const{
        PSSM_2D_View<T> view;
        view.numRows = alphabetSize;
        view.numColumns = columnstride;
        view.stride = columnstride;
        view.data = reinterpret_cast<const T*>(data.data());

        return view;
    }

    template<class T>
    PSSM_2D_View<T> makeVec2View() const{
        
        PSSM_2D_View<T> view;
        view.numRows = alphabetSize;
        if(kindOfPssm == Kind::Gapless){
            assert(columnstride % 2 == 0);
            view.numColumns = columnstride / 2;
            view.stride = columnstride / 2;
        }else{
            view.numColumns = columnstride;
            view.stride = columnstride;
        }
        view.data = reinterpret_cast<const T*>(data.data());

        return view;
    }

    template<class T>
    PSSM_2D_View<T> makeVec4View() const{

        PSSM_2D_View<T> view;
        view.numRows = alphabetSize;
        if(kindOfPssm == Kind::Gapless){
            assert(columnstride % 4 == 0);
            view.numColumns = columnstride / 4;
            view.stride = columnstride / 4;
        }else{
            view.numColumns = columnstride;
            view.stride = columnstride;
        }
        view.data = reinterpret_cast<const T*>(data.data());

        return view;
    }


    // PSSM_2D_View<half2> makeHalf2View() const{
    //     return makeVec2View<half2>();
    // }

    // PSSM_2D_View<short2> makeShort2View() const{
    //     return makeVec2View<short2>();
    // }

    // PSSM_2D_View<ScoreType_u8x4> makeUint8x4View() const{
    //     return makeVec4View<ScoreType_u8x4>();
    // }



};


template<class ScoreType, class PssmScoreType, class TT, int groupsize, int numRegs>
void makeConvertedGpuPssm_gapless(
    GpuConvertedPSSM& result,
    const cuda::std::mdspan<TT, cuda::std::dextents<int,2>> hostPssmView, //numRows * queryLength elements, score to compare encoded subject letter s and query position q < queryLength must be at pssm[s][q]
    cudaStream_t stream
){
    using T = typename std::remove_cv<TT>::type;
    const int numRows = hostPssmView.extent(0);
    const int queryLength = hostPssmView.extent(1);
 
    int workVectorLanesPerReg = 1;
    if(IsVec2Type<ScoreType>::value){
        workVectorLanesPerReg = 2;
    }else if(IsVec4Type<ScoreType>::value){
        workVectorLanesPerReg = 4;
    }

    //score half2, pssm half2 -> 2
    //score half2, pssm fp8x2 -> 2
    //score int8x4, pssm int8x4 -> 4
    int loadVectorLanesPerReg = 1;
    if(IsVec2Type<ScoreType>::value){
        loadVectorLanesPerReg = 2;
    }else if(IsVec4Type<ScoreType>::value){
        loadVectorLanesPerReg = 4;
    }

    using ScalarPssmScoreType = typename ScalarScoreType<PssmScoreType>::type;


    std::vector<T> h_pssm(numRows * queryLength);
    int substitutionScoreBias = 0;
    if constexpr(IsVec4Type<ScoreType>::value){
        for(int s = 0; s < numRows; s++){
            for(int p = 0; p < queryLength; p++){
                int val = hostPssmView(s,p);
                substitutionScoreBias = std::min(substitutionScoreBias, val);
            }
        }
        substitutionScoreBias = std::abs(substitutionScoreBias);
        result.substitutionScoreBias_optional = substitutionScoreBias;
    }
    for(int s = 0; s < numRows; s++){
        for(int p = 0; p < queryLength; p++){
            h_pssm[s * queryLength + p] = hostPssmView(s,p) + substitutionScoreBias;
        }
    }
    T* d_pssm;
    CUDACHECK(helper_gpuMallocAsync((void**)&d_pssm, sizeof(T) * numRows * queryLength, stream));
    CUDACHECK(cudaMemcpyAsync(d_pssm, h_pssm.data(), sizeof(T) * numRows * queryLength, cudaMemcpyHostToDevice, stream));


    constexpr int accessSizeBytes = 16;
    int relaxChunkSize = accessSizeBytes / sizeof(PssmScoreType);
    const int numRelaxChunks = SDIV(numRegs, relaxChunkSize);

    result.resize_gapless<ScoreType, ScalarPssmScoreType>(groupsize, numRegs, numRows, queryLength, relaxChunkSize * numRelaxChunks, stream);


    PSSM_2D_View<T> unpermutedView;
    unpermutedView.numRows = numRows;
    unpermutedView.numColumns = queryLength;
    unpermutedView.stride = queryLength;
    unpermutedView.data = reinterpret_cast<T*>(d_pssm);

    PSSM_2D_ModifiableView<ScalarPssmScoreType> resultView;
    resultView.numRows = numRows;
    resultView.numColumns = queryLength;
    resultView.stride = result.columnstride;
    resultView.data = reinterpret_cast<ScalarPssmScoreType*>(result.data.data());

    dim3 block(128,1,1);
    dim3 grid(SDIV(unpermutedView.numColumns, block.x), unpermutedView.numRows, 1);

#if !defined(__METAL_BACKEND__) && !defined(__SIMD_BACKEND__)
    permute_PSSM_for_gapless_kernel_withRelaxChunkSize<<<grid, block, 0, stream>>>(
        resultView,
        unpermutedView,
        numRegs,
        groupsize,
        relaxChunkSize, // number of regs per thread for which data is loaded per load instruction
        workVectorLanesPerReg, //number of packed elements per reg for gapless computation
        loadVectorLanesPerReg //number of packed elements per reg for data loading during computation
    );
#endif
    CUDACHECKASYNC;

    CUDACHECK(helper_gpuFreeAsync((void*)d_pssm, stream));
}


template<class ScoreType, class PssmScoreType, class TT, int groupsize, int numRegs>
void makeConvertedGpuPssm_sw(
    GpuConvertedPSSM& result,
    const cuda::std::mdspan<TT, cuda::std::dextents<int,2>> hostPssmView, //numRows * queryLength elements, score to compare encoded subject letter s and query position q < queryLength must be at pssm[s][q]
    cudaStream_t stream
){
    using T = typename std::remove_cv<TT>::type;
    const int numRows = hostPssmView.extent(0);
    const int queryLength = hostPssmView.extent(1);

    std::vector<T> h_pssm(numRows * queryLength);
    int substitutionScoreBias = 0;
    if constexpr(IsVec4Type<ScoreType>::value){
        for(int s = 0; s < numRows; s++){
            for(int p = 0; p < queryLength; p++){
                int val = hostPssmView(s,p);
                substitutionScoreBias = std::min(substitutionScoreBias, val);
            }
        }
        substitutionScoreBias = std::abs(substitutionScoreBias);
        result.substitutionScoreBias_optional = substitutionScoreBias;
    }
    for(int s = 0; s < numRows; s++){
        for(int p = 0; p < queryLength; p++){
            h_pssm[s * queryLength + p] = hostPssmView(s,p) + substitutionScoreBias;
        }
    }
    T* d_pssm;
    CUDACHECK(helper_gpuMallocAsync((void**)&d_pssm, sizeof(T) * numRows * queryLength, stream));
    CUDACHECK(cudaMemcpyAsync(d_pssm, h_pssm.data(), sizeof(T) * numRows * queryLength, cudaMemcpyHostToDevice, stream));

    result.resize_sw<PssmScoreType>(groupsize, numRegs, numRows, queryLength, stream);

    PSSM_2D_View<T> unpermutedView;
    unpermutedView.numRows = numRows;
    unpermutedView.numColumns = queryLength;
    unpermutedView.stride = queryLength;
    unpermutedView.data = reinterpret_cast<T*>(d_pssm);

    PSSM_2D_ModifiableView<PssmScoreType> resultView;
    resultView.numRows = numRows;
    resultView.numColumns = queryLength;
    resultView.stride = result.columnstride;
    resultView.data = reinterpret_cast<PssmScoreType*>(result.data.data());

    dim3 block(128,1,1);
    dim3 grid(SDIV(unpermutedView.numColumns, block.x), unpermutedView.numRows, 1);

    constexpr int accessSizeBytes = 16;
#if !defined(__METAL_BACKEND__) && !defined(__SIMD_BACKEND__)
    permute_PSSM_for_SW_kernel<accessSizeBytes><<<grid, block, 0, stream>>>(
        resultView,
        unpermutedView,
        numRegs,
        groupsize
    );
#endif
    CUDACHECKASYNC;

    CUDACHECK(helper_gpuFreeAsync((void*)d_pssm, stream));
}



struct GpuPlainPSSM{
    int numRows;
    int numColumns;
    helpers::SimpleAllocationDevice<char, 0> buffer;

    template<class T>
    void resize(int numRows_, int numColumns_, cudaStream_t){
        numRows = numRows_;
        numColumns = numColumns_;
        buffer.resize(sizeof(T) * numRows_ * numColumns_);
    }

    template<class T>
    T data(){
        return reinterpret_cast<T>(buffer.data());
    }

};


LIBMARV_NAMESPACE_WITH_NESTING_END


#endif
