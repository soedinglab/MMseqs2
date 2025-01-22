#ifndef PSSM_CUH
#define PSSM_CUH

#include "config.hpp"
#include "types.hpp"
#include "convert.cuh"
#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/simple_allocation.cuh"

#include <vector>
#include <cassert>

namespace cudasw4{

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
    const int group_size
) {
    static_assert(accessSizeBytes == 4 || accessSizeBytes == 8 || accessSizeBytes == 16); //float, float2, float4
    constexpr int numFloatsPerAccess = accessSizeBytes / 4;
    constexpr int numHalfsPerAccess = numFloatsPerAccess * 2;

    const int thid = threadIdx.x + blockIdx.x*blockDim.x;

    const int numColumns = inputView.numColumns;
    const int numRows = inputView.numRows;

    const int tileSize = (numRegs * group_size * 2);

    for(int inputRow = blockIdx.y; inputRow < numRows; inputRow += gridDim.y){

        for (int inputCol = thid; inputCol < numColumns; inputCol += blockDim.x*gridDim.x) {
            const int tileId = inputCol / tileSize;
            const int tileColumnOffset = tileId * tileSize;
            const int columnInTile = inputCol - tileColumnOffset;

            const int l_2 = tileSize/2;
            const int offset = (2*(columnInTile%l_2) + columnInTile/l_2)%numHalfsPerAccess;
            const int thread = (columnInTile%l_2)/numRegs;
            const int part = (columnInTile%numRegs)/numFloatsPerAccess;
            const int resultCol = numHalfsPerAccess*thread + offset + part*numHalfsPerAccess*group_size;

            resultView[inputRow][tileColumnOffset + resultCol] = inputView[inputRow][inputCol];
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
    int group_size;
    int columnstride;
    SequenceLengthT queryLength;
    helpers::SimpleAllocationDevice<char, 0> data;

    template<class T>
    void resize(int group_size_, int numRegs_, int alphabetSize_, int queryLength_, cudaStream_t stream){
        assert(512 % sizeof(T) == 0);

        if((group_size_ * numRegs_) % 2 == 1){
            throw std::runtime_error("GpuPermutedPSSMforGapless resize error. elements per row must be even");
        }
        group_size = group_size_;
        numRegs = numRegs_;
        alphabetSize = alphabetSize_;
        queryLength = queryLength_;

        const int tileSize = (numRegs * group_size * 2);
        const int numTiles = SDIV(queryLength, tileSize);
        columnstride = numTiles * tileSize;
        //numPaddedColumns = (SDIV(groupsize * numItems * sizeof(PssmScoreType), 512) * 512) / sizeof(PssmScoreType);

        data.resize(sizeof(T) * alphabetSize * tileSize * numTiles);
        //init with 0 so oob elements won't contribute to the score
        cudaMemsetAsync(data.data(), 0, sizeof(char) * data.size(), stream);
    }


    PSSM_2D_View<half2> makeHalf2View() const{
        assert(columnstride % 2 == 0);

        PSSM_2D_View<half2> view;
        view.numRows = alphabetSize;
        view.numColumns = columnstride / 2;
        view.stride = columnstride / 2;
        view.data = reinterpret_cast<const half2*>(data.data());

        return view;
    }

    PSSM_2D_View<short2> makeShort2View() const{
        assert(columnstride % 2 == 0);

        PSSM_2D_View<short2> view;
        view.numRows = alphabetSize;
        view.numColumns = columnstride / 2;
        view.stride = columnstride / 2;
        view.data = reinterpret_cast<const short2*>(data.data());

        return view;
    }

    template<class OutputT, int accessSizeBytes, class InputT>
    void fromGpuPSSMView(PSSM_2D_View<InputT> inputView, int group_size_, int numRegs_, cudaStream_t stream){
        resize<OutputT>(group_size_, numRegs_, inputView.numRows, inputView.numColumns, stream);

        PSSM_2D_ModifiableView<OutputT> resultView;
        resultView.numRows = alphabetSize;
        resultView.numColumns = queryLength;
        resultView.stride = columnstride;
        resultView.data = reinterpret_cast<OutputT*>(data.data());

        dim3 block(128,1,1);
        dim3 grid(SDIV(inputView.numColumns, block.x), inputView.numRows, 1);

        permute_PSSM_for_gapless_kernel<accessSizeBytes><<<grid, block, 0, stream>>>(
            resultView,
            inputView,
            numRegs,
            group_size
        ); CUERR;
    }
};

struct GpuPermutedPSSMforSW{
    int alphabetSize;
    int numRegs;
    int group_size;
    int columnstride;
    SequenceLengthT queryLength;
    helpers::SimpleAllocationDevice<char, 0> data;

    template<class T>
    void resize(int group_size_, int numRegs_, int alphabetSize_, int queryLength_, cudaStream_t stream){
        assert(512 % sizeof(T) == 0);

        group_size = group_size_;
        numRegs = numRegs_;
        alphabetSize = alphabetSize_;
        queryLength = queryLength_;

        const int tileSize = (numRegs * group_size);
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
    void fromGpuPSSMView(PSSM_2D_View<InputT> inputView, int group_size_, int numRegs_, cudaStream_t stream){
        resize<OutputT>(group_size_, numRegs_, inputView.numRows, inputView.numColumns, stream);

        PSSM_2D_ModifiableView<OutputT> resultView;
        resultView.numRows = alphabetSize;
        resultView.numColumns = queryLength;
        resultView.stride = columnstride;
        resultView.data = reinterpret_cast<OutputT*>(data.data());

        dim3 block(128,1,1);
        dim3 grid(SDIV(inputView.numColumns, block.x), inputView.numRows, 1);

        permute_PSSM_for_SW_kernel<accessSizeBytes><<<grid, block, 0, stream>>>(
            resultView,
            inputView,
            numRegs,
            group_size
        ); CUERR;
    }
};


}


#endif
