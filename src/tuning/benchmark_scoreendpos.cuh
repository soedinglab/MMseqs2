
#ifndef BENCHMARK_SCORE_ENDPOS_CUH
#define BENCHMARK_SCORE_ENDPOS_CUH

#include "benchmarkdata.hpp"
#include "common.hpp"
#include "configtuples.hpp"


#include <iostream>
#include <vector>
#include <random>
#include <cstddef>
#include <utility>
#include <string>
#include <algorithm>
#include <string>
#include <thread>
#include <chrono>


#include "../cuda_errorcheck.cuh"

#include "../convert.cuh"
#include "../pssm.cuh"
#include "../offset_iterator.cuh"
#include "../util.cuh"

#include "../aligners/score_with_endpos_alignment_interfaces.cuh"
#include "../aligners/score_with_endpos_alignment.cuh"
// #include "../alignment_algorithms/alignment_interface_data.cuh"

using libmarv::Approach;
using libmarv::Datatype;
using libmarv::AlignmentAlgorithmE;
using libmarv::IsScalarType;
using libmarv::IsVec2Type;
using libmarv::IsVec4Type;
using libmarv::ScoreTypeFromDatatypeE;

#include <omp.h>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>


namespace score_endpos{

template<
    int alphabetSize, 
    AlignmentAlgorithmE algorithm, 
    int blocksize,
    int groupsize,
    int numItems,
    bool hasSingleTile,
    bool hasMultiTile,
    class ApproachAndType
>
std::unique_ptr<libmarv::ScoreWithEndPosAlignmentInterface> makeAligner(){
    if constexpr(algorithm == AlignmentAlgorithmE::Gapless){                    
        if constexpr(IsScalarType<typename ApproachAndType::ScoreType>::value){
            using AlignerClass = libmarv::GaplessAlignment_endPos_scalarscore<
                alphabetSize, 
                ApproachAndType,
                blocksize, groupsize, numItems,
                hasSingleTile,
                hasMultiTile
            >;
            return std::make_unique<AlignerClass>();
        }else if constexpr(IsVec2Type<typename ApproachAndType::ScoreType>::value){
            using AlignerClass = libmarv::GaplessAlignment_endPos_vec2score<
                alphabetSize, 
                ApproachAndType,
                blocksize, groupsize, numItems,
                hasSingleTile,
                hasMultiTile
            >;
            return std::make_unique<AlignerClass>();
        }else if constexpr(IsVec4Type<typename ApproachAndType::ScoreType>::value){

            using AlignerClass = libmarv::GaplessAlignment_endPos_vec4score<
                alphabetSize, 
                ApproachAndType,
                blocksize, groupsize, numItems,
                hasSingleTile,
                hasMultiTile
            >;
            return std::make_unique<AlignerClass>();
        }else{
            throw std::runtime_error("unexpected score type");
        }
    }else if constexpr(algorithm == AlignmentAlgorithmE::SW){
        using AlignerClass = libmarv::SWAlignment_endPos_scalarscore<
            alphabetSize, 
            ApproachAndType,
            blocksize, groupsize, numItems,
            hasSingleTile,
            hasMultiTile
        >;
        return std::make_unique<AlignerClass>();
    }else{
        throw std::runtime_error("unexpected config kind");
    }
}

template<
    int blocksize,
    int groupsize,
    int numItems,
    int tilesize,
    class ApproachAndType,
    class SleepDuration
>
BenchmarkData benchmark_score_endpos(
    int timingIterations, 
    SleepDuration sleepAmount,
    bool checkResults,
    Stuff1& stuff,
    int queryLength,
    libmarv::ScoreWithEndPosAlignmentInterface* aligner
){
    constexpr AlignmentResultKind resultkind = AlignmentResultKind::ScoreEndpos;

    std::cerr << blocksize << " " <<  groupsize << " " << numItems << " " << to_string(ApproachAndType::datatype) << " " << to_string(ApproachAndType::approach) 
        << " " << to_string(resultkind)  << ", query length " << queryLength << "\n";


    if (ApproachAndType::datatype != aligner->getScoreDatatype()){
        return makeNotExecutedBenchmarkData<blocksize, groupsize, numItems, tilesize, ApproachAndType, resultkind>();
    }

    int ccMajor = 0;
    cudaDeviceGetAttribute(&ccMajor, cudaDevAttrComputeCapabilityMajor, 0);
    if(ccMajor < 9 && ApproachAndType::datatype == Datatype::Short2){
        return makeNotExecutedBenchmarkData<blocksize, groupsize, numItems, tilesize, ApproachAndType, resultkind>();
    }
    if(ccMajor < 12 && ApproachAndType::datatype == Datatype::UInt8x4){
        return makeNotExecutedBenchmarkData<blocksize, groupsize, numItems, tilesize, ApproachAndType, resultkind>();
    }
    
    cudaStream_t stream = cudaStreamLegacy;
    
    BenchmarkData benchmarkData = makeNotExecutedBenchmarkData<blocksize, groupsize, numItems, tilesize, ApproachAndType, resultkind>();

    
    const int subjectLength = tilesize;
    const int subjectLengthPadded = ((subjectLength+3)/4) * 4;
    const int numAlignments = std::min(stuff.maxNumAlignments, stuff.maxSequenceBytes / subjectLengthPadded);

    size_t cells = size_t(numAlignments) * size_t(subjectLength) * size_t(queryLength);

    std::vector<size_t> h_subjectBeginOffsets(numAlignments);
    for(int i = 0; i < numAlignments; i++){
        h_subjectBeginOffsets[i] = size_t(i) * subjectLengthPadded;
        // h_subjectBeginOffsets[i] = size_t(53) * subjectLengthPadded;
    }
    // std::copy(
    //     stuff.h_encodedSubjects.begin() + size_t(53) * subjectLengthPadded, 
    //     stuff.h_encodedSubjects.begin() + size_t(54) * subjectLengthPadded,
    //     stuff.h_encodedSubjects.begin()
    // );
    // stuff.d_subjects = stuff.h_encodedSubjects;


    thrust::copy(h_subjectBeginOffsets.begin(), h_subjectBeginOffsets.end(), stuff.d_subjectBeginOffsets.begin());
    thrust::fill(stuff.d_subjectLengths.begin(), stuff.d_subjectLengths.end(), subjectLength);
    thrust::fill(stuff.d_scoreOutput.begin(), stuff.d_scoreOutput.end(), 0);

    std::vector<int8_t> h_pssm(alphabetSize * queryLength);
    for(int s = 0; s < alphabetSize; s++){
        for(int p = 0; p < queryLength; p++){
            h_pssm[s * queryLength + p] = stuff.substitutionMatrix2D[s][stuff.h_encodedQueries[p]];            
        }
    }


    auto hostFullQueryPSSM = libmarv::PSSM::fromPSSM(nullptr, queryLength, h_pssm.data(), alphabetSize);


    // std::cout << "hostFullQueryPSSM\n";
    // for(int s = 0; s < alphabetSize; s++){
    //     for(int p = 0; p < queryLength; p++){
    //         std::cout << hostFullQueryPSSM[s][p] << " ";
    //     }
    //     std::cout << "\n";
    // }

    cuda::std::mdspan<std::int8_t, cuda::std::dextents<int,2>> hostPssmView(h_pssm.data(), alphabetSize, queryLength);
    libmarv::GpuConvertedPSSM gpuConvertedPSSM;
    aligner->makeGpuPssm(gpuConvertedPSSM, hostPssmView, stream);

    libmarv::OneToAllInputDataPSSM inputData;
    inputData.d_subjects = reinterpret_cast<const std::int8_t*>(stuff.d_subjects.data().get());
    inputData.d_subjectOffsets = libmarv::makeCustomOffsetIterator(stuff.d_subjectBeginOffsets.data().get(), false);
    inputData.d_subjectLengths = stuff.d_subjectLengths.data().get();
    inputData.queryLength = queryLength;
    inputData.numAlignments = numAlignments;
    inputData.maximumSubjectLength = subjectLength;
    inputData.indexIndirection = libmarv::PositionsIterator::fromCountingIterator(0);


    cudaEvent_t startevent,stopevent;
    float elapsed_ms = 0;
    CUDACHECK(cudaEventCreate(&startevent));
    CUDACHECK(cudaEventCreate(&stopevent));

    thrust::device_vector<char> d_temp;
    if(!aligner->isSingleTile(inputData)){
        size_t bytes = aligner->getMinimumSuggestedTempBytes_multiTile(inputData);

        size_t free,total;
        CUDACHECK(cudaMemGetInfo(&free, &total));
        if(free > (4 << 20)){
            free -= (4 << 20);
        }
        bytes = std::min(free, bytes);
        d_temp.resize(bytes, thrust::no_init);
    }




    std::vector<int> referenceScores(numAlignments);

    if(checkResults){    
        #pragma omp parallel for
        for(int i = 0; i < numAlignments; i++){
            referenceScores[i] = cpu_gapless(
                hostFullQueryPSSM,
                queryLength,
                stuff.h_encodedSubjects.data() + h_subjectBeginOffsets[i],
                subjectLength,
                false
            );
            // std::cout << "cpu h_subjectBeginOffsets = " << h_subjectBeginOffsets[i] << "\n";
            // for(int k = 0; k < subjectLength; k++){
            //     std::cout << int((stuff.h_encodedSubjects.data() + h_subjectBeginOffsets[i])[k]) << " ";
            // }
            // std::cout << "\n";
        }

    }

    try{
        std::vector<double> gcupsVec;
        for(int iteration = 0; iteration < timingIterations; iteration++){
            if(aligner->isSingleTile(inputData)){
                CUDACHECK(cudaEventRecord(startevent));
                aligner->scoreEndPos_singleTile(
                    stuff.d_scoreOutput.data().get(),
                    stuff.d_subjectEndPositionsExcl.data().get(),
                    stuff.d_queryEndPositionsExcl.data().get(),
                    inputData,
                    gpuConvertedPSSM,
                    libmarv::GapScoreArgs{},
                    stream
                );

                CUDACHECK(cudaEventRecord(stopevent));
            }else{
                CUDACHECK(cudaEventRecord(startevent));
                aligner->scoreEndPos_multiTile(
                    d_temp.data().get(),
                    d_temp.size(),
                    stuff.d_scoreOutput.data().get(),
                    stuff.d_subjectEndPositionsExcl.data().get(),
                    stuff.d_queryEndPositionsExcl.data().get(),
                    inputData,
                    gpuConvertedPSSM,
                    libmarv::GapScoreArgs{},
                    stream
                );

                CUDACHECK(cudaEventRecord(stopevent));
            }
            CUDACHECK(cudaDeviceSynchronize());
            CUDACHECK(cudaEventElapsedTime(&elapsed_ms, startevent, stopevent));
            // std::cout << "cells " << cells << ", elapsed_ms " << elapsed_ms << "\n";
            double GCUPS = (double(cells) / 1000./ 1000./ 1000.) / (elapsed_ms / 1000.);
            gcupsVec.push_back(GCUPS);
            std::this_thread::sleep_for(sleepAmount);

            if(checkResults){
                thrust::host_vector<int> h_scoreOutput = stuff.d_scoreOutput;

                int errors = 0;
                // for(int i = 0; i < std::min(10, numAlignments); i++){
                //     std::cout << referenceScores[i] << " " << h_scoreOutput[i] << "\n";
                // }
                for(int i = 0; i < numAlignments; i++){
                    if(referenceScores[i] != h_scoreOutput[i]/* && h_scoreOutput[i] < 2048*/){
                        errors++;
                        if(errors == 1){
                            std::cerr << "Error subject " << i << ", reference score " << referenceScores[i] << ", gpu score " << h_scoreOutput[i];
                            std::cerr << " subject length " << subjectLength << ", query length " << queryLength << "\n";
                        }
                    }
                }
                if(errors > 0){
                    std::cerr << errors << " errors in total\n";
                }else{
                    std::cerr << "gpu scores ok\n";
                }
            }
        }

        std::sort(gcupsVec.begin(), gcupsVec.end());
        //erase slowest and fastest run
        if(gcupsVec.size() > 2){
            gcupsVec.erase(gcupsVec.begin());
            gcupsVec.erase(gcupsVec.begin() + gcupsVec.size()-1);
        }
        double avggcups = std::reduce(gcupsVec.begin(), gcupsVec.end()) / gcupsVec.size();
        double mingcups = 0;
        double maxgcups = 0;
        if(gcupsVec.size() > 0){
            mingcups = *std::min_element(gcupsVec.begin(), gcupsVec.end());
            maxgcups = *std::max_element(gcupsVec.begin(), gcupsVec.end());
        }
        benchmarkData.mingcups = mingcups;
        benchmarkData.maxgcups = maxgcups;
        benchmarkData.avggcups = avggcups;
        benchmarkData.numTiles = ((queryLength + tilesize - 1) / tilesize);
    }catch(...){ 
        //catch any exceptions to allow other configs to run. do not report performance numbers
        benchmarkData = makeNotExecutedBenchmarkData<blocksize, groupsize, numItems, tilesize, ApproachAndType, resultkind>();
    }

    cudaEventDestroy(startevent);
    cudaEventDestroy(stopevent);


    return benchmarkData;
}







template<
    class TupleOfICBlockSizes,
    class TupleOfICGroupSizes,
    class TupleOfICNumItems,
    class TupleOfApproachAndDatatype,
    bool singleTile
>
struct ForEachCombination_call_singleconfigFunction{
    static constexpr int sizeA = std::tuple_size_v<TupleOfICBlockSizes>;
    static constexpr int sizeB = std::tuple_size_v<TupleOfICGroupSizes>;
    static constexpr int sizeC = std::tuple_size_v<TupleOfICNumItems>;
    static constexpr int sizeD = std::tuple_size_v<TupleOfApproachAndDatatype>;
    static constexpr int numCombinations = sizeA * sizeB * sizeC * sizeD;

    static constexpr size_t getIndexA(size_t I){
        return I / (sizeB * sizeC * sizeD);
    };

    static constexpr size_t getIndexB(size_t I){
        const auto remA = I % (sizeB * sizeC * sizeD);
        return remA / (sizeC * sizeD);
    };

    static constexpr size_t getIndexC(size_t I){
        const auto remA = I % (sizeB * sizeC * sizeD);
        const auto remB = remA % (sizeC * sizeD);
        return remB / (sizeD);
    };

    static constexpr size_t getIndexD(size_t I){
        const auto remA = I % (sizeB * sizeC * sizeD);
        const auto remB = remA % (sizeC * sizeD);
        return remB % (sizeD);
    };

    template<size_t index, class SleepDuration>
    BenchmarkData helper2(
        int timingIterations, 
        SleepDuration sleepAmount,
        bool checkResults,
        int tiles,
        Stuff1& stuff
    ){
        constexpr int blocksize = std::tuple_element_t<getIndexA(index), TupleOfICBlockSizes>::value;
        constexpr int groupsize = std::tuple_element_t<getIndexB(index), TupleOfICGroupSizes>::value;
        constexpr int numItems = std::tuple_element_t<getIndexC(index), TupleOfICNumItems>::value;
        //constexpr auto approachAndType = std::tuple_element_t<getIndexD(index), TupleOfApproachAndDatatype>::value;
        using ApproachAndType = std::tuple_element_t<getIndexD(index), TupleOfApproachAndDatatype>;
        using ScoreType = typename ScoreTypeFromDatatypeE<ApproachAndType::datatype>::type;

        assert(!singleTile || tiles == 1);
        if(singleTile && tiles != 1){
            throw std::runtime_error("singleTile && tiles != 1");
        }

        auto aligner = makeAligner<
            21, 
            AlignmentAlgorithmE::Gapless, 
            blocksize,
            groupsize,
            numItems,
            singleTile,
            !singleTile,
            ApproachAndType
        >();
        
        constexpr int valuesPerItem = IsScalarType<ScoreType>::value ? 1 : (IsVec2Type<ScoreType>::value ? 2 : (IsVec4Type<ScoreType>::value ? 4 : 0));
        static_assert(valuesPerItem > 0);

        constexpr int tilesize = valuesPerItem * groupsize * numItems;

        const int queryLength = tilesize * tiles;

        return benchmark_score_endpos<blocksize, groupsize, numItems, tilesize, ApproachAndType, SleepDuration>(
            timingIterations, sleepAmount, checkResults, stuff, queryLength, aligner.get());
        
    }

    template <class SleepDuration, size_t ... Is>    
    std::vector<BenchmarkData> helper1(
        std::index_sequence<Is...>, 
        int timingIterations, 
        SleepDuration sleepAmount,
        bool checkResults,
        int tiles,
        Stuff1& stuff
    ){
        std::vector<BenchmarkData> result;

        ((result.push_back(helper2<Is, SleepDuration>(timingIterations, sleepAmount, checkResults, tiles, stuff))), ...) ;

        return result;
    }

    template<class SleepDuration>
    std::vector<BenchmarkData> operator()(
        int timingIterations, 
        SleepDuration sleepAmount,
        bool checkResults,
        int tiles
    ){
        Stuff1 stuff = makeStuff1();
        auto indexList = std::make_index_sequence<numCombinations>();

        auto results = helper1(indexList, timingIterations, sleepAmount, checkResults, tiles, stuff);

        return results;
    }
};

} //namespace score_endpos




#endif