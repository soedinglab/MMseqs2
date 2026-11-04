#ifndef COMMON_HPP
#define COMMON_HPP


#include <iostream>
#include <random>
#include <vector>
#include <omp.h>

#include <thrust/device_vector.h>

#include "../pssm.cuh"

constexpr int alphabetSize = 21;

constexpr bool subjectIsCaseSensitive = true;



struct Stuff1{
    size_t maxSequenceBytes = 512*1024ull*1024ull;
    size_t maxNumAlignments = 10000000;

    std::vector<std::vector<int>> substitutionMatrix2D;
    std::vector<int> substitutionMatrix;
    std::vector<std::int8_t> h_encodedSubjects;
    std::vector<std::int8_t> h_encodedQueries;
    thrust::device_vector<std::int8_t> d_subjects;
    thrust::device_vector<int> d_subjectLengths;
    thrust::device_vector<size_t> d_subjectBeginOffsets;
    thrust::device_vector<int> d_scoreOutput;
    thrust::device_vector<int> d_subjectEndPositionsExcl;
    thrust::device_vector<int> d_queryEndPositionsExcl;
};

__forceinline__
Stuff1 makeStuff1(){
    Stuff1 stuff;

    stuff.substitutionMatrix2D = std::vector<std::vector<int>>(alphabetSize, std::vector<int>(alphabetSize, -1));
    for(int i = 0; i < alphabetSize; i++){
        stuff.substitutionMatrix2D[i][i] = 3;
    }

    stuff.substitutionMatrix = std::vector<int>(alphabetSize * alphabetSize);
    for(int r = 0; r < alphabetSize; r++){
        for(int c = 0; c < alphabetSize; c++){
            stuff.substitutionMatrix[r * alphabetSize + c] = stuff.substitutionMatrix2D[r][c];
        }
    }


    constexpr size_t maxSequenceBytes = 512*1024ull*1024ull;
    constexpr size_t maxNumAlignments = 10000000;
    stuff.maxSequenceBytes = maxSequenceBytes;
    stuff.maxNumAlignments = maxNumAlignments;

    stuff.h_encodedSubjects = std::vector<std::int8_t>(maxSequenceBytes);
    stuff.h_encodedQueries = std::vector<std::int8_t>(maxSequenceBytes);

    #pragma omp parallel
    {
        std::mt19937 gen(42 + omp_get_thread_num());
        std::uniform_int_distribution<> dist(0, alphabetSize-1);

        #pragma omp for
        for(size_t i = 0; i < stuff.h_encodedSubjects.size(); i++){
            stuff.h_encodedSubjects[i] = dist(gen);
        }
        #pragma omp for
        for(size_t i = 0; i < stuff.h_encodedQueries.size(); i++){
            stuff.h_encodedQueries[i] = dist(gen);
        }
    }

    stuff.d_subjects = stuff.h_encodedSubjects;
    stuff.d_subjectLengths = thrust::device_vector<int>(maxNumAlignments);
    stuff.d_subjectBeginOffsets = thrust::device_vector<size_t>(maxNumAlignments);
    stuff.d_scoreOutput = thrust::device_vector<int>(maxNumAlignments);
    stuff.d_subjectEndPositionsExcl = thrust::device_vector<int>(maxNumAlignments);
    stuff.d_queryEndPositionsExcl = thrust::device_vector<int>(maxNumAlignments);   

    return stuff;
}




__forceinline__
int cpu_gapless(
    const libmarv::PSSM& pssm,
    int queryLength,
    const std::int8_t* const __restrict__ subjectData,
    const int subjectLength,
    bool printMatrix = false
){
    const int numColumns = queryLength+1;
    const int numRows = subjectLength+1;
    std::vector<int> M(numColumns * numRows, 0);

    int maxScore = 0;
    for(int r = 1; r < numRows; r++){
        const int subjectLetter = subjectData[r-1];
        for(int c = 1; c < numColumns; c++){
            const int diag = M[(r-1) * numColumns + (c-1)];
            const int substScore = pssm[subjectLetter][c-1];
            // std::cerr << substScore << " ";
            const int cell = max(diag + substScore, 0);
            M[r * numColumns + c] = cell;
            maxScore = max(maxScore, cell);
        }
        // std::cerr << "\n";
    }

    if(printMatrix){
        std::cerr << "gapless dp matrix\n";
        for(int r = 0; r < numRows; r++){
            for(int c = 0; c < numColumns; c++){
                const int val = M[r * numColumns + c];
                std::cerr << val << " ";
            }
            std::cerr << "\n";
        }
    }

    return maxScore;
}

#endif