//
// Created by lars on 26.05.15.
//

#ifndef MMSEQS_DISTANCECALCULATOR_H
#define MMSEQS_DISTANCECALCULATOR_H

#include <sstream>
#include <cstring>
#include <vector>

#include "simd.h"
#include "MathUtil.h"
#include "BaseMatrix.h"

class DistanceCalculator {
public:


    template<typename T>
    static unsigned int computeSubstituionDistance(const T *seq1,
                                                   const T *seq2,
                                                   const unsigned int length,
                                                   const char ** subMat, bool globalAlignment = false) {
        int max = 0;
        int score = 0;
        if (globalAlignment)
        {
            for(unsigned int pos = 0; pos < length; pos++){
                max += subMat[seq1[pos]][seq2[pos]];
            }
        } else {
            for(unsigned int pos = 0; pos < length; pos++){
                int curr = subMat[seq1[pos]][seq2[pos]];
                score = curr  + score;
                score = (score < 0) ? 0 : score;
                max = (score > max)? score : max;
            }
        }
        if (max<0)
            max = 0;
            
        return max;
    }



    struct LocalAlignment{
        int startPos;
        int endPos;
        unsigned int score;
        LocalAlignment(int startPos, int endPos, int score)
        : startPos(startPos), endPos(endPos), score(score)
        {}
        LocalAlignment() : startPos(-1), endPos(-1), score(-1)  {}

    };

    template<typename T>
    static LocalAlignment computeSubstituionStartEndDistance(const T *seq1,
                                                           const T *seq2,
                                                           const unsigned int length,
                                                           const char ** subMat) {
        int maxScore = 0;
        int maxEndPos = 0;
        int maxStartPos = 0;
        int minPos = 0;
        int maxMinPos = 0;
        int score = 0;
        for(unsigned int pos = 0; pos < length; pos++){
            int curr = subMat[seq1[pos]][seq2[pos]];
            score = curr  + score;
            const bool isMinScore = (score < 0);
            score =  (isMinScore) ? 0 : score;
            minPos = (isMinScore) ? pos : minPos;
            const bool isNewMaxScore = (score > maxScore);
            maxEndPos = (isNewMaxScore) ? pos : maxEndPos;
            maxStartPos = (isNewMaxScore) ? minPos : maxStartPos;
            maxScore = (isNewMaxScore)? score : maxScore;
        }
        return LocalAlignment(maxStartPos, maxEndPos, maxScore);
    }


    static unsigned int computeHammingDistance(const char *seq1, const char *seq2, unsigned int length){
        unsigned int diff = 0;
        unsigned int simdBlock = length/(VECSIZE_INT*4);
        simd_int * simdSeq1 = (simd_int *) seq1;
        simd_int * simdSeq2 = (simd_int *) seq2;
        for (unsigned int pos = 0; pos < simdBlock; pos++ ) {
            simd_int seq1vec = simdi_loadu(simdSeq1+pos);
            simd_int seq2vec = simdi_loadu(simdSeq2+pos);
            // int _mm_movemask_epi8(__m128i a) creates 16-bit mask from most significant bits of
            // the 16 signed or unsigned 8-bit integers in a and zero-extends the upper bits.
            simd_int seqComparision = simdi8_eq(seq1vec, seq2vec);
            int res = simdi8_movemask(seqComparision);
            diff += MathUtil::popCount(res);  // subtract positions that should not contribute to coverage
        }
         // compute missing rest
        for (unsigned int pos = simdBlock*(VECSIZE_INT*4); pos < length; pos++ ) {
            diff += (seq1[pos] == seq2[pos]);
        }
        return length-diff;
    }


    /*
     * Adapted from levenshtein.js (https://gist.github.com/andrei-m/982927)
     * Changed to hold only one row of the dynamic programing matrix
     * Copyright (c) 2011 Andrei Mackenzie
     * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
     * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
     * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
    static size_t levenshteinDistance(const std::string &s1, const std::string &s2) {
        size_t m = s1.length();
        size_t n = s2.length();

        if (m == 0)
            return n;
        if (n == 0)
            return m;

        // we don't need the full matrix just for the distance
        // we only need space for one row + the last value;
        size_t *row  = new size_t[m + 1];
        for (size_t i = 0; i <= m; ++i) {
            row[i] = i;
        }

        // fill in the rest of the matrix
        for (size_t i = 1; i <= n; i++) {
            size_t prev = i;
            for (size_t j = 1; j <= m; j++) {
                size_t val;
                if (s2[i - 1] == s1[j - 1]) {
                    val = row[j - 1]; // match
                } else {
                    val = std::min(row[j - 1] + 1, // substitution
                                   std::min(prev + 1,       // insertion
                                            row[j] + 1));   // deletion

                }
                row[j - 1] = prev;

                prev = val;
            }
            row[m] = prev;
        }

        // the last element of the matrix is the distance
        size_t distance = row[m];
        delete[] row;

        return distance;
    }

    static double keywordDistance(const std::string &s1, const std::string &s2) {
        std::stringstream seq1(s1);
        std::string segment;
        std::vector<std::string> kwlist1;
        while (std::getline(seq1, segment, ';')) {
            if (strcmp(segment.c_str(), "\n") != 0)
                kwlist1.push_back(segment);
        }
        std::vector<std::string> kwlist2;
        std::stringstream seq2(s2);
        while (std::getline(seq2, segment, ';')) {
            if (strcmp(segment.c_str(), "\n") != 0)
                kwlist2.push_back(segment);
        }
        double score = 0.0;
        for (std::string s1 : kwlist1) {
            for (std::string s2 : kwlist2) {
                if (strcmp(s1.c_str(), s2.c_str()) == 0) {
                    score++;
                    break;
                }
            }
        }
        score = ((score / std::max(kwlist1.size(), kwlist2.size())) +
                 (score / std::min(kwlist1.size(), kwlist2.size()))) / 2;

        return score;
    }

    void prepareGlobalAliParam(const BaseMatrix &subMat)
    {
        globalAliMu = 0;
        globalAliSigma = 0;

        for (size_t i = 0; i<subMat.alphabetSize - 1;i++)
        {

            for (size_t j = 0; j<subMat.alphabetSize - 1;j++)
            {
                globalAliMu += subMat.pBack[i] * subMat.pBack[j] * subMat.subMatrix[i][j];
            }
        }

        for (size_t i = 0; i<subMat.alphabetSize - 1;i++)
        {

            for (size_t j = 0; j<subMat.alphabetSize - 1;j++)
            {
                double distToMean = (subMat.subMatrix[i][j] - globalAliMu);
                globalAliSigma += subMat.pBack[i] * subMat.pBack[j] * distToMean*distToMean;
            }
        }
        globalAliSigma = sqrt(globalAliSigma);

    }

    double getPvalGlobalAli(float score,size_t len)
    {
        return 0.5 - 0.5*erf((score/len-globalAliMu)/(sqrt(2.0/sqrt((float)len))*globalAliSigma));
    }
private:
    
    float globalAliMu,globalAliSigma;
    

};

#endif //MMSEQS_DISTANCECALCULATOR_H
