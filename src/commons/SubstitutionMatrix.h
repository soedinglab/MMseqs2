#ifndef SUBSTITUTION_MATRIX_H
#define SUBSTITUTION_MATRIX_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
//  
// Represents a simple amino acid substitution matrix. 
// Can read a straight-forward substitution matrix file (.mat ending of the file) and probabilities file (.out ending of the file).
// If a probabilities file is given, it calculates a biased matrix (produces shorter = more precise alignments).
//

#include <cstddef>
#include "BaseMatrix.h"
#include "ProfileStates.h"

class SubstitutionMatrix: public BaseMatrix {

    public:
        SubstitutionMatrix(const char *filename, float bitFactor, float scoreBias);

        virtual ~SubstitutionMatrix();

        virtual float getBitFactor() {return bitFactor; }
    
        virtual double getBackgroundProb(size_t aa_index) { return pBack[aa_index]; }

        static void calcLocalAaBiasCorrection(const BaseMatrix *m ,const unsigned char *int_sequence, const int N, float *compositionBias);
        static void calcProfileProfileLocalAaBiasCorrection(short *profileScores,
                                                const size_t profileAASize,
                                                const int N,
                                                size_t alphabetSize);
        static void calcProfileProfileLocalAaBiasCorrectionAln(int8_t *profileScores,
                                                             int N,
                                                             size_t alphabetSize,
                                                             BaseMatrix *subMat);
        static void calcGlobalAaBiasCorrection(const BaseMatrix * m,
                                               short *profileScores,
                                               float *pNullBuffer,
                                               const size_t profileAASize,
                                               const int N);
        bool estimateLambdaAndBackground(const double ** mat, int alphabetSize, double * pBack, double & lambda);


        void setupLetterMapping();


        struct FastMatrix{
            const char ** matrix;
            const char * matrixData;
            const size_t asciiStart;
            FastMatrix(const char ** matrix, const char * matrixData, const size_t asciiStart):
                    matrix(matrix), matrixData(matrixData), asciiStart(asciiStart)
            {}
        };

        // build matrix from ~ (=0) to ~(=122)
        static FastMatrix createAsciiSubMat(BaseMatrix & submat){
            const size_t asciiStart = 0;
            const size_t asciiEnd = 'z'+1;
            const size_t range = asciiEnd-asciiStart;
            char ** matrix = new char *[range];
            char * matrixData = new char[range*range];
            for(size_t i = 0; i < range; i++) {
                matrix[i] = matrixData+(i*range);
                int curr_i = static_cast<int>(submat.aa2num[asciiStart+i]);
                for (size_t j = 0; j < range; j++) {
                    int curr_j = static_cast<int>(submat.aa2num[asciiStart+j]);
                    matrix[i][j] = static_cast<char>(submat.subMatrix[curr_i][curr_j]);
                }
            }
            return FastMatrix((const char**) matrix,
                              (const char*) matrixData,
                              asciiStart);
        }

private:
    int parseAlphabet(char * word, char * num2aa, int * aa2num);

    void readProbMatrix(const std::string &matrixData, bool containsX);

    std::pair<int, bool>  setAaMappingDetectAlphSize(std::string &matrixData);

    const float bitFactor;
};

#endif
