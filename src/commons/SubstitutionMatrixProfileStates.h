#ifndef SubstitutionMatrixProfileStates_H
#define SubstitutionMatrixProfileStates_H
#include <fstream>
#include <string>
#include <vector>

#include "BaseMatrix.h"
#include "Debug.h"
#include "ProfileStates.h"

class SubstitutionMatrixProfileStates : public BaseMatrix {
    public:
    SubstitutionMatrixProfileStates(std::string matrixName,
                                    double **probMatrix, double * pBack,
                                    float **rMatrix, float bitFactor, float scoreBias,
                                    int libAlphabetSize) {
//        alphabetSize = 32;
        int2aa[0] = 'A';
        int2aa[1] = 'C';
        int2aa[2] = 'D';
        int2aa[3] = 'E';
        int2aa[4] = 'F';
        int2aa[5] = 'G';
        int2aa[6] = 'H';
        int2aa[7] = 'I';
        int2aa[8] = 'K';
        int2aa[9] = 'L';
        int2aa[10] = 'M';
        int2aa[11] = 'N';
        int2aa[12] = 'P';
        int2aa[13] = 'Q';
        int2aa[14] = 'R';
        int2aa[15] = 'S';
        int2aa[16] = 'T';
        int2aa[17] = 'V';
        int2aa[18] = 'W';
        int2aa[19] = 'Y';
        int2aa[20] = 'X';
        int2aa[21] = 'Z';
        int2aa[22] = '[';
        int2aa[23] = '\\';
        int2aa[24] = ']';
        int2aa[25] = '^';
        int2aa[26] = '_';
        int2aa[27] = '`';
        int2aa[28] = 'a';
        int2aa[29] = 'b';
        int2aa[30] = 'c';
        int2aa[31] = 'd';
        alphabetSize = 21;
        initMatrixMemory(alphabetSize);

        for (int i = 0; i < alphabetSize; ++i){
            aa2int[(int)int2aa[i]] = i;
        }

        this->matrixName = matrixName;
        this->origAlphabetSize = alphabetSize;
        this->scoreBias = scoreBias;
        for(int i = 0; i < this->alphabetSize; i++) {
            for (int j = 0; j < this->alphabetSize; j++) {
                this->probMatrix[i][j] = probMatrix[i][j];
            }
        }

        for(int i = 0; i < this->alphabetSize; i++) {
                this->pBack[i] = pBack[i];
        }

        for(int i = 0; i < this->alphabetSize; i++) {
            for (int j = 0; j < this->alphabetSize; j++) {
                this->subMatrixPseudoCounts[i][j] = rMatrix[i][j];
            }
        }

        ps = new ProfileStates(libAlphabetSize,this->pBack);
        this->scoreNormalization = ps->getScoreNormalization();
        this->bitFactor = bitFactor * scoreNormalization;
        //this->int2aa[toIndex] = this->int2aa[fromIndex];

        this->subMatrix = new short*[alphabetSize];
        for (int i = 0; i<alphabetSize; i++) {
            this->subMatrix[i] = new short[alphabetSize];
        }
        generateSubMatrix(probMatrix, rMatrix, this->subMatrix, alphabetSize, bitFactor, scoreBias);
        // remove X
        this->alphabetSize = ps->getAlphSize();
    }

    ~SubstitutionMatrixProfileStates(){
        delete ps;
        // delete too much memory in the BaseMatrix destructor if I do not set this back to org. alph. size
        alphabetSize = origAlphabetSize;
    };

    float * getProfileVectorForState(size_t k){
        return ps->getProfile(k);
    }
    virtual float getBitFactor() {return bitFactor; }
    float getScoreNormalization(){ return scoreNormalization; }
    virtual float getScoreBias() {return scoreBias; }

    float scoreState(float *profile, float *pav, size_t k) {
        return ps->score(profile, pav, k);
    }

private:
    ProfileStates * ps;
    int origAlphabetSize;
    float bitFactor;
    float scoreBias;
    float scoreNormalization;
};
#endif
