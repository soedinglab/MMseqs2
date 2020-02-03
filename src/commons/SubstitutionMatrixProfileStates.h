#ifndef SubstitutionMatrixProfileStates_H
#define SubstitutionMatrixProfileStates_H
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
        num2aa[0] = 'A';
        num2aa[1] = 'C';
        num2aa[2] = 'D';
        num2aa[3] = 'E';
        num2aa[4] = 'F';
        num2aa[5] = 'G';
        num2aa[6] = 'H';
        num2aa[7] = 'I';
        num2aa[8] = 'K';
        num2aa[9] = 'L';
        num2aa[10] = 'M';
        num2aa[11] = 'N';
        num2aa[12] = 'P';
        num2aa[13] = 'Q';
        num2aa[14] = 'R';
        num2aa[15] = 'S';
        num2aa[16] = 'T';
        num2aa[17] = 'V';
        num2aa[18] = 'W';
        num2aa[19] = 'Y';
        num2aa[20] = 'X';
        num2aa[21] = 'Z';
        num2aa[22] = '[';
        num2aa[23] = '\\';
        num2aa[24] = ']';
        num2aa[25] = '^';
        num2aa[26] = '_';
        num2aa[27] = '`';
        num2aa[28] = 'a';
        num2aa[29] = 'b';
        num2aa[30] = 'c';
        num2aa[31] = 'd';
        alphabetSize = 21;
        initMatrixMemory(alphabetSize);

        for (int i = 0; i < alphabetSize; ++i){
            aa2num[static_cast<int>(num2aa[i])] = static_cast<unsigned char>(i);
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
        //this->num2aa[toIndex] = this->num2aa[fromIndex];

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
