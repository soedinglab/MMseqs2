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
                                    float **rMatrix, short **subMat,
                                    float bitFactor, float scoreBias, size_t maxSeqLen) {

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


        
        ps = new ProfileStates(this->pBack);
        this->scoreNormalization = ps->getScoreNormalization();
        this->bitFactor = bitFactor * scoreNormalization;
        //this->int2aa[toIndex] = this->int2aa[fromIndex];

        this->subMatrix = new short*[alphabetSize];
        for (int i = 0; i<alphabetSize; i++) {
            this->subMatrix[i] = new short[alphabetSize];
        }
        generateSubMatrix(probMatrix, rMatrix, this->subMatrix, this->subMatrix2Bit, alphabetSize, bitFactor, scoreBias);
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
