#ifndef SubstitutionMatrixWithoutX_H
#define SubstitutionMatrixWithoutX_H
#include <fstream>
#include <string>
#include <vector>

#include "BaseMatrix.h"
#include "Debug.h"

class SubstitutionMatrixWithoutX : public BaseMatrix {
    public:
    SubstitutionMatrixWithoutX(double **probMatrix, float **rMatrix, short **subMat, float bitFactor) {
        // copy data
        for(int i = 0; i < this->alphabetSize; i++) {
            for (int j = 0; j < this->alphabetSize; j++) {
                this->probMatrix[i][j] = probMatrix[i][j];
            }
        }
//        for(size_t i = 0; i < this->alphabetSize; i++) {
//            for (size_t j = 0; j < this->alphabetSize; j++) {
//                this->subMatrixPseudoCounts[i][j] = rMatrix[i][j];
//            }
//        }
        for(int i = 0; i < this->alphabetSize; i++) {
            for (int j = 0; j < this->alphabetSize; j++) {
                this->subMatrix[i][j] = subMat[i][j];
            }
        }
        for (int i = 0; i<alphabetSize; i++)
            delete [] subMatrix[i];
        delete [] subMatrix;

        // remove X
        this->alphabetSize = 20;
        //int toIndex =  this->aa2int[(int)'X'];
        //int fromIndex = this->aa2int[(int)'L'];
        this->aa2int[(int)'X'] = this->aa2int[(int)'L'];
        //this->int2aa[toIndex] = this->int2aa[fromIndex];

        this->subMatrix = new short*[alphabetSize];
        for (int i = 0; i<alphabetSize; i++)
            this->subMatrix[i] = new short[alphabetSize];

        generateSubMatrix(probMatrix, rMatrix, this->subMatrix, this->subMatrix2Bit, alphabetSize, bitFactor, 0.0);

    }
    ~SubstitutionMatrixWithoutX(){};

};
#endif
