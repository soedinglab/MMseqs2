#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include <iostream>
#include <cmath>
#include <cstdio>

class BaseMatrix{
    public:
        BaseMatrix();

        virtual ~BaseMatrix();

        /*contains amino acid to int mapping*/
        int*  aa2int;
        
        /*contains int to amino acid mapping*/
        char* int2aa;
        
        /* size of alphabet*/
        int alphabetSize;

        // substitution matrix
        short** subMatrix;

        // joint probability matrix
        double** probMatrix;

        // background probabilities of the amino acids
        double* pBack;

        // print the substitution matrix
        static void print(short** matrix, char* int2aa, int size);

        static void print(double** matrix, char* int2aa, int size);

        // generate the substitution matrix given the probability matrix, background probabilities and the alphabet size
        static void generateSubMatrix(double ** probMatrix, double ** subMatrix, int size, double bitFactor = 1.0, double scoringBias = 0.0);

        // generate a short data type substitution matrix
        static void generateSubMatrix(double ** probMatrix, short ** subMatrix, int size, double bitFactor = 1.0, double scoringBias = 0.0);

        virtual float getBitFactor() {return 1.0; }

    private:

        static inline double _log2 (double x) { return log10(x)/0.301029996; }
};
#endif
