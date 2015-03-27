#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H
#include "Util.h"
#include <float.h>
#include "Debug.h"
#include <iostream>
#include <cmath>
#include <cstdio>
#include <stdint.h>

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

    // substitution matrix for pseudocounts
    float** subMatrixPseudoCounts;

    // joint probability matrix
    double** probMatrix;

    // background probabilities of the amino acids
    double* pBack;

    // print the substitution matrix
    static void print(short** matrix, char* int2aa, int size);

    static void print(double** matrix, char* int2aa, int size);

    // generate the substitution matrix given the probability matrix, background probabilities and the alphabet size
    static void generateSubMatrix(double ** probMatrix, double ** subMatrix, float ** subMatrixPseudoCounts, int size, double bitFactor = 1.0, double scoringBias = 0.0);

    // generate a short data type substitution matrix
    static void generateSubMatrix(double ** probMatrix, float ** subMatrixPseudoCounts, short ** subMatrix, int size, double bitFactor = 1.0, double scoringBias = 0.0);

    virtual double getBackgroundProb(size_t aa_index) {
        Debug(Debug::ERROR) << "getBackground is not Impl. for this type of Matrix \n";
        EXIT(EXIT_FAILURE);
        return 0.0;
    };

    virtual float getBitFactor() {return 1.0; }

    static inline double _log2 (double x) { return log10(x)/0.301029996; }

    static inline float
    flog2(float x)
    {
        if (x <= 0)
            return -128;
        int *px = (int*) (&x);        // store address of float as pointer to long int
        float e = (float) (((*px & 0x7F800000) >> 23) - 0x7f); // shift right by 23 bits and subtract 127 = 0x7f => exponent
        *px = ((*px & 0x007FFFFF) | 0x3f800000);  // set exponent to 127 (i.e., 0)
        x -= 1.0;         // and calculate x-1.0
        x *= (1.441740
                + x * (-0.7077702 + x * (0.4123442 + x * (-0.1903190 + x * 0.0440047)))); // 5'th order polynomial approx. of log(1+x)
        return x + e;
    }

    static inline double fpow2(float x) {
        if (x>=FLT_MAX_EXP) return FLT_MAX;
        if (x<=FLT_MIN_EXP) return 0.0f;

        int *px = (int*) (&x);        // store address of float as pointer to long int
        float tx = (x - 0.5f) + (3 << 22); // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
        // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127),
        // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)
        int lx = *((int*) &tx) - 0x4b400000;   // integer value of x
        float dx = x - (float) (lx);             // float remainder of x
//   x = 1.0f + dx*(0.69606564f           // cubic apporoximation of 2^x for x in the range [0, 1]
//            + dx*(0.22449433f           // Gives relative deviation < 1.5E-4
//            + dx*(0.07944023f)));       // Speed: 1.9E-8s
        x = 1.0f + dx * (0.693019f // polynomial approximation of 2^x for x in the range [0, 1]
                + dx * (0.241404f             // Gives relative deviation < 4.6E-6
                + dx * (0.0520749f            // Speed: 2.1E-8s
                + dx * 0.0134929f)));
//   x = 1.0f + dx*(0.693153f             // polynomial apporoximation of 2^x for x in the range [0, 1]
//            + dx*(0.240153f             // Gives relative deviation < 2.3E-7
//            + dx*(0.0558282f            // Speed: 2.3E-8s
//            + dx*(0.00898898f
//            + dx* 0.00187682f ))));
        *px += (lx << 23);                      // add integer power of 2 to exponent
        return x;
    }

private:

};
#endif
