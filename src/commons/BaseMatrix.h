#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H
#include "Util.h"
#include "Debug.h"
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

        virtual double getBackgroundProb(size_t aa_index) {
                Debug(Debug::ERROR) << "getBackground is not Impl. for this type of Matrix \n";
                EXIT(EXIT_FAILURE);
                return 0.0;
        };

        virtual float getBitFactor() {return 1.0; }
    
        static inline double _log2 (double x) { return log10(x)/0.301029996; }
    
    static inline float
    fastlog2 (float x)
    {
        union { float f; uint32_t i; } vx = { x };
        union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
        float y = vx.i;
        y *= 1.1920928955078125e-7f;
        
        return y - 124.22551499f
        - 1.498030302f * mx.f
        - 1.72587999f / (0.3520887068f + mx.f);
    }
    
        static inline double fastPow(double a, double b) {
            union {
                double d;
                int x[2];
            } u = { a };
            u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
            u.x[0] = 0;
            return u.d;
        }

    private:

};
#endif
