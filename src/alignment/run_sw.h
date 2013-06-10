#ifndef RUN_SW_H
#define RUN_SW_H

#ifdef OPENMP
#include <omp.h>
#endif

#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN

#include <string>
#include <cstdlib>
#include <sys/stat.h>
#include <malloc.h>
#include <xmmintrin.h>
#include <sstream>
#include <list>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include "smith_waterman_sse2.h"
#include "scoring.h"

typedef struct {
      size_t length;
      unsigned char* sequence;
} sequence_t;

typedef struct {
    std::string name;
    int score;
    float qcov;
    float dbcov;
    float eval;
} hit_t;

void printUsage();

void parseArgs(int argc, char** argv, std::string* ffindexSWBase, std::string* ffindexPrefBase, std::string* mFile, std::string* ffindexDbBase, std::string* ffindexDbNewBase, double* evalThr, double* covThr);

// compares two hits based on their scores (for sorting)
bool compareHits (hit_t first, hit_t second);

// inits ffindex for reading
void initFFIndexRead(std::string ffindexBase, char** data, ffindex_index_t** index);

// inits ffindex for writing
void initFFIndexWrite(std::string ffindexBase, FILE** dataFile, FILE** indexFile);

// converts a sequence string into unsigned char array, removes all spaces and stores the sequence and its length in a sequence_t structure
sequence_t seq2sequence_t(char* seq, size_t len, int* aa2int);

// calculates the Smith Waterman scores for a query sequence and the corresponding list of DB sequences stemming from the prefiltering
std::list<hit_t>* getSWScoresForSequence(char* querySeqData, int querySeqLen, char* prefList, char* dbData, ffindex_index_t* dbIndex, void* workspace, int* aa2int, int** scMatrix, unsigned char bias, int dbSize);

// run the parallelized SW scores calculation
void runSWParallel(std::string ffindexPrefBase, std::string ffindexDbBase, std::string ffindexSWBase, int* aa2int, int** scMatrix, double evalThr, double covThr);

/////////////////////////////////////////////////////////////////////////////////////
// fast 2^x
// ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
// Relative deviation < 4.6E-6  (< 2.3E-7 with 5'th order polynomial)
// Speed: 2.1E-8s (2.3E-8s) per call! (exp(): 8.5E-8, pow(): 1.7E-7)
// Internal representation of float number according to IEEE 754: 
//   1bit sign, 8 bits exponent, 23 bits mantissa: seee eeee emmm mmmm mmmm mmmm mmmm mmmm
//                                    0x4b400000 = 0100 1011 0100 0000 0000 0000 0000 0000 
//   In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127)
/////////////////////////////////////////////////////////////////////////////////////
inline float fpow2(float x)
{
    if (x>FLT_MAX_EXP) return FLT_MAX;
        if (x<FLT_MIN_EXP) return 0.0f;
    int *px = (int*)(&x);                 // store address of float as pointer to long int
    float tx = (x-0.5f) + (3<<22);        // temporary value for truncation: x-0.5 is added to a large integer (3<<22),
                                          // 3<<22 = (1.1bin)*2^23 = (1.1bin)*2^(150-127), 
                                          // which, in internal bits, is written 0x4b400000 (since 10010110bin = 150)
    int lx = *((int*)&tx) - 0x4b400000;   // integer value of x 
    float dx = x-(float)(lx);             // float remainder of x
    x = 1.0f + dx*(0.693019f             // polynomial approximation of 2^x for x in the range [0, 1]
             + dx*(0.241404f             // Gives relative deviation < 4.6E-6 
             + dx*(0.0520749f            // Speed: 2.1E-8s
             + dx* 0.0134929f )));
    *px += (lx<<23);                      // add integer power of 2 to exponent
    return x;
}


int main(int argc, char **argv);

#endif
