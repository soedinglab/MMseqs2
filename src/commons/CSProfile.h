//
// Created by mad on 11/29/16.
//

#ifndef MMSEQS_CSBLASTPROFILES_H
#define MMSEQS_CSBLASTPROFILES_H

#include <stdlib.h>
#include <simd/simd.h>
#include "LibraryReader.h"
#include "Debug.h"
#include "Util.h"
#include "ProfileStates.h"
#include "Sequence.h"
class ContextLibrary{

public:

    size_t wlen_;                            // size of context window.
    size_t center;                           // center of window
    size_t libSize;                          // size of library (e.g. 4000 K4000.crf)
    LibraryReader reader;
    std::vector<std::string> names;
    std::vector<float>  bias_weight;
    // keep ContextLibrary probs
    float *** context_weights; // k * states * aa
    float ** pc_weights;
    float ** pc;

    static ContextLibrary* getContextLibraryInstance()
    {
#pragma omp critical
        {
            if (contextLibrary==NULL) {
                contextLibrary = new ContextLibrary();
            }
        }
        return contextLibrary;
    }

private:
    static ContextLibrary * contextLibrary;

    ContextLibrary();
    ~ContextLibrary();

    // read library file
    void read(std::string &libStr);

    // read context profiles from library
    void readContextProfile(std::stringstream &in, LibraryReader &reader,
                            float ** context_weight, float * pc_weight, float * pc);
};

class CSProfile {
    ContextLibrary * ctxLib;
    float * profile;
    float * pp;
    float * maximums;
    float * sums;
public:

    CSProfile(size_t maxSeqLen) {
        ctxLib = ContextLibrary::getContextLibraryInstance();
        this->profile = (float * )mem_align(ALIGN_FLOAT, (Sequence::PROFILE_AA_SIZE + 4) * maxSeqLen * sizeof(float));
        int segmentSize = (maxSeqLen+VECSIZE_FLOAT-1)/VECSIZE_FLOAT;
        this->pp = (float * ) mem_align(ALIGN_FLOAT, 4000 * segmentSize * VECSIZE_FLOAT * sizeof(float));
        this->maximums = (float * ) mem_align(ALIGN_FLOAT,  segmentSize * VECSIZE_FLOAT * sizeof(float));
        this->sums = (float * ) mem_align(ALIGN_FLOAT,  segmentSize * VECSIZE_FLOAT * sizeof(float));
    };

    ~CSProfile(){
        free(profile);
        free(pp);
        free(maximums);
        free(sums);
    }

    float computeSeqContextScore(float ** context_weights,
                              const unsigned char * seq, const int L,
                              size_t idx, size_t center);

    float computeProfileContextScore(float ** context_weights,
                                 const float * counts, const int L,
                                 size_t idx, size_t center);

    float * computeProfileCs(int seqLen, float * count, float * Neff_M, float pca, float pcb);
    float * computeSequenceCs(unsigned char * numSeq, int seqLen, float tau);
private:
    template<int type>
    float * computeProfile(unsigned char * numSeq, int seqLen, float * count, float * Neff_M, float pTau, float pca, float pcb);
};


#endif //MMSEQS_CSBLASTPROFILES_H
