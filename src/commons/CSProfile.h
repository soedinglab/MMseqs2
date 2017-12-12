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
    int read(std::string &libStr);

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
        this->profile = (float * )mem_align(16, Sequence::PROFILE_AA_SIZE * maxSeqLen * sizeof(float));
        this->pp =  (float * ) mem_align(16, 4000 * maxSeqLen * sizeof(float)); //TODO how can I avoid th 4000?
        this->maximums = new float[maxSeqLen];
        this->sums =  new float[maxSeqLen];
    };

    ~CSProfile(){
        free(profile);
        free(pp);
        delete [] maximums;
        delete [] sums;
    }

    float computeContextScore(float ** context_weights,
                              const int * seq, const int L,
                              size_t idx, size_t center);

    float * computeProfile(Sequence * seq, float neff, float tau);

};


#endif //MMSEQS_CSBLASTPROFILES_H
