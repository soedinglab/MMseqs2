#ifndef SEQUENCE_H
#define SEQUENCE_H
// Written by Martin Steinegger & Maria Hauser mhauser@genzentrum.lmu.de, Martin Steinegger martin.steinegger@snu.ac.kr
//
// Represents a database sequence object, holds its representation in the int array form.
//

#include "MathUtil.h"
#include "BaseMatrix.h"
#include "Parameters.h"
#include "ScoreMatrix.h"

#include <cstdint>
#include <cstddef>
#include <utility>
#include <simd/simd.h>

const int8_t seed_4[]        = {1, 1, 1, 1};
const int8_t spaced_seed_4[] = {1, 1, 1, 0, 1};
const int8_t seed_5[]        = {1, 1, 1, 1, 1};
const int8_t spaced_seed_5[]  = {1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1}; // just 0.001 %worse ROC5 but way faster
const int8_t seed_6[]         = {1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_6[]  = {1, 1, 0, 1, 0, 1, 0, 0, 1, 1}; // better than 11101101
const int8_t seed_7[]         = {1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_7[]  = {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1};
const int8_t seed_8[]         = {1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_8[]  = {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1};
const int8_t seed_9[]         = {1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_9[]  = {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1};
const int8_t seed_10[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_10[] = {1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1};
const int8_t seed_11[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_11[] = {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1};
const int8_t seed_12[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_12[] = {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1};
const int8_t seed_13[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_13[] = {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1};
const int8_t seed_14[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_14[] = {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1};
const int8_t seed_15[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_15[] = {1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1};
const int8_t seed_16[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_16[] = {1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1};
const int8_t seed_17[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_17[] = {1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1};
const int8_t seed_18[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_18[] = {1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1};
const int8_t seed_19[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_19[] = {1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1};
const int8_t seed_20[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_20[] = {1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1};
const int8_t seed_21[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_21[] = {1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1};
const int8_t seed_22[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_22[] = {1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1};
const int8_t seed_23[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_23[] = {1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1};
const int8_t seed_24[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_24[] = {1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1};
const int8_t seed_25[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_25[] = {1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1};
const int8_t seed_26[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_26[] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1};
const int8_t seed_27[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_27[] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1};
const int8_t seed_28[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_28[] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_29[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_29[] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_30[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t spaced_seed_30[] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1};




class Sequence {
public:
    Sequence(size_t maxLen, int seqType, const BaseMatrix *subMat,  const unsigned int kmerSize, const bool spaced, const bool aaBiasCorrection,
             bool shouldAddPC = true, const std::string& userSpacedKmerPattern = "");
    ~Sequence();

    // Map char -> int
    void mapSequence(size_t id, unsigned int dbKey, const char *seq, unsigned int seqLen);

    // map sequence from SequenceLookup
    void mapSequence(size_t id, unsigned int dbKey, std::pair<const unsigned char *, const unsigned int> data);

    // map profile HMM, *data points to start position of Profile
    void mapProfile(const char *profileData, unsigned int seqLen);

    // checks if there is still a k-mer left
    bool hasNextKmer() {
        return (((currItPos + 1) + this->spacedPatternSize) <= this->L);
    }

    // k-mer contains x, is only field aftter nextKmer
    inline bool kmerContainsX(){
        return kmerHasX != 0;
    }

    // returns next k-mer
    inline const unsigned char * nextKmer() {
        if (hasNextKmer() == false) {
            return 0;
        }

        currItPos++;
        const unsigned char *posToRead = numSequence + currItPos;
        unsigned char *currWindowPos = kmerWindow;
        switch (this->kmerSize){
            case 6:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                break;
            case 7:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                break;
            case 10:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                break;
            case 11:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                break;
            case 12:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                break;
            case 13:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                break;
            case 14:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                break;
            case 15:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                kmerWindow[14] = posToRead[aaPosInSpacedPattern[14]];

                break;
            case 16:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                kmerWindow[14] = posToRead[aaPosInSpacedPattern[14]];
                kmerWindow[15] = posToRead[aaPosInSpacedPattern[15]];
                break;
            case 17:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                kmerWindow[14] = posToRead[aaPosInSpacedPattern[14]];
                kmerWindow[15] = posToRead[aaPosInSpacedPattern[15]];
                kmerWindow[16] = posToRead[aaPosInSpacedPattern[16]];
                break;
            case 18:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                kmerWindow[14] = posToRead[aaPosInSpacedPattern[14]];
                kmerWindow[15] = posToRead[aaPosInSpacedPattern[15]];
                kmerWindow[16] = posToRead[aaPosInSpacedPattern[16]];
                kmerWindow[17] = posToRead[aaPosInSpacedPattern[17]];
                break;
            case 19:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                kmerWindow[14] = posToRead[aaPosInSpacedPattern[14]];
                kmerWindow[15] = posToRead[aaPosInSpacedPattern[15]];
                kmerWindow[16] = posToRead[aaPosInSpacedPattern[16]];
                kmerWindow[17] = posToRead[aaPosInSpacedPattern[17]];
                kmerWindow[18] = posToRead[aaPosInSpacedPattern[18]];
                break;
            case 20:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                kmerWindow[14] = posToRead[aaPosInSpacedPattern[14]];
                kmerWindow[15] = posToRead[aaPosInSpacedPattern[15]];
                kmerWindow[16] = posToRead[aaPosInSpacedPattern[16]];
                kmerWindow[17] = posToRead[aaPosInSpacedPattern[17]];
                kmerWindow[18] = posToRead[aaPosInSpacedPattern[18]];
                kmerWindow[19] = posToRead[aaPosInSpacedPattern[19]];
                break;
            case 21:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                kmerWindow[14] = posToRead[aaPosInSpacedPattern[14]];
                kmerWindow[15] = posToRead[aaPosInSpacedPattern[15]];
                kmerWindow[16] = posToRead[aaPosInSpacedPattern[16]];
                kmerWindow[17] = posToRead[aaPosInSpacedPattern[17]];
                kmerWindow[18] = posToRead[aaPosInSpacedPattern[18]];
                kmerWindow[19] = posToRead[aaPosInSpacedPattern[19]];
                kmerWindow[20] = posToRead[aaPosInSpacedPattern[20]];
                break;
            case 22:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                kmerWindow[14] = posToRead[aaPosInSpacedPattern[14]];
                kmerWindow[15] = posToRead[aaPosInSpacedPattern[15]];
                kmerWindow[16] = posToRead[aaPosInSpacedPattern[16]];
                kmerWindow[17] = posToRead[aaPosInSpacedPattern[17]];
                kmerWindow[18] = posToRead[aaPosInSpacedPattern[18]];
                kmerWindow[19] = posToRead[aaPosInSpacedPattern[19]];
                kmerWindow[20] = posToRead[aaPosInSpacedPattern[20]];
                kmerWindow[21] = posToRead[aaPosInSpacedPattern[21]];
                break;
            case 23:
                kmerWindow[0] = posToRead[aaPosInSpacedPattern[0]];
                kmerWindow[1] = posToRead[aaPosInSpacedPattern[1]];
                kmerWindow[2] = posToRead[aaPosInSpacedPattern[2]];
                kmerWindow[3] = posToRead[aaPosInSpacedPattern[3]];
                kmerWindow[4] = posToRead[aaPosInSpacedPattern[4]];
                kmerWindow[5] = posToRead[aaPosInSpacedPattern[5]];
                kmerWindow[6] = posToRead[aaPosInSpacedPattern[6]];
                kmerWindow[7] = posToRead[aaPosInSpacedPattern[7]];
                kmerWindow[8] = posToRead[aaPosInSpacedPattern[8]];
                kmerWindow[9] = posToRead[aaPosInSpacedPattern[9]];
                kmerWindow[10] = posToRead[aaPosInSpacedPattern[10]];
                kmerWindow[11] = posToRead[aaPosInSpacedPattern[11]];
                kmerWindow[12] = posToRead[aaPosInSpacedPattern[12]];
                kmerWindow[13] = posToRead[aaPosInSpacedPattern[13]];
                kmerWindow[14] = posToRead[aaPosInSpacedPattern[14]];
                kmerWindow[15] = posToRead[aaPosInSpacedPattern[15]];
                kmerWindow[16] = posToRead[aaPosInSpacedPattern[16]];
                kmerWindow[17] = posToRead[aaPosInSpacedPattern[17]];
                kmerWindow[18] = posToRead[aaPosInSpacedPattern[18]];
                kmerWindow[19] = posToRead[aaPosInSpacedPattern[19]];
                kmerWindow[20] = posToRead[aaPosInSpacedPattern[20]];
                kmerWindow[21] = posToRead[aaPosInSpacedPattern[21]];
                kmerWindow[22] = posToRead[aaPosInSpacedPattern[22]];

                break;
            default:
                for (unsigned int i = 0; i < this->kmerSize; i++) {
                    unsigned char pos = aaPosInSpacedPattern[i];
                    currWindowPos[0] = posToRead[pos];
                    currWindowPos++;
                }
                break;
        }
        kmerHasX = 0;

        const simd_int xChar = simdi8_set(subMat->aa2num[static_cast<int>('X')]);
        for(size_t i = 0; i < simdKmerRegisterCnt; i++){
            simd_int kmer = simdi_load((((simd_int *) kmerWindow) + i));
            kmerHasX |= static_cast<unsigned int>(simdi8_movemask(simdi8_eq(kmer, xChar)));
        }
        if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE)) {
            nextProfileKmer();
            for (unsigned int i = 0; i < this->kmerSize; i++) {
                kmerWindow[i] = 0;
            }
            return kmerWindow;
        }

        return (const unsigned char *)kmerWindow;
    }

    // resets the sequence position pointer to the start of the sequence
    void resetCurrPos() { currItPos = -1; }

    void print(); // for debugging

    static void extractProfileSequence(const char* data, const BaseMatrix &submat, std::string &result);
    static void extractProfileConsensus(const char* data, const BaseMatrix &submat, std::string &result);

    int getId() const { return id; }

    int getCurrentPosition() { return currItPos; }

    unsigned int getDbKey() { return dbKey; }

    int getSeqType() { return seqType; }

    size_t getMaxLen() { return maxLen; }
    unsigned int getKmerSize(){ return kmerSize; }
    bool isSpaced() { return spaced; }

    // reverse the sequence for the match statistics calculation
    void reverse();

    // submat
    BaseMatrix *subMat;

    // length of sequence
    int L;

    // each amino acid coded as integer
    unsigned char *numSequence;

    // each consensus amino acid as integer (PROFILE ONLY)
    unsigned char *numConsensusSequence;

    // Contains profile information
    short           *profile_score;
    unsigned int    *profile_index;
    float           *neffM;
    uint8_t         *gDel;
    uint8_t         *gIns;
    float           *pseudocountsWeight;
    // const size_t PROFILE_ROW_SIZE = (((size_t) PROFILE_AA_SIZE / (VECSIZE_INT * 4)) + 1) * (VECSIZE_INT * 4);
    size_t profile_row_size;
    static const size_t PROFILE_AA_SIZE = 20;
    static const size_t PROFILE_CONSENSUS = 21;     // new
    static const size_t PROFILE_NEFF = 22;          // new
    static const size_t PROFILE_GAP_DEL = 23;       // new
    static const size_t PROFILE_GAP_INS = 24;       // new
    // 20 AA, 1 query, 1 consensus, 1 Neff M, 2 gap penalties
    static const size_t PROFILE_READIN_SIZE = 25;
    ScoreMatrix **profile_matrix;
    // Memory layout of this profile is qL * AA
    //   Query length
    // A  -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
    // C  -1  -4   2   5  -3  -2   0  -3   1  -3  -2   0  -1   2   0   0  -1  -3  -4  -2
    // ...
    // Y -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
    int8_t *profile_for_alignment;

    std::pair<const char *, unsigned int> getSpacedPattern(bool spaced, unsigned int kmerSize);

    std::pair<const char *, unsigned int> parseSpacedPattern(unsigned int kmerSize, bool spaced, const std::string& spacedKmerPattern);

    const unsigned char *getAAPosInSpacedPattern() { return aaPosInSpacedPattern; }

    void printPSSM();

    int8_t const* getAlignmentProfile() const {
        return profile_for_alignment;
    }

    void printProfile() const;

    int getSequenceType() const {
        return seqType;
    }

    unsigned int getEffectiveKmerSize() {
        return spacedPatternSize;
    }

    static unsigned char scoreMask(float prob) {
        unsigned char charProb = MathUtil::convertFloatToChar(prob);
        // avoid 0
        return charProb + 1;
    }

    static float scoreUnmask(unsigned char score) {
        float prob = MathUtil::convertCharToFloat(score-1);
        return prob;
    }

    static float probaToBitScore(double proba, double pBack) {
        // No score bias when profile proba stored in file
        return MathUtil::flog2(proba / pBack);
    }

    static float scoreToProba(short score, double pBack, double bitFactor, double scoreBias) {
        if (score == -127) {
            return 0.0;
        }
        double dblScore = static_cast<double>(score);
        // No score bias when profile probability stored in file
        return MathUtil::fpow2((float)(dblScore - scoreBias) / bitFactor) * pBack;
    }

    const char *getSeqData() {
        return seqData;
    }

    const std::string &getUserSpacedKmerPattern() const {
        return userSpacedKmerPattern;
    }

private:
    void mapSequence(const char *seq, unsigned int dataLen);
    // read next kmer profile in profile_matrix
    void nextProfileKmer();

    size_t id;
    unsigned int dbKey;
    const char *seqData;

    // current iterator position
    int currItPos;

    // DBTYPE_AMINO_ACIDS or DBTYPE_NUCLEOTIDES
    int seqType;

    // maximum possible length of sequence
    size_t maxLen;

    // size of Pattern
    int spacedPatternSize;

    // contains spaced pattern e.g. 1 1 1 1 0 1 0 1 0 1
    const char *spacedPattern;

    // kmer Size
    unsigned int kmerSize;

    // simd kmer size
    unsigned int simdKmerRegisterCnt;

    // sequence window will be filled by newxtKmer (needed for spaced patterns)
    unsigned char *kmerWindow;

    // set if kmer contains X
    unsigned int kmerHasX;

    // stores position of residues in sequence
    unsigned char *aaPosInSpacedPattern;

    // buffer for background null probability for global aa bias correction
    float *pNullBuffer;

    // bias correction in profiles
    bool aaBiasCorrection;

    // spaced pattern
    bool spaced;

    // should add pseudo-counts when loading the profile?
    bool shouldAddPC;

    // user kmer pattern
    std::string userSpacedKmerPattern;
};
#endif
