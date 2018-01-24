#ifndef SEQUENCE_H
#define SEQUENCE_H
#define READ_BUFFER_SIZE 20971520

// Written by Martin Steinegger & Maria Hauser mhauser@genzentrum.lmu.de, Martin Steinegger martin.steinegger@mpibpc.mpg.de
//
// Represents a database sequence object, holds its representation in the int array form.
//
#include "Debug.h"
#include "MathUtil.h"
#include "BaseMatrix.h"
#include <cstdint>
#include <cstddef>
#include <utility>

struct ScoreMatrix;

const int8_t seed_4[]        = {1, 1, 1, 1};
const int8_t seed_4_spaced[] = {1, 1, 1, 0, 1};
const int8_t seed_5[]        = {1, 1, 1, 1, 1};
//const char seed_5_spaced[] = {1, 1, 1, 0, 1, 1};
//const char seed_5_spaced[] = {1, 1, 0, 1, 0, 1, 1};// better than 111011
const int8_t seed_5_spaced[]  = {1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1}; // just 0.001 %worse ROC5 but way faster
const int8_t seed_6[]         = {1, 1, 1, 1, 1, 1};
//const char seed_6_spaced[] = {1, 1, 1, 0, 1, 1, 0, 1};
const int8_t seed_6_spaced[]  = {1, 1, 0, 1, 0, 1, 0, 0, 1, 1}; // better than 11101101
const int8_t seed_7[]         = {1, 1, 1, 1, 1, 1, 1};
//const char seed_7_spaced[] = {1, 1, 1, 1, 0, 1, 0, 1, 0, 1};
const int8_t seed_7_spaced[]  = {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1};
const int8_t seed_9[]         = {1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_9_spaced[]  = {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1};
const int8_t seed_10[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_10_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_11[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_11_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_12[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_12_spaced[] = {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1};
const int8_t seed_13[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_13_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_14[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_14_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_15[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_15_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_16[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_16_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_17[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_17_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};


class Sequence
{
public:
    Sequence(size_t maxLen, int seqType, const BaseMatrix *subMat,
             const unsigned int kmerSize, const bool spaced, const bool aaBiasCorrection);
    ~Sequence();

    // Map char -> int
    void mapSequence(size_t id, unsigned int dbKey, const char *seq);

    // map sequence from SequenceLookup
    void mapSequence(size_t id, unsigned int dbKey, std::pair<const unsigned char *, const unsigned int> data);

    // map profile HMM, *data points to start position of Profile
    void mapProfile(const char *sequence, bool mapScores);

    // mixture of library and profile prob
    void mapProfileState(const char *sequence);

    // map the profile state sequence
    void mapProfileStateSequence(const char *sequence);

    // checks if there is still a k-mer left
    bool hasNextKmer() {
            return (((currItPos + 1) + this->spacedPatternSize) <= this->L);
    }

    // returns next k-mer
    inline const int * nextKmer() {
            if (hasNextKmer()) {
                    currItPos++;
                    const int * posToRead = int_sequence + currItPos;
                    int * currWindowPos = kmerWindow;
                    switch(this->kmerSize){
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
                        default:
                            for(unsigned int i = 0; i < this->kmerSize; i++) {
                                unsigned char pos = aaPosInSpacedPattern[i];
                                currWindowPos[0] = posToRead[pos];
                                currWindowPos++;
                            }
                            break;
                    }

                    if(seqType == HMM_PROFILE||seqType==PROFILE_STATE_PROFILE) {
                            nextProfileKmer();
                            for(unsigned int i = 0; i < this->kmerSize; i++) {
                                    kmerWindow[i] = 0;
                            }
                            return kmerWindow;
                    }


                    return (const int *) kmerWindow;
            }
            return 0;
    }

    // resets the sequence position pointer to the start of the sequence
    void resetCurrPos() { currItPos = -1; }

    void print(); // for debugging

    int getId() { return id; }

    int getCurrentPosition() { return currItPos; }

    unsigned int getDbKey() { return dbKey; }

    int getSeqType() { return seqType; }

    size_t getMaxLen(){ return maxLen; }
    unsigned int getKmerSize(){ return kmerSize; }
    bool isSpaced(){ return spaced; }

    // reverse the sequence for the match statistics calculation
    void reverse();

    static const int AMINO_ACIDS = 0;
    static const int NUCLEOTIDES = 1;
    static const int HMM_PROFILE = 2;
    static const int PROFILE_STATE_SEQ = 3;
    static const int PROFILE_STATE_PROFILE = 4;


    // submat
    BaseMatrix * subMat;

    // length of sequence
    int L;

    // each amino acid coded as integer
    int * int_sequence;

    // each consensus amino acid as integer (PROFILE ONLY)
    int * int_consensus_sequence;

    // Contains profile information
    short           * profile_score;
    unsigned int    * profile_index;
    float           * profile;
    float           * neffM;
    float           * pseudocountsWeight;
    size_t profile_row_size; // (PROFILE_AA_SIZE / SIMD_SIZE) + 1 * SIMD_SIZE

    static const size_t PROFILE_AA_SIZE = 20;
    static const size_t PROFILE_READIN_SIZE = 23; // 20 AA, 1 query, 1 consensus, 2 for Neff M,
    ScoreMatrix ** profile_matrix;
    // Memory layout of this profile is qL * AA
    //   Query lenght
    // A  -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
    // C  -1  -4   2   5  -3  -2   0  -3   1  -3  -2   0  -1   2   0   0  -1  -3  -4  -2
    // ...
    // Y -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3

    int8_t * profile_for_alignment;

    std::pair<const char *, unsigned int> getSpacedPattern(bool spaced, unsigned int kmerSize);

    const unsigned char *getAAPosInSpacedPattern() {     return aaPosInSpacedPattern;  }

    void printPSSM();

    void printProfileStatePSSM();

    void printProfile();

    int8_t const * getAlignmentProfile()const;

    int getSequenceType()const;

    unsigned int getEffectiveKmerSize();

    static unsigned char scoreMask(float prob)
    {
        unsigned char charProb = MathUtil::convertFloatToChar(prob);
        // avoid 0
        return charProb + 1;
    }

    static float scoreUnmask(unsigned char score)
    {
        float prob = MathUtil::convertCharToFloat(score-1);
        return prob;
    }

    static float probaToBitScore(double proba, double pBack)
    {
        // No score bias when profile proba stored in file
        return MathUtil::flog2(proba / pBack);
    }


    static float scoreToProba(short score, double pBack, double bitFactor, double scoreBias)
    {
        if (score == -127)
            return 0.0;
        double dblScore = static_cast<double>(score);
        // No score bias when profile proba stored in file
        float prob = MathUtil::fpow2( (float) (dblScore - scoreBias) / bitFactor )*pBack;
        return prob;
    }


    const float *getProfile();

    const char * getSeqData(){
        return seqData;
    }

private:
    void mapSequence(const char *seq);
    size_t id;
    unsigned int dbKey;
    const char * seqData;

    // current iterator position
    int currItPos;

    // AMINO_ACIDS or NUCLEOTIDES
    int seqType;

    // maximum possible length of sequence
    size_t maxLen;

    // read next kmer profile in profile_matrix
    void nextProfileKmer();

    // size of Pattern
    int spacedPatternSize;

    // contains spaced pattern e.g. 1 1 1 1 0 1 0 1 0 1
    const char * spacedPattern;

    // kmer Size
    unsigned int kmerSize;

    // sequence window will be filled by newxtKmer (needed for spaced patterns)
    int * kmerWindow;

    // stores position of residues in sequence
    unsigned char *aaPosInSpacedPattern;

    // bias correction in profiles
    bool aaBiasCorrection;

    // spaced pattern
    bool spaced;
};
#endif

