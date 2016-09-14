#ifndef SEQUENCE_H
#define SEQUENCE_H
#define READ_BUFFER_SIZE 20971520

// Written by Martin Steinegger & Maria Hauser mhauser@genzentrum.lmu.de, Martin Steinegger martin.steinegger@mpibpc.mpg.de
// 
// Represents a database sequence object, holds its representation in the int array form.
//

#include <cstdint>
#include <cstddef>
#include <utility>

struct ScoreMatrix;

const int8_t seed_4[]        = {1, 1, 1, 1};
const int8_t seed_4_spaced[] = {1, 1, 1, 0, 1};
const int8_t seed_5[]        = {1, 1, 1, 1, 1};
//const char seed_5_spaced[] = {1, 1, 1, 0, 1, 1};
//const char seed_5_spaced[] = {1, 1, 0, 1, 0, 1, 1};// better than 111011
const int8_t seed_5_spaced[] = {1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1}; // just 0.001 %worse ROC5 but way faster
const int8_t seed_6[]        = {1, 1, 1, 1, 1, 1};
//const char seed_6_spaced[] = {1, 1, 1, 0, 1, 1, 0, 1};
const int8_t seed_6_spaced[] = {1, 1, 0, 1, 0, 1, 0, 0, 1, 1}; // better than 11101101
const int8_t seed_7[]        = {1, 1, 1, 1, 1, 1, 1};
//const char seed_7_spaced[] = {1, 1, 1, 1, 0, 1, 0, 1, 0, 1};
const int8_t seed_7_spaced[] = {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1};


const int8_t seed_9[]         = {1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_9_spaced[]  = {1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_10[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_10_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_11[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_11_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_12[]        = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int8_t seed_12_spaced[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
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
        Sequence(size_t maxLen, int *aa2int, char *int2aa, int seqType,
                 const unsigned int kmerSize, const bool spaced, const bool aaBiasCorrection);
        ~Sequence();

        // Map char -> int
        void mapSequence(size_t id, unsigned int dbKey, const char *seq);

        // map sequence from SequenceLookup
        void mapSequence(size_t id, unsigned int dbKey, std::pair<const unsigned char *, const unsigned int> data);

        // map profile HMM, *data points to start position of Profile
        void mapProfile(const char *data);

        // checks if there is still a k-mer left
        bool hasNextKmer();

        // returns next k-mer
        const int* nextKmer();

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
    
        // length of sequence
        int L;

        // each amino acid coded as integer
        int * int_sequence;  

        // Contains profile information
        short           * profile_score;
        unsigned int    * profile_index;
        size_t profile_row_size; // (PROFILE_AA_SIZE / SIMD_SIZE) + 1 * SIMD_SIZE

        static const size_t PROFILE_AA_SIZE = 20;
        ScoreMatrix ** profile_matrix;
        // Memory layout of this profile is qL * AA
        //   Query lenght
        // A  -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
        // C  -1  -4   2   5  -3  -2   0  -3   1  -3  -2   0  -1   2   0   0  -1  -3  -4  -2
        // ...
        // Y -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3

        int8_t * profile_for_alignment;
    
        int  * aa2int; // ref to mapping from aa -> int
        char * int2aa; // ref mapping from int -> aa

        std::pair<const char *, unsigned int> getSpacedPattern(bool spaced, unsigned int kmerSize);

        const int *getKmerPositons();

        void printProfile();

        int8_t const * getAlignmentProfile()const;

        int getSequenceType()const;

        unsigned int getEffectiveKmerSize();


private:
        void mapProteinSequence(const char *seq);
        void mapNucleotideSequence(const char *seq);
        size_t id;
        unsigned int dbKey;

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

        // contains sequence positions for current kmer
        int *kmerPos;

        // bias correction in profiles
        bool aaBiasCorrection;

        // spaced pattern
        bool spaced;
};
#endif
