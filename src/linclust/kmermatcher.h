#ifndef MMSEQS_KMERMATCHER_H
#define MMSEQS_KMERMATCHER_H
#include "DBReader.h"
#include "DBWriter.h"
#include "Parameters.h"
#include "BaseMatrix.h"

#include <queue>
#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif

using AssignGroupMask = uint32_t;

struct AssignGroupFeature {
    static constexpr AssignGroupMask Default        = 0;
    static constexpr AssignGroupMask AdjacentSeq = 1u << 0;
    static constexpr AssignGroupMask CountTable  = 1u << 1;
};

enum class ComputationPhase {
    Main,
    SetupCountTable
};

inline bool hasFeature(AssignGroupMask mask, AssignGroupMask feature) {
    return (mask & feature) != 0;
}

struct SequencePosition{
    size_t kmer;
    unsigned int pos;
    unsigned short score;
    static bool compareByScore(const SequencePosition &first, const SequencePosition &second){
        if(first.score < second.score)
            return true;
        if(second.score < first.score)
            return false;
        if(first.kmer < second.kmer)
            return true;
        if(second.kmer < first.kmer)
            return false;
        if(first.pos < second.pos)
            return true;
        if(second.pos < first.pos)
            return false;
        return false;
    }
    static bool compareByScoreReverse(const SequencePosition &first, const SequencePosition &second){
        if(first.score < second.score)
            return true;
        if(second.score < first.score)
            return false;

        size_t firstKmer  = BIT_SET(first.kmer, 63);
        size_t secondKmer = BIT_SET(second.kmer, 63);
        if(firstKmer < secondKmer)
            return true;
        if(secondKmer < firstKmer)
            return false;
        if(first.pos < second.pos)
            return true;
        if(second.pos < first.pos)
            return false;
        return false;
    }
};

template <bool Include>
struct AdjacentData {
    unsigned char data[6];
    void set(int i, unsigned char v) { data[i] = v; }
    unsigned char get(int i) const { return data[i]; }
};

template <>
struct __attribute__((__packed__)) AdjacentData<false> {
    void set(int, unsigned char) { }
    unsigned char get(int) const { return '\0'; }
};

// [ isWholeSequence | k-mer | adjKey | maxSeqLen - seqLen | id | pos ] in integer order; a whole sequence holds its hash from adjKeyShift up
struct KmerPositionBit {
    static const unsigned int WHOLE_SEQUENCE_SHIFT = 127;
    static const unsigned int RESIDUE_BITS = 5;
    static unsigned int kmerShift;
    static unsigned int adjKeyShift;
    static unsigned int seqLenShift;
    static unsigned int idShift;
    static size_t maxSeqLen;
    __extension__ typedef unsigned __int128 uint128;
    uint128 value;

    void set(size_t kmer, unsigned int adjKey, size_t seqLen, DBKeyType id, size_t pos, bool isWholeSequence) {
        value = (static_cast<uint128>(isWholeSequence) << WHOLE_SEQUENCE_SHIFT)
                | (static_cast<uint128>(id) << idShift);
        // --createdb-mode 3 records have no length and no position field
        if (adjKeyShift > seqLenShift) {
            value |= static_cast<uint128>(maxSeqLen - seqLen) << seqLenShift;
        }
        if (isWholeSequence) {
            value |= static_cast<uint128>(kmer) << adjKeyShift;
        } else {
            value |= (static_cast<uint128>(kmer) << kmerShift) | (static_cast<uint128>(adjKey) << adjKeyShift);
            if (idShift > 0) {
                value |= pos;
            }
        }
    }

    bool isWholeSequence() const { return (value >> WHOLE_SEQUENCE_SHIFT) != 0; }
    uint128 getKmer() const { return value >> (isWholeSequence() ? adjKeyShift : kmerShift); }
    unsigned int getAdjKey() const { return getBits(adjKeyShift, kmerShift); }
    static unsigned char getAdjacentSeq(unsigned int adjKey, int index) { return (adjKey >> (RESIDUE_BITS * index)) & ((1u << RESIDUE_BITS) - 1); }
    unsigned char getAdjacentSeq(int index) const { return getAdjacentSeq(getAdjKey(), index); }
    size_t getSeqLen() const { return maxSeqLen - getBits(seqLenShift, adjKeyShift); }
    DBKeyType getId() const { return getBits(idShift, seqLenShift); }
    size_t getPos() const { return getBits(0, idShift); }

    bool operator<(const KmerPositionBit &other) const { return value < other.value; }

private:
    size_t getBits(unsigned int shift, unsigned int nextShift) const {
        return (value >> shift) & ((static_cast<uint128>(1) << (nextShift - shift)) - 1);
    }
};

template <typename T, bool includeAdjacency = false, typename L = T>
struct __attribute__((__packed__)) KmerPosition {
    size_t kmer;
    DBKeyType id;
    T pos;
    L seqLen;
    AdjacentData<includeAdjacency> adj;

    void setAdjacentSeq(int index, unsigned char val) {
        adj.set(index, val);
    }
    unsigned char getAdjacentSeq(int index) {
        return adj.get(index);
    }

    bool isWholeSequence() const { return false; }
    size_t getKmer() const { return kmer; }
    size_t getAdjKey() const {
        size_t adjKey = 0;
        for (int i = 5; i >= 0; i--) {
            adjKey = (adjKey << KmerPositionBit::RESIDUE_BITS) | adj.get(i);
        }
        return adjKey;
    }
    size_t getSeqLen() const { return seqLen; }
    DBKeyType getId() const { return id; }
    size_t getPos() const { return pos; }

    static bool compareKmer(const KmerPosition &first, const KmerPosition &second) {
        return first.kmer < second.kmer;
    }

    static bool compareKmerAndAdjKeyAndSeqLenAndIdAndPos(const KmerPosition &first, const KmerPosition &second) {
        if(first.kmer < second.kmer) return true;
        if(second.kmer < first.kmer) return false;
        if(first.getAdjKey() < second.getAdjKey()) return true;
        if(second.getAdjKey() < first.getAdjKey()) return false;
        if(first.seqLen > second.seqLen) return true;
        if(second.seqLen > first.seqLen) return false;
        if(first.id < second.id) return true;
        if(second.id < first.id) return false;
        if(first.pos < second.pos) return true;
        if(second.pos < first.pos) return false;
        return false;
    }

    static bool compareRepSequenceAndIdAndPos(
            const KmerPosition &first, const KmerPosition &second) {
        if(first.kmer < second.kmer) return true;
        if(second.kmer < first.kmer) return false;
        if(first.id != DB_KEY_INVALID && second.id != DB_KEY_INVALID) {
            L len1 = first.seqLen;
            L len2 = second.seqLen;
            if(len1 > len2) return true;
            if(len2 > len1) return false;
        }
        if(first.id < second.id) return true;
        if(second.id < first.id) return false;
        if(first.pos < second.pos) return true;
        if(second.pos < first.pos) return false;
        return false;
    }

    static bool compareRepSequenceAndIdAndPosReverse(
            const KmerPosition &first, const KmerPosition &second) {
        size_t firstKmer  = BIT_SET(first.kmer, 63);
        size_t secondKmer = BIT_SET(second.kmer, 63);
        if(firstKmer < secondKmer) return true;
        if(secondKmer < firstKmer) return false;
        if(first.id != DB_KEY_INVALID && second.id != DB_KEY_INVALID) {
            L len1 = first.seqLen;
            L len2 = second.seqLen;
            if(len1 > len2) return true;
            if(len2 > len1) return false;
        }
        if(first.id < second.id) return true;
        if(second.id < first.id) return false;
        if(first.pos < second.pos) return true;
        if(second.pos < first.pos) return false;
        return false;
    }

    static bool compareRepSequenceAndIdAndDiagReverse(
            const KmerPosition &first, const KmerPosition &second) {
        size_t firstKmer  = BIT_SET(first.kmer, 63);
        size_t secondKmer = BIT_SET(second.kmer, 63);
        if(firstKmer < secondKmer) return true;
        if(secondKmer < firstKmer) return false;
        if(first.id < second.id) return true;
        if(second.id < first.id) return false;
        if(first.pos < second.pos) return true;
        if(second.pos < first.pos) return false;
        return false;
    }

    static bool compareRepSequenceAndIdAndDiag(
            const KmerPosition &first, const KmerPosition &second) {
        if(first.kmer < second.kmer) return true;
        if(second.kmer < first.kmer) return false;
        if(first.id < second.id) return true;
        if(second.id < first.id) return false;
        if(first.pos < second.pos) return true;
        if(second.pos < first.pos) return false;
        return false;
    }
};

struct KmerGroupRepSeq {
    DBKeyType repSeqId;
    size_t groupStart;

    static bool compareRepSequence(const KmerGroupRepSeq &first, const KmerGroupRepSeq &second) {
        return first.repSeqId < second.repSeqId;
    }
};

struct __attribute__((__packed__)) KmerEntry {
    DBKeyType seqId;
    short diagonal;
    unsigned char score;
    void setReverse(bool ){
        ;
    }
    unsigned char getRev(){
        return 0;
    }
};

struct __attribute__((__packed__)) KmerEntryRev {
    DBKeyType seqId;
    short diagonal;
    unsigned char score;
    unsigned char rev;
    void setReverse(bool rev){
        this->rev = rev;
    }
    unsigned char getRev(){
        return this->rev;
    }
};

struct FileKmerPosition {
    DBKeyType repSeq;
    DBKeyType id;
    short pos;
    unsigned char score;
    unsigned int file;
    char reverse;
    FileKmerPosition(){}
    FileKmerPosition(DBKeyType repSeq, DBKeyType id,short pos, unsigned char score, unsigned int file):
            repSeq(repSeq), id(id), pos(pos), score(score), file(file), reverse(0) {}
    FileKmerPosition(DBKeyType repSeq, DBKeyType id,short pos, unsigned char score, char reverse, unsigned int file):
            repSeq(repSeq), id(id), pos(pos), score(score), file(file), reverse(reverse) {}
};

class CompareResultBySeqId {
public:
    bool operator() (FileKmerPosition & first, FileKmerPosition & second) const {
        //return (first.eval < second.eval);
        if(first.repSeq > second.repSeq )
            return true;
        if(second.repSeq > first.repSeq )
            return false;
        if(first.id > second.id )
            return true;
        if(second.id > first.id )
            return false;
        if(first.pos > second.pos )
            return true;
        if(second.pos > first.pos )
            return false;
        return false;
    }
};

template <int TYPE, typename T, bool includeAdjacency = false>
void mergeKmerFilesAndOutput(DBWriter & dbw, std::vector<std::string> tmpFiles, std::vector<char> &repSequence, int numThreads = 1, int maxIter = 1);

typedef std::priority_queue<FileKmerPosition, std::vector<FileKmerPosition>, CompareResultBySeqId> KmerPositionQueue;

template <int TYPE, typename T>
size_t queueNextEntry(KmerPositionQueue &queue, int file, size_t offsetPos, T *entries, size_t entrySize);

void setKmerLengthAndAlphabet(Parameters &parameters, size_t aaDbSize, int seqType);

template <int TYPE, typename T, typename posType, bool includeAdjacency, typename L>
void writeKmersToDisk(std::string tmpFile, KmerPosition<posType, includeAdjacency, L> *kmers, size_t totalKmers, int numThreads = 1, std::vector<size_t> *threadQueryOffsets = NULL, int iteration = 0);

template <int TYPE, typename T, bool includeAdjacency = false, typename L = T>
void writeKmerMatcherResult(DBWriter & dbw, KmerPosition<T, includeAdjacency, L> *hashSeqPair, size_t totalKmers,
                            std::vector<char> &repSequence, size_t threads);


template <typename T, bool includeAdjacency = false, typename L = T>
KmerPosition<T, includeAdjacency, L> * doComputation(size_t totalKmers, size_t split, size_t splits, std::string splitFile,
                                DBReader<DBKeyType> & seqDbr, Parameters & par, BaseMatrix  * subMat,
                                size_t KMER_SIZE, size_t chooseTopKmer, float chooseTopKmerScale = 0.0);

template <typename T, bool includeAdjacency = false, typename L = T>
KmerPosition<T, includeAdjacency, L> *initKmerPositionMemory(size_t size);

template <int TYPE, typename T, bool includeAdjacency = false, typename L = T, typename KmerPositionType>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPositionType * kmerArray, size_t kmerArraySize, DBReader<DBKeyType> &seqDbr,
                                                 Parameters & par, BaseMatrix * subMat, bool hashWholeSequence,
                                                 size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);


void maskSequence(int maskMode, int maskLowerCase,
                  Sequence &seq, int maskLetter, ProbabilityMatrix * probMatrix);

template <typename T, bool includeAdjacency = false, typename L = T>
size_t computeMemoryNeededLinearfilter(size_t totalKmer);

template <typename T, bool includeAdjacency = false, typename L = T>
std::vector<std::pair<size_t, size_t>> setupKmerSplits(Parameters &par, BaseMatrix * subMat, DBReader<DBKeyType> &seqDbr, size_t totalKmers, size_t splits,
                                                       size_t kmerMemory = sizeof(KmerPosition<T, includeAdjacency, L>));
size_t computeKmerCount(DBReader<DBKeyType> &reader, size_t KMER_SIZE, size_t chooseTopKmer,
                        float chooseTopKmerScale = 0.0);

void setLinearFilterDefault(Parameters *p);

size_t computeMemoryNeededLinearfilter(size_t totalKmer);


#undef SIZE_T_MAX


#endif
