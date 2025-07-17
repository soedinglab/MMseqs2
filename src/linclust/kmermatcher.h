#ifndef MMSEQS_KMERMATCHER_H
#define MMSEQS_KMERMATCHER_H
#include "DBReader.h"
#include "DBWriter.h"
#include "Parameters.h"
#include "BaseMatrix.h"
#include "SequenceWeights.h"

#include <queue>

struct SequencePosition{
    unsigned short score;
    size_t kmer;
    unsigned int pos;
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
struct AdjacentSeqArray {
    void setAdjacentSeq(int index, const unsigned char val) {
        _adjacentSeq[index] = val;
    }
    unsigned char getAdjacentSeq(int index) {
        return _adjacentSeq[index];
    }

    unsigned char _adjacentSeq[6];
};

// save memory when adjacent sequence is unused
template <>
struct AdjacentSeqArray<false> {
    void setAdjacentSeq(const int, const unsigned char) {};
    unsigned char getAdjacentSeq(int) {
        Debug(Debug::ERROR) << "Invalid read attempt at adjacent sequence array";
        return '\0';
    }
};

template <typename T, bool IncludeAdjacentSeq = false>
struct __attribute__((__packed__))KmerPosition : public AdjacentSeqArray<IncludeAdjacentSeq> {
    size_t kmer;
    unsigned int id;
    T seqLen;
    T pos;

    static bool compareRepSequenceAndIdAndPos(const KmerPosition<T, IncludeAdjacentSeq> &first, const KmerPosition<T, IncludeAdjacentSeq> &second){
        if(first.kmer < second.kmer )
            return true;
        if(second.kmer < first.kmer )
            return false;
        if(first.seqLen > second.seqLen )
            return true;
        if(second.seqLen > first.seqLen )
            return false;
        if(first.id < second.id )
            return true;
        if(second.id < first.id )
            return false;
        if(first.pos < second.pos )
            return true;
        if(second.pos < first.pos )
            return false;
        return false;
    }

    static bool compareRepSequenceAndIdAndPosReverse(const KmerPosition<T, IncludeAdjacentSeq> &first, const KmerPosition<T, IncludeAdjacentSeq> &second){
        size_t firstKmer  = BIT_SET(first.kmer, 63);
        size_t secondKmer = BIT_SET(second.kmer, 63);
        if(firstKmer < secondKmer )
            return true;
        if(secondKmer < firstKmer )
            return false;
        if(first.seqLen > second.seqLen )
            return true;
        if(second.seqLen > first.seqLen )
            return false;
        if(first.id < second.id )
            return true;
        if(second.id < first.id )
            return false;
        if(first.pos < second.pos )
            return true;
        if(second.pos < first.pos )
            return false;
        return false;
    }

    static bool compareRepSequenceAndIdAndDiagReverse(const KmerPosition<T, IncludeAdjacentSeq> &first, const KmerPosition<T, IncludeAdjacentSeq> &second){
        size_t firstKmer  = BIT_SET(first.kmer, 63);
        size_t secondKmer = BIT_SET(second.kmer, 63);
        if(firstKmer < secondKmer)
            return true;
        if(secondKmer < firstKmer)
            return false;
        if(first.id < second.id)
            return true;
        if(second.id < first.id)
            return false;
        if(first.pos < second.pos)
            return true;
        if(second.pos < first.pos)
            return false;
        return false;
    }

    static bool compareRepSequenceAndIdAndDiag(const KmerPosition<T, IncludeAdjacentSeq> &first, const KmerPosition<T, IncludeAdjacentSeq> &second){
        if(first.kmer < second.kmer)
            return true;
        if(second.kmer < first.kmer)
            return false;
        if(first.id < second.id)
            return true;
        if(second.id < first.id)
            return false;
        if(first.pos < second.pos)
            return true;
        if(second.pos < first.pos)
            return false;
        return false;
    }
};



struct __attribute__((__packed__)) KmerEntry {
    unsigned int seqId;
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
    unsigned int seqId;
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
    size_t repSeq;
    unsigned int id;
    short pos;
    unsigned char score;
    unsigned int file;
    char reverse;
    FileKmerPosition(){}
    FileKmerPosition(size_t repSeq, unsigned int id,short pos, unsigned char score, unsigned int file):
            repSeq(repSeq), id(id), pos(pos), score(score), file(file), reverse(0) {}
    FileKmerPosition(size_t repSeq, unsigned int id,short pos, unsigned char score, char reverse, unsigned int file):
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


template  <int TYPE, typename T>
size_t assignGroup(KmerPosition<T, false> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr,
                   SequenceWeights * sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);
template  <int TYPE, typename T>
size_t assignGroup(KmerPosition<T, true> *kmers, size_t splitKmerCount, bool includeOnlyExtendable, int covMode, float covThr,
                   SequenceWeights * sequenceWeights, float weightThr, BaseMatrix *subMat, float &hashSeqBuffer, std::string tmpFile, int &numDiskBuffer);

template <int TYPE, typename T>
void mergeKmerFilesAndOutput(DBWriter & dbw, std::vector<std::string> tmpFiles, std::vector<char> &repSequence);

typedef std::priority_queue<FileKmerPosition, std::vector<FileKmerPosition>, CompareResultBySeqId> KmerPositionQueue;

template <int TYPE, typename T>
size_t queueNextEntry(KmerPositionQueue &queue, int file, size_t offsetPos, T *entries, size_t entrySize);

void setKmerLengthAndAlphabet(Parameters &parameters, size_t aaDbSize, int seqType);

template <int TYPE, typename T, typename seqLenType, bool IncludeAdjacentSeq>
void writeKmersToDisk(std::string tmpFile, KmerPosition<seqLenType, IncludeAdjacentSeq> *kmers, size_t totalKmers);

template <int TYPE, typename T, bool IncludeAdjacentSeq>
void writeKmerMatcherResult(DBWriter & dbw, KmerPosition<T, IncludeAdjacentSeq> *hashSeqPair, size_t totalKmers,
                            std::vector<char> &repSequence, size_t threads);

template <typename T, bool IncludeAdjacentSeq>
void resizeBuffer(size_t totalKmers, size_t hashStartRange, size_t hashEndRange, DBReader<unsigned int> & seqDbr, 
                   Parameters & par, BaseMatrix * subMat);
template <typename T, bool IncludeAdjacentSeq>
KmerPosition<T, IncludeAdjacentSeq> * doComputation(size_t &totalKmers, size_t split, size_t splits, std::string splitFile,
                                                    DBReader<unsigned int> & seqDbr, Parameters & par, BaseMatrix * subMat,
                                                    size_t KMER_SIZE, size_t chooseTopKmer, float chooseTopKmerScale = 0.0);
template <typename T, bool IncludeAdjacentSeq = false>
KmerPosition<T, IncludeAdjacentSeq> *initKmerPositionMemory(size_t size);

template <int TYPE, typename T, bool IncludeAdjacentSeq>
std::pair<size_t, size_t>  fillKmerPositionArray(KmerPosition<T, IncludeAdjacentSeq> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                 Parameters & par, BaseMatrix * subMat, bool hashWholeSequence,
                                                 size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);


void maskSequence(int maskMode, int maskLowerCase,
                  Sequence &seq, int maskLetter, ProbabilityMatrix * probMatrix);

template <typename T, bool IncludeAdjacentSeq = false>
size_t computeMemoryNeededLinearfilter(size_t totalKmer);

template <typename T, bool IncludeAdjacentSeq = false>
std::vector<std::pair<size_t, size_t>> setupKmerSplits(Parameters &par, BaseMatrix * subMat, DBReader<unsigned int> &seqDbr, size_t totalKmers, size_t splits);

size_t computeKmerCount(DBReader<unsigned int> &reader, size_t KMER_SIZE, size_t chooseTopKmer,
                        float chooseTopKmerScale = 0.0);

void setLinearFilterDefault(Parameters *p);

size_t computeMemoryNeededLinearfilter(size_t totalKmer);


#undef SIZE_T_MAX


#endif
