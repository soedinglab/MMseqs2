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
struct AdjacentData {
    unsigned char data[6];
    void set(int i, unsigned char v) { data[i] = v; }
    unsigned char get(int i) { return data[i]; }
};

template <>
struct __attribute__((__packed__)) AdjacentData<false> {
    void set(int, unsigned char) { }
    unsigned char get(int) { return '\0'; }
};
template <typename T, bool IncludeSeqLen>
struct SeqLenData {};

template <typename T>
struct SeqLenData<T, true> {
    T seqlen;

    T getSeqLen(unsigned int) const { return seqlen; }
    void setSeqLen(T len) { seqlen = len; }
};

template <typename T>
struct SeqLenData<T, false> {
    static T* seqkey_to_len;

    T getSeqLen(unsigned int id) const { return seqkey_to_len[id]; }
    void setSeqLen(T) {}
};

template <typename T>
T* SeqLenData<T, false>::seqkey_to_len = NULL;
template <typename T, bool IncludeAdjacentSeq = false, bool IncludeSeqLen = false>
struct __attribute__((__packed__)) KmerPosition {
    size_t kmer;
    unsigned int id;
    T pos;
    SeqLenData<T, IncludeSeqLen> sl;
    AdjacentData<IncludeAdjacentSeq> adj;

    T getSeqLen() const {
        return sl.getSeqLen(id);
    }

    void setAdjacentSeq(int index, unsigned char val) {
        adj.set(index, val);
    }
    unsigned char getAdjacentSeq(int index) {
        return adj.get(index);
    }

    static bool compareRepSequenceAndIdAndPos(
            const KmerPosition &first, const KmerPosition &second) {
        if(first.kmer < second.kmer) return true;
        if(second.kmer < first.kmer) return false;
        if(first.id != SIZE_T_MAX && second.id != SIZE_T_MAX) {
            T len1 = first.getSeqLen();
            T len2 = second.getSeqLen();
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
        if(first.id != SIZE_T_MAX && second.id != SIZE_T_MAX) {
            T len1 = first.getSeqLen();
            T len2 = second.getSeqLen(); 
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


template <int TYPE, typename T, bool IncludeAdjacentSeq = false, bool IncludeSeqLen = false> 
size_t assignGroup(KmerPosition<T, IncludeAdjacentSeq, IncludeSeqLen> *kmers, KmerPosition<T, IncludeAdjacentSeq, IncludeSeqLen> *writeSeqPair, bool includeOnlyExtendable, int covMode, float covThr);

template <int TYPE, typename T, bool IncludeAdjacentSeq = false>
void mergeKmerFilesAndOutput(DBWriter & dbw, std::vector<std::string> tmpFiles, std::vector<char> &repSequence, int numThreads = 1, int maxIter = 1);

typedef std::priority_queue<FileKmerPosition, std::vector<FileKmerPosition>, CompareResultBySeqId> KmerPositionQueue;

template <int TYPE, typename T>
size_t queueNextEntry(KmerPositionQueue &queue, int file, size_t offsetPos, T *entries, size_t entrySize);

void setKmerLengthAndAlphabet(Parameters &parameters, size_t aaDbSize, int seqType);

template <int TYPE, typename T, typename seqLenType, bool IncludeAdjacentSeq = false, bool IncludeSeqLen = false>
void writeKmersToDisk(std::string tmpFile, KmerPosition<seqLenType, IncludeAdjacentSeq, IncludeSeqLen> *kmers, size_t totalKmers, int numThreads = 1, std::vector<size_t> *threadQueryOffsets = NULL, int iteration = 0);

template <int TYPE, typename T, bool IncludeAdjacentSeq = false, bool IncludeSeqLen = false>
void writeKmerMatcherResult(DBWriter & dbw, KmerPosition<T, IncludeAdjacentSeq, IncludeSeqLen> *hashSeqPair, size_t totalKmers,
                            std::vector<char> &repSequence, size_t threads);


template <typename T, bool IncludeAdjacentSeq = false, bool IncludeSeqLen = false>
KmerPosition<T, IncludeAdjacentSeq, IncludeSeqLen> * doComputation(size_t totalKmers, size_t split, size_t splits, std::string splitFile,
                                DBReader<unsigned int> & seqDbr, Parameters & par, BaseMatrix  * subMat,
                                size_t KMER_SIZE, size_t chooseTopKmer, float chooseTopKmerScale = 0.0);

template <typename T, bool IncludeAdjacentSeq = false, bool IncludeSeqLen = false>
KmerPosition<T, IncludeAdjacentSeq, IncludeSeqLen> *initKmerPositionMemory(size_t size);

template <int TYPE, typename T, bool IncludeAdjacentSeq = false, bool IncludeSeqLen = false>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPosition<T, IncludeAdjacentSeq, IncludeSeqLen> * kmerArray, size_t kmerArraySize, DBReader<unsigned int> &seqDbr,
                                                 Parameters & par, BaseMatrix * subMat, bool hashWholeSequence,
                                                 size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution);


void maskSequence(int maskMode, int maskLowerCase,
                  Sequence &seq, int maskLetter, ProbabilityMatrix * probMatrix);

template <typename T, bool IncludeAdjacentSeq = false, bool IncludeSeqLen = false>
size_t computeMemoryNeededLinearfilter(size_t totalKmer);

template <typename T, bool IncludeAdjacentSeq = false, bool IncludeSeqLen = false>
std::vector<std::pair<size_t, size_t>> setupKmerSplits(Parameters &par, BaseMatrix * subMat, DBReader<unsigned int> &seqDbr, size_t totalKmers, size_t splits);
size_t computeKmerCount(DBReader<unsigned int> &reader, size_t KMER_SIZE, size_t chooseTopKmer,
                        float chooseTopKmerScale = 0.0);

void setLinearFilterDefault(Parameters *p);

size_t computeMemoryNeededLinearfilter(size_t totalKmer);


#undef SIZE_T_MAX


#endif
