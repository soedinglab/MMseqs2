#ifndef MMSEQS_KMERMATCHER_H
#define MMSEQS_KMERMATCHER_H
#include <queue>
#include "DBWriter.h"
#include "Util.h"
#include "DBReader.h"
#include "Parameters.h"
#include "BaseMatrix.h"

struct KmerPosition {
    size_t kmer;
    unsigned int id;
    unsigned short seqLen;
    short pos;
    KmerPosition(){}
    KmerPosition(size_t kmer, unsigned int id, unsigned short seqLen, short pos):
            kmer(kmer), id(id), seqLen(seqLen), pos(pos) {}
    static bool compareRepSequenceAndIdAndPos(const KmerPosition &first, const KmerPosition &second){
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

    static bool compareRepSequenceAndIdAndPosReverse(const KmerPosition &first, const KmerPosition &second){
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

    static bool compareRepSequenceAndIdAndDiagReverse(const KmerPosition &first, const KmerPosition &second){
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


    static bool compareRepSequenceAndIdAndDiag(const KmerPosition &first, const KmerPosition &second){
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
    unsigned int file;
    char reverse;
    FileKmerPosition(){}
    FileKmerPosition(size_t repSeq, unsigned int id,short pos, unsigned int file):
            repSeq(repSeq), id(id), pos(pos), file(file), reverse(0) {}
    FileKmerPosition(size_t repSeq, unsigned int id,short pos, char reverse, unsigned int file):
            repSeq(repSeq), id(id), pos(pos), file(file), reverse(reverse) {}
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


template  <int TYPE>
size_t assignGroup(KmerPosition *kmers, size_t splitKmerCount, bool includeOnlyExtendable);

template <int TYPE, typename T>
void mergeKmerFilesAndOutput(DBReader<unsigned int> & seqDbr, DBWriter & dbw,
                             std::vector<std::string> tmpFiles, std::vector<char> &repSequence,
                             int covMode, float covThr) ;

typedef std::priority_queue<FileKmerPosition, std::vector<FileKmerPosition>, CompareResultBySeqId> KmerPositionQueue;

template <int TYPE, typename T>
size_t queueNextEntry(KmerPositionQueue &queue, int file, size_t offsetPos, T *entries, size_t entrySize);

void setKmerLengthAndAlphabet(Parameters &parameters, size_t aaDbSize, int seqType);

template <int TYPE, typename T>
void writeKmersToDisk(std::string tmpFile, KmerPosition *kmers, size_t totalKmers);

template <int TYPE>
void writeKmerMatcherResult(DBReader<unsigned int> & seqDbr, DBWriter & dbw,
                            KmerPosition *hashSeqPair, size_t totalKmers,
                            std::vector<char> &repSequence, int covMode, float covThr,
                            size_t threads);

KmerPosition * doComputation(size_t totalKmers, size_t split, size_t splits, std::string splitFile,
                             DBReader<unsigned int> & seqDbr, Parameters & par, BaseMatrix  * subMat,
                             size_t KMER_SIZE, size_t chooseTopKmer);
template  <int TYPE>
size_t fillKmerPositionArray(KmerPosition * hashSeqPair, DBReader<unsigned int> &seqDbr,
                             Parameters & par, BaseMatrix * subMat,
                             size_t KMER_SIZE, size_t chooseTopKmer,
                             size_t splits, size_t split);

size_t computeMemoryNeededLinearfilter(size_t totalKmer);

size_t computeKmerCount(DBReader<unsigned int> &reader, size_t KMER_SIZE, size_t chooseTopKmer);

size_t computeMemoryNeededLinearfilter(size_t totalKmer);

unsigned circ_hash(const int * x, unsigned length, const unsigned rol);

unsigned circ_hash_next(const int * x, unsigned length, int x_first, short unsigned h, const unsigned rol);



#undef SIZE_T_MAX


#endif
