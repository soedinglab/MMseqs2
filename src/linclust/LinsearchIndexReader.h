#ifndef MMSEQS_INDEXKMERDB_H
#define MMSEQS_INDEXKMERDB_H

#include "kmermatcher.h"


struct FileKmer {
    size_t kmer;
    unsigned int id;
    short pos;
    unsigned short seqLen;
    unsigned int file;
    FileKmer(){}
    FileKmer(size_t kmer, unsigned int id, short pos, short seqLen, unsigned int file):
            kmer(kmer), id(id), pos(pos), seqLen(seqLen), file(file) {}

};

class CompareRepSequenceAndIdAndDiag {
public:
    bool operator() (FileKmer & first, FileKmer & second) const {
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

class CompareRepSequenceAndIdAndDiagReverse {
public:
    bool operator() (FileKmer & first, FileKmer & second) const {
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
};

class LinsearchIndexReader {
public:


    template<int TYPE>
    static size_t pickCenterKmer(KmerPosition *kmers, size_t splitKmerCount);

    template<int TYPE, class Comp>
    static void mergeAndWriteIndex(DBWriter &dbw, std::vector<std::string> tmpFiles, int alphSize, int kmerSize);

    template<int TYPE>
    static void writeIndex(DBWriter &dbw,
                           KmerPosition *hashSeqPair, size_t totalKmers,
                           int alphSize, int kmerSize);

    static std::string indexName(std::string baseName);

    static void writeKmerIndexToDisk(std::string fileName, KmerPosition *kmers, size_t kmerCnt);

    static bool checkIfIndexFile(DBReader<unsigned int> *pReader);

    static bool isIndexCompatible(DBReader<unsigned int> & index, Parameters &parameters, int dbtype);

    static std::string searchForIndex(std::string dbName);
};
#endif
