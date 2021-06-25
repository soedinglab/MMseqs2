#ifndef MMSEQS_KMERINDEX_H
#define MMSEQS_KMERINDEX_H

// Written by Martin Steinegger martin.steinegger@snu.ac.kr
// storage for k-mers
#include "MathUtil.h"
#include <string>
#include <algorithm>
#include <fcntl.h>
#include <sys/mman.h>

class KmerIndex{

private:
    // index data
    const size_t indexGridResolution = 32768;
    size_t indexGridSize;
    size_t * entryOffsets;
    size_t prevKmerStartRange;
    long long iteratorPos;
    size_t entryOffsetPos;
    size_t writingPosition;
    // total entries count
    size_t entryCount;
    size_t maxWriteEntries;

    bool isMmaped;
    int alphabetSize;
    int kmerSize;

    struct __attribute__((__packed__)) KmerEntryRelative {
        unsigned int id;
        unsigned short kmerOffset;
        unsigned short pos;
        unsigned short seqLen;
    };
    KmerEntryRelative * entries;

public:
    struct KmerEntry {
        size_t kmer;
        unsigned int id;
        unsigned short pos;
        unsigned short seqLen;
        KmerEntry(size_t kmer, unsigned int id, unsigned short pos, unsigned short seqLen) :
                kmer(kmer), id(id), pos(pos), seqLen(seqLen) {}

        KmerEntry() {}

    };


    KmerIndex(){}

    void init(const size_t alphabetSize,
              const size_t kmerSize,
              const size_t entryCount){
        this->isMmaped = false;
        this->entries = new KmerEntryRelative[entryCount];
        this->maxWriteEntries = entryCount;
        this->entryCount = 0;
        this->indexGridSize = MathUtil::ceilIntDivision( MathUtil::ipow<size_t>(alphabetSize, kmerSize),indexGridResolution);
        this->entryOffsets = new size_t[indexGridSize+1];
        memset(entryOffsets, 0, sizeof(size_t)*(indexGridSize + 1));
        this->prevKmerStartRange = 0;
        this->writingPosition = 0;
    }
    // use for writing the data to file
    KmerIndex(const size_t alphabetSize,
              const size_t kmerSize) : alphabetSize(alphabetSize), kmerSize(kmerSize)  {
        init(alphabetSize, kmerSize, indexGridResolution);
    }


    bool hasNextEntry(){
        return (iteratorPos + 1 < static_cast<long long>(entryCount));
    }

    template <int TYPE>
    KmerEntry getNextEntry(){
        iteratorPos++;
        while(iteratorPos >= static_cast<long long>(entryOffsets[entryOffsetPos+1])){
            entryOffsetPos++;
        }
        size_t kmer;

        if(TYPE==Parameters::DBTYPE_NUCLEOTIDES){
            bool isReverse = BIT_CHECK(entries[iteratorPos].kmerOffset, 15);
            kmer = BIT_CLEAR(entries[iteratorPos].kmerOffset, 15);
            kmer = entryOffsetPos*indexGridResolution + kmer;
            kmer = (isReverse) ? kmer :  BIT_SET(kmer, 63);
        }else{
//            std::cout << entryOffsetPos <<  "\t" << iteratorPos << "\t" << entryOffsets[entryOffsetPos] <<  "\t" << entries[iteratorPos].kmerOffset << std::endl;
            kmer= entryOffsetPos*indexGridResolution + static_cast<size_t >(entries[iteratorPos].kmerOffset);
        }
        return KmerEntry(kmer, entries[iteratorPos].id, entries[iteratorPos].pos, entries[iteratorPos].seqLen);
    }


    size_t getGridPosition(size_t kmer){
        return kmer / indexGridResolution;
    }

    size_t getKmerStartRange(size_t kmer) {
        return (kmer / indexGridResolution ) * indexGridResolution;
    }

    size_t getGridResolution(){
        return indexGridResolution;
    }

    void setupOffsetTable(){
        size_t prevElementLength = entryOffsets[0];
        entryOffsets[0] = 0;
        for(size_t i = 0; i < indexGridSize; i++) {
            const size_t currElementLength = entryOffsets[i + 1];
            entryOffsets[i + 1] = entryOffsets[i] + prevElementLength;
//            if(currElementLength  > 0){
//                std::cout  << i << "\t" << entryOffsets[i + 1] << std::endl;
//            }
            prevElementLength = currElementLength;
        }
    }

    bool needsFlush(const size_t kmer){
        bool needsFlash = false;
        //const size_t gridPosition = getGridPosition(kmer);
        size_t kmerStartRange = getKmerStartRange(kmer);
        if(kmerStartRange != prevKmerStartRange){
            needsFlash = true;
        }
        return needsFlash;
    }

    void flush(DBWriter & writer){
        writer.writeAdd((char*)entries, writingPosition* sizeof(KmerEntryRelative), 0);
        writingPosition = 0;
    }

    void addElementSorted(const size_t kmer, unsigned int id, unsigned short pos, unsigned short seqLen, bool isReverse){
        const size_t gridPosition = getGridPosition(kmer);
        size_t kmerStartRange = getKmerStartRange(kmer);
        prevKmerStartRange = kmerStartRange;
        entryOffsets[gridPosition]++;
        if(writingPosition >= maxWriteEntries){
            Debug(Debug::ERROR) << "addElement overflows. Current write position is " << writingPosition << "\n";
            EXIT(EXIT_FAILURE);
        }
        entries[writingPosition].id = id;
        entries[writingPosition].kmerOffset = kmer - kmerStartRange;
        entries[writingPosition].kmerOffset = (isReverse) ? BIT_SET(entries[writingPosition].kmerOffset, 15) :  entries[writingPosition].kmerOffset;
        entries[writingPosition].pos = pos;
        entries[writingPosition].seqLen = seqLen;
        writingPosition++;
        entryCount++;
    }
//
//    std::pair<size_t, size_t> getEntryRange(size_t kmer) {
//        size_t gridPosition = getGridPosition(kmer);
//        size_t abundanceEntrySize = entryOffsets[gridPosition + 1] - entryOffsets[gridPosition];
//        return std::make_pair(entryOffsets[gridPosition], abundanceEntrySize);
//    }

    const size_t  *getOffsets() {
        return entryOffsets;
    }

    size_t getOffsetsSize() {
        return indexGridSize;
    }

    uint64_t getTableEntriesNum() {
        return entryCount;
    }

    KmerIndex(int alphabetSize, int kmerSize, char *entriesData, char *entriesOffetData,
              size_t entryCount, size_t gridResolution) {
        this->alphabetSize = alphabetSize;
        this->kmerSize = kmerSize;
        this->isMmaped = true;
        this->entries =  (KmerEntryRelative * )entriesData;
        this->entryCount = entryCount;
        this->indexGridSize = MathUtil::ceilIntDivision( MathUtil::ipow<size_t>(alphabetSize, kmerSize), gridResolution );
        this->entryOffsets = (size_t *) entriesOffetData;
#if HAVE_POSIX_MADVISE
        if (entryCount > 0 && posix_madvise (entriesData, entryCount* sizeof(KmerEntryRelative), POSIX_MADV_SEQUENTIAL) != 0){
            Debug(Debug::ERROR) << "KmerIndex posix_madvise returned an error\n";
        }
#endif

        this->prevKmerStartRange = 0;
        this->iteratorPos = -1;
        this->entryOffsetPos = 0;
    }

    ~KmerIndex(){
        if(isMmaped==false){
            if(entries){
                delete[] this->entries;
            }
            if(entryOffsets){
                delete[] this->entryOffsets;
            }
        }
    }

    template <int TYPE>
    void printIndex(BaseMatrix * mat){
        Indexer indexer(alphabetSize, kmerSize);
        reset();
        Debug(Debug::INFO)  << "Entry Count: " << entryCount << "\n";
        size_t id = 0;
        while(hasNextEntry()){
            KmerEntry kmer = getNextEntry<TYPE>();
            Debug(Debug::INFO) << id++ << "\t";
            size_t kmerIdx = kmer.kmer;
            if(TYPE==Parameters::DBTYPE_NUCLEOTIDES){
                kmerIdx = BIT_CLEAR(kmerIdx, 15);
                Indexer::printKmer(kmerIdx, kmerSize);
//                indexer.printKmer(kmer.kmer, kmerSize, mat->num2aa);
            }else{
                indexer.printKmer(kmerIdx, kmerSize, mat->num2aa);
            }
            Debug(Debug::INFO) << "\t";
            Debug(Debug::INFO) << kmerIdx << "\t";
            Debug(Debug::INFO) << kmer.id << "\t";
            Debug(Debug::INFO) << kmer.pos << "\t";
            Debug(Debug::INFO) << kmer.seqLen << "\t";
            Debug(Debug::INFO) << "\n";
        }
    }

    void reset() {
        this->iteratorPos = -1;
        this->entryOffsetPos = 0;
    }
};

#endif //MMSEQS_KMERINDEX_H
