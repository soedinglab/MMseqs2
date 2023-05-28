#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

//
// Written by Martin Steinegger martin.steinegger@snu.ac.kr and Maria Hauser mhauser@genzentrum.lmu.de
//
// Abstract: Index table stores the list of DB sequences containing a certain k-mer, for each k-mer.
//

#include "DBReader.h"
#include "Sequence.h"
#include "Indexer.h"
#include "Debug.h"
#include "Util.h"
#include "SequenceLookup.h"
#include "MathUtil.h"
#include "KmerGenerator.h"
#include "Parameters.h"
#include "FastSort.h"
#include <stdlib.h>
#include <algorithm>

// IndexEntryLocal is an entry with position and seqId for a kmer
// structure needs to be packed or it will need 8 bytes instead of 6
struct __attribute__((__packed__)) IndexEntryLocal {
    unsigned int seqId;
    unsigned short position_j;
    static bool comapreByIdAndPos(IndexEntryLocal first, IndexEntryLocal second){
        if(first.seqId < second.seqId )
            return true;
        if(second.seqId < first.seqId )
            return false;
        if(first.position_j < second.position_j )
            return true;
        if(second.position_j < first.position_j )
            return false;
        return false;
    }
};

struct __attribute__((__packed__)) IndexEntryLocalTmp {
    unsigned int kmer;
    unsigned int seqId;
    unsigned short position_j;

    IndexEntryLocalTmp(unsigned int kmer, unsigned int seqId, unsigned short position_j)
            :kmer(kmer),seqId(seqId), position_j(position_j)
    {}

    IndexEntryLocalTmp() {}

    static bool comapreByIdAndPos(IndexEntryLocalTmp first, IndexEntryLocalTmp second){
        if(first.kmer < second.kmer )
            return true;
        if(second.kmer < first.kmer )
            return false;
        if(first.position_j < second.position_j )
            return true;
        if(second.position_j < first.position_j )
            return false;
        return false;
    }
};

class IndexTable {
public:
    IndexTable(int alphabetSize, int kmerSize, bool externalData)
            : tableSize(MathUtil::ipow<size_t>(alphabetSize, kmerSize)), alphabetSize(alphabetSize),
              kmerSize(kmerSize), externalData(externalData), tableEntriesNum(0), size(0),
              indexer(new Indexer(alphabetSize, kmerSize)), entries(NULL), offsets(NULL) {
        if (externalData == false) {
            offsets = new(std::nothrow) size_t[tableSize + 1];
            Util::checkAllocation(offsets, "Can not allocate entries memory in IndexTable");
            memset(offsets, 0, (tableSize + 1) * sizeof(size_t));
        }
    }

    virtual ~IndexTable() {
        deleteEntries();
        delete indexer;
    }

    void deleteEntries() {
        if (externalData == false) {
            if (entries != NULL) {
                delete[] entries;
                entries = NULL;
            }
            if (offsets != NULL) {
                delete[] offsets;
                offsets = NULL;
            }
        }
    }

    // count k-mers in the sequence, so enough memory for the sequence lists can be allocated in the end
    size_t addSimilarKmerCount(Sequence* s, KmerGenerator* kmerGenerator){
        s->resetCurrPos();
        std::vector<unsigned int> seqKmerPosBuffer;

        //idxer->reset();
        while(s->hasNextKmer()){
            const unsigned char * kmer = s->nextKmer();
            if(s->kmerContainsX()){
                continue;
            }
            const std::pair<size_t *, size_t> kmerList = kmerGenerator->generateKmerList(kmer);

            //unsigned int kmerIdx = idxer->int2index(kmer, 0, kmerSize);
            for(size_t i = 0; i < kmerList.second; i++){
                seqKmerPosBuffer.push_back(kmerList.first[i]);
            }
        }
        if(seqKmerPosBuffer.size() > 1){
            SORT_SERIAL(seqKmerPosBuffer.begin(), seqKmerPosBuffer.end());
        }
        size_t countUniqKmer = 0;
        unsigned int prevKmerIdx = UINT_MAX;
        for(size_t i = 0; i < seqKmerPosBuffer.size(); i++){
            unsigned int kmerIdx = seqKmerPosBuffer[i];
            if(prevKmerIdx != kmerIdx){
                //table[kmerIdx] += 1;
                // size increases by one
                __sync_fetch_and_add(&(offsets[kmerIdx]), 1);
                countUniqKmer++;
            }
            prevKmerIdx = kmerIdx;
        }
        return countUniqKmer;
    }

    // count k-mers in the sequence, so enough memory for the sequence lists can be allocated in the end
    size_t addKmerCount(Sequence *s, Indexer *idxer, unsigned int *seqKmerPosBuffer,
                        int threshold, char *diagonalScore) {
        s->resetCurrPos();
        size_t countKmer = 0;
        bool removeX = (Parameters::isEqualDbtype(s->getSequenceType(), Parameters::DBTYPE_NUCLEOTIDES) ||
                        Parameters::isEqualDbtype(s->getSequenceType(), Parameters::DBTYPE_AMINO_ACIDS));
        while(s->hasNextKmer()){
            const unsigned char * kmer = s->nextKmer();
            if(removeX && s->kmerContainsX()){
                continue;
            }
            if(threshold > 0){
                int score = 0;
                for(int pos = 0; pos < kmerSize; pos++){
                    score += diagonalScore[kmer[pos]];
                }
                if(score < threshold){
                    continue;
                }
            }
            unsigned int kmerIdx = idxer->int2index(kmer, 0, kmerSize);
            seqKmerPosBuffer[countKmer] = kmerIdx;
            countKmer++;
        }
        if(countKmer > 1){
            SORT_SERIAL(seqKmerPosBuffer, seqKmerPosBuffer + countKmer);
        }
        size_t countUniqKmer = 0;
        unsigned int prevKmerIdx = UINT_MAX;
        for(size_t i = 0; i < countKmer; i++){
            unsigned int kmerIdx = seqKmerPosBuffer[i];
            if(prevKmerIdx != kmerIdx){
                //table[kmerIdx] += 1;
                // size increases by one
                __sync_fetch_and_add(&(offsets[kmerIdx]), 1);
                countUniqKmer++;
            }
            prevKmerIdx = kmerIdx;
        }
        return countUniqKmer;
    }

    // get list of DB sequences containing this k-mer
    inline IndexEntryLocal *getDBSeqList(size_t kmer, size_t *matchedListSize) {
        const ptrdiff_t diff = offsets[kmer + 1] - offsets[kmer];
        *matchedListSize = static_cast<size_t>(diff);
        return (entries + offsets[kmer]);
    }

    void sortDBSeqLists() {
        #pragma omp parallel for
        for (size_t i = 0; i < tableSize; i++) {
            size_t entrySize;
            IndexEntryLocal *entries = getDBSeqList(i, &entrySize);
            SORT_SERIAL(entries, entries + entrySize, IndexEntryLocal::comapreByIdAndPos);
        }
    }

    // get pointer to entries array
    IndexEntryLocal *getEntries() {
        return entries;
    }

    inline size_t getOffset(size_t kmer) {
        return offsets[kmer];
    }

    size_t *getOffsets() {
        return offsets;
    }

    // init the arrays for the sequence lists
    void initMemory(size_t dbSize) {
        size_t tableEntriesNum = 0;
        for (size_t i = 0; i < getTableSize(); i++) {
            tableEntriesNum += getOffset(i);
        }

        this->tableEntriesNum = tableEntriesNum;
        this->size = dbSize; // amount of sequences added

        // allocate memory for the sequence id lists
        entries = new(std::nothrow) IndexEntryLocal[tableEntriesNum];
        Util::checkAllocation(entries, "Can not allocate entries memory in IndexTable::initMemory");
    }

    // allocates memory for index tables
    void init() {
        // set the pointers in the index table to the start of the list for a certain k-mer
        size_t offset = 0;
        for (size_t i = 0; i < tableSize; i++) {
            const size_t currentOffset = offsets[i];
            offsets[i] = offset;
            offset += currentOffset;
        }
        offsets[tableSize] = offset;
    }

    // init index table with external data (needed for index readin)
    void initTableByExternalData(size_t sequenceCount, size_t tableEntriesNum, IndexEntryLocal *entries, size_t *entryOffsets) {
        this->tableEntriesNum = tableEntriesNum;
        this->size = sequenceCount;

        this->entries = entries;
        this->offsets = entryOffsets;
    }

    void initTableByExternalDataCopy(size_t sequenceCount, size_t tableEntriesNum, IndexEntryLocal *entries, size_t *entryOffsets) {
        this->tableEntriesNum = tableEntriesNum;
        this->size = sequenceCount;

        this->entries = new(std::nothrow) IndexEntryLocal[tableEntriesNum];
        Util::checkAllocation(entries, "Can not allocate " + SSTR(tableEntriesNum * sizeof(IndexEntryLocal)) + " bytes for entries in IndexTable::initMemory");
        memcpy(this->entries, entries, tableEntriesNum * sizeof(IndexEntryLocal));

        memcpy(this->offsets, entryOffsets, (tableSize + 1) * sizeof(size_t));
    }

    void revertPointer() {
        for (size_t i = tableSize; i > 0; i--) {
            offsets[i] = offsets[i - 1];
        }
        offsets[0] = 0;
    }

    void printStatistics(char *num2aa) {
        const size_t top_N = 10;
        std::pair<size_t, size_t> topElements[top_N];
        for (size_t j = 0; j < top_N; j++) {
            topElements[j].first = 0;
        }

        size_t entrySize = 0;
        size_t minKmer = 0;
        size_t emptyKmer = 0;
        for (size_t i = 0; i < tableSize; i++) {
            const ptrdiff_t size = offsets[i + 1] - offsets[i];
            minKmer = std::min(minKmer, (size_t) size);
            entrySize += size;
            if (size == 0) {
                emptyKmer++;
            }
            if (((size_t) size) < topElements[top_N - 1].first)
                continue;
            for (size_t j = 0; j < top_N; j++) {
                if (topElements[j].first < ((size_t) size)) {
                    topElements[j].first = static_cast<unsigned long>(size);
                    topElements[j].second = i;
                    break;
                }
            }
        }

        double avgKmer = ((double) entrySize) / ((double) tableSize);
        Debug(Debug::INFO) << "Index statistics\n";
        Debug(Debug::INFO) << "Entries:          " << entrySize << "\n";
        Debug(Debug::INFO) << "DB size:          " << (entrySize * sizeof(IndexEntryLocal) + tableSize * sizeof(size_t))/1024/1024 << " MB\n";
        Debug(Debug::INFO) << "Avg k-mer size:   " << avgKmer << "\n";
        Debug(Debug::INFO) << "Top " << top_N << " k-mers\n";
        for (size_t j = 0; j < top_N; j++) {
            Debug(Debug::INFO) << "    ";
            indexer->printKmer(topElements[j].second, kmerSize, num2aa);
            Debug(Debug::INFO) << "\t" << topElements[j].first << "\n";
        }
    }

    // FUNCTIONS TO OVERWRITE
    // add k-mers of the sequence to the index table
    void addSimilarSequence(Sequence* s, KmerGenerator* kmerGenerator, IndexEntryLocalTmp ** buffer, size_t &bufferSize, Indexer * idxer) {
        // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
        s->resetCurrPos();
        idxer->reset();
        size_t kmerPos = 0;
        while(s->hasNextKmer()){
            const unsigned char * kmer = s->nextKmer();
            if(s->kmerContainsX()){
                continue;
            }
            std::pair<size_t *, size_t> scoreMatrix = kmerGenerator->generateKmerList(kmer);
            if(kmerPos+scoreMatrix.second >= bufferSize){
                *buffer = static_cast<IndexEntryLocalTmp*>(realloc(*buffer, sizeof(IndexEntryLocalTmp) * bufferSize*2));
                bufferSize = bufferSize*2;
            }
            for(size_t i = 0; i < scoreMatrix.second; i++) {
                unsigned int kmerIdx = scoreMatrix.first[i];

                // if region got masked do not add kmer
                if (offsets[kmerIdx + 1] - offsets[kmerIdx] == 0)
                    continue;
                (*buffer)[kmerPos].kmer = kmerIdx;
                (*buffer)[kmerPos].seqId = s->getId();
                (*buffer)[kmerPos].position_j = s->getCurrentPosition();
                kmerPos++;
            }

        }

        if(kmerPos>1){
            SORT_SERIAL(*buffer, *buffer+kmerPos, IndexEntryLocalTmp::comapreByIdAndPos);
        }
        unsigned int prevKmer = UINT_MAX;
        for(size_t pos = 0; pos < kmerPos; pos++){
            unsigned int kmerIdx = (*buffer)[pos].kmer;
            if(kmerIdx != prevKmer){
                size_t offset = __sync_fetch_and_add(&(offsets[kmerIdx]), 1);
                IndexEntryLocal *entry = &entries[offset];
                entry->seqId      = (*buffer)[pos].seqId;
                entry->position_j = (*buffer)[pos].position_j;
            }
            prevKmer = kmerIdx;
        }
    }

    // add k-mers of the sequence to the index table
    void addSequence (Sequence* s, Indexer * idxer,
                      IndexEntryLocalTmp ** buffer, size_t bufferSize,
                      int threshold, char * diagonalScore){
        // iterate over all k-mers of the sequence and add the id of s to the sequence list of the k-mer (tableDummy)
        s->resetCurrPos();
        idxer->reset();
        size_t kmerPos = 0;
        bool removeX = (Parameters::isEqualDbtype(s->getSequenceType(), Parameters::DBTYPE_NUCLEOTIDES) ||
                        Parameters::isEqualDbtype(s->getSequenceType(), Parameters::DBTYPE_AMINO_ACIDS));
        while (s->hasNextKmer()){
            const unsigned char * kmer = s->nextKmer();
            if(removeX && s->kmerContainsX()){
                continue;
            }
            if(threshold > 0) {
                int score = 0;
                for (int pos = 0; pos < kmerSize; pos++) {
                    score += diagonalScore[kmer[pos]];
                }
                if (score < threshold) {
                    continue;
                }
            }
            unsigned int kmerIdx = idxer->int2index(kmer, 0, kmerSize);
            // if region got masked do not add kmer
            if (offsets[kmerIdx + 1] - offsets[kmerIdx] == 0)
                continue;

            (*buffer)[kmerPos].kmer = kmerIdx;
            (*buffer)[kmerPos].seqId      = s->getId();
            (*buffer)[kmerPos].position_j = s->getCurrentPosition();
            kmerPos++;
            if(kmerPos >= bufferSize){
                *buffer = static_cast<IndexEntryLocalTmp*>(realloc(*buffer, sizeof(IndexEntryLocalTmp) * bufferSize*2));
                bufferSize = bufferSize*2;
            }
        }

        if(kmerPos>1){
            SORT_SERIAL(*buffer, *buffer+kmerPos, IndexEntryLocalTmp::comapreByIdAndPos);
        }

        unsigned int prevKmer = UINT_MAX;
        for(size_t pos = 0; pos < kmerPos; pos++){
            unsigned int kmerIdx = (*buffer)[pos].kmer;
            if(kmerIdx != prevKmer){
                size_t offset = __sync_fetch_and_add(&(offsets[kmerIdx]), 1);
                IndexEntryLocal *entry = &entries[offset];
                entry->seqId      = (*buffer)[pos].seqId;
                entry->position_j = (*buffer)[pos].position_j;
            }
            prevKmer = kmerIdx;
        }
    }

    // prints the IndexTable
    void print(char *num2aa) {
        for (size_t i = 0; i < tableSize; i++) {
            ptrdiff_t entrySize = offsets[i + 1] - offsets[i];
            if (entrySize > 0) {
                indexer->printKmer(i, kmerSize, num2aa);

                Debug(Debug::INFO) << "\n";
                IndexEntryLocal *e = &entries[offsets[i]];
                for (ptrdiff_t j = 0; j < entrySize; j++) {
                    Debug(Debug::INFO) << "\t(" << e[j].seqId << ", " << e[j].position_j << ")\n";
                }
            }
        }
    };

    // get amount of sequences in Index
    size_t getSize() { return size; };

    // returns the size of  table entries
    uint64_t getTableEntriesNum() { return tableEntriesNum; };

    // returns table size
    size_t getTableSize() { return tableSize; };

    // returns the size of the entry (int for global) (IndexEntryLocal for local)
    size_t getSizeOfEntry() { return sizeof(IndexEntryLocal); }

    int getKmerSize() {
        return kmerSize;
    }

    int getAlphabetSize() {
        return alphabetSize;
    }

    static int computeKmerSize(size_t aaSize) {
        return aaSize < getUpperBoundAACountForKmerSize(6) ? 6 : 7;
    }


    static size_t getUpperBoundAACountForKmerSize(int kmerSize) {
        switch (kmerSize) {
            case 6:
                return 3350000000;
            case 7:
                return (SIZE_MAX - 1); // SIZE_MAX is often reserved as safe flag
            default:
                Debug(Debug::ERROR) << "Invalid kmer size of " << kmerSize << "!\n";
                EXIT(EXIT_FAILURE);
        }
    }

    static size_t getUpperBoundNucCountForKmerSize(int kmerSize) {
        switch (kmerSize) {
            case 14:
                return 3350000000;
            case 15:
                return (SIZE_MAX - 1); // SIZE_MAX is often reserved as safe flag
            default:
                Debug(Debug::ERROR) << "Invalid kmer size of " << kmerSize << "!\n";
                EXIT(EXIT_FAILURE);
        }
    }

protected:
    // alphabetSize**kmerSize
    const size_t tableSize;
    const int alphabetSize;
    const int kmerSize;

    // external data from mmap
    const bool externalData;

    // number of entries in all sequence lists - must be 64bit
    uint64_t tableEntriesNum;
    // number of sequences in Index
    size_t size;

    Indexer *indexer;

    // Index table entries: ids of sequences containing a certain k-mer, stored sequentially in the memory
    IndexEntryLocal *entries;
    size_t *offsets;

    // sequence lookup
    SequenceLookup *sequenceLookup;
};
#endif
