//
// Created by Martin Steinegger on 2/25/16.
//

#include <limits>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <Indexer.h>
#include "ReducedMatrix.h"
#include "DBWriter.h"
#include "SubstitutionMatrix.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "omptl/omptl_algorithm"
#include <sys/time.h>
#include <MathUtil.h>

#ifdef OPENMP
#include <omp.h>
#endif

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif
struct KmerPosition {
    size_t kmer;
    unsigned int id;
    unsigned short size;
    unsigned short pos;
    
    static bool compareKmerPositionByKmer(KmerPosition first, KmerPosition second){
        return (first.kmer < second.kmer) ? true : false;
    }
    
    static bool compareKmerPositionByKmerAndSize(KmerPosition first, KmerPosition second){
        if(first.size > second.size )
            return true;
        if(second.size > first.size )
            return false;
        if(first.kmer < second.kmer )
            return true;
        if(second.kmer < first.kmer )
            return false;
        return false;
    }
    static bool compareRepSequenceAndId(KmerPosition first, KmerPosition second){
        if(first.kmer < second.kmer )
            return true;
        if(second.kmer < first.kmer )
            return false;
        if(first.id < second.id )
            return true;
        if(second.id < first.id )
            return false;
        return false;
    }
};

struct SequencePosition{
    size_t kmer;
    unsigned int pos;
    float score;
    static bool compareByScore(SequencePosition first, SequencePosition second){
        return (first.score < second.score) ? true : false;
    }
};

void setLinearFilterDefault(Parameters *p) {
    p->spacedKmer = false;
    p->covThr = 0.8;
    p->kmerSize = 10;
    p->alphabetSize = 13;
    p->kmersPerSequence = 20;
}


size_t computeKmerCount(DBReader<unsigned int> &reader, size_t KMER_SIZE, size_t chooseTopKmer) {
    size_t totalKmers = 0;
    for(size_t id = 0; id < reader.getSize(); id++ ){
        int kmerAdjustedSeqLen = std::max(0, static_cast<int>(reader.getSeqLens(id)) - static_cast<int>(KMER_SIZE - 1)-2) ;
        totalKmers += std::min(kmerAdjustedSeqLen, static_cast<int>( chooseTopKmer ) );
    }
    return totalKmers;
}

size_t computeMemoryNeededLinearfilter(size_t totalKmer) {
    return sizeof(KmerPosition) * totalKmer;
}


#define RoL(val, numbits) (val << numbits) ^ (val >> (32 - numbits))
unsigned circ_hash(const int * x, unsigned length){
    short unsigned RAND[21] = {0x4567, 0x23c6, 0x9869, 0x4873, 0xdc51, 0x5cff, 0x944a, 0x58ec, 0x1f29, 0x7ccd, 0x58ba, 0xd7ab, 0x41f2, 0x1efb, 0xa9e3, 0xe146, 0x007c, 0x62c2, 0x0854, 0x27f8, 0x231b};
    short unsigned h = 0x0;
    h = h^ RAND[x[0]];                  // XOR h and ki      
    for (int i = 1; i < length; ++i){
        h = RoL(h, 5); 
        h ^= RAND[x[i]];                   // XOR h and ki      
    }
    return h;
}

// Rolling hash for CRC variant: compute hash value for next key x[0:length-1] from previous hash value hash( x[-1:length-2] ) and x_first = x[-1]
unsigned circ_hash_next(const int * x, unsigned length, int x_first, short unsigned h){
    short unsigned RAND[21] = {0x4567, 0x23c6, 0x9869, 0x4873, 0xdc51, 0x5cff, 0x944a, 0x58ec, 0x1f29, 0x7ccd, 0x58ba, 0xd7ab, 0x41f2, 0x1efb, 0xa9e3, 0xe146, 0x007c, 0x62c2, 0x0854, 0x27f8, 0x231b};
    h ^= RoL(RAND[x_first], (5*(length-1)) % 16); // undo INITIAL_VALUE and first letter x[0] of old key
    h =  RoL(h, 5); // circularly permute all letters x[1:length-1] to 5 positions to left
    h ^= RAND[x[length-1]]; // add new, last letter of new key x[1:length]
    return h;
}
#undef RoL

size_t fillKmerPositionArray(KmerPosition * hashSeqPair, DBReader<unsigned int> &seqDbr,
                             Parameters & par, BaseMatrix * subMat,
                             size_t KMER_SIZE, size_t chooseTopKmer){
    size_t kmerCounter = 0;
#pragma omp parallel reduction (+: kmerCounter)
    {
        Sequence seq(par.maxSeqLen, subMat->aa2int, subMat->int2aa, Sequence::AMINO_ACIDS, KMER_SIZE, false, false);
        Indexer idxer(subMat->alphabetSize, KMER_SIZE);

        size_t dbFrom, dbSize;
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        Util::decomposeDomainByAminoAcid(seqDbr.getAminoAcidDBSize(), seqDbr.getSeqLens(), seqDbr.getSize(),
                                         thread_idx, par.threads, &dbFrom, &dbSize);
        size_t kmerStartPos = 0;
        for(size_t id = 0; id < dbFrom; id++ ){
            int kmerAdjustedSeqLen = std::max(0, static_cast<int>(seqDbr.getSeqLens(id)) - static_cast<int>(KMER_SIZE)  - 2);
            kmerStartPos += std::min(kmerAdjustedSeqLen, static_cast<int>(chooseTopKmer));
        }
        KmerPosition * tmpHashSeqPair=hashSeqPair + kmerStartPos;
        SequencePosition * kmers = new SequencePosition[par.maxSeqLen];
        for(size_t id = dbFrom; id < (dbFrom + dbSize); id++){
            Debug::printProgress(id);
            seq.mapSequence(id, id, seqDbr.getData(id));
            Util::maskLowComplexity(subMat, &seq, seq.L, 12, 3,
                                    par.alphabetSize, seq.aa2int[(unsigned char) 'X'], true, false, false, true);
            int seqKmerCount = 0;
            unsigned int seqId = seq.getId();
            unsigned int prevHash=0;
            unsigned int prevFirstRes = 0;
            if(seq.hasNextKmer()){
                    const int *kmer = seq.nextKmer();
                    prevHash = circ_hash(kmer, KMER_SIZE);
                    prevFirstRes = kmer[0];
            }
            while(seq.hasNextKmer()) {
                const int *kmer = seq.nextKmer();
                //float kmerScore = 1.0;
                prevHash = circ_hash_next(kmer, KMER_SIZE, prevFirstRes, prevHash);
                prevFirstRes = kmer[0];
                float kmerScore = prevHash;
                size_t xCount = 0;
                for (size_t kpos = 0; kpos < KMER_SIZE; kpos++) {
                    xCount += (kmer[kpos] == subMat->aa2int[(int)'X']);
                }
                if(xCount > 0){
                    continue;
                }
                
                (kmers+seqKmerCount)->score = kmerScore;
                size_t kmerIdx = idxer.int2index(kmer, 0, KMER_SIZE);
                (kmers+seqKmerCount)->kmer = kmerIdx;
                (kmers+seqKmerCount)->pos = seq.getCurrentPosition();
                seqKmerCount++;
            }
            if(seqKmerCount > 1){
                std::sort(kmers, kmers+seqKmerCount, SequencePosition::compareByScore);
            }
            size_t kmerConsidered = std::min(static_cast<int>(chooseTopKmer), seqKmerCount);
            
            for(size_t topKmer = 0; topKmer < kmerConsidered; topKmer++){
                tmpHashSeqPair->kmer  = (kmers + topKmer)->kmer;
                tmpHashSeqPair->id    = seqId;
                tmpHashSeqPair->pos   = (kmers + topKmer)->pos;
                //idxer.printKmer(tmpHashSeqPair->kmer, KMER_SIZE, subMat->int2aa);
                //std::cout << std::endl;
                tmpHashSeqPair++;
            }
            kmerCounter += kmerConsidered;
        }
        delete [] kmers;
    }
    return kmerCounter;
}


int kmermatcher(int argc, const char **argv, const Command &command) {
    Parameters& par = Parameters::getInstance();
    setLinearFilterDefault(&par);
    par.parseParameters(argc, argv, command, 2, true, false, MMseqsParameter::COMMAND_CLUSTLINEAR);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    struct timeval start, end;
    gettimeofday(&start, NULL);

    BaseMatrix * subMat;
    if(par.alphabetSize == 21){
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 8.0, -0.2);
    }else{
        SubstitutionMatrix sMat(par.scoringMatrixFile.c_str(), 8.0, -0.2);
        subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, par.alphabetSize, 8.0);
    }


    DBReader<unsigned int> seqDbr(par.db1.c_str(), (par.db1 + ".index").c_str());
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    seqDbr.readMmapedDataInMemory();
    const size_t KMER_SIZE = par.kmerSize;
    size_t chooseTopKmer = par.kmersPerSequence;
    DBWriter dbw(par.db2.c_str(), std::string(par.db2 + ".index").c_str(), 1);
    dbw.open();
    size_t totalMemoryInByte =  Util::getTotalSystemMemory();
    Debug(Debug::WARNING) << "Check requirements\n";
    size_t totalKmers = computeKmerCount(seqDbr, KMER_SIZE, chooseTopKmer);
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter(totalKmers);
    Debug(Debug::INFO) << "Needed memory (" << totalSizeNeeded << " byte) of total memory (" << totalMemoryInByte << " byte)\n";
    while(totalSizeNeeded > totalMemoryInByte){
        Debug(Debug::INFO) << "Adjust memory demand by reducing k-mers per sequence\n";
        chooseTopKmer = chooseTopKmer / 2;
        totalKmers = computeKmerCount(seqDbr, KMER_SIZE, chooseTopKmer);
        totalSizeNeeded = computeMemoryNeededLinearfilter(totalKmers);
        if(chooseTopKmer == 1 && totalSizeNeeded > totalMemoryInByte){
            Debug(Debug::ERROR) << "There is not enough memory on this machine\n";
            EXIT(EXIT_FAILURE);
        }
        Debug(Debug::INFO) << chooseTopKmer << " k-mer per sequence needs " << totalSizeNeeded << " byte of memory.\n";
    }
    Debug(Debug::WARNING) << "Generate k-mers list ... \n";
    KmerPosition * hashSeqPair = new KmerPosition[totalKmers+1];
#pragma omp parallel for
    for(size_t i = 0; i < totalKmers+1; i++){
        hashSeqPair[i].kmer = SIZE_T_MAX;
    }


    size_t kmerCounter = fillKmerPositionArray(hashSeqPair, seqDbr, par, subMat, KMER_SIZE, chooseTopKmer);
    Debug(Debug::WARNING) << "Done." << "\n";
    Debug(Debug::WARNING) << "Sort kmer ... ";
    omptl::sort(hashSeqPair, hashSeqPair + kmerCounter, KmerPosition::compareRepSequenceAndId);
    Debug(Debug::WARNING) << "Done." << "\n";
    // sort by hash and set up the pointer for parallel processing
    size_t prevHash = hashSeqPair[0].kmer;
    size_t prevHashStart = 0;
    size_t prevSetSize = 0;
    for(size_t pos = 0; pos < kmerCounter + 1; pos++) {
        if(prevHash !=  hashSeqPair[pos].kmer){
            for(size_t i = prevHashStart; i < pos; i++){
                hashSeqPair[i].size = (prevSetSize < SHRT_MAX) ? prevSetSize : SHRT_MAX;
            }
            prevSetSize = 0;
            prevHashStart = pos;
        }
        prevSetSize++;
        prevHash = hashSeqPair[pos].kmer;
    }
    Debug(Debug::WARNING) << "Sort again ... ";

    omptl::sort(hashSeqPair, hashSeqPair + kmerCounter, KmerPosition::compareKmerPositionByKmerAndSize);
    Debug(Debug::WARNING) << "Done\n";

    std::vector<bool> repSequence(seqDbr.getSize());
    std::fill(repSequence.begin(), repSequence.end(), false);

    // set rep. sequences for each kmer clusters
    size_t pos;
//    Indexer idxer(subMat->alphabetSize, KMER_SIZE);

    for(pos = 0; pos < kmerCounter; pos++) {
        if (hashSeqPair[pos].size <= 1) { // stop when there is just one member per cluster
            // speed up computation
            break;
        }

        size_t initKmer = hashSeqPair[pos].kmer;
//        idxer.printKmer(initKmer, KMER_SIZE, subMat->int2aa);
//        std::cout << std::endl;
        unsigned int queryId = UINT_MAX;
        unsigned short max_i_Pos = 0;
        size_t maxSeqLen = 0;
        const size_t startPos = pos;
        while (hashSeqPair[pos].kmer == initKmer) {
            const unsigned int id = hashSeqPair[pos].id;
            const unsigned int i_pos = hashSeqPair[pos].pos;
            size_t currSeqLen = seqDbr.getSeqLens(id);
            // if it is not yet a representative sequence and if it is the longest sequence
            if (currSeqLen > maxSeqLen ) {
                queryId = id;
                maxSeqLen = currSeqLen;
                max_i_Pos = i_pos;
            }
            pos++;
        }
        repSequence[queryId] = true;

        // set rep. sequence for all kmer hits
        pos = startPos;
        while (hashSeqPair[pos].kmer == initKmer) {
            const unsigned int j_pos = hashSeqPair[pos].pos;
            hashSeqPair[pos].kmer = queryId;
            short diagonal = max_i_Pos - j_pos;
            hashSeqPair[pos].pos = diagonal;
            pos++;
        }
        pos--;
    }

    // sort by rep. sequence (stored in kmer) and sequence id
    Debug(Debug::WARNING) << "Sort by rep. sequence ... ";
    omptl::sort(hashSeqPair, hashSeqPair + pos, KmerPosition::compareRepSequenceAndId);
    Debug(Debug::WARNING) << "Done\n";

    unsigned int repSeqId = UINT_MAX;
    unsigned int lastTargetId = UINT_MAX;
    unsigned int writeSets = 0;
    unsigned int queryLength = 0;
    // write result
    std::stringstream swResultsSs;
    for(size_t kmerPos = 0; kmerPos < pos; kmerPos++){
        if(repSeqId != hashSeqPair[kmerPos].kmer) {
            if (writeSets > 0) {
                repSequence[repSeqId] = true;
                std::string swResultsString = swResultsSs.str();
                const char *swResultsStringData = swResultsString.c_str();
                dbw.writeData(swResultsStringData, swResultsString.length(), seqDbr.getDbKey(repSeqId), 0);
            }else{
                if(repSeqId != UINT_MAX) {
                    repSequence[repSeqId] = false;
                }
            }
            lastTargetId = UINT_MAX;
            swResultsSs.str(std::string());
            repSeqId = hashSeqPair[kmerPos].kmer;
            queryLength = std::max(seqDbr.getSeqLens(repSeqId), 3ul) - 2;
            swResultsSs << seqDbr.getDbKey(repSeqId) << "\t";
            swResultsSs << 0 << "\t";
            swResultsSs << 0 << "\n";
        }
        unsigned int targetId = hashSeqPair[kmerPos].id;
        unsigned short diagonal = hashSeqPair[kmerPos].pos;
        unsigned int targetLength = std::max(seqDbr.getSeqLens(targetId), 3ul) - 2;
        if(targetId != repSeqId && lastTargetId != targetId ){
            if ( (((float) queryLength) / ((float) targetLength) < par.covThr) ||
                 (((float) targetLength) / ((float) queryLength) < par.covThr) ) {
                //                    foundAlready[targetId] = 0;
                lastTargetId = targetId;
                continue;
            }
        }else{
            lastTargetId = targetId;
            continue;
        }
        swResultsSs << seqDbr.getDbKey(targetId) << "\t";
        swResultsSs << 0 << "\t";
        swResultsSs << static_cast<short>(diagonal) << "\n";
        lastTargetId = targetId;
        writeSets++;
    }
    if (writeSets > 0) {
        repSequence[repSeqId] = true;
        std::string swResultsString = swResultsSs.str();
        const char *swResultsStringData = swResultsString.c_str();
        dbw.writeData(swResultsStringData, swResultsString.length(), seqDbr.getDbKey(repSeqId), 0);
    }else{
        if(repSeqId != UINT_MAX) {
            repSequence[repSeqId] = false;
        }
    }
    // add missing entries to the result (needed for clustering)
    for(size_t id = 0; id < seqDbr.getSize(); id++){
        if(repSequence[id] == false){
            std::stringstream swResultsSs;
            swResultsSs << seqDbr.getDbKey(id) << "\t";
            swResultsSs << 0 << "\t";
            swResultsSs << 0 << "\n";
            std::string swResultsString = swResultsSs.str();
            const char* swResultsStringData = swResultsString.c_str();
            dbw.writeData(swResultsStringData, swResultsString.length(), seqDbr.getDbKey(id), 0);
        }
    }
    // free memory
    delete subMat;
    delete [] hashSeqPair;
    seqDbr.close();
    dbw.close();

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for processing: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    return 0;
}

#undef SIZE_T_MAX
