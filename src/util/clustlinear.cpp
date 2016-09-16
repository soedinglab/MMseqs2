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

#ifdef OPENMP
#include <omp.h>
#endif

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif
struct KmerPosition {
    size_t kmer;
    unsigned int id;
    unsigned int size;
};

void setLinearFilterDefault(Parameters *p) {
    p->spacedKmer = false;
    p->covThr = 0.8;
    p->kmerSize = 14;
    p->alphabetSize = 14;
    p->kmersPerSequence = 20;
}

bool compareKmerPositionByKmer(KmerPosition first, KmerPosition second){
    return (first.kmer < second.kmer) ? true : false;
}

float computeMutualInformation(int x, BaseMatrix *subMat){
    float mutalInformation=0.0;
    for(int aa = 0; aa < subMat->alphabetSize; aa++){
        mutalInformation += subMat->probMatrix[x][aa] * subMat->subMatrix[x][aa];
    }
    return mutalInformation;
}

bool compareByMutalInformation(std::pair<float, size_t > first, std::pair<float, size_t > second){
    return (first.first > second.first) ? true : false;
}

bool compareKmerPositionByKmerAndSize(KmerPosition first, KmerPosition second){
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

int clustlinear(int argc, const char **argv, const Command& command) {
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
    for(size_t i = 0; i < totalKmers+1; i++){
        hashSeqPair[i].kmer = SIZE_T_MAX;
    }
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
        std::pair<float, size_t > * kmers = new std::pair<float, size_t >[par.maxSeqLen];
        for(size_t id = dbFrom; id < (dbFrom + dbSize); id++){
            Debug::printProgress(id);
            seq.mapSequence(id, id, seqDbr.getData(id));
            int seqKmerCount = 0;
            unsigned int seqId = seq.getId();
            while(seq.hasNextKmer()) {
                const int *kmer = seq.nextKmer();
                float kmerMutualInformation = 0.0;
                for (size_t kpos = 0; kpos < KMER_SIZE; kpos++) {
                    kmerMutualInformation += computeMutualInformation(kmer[kpos], subMat);
                }
                (kmers+seqKmerCount)->first = kmerMutualInformation;
                size_t kmerIdx = idxer.int2index(kmer, 0, KMER_SIZE);

                (kmers+seqKmerCount)->second = kmerIdx;
                seqKmerCount++;
            }
            if(seqKmerCount > 0){
                std::sort(kmers, kmers+seqKmerCount, compareByMutalInformation);
            }
            size_t kmerConsidered = std::min(static_cast<int>(chooseTopKmer), seqKmerCount);

            for(size_t topKmer = 0; topKmer < kmerConsidered; topKmer++){
                tmpHashSeqPair->kmer  = (kmers + topKmer)->second;
                tmpHashSeqPair->id = seqId;
                tmpHashSeqPair++;
            }
            kmerCounter += kmerConsidered;
        }
        delete [] kmers;
    }
    Debug(Debug::WARNING) << "Done." << "\n";
    if(totalKmers != kmerCounter ){
        Debug(Debug::WARNING) << "Problem totalKmers(" << totalKmers << ") != kmerCounter("<<kmerCounter << ") \n";
    }
    Debug(Debug::WARNING) << "Sort kmer ... ";
    omptl::sort(hashSeqPair, hashSeqPair + kmerCounter, compareKmerPositionByKmer);
    Debug(Debug::WARNING) << "Done." << "\n";
    // sort by hash and set up the pointer for parallel processing
    size_t uniqHashes = 1;
    size_t prevHash = hashSeqPair[0].kmer;
    size_t prevHashStart = 0;
    size_t prevSetSize = 0;
    for(size_t pos = 0; pos < kmerCounter + 1; pos++) {
        if(prevHash !=  hashSeqPair[pos].kmer){
            uniqHashes++;
            for(size_t i = prevHashStart; i < pos; i++){
                hashSeqPair[i].size = prevSetSize;
            }
            prevSetSize = 0;
            prevHashStart = pos;
        }
        prevSetSize++;
        prevHash = hashSeqPair[pos].kmer;
    }
    Debug(Debug::WARNING) << "Uniq kmer " << uniqHashes << "\n";
    Debug(Debug::WARNING) << "Sort again ... ";

    omptl::sort(hashSeqPair, hashSeqPair + kmerCounter, compareKmerPositionByKmerAndSize);
    Debug(Debug::WARNING) << "Done\n";

    char * foundAlready = new char[seqDbr.getSize()];
    memset(foundAlready, 0, seqDbr.getSize() * sizeof(char));

    for(size_t pos = 0; pos < kmerCounter; pos++) {
        if(hashSeqPair[pos].size <= 1){ // stop when there is just one member per cluster
                                        // speed up computation
            break;
        }

        size_t initKmer = hashSeqPair[pos].kmer;
        std::vector<unsigned int> setIds;
        unsigned int queryId = 0;
        size_t maxSeqLen = 0;
        while(hashSeqPair[pos].kmer == initKmer){
            const unsigned int id = hashSeqPair[pos].id;
            if(foundAlready[id] == 0){
                setIds.push_back(id);
                foundAlready[id] = 1;
                size_t currSeqLen = seqDbr.getSeqLens(id);
                if( currSeqLen > maxSeqLen ){
                    queryId = id;
                    maxSeqLen = currSeqLen;
                }
            }
            pos++;
        }
        pos--;
        
        std::stringstream swResultsSs;
        unsigned int queryLength = std::max(seqDbr.getSeqLens(queryId), 3ul) - 2;
        unsigned int writeSets = 0;
        for(size_t i = 0; i < setIds.size(); i++) {
            unsigned int targetLength = std::max(seqDbr.getSeqLens(setIds[i]), 3ul) - 2;
            if(par.fragmentMerge == false && setIds[i] != queryId ){
                if ( (((float) queryLength) / ((float) targetLength) < par.covThr) ||
                     (((float) targetLength) / ((float) queryLength) < par.covThr) ) {
                    foundAlready[setIds[i]] = 0;
                    continue;
                }
            }
            swResultsSs << SSTR(seqDbr.getDbKey(setIds[i])).c_str() << "\t";
            swResultsSs << 255 << "\t";
            swResultsSs << std::fixed << std::setprecision(3) << 1.0f << "\t";
            swResultsSs << std::scientific << 0 << "\t";
            swResultsSs << 0 << "\t";
            swResultsSs << queryLength - 1 << "\t";
            swResultsSs << queryLength << "\t";
            swResultsSs << 0 << "\t";
            swResultsSs << targetLength - 1 << "\t";
            swResultsSs << targetLength << "\n";
            writeSets++;
        }
        if(writeSets > 0){
            foundAlready[queryId] = 2;
            std::string swResultsString = swResultsSs.str();
            const char* swResultsStringData = swResultsString.c_str();
            dbw.writeData(swResultsStringData, swResultsString.length(), SSTR(seqDbr.getDbKey(queryId)).c_str(), 0);
        }
        setIds.clear();
        swResultsSs.clear();
    }
    // add missing entries to the result (needed for clustering)
    for(size_t id = 0; id < seqDbr.getSize(); id++){
        if(foundAlready[id] != 2){
            std::stringstream swResultsSs;
            unsigned int queryLength = std::max(seqDbr.getSeqLens(id), 3ul) - 2;
            swResultsSs << SSTR(seqDbr.getDbKey(id)).c_str() << "\t";
            swResultsSs << 255 << "\t";
            swResultsSs << std::fixed << std::setprecision(3) << 1.0f << "\t";
            swResultsSs << std::scientific << 0 << "\t";
            swResultsSs << 0 << "\t";
            swResultsSs << queryLength - 1 << "\t";
            swResultsSs << queryLength << "\t";
            swResultsSs << 0 << "\t";
            swResultsSs << queryLength - 1 << "\t";
            swResultsSs << queryLength << "\n";
            std::string swResultsString = swResultsSs.str();
            const char* swResultsStringData = swResultsString.c_str();
            dbw.writeData(swResultsStringData, swResultsString.length(), SSTR(seqDbr.getDbKey(id)).c_str(), 0);
        }
    }
    // free memory
    delete subMat;
    delete [] foundAlready;
    delete [] hashSeqPair;
    seqDbr.close();
    dbw.close();

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for processing: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    return 0;
}

#undef SIZE_T_MAX
