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
#include "Log.h"
#include "SubstitutionMatrix.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "omptl/omptl_algorithm"
#ifdef OPENMP
#include <omp.h>
#endif

struct KmerPosition {
    size_t kmer;
    unsigned int id;
    unsigned int size;
};

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


int linearfilter (int argc, const char * argv[])
{
    std::string usage;
    usage.append("Detects redundant sequences based on reduced alphabet and k-mer sorting.\n");
    usage.append("USAGE: <sequenceDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>.\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.linearfilter, 2);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    BaseMatrix * subMat;
    if(par.alphabetSize == 21){
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    }else{
        SubstitutionMatrix sMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
        subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, par.alphabetSize, 2.0);
    }

    DBReader<unsigned int> seqDbr(par.db1.c_str(), (par.db1 + ".index").c_str());
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    seqDbr.readMmapedDataInMemory();

//  Sensitivity: For good sensitivity we will need a reduced alphabet, e.g. 8 states,
//  while for maintaining good specificity we need a sufficient k-mer length, as shown now.
//  Say in reduced alphabet the sequence identity is seqid, e.g. 90% for 70% sequence identity in original alphabet.
//  Then the probability for two homologous k-mers to be identical is 1 - (1-seqid)^k.
//  The probability that two full-length homologous sequences of length L will have a common k-mer is
//  therefore approximately 1 - (1 - (1-seqid)^k)^L. This number should be near 1.

//  Specificity: Each sequence should have a probability much less than 1
//  that it will have a chance k-mer match with a non-homologous sequence in the database:
//  L * NL * (1/A)^k «1. Let's say we demand a probability 0.001, then we get k ≥ 56 / log2(A).
//  At A=8 this gives k=19, at A=4 it gives 28.
//  That is probably already too long since spacing of the pattern will require no gap within ~50 residues.
//  That will reduce sensitivity.
//  Combined
//  These two conditions can be combined into 56 / log2(A) ≤ k « log2(1/L)/Log2(1-seqid).
//  The best sensitivity at sufficient specificity can be obtained by choosing the alphabet size with
//  largest log2(A) / |log2 (1-seqid)| and then choosing the k = 56 / log2(A).
    const size_t KMER_SIZE = par.kmerSize;
    size_t chooseTopKmer = 10;
    DBWriter dbw(par.db2.c_str(), std::string(par.db2 + ".index").c_str(), 1);
    dbw.open();
    Debug(Debug::WARNING) << "Generate k-mers list ... \n";
    size_t totalKmers = 0;
    for(size_t id = 0; id < seqDbr.getSize(); id++ ){
        int kmerAdjustedSeqLen = std::max(0, static_cast<int>(seqDbr.getSeqLens(id)) - static_cast<int>(KMER_SIZE - 1)-2) ;
        totalKmers += std::min(kmerAdjustedSeqLen, static_cast<int>( chooseTopKmer ) );
    }
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
            Log::printProgress(id);
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

    Debug(Debug::WARNING) << "Kmer considered " << kmerCounter << "\n";
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

    //DEBUG code
//    Indexer idxer(subMat->alphabetSize, KMER_SIZE);
//    for(size_t id = 0; id < kmerCounter; id++){
//        idxer.printKmer(hashSeqPair[id].kmer, KMER_SIZE, subMat->int2aa);
//
//        Debug(Debug::WARNING) << "\t" << hashSeqPair[id].id;
//        Debug(Debug::WARNING) << "\t" << hashSeqPair[id].size;
//        Debug(Debug::WARNING) << "\n";
//    }


    char * foundAlready = new char[seqDbr.getSize()];
    memset(foundAlready, 0, seqDbr.getSize() * sizeof(char));

    for(size_t pos = 0; pos < kmerCounter; pos++) {
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
        for(size_t i = 0; i < setIds.size(); i++) {
            unsigned int targetLength = std::max(seqDbr.getSeqLens(setIds[i]), 3ul) - 2;
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
        }
        if(setIds.size() > 0){
            foundAlready[queryId] = 2;
            std::string swResultsString = swResultsSs.str();
            const char* swResultsStringData = swResultsString.c_str();
            dbw.write(swResultsStringData, swResultsString.length(), SSTR(seqDbr.getDbKey(queryId)).c_str(), 0);
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
            dbw.write(swResultsStringData, swResultsString.length(), SSTR(seqDbr.getDbKey(id)).c_str(), 0);
        }
    }

    delete subMat;
    delete [] foundAlready;
    delete [] hashSeqPair;
    seqDbr.close();
    dbw.close();
    return 0;
}
