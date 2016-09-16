//
// Created by mad on 2/25/16.
//

//
// Created by mad on 10/21/15.
//
#include <limits>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

#include "ReducedMatrix.h"
#include "DBWriter.h"
#include "SubstitutionMatrix.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

unsigned int computeHammingDistance(const char *seq1, const char *seq2, unsigned int length){
    unsigned int diff = 0;
    for (unsigned int pos = 0; pos < length; pos++ ) {
        diff += (seq1[pos] != seq2[pos]);
    }
    return diff;
}

size_t hash(int * x, size_t length){
    const size_t INITIAL_VALUE = 0;
    const size_t A = 31;
    size_t h = INITIAL_VALUE;
    for (size_t i = 0; i < length; ++i){
        h = ((h*A) + x[i]);
    }
    return h;
}

void setClustHashDefaults(Parameters *p) {
    p->alphabetSize = Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE;

}

int clusthash(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setClustHashDefaults(&par);
    par.parseParameters(argc, argv, command, 2);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, -0.2);
    ReducedMatrix redSubMat(subMat.probMatrix, subMat.subMatrixPseudoCounts, par.alphabetSize, 2.0);

    DBReader<unsigned int> seqDbr(par.db1.c_str(), (par.db1 + ".index").c_str());
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    seqDbr.readMmapedDataInMemory();

    DBWriter dbw(par.db2.c_str(), std::string(par.db2 + ".index").c_str(), par.threads);
    dbw.open();
    Debug(Debug::WARNING) << "Hashing sequences ... \n";
    std::pair<size_t, unsigned int> * hashSeqPair = new  std::pair<size_t, unsigned int>[seqDbr.getSize()+1];
    hashSeqPair[seqDbr.getSize()] = std::make_pair(UINT_MAX, 0); // needed later to check if one of array
#pragma omp parallel
    {
        Sequence seq(par.maxSeqLen, redSubMat.aa2int, redSubMat.int2aa, Sequence::AMINO_ACIDS, 0, false, false);
#pragma omp for schedule(dynamic, 10000)
        for(size_t id = 0; id < seqDbr.getSize(); id++){
            Debug::printProgress(id);
            unsigned int queryKey = seqDbr.getDbKey(id);
            char * data = seqDbr.getData(id);
            seq.mapSequence(id, queryKey, data);
            size_t seqHash = hash(seq.int_sequence, seq.L);
            hashSeqPair[id] = std::make_pair(seqHash, id);
        }
    }
    Debug(Debug::WARNING) << "Done." << "\n";


    // sort by hash and set up the pointer for parallel processing
    std::sort(hashSeqPair, hashSeqPair + seqDbr.getSize());
    size_t uniqHashes = 1;
    size_t prevHash = hashSeqPair[0].first;
    for(size_t id = 0; id < seqDbr.getSize(); id++) {
        if(prevHash !=  hashSeqPair[id].first){
            uniqHashes++;
        }
        prevHash =  hashSeqPair[id].first;
    }
    std::pair<size_t, unsigned int> ** hashLookup = new std::pair<size_t, unsigned int> *[uniqHashes];
    hashLookup[0] = hashSeqPair;
    size_t currKey = 1;
    prevHash = hashSeqPair[0].first;
    for(size_t id = 0; id < seqDbr.getSize(); id++) {
        if (prevHash != hashSeqPair[id].first) {
            hashLookup[currKey] = (hashSeqPair + id);
            currKey++;
        }
        prevHash = hashSeqPair[id].first;
    }
    Debug(Debug::WARNING) << "Compute "<< uniqHashes <<" unique hashes.\n";

#pragma omp parallel
    {
        std::vector<unsigned int> setIds;
        std::vector<bool> found;

#pragma omp for schedule(dynamic, 2)
        for(size_t hashId = 0; hashId < uniqHashes; hashId++) {
            size_t initHash = hashLookup[hashId]->first;
            size_t pos = 0;
            Debug::printProgress(hashId);

            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            while(hashLookup[hashId][pos].first == initHash ){
                setIds.push_back(hashLookup[hashId][pos].second);
                found.push_back(false);
                pos++;
            }
            for(size_t i = 0; i < setIds.size(); i++) {
                unsigned int queryLength = std::max(seqDbr.getSeqLens(setIds[i]), 3ul) - 2;
                const char * querySeq =  seqDbr.getData(setIds[i]);
                std::stringstream swResultsSs;
                swResultsSs << SSTR(seqDbr.getDbKey(setIds[i])).c_str() << "\t";
                swResultsSs << 255 << "\t";
                swResultsSs << std::fixed << std::setprecision(3) << 1.0f << "\t";
                swResultsSs << std::scientific << 0 << "\t";
                swResultsSs << 0 << "\t";
                swResultsSs << queryLength - 1 << "\t";
                swResultsSs << queryLength << "\t";
                swResultsSs << 0 << "\t";
                swResultsSs << queryLength - 1 << "\t";
                swResultsSs << queryLength << "\n";
                if(found[i] == true){
                    goto outer;
                }

                for (size_t j = 0; j < setIds.size(); j++) {
                    if(found[j] == true)
                        continue;
                    unsigned int targetLength = std::max(seqDbr.getSeqLens(setIds[j]), 3ul) - 2;
                    if(i != j && queryLength == targetLength){
                        const char * targetSeq = seqDbr.getData(setIds[j]);
                        unsigned int distance = computeHammingDistance(querySeq, targetSeq, queryLength);
                        float seqId = (static_cast<float>(queryLength) - static_cast<float>(distance))/static_cast<float>(queryLength);
                        if(seqId > par.seqIdThr) {
                            swResultsSs << SSTR(seqDbr.getDbKey(setIds[j])).c_str() << "\t";
                            swResultsSs << 255 << "\t";
                            swResultsSs << std::fixed << std::setprecision(3) << seqId << "\t";
                            swResultsSs << std::scientific << 0 << "\t";
                            swResultsSs << 0 << "\t";
                            swResultsSs << queryLength - 1 << "\t";
                            swResultsSs << queryLength << "\t";
                            swResultsSs << 0 << "\t";
                            swResultsSs << queryLength - 1 << "\t";
                            swResultsSs << queryLength << "\n";
                            found[j] = true;
                        }
                    }
                }
                outer:
                std::string swResultsString = swResultsSs.str();
                const char* swResultsStringData = swResultsString.c_str();
                dbw.writeData(swResultsStringData, swResultsString.length(), SSTR(seqDbr.getDbKey(setIds[i])).c_str(),
                              thread_idx);
            }
            setIds.clear();
            found.clear();
        }
    }
    delete [] hashLookup;
    delete [] hashSeqPair;
    seqDbr.close();
    dbw.close();
    return 0;
}
