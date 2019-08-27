//
// Created by mad on 2/25/16.
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
#include "DistanceCalculator.h"

#ifdef OPENMP
#include <omp.h>
#endif

void setClustHashDefaults(Parameters *p) {
    p->alphabetSize = Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE;

}

int clusthash(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setClustHashDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, -0.2);
    ReducedMatrix redSubMat(subMat.probMatrix, subMat.subMatrixPseudoCounts, subMat.aa2int, subMat.int2aa, subMat.alphabetSize, par.alphabetSize, 2.0);

    DBReader<unsigned int> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    seqDbr.open(DBReader<unsigned int>::NOSORT);
    seqDbr.readMmapedDataInMemory();

    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_PREFILTER_RES);
    dbw.open();
    Debug(Debug::INFO) << "Hashing sequences ... \n";
    std::pair<size_t, unsigned int> * hashSeqPair = new  std::pair<size_t, unsigned int>[seqDbr.getSize()+1];
    hashSeqPair[seqDbr.getSize()] = std::make_pair(UINT_MAX, 0); // needed later to check if one of array
    Debug::Progress progress(seqDbr.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence seq(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, &redSubMat, 0, false, false);
#pragma omp for schedule(dynamic, 10000)
        for(size_t id = 0; id < seqDbr.getSize(); id++){
            progress.updateProgress();
            unsigned int queryKey = seqDbr.getDbKey(id);
            char * data = seqDbr.getData(id, thread_idx);
            seq.mapSequence(id, queryKey, data, seqDbr.getSeqLen(id));
            size_t seqHash = Util::hash(seq.int_sequence, seq.L);
            hashSeqPair[id] = std::make_pair(seqHash, id);
        }
    }



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
    Debug(Debug::INFO) << "Compute " << uniqHashes << " unique hashes.\n";

#pragma omp parallel
    {
        std::vector<unsigned int> setIds;
        std::vector<bool> found;

#pragma omp for schedule(dynamic, 2)
        for(size_t hashId = 0; hashId < uniqHashes; hashId++) {
            size_t initHash = hashLookup[hashId]->first;
            size_t pos = 0;
            progress.updateProgress();

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
                unsigned int queryLength = seqDbr.getSeqLen(setIds[i]);
                const char * querySeq =  seqDbr.getData(setIds[i], thread_idx);
                std::stringstream swResultsSs;
                swResultsSs << seqDbr.getDbKey(setIds[i]) << "\t";
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
                    unsigned int targetLength = seqDbr.getSeqLen(setIds[j]);
                    if(i != j && queryLength == targetLength){
                        const char * targetSeq = seqDbr.getData(setIds[j], thread_idx);
                        unsigned int distance = DistanceCalculator::computeInverseHammingDistance(querySeq, targetSeq,
                                                                                                  queryLength);
                        float seqId = (static_cast<float>(distance))/static_cast<float>(queryLength);
                        if(seqId >= par.seqIdThr) {
                            swResultsSs << seqDbr.getDbKey(setIds[j]) << "\t";
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
                dbw.writeData(swResultsStringData, swResultsString.length(), seqDbr.getDbKey(setIds[i]), thread_idx);
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
