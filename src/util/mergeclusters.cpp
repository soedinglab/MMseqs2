#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "itoa.h"

#include <list>

#ifdef OPENMP
#include <omp.h>
#endif

void mergeClusteringResults(std::string seqDB, std::string outDB, std::list<std::string> cluSteps, int threads, int compressed){
    // open the sequence database
    // it will serve as the reference for sequence indexes
    std::string seqDBIndex = seqDB + ".index";
    DBReader<unsigned int> dbr(seqDB.c_str(), seqDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::NOSORT);

    // init the structure for cluster merging
    // it has the size of all possible cluster (sequence amount)
    std::list<unsigned int>** mergedClustering = new std::list<unsigned int>*[dbr.getSize()];
    Debug(Debug::INFO) << "List amount "<< dbr.getSize() << "\n";
#pragma omp parallel for
    for (size_t i = 0; i < dbr.getSize(); i++){
        mergedClustering[i] = new std::list<unsigned int>();
    }

    // read the clustering from the first clustering step
    std::string firstCluStep = cluSteps.front();
    cluSteps.pop_front();
    std::string firstCluStepIndex = firstCluStep + ".index";
    DBReader<unsigned int>* cluStepDbr = new DBReader<unsigned int>(firstCluStep.c_str(), firstCluStepIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    cluStepDbr->open(DBReader<unsigned int>::NOSORT);
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
#pragma omp for
        for (size_t i = 0; i < cluStepDbr->getSize(); i++){
            unsigned int clusterId = cluStepDbr->getDbKey(i);
            size_t cluId = dbr.getId(clusterId);
            char * data = cluStepDbr->getData(i, thread_idx);
            // go through the sequences in the cluster and add them to the initial clustering
            char keyBuffer[255];
            while (*data != '\0'){
                Util::parseKey(data, keyBuffer);
                unsigned int key = Util::fast_atoi<unsigned int>(keyBuffer);
                size_t seqId = dbr.getId(key);
                mergedClustering[cluId]->push_back(seqId);
                data = Util::skipLine(data);
            }
        }
    };
    cluStepDbr->close();
    delete cluStepDbr;
    Debug(Debug::INFO) << "Clustering step 1\n";

    // merge later clustering steps into the initial clustering step
    int cnt = 2;
    while(!cluSteps.empty()){
        // open the next clustering database
        std::string cluStep = cluSteps.front();
        std::string cluStepIndex = cluStep + ".index";
        cluSteps.pop_front();

        cluStepDbr = new DBReader<unsigned int>(cluStep.c_str(), cluStepIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        cluStepDbr->open(DBReader<unsigned int>::NOSORT);

        // go through the clusters and merge them into the clusters from the previous clustering step
#pragma omp parallel
        {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
#pragma omp for
            for (size_t i = 0; i < cluStepDbr->getSize(); i++){
                size_t cluId = dbr.getId(cluStepDbr->getDbKey(i));
                char* cluData = cluStepDbr->getData(i, thread_idx);
                // go through the sequences in the cluster and add them and their clusters to the cluster of cluId
                // afterwards, delete the added cluster from the clustering
                char * data = cluData;
                char keyBuffer[255];
                while (*data != '\0') {
                    Util::parseKey(data, keyBuffer);
                    unsigned int key = Util::fast_atoi<unsigned int>(keyBuffer);
                    size_t seqId = dbr.getId(key);
                    if(seqId != cluId) { // to avoid copies of the same cluster list
                        mergedClustering[cluId]->splice(mergedClustering[cluId]->end(), *mergedClustering[seqId]);
                    }
                    data = Util::skipLine(data);
                }
            }
        }
        cluStepDbr->close();
        delete cluStepDbr;
        Debug(Debug::INFO) << "Clustering step " << cnt << "\n";
        cnt++;
    }

    Debug(Debug::INFO) << "Writing the results\n";

    std::string outDBIndex = outDB + ".index";
    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), threads, compressed, Parameters::DBTYPE_CLUSTER_RES);
    dbw->open();

#pragma omp parallel
    {
        std::string res;
        res.reserve(1024*1024);

        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        // go through all sequences in the database
#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < dbr.getSize(); i++){
            // no cluster for this representative
            if (mergedClustering[i]->size() == 0)
                continue;

            // representative
            unsigned int dbKey = dbr.getDbKey(i);
            char buffer[32];
            for(std::list<unsigned int>::iterator it = mergedClustering[i]->begin();
                it != mergedClustering[i]->end(); ++it){
                char * tmpBuff = Itoa::u32toa_sse2(dbr.getDbKey(*it), buffer);
                size_t length = tmpBuff - buffer - 1;
                res.append(buffer, length);
                res.push_back('\n');
            }

            dbw->writeData(res.c_str(), res.length(), dbKey, thread_idx);
            res.clear();
        }
    }
    dbw->close();
    delete dbw;

    // delete the clustering data structure
    for (unsigned int i = 0; i < dbr.getSize(); i++){
        delete mergedClustering[i];
    }
    delete[] mergedClustering;
}

int mergeclusters(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3, true, true);

    std::list<std::string> clusterings;
    for (size_t i = 2; i < par.filenames.size(); i++) {
        clusterings.push_back(par.filenames[i]);
    }

    mergeClusteringResults(par.db1, par.db2, clusterings, par.threads, par.compressed);

    return 0;
}
