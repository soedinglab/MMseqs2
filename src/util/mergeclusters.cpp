#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "itoa.h"

#include <list>

#ifdef OPENMP
#include <omp.h>
#endif

int mergeclusters(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, true, 0);

    std::list<std::string> clusterings;
    for (size_t i = 2; i < par.filenames.size(); i++) {
        clusterings.push_back(par.filenames[i]);
    }

    // the sequence database will serve as the reference for sequence indexes
    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    dbr.open(DBReader<unsigned int>::NOSORT);

    // init the structure for cluster merging
    // it has the size of all possible cluster (sequence amount)
    std::list<unsigned int> *mergedClustering = new std::list<unsigned int>[dbr.getSize()];

    // read the clustering from the first clustering step
    std::string firstClu = clusterings.front();
    std::string firstCluStepIndex = firstClu + ".index";
    clusterings.pop_front();

    Debug(Debug::INFO) << "Clustering step 1\n";
    DBReader<unsigned int> cluDb(firstClu.c_str(), firstCluStepIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    cluDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Debug::Progress progress(cluDb.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        char keyBuffer[255];
#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < cluDb.getSize(); i++) {
            progress.updateProgress();
            unsigned int clusterId = cluDb.getDbKey(i);
            size_t cluId = dbr.getId(clusterId);
            char *data = cluDb.getData(i, thread_idx);
            // go through the sequences in the cluster and add them to the initial clustering
            while (*data != '\0') {
                Util::parseKey(data, keyBuffer);
                unsigned int key = Util::fast_atoi<unsigned int>(keyBuffer);
                size_t seqId = dbr.getId(key);
                mergedClustering[cluId].push_back(seqId);
                data = Util::skipLine(data);
            }
        }
    }
    cluDb.close();

    // merge later clustering steps into the initial clustering step
    int cnt = 2;
    while (!clusterings.empty()) {
        Debug(Debug::INFO) << "Clustering step " << cnt << "\n";

        std::string cluStep = clusterings.front();
        std::string cluStepIndex = cluStep + ".index";
        clusterings.pop_front();

        DBReader<unsigned int> cluDb(cluStep.c_str(), cluStepIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        cluDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

        progress.reset(cluDb.getSize());
        // go through the clusters and merge them into the clusters from the previous clustering step
#pragma omp parallel
        {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char keyBuffer[255];
#pragma omp for schedule(dynamic, 100)
            for (size_t i = 0; i < cluDb.getSize(); i++) {
                progress.updateProgress();
                // go through the sequences in the cluster and add them and their clusters to the cluster of cluId
                // afterwards, delete the added cluster from the clustering
                size_t cluId = dbr.getId(cluDb.getDbKey(i));
                char *data = cluDb.getData(i, thread_idx);
                while (*data != '\0') {
                    Util::parseKey(data, keyBuffer);
                    unsigned int key = Util::fast_atoi<unsigned int>(keyBuffer);
                    size_t seqId = dbr.getId(key);
                    if (seqId != cluId) { // to avoid copies of the same cluster list
                        mergedClustering[cluId].splice(mergedClustering[cluId].end(), mergedClustering[seqId]);
                    }
                    data = Util::skipLine(data);
                }
            }
        }
        cluDb.close();
        cnt++;
    }

    Debug(Debug::INFO) << "Write merged clustering\n";
    DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_CLUSTER_RES);
    dbw.open();
    progress.reset(dbr.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::string res;
        res.reserve(1024 * 1024);

        char buffer[32];

        // go through all sequences in the database
#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < dbr.getSize(); i++) {
            progress.updateProgress();

            // no cluster for this representative
            if (mergedClustering[i].empty())
                continue;

            // representative
            unsigned int dbKey = dbr.getDbKey(i);
            for (std::list<unsigned int>::iterator it = mergedClustering[i].begin(); it != mergedClustering[i].end(); ++it) {
                char *tmpBuff = Itoa::u32toa_sse2(dbr.getDbKey(*it), buffer);
                size_t length = tmpBuff - buffer - 1;
                res.append(buffer, length);
                res.push_back('\n');
            }

            dbw.writeData(res.c_str(), res.length(), dbKey, thread_idx);
            res.clear();
        }
    }
    dbw.close();
    dbr.close();

    delete[] mergedClustering;

    return EXIT_SUCCESS;
}
