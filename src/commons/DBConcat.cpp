#include "DBConcat.h"
#include "DBWriter.h"

#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Parameters.h"

#include <algorithm>
#include <sys/time.h>

#ifdef OPENMP
#include <omp.h>
#endif

DBConcat::DBConcat(const std::string &dataFileNameA, const std::string &indexFileNameA,
                   const std::string &dataFileNameB, const std::string &indexFileNameB,
                   const std::string &dataFileNameC, const std::string &indexFileNameC,
                   unsigned int threads, int dataMode, bool preserveKeysA, bool preserveKeysB)
        : DBReader((dataFileNameA == dataFileNameB ? dataFileNameA : dataFileNameC).c_str(),
                   (indexFileNameA == indexFileNameB ? indexFileNameA : indexFileNameC).c_str(), dataMode),
          dataFileNameA(dataFileNameA), indexFileNameA(indexFileNameA),
          dataFileNameB(dataFileNameB), indexFileNameB(indexFileNameB),
          dataFileNameC(dataFileNameC), indexFileNameC(indexFileNameC),
          threads(threads), preserveKeysA(preserveKeysA),preserveKeysB(preserveKeysB) {
    sameDatabase = dataFileNameA == dataFileNameB;
}

// If dbA != dbB, then Concatenate dbA and dbB in concatWriter ("dataFileNameC")
// and "this" will be a reader on "dataFileNameC" after calling open()
// otherwise, do nothing and "this"  will be a reader on "dataFileNameA"
void DBConcat::concat(bool write) {
    if (sameDatabase) {
        return;
    }

    DBReader<unsigned int> dbA(dataFileNameA.c_str(), indexFileNameA.c_str());
    DBReader<unsigned int> dbB(dataFileNameB.c_str(), indexFileNameB.c_str());


    dbA.open(DBReader<unsigned int>::NOSORT);
    dbB.open(DBReader<unsigned int>::NOSORT);

    indexSizeA = dbA.getSize();
    indexSizeB = dbB.getSize();

    // keys paris are like : (key,i) where key is the ith key in the ffindex
    keysA = new std::pair<unsigned int, unsigned int>[indexSizeA];
    keysB = new std::pair<unsigned int, unsigned int>[indexSizeB];

    DBWriter* concatWriter = NULL;
    if (write) {
        concatWriter = new DBWriter(dataFileNameC.c_str(), indexFileNameC.c_str(), threads, DBWriter::BINARY_MODE);
        concatWriter->open();
    }

    // where the new key numbering of B should start
    unsigned int maxKeyA = 0;
#pragma omp parallel for schedule(static) num_threads(threads) reduction(max:maxKeyA)
    for (size_t id = 0; id < indexSizeA; id++) {
        Debug::printProgress(id);
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        unsigned int newKey;

        if (preserveKeysA) {
            newKey = dbA.getDbKey(id);
        } else {
            newKey = static_cast<unsigned int>(id);
        }

        if (write) {
            char *data = dbA.getData(id);
            concatWriter->writeData(data, dbA.getSeqLens(id) - 1, SSTR(newKey).c_str(), thread_idx);
        }

        // need to store the index, because it'll be sorted out by keys later
        keysA[id] = std::make_pair(dbA.getDbKey(id), newKey);
        maxKeyA = std::max(maxKeyA, newKey);
    }
    maxKeyA++;

#pragma omp parallel for schedule(static) num_threads(threads)
    for (size_t id = 0; id < indexSizeB; id++) {
        Debug::printProgress(id);
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        unsigned int newKey;
        
        if (preserveKeysB) {
            newKey = dbB.getDbKey(id);
        } else {
            newKey = static_cast<unsigned int>(id) + maxKeyA;
        }
        
        if (write) {
            char *data = dbB.getData(id);
            concatWriter->writeData(data, dbB.getSeqLens(id) - 1, SSTR(newKey).c_str(), thread_idx);
        }

        // need to store the index, because it'll be sorted out by keys later
        keysB[id] = std::make_pair(dbB.getDbKey(id), id + maxKeyA);
    }

    //sort by key
    std::stable_sort(keysA, keysA + indexSizeA, compareFirstEntry());
    std::stable_sort(keysB, keysB + indexSizeB, compareFirstEntry());

    if (write) {
        concatWriter->close();
        delete concatWriter;
    }
    dbA.close();
    dbB.close();
}

unsigned int DBConcat::dbAKeyMap(unsigned int key) {
    if (sameDatabase)
        return key;

    std::pair<unsigned int, unsigned int> *originalMap = std::upper_bound(keysA, keysA + indexSizeA, key,
                                                                          compareKeyToFirstEntry());
    return (*originalMap).second;
}

unsigned int DBConcat::dbBKeyMap(unsigned int key) {
    if (sameDatabase)
        return key;

    std::pair<unsigned int, unsigned int> *originalMap = std::upper_bound(keysB, keysB + indexSizeB, key,
                                                                          compareKeyToFirstEntry());
    return (*originalMap).second;
}

DBConcat::~DBConcat() {
    if (!sameDatabase) {
        delete keysA;
        delete keysB;
    }
}

void setDbConcatDefault(Parameters *par) {
    par->threads = 1;
}

int concatdbs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setDbConcatDefault(&par);
    par.parseParameters(argc, argv, command, 3);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    struct timeval start, end;
    gettimeofday(&start, NULL);

    int datamode = DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX;
    DBConcat outDB(par.db1.c_str(), par.db1Index.c_str(),
                   par.db2.c_str(), par.db2Index.c_str(),
                   par.db3.c_str(), par.db3Index.c_str(),
                   static_cast<unsigned int>(par.threads), datamode, true, par.preserveKeysB);
    outDB.concat();

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for concatenating DBs: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    return EXIT_SUCCESS;
}
