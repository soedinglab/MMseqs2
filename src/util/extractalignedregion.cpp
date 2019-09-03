#include "DBWriter.h"
#include "Debug.h"
#include "DBReader.h"
#include "Util.h"
#include "Matcher.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

int extractalignedregion(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // never allow deletions
    par.allowDeletion = false;

    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qdbr.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        qdbr.readMmapedDataInMemory();
    }

    bool sameDB = false;
    DBReader<unsigned int> *tdbr = NULL;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            tdbr->readMmapedDataInMemory();
        }
    }

    DBReader<unsigned int> alndbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alndbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, tdbr->getDbtype());
    dbw.open();
    Debug::Progress progress(alndbr.getSize());

    const char newline = '\n';
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> results;
        results.reserve(300);
#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < alndbr.getSize(); i++) {
            progress.updateProgress();

            unsigned int queryKey = alndbr.getDbKey(i);
            char *qSeq = NULL;
            if (par.extractMode == Parameters::EXTRACT_QUERY) {
                qSeq = qdbr.getDataByDBKey(queryKey, thread_idx);
            }

            char *data = alndbr.getData(i, thread_idx);
            Matcher::readAlignmentResults(results, data);
            for (size_t j = 0; j < results.size(); j++) {
                Matcher::result_t& res = results[j];
                size_t length = 0;
                char* seq = NULL;
                if (qSeq) {
                    seq = qSeq + res.qStartPos;
                    length = res.qEndPos - res.qStartPos+1;
                } else if (par.extractMode == Parameters::EXTRACT_TARGET) {
                    seq = tdbr->getDataByDBKey(res.dbKey, thread_idx) + res.dbStartPos;
                    length = res.dbEndPos - res.dbStartPos+1;
                } else {
                    Debug(Debug::ERROR) << "Missing extraction type!\n";
                    EXIT(EXIT_FAILURE);
                }

                dbw.writeStart(thread_idx);
                dbw.writeAdd(seq, length, thread_idx);
                dbw.writeAdd(&newline, 1, thread_idx);
                dbw.writeEnd(queryKey, thread_idx);
            }
            results.clear();
        }
    }
    dbw.close();

    if (par.extractMode == Parameters::EXTRACT_QUERY) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.db4, DBFiles::SEQUENCE_ANCILLARY);
    } else {
        DBReader<unsigned int>::softlinkDb(par.db2, par.db4, DBFiles::SEQUENCE_ANCILLARY);
    }

    alndbr.close();
    qdbr.close();
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }

    return EXIT_SUCCESS;
}
