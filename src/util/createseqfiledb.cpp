#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"

#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int createseqfiledb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> headerDb(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    headerDb.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        headerDb.readMmapedDataInMemory();
    }

    DBReader<unsigned int> seqDb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    seqDb.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        seqDb.readMmapedDataInMemory();
    }

    DBReader<unsigned int> resultDb(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    resultDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer.open();

    Debug::Progress progress(resultDb.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string result;
        result.reserve(1024);

        char dbKey[255];
#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < resultDb.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = resultDb.getDbKey(i);
            char *data = resultDb.getData(i, thread_idx);

            size_t entries = Util::countLines(data, resultDb.getEntryLen(i) - 1);
            if (entries < (unsigned int) par.minSequences || entries > (unsigned int) par.maxSequences) {
                continue;
            }

            size_t entries_num = 0;
            while (*data != '\0') {
                entries_num++;
                Util::parseKey(data, dbKey);
                data = Util::skipLine(data);

                const unsigned int memberKey = (unsigned int) strtoul(dbKey, NULL, 10);
                size_t headerId = headerDb.getId(memberKey);
                if (headerId == UINT_MAX) {
                    Debug(Debug::ERROR) << "Entry " << key << " does not contain a sequence!" << "\n";
                    EXIT(EXIT_FAILURE);
                }
                size_t seqId = seqDb.getId(memberKey);
                if (seqId == UINT_MAX) {
                    Debug(Debug::ERROR) << "Entry " << key << " does not contain a sequence!" << "\n";
                    EXIT(EXIT_FAILURE);
                }
                if (entries_num == 1 && par.hhFormat) {
                    char *header = headerDb.getData(headerId, thread_idx);
                    size_t headerLen = headerDb.getEntryLen(headerId) - 1;
                    size_t accessionLen = Util::skipNoneWhitespace(header);
                    char *sequence = seqDb.getData(headerId, thread_idx);
                    size_t sequenceLen = seqDb.getEntryLen(headerId) - 1;
                    result.append(1, '#');
                    result.append(header, headerLen);
                    result.append(1, '>');
                    result.append(header, accessionLen);
                    result.append("_consensus\n");
                    result.append(sequence, seqDb.getEntryLen(headerId) - 1);
                    result.append(1, '>');
                    result.append(header, headerLen);
                    result.append(sequence, sequenceLen);
                } else {
                    result.append(1, '>');
                    result.append(headerDb.getData(headerId, thread_idx), headerDb.getEntryLen(headerId) - 1);
                    result.append(seqDb.getData(headerId, thread_idx), seqDb.getEntryLen(headerId) - 1);
                }
            }
            writer.writeData(result.c_str(), result.length(), key, thread_idx);
            result.clear();
        }
    }

    writer.close();
    resultDb.close();
    seqDb.close();
    headerDb.close();

    return EXIT_SUCCESS;
}
