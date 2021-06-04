#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include "QueryMatcher.h"
#include "FastSort.h"

#ifdef OPENMP
#include <omp.h>
#endif

int sortresult(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();
    Debug::Progress progress(reader.getSize());

   #pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        const char *entry[255];
        char buffer[1024 + 32768*4];

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);

        std::vector<hit_t> prefResults;
        prefResults.reserve(300);

#pragma omp for schedule(dynamic, 5)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);

            int format = -1;
            while (*data != '\0') {
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns >= Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    alnResults.emplace_back(Matcher::parseAlignmentRecord(data, true));
                    format = columns >= Matcher::ALN_RES_WITH_BT_COL_CNT ? 1 : 0;
                } else if (columns == 3) {
                    prefResults.emplace_back(QueryMatcher::parsePrefilterHit(data));
                    format = 2;
                } else {
                    Debug(Debug::ERROR) << "Invalid input result format ("<<columns<<" columns).\n";
                    EXIT(EXIT_FAILURE);
                }
                data = Util::skipLine(data);
            }

            writer.writeStart(thread_idx);
            if (format == 0 || format == 1) {
                SORT_SERIAL(alnResults.begin(), alnResults.end(), Matcher::compareHits);
                for (size_t i = 0; i < alnResults.size(); ++i) {
                    size_t length = Matcher::resultToBuffer(buffer, alnResults[i], format == 1, false);
                    writer.writeAdd(buffer, length, thread_idx);
                }
            } else if (format == 2) {
                SORT_SERIAL(prefResults.begin(), prefResults.end(), hit_t::compareHitsByScoreAndId);
                for (size_t i = 0; i < prefResults.size(); ++i) {
                    size_t length = QueryMatcher::prefilterHitToBuffer(buffer, prefResults[i]);
                    writer.writeAdd(buffer, length, thread_idx);
                }
            }
            writer.writeEnd(key, thread_idx);

            alnResults.clear();
            prefResults.clear();
        }
    }

    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}


