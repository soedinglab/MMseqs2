#include "Matcher.h"
#include "DBReader.h"
#include "Debug.h"
#include "DBWriter.h"
#include "Util.h"
#include "Parameters.h"

#include <map>

#ifdef OPENMP
#include <omp.h>
#endif

int subtractdbs(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, 0, 0);
    par.evalProfile = (par.evalThr < par.evalProfile) ? par.evalThr : par.evalProfile;
    par.printParameters(command.cmd, argc, argv, *command.params);
    const double evalThreshold = par.evalProfile;

    DBReader<unsigned int> leftDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    leftDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> rightDbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    rightDbr.open(DBReader<unsigned int>::NOSORT);

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, leftDbr.getSize()), (size_t)1);
#endif

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), localThreads, par.compressed, leftDbr.getDbtype());
    writer.open();

    Debug::Progress progress(leftDbr.getSize());
#pragma omp parallel num_threads(localThreads)
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        const char *entry[255];
        char key[255];
        std::string result;
        result.reserve(100000);

#pragma omp  for schedule(dynamic, 10)
        for (size_t id = 0; id < leftDbr.getSize(); id++) {
            progress.updateProgress();
            std::map<unsigned int, bool> elementLookup;
            const char *leftData = leftDbr.getData(id, thread_idx);
            unsigned int leftDbKey = leftDbr.getDbKey(id);

            // fill element id look up with left side elementLookup
            {
                char *data = (char *) leftData;
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    unsigned int dbKey = std::strtoul(key, NULL, 10);
                    double evalue = 0.0;
                    const size_t columns = Util::getWordsOfLine(data, entry, 255);
                    // its an aln result (parse e-value)
                    if (columns >= Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                        evalue = strtod(entry[3], NULL);
                    }
                    if (evalue <= evalThreshold) {
                        elementLookup[dbKey] = true;
                    }
                    data = Util::skipLine(data);
                }
            }
            // get all data for the leftDbkey from rightDbr
            // check if right ids are in elementsId
            char *data = rightDbr.getDataByDBKey(leftDbKey, thread_idx);

            if (data != NULL) {
                while (*data != '\0') {
                    Util::parseKey(data, key);
                    unsigned int element = std::strtoul(key, NULL, 10);
                    double evalue = 0.0;
                    const size_t columns = Util::getWordsOfLine(data, entry, 255);
                    if (columns >= Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                        evalue = strtod(entry[3], NULL);
                    }
                    if (evalue <= evalThreshold) {
                        elementLookup[element] = false;
                    }
                    data = Util::skipLine(data);
                }
            }
            // write only elementLookup that are not found in rightDbr (id != UINT_MAX)
            {
                char *data = (char *) leftData;
                while (*data != '\0') {
                    char *start = data;
                    data = Util::skipLine(data);
                    Util::parseKey(start, key);
                    unsigned int elementIdx = std::strtoul(key, NULL, 10);
                    if (elementLookup[elementIdx]) {
                        result.append(start, data - start);
                    }
                }
            }

            writer.writeData(result.c_str(), result.length(), leftDbKey, thread_idx);
            result.clear();
        }
    }
    writer.close();

    leftDbr.close();
    rightDbr.close();

    return EXIT_SUCCESS;
}
