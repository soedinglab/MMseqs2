#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "MathUtil.h"
#include "Matcher.h"

#ifdef OPENMP
#include <omp.h>
#endif


int summarizealis(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);


    DBReader<unsigned int> resultReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer.open();
    Debug::Progress progress(resultReader.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);

        std::string annotation;
        annotation.reserve(1024*1024);

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int dbKey = resultReader.getDbKey(i);
            char *tabData = resultReader.getData(i,  thread_idx);
            Matcher::readAlignmentResults(alnResults, tabData);
            std::stable_sort(alnResults.begin(), alnResults.end(), Matcher::compareHitByPos);

            if (alnResults.size() == 0) {
                continue;
            }
            float resCov = 0;
            float avgSeqId = 0.0f;
            unsigned int seqLen = 1;
            float uniqCov = 0;
            std::vector<bool> covered(alnResults[0].qLen, false);
            int prevQEndPos = -1;

            for (size_t i = 0; i < alnResults.size(); i++) {
                Matcher::result_t res = alnResults[i];
                seqLen = res.qLen;
                int qStartPos  = std::min( res.qStartPos, res.qEndPos);
                int qEndPos = std::max( res.qStartPos, res.qEndPos);
                uniqCov += std::max(prevQEndPos, res.qEndPos) - std::max(prevQEndPos, res.qStartPos);
                resCov += static_cast<float>(qEndPos - qStartPos);
                avgSeqId += res.seqId;
                prevQEndPos = std::max(prevQEndPos, res.qEndPos);
            }
            avgSeqId = avgSeqId / static_cast<float>(alnResults.size());
            resCov = resCov / static_cast<float>(seqLen);
            uniqCov = uniqCov / static_cast<float>(seqLen);
            annotation.append(SSTR(alnResults.size()));
            annotation.append("\t");
            annotation.append(SSTR(uniqCov));
            annotation.append("\t");
            annotation.append(SSTR(resCov));
            annotation.append("\t");
            annotation.append(SSTR(avgSeqId));
            annotation.append("\n");
            writer.writeData(annotation.c_str(), annotation.length(), dbKey, thread_idx);
            alnResults.clear();
            annotation.clear();
        }
    }
    writer.close();

    return EXIT_SUCCESS;
}


