#include "Debug.h"
#include "Parameters.h"
#include "Aggregation.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

class BestHitBySetFilter : public Aggregation {
public :
    BestHitBySetFilter(const std::string &targetDbName, const std::string &resultDbName,
                       const std::string &outputDbName, bool simpleBestHitMode, unsigned int threads) :
            Aggregation(targetDbName, resultDbName, outputDbName, threads), simpleBestHitMode(simpleBestHitMode) {
        std::string sizeDbName = targetDbName + "_set_size";
        std::string sizeDbIndex = targetDbName + "_set_size.index";
        targetSizeReader = new DBReader<unsigned int>(sizeDbName.c_str(), sizeDbIndex.c_str());
        targetSizeReader->open(DBReader<unsigned int>::NOSORT);
    }

    ~BestHitBySetFilter() {
        targetSizeReader->close();
        delete targetSizeReader;
    }


    void prepareInput(unsigned int, unsigned int) {}

    std::string aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, unsigned int, unsigned int targetSetKey, unsigned int)  {
        std::string buffer;
        buffer.reserve(1024);

        double bestScore = 0;
        double secondBestScore = 0;
        double bestEval = DBL_MAX;

        double correctedPval = 0;

        // Look for the lowest p-value and retain only this line
        // dataToAggregate = [nbrTargetGene][Field of result]
        std::vector<std::string> *bestEntry;
        for (size_t i = 0; i < dataToAggregate.size(); i++) {
            double score = strtod(dataToAggregate[i][1].c_str(), NULL);
            double eval = strtod(dataToAggregate[i][3].c_str(), NULL);

            if (score > bestScore) {
                secondBestScore = bestScore;
                bestScore = score;
                if (simpleBestHitMode == false) {
                    bestEntry = &dataToAggregate[i];
                }
            }

            if (simpleBestHitMode == true && bestEval > eval) {
                bestEval = eval;
                bestEntry = &dataToAggregate[i];
            }
        }

        size_t targetId = targetSizeReader->getId(targetSetKey);
        if (targetId == UINT_MAX) {
            Debug(Debug::ERROR) << "Invalid target size database key " << targetSetKey << ".\n";
            EXIT(EXIT_FAILURE);
        }
        char *data = targetSizeReader->getData(targetId);
        unsigned int nbrGenes = Util::fast_atoi<unsigned int>(data);

        if (simpleBestHitMode) {
            correctedPval = bestEval / nbrGenes;
        } else {
            // if no second hit is available, update pvalue with fake hit
            if (dataToAggregate.size() < 2) {
                secondBestScore = 2.0 * log((nbrGenes + 1) / 2) / log(2.0);
            }
            correctedPval = pow(2.0, secondBestScore / 2 - bestScore / 2);
        }

        // Aggregate the full line into string
        for (size_t i = 0; i < bestEntry->size(); ++i) {
            if (i == 1) {
                char tmpBuf[15];
                sprintf(tmpBuf, "%.3E", correctedPval);
                buffer.append(tmpBuf);
            } else {
                buffer.append(bestEntry->at(i));
            }
            buffer.append("\t");
        }

        return buffer;
    }

private:
    DBReader<unsigned int> *targetSizeReader;
    bool simpleBestHitMode;
};


int besthitperset(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3, true);

    BestHitBySetFilter aggregation(par.db2, par.db3, par.db4, par.simpleBestHit, (unsigned int) par.threads);
    return aggregation.run();
}
