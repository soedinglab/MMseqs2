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
                       const std::string &outputDbName, bool simpleBestHitMode, unsigned int threads, unsigned int compressed) :
            Aggregation(targetDbName, resultDbName, outputDbName, threads, compressed), simpleBestHitMode(simpleBestHitMode) {
        std::string sizeDbName = targetDbName + "_set_size";
        std::string sizeDbIndex = targetDbName + "_set_size.index";
        targetSizeReader = new DBReader<unsigned int>(sizeDbName.c_str(), sizeDbIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        targetSizeReader->open(DBReader<unsigned int>::NOSORT);
    }

    ~BestHitBySetFilter() {
        targetSizeReader->close();
        delete targetSizeReader;
    }


    void prepareInput(unsigned int, unsigned int) {}

    std::string aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, unsigned int, unsigned int targetSetKey, unsigned int thread_idx)  {
        double bestScore = -DBL_MAX;
        double secondBestScore = -DBL_MAX;
        double bestEval = DBL_MAX;
        double logBestHitCalibration = log(1);

        double logCorrectedPval = 0;

        // Look for the lowest p-value and retain only this line
        // dataToAggregate = [nbrTargetGene][Field of result]
        size_t targetId = targetSizeReader->getId(targetSetKey);
        if (targetId == UINT_MAX) {
            Debug(Debug::ERROR) << "Invalid target size database key " << targetSetKey << ".\n";
            EXIT(EXIT_FAILURE);
        }
        char *data = targetSizeReader->getData(targetId, thread_idx);
        unsigned int nbrGenes = Util::fast_atoi<unsigned int>(data);

        std::vector<std::string> *bestEntry = NULL;
        for (size_t i = 0; i < dataToAggregate.size(); i++) {
            double eval = strtod(dataToAggregate[i][3].c_str(), NULL);
            double pval = eval/nbrGenes;
            //prevent log(0)
            if (pval == 0) {
                pval = DBL_MIN;
            }
            double score = -log(pval);

            //if only one hit use simple best hit
            if(simpleBestHitMode ||dataToAggregate.size() < 2) {
                if(bestEval > eval){
                    bestEval = eval;
                    bestEntry = &dataToAggregate[i];
                }
            }
            else {
                if (score >= bestScore) {
                    secondBestScore = bestScore;
                    bestScore = score;
                    bestEntry = &dataToAggregate[i];
                } 
                else if (score > secondBestScore) {
                    secondBestScore = score;
                }
            }
        }


        if (simpleBestHitMode ||dataToAggregate.size() < 2) {
            if(bestEval == 0) {
                logCorrectedPval = log(DBL_MIN)-logBestHitCalibration;
            }
            else if (bestEval > 0 && bestEval < 10e-4){
                logCorrectedPval = log(bestEval)-logBestHitCalibration;
            }
            else {
                logCorrectedPval = log(1 - exp(-bestEval))-logBestHitCalibration;
            }
            
        } 
        else {
            logCorrectedPval = secondBestScore - bestScore;
        }

        if (bestEntry == NULL) {
            return "";
        }

        std::string buffer;
        buffer.reserve(1024);

        // Aggregate the full line into string
        for (size_t i = 0; i < bestEntry->size(); ++i) {
            if (i == 1) {
                char tmpBuf[15];
                snprintf(tmpBuf, sizeof(tmpBuf), "%.3E", logCorrectedPval);
                buffer.append(tmpBuf);
            } else {
                buffer.append(bestEntry->at(i));
            }
            if (i != (bestEntry->size() - 1)) {
                buffer.append(1, '\t');
            }
        }

        return buffer;
    }

private:
    DBReader<unsigned int> *targetSizeReader;
    bool simpleBestHitMode;
};


int besthitperset(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    BestHitBySetFilter aggregation(par.db2, par.db3, par.db4, par.simpleBestHit, (unsigned int) par.threads, par.compressed);
    return aggregation.run();
}
