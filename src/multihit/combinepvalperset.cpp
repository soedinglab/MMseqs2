#include "Debug.h"
#include "Parameters.h"
#include "Aggregation.h"
#include "itoa.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

double LBinCoeff(int M, int k) {
    return lgamma(M + 1) - lgamma(M - k + 1) - lgamma(k + 1);
}

class PvalueAggregator : public Aggregation {
public:
    PvalueAggregator(std::string queryDbName, std::string targetDbName, const std::string &resultDbName,
                     const std::string &outputDbName, float alpha, unsigned int threads) :
            Aggregation(targetDbName, resultDbName, outputDbName, threads), alpha(alpha) {

        std::string sizeDBName = queryDbName + "_set_size";
        std::string sizeDBIndex = queryDbName + "_set_size.index";
        querySizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str());
        querySizeReader->open(DBReader<unsigned int>::NOSORT);

        sizeDBName = targetDbName + "_set_size";
        sizeDBIndex = targetDbName + "_set_size.index";
        targetSizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str());
        targetSizeReader->open(DBReader<unsigned int>::NOSORT);
    }

    ~PvalueAggregator() {
        targetSizeReader->close();
        delete targetSizeReader;

        querySizeReader->close();
        delete querySizeReader;
    }

    //Get all result of a single Query Set VS a Single Target Set and return the multiple-match p-value for it
    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, unsigned int querySetKey,
                               unsigned int targetSetKey) {
        unsigned int targetGeneCount = (unsigned int) strtoull(targetSizeReader->getDataByDBKey(targetSetKey), NULL, 10);

        double pvalThreshold = alpha / targetGeneCount;
        const size_t numSets = targetSizeReader->getSize();

        std::string buffer;
        char keyBuffer[255];
        char *tmpBuff = Itoa::u32toa_sse2(targetSetKey, keyBuffer);
        buffer.append(keyBuffer, tmpBuff - keyBuffer - 1);
        buffer.append("\t");

        //edge case p0 = 0
        if (pvalThreshold == 0.0) {
            buffer.append(SSTR(numSets));
            return buffer;
        }

        size_t k = 0;
        double r = 0;
        const double logPvalThr = log(pvalThreshold);
        for (size_t i = 0; i < dataToAggregate.size(); ++i) {
            double pvalue = std::strtod(dataToAggregate[i][1].c_str(), NULL);
            if (pvalue < pvalThreshold) {
                k++;
                r -= log(pvalue) - logPvalThr;
            }
        }

        if (r == 0) {
            buffer.append(SSTR(numSets));
            return buffer;
        }

        if (std::isinf(r)) {
            buffer.append("0");
            return buffer;
        }        

        const double expMinusR = exp(-r);
        
        //edge case p0 = 1
        if (pvalThreshold == 1.0) {
            buffer.append(SSTR(expMinusR * numSets));
            return buffer;
        }


        unsigned int orfCount = (unsigned int) strtoull(querySizeReader->getDataByDBKey(querySetKey), NULL, 10); 
        double truncatedFisherPval = 0;
        double logRatioRToFactorial = 0;
        const double logR = log(r);
        const double log1MinusPvalThr = log(1.0 - pvalThreshold);
        for (size_t i = 0; i < orfCount; ++i) { 
            double logProbKGoodHits = 1;
            double firstFactor = (LBinCoeff(orfCount, (i + 1)) + (i + 1) * logPvalThr + (orfCount - (i +1)) * log1MinusPvalThr); 
            for (size_t j = i + 2; j < orfCount + 1; ++j) { 
                logProbKGoodHits += exp(LBinCoeff(orfCount, j) + j * logPvalThr + (orfCount - j) * log1MinusPvalThr - firstFactor);
            }
            logProbKGoodHits = firstFactor + log(logProbKGoodHits);
            if (i > 0) {
                logRatioRToFactorial = logRatioRToFactorial + logR - log(i);
            }
            truncatedFisherPval += exp(logRatioRToFactorial + logProbKGoodHits);
        }


        
        double updatedPval = expMinusR * truncatedFisherPval;
        double updatedEval = updatedPval * numSets;
        buffer.append(SSTR(updatedEval));
        return buffer;
    }

private:
    double alpha;
    DBReader<unsigned int> *querySizeReader;
    DBReader<unsigned int> *targetSizeReader;
};

int combinepvalperset(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4, true);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    PvalueAggregator aggregation(par.db1, par.db2, par.db3, par.db4, par.alpha, (unsigned int) par.threads);
    return aggregation.run();
}
