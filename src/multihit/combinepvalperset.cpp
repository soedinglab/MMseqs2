#include "Debug.h"
#include "Parameters.h"
#include "Aggregation.h"
#include "itoa.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif


double LBinCoeff(double* lookup, int M, int k) {
    return lookup[M + 1] - lookup[M - k + 1] - lookup[k + 1];
}

// Precompute coefficients logB[i] = log(B[i])
void precomputeLogB(const unsigned int orfCount, const double pvalThreshold, double* lGammaLookup, double *logB) { 
    double logPvalThr = log(pvalThreshold);
    double log1MinusPvalThr = log(1 - pvalThreshold);
    logB[orfCount - 1] = orfCount * logPvalThr;
    for (int i = (orfCount - 2); i >= 0; i--){
        int k = i + 1;
        double log_newTerm = LBinCoeff(lGammaLookup, orfCount, k) + k * logPvalThr + (orfCount - k) * log1MinusPvalThr;
        logB[i] = logB[i + 1] + log(1 + exp(log_newTerm - logB[i+1]));
    }
}


class PvalueAggregator : public Aggregation {
public:
    PvalueAggregator(std::string queryDbName, std::string targetDbName, const std::string &resultDbName,
                     const std::string &outputDbName, float alpha, unsigned int threads, unsigned int compressed) :
            Aggregation(targetDbName, resultDbName, outputDbName, threads, compressed), alpha(alpha) {

        std::string sizeDBName = queryDbName + "_set_size";
        std::string sizeDBIndex = queryDbName + "_set_size.index";
        querySizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        querySizeReader->open(DBReader<unsigned int>::NOSORT);

        sizeDBName = targetDbName + "_set_size";
        sizeDBIndex = targetDbName + "_set_size.index";
        targetSizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        targetSizeReader->open(DBReader<unsigned int>::NOSORT);

        unsigned int maxOrfCount = 0;
        for (size_t i = 0; i < querySizeReader->getSize(); ++i) { 
            unsigned int currentCount = Util::fast_atoi<unsigned int>(querySizeReader->getData(i, 0));
            if (currentCount > maxOrfCount) {
                maxOrfCount = currentCount;
            };
        }

        lGammaLookup = new double[maxOrfCount + 2];
        for (size_t i = 0; i < maxOrfCount + 2; ++i) { 
            lGammaLookup[i] = lgamma(i);
        }

        logBiLookup = new double*[threads];
        for (size_t i = 0; i < threads; ++i) {
            logBiLookup[i] = new double[maxOrfCount];
        }
    }

    ~PvalueAggregator() {
        for (size_t i = 0; i < threads; ++i) {
            delete[] logBiLookup[i];
        }
        delete[] logBiLookup;

        delete[] lGammaLookup;

        targetSizeReader->close();
        delete targetSizeReader;

        querySizeReader->close();
        delete querySizeReader;
    }

    void prepareInput(unsigned int querySetKey, unsigned int thread_idx) {
        unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
        precomputeLogB(orfCount, alpha/(orfCount + 1), lGammaLookup, logBiLookup[thread_idx]);
    }

    //Get all result of a single Query Set VS a Single Target Set and return the multiple-match p-value for it
    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, unsigned int querySetKey,
                               unsigned int targetSetKey, unsigned int thread_idx) {
        unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
        double pvalThreshold = alpha / (orfCount + 1);
        const size_t numTargetSets = targetSizeReader->getSize();

        std::string buffer;
        char keyBuffer[255];
        char *tmpBuff = Itoa::u32toa_sse2(targetSetKey, keyBuffer);
        buffer.append(keyBuffer, tmpBuff - keyBuffer - 1);
        buffer.append("\t");

        //edge case p0 = 0
        if (pvalThreshold == 0.0) {
            buffer.append(SSTR(numTargetSets));
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
            buffer.append(SSTR(numTargetSets));
            return buffer;
        }

        if (std::isinf(r)) {
            buffer.append("0");
            return buffer;
        }        

        const double expMinusR = exp(-r);
        
        //edge case p0 = 1
        if (pvalThreshold == 1.0) {
            buffer.append(SSTR(expMinusR * numTargetSets));
            return buffer;
        }


        double truncatedFisherPval = 0;
        const double logR = log(r);
        for (size_t i = 0; i < orfCount; ++i) { 
            truncatedFisherPval += exp(i*logR - lGammaLookup[i+1] + logBiLookup[thread_idx][i]);
        }
        
        double updatedPval = expMinusR * truncatedFisherPval;
        double updatedEval = updatedPval * numTargetSets;
        buffer.append(SSTR(updatedEval));
        return buffer;
    }

private:
    double alpha;
    DBReader<unsigned int> *querySizeReader;
    DBReader<unsigned int> *targetSizeReader;
    double* lGammaLookup;
    double** logBiLookup;
};

int combinepvalperset(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4, true);

    PvalueAggregator aggregation(par.db1, par.db2, par.db3, par.db4, par.alpha, (unsigned int) par.threads, par.compressed);
    return aggregation.run();
}
