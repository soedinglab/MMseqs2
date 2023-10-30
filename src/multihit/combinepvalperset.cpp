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
                     const std::string &outputDbName, float alpha, unsigned int threads, unsigned int compressed, int aggregationMode) :
            Aggregation(targetDbName, resultDbName, outputDbName, threads, compressed), alpha(alpha), aggregationMode(aggregationMode) {

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
        
        const size_t numTargetSets = targetSizeReader->getSize();  
        double updatedPval;

        std::string buffer;
        char keyBuffer[255];
        char *tmpBuff = Itoa::u32toa_sse2(targetSetKey, keyBuffer);
        buffer.append(keyBuffer, tmpBuff - keyBuffer - 1);
        buffer.append("\t");

        //0) multihit P-values
        if(aggregationMode == Parameters::AGGREGATION_MODE_MULTIHIT){
            unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx)); 
            double pvalThreshold = alpha / (orfCount + 1);

            //multihit edge case p0 = 0
            if (pvalThreshold == 0.0) {
                buffer.append(SSTR(numTargetSets));
                return buffer;
            }

            // size_t k = 0;
            double r = 0;
            const double logPvalThr = log(pvalThreshold);
            for (size_t i = 0; i < dataToAggregate.size(); ++i) {
                double logPvalue = std::strtod(dataToAggregate[i][1].c_str(), NULL);
                if (logPvalue < logPvalThr) {
                    // k++;
                    r -= logPvalue - logPvalThr;
                }
            }
            //multihit edge case r = 0
            if (r == 0) {
                buffer.append(SSTR(numTargetSets));
                return buffer;
            }

            if (std::isinf(r)) {
                buffer.append("0");
                return buffer;
            }        

            const double expMinusR = exp(-r);
            
            //multihit edge case p0 = 1
            if (pvalThreshold == 1.0) {
                buffer.append(SSTR(expMinusR * numTargetSets));
                return buffer;
            }


            double truncatedFisherPval = 0;
            const double logR = log(r);
            for (size_t i = 0; i < orfCount; ++i) { 
                truncatedFisherPval += exp(i*logR - lGammaLookup[i+1] + logBiLookup[thread_idx][i]);
            }            
            updatedPval = expMinusR * truncatedFisherPval;       
        } 

        //1) the minimum of all P-values(as a baseline)
        else if(aggregationMode == Parameters::AGGREGATION_MODE_MIN_PVAL){
            unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
            double minLogPval = 0;
            for (size_t i = 0; i < dataToAggregate.size(); ++i) { 
                double currentLogPval = std::strtod(dataToAggregate[i][1].c_str(), NULL);
                if (currentLogPval < minLogPval) {
                    minLogPval = currentLogPval;
                };
            }
            updatedPval = 1 -  exp( - exp(minLogPval) * orfCount);;    
        }      

        //2) the P-value for the product-of-P-values
        else if (aggregationMode == Parameters::AGGREGATION_MODE_PRODUCT)    {
            double  sumLogPval= 0;
            for (size_t i = 0; i < dataToAggregate.size(); ++i) {
                double logPvalue = std::strtod(dataToAggregate[i][1].c_str(), NULL);
                sumLogPval += logPvalue;
            }
            updatedPval = exp(sumLogPval);   
        }

        //3) the P-values of the (modified) truncated product method 
        else if(aggregationMode == Parameters::AGGREGATION_MODE_TRUNCATED_PRODUCT){
            //new theory: taking the best hit regardless of threshold and (from second hit on)sum of how much it surpassed threshold
            unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
            double logPvalThreshold = log(alpha / (orfCount + 1));
            double minLogPval = 0;
            double sumLogPval = 0; 
            size_t k = 0;
            for (size_t i = 0; i < dataToAggregate.size(); ++i) {
                double logPvalue = std::strtod(dataToAggregate[i][1].c_str(), NULL);
                if (logPvalue < minLogPval) {
                    if (logPvalue == 0) {
                        //to avoid -0.0
                        minLogPval = logPvalue;
                    }
                    else {minLogPval = -logPvalue;}
                }
                if (logPvalue < logPvalThreshold) {
                    //sum up the part exceeding logThreshold, add a minus to make score positive
                    sumLogPval -= logPvalue - logPvalThreshold;
                    k++;
                }
            }
            if(k == 0){
                //if no hit passed thr, take the -log of best hit pval as score
                buffer.append(SSTR(minLogPval));
                return buffer;
            }
            else {
                //if one or more hits passed thr
                buffer.append(SSTR(sumLogPval - logPvalThreshold));
                return buffer;
            }
        }

        
        else {
            Debug(Debug::ERROR) << "Invalid aggregation function!\n";
            EXIT(EXIT_FAILURE);
        }
        double updatedEval = updatedPval * numTargetSets;
        buffer.append(SSTR(updatedEval));
        return buffer;
    }

private:
    double alpha;
    int aggregationMode;
    DBReader<unsigned int> *querySizeReader;
    DBReader<unsigned int> *targetSizeReader;
    double* lGammaLookup;
    double** logBiLookup;
};

int combinepvalperset(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    PvalueAggregator aggregation(par.db1, par.db2, par.db3, par.db4, par.alpha, (unsigned int) par.threads, par.compressed, par.aggregationMode);
    return aggregation.run();
}
