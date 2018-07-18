#include "Aggregation.h"
#include "Util.h"
#include "Parameters.h"
#include "Debug.h"

#include <cfloat>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <cmath>

#ifdef OPENMP

#include <omp.h>

#endif

#define EVAL_COLUMN 3

struct compareByStart {
    bool operator()(const std::pair<long, long> &lhs,
                    const std::pair<long, long> &rhs) const {
        return (lhs.first < rhs.first);
    }
};

Aggregation::Aggregation(const std::string &resultDbName, const std::string &outputDbName, size_t setColumn, unsigned int threads)
        : resultDbName(resultDbName), outputDbName(outputDbName), setColumn(setColumn), threads(threads) {}

// build a map with the value in [target column] field as a key and the rest of the line, cut in fields, as values
bool Aggregation::buildMap(std::stringstream &dataStream, std::map<unsigned int, std::vector<std::vector<std::string> > > &dataToAggregate) {
    dataToAggregate.clear();
    std::string line;
    while (std::getline(dataStream, line)) {
        if (!line.empty()) {
            std::vector<std::string> columns = Util::split(line, "\t");
            unsigned int dbKey = (unsigned int) strtoull(columns[setColumn].c_str(), NULL, 10);
            dataToAggregate[dbKey].push_back(columns);
        }
    }
    return true;
}

// General workflow for search Result processing
int Aggregation::run() {
    std::string inputDBIndex = resultDbName + ".index";
    DBReader<unsigned int> reader(resultDbName.c_str(), inputDBIndex.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    std::string outputDBIndex = outputDbName + ".index";
    DBWriter writer(outputDbName.c_str(), outputDBIndex.c_str(), threads);
    writer.open();

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        std::map<unsigned int, std::vector<std::vector<std::string>>> dataToMerge;

#pragma omp for
        for (size_t i = 0; i < reader.getSize(); i++) {
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i);
            std::stringstream dataStream(data);
            if (buildMap(dataStream, dataToMerge)) {
                std::string outputBuffer;
                outputBuffer.reserve(1048576);
                for (std::map<unsigned int, std::vector<std::vector<std::string>>>::const_iterator it = dataToMerge.begin();
                     it != dataToMerge.end(); ++it) {
                    unsigned int targetKey = it->first;
                    std::vector<std::vector<std::string>> columns = it->second;
                    outputBuffer.append(aggregateEntry(columns, key, targetKey));
                    outputBuffer.append("\n");
                }

                writer.writeData(outputBuffer.c_str(), outputBuffer.length(), key, thread_idx);
            } else {
                Debug(Debug::ERROR) << "buildMap failed \n";
                EXIT(EXIT_FAILURE);
            }
        }
    };
    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}


BestHitAggregator::BestHitAggregator(const std::string &targetDbName, const std::string &resultDbName,
                                     const std::string &outputDbName, bool simpleBestHitMode, unsigned int threads) :
        Aggregation(resultDbName, outputDbName, 10, threads), simpleBestHitMode(simpleBestHitMode) {
    std::string sizeDbName = targetDbName + "_orfs_size";
    std::string sizeDbIndex = targetDbName + "_orfs_size.index";
    targetSizeReader = new DBReader<unsigned int>(sizeDbName.c_str(), sizeDbIndex.c_str());
    targetSizeReader->open(DBReader<unsigned int>::NOSORT);
}

// Only Keep the best hits of a protein against each Target Set
std::string BestHitAggregator::aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, unsigned int querySetKey,
                                  unsigned int targetSetKey) {
    std::string buffer;
    buffer.reserve(1024);

    double bestScore = 0;
    size_t maxId = 0;
    size_t bestEvalId = 0;
    double secondBestScore = 0;
    double correctedPval = 0;

    double eval = 0;
    double bestEval = DBL_MAX;
    // Look for the lowest p-value and retain only this line
    // dataToAggregate = [nbrTargetGene][Field of result]
    for (size_t i = 0; i < dataToAggregate.size(); i++) {
        double score = strtod(dataToAggregate[i][1].c_str(), NULL);
        eval = strtod(dataToAggregate[i][EVAL_COLUMN].c_str(), NULL);

        if (score > bestScore) {
            secondBestScore = bestScore;
            bestScore = score;
            maxId = i;
        }

        if (bestEval > eval) {
            bestEval = eval;
            bestEvalId = i;
        }

    }

    unsigned int targetKey = static_cast<unsigned int>(std::strtoull(dataToAggregate[maxId][10].c_str(), NULL, 10));
    char *data = targetSizeReader->getDataByDBKey(targetKey);
    unsigned int nbrGenes = (unsigned int) std::strtoull(data, NULL, 10);

    if (simpleBestHitMode) {
        maxId = bestEvalId;
        correctedPval = bestEval / nbrGenes;
    } else {
        // ortholog mode
        // exp(log(secondBestEval / nbrGenes) - log(secondBestEval / nbrGenes));
        // TODO what happens with == 1 ????
        if (dataToAggregate.size() < 2) {
            // (2.00/(nbrGenes+1.00))
            // (2.00*nbrGenes/(nbrGenes+1.00));
            secondBestScore = 2.0 * log((nbrGenes + 1) / 2) / log(2.0);
        }
        correctedPval = pow(2.0, secondBestScore / 2 - bestScore / 2);
    }


    char tmpBuf[15];
    sprintf(tmpBuf, "%.3E", correctedPval);
    dataToAggregate[maxId][1] = std::string(tmpBuf);

    // Aggregate the full line in one string
    for (std::vector<std::string>::const_iterator it = dataToAggregate[maxId].begin();
         it != dataToAggregate[maxId].end(); ++it) {
        buffer.append(*it);
        buffer.append("\t");
    }
    return buffer;
}


PvalueAggregator::PvalueAggregator(std::string queryDbName, std::string targetDbName, const std::string &argInputDBname,
                                   const std::string &argOutputDBname, float alpha, unsigned int threads) :
        Aggregation(argInputDBname, argOutputDBname, 10, threads), alpha(alpha) {

    std::string sizeDBName  = queryDbName + "_orfs_size";
    std::string sizeDBIndex = queryDbName + "_orfs_size.index";
    querySizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str());
    querySizeReader->open(DBReader<unsigned int>::NOSORT);

    sizeDBName = targetDbName + "_orfs_size";
    sizeDBIndex = targetDbName + "_orfs_size.index";
    targetSizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str());
    targetSizeReader->open(DBReader<unsigned int>::NOSORT);
}

double LBinCoeff(int M, int k) {
    return lgamma(M + 1) - lgamma(M - k + 1) - lgamma(k + 1);
}

double factorial(size_t i) {
    double fac = 1;
    for (size_t j = 1; j < i + 1; j++) {
        fac *= j;
    }
    return fac;
}

//Get all result of a single Query Set VS a Single Target Set and return the multiple-match p-value for it
std::string PvalueAggregator::aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, unsigned int querySetKey, unsigned int targetSetKey) {
    unsigned int orfCount = (unsigned int) std::stoul(querySizeReader->getDataByDBKey(querySetKey));
    unsigned targetGeneCount = (unsigned int) std::stoul(targetSizeReader->getDataByDBKey(targetSetKey));

    double pvalThreshold = alpha / targetGeneCount;
    double r = 0;
    double pvalue;
    size_t k = 0;
    for (size_t i = 0; i < dataToAggregate.size(); ++i) {
        pvalue = std::strtod(dataToAggregate[i][1].c_str(), NULL);
        if (pvalue < pvalThreshold) {
            k++;
            r -= log(pvalue / pvalThreshold);
        }
    }

    double totalSum = 0;
    double rightSum = 0;
    for (size_t i = 0; i < k; i++) {
        // LeftSum
        rightSum = 0.0;
        for (size_t j = i + 1; j < k + 1; j++) {
            // RightSum
            //BinCoeff(orfCount, j) * pow(P0, j) * pow((1.0 - P0), (orfCount - j));
            rightSum += exp(LBinCoeff(orfCount, j) + j * log(pvalThreshold) +
                            (orfCount - j) * log(1.0 - pvalThreshold));
        }
        // FIXME: faster
        totalSum += (pow(r, i) / factorial(i)) * rightSum;
    }

    double I = 0;
    if (r == 0) {
        I = 1;
    }

    double updatedPval = (1.0 - pow((1.0 - pvalThreshold), orfCount)) * I + exp(-r) * totalSum;
    double updatedEval = updatedPval * orfCount;

    if (std::isinf(r)) {
        Debug(Debug::WARNING) << "R is infinity !\n";
        updatedEval = 0;
    }

    if (updatedEval == 0) {
        Debug(Debug::WARNING) << "Eval is 0 !\n";
        pvalue = std::strtod(dataToAggregate.back()[1].c_str(), NULL);
    }

    std::string buffer;
    // FIXME: Itoa
    buffer.append(SSTR(targetSetKey));
    buffer.append("\t");
    buffer.append(SSTR(updatedEval));
    return buffer;
}

HitDistanceAggregator::HitDistanceAggregator(const std::string &queryDbName, const std::string &targetDbName,
                                             const std::string &resultDbName, const std::string &outputDbName, bool shortOutput,
                                             float alpha, unsigned int threads)
        : Aggregation(resultDbName, targetDbName, 12, threads), alpha(alpha), shortOutput(shortOutput) {
\
    std::string data = queryDbName + "_orfs_size";
    std::string index = queryDbName + "_orfs_size.index";
    querySizeReader = new DBReader<unsigned int>(data.c_str(), index.c_str());
    querySizeReader->open(DBReader<unsigned int>::NOSORT);

    data = targetDbName + "_orfs_size";
    index = targetDbName + "_orfs_size.index";
    targetSizeReader = new DBReader<unsigned int>(data.c_str(), index.c_str());
    targetSizeReader->open(DBReader<unsigned int>::NOSORT);

    data = targetDbName + "_nucl_size";
    index = targetDbName + "_nucl_size.index";
    targetSourceReader = new DBReader<unsigned int>(data.c_str(), index.c_str());
    targetSourceReader->open(DBReader<unsigned int>::USE_INDEX);
}

HitDistanceAggregator::~HitDistanceAggregator() {
    targetSourceReader->close();
    delete targetSourceReader;

    targetSizeReader->close();
    delete targetSizeReader;

    querySizeReader->close();
    delete querySizeReader;
}

double geomMedianProba(size_t i, size_t K, double rate) {
    return pow((1.0 - pow(1.0 - rate, i)) * pow(1.0 - rate, (i + 1)), (K / 2));
}

double normalizedSpreadPval(size_t median, size_t K, double rate) {
    double s = 0.0;

    if (median >= 1.0 / rate) {
        return 1.0;
    }

    for (size_t i = 0; i <= median; i++) {
        s += geomMedianProba(i, K, rate);
    }

    double sNorm = s;
    for (size_t i = median + 1; i <= 1.0 / rate; i++) {
        sNorm += geomMedianProba(i, K, rate);
    }

    // would need *exp(LBinCoeff(K,K/2)); but does not matter since we normalize
    return s / sNorm;
}

double spreadPval(size_t median, double rate, size_t K) {
    if (median < 1) {
        median = 1;
    }

    if (median > 100000) {
        return 1.0;
    }

    return normalizedSpreadPval(median, K, rate);
}


// return the median of the distance between genes of dataToAggregate
std::string HitDistanceAggregator::aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, unsigned int querySetKey, unsigned int targetSetKey) {
    double targetGeneCount = std::strtod(targetSizeReader->getDataByDBKey(targetSetKey), NULL);
    double pvalThreshold = this->alpha / targetGeneCount;
    std::vector<std::pair<long, long>> genesPositions;
    size_t hitsUnderThreshold = 0;
    double meanEval = 0;
    std::string eVals;
    std::string genesID;
    std::string positionsStr;
    unsigned int nbrGoodEvals = 0;
    for (auto &it : dataToAggregate) {
        double Pval = std::strtod(it[3].c_str(), nullptr);
        if (Pval < pvalThreshold) {
            unsigned long start = static_cast<unsigned long>(std::stol(it[8]));
            unsigned long stop = static_cast<unsigned long>(std::stol(it[10]));
            genesPositions.emplace_back(std::make_pair(start, stop));
            hitsUnderThreshold++;
            if (shortOutput == false) {
                meanEval += log10(std::strtod(it[3].c_str(), nullptr));
                eVals += it[3] + ",";
                genesID += it[0] + ",";
                positionsStr += std::to_string(start) + "," + std::to_string(stop) + ",";
                if (std::strtod(it[3].c_str(), nullptr) < 1e-10) {
                    nbrGoodEvals++;
                }
            }
        }
    }

    std::sort(genesPositions.begin(), genesPositions.end(), compareByStart());
    double genomeSize = (targetSourceReader->getSeqLens(targetSourceReader->getId(targetSetKey)) - 2);
    double rate = ((double) hitsUnderThreshold) / genomeSize;

    std::string buffer;
    if (hitsUnderThreshold > 1) {
        std::vector<long> interGeneSpaces;
        for (size_t i = 0; i < hitsUnderThreshold - 1; i++) {
            size_t currentInterGenePosition = genesPositions[i + 1].first - genesPositions[i].second;
            interGeneSpaces.push_back(currentInterGenePosition);
        }
        std::sort(interGeneSpaces.begin(), interGeneSpaces.end());

        //if odd number
        unsigned int interSpaceIndex = 0;
        if (interGeneSpaces.size() == 1) {
            interSpaceIndex = 0;
        } else if (interGeneSpaces.size() % 2) {
            interSpaceIndex = static_cast<unsigned int>((interGeneSpaces.size() + 1) / 2);
        } else {
            interSpaceIndex = static_cast<unsigned int>(interGeneSpaces.size() / 2);
        }

        if (shortOutput == false) {
            // T    HitDist     nbrGeneT        nbrHits     nbrGeneQ    spreadPval      genomeSize      nbrGoodEvals
            buffer.append(std::to_string(targetSetKey)); // TODO ito
            buffer.append("\t");
            buffer.append(std::to_string(interGeneSpaces[interSpaceIndex]));
            buffer.append("\t");
            buffer.append(targetSizeReader->getDataByDBKey(targetSetKey));
            buffer.append("\t");
            buffer.append(std::to_string(hitsUnderThreshold));
            buffer.append("\t");
            buffer.append(std::to_string(meanEval / hitsUnderThreshold));
            buffer.append("\t");
            buffer.append(querySizeReader->getDataByDBKey(querySetKey));
            buffer.append("\t");
            buffer.append(std::to_string(spreadPval(interGeneSpaces[interSpaceIndex], rate, hitsUnderThreshold)));
            buffer.append("\t");
            buffer.append(std::to_string(genomeSize));
            buffer.append("\t");
            buffer.append(std::to_string(nbrGoodEvals));
        } else {
            buffer.append(std::to_string(targetSetKey));
            buffer.append("\t");
            buffer.append(std::to_string(spreadPval(interGeneSpaces[interSpaceIndex], rate, hitsUnderThreshold)));
        }
    } else {
        if (shortOutput == false) {
            buffer.append(std::to_string(targetSetKey));
            buffer.append("\t0\t");
            buffer.append(targetSizeReader->getDataByDBKey(targetSetKey));
            buffer.append("\t");
            buffer.append(std::to_string(hitsUnderThreshold));
            buffer.append("\t");
            buffer.append(std::to_string(meanEval / hitsUnderThreshold ));
            buffer.append("\t");
            buffer.append(querySizeReader->getDataByDBKey(querySetKey));
            buffer.append("\t1.0\t");
            buffer.append(std::to_string(genomeSize));
            buffer.append("\t");
            buffer.append(std::to_string(nbrGoodEvals)); //NaMM

        } else {
            buffer.append(std::to_string(targetSetKey));
            buffer.append("\t0");
        }
    }
    if (shortOutput == false) {
        buffer.append("\t");
        buffer.append(positionsStr);
        buffer.append("\t");
        buffer.append(genesID);
        buffer.append("\t");
        buffer.append(eVals);
    }
    return buffer;
}
