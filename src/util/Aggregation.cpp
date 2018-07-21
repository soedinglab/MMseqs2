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

Aggregation::Aggregation(const std::string &targetDbName, const std::string &resultDbName, const std::string &outputDbName, unsigned int threads)
        : resultDbName(resultDbName), outputDbName(outputDbName), threads(threads) {
    std::string sizeDbName = targetDbName + "_set_lookup";
    std::string sizeDbIndex = targetDbName + "_set_lookup.index";
    targetSetReader = new DBReader<unsigned int>(sizeDbName.c_str(), sizeDbIndex.c_str());
    targetSetReader->open(DBReader<unsigned int>::NOSORT);
}

// build a map with the value in [target column] field as a key and the rest of the line, cut in fields, as values
bool Aggregation::buildMap(std::stringstream &dataStream, std::map<unsigned int, std::vector<std::vector<std::string>>> &dataToAggregate) {
    dataToAggregate.clear();
    std::string line;
    while (std::getline(dataStream, line)) {
        if (!line.empty()) {
            std::vector<std::string> columns = Util::split(line, "\t");
            unsigned int targetKey = (unsigned int) strtoull(columns[0].c_str(), NULL, 10);
            size_t setId = targetSetReader->getId(targetKey);
            if (setId == UINT_MAX) {
                Debug(Debug::ERROR) << "Invalid target database key " << targetKey << ".\n";
                EXIT(EXIT_FAILURE);
            }
            char *data = targetSetReader->getData(setId);
            unsigned int setKey = (unsigned int)strtoull(data, NULL, 10);
            dataToAggregate[setKey].push_back(columns);
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


BestHitBySetFilter::BestHitBySetFilter(const std::string &targetDbName, const std::string &resultDbName,
                                     const std::string &outputDbName, bool simpleBestHitMode, unsigned int threads) :
        Aggregation(targetDbName, resultDbName, outputDbName, threads), simpleBestHitMode(simpleBestHitMode) {
    std::string sizeDbName = targetDbName + "_set_size";
    std::string sizeDbIndex = targetDbName + "_set_size.index";
    targetSizeReader = new DBReader<unsigned int>(sizeDbName.c_str(), sizeDbIndex.c_str());
    targetSizeReader->open(DBReader<unsigned int>::NOSORT);
}

// Only Keep the best hits of a protein against each Target Set
std::string BestHitBySetFilter::aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, unsigned int, unsigned int targetSetKey) {
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
        double eval = strtod(dataToAggregate[i][EVAL_COLUMN].c_str(), NULL);

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
    unsigned int nbrGenes = (unsigned int) std::strtoull(data, NULL, 10);

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


PvalueAggregator::PvalueAggregator(std::string queryDbName, std::string targetDbName, const std::string &resultDbName,
                                   const std::string &outputDbName, float alpha, unsigned int threads) :
        Aggregation(targetDbName, resultDbName, outputDbName, threads), alpha(alpha) {

    std::string sizeDBName  = queryDbName + "_set_size";
    std::string sizeDBIndex = queryDbName + "_set_size.index";
    querySizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str());
    querySizeReader->open(DBReader<unsigned int>::NOSORT);

    sizeDBName = targetDbName + "_set_size";
    sizeDBIndex = targetDbName + "_set_size.index";
    targetSizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str());
    targetSizeReader->open(DBReader<unsigned int>::NOSORT);
}

double LBinCoeff(int M, int k) {
    return lgamma(M + 1) - lgamma(M - k + 1) - lgamma(k + 1);
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

    double updatedEval = 0;
    if (std::isinf(r) == false) {
        double truncatedFisherPval = 0;
        double powR = 1;
        size_t factorial = 1;
        for (size_t i = 0; i < k; i++) {
            // K hits with p-val < p0
            double probKGoodHits = 0.0;
            for (size_t j = i + 1; j < k + 1; j++) {
                probKGoodHits += exp(LBinCoeff(orfCount, j) + j * log(pvalThreshold) +
                                (orfCount - j) * log(1.0 - pvalThreshold));
            }
            if (i > 0) {
                factorial *= i;
            }
            truncatedFisherPval += (powR / factorial) * probKGoodHits;
            powR *= r;
        }

        double I = 0;
        if (r == 0) {
            I = 1;
        }

        double updatedPval = (1.0 - pow((1.0 - pvalThreshold), orfCount)) * I + exp(-r) * truncatedFisherPval;
        updatedEval = updatedPval * targetSizeReader->getSize();
    }

    std::string buffer;
    // FIXME: Itoa
    buffer.append(SSTR(targetSetKey));
    buffer.append("\t");
    buffer.append(SSTR(updatedEval));
    return buffer;
}

SetSummaryAggregator::SetSummaryAggregator(const std::string &queryDbName, const std::string &targetDbName,
                                             const std::string &resultDbName, const std::string &outputDbName, bool shortOutput,
                                             float alpha, unsigned int threads)
        : Aggregation(targetDbName, resultDbName, outputDbName, threads), alpha(alpha), shortOutput(shortOutput) {
    std::string data = queryDbName + "_set_size";
    std::string index = queryDbName + "_set_size.index";
    querySizeReader = new DBReader<unsigned int>(data.c_str(), index.c_str());
    querySizeReader->open(DBReader<unsigned int>::NOSORT);

    data = targetDbName + "_set_size";
    index = targetDbName + "_set_size.index";
    targetSizeReader = new DBReader<unsigned int>(data.c_str(), index.c_str());
    targetSizeReader->open(DBReader<unsigned int>::NOSORT);

    data = targetDbName + "_nucl";
    index = targetDbName + "_nucl.index";
    targetSourceReader = new DBReader<unsigned int>(data.c_str(), index.c_str());
    targetSourceReader->open(DBReader<unsigned int>::USE_INDEX);
}

SetSummaryAggregator::~SetSummaryAggregator() {
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
std::string SetSummaryAggregator::aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate, unsigned int querySetKey, unsigned int targetSetKey) {
    double targetGeneCount = std::strtod(targetSizeReader->getDataByDBKey(targetSetKey), NULL);
    double pvalThreshold = this->alpha / targetGeneCount;
    std::vector<std::pair<long, long>> genesPositions;
    size_t hitsUnderThreshold = 0;
    double meanEval = 0;
    std::string eVals;
    std::string genesID;
    std::string positionsStr;
    unsigned int nbrGoodEvals = 0;
    for (std::vector<std::vector<std::string>>::const_iterator it = dataToAggregate.begin(); it != dataToAggregate.end(); ++it) {
        double Pval = std::strtod((*it)[3].c_str(), NULL);
        if (Pval >= pvalThreshold) {
            continue;
        }

        unsigned long start = static_cast<unsigned long>(std::stol((*it)[8]));
        unsigned long stop = static_cast<unsigned long>(std::stol((*it)[10]));
        genesPositions.emplace_back(std::make_pair(start, stop));
        hitsUnderThreshold++;

        if (shortOutput) {
            continue;
        }
        meanEval += log10(std::strtod((*it)[3].c_str(), NULL));
        eVals += (*it)[3] + ",";
        genesID += (*it)[0] + ",";
        positionsStr += std::to_string(start) + "," + std::to_string(stop) + ",";
        if (std::strtod((*it)[3].c_str(), NULL) < 1e-10) {
            nbrGoodEvals++;
        }
    }

    std::sort(genesPositions.begin(), genesPositions.end(), compareByStart());

    // TODO: Get size for whole genome (multiple contigs, proteines etc.)
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
