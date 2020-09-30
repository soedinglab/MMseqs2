#include "Debug.h"
#include "Parameters.h"
#include "FastSort.h"
#include "Aggregation.h"

#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

struct compareByStart {
    bool operator()(const std::pair<long, long> &lhs,
                    const std::pair<long, long> &rhs) const {
        return (lhs.first < rhs.first);
    }
};

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

class SetSummaryAggregator : public Aggregation {
public:
    SetSummaryAggregator(const std::string &queryDbName, const std::string &targetDbName,
                         const std::string &resultDbName, const std::string &outputDbName, bool shortOutput,
                         float alpha, unsigned int threads, unsigned int compressed)
            : Aggregation(targetDbName, resultDbName, outputDbName, threads, compressed), alpha(alpha), shortOutput(shortOutput) {
        std::string data = queryDbName + "_set_size";
        std::string index = queryDbName + "_set_size.index";
        querySizeReader = new DBReader<unsigned int>(data.c_str(), index.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        querySizeReader->open(DBReader<unsigned int>::NOSORT);

        data = targetDbName + "_set_size";
        index = targetDbName + "_set_size.index";
        targetSizeReader = new DBReader<unsigned int>(data.c_str(), index.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        targetSizeReader->open(DBReader<unsigned int>::NOSORT);

        data = targetDbName + "_nucl";
        index = targetDbName + "_nucl.index";
        targetSourceReader = new DBReader<unsigned int>(data.c_str(), index.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        targetSourceReader->open(DBReader<unsigned int>::USE_INDEX);
    }

    ~SetSummaryAggregator() {
        targetSourceReader->close();
        delete targetSourceReader;

        targetSizeReader->close();
        delete targetSizeReader;

        querySizeReader->close();
        delete querySizeReader;
    }


    void prepareInput(unsigned int, unsigned int) {}

    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, unsigned int querySetKey,
                               unsigned int targetSetKey, unsigned int thread_idx) {
        double targetGeneCount = std::strtod(targetSizeReader->getDataByDBKey(targetSetKey, thread_idx), NULL);
        double pvalThreshold = this->alpha / targetGeneCount;
        std::vector<std::pair<long, long>> genesPositions;
        size_t hitsUnderThreshold = 0;
        double meanEval = 0;
        std::string eVals;
        std::string genesID;
        std::string positionsStr;
        unsigned int nbrGoodEvals = 0;
        for (std::vector<std::vector<std::string>>::const_iterator it = dataToAggregate.begin();
             it != dataToAggregate.end(); ++it) {
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

        SORT_SERIAL(genesPositions.begin(), genesPositions.end(), compareByStart());

        // TODO: Get size for whole genome (multiple contigs, proteines etc.)
        double genomeSize = (targetSourceReader->getSeqLen(targetSourceReader->getId(targetSetKey)));
        double rate = ((double) hitsUnderThreshold) / genomeSize;

        std::string buffer;
        if (hitsUnderThreshold > 1) {
            std::vector<long> interGeneSpaces;
            for (size_t i = 0; i < hitsUnderThreshold - 1; i++) {
                size_t currentInterGenePosition = genesPositions[i + 1].first - genesPositions[i].second;
                interGeneSpaces.push_back(currentInterGenePosition);
            }
            SORT_SERIAL(interGeneSpaces.begin(), interGeneSpaces.end());

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
                buffer.append(targetSizeReader->getDataByDBKey(targetSetKey, thread_idx));
                buffer.append("\t");
                buffer.append(std::to_string(hitsUnderThreshold));
                buffer.append("\t");
                buffer.append(std::to_string(meanEval / hitsUnderThreshold));
                buffer.append("\t");
                buffer.append(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
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
                buffer.append(targetSizeReader->getDataByDBKey(targetSetKey, thread_idx));
                buffer.append("\t");
                buffer.append(std::to_string(hitsUnderThreshold));
                buffer.append("\t");
                buffer.append(std::to_string(meanEval / hitsUnderThreshold));
                buffer.append("\t");
                buffer.append(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
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

private:
    DBReader<unsigned int> *querySizeReader;
    DBReader<unsigned int> *targetSourceReader;
    DBReader<unsigned int> *targetSizeReader;
    float alpha;
    bool shortOutput;
};

int resultsbyset(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    SetSummaryAggregator aggregation(par.db1, par.db2, par.db3, par.db4, par.shortOutput, par.alpha, (unsigned int) par.threads, par.compressed);
    return aggregation.run();
}
