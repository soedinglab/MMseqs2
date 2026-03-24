#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include "FastSort.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
struct ReclassEntry {
    Matcher::result_t result;
    double abundance;
    double posterior;
    double entropyPenalty;
};

typedef std::unordered_map<unsigned int, std::vector<ReclassEntry> > MappingTable;

struct ReclassContext {
    MappingTable mappingTable;
    std::unordered_set<unsigned int> targetSet;
    size_t queryCount;
    bool hasBacktrace;
    bool hasOrfPosition;

    ReclassContext() : queryCount(0), hasBacktrace(false), hasOrfPosition(false) {}
};

static const double DEFAULT_LAMBDA = 0.02;
static const double DEFAULT_ALPHA = 1.0;
static const double DEFAULT_BETA = 1.0;
static const double DEFAULT_GAMMA = 1.0;
static const int DEFAULT_MAX_ITER = 100;
static const double DEFAULT_TOL = 1e-5;
static const double STEP_MIN = -1.0;
static const double STEP_MAX = 1.0;
static const double EPS = 1e-12;

static std::vector<unsigned int> targetListFromSet(const std::unordered_set<unsigned int> &targets) {
    std::vector<unsigned int> out(targets.begin(), targets.end());
    SORT_SERIAL(out.begin(), out.end());
    return out;
}

static void loadAlignmentDb(DBReader<unsigned int> &reader, ReclassContext &ctx) {
    Debug::Progress progress(reader.getSize());
    const char *entry[255];

    for (size_t i = 0; i < reader.getSize(); ++i) {
        progress.updateProgress();
        const unsigned int queryKey = reader.getDbKey(i);
        char *data = reader.getData(i, 0);

        if (reader.getEntryLen(i) <= 1) {
            continue;
        }

        std::vector<ReclassEntry> &records = ctx.mappingTable[queryKey];
        while (*data != '\0') {
            const size_t columns = Util::getWordsOfLine(data, entry, 255);
            if (columns < Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                Debug(Debug::ERROR) << "Invalid alignment result record in query " << queryKey << ".\n";
                EXIT(EXIT_FAILURE);
            }

            if (columns == Matcher::ALN_RES_WITH_BT_COL_CNT || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT) {
                ctx.hasBacktrace = true;
            }
            if (columns == Matcher::ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT) {
                ctx.hasOrfPosition = true;
            }

            Matcher::result_t result = Matcher::parseAlignmentRecord(data, true);
            records.push_back(ReclassEntry{result, 0.0, 0.0, 0.0});
            ctx.targetSet.insert(result.dbKey);
            data = Util::skipLine(data);
        }
    }

    ctx.queryCount = ctx.mappingTable.size();
}

static void initAbundance(MappingTable &mappingTable, const std::unordered_set<unsigned int> &targetSet, size_t queryCount) {
    std::unordered_map<unsigned int, double> initAbundance;
    initAbundance.reserve(targetSet.size());
    for (std::unordered_set<unsigned int>::const_iterator it = targetSet.begin(); it != targetSet.end(); ++it) {
        initAbundance[*it] = 0.0;
    }

    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        double scoreSum = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            scoreSum += it->second[j].result.score;
        }
        if (scoreSum <= 0.0) {
            continue;
        }
        for (size_t j = 0; j < it->second.size(); ++j) {
            initAbundance[it->second[j].result.dbKey] += it->second[j].result.score / scoreSum;
        }
    }

    if (queryCount > 0) {
        const double denom = static_cast<double>(queryCount);
        for (std::unordered_map<unsigned int, double>::iterator it = initAbundance.begin(); it != initAbundance.end(); ++it) {
            it->second /= denom;
        }
    }

    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            it->second[j].abundance = initAbundance[it->second[j].result.dbKey];
        }
    }
}

static void initEntropy(MappingTable &mappingTable, const std::unordered_set<unsigned int> &targetSet, double lambda) {
    std::unordered_map<unsigned int, int> targetMin;
    std::unordered_map<unsigned int, int> targetMax;

    for (std::unordered_set<unsigned int>::const_iterator it = targetSet.begin(); it != targetSet.end(); ++it) {
        targetMin[*it] = std::numeric_limits<int>::max();
        targetMax[*it] = std::numeric_limits<int>::min();
    }

    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            const unsigned int target = it->second[j].result.dbKey;
            if (it->second[j].result.dbStartPos < targetMin[target]) {
                targetMin[target] = it->second[j].result.dbStartPos;
            }
            if (it->second[j].result.dbEndPos > targetMax[target]) {
                targetMax[target] = it->second[j].result.dbEndPos;
            }
        }
    }

    std::unordered_map<unsigned int, std::vector<double> > coverage;
    coverage.reserve(targetSet.size());
    for (std::unordered_set<unsigned int>::const_iterator it = targetSet.begin(); it != targetSet.end(); ++it) {
        const int start = targetMin[*it];
        const int end = targetMax[*it];
        const int len = (end >= start) ? (end - start + 1) : 1;
        coverage[*it] = std::vector<double>(len, 0.0);
    }

    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            const Matcher::result_t &result = it->second[j].result;
            const int targetLen = result.dbEndPos - result.dbStartPos + 1;
            if (targetLen <= 0) {
                continue;
            }

            const double mq = std::exp(lambda * static_cast<double>(result.score)) / static_cast<double>(targetLen);
            std::vector<double> &cov = coverage[result.dbKey];
            const int start = std::max(0, result.dbStartPos - targetMin[result.dbKey]);
            const int end = std::min(static_cast<int>(cov.size()) - 1, result.dbEndPos - targetMin[result.dbKey]);
            for (int pos = start; pos <= end; ++pos) {
                cov[pos] += mq;
            }
        }
    }

    std::unordered_map<unsigned int, double> entropy;
    entropy.reserve(targetSet.size());
    for (std::unordered_set<unsigned int>::const_iterator it = targetSet.begin(); it != targetSet.end(); ++it) {
        const std::vector<double> &cov = coverage[*it];
        const double covSum = std::accumulate(cov.begin(), cov.end(), 0.0);
        if (covSum <= 0.0) {
            entropy[*it] = 0.0;
            continue;
        }

        double ent = 0.0;
        for (size_t pos = 0; pos < cov.size(); ++pos) {
            if (cov[pos] <= 0.0) {
                continue;
            }
            const double p = cov[pos] / covSum;
            ent -= p * std::log2(p);
        }
        entropy[*it] = ent;
    }

    double entropySum = 0.0;
    for (std::unordered_map<unsigned int, double>::const_iterator it = entropy.begin(); it != entropy.end(); ++it) {
        entropySum += it->second;
    }

    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            const unsigned int target = it->second[j].result.dbKey;
            it->second[j].entropyPenalty = (entropySum > 0.0) ? (1.0 - (entropy[target] / entropySum)) : 0.0;
        }
    }
}

static double scoreTerm(const ReclassEntry &entry, double lambda, double alpha, double beta, double gamma) {
    if (entry.result.alnLength == 0) {
        return 0.0;
    }

    const double bitPerLen = static_cast<double>(entry.result.score) / static_cast<double>(entry.result.alnLength);
    const double seqId = std::max(static_cast<double>(entry.result.seqId) / 100.0, EPS);
    const double abundance = std::max(entry.abundance, EPS);
    const double entropyPenalty = std::max(entry.entropyPenalty, EPS);
    return std::exp(lambda * bitPerLen) * std::pow(seqId, beta) * std::pow(abundance, alpha) * std::pow(entropyPenalty, gamma);
}

static void computePosterior(MappingTable &mappingTable, double lambda, double alpha, double beta, double gamma) {
    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        double denom = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            denom += scoreTerm(it->second[j], lambda, alpha, beta, gamma);
        }
        for (size_t j = 0; j < it->second.size(); ++j) {
            const double value = scoreTerm(it->second[j], lambda, alpha, beta, gamma);
            it->second[j].posterior = (denom > 0.0) ? (value / denom) : 0.0;
        }
    }
}

static double logLikelihood(const MappingTable &mappingTable, double lambda, double alpha, double beta, double gamma, size_t queryCount) {
    if (queryCount == 0) {
        return 0.0;
    }

    double ll = 0.0;
    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        double denom = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            const double value = scoreTerm(it->second[j], lambda, alpha, beta, gamma);
            denom += value;
            if (it->second[j].posterior > 0.0) {
                ll += it->second[j].posterior * std::log(value > 0.0 ? value : 1e-300);
            }
        }
        ll -= std::log(denom > 0.0 ? denom : 1e-300);
    }
    return ll / static_cast<double>(queryCount);
}

static std::vector<double> abundanceVectorFromTable(const MappingTable &mappingTable, const std::vector<unsigned int> &targetList) {
    std::unordered_map<unsigned int, double> abundance;
    abundance.reserve(targetList.size());
    for (size_t i = 0; i < targetList.size(); ++i) {
        abundance[targetList[i]] = 0.0;
    }

    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            abundance[it->second[j].result.dbKey] = it->second[j].abundance;
        }
    }

    std::vector<double> out;
    out.reserve(targetList.size());
    for (size_t i = 0; i < targetList.size(); ++i) {
        out.push_back(abundance[targetList[i]]);
    }
    return out;
}

static void setAbundance(MappingTable &mappingTable, const std::vector<unsigned int> &targetList, const std::vector<double> &abundanceVector) {
    std::unordered_map<unsigned int, double> abundance;
    abundance.reserve(targetList.size());
    for (size_t i = 0; i < targetList.size(); ++i) {
        abundance[targetList[i]] = abundanceVector[i];
    }

    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            it->second[j].abundance = abundance[it->second[j].result.dbKey];
        }
    }
}

static std::vector<double> emUpdate(MappingTable &mappingTable,
                                    double lambda,
                                    const std::unordered_map<unsigned int, double> &fixedEntropy,
                                    const std::vector<unsigned int> &targetList,
                                    size_t queryCount,
                                    double alpha,
                                    double beta,
                                    double gamma) {
    computePosterior(mappingTable, lambda, alpha, beta, gamma);

    std::unordered_map<unsigned int, double> nextAbundance;
    nextAbundance.reserve(targetList.size());
    for (size_t i = 0; i < targetList.size(); ++i) {
        nextAbundance[targetList[i]] = 0.0;
    }

    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            nextAbundance[it->second[j].result.dbKey] += it->second[j].posterior;
        }
    }

    if (queryCount > 0) {
        const double denom = static_cast<double>(queryCount);
        for (std::unordered_map<unsigned int, double>::iterator it = nextAbundance.begin(); it != nextAbundance.end(); ++it) {
            it->second /= denom;
        }
    }

    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            const unsigned int target = it->second[j].result.dbKey;
            it->second[j].abundance = nextAbundance[target];
            std::unordered_map<unsigned int, double>::const_iterator fixed = fixedEntropy.find(target);
            it->second[j].entropyPenalty = (fixed != fixedEntropy.end()) ? fixed->second : 0.0;
        }
    }

    std::vector<double> out;
    out.reserve(targetList.size());
    for (size_t i = 0; i < targetList.size(); ++i) {
        out.push_back(nextAbundance[targetList[i]]);
    }
    return out;
}

static std::vector<double> projectSimplex(const std::vector<double> &x) {
    std::vector<double> projected = x;
    for (size_t i = 0; i < projected.size(); ++i) {
        if (projected[i] < 0.0) {
            projected[i] = 0.0;
        }
    }

    const double sum = std::accumulate(projected.begin(), projected.end(), 0.0);
    if (sum > 0.0) {
        for (size_t i = 0; i < projected.size(); ++i) {
            projected[i] /= sum;
        }
    }
    return projected;
}

static void squarem(ReclassContext &ctx, double lambda, int maxIter, double tol, double alpha, double beta, double gamma) {
    if (ctx.queryCount == 0 || ctx.targetSet.empty()) {
        return;
    }

    initAbundance(ctx.mappingTable, ctx.targetSet, ctx.queryCount);
    initEntropy(ctx.mappingTable, ctx.targetSet, lambda);

    std::unordered_map<unsigned int, double> fixedEntropy;
    fixedEntropy.reserve(ctx.targetSet.size());
    for (MappingTable::const_iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            fixedEntropy[it->second[j].result.dbKey] = it->second[j].entropyPenalty;
        }
    }

    const std::vector<unsigned int> targetList = targetListFromSet(ctx.targetSet);
    std::vector<double> x0 = abundanceVectorFromTable(ctx.mappingTable, targetList);
    std::vector<double> logLikelihoods;

    for (int iter = 0; iter < maxIter; ++iter) {
        const std::vector<double> x1 = emUpdate(ctx.mappingTable, lambda, fixedEntropy, targetList, ctx.queryCount, alpha, beta, gamma);
        const std::vector<double> x2 = emUpdate(ctx.mappingTable, lambda, fixedEntropy, targetList, ctx.queryCount, alpha, beta, gamma);

        std::vector<double> r(x0.size(), 0.0);
        std::vector<double> v(x0.size(), 0.0);
        for (size_t i = 0; i < x0.size(); ++i) {
            r[i] = x1[i] - x0[i];
            v[i] = x2[i] - x1[i] - r[i];
        }

        double normR = 0.0;
        double normV = 0.0;
        for (size_t i = 0; i < x0.size(); ++i) {
            normR += r[i] * r[i];
            normV += v[i] * v[i];
        }
        normR = std::sqrt(normR);
        normV = std::sqrt(normV);

        double accel = (normV == 0.0) ? -1.0 : -(normR / normV);
        accel = std::max(STEP_MIN, std::min(STEP_MAX, accel));

        std::vector<double> xNew(x0.size(), 0.0);
        for (size_t i = 0; i < x0.size(); ++i) {
            xNew[i] = x0[i] - 2.0 * accel * r[i] + accel * accel * v[i];
        }
        xNew = projectSimplex(xNew);

        setAbundance(ctx.mappingTable, targetList, xNew);
        computePosterior(ctx.mappingTable, lambda, alpha, beta, gamma);
        double currentLl = logLikelihood(ctx.mappingTable, lambda, alpha, beta, gamma, ctx.queryCount);

        if (!logLikelihoods.empty() && currentLl < logLikelihoods.back() - 1e-9) {
            setAbundance(ctx.mappingTable, targetList, x2);
            computePosterior(ctx.mappingTable, lambda, alpha, beta, gamma);
            currentLl = logLikelihood(ctx.mappingTable, lambda, alpha, beta, gamma, ctx.queryCount);
            xNew = x2;
        }

        logLikelihoods.push_back(currentLl);

        double parameterChange = 0.0;
        for (size_t i = 0; i < x0.size(); ++i) {
            parameterChange = std::max(parameterChange, std::fabs(xNew[i] - x0[i]));
        }

        Debug(Debug::INFO) << "Reclassify iteration " << iter << ": LL=" << currentLl << " delta=" << parameterChange << "\n";
        x0 = xNew;
        if (parameterChange < tol && iter > 5) {
            Debug(Debug::INFO) << "Reclassify converged after " << (iter + 1) << " iterations.\n";
            break;
        }
    }
}

static bool compareByPosterior(const ReclassEntry &a, const ReclassEntry &b) {
    if (a.posterior != b.posterior) {
        return a.posterior > b.posterior;
    }
    return Matcher::compareHits(a.result, b.result);
}
}

int reclassify(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    ReclassContext ctx;
    loadAlignmentDb(reader, ctx);
    Debug(Debug::INFO) << "Loaded " << ctx.queryCount << " queries with hits and " << ctx.targetSet.size() << " unique targets.\n";

    squarem(ctx, DEFAULT_LAMBDA, DEFAULT_MAX_ITER, DEFAULT_TOL, DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_GAMMA);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    Debug::Progress progress(reader.getSize());
    char buffer[1024 + 32768 * 4];
    for (size_t i = 0; i < reader.getSize(); ++i) {
        progress.updateProgress();
        const unsigned int queryKey = reader.getDbKey(i);
        MappingTable::iterator it = ctx.mappingTable.find(queryKey);
        if (it == ctx.mappingTable.end() || it->second.empty()) {
            writer.writeData("", 0, queryKey, 0);
            continue;
        }

        SORT_SERIAL(it->second.begin(), it->second.end(), compareByPosterior);
        writer.writeStart(0);
        for (size_t j = 0; j < it->second.size(); ++j) {
            const size_t len = Matcher::resultToBuffer(buffer, it->second[j].result, ctx.hasBacktrace, false, ctx.hasOrfPosition);
            writer.writeAdd(buffer, len, 0);
        }
        writer.writeEnd(queryKey, 0);
    }

    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}
