#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include "FastSort.h"
#include "FileUtil.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cctype>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
struct ReclassTaxEntry {
    Matcher::result_t result;
    double abundance;
    double posterior;
    double entropyValue;
    double entropyPenalty;
};

typedef std::unordered_map<unsigned int, std::vector<ReclassTaxEntry> > MappingTable;

struct Interval {
    int start;
    int end;
};

struct TargetStats {
    unsigned int key;
    unsigned int targetLength;
    double abundance;
    double entropyValue;
    double entropyPenalty;
    bool dropped;
    std::vector<Interval> intervals;
};

struct TaxonomyStats {
    unsigned int taxId;
    double abundance;
    double entropySum;
    double entropyPenaltySum;
    size_t proteinCount;

    TaxonomyStats() : taxId(0), abundance(0.0), entropySum(0.0), entropyPenaltySum(0.0), proteinCount(0) {}
};

struct ReclassTaxContext {
    MappingTable mappingTable;
    std::vector<unsigned int> queryOrder;
    std::unordered_set<unsigned int> targetSet;
    size_t queryCount;
    bool hasBacktrace;
    bool hasOrfPosition;

    ReclassTaxContext() : queryCount(0), hasBacktrace(false), hasOrfPosition(false) {}
};

static const double STEP_MIN = -1.0;
static const double STEP_MAX = 1.0;
static const double EPS = 1e-12;

static std::vector<unsigned int> targetListFromSet(const std::unordered_set<unsigned int> &targets) {
    std::vector<unsigned int> out(targets.begin(), targets.end());
    SORT_SERIAL(out.begin(), out.end());
    return out;
}

static void loadAlignmentDb(DBReader<unsigned int> &reader, ReclassTaxContext &ctx) {
    Debug::Progress progress(reader.getSize());
    const char *entry[255];

    for (size_t i = 0; i < reader.getSize(); ++i) {
        progress.updateProgress();
        const unsigned int queryKey = reader.getDbKey(i);
        char *data = reader.getData(i, 0);

        if (reader.getEntryLen(i) <= 1) {
            continue;
        }

        std::vector<ReclassTaxEntry> &records = ctx.mappingTable[queryKey];
        if (records.empty()) {
            ctx.queryOrder.push_back(queryKey);
        }
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
            records.push_back(ReclassTaxEntry{result, 0.0, 0.0, 0.0, 0.0});
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
            it->second[j].entropyValue = entropy[target];
            it->second[j].entropyPenalty = (entropySum > 0.0) ? (1.0 - (entropy[target] / entropySum)) : 0.0;
        }
    }
}

static double maxQueryBitScore(const std::vector<ReclassTaxEntry> &entries) {
    double maxScore = 0.0;
    for (size_t i = 0; i < entries.size(); ++i) {
        maxScore = std::max(maxScore, static_cast<double>(entries[i].result.score));
    }
    return maxScore;
}

static double scoreTerm(const ReclassTaxEntry &entry, double queryMaxScore, double lambda, double alpha, double gamma) {
    if (queryMaxScore <= 0.0) {
        return 0.0;
    }

    const double normalizedScore = static_cast<double>(entry.result.score) / queryMaxScore;
    const double abundance = std::max(entry.abundance, EPS);
    const double entropyPenalty = std::max(entry.entropyPenalty, EPS);
    return std::exp(lambda * normalizedScore) * std::pow(abundance, alpha) * std::pow(entropyPenalty, gamma);
}

static void computePosterior(MappingTable &mappingTable, double lambda, double alpha, double gamma) {
    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        const double queryMaxScore = maxQueryBitScore(it->second);
        double denom = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            denom += scoreTerm(it->second[j], queryMaxScore, lambda, alpha, gamma);
        }
        for (size_t j = 0; j < it->second.size(); ++j) {
            const double value = scoreTerm(it->second[j], queryMaxScore, lambda, alpha, gamma);
            it->second[j].posterior = (denom > 0.0) ? (value / denom) : 0.0;
        }
    }
}

static double logLikelihood(const MappingTable &mappingTable, double lambda, double alpha, double gamma, size_t queryCount) {
    if (queryCount == 0) {
        return 0.0;
    }

    double ll = 0.0;
    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        const double queryMaxScore = maxQueryBitScore(it->second);
        double denom = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            const double value = scoreTerm(it->second[j], queryMaxScore, lambda, alpha, gamma);
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
                                    const std::unordered_map<unsigned int, std::pair<double, double> > &fixedEntropy,
                                    const std::vector<unsigned int> &targetList,
                                    size_t queryCount,
                                    double alpha,
                                    double gamma) {
    computePosterior(mappingTable, lambda, alpha, gamma);

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
            std::unordered_map<unsigned int, std::pair<double, double> >::const_iterator fixed = fixedEntropy.find(target);
            if (fixed != fixedEntropy.end()) {
                it->second[j].entropyValue = fixed->second.first;
                it->second[j].entropyPenalty = fixed->second.second;
            } else {
                it->second[j].entropyValue = 0.0;
                it->second[j].entropyPenalty = 0.0;
            }
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

static void squarem(ReclassTaxContext &ctx, double lambda, int maxIter, double tol, double alpha, double gamma) {
    if (ctx.queryCount == 0 || ctx.targetSet.empty()) {
        return;
    }

    initAbundance(ctx.mappingTable, ctx.targetSet, ctx.queryCount);
    initEntropy(ctx.mappingTable, ctx.targetSet, lambda);

    std::unordered_map<unsigned int, std::pair<double, double> > fixedEntropy;
    fixedEntropy.reserve(ctx.targetSet.size());
    for (MappingTable::const_iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            fixedEntropy[it->second[j].result.dbKey] = std::make_pair(it->second[j].entropyValue, it->second[j].entropyPenalty);
        }
    }

    const std::vector<unsigned int> targetList = targetListFromSet(ctx.targetSet);
    std::vector<double> x0 = abundanceVectorFromTable(ctx.mappingTable, targetList);
    std::vector<double> logLikelihoods;

    for (int iter = 0; iter < maxIter; ++iter) {
        const std::vector<double> x1 = emUpdate(ctx.mappingTable, lambda, fixedEntropy, targetList, ctx.queryCount, alpha, gamma);
        const std::vector<double> x2 = emUpdate(ctx.mappingTable, lambda, fixedEntropy, targetList, ctx.queryCount, alpha, gamma);

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
        computePosterior(ctx.mappingTable, lambda, alpha, gamma);
        double currentLl = logLikelihood(ctx.mappingTable, lambda, alpha, gamma, ctx.queryCount);

        if (!logLikelihoods.empty() && currentLl < logLikelihoods.back() - 1e-9) {
            setAbundance(ctx.mappingTable, targetList, x2);
            computePosterior(ctx.mappingTable, lambda, alpha, gamma);
            currentLl = logLikelihood(ctx.mappingTable, lambda, alpha, gamma, ctx.queryCount);
            xNew = x2;
        }

        logLikelihoods.push_back(currentLl);

        double parameterChange = 0.0;
        for (size_t i = 0; i < x0.size(); ++i) {
            parameterChange = std::max(parameterChange, std::fabs(xNew[i] - x0[i]));
        }

        Debug(Debug::INFO) << "Reclassify-taxonomy iteration " << iter << ": LL=" << currentLl << " delta=" << parameterChange << "\n";
        x0 = xNew;
        if (parameterChange < tol && iter > 5) {
            Debug(Debug::INFO) << "Reclassify-taxonomy converged after " << (iter + 1) << " iterations.\n";
            break;
        }
    }
}

static bool compareByPosterior(const ReclassTaxEntry &a, const ReclassTaxEntry &b) {
    if (a.posterior != b.posterior) {
        return a.posterior > b.posterior;
    }
    return Matcher::compareHits(a.result, b.result);
}

static const char *headerForKey(DBReader<unsigned int> &headerReader, unsigned int key, unsigned int threadIdx) {
    size_t id = headerReader.getId(key);
    if (id == UINT_MAX) {
        return NULL;
    }
    return headerReader.getData(id, threadIdx);
}

static std::string identifierForKey(DBReader<unsigned int> &headerReader, unsigned int key, unsigned int threadIdx) {
    const char *header = headerForKey(headerReader, key, threadIdx);
    if (header == NULL) {
        return SSTR(key);
    }
    std::string parsed = Util::parseFastaHeader(header);
    return parsed.empty() ? SSTR(key) : parsed;
}

static void computeAlignmentCounts(const Matcher::result_t &res, unsigned int &alnLen, unsigned int &mismatchCount, unsigned int &gapOpenCount) {
    gapOpenCount = 0;
    alnLen = res.alnLength;
    mismatchCount = 0;

    if (!res.backtrace.empty()) {
        size_t matchCount = 0;
        alnLen = 0;
        for (size_t pos = 0; pos < res.backtrace.size(); ++pos) {
            int cnt = 0;
            if (std::isdigit(static_cast<unsigned char>(res.backtrace[pos]))) {
                cnt += Util::fast_atoi<int>(res.backtrace.c_str() + pos);
                while (std::isdigit(static_cast<unsigned char>(res.backtrace[pos]))) {
                    pos++;
                }
            }
            alnLen += cnt;

            switch (res.backtrace[pos]) {
                case 'M':
                    matchCount += cnt;
                    break;
                case 'D':
                case 'I':
                    gapOpenCount += 1;
                    break;
            }
        }
        const unsigned int identical = static_cast<unsigned int>(res.seqId * static_cast<float>(alnLen) + 0.5f);
        mismatchCount = static_cast<unsigned int>(matchCount - identical);
    } else {
        const int adjustQstart = (res.qStartPos == -1) ? 0 : res.qStartPos;
        const int adjustDBstart = (res.dbStartPos == -1) ? 0 : res.dbStartPos;
        const float bestMatchEstimate = static_cast<float>(std::min(abs(res.qEndPos - adjustQstart), abs(res.dbEndPos - adjustDBstart)));
        mismatchCount = static_cast<unsigned int>(bestMatchEstimate * (1.0f - res.seqId) + 0.5f);
    }
}

static void addInterval(std::vector<Interval> &intervals, int start, int end) {
    Interval interval;
    interval.start = std::min(start, end);
    interval.end = std::max(start, end);
    intervals.push_back(interval);
}

static std::vector<Interval> mergeIntervals(std::vector<Interval> intervals) {
    if (intervals.empty()) {
        return intervals;
    }

    std::sort(intervals.begin(), intervals.end(), [](const Interval &lhs, const Interval &rhs) {
        if (lhs.start != rhs.start) {
            return lhs.start < rhs.start;
        }
        return lhs.end < rhs.end;
    });

    std::vector<Interval> merged;
    merged.push_back(intervals[0]);
    for (size_t i = 1; i < intervals.size(); ++i) {
        if (intervals[i].start <= merged.back().end + 1) {
            merged.back().end = std::max(merged.back().end, intervals[i].end);
        } else {
            merged.push_back(intervals[i]);
        }
    }
    return merged;
}

static std::string intervalsToString(const std::vector<Interval> &intervals) {
    std::string out;
    for (size_t i = 0; i < intervals.size(); ++i) {
        if (i > 0) {
            out.append(",");
        }
        out.append(SSTR(intervals[i].start + 1));
        out.append(":");
        out.append(SSTR(intervals[i].end + 1));
    }
    return out;
}

static unsigned int intervalCoverage(const std::vector<Interval> &intervals) {
    unsigned int covered = 0;
    for (size_t i = 0; i < intervals.size(); ++i) {
        covered += static_cast<unsigned int>(intervals[i].end - intervals[i].start + 1);
    }
    return covered;
}

static std::vector<TargetStats> collectTargetStats(const ReclassTaxContext &ctx) {
    std::unordered_map<unsigned int, TargetStats> statsByTarget;

    for (MappingTable::const_iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            const ReclassTaxEntry &entry = it->second[j];
            TargetStats &stats = statsByTarget[entry.result.dbKey];
            stats.key = entry.result.dbKey;
            stats.targetLength = entry.result.dbLen;
            stats.abundance = entry.abundance;
            stats.entropyValue = entry.entropyValue;
            stats.entropyPenalty = entry.entropyPenalty;
            stats.dropped = false;
            addInterval(stats.intervals, entry.result.dbStartPos, entry.result.dbEndPos);
        }
    }

    std::vector<TargetStats> out;
    out.reserve(statsByTarget.size());
    for (std::unordered_map<unsigned int, TargetStats>::iterator it = statsByTarget.begin(); it != statsByTarget.end(); ++it) {
        it->second.intervals = mergeIntervals(it->second.intervals);
        out.push_back(it->second);
    }

    std::sort(out.begin(), out.end(), [](const TargetStats &lhs, const TargetStats &rhs) {
        if (lhs.abundance != rhs.abundance) {
            return lhs.abundance > rhs.abundance;
        }
        return lhs.key < rhs.key;
    });
    return out;
}

static bool largestJumpCutoff(std::vector<double> values, double &cutoff) {
    cutoff = 0.0;
    if (values.size() < 4) {
        return false;
    }

    std::sort(values.begin(), values.end());
    double bestGap = 0.0;
    size_t bestIdx = 0;
    for (size_t i = 0; i + 1 < values.size(); ++i) {
        const double gap = values[i + 1] - values[i];
        if (gap > bestGap) {
            bestGap = gap;
            bestIdx = i;
        }
    }

    if (bestGap <= EPS) {
        return false;
    }

    cutoff = 0.5 * (values[bestIdx] + values[bestIdx + 1]);
    return true;
}

static std::unordered_set<unsigned int> selectDroppedTargets(const std::vector<TargetStats> &stats,
                                                             double &abundanceCutoff,
                                                             double &entropyCutoff) {
    std::unordered_set<unsigned int> dropped;
    if (stats.empty()) {
        abundanceCutoff = 0.0;
        entropyCutoff = 0.0;
        return dropped;
    }
    if (stats.size() < 4) {
        abundanceCutoff = 0.0;
        entropyCutoff = 0.0;
        return dropped;
    }

    std::vector<double> abundances;
    std::vector<double> entropies;
    abundances.reserve(stats.size());
    entropies.reserve(stats.size());
    for (size_t i = 0; i < stats.size(); ++i) {
        abundances.push_back(stats[i].abundance);
        entropies.push_back(stats[i].entropyValue);
    }

    const bool hasAbundanceCutoff = largestJumpCutoff(abundances, abundanceCutoff);
    const bool hasEntropyCutoff = largestJumpCutoff(entropies, entropyCutoff);
    if (hasAbundanceCutoff == false || hasEntropyCutoff == false) {
        abundanceCutoff = 0.0;
        entropyCutoff = 0.0;
        return dropped;
    }

    for (size_t i = 0; i < stats.size(); ++i) {
        if (stats[i].abundance <= abundanceCutoff && stats[i].entropyValue >= entropyCutoff) {
            dropped.insert(stats[i].key);
        }
    }
    if (dropped.size() == stats.size()) {
        dropped.clear();
    }
    return dropped;
}

static void applyDroppedTargets(ReclassTaxContext &ctx,
                                const std::unordered_set<unsigned int> &dropped,
                                size_t totalTargets,
                                double abundanceCutoff,
                                double entropyCutoff) {
    if (dropped.empty()) {
        Debug(Debug::INFO) << "Reclassify-taxonomy target filter kept all targets. abundance cutoff="
                           << abundanceCutoff << " entropy cutoff=" << entropyCutoff << "\n";
        return;
    }

    for (MappingTable::iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end();) {
        std::vector<ReclassTaxEntry> &records = it->second;
        records.erase(std::remove_if(records.begin(), records.end(), [&dropped](const ReclassTaxEntry &entry) {
            return dropped.find(entry.result.dbKey) != dropped.end();
        }), records.end());

        if (records.empty()) {
            it = ctx.mappingTable.erase(it);
        } else {
            ++it;
        }
    }

    ctx.queryOrder.erase(std::remove_if(ctx.queryOrder.begin(), ctx.queryOrder.end(), [&ctx](unsigned int queryKey) {
        return ctx.mappingTable.find(queryKey) == ctx.mappingTable.end();
    }), ctx.queryOrder.end());

    for (std::unordered_set<unsigned int>::const_iterator it = dropped.begin(); it != dropped.end(); ++it) {
        ctx.targetSet.erase(*it);
    }

    const double removedPct = (totalTargets > 0)
                                ? (100.0 * static_cast<double>(dropped.size()) / static_cast<double>(totalTargets))
                                : 0.0;
    Debug(Debug::INFO) << "Reclassify-taxonomy dropped " << dropped.size()
                       << " of " << totalTargets
                       << " targets (" << removedPct << "%)"
                       << " using abundance <= " << abundanceCutoff
                       << " and entropy >= " << entropyCutoff << ".\n";
}

static void markDroppedTargets(std::vector<TargetStats> &stats, const std::unordered_set<unsigned int> &dropped) {
    for (size_t i = 0; i < stats.size(); ++i) {
        stats[i].dropped = (dropped.find(stats[i].key) != dropped.end());
    }
}

static void convertAbundanceToPercent(std::vector<TargetStats> &stats) {
    double total = 0.0;
    for (size_t i = 0; i < stats.size(); ++i) {
        total += stats[i].abundance;
    }

    if (total <= 0.0) {
        for (size_t i = 0; i < stats.size(); ++i) {
            stats[i].abundance = 0.0;
        }
        return;
    }

    for (size_t i = 0; i < stats.size(); ++i) {
        stats[i].abundance = 100.0 * (stats[i].abundance / total);
    }
}

static void writeReclassifiedM8(const ReclassTaxContext &ctx,
                                DBReader<unsigned int> &queryHeaderReader,
                                DBReader<unsigned int> &targetHeaderReader,
                                const std::string &path) {
    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    char line[4096];

    for (size_t i = 0; i < ctx.queryOrder.size(); ++i) {
        const unsigned int queryKey = ctx.queryOrder[i];
        MappingTable::const_iterator recordsIt = ctx.mappingTable.find(queryKey);
        if (recordsIt == ctx.mappingTable.end()) {
            continue;
        }

        std::string queryId = identifierForKey(queryHeaderReader, queryKey, 0);
        std::vector<ReclassTaxEntry> records = recordsIt->second;
        SORT_SERIAL(records.begin(), records.end(), compareByPosterior);

        for (size_t j = 0; j < records.size(); ++j) {
            const Matcher::result_t &res = records[j].result;
            const std::string targetId = identifierForKey(targetHeaderReader, res.dbKey, 0);

            unsigned int alnLen = 0;
            unsigned int mismatchCount = 0;
            unsigned int gapOpenCount = 0;
            computeAlignmentCounts(res, alnLen, mismatchCount, gapOpenCount);

            const int written = snprintf(line, sizeof(line),
                                         "%s\t%s\t%1.3f\t%u\t%u\t%u\t%d\t%d\t%d\t%d\t%.2E\t%d\n",
                                         queryId.c_str(), targetId.c_str(), res.seqId, alnLen,
                                         mismatchCount, gapOpenCount,
                                         res.qStartPos + 1, res.qEndPos + 1,
                                         res.dbStartPos + 1, res.dbEndPos + 1,
                                         res.eval, res.score);
            if (written < 0 || static_cast<size_t>(written) >= sizeof(line)) {
                Debug(Debug::WARNING) << "Truncated M8 line for query " << queryKey << " and target " << res.dbKey << ".\n";
                continue;
            }
            fputs(line, handle);
        }
    }

    fclose(handle);
}

static void writeProteinStats(const std::vector<TargetStats> &stats,
                              DBReader<unsigned int> &targetHeaderReader,
                              MappingReader &mapping,
                              NcbiTaxonomy *taxonomy,
                              const std::string &path) {
    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    fputs("target_key\ttarget_id\tabundance_pct\tentropy\tDrop(y/n)\tmapping_parts\tmapped_length\ttarget_length\ttaxid\trank\ttaxname\ttaxlineage\n", handle);

    for (size_t i = 0; i < stats.size(); ++i) {
        const unsigned int key = stats[i].key;
        const std::string targetId = identifierForKey(targetHeaderReader, key, 0);
        const unsigned int taxId = mapping.lookup(key);
        const TaxonNode *node = (taxId != 0) ? taxonomy->taxonNode(taxId, false) : NULL;
        const std::string lineage = (node != NULL) ? taxonomy->taxLineage(node, true) : "unclassified";
        const std::string parts = intervalsToString(stats[i].intervals);
        const unsigned int mappedLength = intervalCoverage(stats[i].intervals);

        fprintf(handle, "%u\t%s\t%.12g\t%.12g\t%s\t%s\t%u\t%u\t%u\t%s\t%s\t%s\n",
                key,
                targetId.c_str(),
                stats[i].abundance,
                stats[i].entropyValue,
                stats[i].dropped ? "y" : "n",
                parts.c_str(),
                mappedLength,
                stats[i].targetLength,
                taxId,
                (node != NULL) ? taxonomy->getString(node->rankIdx) : "unclassified",
                (node != NULL) ? taxonomy->getString(node->nameIdx) : "unclassified",
                lineage.c_str());
    }

    fclose(handle);
}

static void writeTaxonomyStats(const std::vector<TargetStats> &stats,
                               MappingReader &mapping,
                               NcbiTaxonomy *taxonomy,
                               const std::string &path) {
    std::unordered_map<unsigned int, TaxonomyStats> aggregated;
    for (size_t i = 0; i < stats.size(); ++i) {
        const unsigned int taxId = mapping.lookup(stats[i].key);
        TaxonomyStats &entry = aggregated[taxId];
        entry.taxId = taxId;
        entry.abundance += stats[i].abundance;
        entry.entropySum += stats[i].entropyValue;
        entry.entropyPenaltySum += stats[i].entropyPenalty;
        entry.proteinCount += 1;
    }

    std::vector<TaxonomyStats> rows;
    rows.reserve(aggregated.size());
    for (std::unordered_map<unsigned int, TaxonomyStats>::const_iterator it = aggregated.begin(); it != aggregated.end(); ++it) {
        rows.push_back(it->second);
    }
    std::sort(rows.begin(), rows.end(), [](const TaxonomyStats &lhs, const TaxonomyStats &rhs) {
        if (lhs.abundance != rhs.abundance) {
            return lhs.abundance > rhs.abundance;
        }
        return lhs.taxId < rhs.taxId;
    });

    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    fputs("taxid\trank\ttaxname\ttaxlineage\tprotein_abundance_pct\tprotein_count\tmean_entropy\tmean_entropy_penalty\n", handle);

    for (size_t i = 0; i < rows.size(); ++i) {
        const TaxonNode *node = (rows[i].taxId != 0) ? taxonomy->taxonNode(rows[i].taxId, false) : NULL;
        const std::string lineage = (node != NULL) ? taxonomy->taxLineage(node, true) : "unclassified";
        const double denom = (rows[i].proteinCount > 0) ? static_cast<double>(rows[i].proteinCount) : 1.0;
        fprintf(handle, "%u\t%s\t%s\t%s\t%.12g\t%zu\t%.12g\t%.12g\n",
                rows[i].taxId,
                (node != NULL) ? taxonomy->getString(node->rankIdx) : "unclassified",
                (node != NULL) ? taxonomy->getString(node->nameIdx) : "unclassified",
                lineage.c_str(),
                rows[i].abundance,
                rows[i].proteinCount,
                rows[i].entropySum / denom,
                rows[i].entropyPenaltySum / denom);
    }

    fclose(handle);
}
}

int reclassifytaxonomy(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> queryHeaderReader((par.db1 + "_h").c_str(), (par.db1 + "_h.index").c_str(), par.threads,
                                             DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    queryHeaderReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> targetHeaderReader((par.db2 + "_h").c_str(), (par.db2 + "_h.index").c_str(), par.threads,
                                              DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    targetHeaderReader.open(DBReader<unsigned int>::NOSORT);

    NcbiTaxonomy *taxonomy = NcbiTaxonomy::openTaxonomy(par.db2);
    MappingReader mapping(par.db2);

    ReclassTaxContext ctx;
    loadAlignmentDb(reader, ctx);
    Debug(Debug::INFO) << "Loaded " << ctx.queryCount << " queries with hits and " << ctx.targetSet.size() << " unique targets.\n";

    squarem(ctx,
            par.reclassifyLambda,
            par.reclassifyMaxIterations,
            par.reclassifyTolerance,
            par.reclassifyAlpha,
            par.reclassifyGamma);

    const std::string outDir = par.db4;
    const std::string m8Path = outDir + "/new_alignment_result.m8";
    const std::string proteinPath = outDir + "/protein_abundance.tsv";
    const std::string taxonomyPath = outDir + "/taxonomy_abundance.tsv";

    std::vector<TargetStats> allTargetStats = collectTargetStats(ctx);
    double abundanceCutoff = 0.0;
    double entropyCutoff = 0.0;
    const std::unordered_set<unsigned int> dropped = selectDroppedTargets(allTargetStats,
                                                                          abundanceCutoff,
                                                                          entropyCutoff);
    markDroppedTargets(allTargetStats, dropped);
    convertAbundanceToPercent(allTargetStats);

    std::vector<TargetStats> targetStats = allTargetStats;
    targetStats.erase(std::remove_if(targetStats.begin(), targetStats.end(), [](const TargetStats &entry) {
        return entry.dropped;
    }), targetStats.end());
    applyDroppedTargets(ctx, dropped, allTargetStats.size(), abundanceCutoff, entropyCutoff);
    convertAbundanceToPercent(targetStats);

    writeReclassifiedM8(ctx, queryHeaderReader, targetHeaderReader, m8Path);
    writeProteinStats(allTargetStats, targetHeaderReader, mapping, taxonomy, proteinPath);
    if (par.reclassifyTaxonomy == 1) {
        writeTaxonomyStats(targetStats, mapping, taxonomy, taxonomyPath);
    }

    delete taxonomy;
    targetHeaderReader.close();
    queryHeaderReader.close();
    reader.close();
    return EXIT_SUCCESS;
}
