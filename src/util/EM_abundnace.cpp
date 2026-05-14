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
#include <cctype>
#include <cstdio>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef OPENMP
#include <omp.h>
#endif

namespace {
struct ReclassTaxEntry {
    Matcher::result_t result;
    double abundance;
    double posterior;
    double coverageConfidence;
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
    double coverageConfidence;
    bool dropped;
    std::vector<Interval> intervals;
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

static const double EPS = 1e-12;
static const size_t MIN_FILTER_TARGETS = 20;
static const size_t MIN_TAIL_TARGETS = 2;

static double clamp01(double value);

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

            if (columns == Matcher::ALN_RES_WITH_BT_COL_CNT || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT
                || columns == Matcher::ALN_RES_WITH_BT_COL_CNT + 1 || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT + 1) {
                ctx.hasBacktrace = true;
            }
            if (columns == Matcher::ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT
                || columns == Matcher::ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT + 1 || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT + 1) {
                ctx.hasOrfPosition = true;
            }

            Matcher::result_t result = Matcher::parseAlignmentRecord(data, true);

            // Always use the 3rd alignment column (seqId) as posterior.
            // Do not consume optional trailing columns as posterior.
            double posterior = static_cast<double>(result.seqId);

            records.push_back(ReclassTaxEntry{result, 0.0, posterior, 0.0});
            ctx.targetSet.insert(result.dbKey);
            data = Util::skipLine(data);
        }
    }

    ctx.queryCount = ctx.mappingTable.size();
}

struct TargetHitRef {
    const ReclassTaxEntry *entry;
    double score;
    double weight;
};

static void initCoverageConfidence(MappingTable &mappingTable,
                                   const std::unordered_set<unsigned int> &targetSet,
                                   int threads) {
    (void)threads;
    std::unordered_map<unsigned int, int> targetMin;
    std::unordered_map<unsigned int, int> targetMax;
    std::unordered_map<unsigned int, std::vector<TargetHitRef> > hitsByTarget;

    for (std::unordered_set<unsigned int>::const_iterator it = targetSet.begin(); it != targetSet.end(); ++it) {
        targetMin[*it] = std::numeric_limits<int>::max();
        targetMax[*it] = std::numeric_limits<int>::min();
        hitsByTarget.emplace(*it, std::vector<TargetHitRef>());
    }

    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        double scoreSum = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            scoreSum += static_cast<double>(it->second[j].result.score);
        }
        for (size_t j = 0; j < it->second.size(); ++j) {
            const unsigned int target = it->second[j].result.dbKey;
            if (it->second[j].result.dbStartPos < targetMin[target]) {
                targetMin[target] = it->second[j].result.dbStartPos;
            }
            if (it->second[j].result.dbEndPos > targetMax[target]) {
                targetMax[target] = it->second[j].result.dbEndPos;
            }
            const double score = static_cast<double>(it->second[j].result.score);
            const double weight = (scoreSum > 0.0) ? (score / scoreSum) : 0.0;
            hitsByTarget[target].push_back(TargetHitRef{&it->second[j], score, weight});
        }
    }

    std::unordered_map<unsigned int, double> coverageFraction;
    coverageFraction.reserve(targetSet.size());
    const std::vector<unsigned int> targetList = targetListFromSet(targetSet);
    std::vector<double> coverageFractionByIndex(targetList.size(), 0.0);

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (size_t i = 0; i < targetList.size(); ++i) {
        const unsigned int target = targetList[i];
        const int startPos = targetMin[target];
        const int endPos = targetMax[target];
        const int len = (endPos >= startPos) ? (endPos - startPos + 1) : 1;
        std::vector<double> cov(static_cast<size_t>(len), 0.0);
        std::vector<double> covConf(static_cast<size_t>(len), 0.0);

        std::unordered_map<unsigned int, std::vector<TargetHitRef> >::const_iterator hitIt = hitsByTarget.find(target);
        if (hitIt != hitsByTarget.end()) {
            const std::vector<TargetHitRef> &hits = hitIt->second;
            for (size_t h = 0; h < hits.size(); ++h) {
                const Matcher::result_t &result = hits[h].entry->result;
                const int targetLen = result.dbEndPos - result.dbStartPos + 1;
                if (targetLen <= 0) {
                    continue;
                }

                const double mq = hits[h].score / static_cast<double>(targetLen);
                const int start = std::max(0, result.dbStartPos - startPos);
                const int end = std::min(len - 1, result.dbEndPos - startPos);
                for (int pos = start; pos <= end; ++pos) {
                    cov[static_cast<size_t>(pos)] += mq;
                    covConf[static_cast<size_t>(pos)] += hits[h].weight;
                }
            }
        }

        double covered = 0.0;
        double squaredCovered = 0.0;
        for (size_t pos = 0; pos < covConf.size(); ++pos) {
            const double clipped = std::min(1.0, covConf[pos]);
            covered += clipped;
            squaredCovered += clipped * clipped;
        }
        const double fraction = covered / static_cast<double>(len);
        const double hhi = (covered > 0.0) ? (squaredCovered / (covered * covered)) : 1.0;
        const double concentrationPenalty = 1.0 - hhi;
        coverageFractionByIndex[i] = clamp01(fraction * concentrationPenalty);
    }

    for (size_t i = 0; i < targetList.size(); ++i) {
        coverageFraction[targetList[i]] = coverageFractionByIndex[i];
    }

    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            const unsigned int target = it->second[j].result.dbKey;
            std::unordered_map<unsigned int, double>::const_iterator cf = coverageFraction.find(target);
            it->second[j].coverageConfidence = (cf != coverageFraction.end()) ? cf->second : 0.0;
        }
    }
}

static void computeAbundanceFromPosterior(MappingTable &mappingTable,
                                          const std::unordered_set<unsigned int> &targetSet,
                                          size_t queryCount) {
    std::unordered_map<unsigned int, double> abundance;
    abundance.reserve(targetSet.size());
    for (std::unordered_set<unsigned int>::const_iterator it = targetSet.begin(); it != targetSet.end(); ++it) {
        abundance[*it] = 0.0;
    }

    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            abundance[it->second[j].result.dbKey] += it->second[j].posterior;
        }
    }

    if (queryCount > 0) {
        const double denom = static_cast<double>(queryCount);
        for (std::unordered_map<unsigned int, double>::iterator it = abundance.begin(); it != abundance.end(); ++it) {
            it->second /= denom;
        }
    }

    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            it->second[j].abundance = abundance[it->second[j].result.dbKey];
        }
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
            stats.coverageConfidence = entry.coverageConfidence;
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

static double clamp01(double value) {
    return std::max(0.0, std::min(1.0, value));
}

static bool largestJumpCutoff(std::vector<double> values,
                              bool useLowTail,
                              double maxTailFraction,
                              double &cutoff,
                              size_t &tailCount) {
    cutoff = 0.0;
    tailCount = 0;
    if (values.size() < MIN_FILTER_TARGETS) {
        return false;
    }

    std::sort(values.begin(), values.end());
    maxTailFraction = clamp01(maxTailFraction);
    const double totalMass = std::accumulate(values.begin(), values.end(), 0.0);
    if (totalMass <= EPS || maxTailFraction <= 0.0) {
        return false;
    }
    const double maxTailMass = maxTailFraction * totalMass;

    double bestGap = 0.0;
    size_t bestIdx = 0;
    double lowTailMass = 0.0;
    for (size_t i = 0; i + 1 < values.size(); ++i) {
        const size_t lowTailCount = i + 1;
        const size_t highTailCount = values.size() - lowTailCount;
        const size_t candidateTailCount = useLowTail ? lowTailCount : highTailCount;
        lowTailMass += values[i];
        const double highTailMass = totalMass - lowTailMass;
        const double candidateTailMass = useLowTail ? lowTailMass : highTailMass;
        if (candidateTailCount < MIN_TAIL_TARGETS || candidateTailMass > (maxTailMass + EPS)) {
            continue;
        }

        const double gap = values[i + 1] - values[i];
        if (gap > bestGap) {
            bestGap = gap;
            bestIdx = i;
            tailCount = candidateTailCount;
        }
    }

    if (bestGap <= EPS) {
        return false;
    }

    cutoff = 0.5 * (values[bestIdx] + values[bestIdx + 1]);
    return true;
}

static bool tailQuantileCutoff(std::vector<double> values,
                               bool useLowTail,
                               double maxTailFraction,
                               double &cutoff,
                               size_t &tailCount) {
    cutoff = 0.0;
    tailCount = 0;
    if (values.size() < MIN_FILTER_TARGETS) {
        return false;
    }

    std::sort(values.begin(), values.end());
    maxTailFraction = clamp01(maxTailFraction);
    const double totalMass = std::accumulate(values.begin(), values.end(), 0.0);
    if (totalMass <= EPS || maxTailFraction <= 0.0) {
        return false;
    }
    const double maxTailMass = maxTailFraction * totalMass;

    double accumulatedMass = 0.0;
    size_t maxTailCount = 0;
    for (size_t i = 0; i < values.size(); ++i) {
        const double candidate = accumulatedMass + values[i];
        if (candidate > (maxTailMass + EPS)) {
            break;
        }
        accumulatedMass = candidate;
        ++maxTailCount;
    }
    if (maxTailCount < MIN_TAIL_TARGETS || maxTailCount >= values.size()) {
        return false;
    }

    tailCount = maxTailCount;
    if (useLowTail) {
        cutoff = values[tailCount - 1];
    } else {
        cutoff = values[values.size() - tailCount];
    }
    return true;
}

static std::unordered_set<unsigned int> selectTailTargets(const std::vector<TargetStats> &stats,
                                                          bool useLowTail,
                                                          size_t tailCount,
                                                          double maxTailFraction) {
    std::vector<const TargetStats *> ordered;
    ordered.reserve(stats.size());
    for (size_t i = 0; i < stats.size(); ++i) {
        ordered.push_back(&stats[i]);
    }

    std::sort(ordered.begin(), ordered.end(), [useLowTail](const TargetStats *lhs, const TargetStats *rhs) {
        const double lhsValue = useLowTail ? lhs->abundance : lhs->coverageConfidence;
        const double rhsValue = useLowTail ? rhs->abundance : rhs->coverageConfidence;
        if (lhsValue != rhsValue) {
            return useLowTail ? (lhsValue < rhsValue) : (lhsValue > rhsValue);
        }
        return lhs->key < rhs->key;
    });

    double totalMass = 0.0;
    for (size_t i = 0; i < ordered.size(); ++i) {
        const double value = useLowTail ? ordered[i]->abundance : ordered[i]->coverageConfidence;
        totalMass += value;
    }
    const double maxTailMass = clamp01(maxTailFraction) * totalMass;

    std::unordered_set<unsigned int> selected;
    const size_t limit = std::min(tailCount, ordered.size());
    double selectedMass = 0.0;
    selected.reserve(limit);
    for (size_t i = 0; i < limit; ++i) {
        const double value = useLowTail ? ordered[i]->abundance : ordered[i]->coverageConfidence;
        if (selected.size() >= MIN_TAIL_TARGETS && (selectedMass + value) > (maxTailMass + EPS)) {
            break;
        }
        selectedMass += value;
        selected.insert(ordered[i]->key);
    }
    return selected;
}

static std::unordered_set<unsigned int> selectDroppedTargets(const std::vector<TargetStats> &stats,
                                                             double maxDropPercentage,
                                                             double &abundanceCutoff) {
    std::unordered_set<unsigned int> dropped;
    if (stats.empty()) {
        abundanceCutoff = 0.0;
        return dropped;
    }
    if (stats.size() < MIN_FILTER_TARGETS) {
        abundanceCutoff = 0.0;
        return dropped;
    }

    std::vector<double> abundances;
    abundances.reserve(stats.size());
    for (size_t i = 0; i < stats.size(); ++i) {
        abundances.push_back(stats[i].abundance);
    }

    const double maxTailFraction = clamp01(maxDropPercentage / 100.0);
    size_t abundanceTailCount = 0;
    bool hasAbundanceCutoff = largestJumpCutoff(abundances, true, maxTailFraction, abundanceCutoff, abundanceTailCount);
    if (hasAbundanceCutoff == false) {
        hasAbundanceCutoff = tailQuantileCutoff(abundances, true, maxTailFraction, abundanceCutoff, abundanceTailCount);
    }
    if (hasAbundanceCutoff == false) {
        abundanceCutoff = 0.0;
        return dropped;
    }

    const std::unordered_set<unsigned int> lowAbundanceTargets = selectTailTargets(stats, true, abundanceTailCount, maxTailFraction);
    for (std::unordered_set<unsigned int>::const_iterator it = lowAbundanceTargets.begin(); it != lowAbundanceTargets.end(); ++it) {
        dropped.insert(*it);
    }
    if (dropped.size() == stats.size()) {
        dropped.clear();
    }
    return dropped;
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

static void writeProteinStats(const std::vector<TargetStats> &stats,
                              DBReader<unsigned int> &targetHeaderReader,
                              const std::string &path) {
    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    fputs("target_key\ttarget_id\tabundance_pct\tcoverage_confidence\tDrop(y/n)\tmapped_length\ttarget_length\n", handle);

    for (size_t i = 0; i < stats.size(); ++i) {
        const unsigned int key = stats[i].key;
        const std::string targetId = identifierForKey(targetHeaderReader, key, 0);
        const unsigned int mappedLength = intervalCoverage(stats[i].intervals);

        fprintf(handle, "%u\t%s\t%.12g\t%.12g\t%s\t%u\t%u\n",
                key,
                targetId.c_str(),
                stats[i].abundance,
                stats[i].coverageConfidence,
                stats[i].dropped ? "y" : "n",
                mappedLength,
                stats[i].targetLength);
    }

    fclose(handle);
}

static void writeKrakenReport(const std::vector<TargetStats> &stats,
                              MappingReader &mapping,
                              NcbiTaxonomy *taxonomy,
                              size_t queryCount,
                              const std::string &path) {
    std::unordered_map<TaxID, unsigned int> directCounts;
    directCounts.reserve(stats.size());

    for (size_t i = 0; i < stats.size(); ++i) {
        const TaxID taxId = mapping.lookup(stats[i].key);
        if (taxId == 0) {
            continue;
        }
        const double expectedReads = stats[i].abundance * static_cast<double>(queryCount) / 100.0;
        directCounts[taxId] += static_cast<unsigned int>(std::floor(expectedReads + 0.5));
    }

    const std::unordered_map<TaxID, std::vector<TaxID> > parentToChildren = taxonomy->getParentToChildren();
    const std::unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(directCounts, parentToChildren);

    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    const double totalReads = (queryCount > 0) ? static_cast<double>(queryCount) : 1.0;

    std::vector<TaxID> stack;
    std::vector<int> depthStack;
    stack.push_back(1);
    depthStack.push_back(0);

    while (!stack.empty()) {
        TaxID taxId = stack.back();
        stack.pop_back();
        int depth = depthStack.back();
        depthStack.pop_back();

        unsigned int cladeCount = 0;
        unsigned int directCount = 0;
        std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxId);
        if (it != cladeCounts.end()) {
            cladeCount = it->second.cladeCount;
            directCount = it->second.taxCount;
        }

        if (cladeCount > 0) {
            const TaxonNode *node = taxonomy->taxonNode(taxId, false);
            const char *rankStr = (node != NULL) ? taxonomy->getString(node->rankIdx) : NULL;
            char rankCode = '-';
            if (rankStr != NULL) {
                std::map<std::string, char>::const_iterator rankIt = NcbiShortRanks.find(std::string(rankStr));
                if (rankIt != NcbiShortRanks.end()) {
                    rankCode = rankIt->second;
                }
            }
            const char *name = (node != NULL) ? taxonomy->getString(node->nameIdx) : "unclassified";
            const double pct = 100.0 * static_cast<double>(cladeCount) / totalReads;

            for (int i = 0; i < depth; ++i) {
                fputs("  ", handle);
            }
            fprintf(handle, "%.4f\t%u\t%u\t%c\t%u\t%s\n",
                    pct, cladeCount, directCount, rankCode, static_cast<unsigned int>(taxId), name);
        }

        std::unordered_map<TaxID, TaxonCounts>::const_iterator childIt = cladeCounts.find(taxId);
        if (childIt != cladeCounts.end()) {
            const std::vector<TaxID> &children = childIt->second.children;
            for (size_t i = 0; i < children.size(); ++i) {
                stack.push_back(children[i]);
                depthStack.push_back(depth + 1);
            }
        }
    }

    fclose(handle);
}

static void writeBrackenReport(const std::vector<TargetStats> &stats,
                               MappingReader &mapping,
                               NcbiTaxonomy *taxonomy,
                               size_t queryCount,
                               const std::string &path) {
    std::unordered_map<TaxID, unsigned int> directCounts;
    directCounts.reserve(stats.size());

    for (size_t i = 0; i < stats.size(); ++i) {
        const TaxID taxId = mapping.lookup(stats[i].key);
        if (taxId == 0) {
            continue;
        }
        const double expectedReads = stats[i].abundance * static_cast<double>(queryCount) / 100.0;
        directCounts[taxId] += static_cast<unsigned int>(std::floor(expectedReads + 0.5));
    }

    const std::unordered_map<TaxID, std::vector<TaxID> > parentToChildren = taxonomy->getParentToChildren();
    const std::unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(directCounts, parentToChildren);

    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    fputs("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n", handle);

    const double totalReads = (queryCount > 0) ? static_cast<double>(queryCount) : 1.0;
    for (std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.begin(); it != cladeCounts.end(); ++it) {
        const TaxID taxId = it->first;
        const TaxonCounts &counts = it->second;
        if (counts.cladeCount == 0) {
            continue;
        }

        const TaxonNode *node = taxonomy->taxonNode(taxId, false);
        const char *rankStr = (node != NULL) ? taxonomy->getString(node->rankIdx) : "-";
        const char *name = (node != NULL) ? taxonomy->getString(node->nameIdx) : "unclassified";
        const unsigned int cladeCount = counts.cladeCount;
        const unsigned int directCount = counts.taxCount;
        const unsigned int addedReads = (cladeCount >= directCount) ? (cladeCount - directCount) : 0;
        const double fraction = static_cast<double>(cladeCount) / totalReads;

        fprintf(handle, "%s\t%u\t%s\t%u\t%u\t%u\t%.12g\n",
                name,
                static_cast<unsigned int>(taxId),
                rankStr,
                directCount,
                addedReads,
                cladeCount,
                fraction);
    }

    fclose(handle);
}
}

int emabundance(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> targetHeaderReader((par.db2 + "_h").c_str(), (par.db2 + "_h.index").c_str(), par.threads,
                                              DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    targetHeaderReader.open(DBReader<unsigned int>::NOSORT);

    const bool withTaxonomy = (par.reclassifyTaxonomy == 1);
    NcbiTaxonomy *taxonomy = NULL;
    MappingReader *mapping = NULL;
    if (withTaxonomy) {
        taxonomy = NcbiTaxonomy::openTaxonomy(par.db2);
        mapping = new MappingReader(par.db2);
    }

    ReclassTaxContext ctx;
    loadAlignmentDb(reader, ctx);
    Debug(Debug::INFO) << "Loaded " << ctx.queryCount << " queries with hits and " << ctx.targetSet.size() << " unique targets.\n";

    initCoverageConfidence(ctx.mappingTable, ctx.targetSet, par.threads);
    computeAbundanceFromPosterior(ctx.mappingTable, ctx.targetSet, ctx.queryCount);

    std::vector<TargetStats> allTargetStats = collectTargetStats(ctx);
    double abundanceCutoff = 0.0;
    const std::unordered_set<unsigned int> dropped = selectDroppedTargets(allTargetStats,
                                                                          par.reclassifyMaxDropPercentage,
                                                                          abundanceCutoff);
    markDroppedTargets(allTargetStats, dropped);
    convertAbundanceToPercent(allTargetStats);

    if (withTaxonomy) {
        std::vector<TargetStats> targetStats = allTargetStats;
        targetStats.erase(std::remove_if(targetStats.begin(), targetStats.end(), [](const TargetStats &entry) {
            return entry.dropped;
        }), targetStats.end());
        writeKrakenReport(targetStats, *mapping, taxonomy, ctx.queryCount, par.db4);
        writeBrackenReport(targetStats, *mapping, taxonomy, ctx.queryCount, par.db4 + ".bracken");
    } else {
        writeProteinStats(allTargetStats, targetHeaderReader, par.db4);
    }

    delete mapping;
    delete taxonomy;
    targetHeaderReader.close();
    reader.close();
    return EXIT_SUCCESS;
}
