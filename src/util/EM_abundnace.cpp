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
    double posterior;   // EM posterior gamma_qt (read from the reclassified alignment's column 3)
};

typedef std::unordered_map<unsigned int, std::vector<ReclassTaxEntry> > MappingTable;

struct TargetStats {
    unsigned int key;
    unsigned int targetLength;
    double abundance;
    bool dropped;
    // Posterior-weighted coverage, using the EM assignment probabilities gamma_qt.
    //   coverageDepth  = mean expected read depth over the target = (1/L) * sum_hits gamma * (#covered M cols)
    //   coveredLength  = expected number of covered residues = sum_p [1 - prod_hits(1 - gamma)]
    //   breadth        = coveredLength / L  (derived at write time)
    double coverageDepth;
    double coveredLength;
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

            records.push_back(ReclassTaxEntry{result, 0.0, posterior});
            ctx.targetSet.insert(result.dbKey);
            data = Util::skipLine(data);
        }
    }

    ctx.queryCount = ctx.mappingTable.size();
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

// Mark the target residues covered by a single hit onto `covered` (indexed by target position).
// The backtrace is a run-length-encoded CIGAR (e.g. "153M2D10M"; a missing count means 1).
// Only 'M' columns place a read residue on a target position, so only they mark coverage.
// 'M' and 'D' both advance the target position (a residue is consumed on the target);
// 'I' consumes the query only. Positions outside [0, targetLength) are ignored, so a hit
// can never contribute coverage beyond the target itself.
// Accumulate one hit's soft (posterior-weighted) coverage. For each target position covered by an
// 'M' column (walking the CIGAR from dbStartPos, clipped to [0, L)):
//   - multiply notCov[p] by (1 - gamma), so 1 - notCov[p] becomes P[position covered by >=1 read]
//   - count the position toward the hit's covered-M length (returned, for the depth sum)
// gamma is the EM posterior gamma_qt of this read->target assignment.
static double accumulateBacktrace(std::vector<double> &notCov,
                                  const std::string &compressedBacktrace,
                                  int dbStartPos,
                                  double gamma) {
    const int targetLength = static_cast<int>(notCov.size());
    const double keep = 1.0 - gamma;
    int targetPos = dbStartPos;
    size_t count = 0;
    double coveredM = 0.0;
    for (size_t i = 0; i < compressedBacktrace.size(); ++i) {
        const char c = compressedBacktrace[i];
        if (c >= '0' && c <= '9') {
            count = count * 10 + static_cast<size_t>(c - '0');
            continue;
        }
        const int n = (count == 0) ? 1 : static_cast<int>(count);
        count = 0;
        if (c == 'M') {
            for (int k = 0; k < n; ++k, ++targetPos) {
                if (targetPos >= 0 && targetPos < targetLength) {
                    notCov[static_cast<size_t>(targetPos)] *= keep;
                    coveredM += 1.0;
                }
            }
        } else if (c == 'D') {
            targetPos += n;
        }
        // 'I' consumes the query only: the target position and coverage are unchanged.
    }
    return coveredM;
}

// Fallback (no backtrace): treat the whole aligned span [dbStartPos, dbEndPos] as covered,
// clipped to [0, L). Same posterior weighting as accumulateBacktrace.
static double accumulateSpan(std::vector<double> &notCov, int dbStartPos, int dbEndPos, double gamma) {
    const int targetLength = static_cast<int>(notCov.size());
    const double keep = 1.0 - gamma;
    const int start = std::max(0, std::min(dbStartPos, dbEndPos));
    const int end = std::min(targetLength - 1, std::max(dbStartPos, dbEndPos));
    double coveredM = 0.0;
    for (int pos = start; pos <= end; ++pos) {
        notCov[static_cast<size_t>(pos)] *= keep;
        coveredM += 1.0;
    }
    return coveredM;
}

static std::vector<TargetStats> collectTargetStats(const ReclassTaxContext &ctx, int threads) {
    (void)threads;
    // Group the hits by target and record the per-target metadata in one pass. We keep a pointer to
    // the whole entry (not just result) because the coverage weighting needs its posterior.
    std::unordered_map<unsigned int, TargetStats> statsByTarget;
    std::unordered_map<unsigned int, std::vector<const ReclassTaxEntry *> > hitsByTarget;

    for (MappingTable::const_iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            const ReclassTaxEntry &entry = it->second[j];
            TargetStats &stats = statsByTarget[entry.result.dbKey];
            stats.key = entry.result.dbKey;
            stats.targetLength = entry.result.dbLen;
            stats.abundance = entry.abundance;
            stats.dropped = false;
            stats.coverageDepth = 0.0;
            stats.coveredLength = 0.0;
            hitsByTarget[entry.result.dbKey].push_back(&entry);
        }
    }

    std::vector<TargetStats> out;
    out.reserve(statsByTarget.size());
    for (std::unordered_map<unsigned int, TargetStats>::const_iterator it = statsByTarget.begin(); it != statsByTarget.end(); ++it) {
        out.push_back(it->second);
    }

    // Posterior-weighted coverage per target, on a map sized to the target length. Every position is
    // confined to [0, targetLength), so coveredLength <= targetLength by construction (no clamping).
    Debug::Progress covProgress(out.size());
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (size_t i = 0; i < out.size(); ++i) {
        covProgress.updateProgress();
        TargetStats &stats = out[i];
        if (stats.targetLength == 0) {
            stats.coverageDepth = 0.0;
            stats.coveredLength = 0.0;
            continue;
        }
        // notCov[p] = prod over hits covering p of (1 - gamma); starts at 1 (uncovered).
        std::vector<double> notCov(stats.targetLength, 1.0);
        double depthSum = 0.0;   // sum over hits of gamma * (#covered M columns)
        const std::vector<const ReclassTaxEntry *> &hits = hitsByTarget[stats.key];
        for (size_t h = 0; h < hits.size(); ++h) {
            const ReclassTaxEntry *entry = hits[h];
            const double gamma = entry->posterior;
            double coveredM;
            if (entry->result.backtrace.empty() == false) {
                coveredM = accumulateBacktrace(notCov, entry->result.backtrace, entry->result.dbStartPos, gamma);
            } else {
                coveredM = accumulateSpan(notCov, entry->result.dbStartPos, entry->result.dbEndPos, gamma);
            }
            depthSum += gamma * coveredM;
        }
        // Expected covered residues = sum_p P[covered] = sum_p (1 - notCov[p]).
        double coveredLength = 0.0;
        for (size_t pos = 0; pos < notCov.size(); ++pos) {
            coveredLength += (1.0 - notCov[pos]);
        }
        stats.coveredLength = coveredLength;
        stats.coverageDepth = depthSum / static_cast<double>(stats.targetLength);
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

    // Dropping ranks by abundance only (low tail); the function is never called on any other key.
    std::sort(ordered.begin(), ordered.end(), [useLowTail](const TargetStats *lhs, const TargetStats *rhs) {
        const double lhsValue = lhs->abundance;
        const double rhsValue = rhs->abundance;
        if (lhsValue != rhsValue) {
            return useLowTail ? (lhsValue < rhsValue) : (lhsValue > rhsValue);
        }
        return lhs->key < rhs->key;
    });

    double totalMass = 0.0;
    for (size_t i = 0; i < ordered.size(); ++i) {
        totalMass += ordered[i]->abundance;
    }
    const double maxTailMass = clamp01(maxTailFraction) * totalMass;

    std::unordered_set<unsigned int> selected;
    const size_t limit = std::min(tailCount, ordered.size());
    double selectedMass = 0.0;
    selected.reserve(limit);
    for (size_t i = 0; i < limit; ++i) {
        const double value = ordered[i]->abundance;
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
    fputs("target_key\ttarget_id\tabundance_pct\tcoverage_depth\tbreadth_of_coverage\tDrop(y/n)\tcovered_length\ttarget_length\n", handle);

    for (size_t i = 0; i < stats.size(); ++i) {
        const unsigned int key = stats[i].key;
        const std::string targetId = identifierForKey(targetHeaderReader, key, 0);
        // Posterior-weighted coverage (see collectTargetStats). coveredLength is the expected number
        // of covered residues and is bounded by targetLength, so breadth is naturally in [0, 1].
        const double coveredLength = stats[i].coveredLength;
        const double breadthOfCoverage = (stats[i].targetLength > 0)
                ? coveredLength / static_cast<double>(stats[i].targetLength)
                : 0.0;

        fprintf(handle, "%u\t%s\t%.12g\t%.12g\t%.12g\t%s\t%.12g\t%u\n",
                key,
                targetId.c_str(),
                stats[i].abundance,
                stats[i].coverageDepth,
                breadthOfCoverage,
                stats[i].dropped ? "y" : "n",
                coveredLength,
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

    Debug(Debug::INFO) << "Computing abundance from posterior...\n";
    computeAbundanceFromPosterior(ctx.mappingTable, ctx.targetSet, ctx.queryCount);

    Debug(Debug::INFO) << "Collecting per-target statistics...\n";
    std::vector<TargetStats> allTargetStats = collectTargetStats(ctx, par.threads);
    double abundanceCutoff = 0.0;
    const std::unordered_set<unsigned int> dropped = selectDroppedTargets(allTargetStats,
                                                                          par.reclassifyMaxDropPercentage,
                                                                          abundanceCutoff);
    markDroppedTargets(allTargetStats, dropped);
    convertAbundanceToPercent(allTargetStats);
    Debug(Debug::INFO) << "Writing output...\n";

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
