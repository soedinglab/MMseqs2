#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include "FileUtil.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"
#include "EMCoverage.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
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
    //   coverageDepth     = mean expected read depth over the target = (1/L) * sum_p depth_p
    //   coveredLength     = expected number of covered residues = sum_p [1 - prod_hits(1 - gamma)]
    //   breadth           = coveredLength / L  (derived at write time)
    //   coverageEvenness  = how uniformly that depth is spread over the target (EMCoverage::evenness)
    //
    // breadth and evenness answer different questions and neither replaces the other. breadth
    // saturates: 1 - prod(1 - gamma) -> 1 as soon as enough reads touch a position, so with
    // --max-seqs 300 on a deep sample even a spurious target reaches breadth ~= 1 and the column
    // stops discriminating. Evenness stays informative there, because it drops when the coverage
    // piles up on a small part of the target instead of spreading over it.
    double coverageDepth;
    double coveredLength;
    double coverageEvenness;
};

// No queryOrder / hasOrfPosition here: unlike `reclassify`, `abundance` writes no alignment DB,
// so nothing ever read them. The counters below are diagnostics only -- they are what tells you
// which coverage path actually ran (see reportInputDiagnostics).
struct ReclassTaxContext {
    MappingTable mappingTable;
    std::unordered_set<unsigned int> targetSet;
    size_t queryCount;
    size_t hitCount;
    size_t backtraceHitCount;   // hits whose coverage came from the CIGAR
    size_t reverseStrandHits;   // hits with dbStartPos > dbEndPos (reverse-strand target)

    ReclassTaxContext() : queryCount(0), hitCount(0), backtraceHitCount(0), reverseStrandHits(0) {}
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
        while (*data != '\0') {
            const size_t columns = Util::getWordsOfLine(data, entry, 255);
            if (columns < Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                Debug(Debug::ERROR) << "Invalid alignment result record in query " << queryKey << ".\n";
                EXIT(EXIT_FAILURE);
            }

            Matcher::result_t result = Matcher::parseAlignmentRecord(data, true);

            // Always use the 3rd alignment column (seqId) as posterior.
            // Do not consume optional trailing columns as posterior.
            // `mmseqs reclassify` writes that column with resultToBuffer(..., highPrecisionSeqId=true),
            // and parseAlignmentRecord reads it through fast_float, so the full float range arrives
            // intact. Clamp anyway: result_t::seqId is a float, so a round-tripped 1.0 can come back
            // as 1.0000001, and accumulateBacktrace/accumulateSpan would then multiply by a negative
            // (1 - gamma) and report a covered length above the target length.
            double posterior = clamp01(static_cast<double>(result.seqId));

            ++ctx.hitCount;
            if (result.backtrace.empty() == false) {
                ++ctx.backtraceHitCount;
            }
            if (result.dbStartPos > result.dbEndPos) {
                ++ctx.reverseStrandHits;
            }

            records.push_back(ReclassTaxEntry{result, 0.0, posterior});
            ctx.targetSet.insert(result.dbKey);
            data = Util::skipLine(data);
        }
    }

    ctx.queryCount = ctx.mappingTable.size();
}

// `abundance` trusts alignment column 3 to hold the EM posterior gamma_qt that `reclassify` wrote
// there. Nothing in the file format distinguishes that from an ordinary sequence identity, so an
// alignment DB that never went through `reclassify` produces a plausible-looking but meaningless
// abundance table. The one invariant that separates the two is that EM posteriors sum to 1 over
// each query's hits (computePosterior() normalises by construction, and denom > 0 always holds
// because phi >= exp(-60) and pi >= 1e-8), whereas sequence identities do not.
//
// Also report which coverage path the hits will take: without a backtrace the coverage falls back
// to the whole aligned span, which ignores indels and over-counts by the gap length. That happens
// silently when the upstream search ran without -a.
static const char *WARNING_RULE =
        "========================================================================\n";

// Both conditions below mean the numbers `abundance` is about to print are wrong, not merely
// imprecise, so they get a banner rather than a line that scrolls past in the middle of the
// progress output.
static void warningBanner(const std::string &headline, const std::string &detail) {
    Debug(Debug::WARNING) << "\n" << WARNING_RULE
                          << "  WARNING: " << headline << "\n"
                          << detail
                          << WARNING_RULE << "\n";
}

static void reportInputDiagnostics(const ReclassTaxContext &ctx) {
    Debug(Debug::INFO) << ctx.hitCount << " hits; coverage source: "
                       << ctx.backtraceHitCount << " from backtrace, "
                       << (ctx.hitCount - ctx.backtraceHitCount) << " from aligned span";
    if (ctx.reverseStrandHits > 0) {
        Debug(Debug::INFO) << "; " << ctx.reverseStrandHits << " reverse-strand";
    }
    Debug(Debug::INFO) << ".\n";
    if (ctx.backtraceHitCount < ctx.hitCount) {
        warningBanner("Some hits carry no backtrace.",
                      "           You should run 'mmseqs search' with the '-a' parameter.\n"
                      "\n"
                      "           " + SSTR(ctx.hitCount - ctx.backtraceHitCount) + " of "
                      + SSTR(ctx.hitCount) + " hits have no backtrace, so their coverage is\n"
                      "           measured over the whole aligned span instead of the CIGAR.\n"
                      "           Indels are ignored and coverage is over-counted by the gap\n"
                      "           length, so coverage_depth and breadth_of_coverage are too high.\n");
    }

    if (ctx.queryCount == 0) {
        return;
    }
    double deviationSum = 0.0;
    double maxDeviation = 0.0;
    size_t offQueries = 0;
    for (MappingTable::const_iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end(); ++it) {
        double posteriorSum = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            posteriorSum += it->second[j].posterior;
        }
        const double deviation = std::fabs(posteriorSum - 1.0);
        deviationSum += deviation;
        maxDeviation = std::max(maxDeviation, deviation);
        if (deviation > 1e-3) {
            ++offQueries;
        }
    }
    const double meanDeviation = deviationSum / static_cast<double>(ctx.queryCount);
    Debug(Debug::INFO) << "Posterior check: mean |sum(gamma) - 1| = " << meanDeviation
                       << ", max = " << maxDeviation << ", "
                       << offQueries << "/" << ctx.queryCount << " queries off by > 1e-3.\n";
    if (meanDeviation > 0.01) {
        warningBanner("You should run 'mmseqs reclassify' first.",
                      "\n"
                      "           Alignment column 3 does not sum to 1 per query, so it does not\n"
                      "           hold EM posteriors -- it still holds sequence identities.\n"
                      "           mean |sum(gamma) - 1| = " + SSTR(meanDeviation)
                      + " over " + SSTR(ctx.queryCount) + " queries.\n"
                      "\n"
                      "           Every number below is meaningless until this DB has been\n"
                      "           through 'mmseqs reclassify'.\n");
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
            stats.coverageEvenness = 0.0;
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
            stats.coverageEvenness = 0.0;
            continue;
        }
        // notCov[p] = prod over hits covering p of (1 - gamma); starts at 1 (uncovered).
        // depth[p]  = sum over hits covering p of gamma; the expected read depth at p.
        std::vector<double> notCov(stats.targetLength, 1.0);
        std::vector<double> depth(stats.targetLength, 0.0);
        // find() rather than operator[]: this runs inside a parallel region, and operator[] is a
        // non-const member, so calling it concurrently is a data race even when every key already
        // exists and no rehash can happen. Same reason EM_reclassify.cpp uses find() here.
        std::unordered_map<unsigned int, std::vector<const ReclassTaxEntry *> >::const_iterator hitIt =
                hitsByTarget.find(stats.key);
        if (hitIt != hitsByTarget.end()) {
            const std::vector<const ReclassTaxEntry *> &hits = hitIt->second;
            for (size_t h = 0; h < hits.size(); ++h) {
                const ReclassTaxEntry *entry = hits[h];
                if (entry->result.backtrace.empty() == false) {
                    EMCoverage::paintBacktrace(notCov, depth, entry->result.backtrace,
                                               entry->result.dbStartPos, entry->result.dbEndPos,
                                               entry->posterior);
                } else {
                    EMCoverage::paintSpan(notCov, depth,
                                          entry->result.dbStartPos, entry->result.dbEndPos,
                                          entry->posterior);
                }
            }
        }
        stats.coveredLength = EMCoverage::expectedCoveredLength(notCov);
        stats.coverageDepth = EMCoverage::meanDepth(depth);
        stats.coverageEvenness = EMCoverage::evenness(depth);
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

// `--drop P` means exactly "drop the lowest-abundance targets that together hold at most P% of the
// abundance mass" -- nothing else. The cut point is the user's percentage and only that.
//
// A largest-gap heuristic used to run first and pick the cut wherever the sorted abundances showed
// their biggest jump (within the same mass budget). It was removed on request: it made the number
// of dropped targets a non-monotone, data-dependent function of P that the user could neither
// predict nor control, and the gap it keyed on is meaningless for an exponential-looking abundance
// distribution, where the largest absolute gap sits at the top of the allowed range regardless of
// whether anything there is spurious.
static bool tailQuantileCutoff(std::vector<double> values,
                               double maxTailFraction,
                               size_t &tailCount) {
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
                                                             double maxDropPercentage) {
    std::unordered_set<unsigned int> dropped;
    if (stats.empty() || stats.size() < MIN_FILTER_TARGETS) {
        return dropped;
    }

    std::vector<double> abundances;
    abundances.reserve(stats.size());
    for (size_t i = 0; i < stats.size(); ++i) {
        abundances.push_back(stats[i].abundance);
    }

    const double maxTailFraction = clamp01(maxDropPercentage / 100.0);
    size_t abundanceTailCount = 0;
    if (tailQuantileCutoff(abundances, maxTailFraction, abundanceTailCount) == false) {
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
    fputs("target_key\ttarget_id\tabundance_pct\tcoverage_depth\tbreadth_of_coverage\tcoverage_evenness\tDrop(y/n)\tcovered_length\ttarget_length\n", handle);

    for (size_t i = 0; i < stats.size(); ++i) {
        const unsigned int key = stats[i].key;
        const std::string targetId = identifierForKey(targetHeaderReader, key, 0);
        // Posterior-weighted coverage (see collectTargetStats). coveredLength is the expected number
        // of covered residues and is bounded by targetLength, so breadth is naturally in [0, 1].
        const double coveredLength = stats[i].coveredLength;
        const double breadthOfCoverage = (stats[i].targetLength > 0)
                ? coveredLength / static_cast<double>(stats[i].targetLength)
                : 0.0;

        fprintf(handle, "%u\t%s\t%.12g\t%.12g\t%.12g\t%.12g\t%s\t%.12g\t%u\n",
                key,
                targetId.c_str(),
                stats[i].abundance,
                stats[i].coverageDepth,
                breadthOfCoverage,
                stats[i].coverageEvenness,
                stats[i].dropped ? "y" : "n",
                coveredLength,
                stats[i].targetLength);
    }

    fclose(handle);
}

// Expected reads assigned directly to each taxon, kept fractional.
//
// Rounding these to integers (the previous behaviour, floor(x + 0.5)) deleted every
// taxon expected to hold less than half a read, and made the clade counts a sum of
// rounded values so the reported percentages did not add up to 100. EM produces a
// continuous abundance; there is no reason to quantise it on the way out.
static std::unordered_map<TaxID, double> directCountsFromStats(const std::vector<TargetStats> &stats,
                                                               MappingReader &mapping,
                                                               size_t assignedReads) {
    std::unordered_map<TaxID, double> directCounts;
    directCounts.reserve(stats.size());
    for (size_t i = 0; i < stats.size(); ++i) {
        const TaxID taxId = mapping.lookup(stats[i].key);
        if (taxId == 0) {
            continue;
        }
        directCounts[taxId] += stats[i].abundance * static_cast<double>(assignedReads) / 100.0;
    }
    return directCounts;
}

// Clade counts in double precision. Mirrors NcbiTaxonomy::getCladeCounts, which is
// hard-coded to `unsigned int` in TaxonCounts and so cannot carry fractional reads.
// Kept local rather than templating the shared taxonomy code.
static std::unordered_map<TaxID, double> cladeCountsFromDirect(const std::unordered_map<TaxID, double> &directCounts,
                                                              NcbiTaxonomy *taxonomy) {
    std::unordered_map<TaxID, double> cladeCounts;
    for (std::unordered_map<TaxID, double>::const_iterator it = directCounts.begin(); it != directCounts.end(); ++it) {
        cladeCounts[it->first] += it->second;
        if (taxonomy->nodeExists(it->first) == false) {
            continue;
        }
        const TaxonNode *node = taxonomy->taxonNode(it->first);
        while (node->parentTaxId != node->taxId && taxonomy->nodeExists(node->parentTaxId)) {
            node = taxonomy->taxonNode(node->parentTaxId);
            cladeCounts[node->taxId] += it->second;
        }
    }
    return cladeCounts;
}

static char rankCodeOf(NcbiTaxonomy *taxonomy, const TaxonNode *node) {
    if (node == NULL) {
        return '-';
    }
    const char *rankStr = taxonomy->getString(node->rankIdx);
    if (rankStr == NULL) {
        return '-';
    }
    std::map<std::string, char>::const_iterator rankIt = NcbiShortRanks.find(std::string(rankStr));
    return (rankIt != NcbiShortRanks.end()) ? rankIt->second : '-';
}

static double lookupCount(const std::unordered_map<TaxID, double> &counts, TaxID taxId) {
    std::unordered_map<TaxID, double>::const_iterator it = counts.find(taxId);
    return (it != counts.end()) ? it->second : 0.0;
}

// `totalReads` is every read in the query DB, not just the reads that got a hit, so the
// percentages are on the same basis Kraken2 uses and the `U` row accounts for the rest
// (reads with no hit, plus mass that went to targets dropped by --drop).
static void writeKrakenReport(const std::unordered_map<TaxID, double> &directCounts,
                             const std::unordered_map<TaxID, double> &cladeCounts,
                             const std::unordered_map<TaxID, std::vector<TaxID> > &parentToChildren,
                             NcbiTaxonomy *taxonomy,
                             size_t totalReads,
                             const std::string &path) {
    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    const double denom = (totalReads > 0) ? static_cast<double>(totalReads) : 1.0;

    // Unclassified row, so the percentage column sums to 100.
    const double classified = lookupCount(cladeCounts, 1);
    const double unclassified = std::max(0.0, denom - classified);
    fprintf(handle, "%.4f\t%.2f\t%.2f\tU\t0\tunclassified\n",
            100.0 * unclassified / denom, unclassified, unclassified);

    std::vector<TaxID> stack;
    std::vector<int> depthStack;
    stack.push_back(1);
    depthStack.push_back(0);

    while (stack.empty() == false) {
        const TaxID taxId = stack.back();
        stack.pop_back();
        const int depth = depthStack.back();
        depthStack.pop_back();

        const double cladeCount = lookupCount(cladeCounts, taxId);
        if (cladeCount <= 0.0) {
            continue;
        }
        const double directCount = lookupCount(directCounts, taxId);
        const TaxonNode *node = taxonomy->taxonNode(taxId, false);
        const char *name = (node != NULL) ? taxonomy->getString(node->nameIdx) : "unclassified";

        for (int i = 0; i < depth; ++i) {
            fputs("  ", handle);
        }
        fprintf(handle, "%.4f\t%.2f\t%.2f\t%c\t%u\t%s\n",
                100.0 * cladeCount / denom, cladeCount, directCount,
                rankCodeOf(taxonomy, node), static_cast<unsigned int>(taxId), name);

        // Kraken orders siblings by descending clade count. Sort ascending and push in
        // that order so the stack pops them descending; ties break on taxId so the
        // output is byte-stable across runs.
        std::unordered_map<TaxID, std::vector<TaxID> >::const_iterator childIt = parentToChildren.find(taxId);
        if (childIt == parentToChildren.end()) {
            continue;
        }
        std::vector<std::pair<double, TaxID> > children;
        children.reserve(childIt->second.size());
        for (size_t i = 0; i < childIt->second.size(); ++i) {
            const TaxID child = childIt->second[i];
            const double count = lookupCount(cladeCounts, child);
            if (count > 0.0) {
                children.push_back(std::make_pair(count, child));
            }
        }
        std::sort(children.begin(), children.end(), [](const std::pair<double, TaxID> &lhs,
                                                      const std::pair<double, TaxID> &rhs) {
            if (lhs.first != rhs.first) {
                return lhs.first < rhs.first;
            }
            return lhs.second > rhs.second;
        });
        for (size_t i = 0; i < children.size(); ++i) {
            stack.push_back(children[i].second);
            depthStack.push_back(depth + 1);
        }
    }

    fclose(handle);
}

static void writeBrackenReport(const std::unordered_map<TaxID, double> &directCounts,
                              const std::unordered_map<TaxID, double> &cladeCounts,
                              NcbiTaxonomy *taxonomy,
                              size_t totalReads,
                              const std::string &path) {
    FILE *handle = FileUtil::openFileOrDie(path.c_str(), "w", false);
    fputs("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n", handle);

    const double denom = (totalReads > 0) ? static_cast<double>(totalReads) : 1.0;

    // Sorted output: iterating the hash map directly made the row order differ between
    // runs on identical input, which makes results impossible to diff.
    std::vector<std::pair<double, TaxID> > rows;
    rows.reserve(cladeCounts.size());
    for (std::unordered_map<TaxID, double>::const_iterator it = cladeCounts.begin(); it != cladeCounts.end(); ++it) {
        if (it->second > 0.0) {
            rows.push_back(std::make_pair(it->second, it->first));
        }
    }
    std::sort(rows.begin(), rows.end(), [](const std::pair<double, TaxID> &lhs,
                                          const std::pair<double, TaxID> &rhs) {
        if (lhs.first != rhs.first) {
            return lhs.first > rhs.first;
        }
        return lhs.second < rhs.second;
    });

    for (size_t i = 0; i < rows.size(); ++i) {
        const TaxID taxId = rows[i].second;
        const double cladeCount = rows[i].first;
        const double directCount = lookupCount(directCounts, taxId);
        const TaxonNode *node = taxonomy->taxonNode(taxId, false);
        const char *name = (node != NULL) ? taxonomy->getString(node->nameIdx) : "unclassified";
        // Bracken's taxonomy_lvl is a single-letter rank code, not the full rank name.
        fprintf(handle, "%s\t%u\t%c\t%.2f\t%.2f\t%.2f\t%.12g\n",
                name,
                static_cast<unsigned int>(taxId),
                rankCodeOf(taxonomy, node),
                directCount,
                std::max(0.0, cladeCount - directCount),
                cladeCount,
                cladeCount / denom);
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
    reportInputDiagnostics(ctx);

    Debug(Debug::INFO) << "Computing abundance from posterior...\n";
    computeAbundanceFromPosterior(ctx.mappingTable, ctx.targetSet, ctx.queryCount);

    Debug(Debug::INFO) << "Collecting per-target statistics...\n";
    std::vector<TargetStats> allTargetStats = collectTargetStats(ctx, par.threads);
    const std::unordered_set<unsigned int> dropped = selectDroppedTargets(allTargetStats,
                                                                          par.reclassifyMaxDropPercentage);
    markDroppedTargets(allTargetStats, dropped);
    convertAbundanceToPercent(allTargetStats);
    Debug(Debug::INFO) << "Writing output...\n";

    if (withTaxonomy) {
        // The percentage basis is every read in the query DB, not just the reads that got a
        // hit (ctx.queryCount). Reads with no hit at all have to be counted somewhere for the
        // report to be comparable with Kraken2/Bracken, and they land in the `U` row.
        DBReader<unsigned int> queryReader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                           DBReader<unsigned int>::USE_INDEX);
        queryReader.open(DBReader<unsigned int>::NOSORT);
        const size_t totalReads = queryReader.getSize();
        queryReader.close();
        Debug(Debug::INFO) << totalReads << " reads in the query DB, " << ctx.queryCount
                           << " with at least one hit.\n";

        std::vector<TargetStats> targetStats = allTargetStats;
        targetStats.erase(std::remove_if(targetStats.begin(), targetStats.end(), [](const TargetStats &entry) {
            return entry.dropped;
        }), targetStats.end());

        // Built once and shared by both writers; getParentToChildren() walks the whole
        // taxonomy, so calling it per writer was pure duplicated work.
        const std::unordered_map<TaxID, double> directCounts =
                directCountsFromStats(targetStats, *mapping, ctx.queryCount);
        const std::unordered_map<TaxID, double> cladeCounts =
                cladeCountsFromDirect(directCounts, taxonomy);
        const std::unordered_map<TaxID, std::vector<TaxID> > parentToChildren = taxonomy->getParentToChildren();

        writeKrakenReport(directCounts, cladeCounts, parentToChildren, taxonomy, totalReads, par.db4);
        writeBrackenReport(directCounts, cladeCounts, taxonomy, totalReads, par.db4 + ".bracken");
    } else {
        DBReader<unsigned int> targetHeaderReader((par.db2 + "_h").c_str(), (par.db2 + "_h.index").c_str(), par.threads,
                                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        targetHeaderReader.open(DBReader<unsigned int>::NOSORT);
        writeProteinStats(allTargetStats, targetHeaderReader, par.db4);
        targetHeaderReader.close();
    }

    delete mapping;
    delete taxonomy;
    reader.close();
    return EXIT_SUCCESS;
}
