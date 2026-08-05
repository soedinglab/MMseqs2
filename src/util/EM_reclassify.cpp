#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include "FastSort.h"
#include "EMCoverage.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
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

struct TargetStats {
    unsigned int key;
    unsigned int targetLength;
    double abundance;
    bool dropped;
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

// SQUAREM step length. The accelerated step is
//     x_new = x0 - 2*alpha*r + alpha^2*v,   alpha = -||r|| / ||v||
// and Varadhan & Roland's global-convergence argument requires alpha <= -1, i.e. |alpha| >= 1:
// alpha == -1 reproduces the plain EM double step (x_new == x2) and |alpha| > 1 is the actual
// acceleration. Clamping |alpha| to <= 1 instead would remove the acceleration entirely, and
// letting alpha approach 0 would make x_new == x0, i.e. stall the iteration.
// The magnitude bound is adaptive: it grows while the cap keeps binding and the accelerated step
// is accepted, and shrinks whenever the step is rejected.
static const double STEP_LENGTH_MIN = 1.0;
static const double STEP_LENGTH_INIT = 1.0;
static const double STEP_LENGTH_CAP = 1.0e4;
static const double STEP_LENGTH_GROWTH = 4.0;
static const double LOG_COMPATIBILITY_MIN = -60.0;
static const double LOG_COMPATIBILITY_MAX = 60.0;
static const double ABUNDANCE_SMOOTH_EPS = 1e-8;

// Alpha annealing. alpha is ramped in over the first iterations instead of starting at alpha_max:
//
//     alpha_k = alpha_max * (1 - exp(-(k+1)/tau)),   tau = `--alpha-tau`
//
// Rationale: theta^(0) is a crude score-share estimate, and letting a bad theta drive the E-step at
// full strength from iteration 1 can lock in a self-reinforcing (rich-get-richer) misassignment.
// alpha ~ 0 makes the first E-step almost pure alignment score, which is the more trustworthy
// signal, and abundance is phased in as theta improves. This is deterministic annealing.
//
// The cost is that the objective changes every iteration, so this is not EM on a fixed objective:
// there is no monotone-ascent guarantee, and SQUAREM's convergence argument assumes a fixed
// fixed-point map. squarem() already works around the first point by scoring the accelerated step
// against x2 under the *same* alpha_k rather than across iterations.
//
// Whether annealing actually helps has never been measured -- there was no control arm while tau
// was hard-coded. `--alpha-tau 0` is that control: it disables annealing entirely (alpha_k =
// alpha_max for every k) so the two can be swept against each other.
static double annealedAbundanceExponent(double alphaMax, double tau, int iter) {
    if (tau <= 0.0) {
        return alphaMax;
    }
    return alphaMax * (1.0 - std::exp(-static_cast<double>(iter + 1) / tau));
}
// The coverage-prior weight w is `--cov-prior` (par.reclassifyCoveragePrior). It is a scale-free
// prior strength: the prior's share of the abundance mass is w/(1+w), independent of read and
// target count. See the kappa derivation in emUpdate().

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

            // No `+ 1` variants: parseAlignmentRecord() has no case for 12 or 16 columns
            // and exits before these flags could be consulted.
            if (columns == Matcher::ALN_RES_WITH_BT_COL_CNT
                || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT) {
                ctx.hasBacktrace = true;
            }
            if (columns == Matcher::ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT
                || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT) {
                ctx.hasOrfPosition = true;
            }

            Matcher::result_t result = Matcher::parseAlignmentRecord(data, true);
            records.push_back(ReclassTaxEntry{result, 0.0, static_cast<double>(result.seqId), 0.0});
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
            scoreSum += std::max(static_cast<double>(it->second[j].result.score), 0.0);
        }
        if (scoreSum <= 0.0) {
            continue;
        }
        for (size_t j = 0; j < it->second.size(); ++j) {
            const double nonNegativeScore = std::max(static_cast<double>(it->second[j].result.score), 0.0);
            initAbundance[it->second[j].result.dbKey] += nonNegativeScore / scoreSum;
        }
    }

    if (queryCount > 0) {
        const double denom = static_cast<double>(queryCount);
        for (std::unordered_map<unsigned int, double>::iterator it = initAbundance.begin(); it != initAbundance.end(); ++it) {
            it->second /= denom;
        }
    }
    double totalAbundance = 0.0;
    for (std::unordered_map<unsigned int, double>::const_iterator it = initAbundance.begin(); it != initAbundance.end(); ++it) {
        totalAbundance += it->second;
    }
    if (totalAbundance > 0.0) {
        for (std::unordered_map<unsigned int, double>::iterator it = initAbundance.begin(); it != initAbundance.end(); ++it) {
            it->second /= totalAbundance;
        }
    }

    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            it->second[j].abundance = initAbundance[it->second[j].result.dbKey];
        }
    }
}

struct TargetHitRef {
    const ReclassTaxEntry *entry;
    double weight;
};

// Coverage confidence C_t: how broadly and how evenly a target is covered by its hits. C_t = F_t *
// E_t, both in [0, 1] and both defined over the full target length L_t (see EMCoverage.h for the
// formulas and for why each denominator is L_t rather than the hit span or the covered-position
// count):
//
//   F_t = breadth of the *score-share-weighted* depth -- how much of the target is covered, and how
//         confidently, since a hit contributes its within-query bit-score share w_qt rather than 1.
//   E_t = evenness of the *raw read* depth -- how uniformly that coverage is spread, so a repeat
//         pileup on a small part of the target is penalised.
//
// The two factors read different depth arrays on purpose. Weighting the evenness term as well made
// it respond to hit ambiguity instead of to depth raggedness, double-counting something F_t already
// prices in. Feeding the evenness a covered-positions denominator made it score a pure repeat
// pileup *above* an end-to-end-covered target; both regressions are documented with numbers at
static void initCoverageConfidence(MappingTable &mappingTable,
                                   const std::unordered_set<unsigned int> &targetSet,
                                   int threads) {
    (void)threads;  // only read by the num_threads() clause, which vanishes without OPENMP
    std::unordered_map<unsigned int, unsigned int> targetLength;
    std::unordered_map<unsigned int, std::vector<TargetHitRef> > hitsByTarget;
    targetLength.reserve(targetSet.size());
    hitsByTarget.reserve(targetSet.size());

    for (std::unordered_set<unsigned int>::const_iterator it = targetSet.begin(); it != targetSet.end(); ++it) {
        targetLength[*it] = 0;
        hitsByTarget.emplace(*it, std::vector<TargetHitRef>());
    }

    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        // Clamp at 0 to match initAbundance: a negative score must not make the share negative or
        // cancel other hits out of the denominator.
        double scoreSum = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            scoreSum += std::max(static_cast<double>(it->second[j].result.score), 0.0);
        }
        for (size_t j = 0; j < it->second.size(); ++j) {
            const unsigned int target = it->second[j].result.dbKey;
            targetLength[target] = it->second[j].result.dbLen;
            const double score = std::max(static_cast<double>(it->second[j].result.score), 0.0);
            const double weight = (scoreSum > 0.0) ? (score / scoreSum) : 0.0;
            hitsByTarget[target].push_back(TargetHitRef{&it->second[j], weight});
        }
    }

    std::unordered_map<unsigned int, double> coverageFraction;
    coverageFraction.reserve(targetSet.size());
    const std::vector<unsigned int> targetList = targetListFromSet(targetSet);
    std::vector<double> coverageFractionByIndex(targetList.size(), 0.0);

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (size_t i = 0; i < targetList.size(); ++i) {
        const unsigned int target = targetList[i];
        // find() rather than operator[]: this runs inside a parallel region, where an accidental
        // insertion would rehash the map under the other threads.
        std::unordered_map<unsigned int, unsigned int>::const_iterator lenIt = targetLength.find(target);
        const int len = (lenIt != targetLength.end()) ? static_cast<int>(lenIt->second) : 0;
        if (len <= 0) {
            coverageFractionByIndex[i] = 0.0;
            continue;
        }
        // weightedDepth carries the score-share weights (breadth), readDepth the raw hit counts
        // (evenness). See the C_t comment above for why the two terms must not share an array.
        std::vector<double> weightedDepth(static_cast<size_t>(len), 0.0);
        std::vector<double> readDepth(static_cast<size_t>(len), 0.0);

        std::unordered_map<unsigned int, std::vector<TargetHitRef> >::const_iterator hitIt = hitsByTarget.find(target);
        if (hitIt != hitsByTarget.end()) {
            const std::vector<TargetHitRef> &hits = hitIt->second;
            for (size_t h = 0; h < hits.size(); ++h) {
                const Matcher::result_t &result = hits[h].entry->result;
                // Coverage confidence uses the full aligned span (it is not CIGAR-aware); only the
                // breadth/depth columns in `abundance` refine this to M columns. paintSpanWeighted
                // clips to [0, L) and is strand-safe.
                EMCoverage::paintSpanWeighted(weightedDepth, readDepth,
                                              result.dbStartPos, result.dbEndPos, hits[h].weight);
            }
        }

        coverageFractionByIndex[i] = clamp01(EMCoverage::breadth(weightedDepth)
                                             * EMCoverage::evenness(readDepth));
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

static double maxQueryBitScore(const std::vector<ReclassTaxEntry> &entries) {
    double maxScore = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i < entries.size(); ++i) {
        maxScore = std::max(maxScore, static_cast<double>(entries[i].result.score));
    }
    return maxScore;
}

// phi_qt = exp(beta_bit * (score - maxScore)). The bit score is the only alignment evidence
// the E-step uses.
//
// A `+ beta_seqid * seqId` term was removed on 2026-08-04, following the same fate as the
// earlier `+ beta_cov * cov` term. Swept over beta_seqid in {0, 0.25, 0.5, 1, 2, 4} at three
// taxonomic levels x three read lengths, the whole 16x range moved Spearman by 0.0003-0.0023
// and Bray-Curtis by <=0.0005 -- against 0.0132 for merely doubling alpha_max -- and the
// optimum pointed in opposite directions per level (family preferred 0, order preferred 4).
// Deleting the term (beta_seqid = 0) changed every metric by <=0.0007 with mixed sign, so it
// costs nothing measurable. The one place the term did measurably act was family L80 at
// beta_seqid = 4, where it *degraded* Spearman by 0.0086 by letting a bounded [0,1] feature
// override 16 bits of alignment score.
//
// Removing it also drops two hazards: seqId is quantised to 1e-3 by the aligner's
// fastSeqIdToBuffer, and if `reclassify` were ever run on its own output, column 3 holds a
// posterior rather than a sequence identity, which would have silently changed the model.
static double compatibilityLogTerm(const ReclassTaxEntry &entry,
                                   double queryMaxScore,
                                   double betaBit) {
    const double r = betaBit * (static_cast<double>(entry.result.score) - queryMaxScore);
    return std::max(LOG_COMPATIBILITY_MIN, std::min(LOG_COMPATIBILITY_MAX, r));
}

static void computePosterior(MappingTable &mappingTable,
                             double betaBit,
                             double abundanceExponent) {
    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        const double queryMaxScore = maxQueryBitScore(it->second);
        std::vector<double> numerators(it->second.size(), 0.0);
        double denom = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            const double phi = std::exp(compatibilityLogTerm(it->second[j], queryMaxScore, betaBit));
            const double abundance = std::max(it->second[j].abundance, ABUNDANCE_SMOOTH_EPS);
            const double weighted = phi * std::pow(abundance, abundanceExponent);
            numerators[j] = weighted;
            denom += weighted;
        }
        for (size_t j = 0; j < it->second.size(); ++j) {
            it->second[j].posterior = (denom > 0.0) ? (numerators[j] / denom) : 0.0;
        }
    }
}

static double logLikelihood(const MappingTable &mappingTable,
                            double betaBit,
                            double abundanceExponent,
                            size_t queryCount) {
    if (queryCount == 0) {
        return 0.0;
    }

    double ll = 0.0;
    for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        const double queryMaxScore = maxQueryBitScore(it->second);
        double mixture = 0.0;
        for (size_t j = 0; j < it->second.size(); ++j) {
            const double phi = std::exp(compatibilityLogTerm(it->second[j], queryMaxScore, betaBit));
            const double abundance = std::max(it->second[j].abundance, ABUNDANCE_SMOOTH_EPS);
            mixture += phi * std::pow(abundance, abundanceExponent);
        }
        ll += std::log(mixture > 0.0 ? mixture : 1e-300);
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
                                    double betaBit,
                                    const std::unordered_map<unsigned int, double> &fixedCoverageConfidence,
                                    double confidenceSum,
                                    double coveragePriorWeight,
                                    const std::vector<unsigned int> &targetList,
                                    size_t queryCount,
                                    double abundanceExponent) {
    computePosterior(mappingTable, betaBit, abundanceExponent);

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

    // Dirichlet-style coverage prior with a scale-free strength. The data term sums to the number of
    // queries (sum_t R_t = N_q) while the raw confidences sum to O(|T|), so adding C_t directly made
    // the prior's influence depend on the read-to-target ratio: on a small run (few hundred reads,
    // thousands of targets) the prior outweighed the likelihood by an order of magnitude, and on a
    // deep run it was negligible. Normalising C_t to a distribution and scaling the total pseudo-count
    // mass to kappa = w * N_q pins the prior's share of the posterior mass at exactly w/(1+w),
    // independent of N_q and |T|.
    const double kappa = coveragePriorWeight * static_cast<double>(queryCount);
    double denom = 0.0;
    for (std::unordered_map<unsigned int, double>::iterator it = nextAbundance.begin(); it != nextAbundance.end(); ++it) {
        const std::unordered_map<unsigned int, double>::const_iterator fixed = fixedCoverageConfidence.find(it->first);
        const double confidence = (fixed != fixedCoverageConfidence.end()) ? fixed->second : 0.0;
        // No constant floor is added here on purpose. A per-target epsilon can never be zeroed out,
        // which destroys the sparsity EM is supposed to produce; log(0) is already guarded by the
        // max(pi, ABUNDANCE_SMOOTH_EPS) inside computePosterior()/logLikelihood().
        if (confidenceSum > 0.0) {
            it->second += kappa * (confidence / confidenceSum);
        }
        denom += it->second;
    }
    if (denom > 0.0) {
        for (std::unordered_map<unsigned int, double>::iterator it = nextAbundance.begin(); it != nextAbundance.end(); ++it) {
            it->second /= denom;
        }
    }

    for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
        for (size_t j = 0; j < it->second.size(); ++j) {
            it->second[j].abundance = nextAbundance[it->second[j].result.dbKey];
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

static void squarem(ReclassTaxContext &ctx,
                    double betaBit,
                    int maxIter,
                    double tol,
                    double alphaMax,
                    double alphaTau,
                    double coveragePriorWeight,
                    int threads) {
    if (ctx.queryCount == 0 || ctx.targetSet.empty()) {
        return;
    }

    initAbundance(ctx.mappingTable, ctx.targetSet, ctx.queryCount);

    std::unordered_map<unsigned int, double> fixedCoverageConfidence;
    double confidenceSum = 0.0;
    if (coveragePriorWeight > 0.0) {
        initCoverageConfidence(ctx.mappingTable, ctx.targetSet, threads);
        fixedCoverageConfidence.reserve(ctx.targetSet.size());
        for (MappingTable::const_iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end(); ++it) {
            for (size_t j = 0; j < it->second.size(); ++j) {
                fixedCoverageConfidence[it->second[j].result.dbKey] = it->second[j].coverageConfidence;
            }
        }
        // fixedCoverageConfidence never changes during EM, so the normaliser is computed once.
        for (std::unordered_map<unsigned int, double>::const_iterator it = fixedCoverageConfidence.begin();
             it != fixedCoverageConfidence.end(); ++it) {
            confidenceSum += it->second;
        }
        // Report the balance the prior scaling exists to control. The data term contributes N_q of mass
        // per M-step and the prior contributes w*N_q, so the prior's share is w/(1+w) regardless of the
        // read-to-target ratio. Before the rescaling the prior added a raw sum(C_t) of mass against the
        // data's N_q, which meant its influence silently tracked N_q/|T|: ~0.05% on a 1M-read run over a
        // few hundred targets, but dominant on a few-hundred-read run over thousands of targets.
        Debug(Debug::INFO) << "Coverage prior: weight=" << coveragePriorWeight
                           << " sum(C_t)=" << confidenceSum
                           << " over " << ctx.targetSet.size() << " targets, "
                           << ctx.queryCount << " queries -> prior mass share "
                           << (coveragePriorWeight / (1.0 + coveragePriorWeight)) << "\n";
    } else {
        // confidenceSum stays 0, which makes emUpdate() skip the prior term entirely. Skipping
        // initCoverageConfidence() also skips the most memory-hungry phase of reclassify.
        Debug(Debug::INFO) << "Coverage prior disabled (--cov-prior 0); skipping coverage confidence.\n";
    }

    const std::vector<unsigned int> targetList = targetListFromSet(ctx.targetSet);
    std::vector<double> x0 = abundanceVectorFromTable(ctx.mappingTable, targetList);

    // abundanceExponent is annealed over the first iterations, so the parameter-change test must
    // not be allowed to declare convergence while the exponent is still moving. After 4*tau
    // iterations the exponent is within 2% of alphaMax (1 - e^-4 = 0.9817). With annealing off
    // (--alpha-tau 0) there is nothing to wait for and this falls back to the plain floor of 6.
    const int minIterations = std::max(6, static_cast<int>(std::ceil(4.0 * std::max(0.0, alphaTau))));
    double stepLengthMax = STEP_LENGTH_INIT;
    double acceptedLl = 0.0;
    bool converged = false;
    int lastIter = 0;

    for (int iter = 0; iter < maxIter; ++iter) {
        lastIter = iter;
        const double abundanceExponent = annealedAbundanceExponent(alphaMax, alphaTau, iter);
        const std::vector<double> x1 = emUpdate(ctx.mappingTable,
                                                betaBit,
                                                fixedCoverageConfidence,
                                                confidenceSum,
                                                coveragePriorWeight,
                                                targetList,
                                                ctx.queryCount,
                                                abundanceExponent);
        const std::vector<double> x2 = emUpdate(ctx.mappingTable,
                                                betaBit,
                                                fixedCoverageConfidence,
                                                confidenceSum,
                                                coveragePriorWeight,
                                                targetList,
                                                ctx.queryCount,
                                                abundanceExponent);

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

        // alpha is confined to [-stepLengthMax, -1] (see the STEP_LENGTH_* comment): |alpha| >= 1
        // guarantees at least the plain EM double step, and the upper magnitude bound is what keeps
        // the extrapolation from leaving the neighbourhood where the local linear model holds.
        double accel = (normV <= 0.0) ? -STEP_LENGTH_MIN : -(normR / normV);
        if (std::isfinite(accel) == false) {
            accel = -STEP_LENGTH_MIN;
        }
        accel = std::max(-stepLengthMax, std::min(-STEP_LENGTH_MIN, accel));
        const bool stepCapBinding = (accel <= (-stepLengthMax + 1e-12));

        std::vector<double> xNew(x0.size(), 0.0);
        for (size_t i = 0; i < x0.size(); ++i) {
            xNew[i] = x0[i] - 2.0 * accel * r[i] + accel * accel * v[i];
        }
        xNew = projectSimplex(xNew);

        // Score the accelerated point against its own fallback candidate x2 (the plain EM double
        // step) under the *same* abundanceExponent. Comparing against the previous iteration's
        // log-likelihood would be meaningless: the objective is sum_q log sum_t phi_qt*theta_t^alpha,
        // alpha is annealed upward every iteration, and theta_t^alpha decreases monotonically in
        // alpha for theta_t < 1. The objective therefore drifts downward on its own -- for 1e5
        // targets an alpha step of 0.11 costs about 1.3 nats per query, dwarfing any 1e-9 tolerance
        // -- so a cross-iteration test would reject essentially every accelerated step.
        // logLikelihood() reads only `abundance`, so no computePosterior() is needed to evaluate it.
        setAbundance(ctx.mappingTable, targetList, x2);
        const double emLl = logLikelihood(ctx.mappingTable, betaBit, abundanceExponent, ctx.queryCount);
        setAbundance(ctx.mappingTable, targetList, xNew);
        const double accelLl = logLikelihood(ctx.mappingTable, betaBit, abundanceExponent, ctx.queryCount);

        // A NaN accelLl fails this test, so a degenerate extrapolation falls back to x2.
        const bool accelAccepted = std::isfinite(accelLl) && (accelLl >= (emLl - 1e-12));
        if (accelAccepted) {
            acceptedLl = accelLl;
            if (stepCapBinding) {
                stepLengthMax = std::min(STEP_LENGTH_CAP, stepLengthMax * STEP_LENGTH_GROWTH);
            }
        } else {
            xNew = x2;
            setAbundance(ctx.mappingTable, targetList, x2);
            acceptedLl = emLl;
            stepLengthMax = std::max(STEP_LENGTH_INIT, stepLengthMax / STEP_LENGTH_GROWTH);
        }
        // The posteriors that get written out have to match the accepted abundance vector.
        computePosterior(ctx.mappingTable, betaBit, abundanceExponent);

        double parameterChange = 0.0;
        for (size_t i = 0; i < x0.size(); ++i) {
            parameterChange = std::max(parameterChange, std::fabs(xNew[i] - x0[i]));
        }

        // '-' the accelerated step was accepted, '.' it fell back to the plain EM double step.
        Debug(Debug::INFO) << (accelAccepted ? "-" : ".");
        std::cout << std::flush;
        x0 = xNew;
        if (parameterChange < tol && (iter + 1) >= minIterations) {
            converged = true;
            break;
        }
    }

    Debug(Debug::INFO) << "\n";
    if (converged) {
        Debug(Debug::INFO) << "Converged after " << (lastIter + 1) << " iterations (logLikelihood="
                           << acceptedLl << ").\n";
    } else {
        Debug(Debug::INFO) << "Reached the iteration limit (" << maxIter << ") without meeting tol="
                           << tol << " (logLikelihood=" << acceptedLl << ").\n";
    }
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
            stats.dropped = false;
        }
    }

    std::vector<TargetStats> out;
    out.reserve(statsByTarget.size());
    for (std::unordered_map<unsigned int, TargetStats>::const_iterator it = statsByTarget.begin(); it != statsByTarget.end(); ++it) {
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

static void printAbundanceDistribution(const std::vector<TargetStats> &stats) {
    if (stats.empty()) {
        Debug(Debug::INFO) << "Abundance distribution: no targets.\n";
        return;
    }

    std::vector<double> values;
    values.reserve(stats.size());
    for (size_t i = 0; i < stats.size(); ++i) {
        values.push_back(stats[i].abundance);
    }
    std::sort(values.begin(), values.end());

    const auto quantile = [&values](double q) -> double {
        const size_t idx = static_cast<size_t>(q * static_cast<double>(values.size() - 1));
        return values[idx];
    };

    const double totalMass = std::accumulate(values.begin(), values.end(), 0.0);
    const auto cumulativeCountAtMassFraction = [&values, totalMass](double massFraction, size_t &count, double &countFrac) {
        count = 0;
        countFrac = 0.0;
        if (values.empty() || totalMass <= 0.0 || massFraction <= 0.0) {
            return;
        }

        const double targetMass = massFraction * totalMass;
        double cumulativeMass = 0.0;
        for (size_t i = 0; i < values.size(); ++i) {
            if ((cumulativeMass + values[i]) <= targetMass) {
                cumulativeMass += values[i];
                ++count;
            } else {
                break;
            }
        }
        countFrac = static_cast<double>(count) / static_cast<double>(values.size());
    };

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(8);
    oss << "Abundance distribution (targets=" << values.size() << "):"
        << " min=" << values.front()
        << " p25=" << quantile(0.25)
        << " p50=" << quantile(0.50)
        << " p75=" << quantile(0.75)
        << " p90=" << quantile(0.90)
        << " p95=" << quantile(0.95)
        << " p99=" << quantile(0.99)
        << " max=" << values.back();
    Debug(Debug::INFO) << oss.str() << "\n";

    const double cutoffs[] = {0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99};
    std::ostringstream cumulativeOss;
    cumulativeOss << std::fixed << std::setprecision(8);
    cumulativeOss << "Abundance cumulative:";
    for (size_t i = 0; i < sizeof(cutoffs) / sizeof(cutoffs[0]); ++i) {
        size_t count = 0;
        double countFrac = 0.0;
        cumulativeCountAtMassFraction(cutoffs[i], count, countFrac);
        cumulativeOss << " <=" << cutoffs[i]
                      << "[count=" << count
                      << ",countFrac=" << countFrac << "]";
    }
    Debug(Debug::INFO) << cumulativeOss.str() << "\n";
}

static double clamp01(double value) {
    return std::max(0.0, std::min(1.0, value));
}

static bool compareByPosteriorThenBitScore(const ReclassTaxEntry &a, const ReclassTaxEntry &b) {
    if (a.posterior != b.posterior) {
        return a.posterior > b.posterior;
    }
    if (a.result.score != b.result.score) {
        return a.result.score > b.result.score;
    }
    return Matcher::compareHits(a.result, b.result);
}

static void writeReclassifiedDb(const ReclassTaxContext &ctx,
                                int dbType,
                                const std::string &outDb,
                                const std::string &outIndex,
                                int threads,
                                bool compress) {
    DBWriter writer(outDb.c_str(), outIndex.c_str(), threads, compress, dbType);
    writer.open();

    Debug::Progress progress(ctx.queryOrder.size());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[1024 + 32768 * 4];

#pragma omp for schedule(dynamic, 5)
        for (size_t i = 0; i < ctx.queryOrder.size(); ++i) {
            progress.updateProgress();
            const unsigned int queryKey = ctx.queryOrder[i];
            MappingTable::const_iterator recordsIt = ctx.mappingTable.find(queryKey);
            if (recordsIt == ctx.mappingTable.end()) {
                continue;
            }

            std::vector<ReclassTaxEntry> records = recordsIt->second;
            SORT_SERIAL(records.begin(), records.end(), compareByPosteriorThenBitScore);

            writer.writeStart(thread_idx);
            for (size_t j = 0; j < records.size(); ++j) {
                Matcher::result_t res = records[j].result;
                res.seqId = static_cast<float>(records[j].posterior);
                // The backtrace was parsed with readCompressed=true, so res.backtrace already holds
                // the compressed CIGAR (e.g. "55M"). resultToBuffer(..., compress, addOrfPosition):
                // write it verbatim (compress=false) instead of re-compressing it into garbage
                // (e.g. "46M" -> "0M14161M"), and forward hasOrfPosition to the addOrfPosition slot.
                // highPrecisionSeqId=true is required: column 3 carries the EM posterior here, and
                // the default writer is a 3-decimal fixed-point encoder that truncates. With it,
                // every gamma < 0.001 would become exactly 0 and the per-read posteriors would no
                // longer sum to 1 -- with an error proportional to the read's hit count, which
                // biases abundance toward reads that map uniquely.
                size_t len = Matcher::resultToBuffer(buffer, res, ctx.hasBacktrace, false, ctx.hasOrfPosition, true);
                writer.writeAdd(buffer, len, thread_idx);
            }
            writer.writeEnd(queryKey, thread_idx);
        }
    }

    writer.close();
}
}

int emreclassify(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    ReclassTaxContext ctx;
    loadAlignmentDb(reader, ctx);
    Debug(Debug::INFO) << "Loaded " << ctx.queryCount << " queries with hits and " << ctx.targetSet.size() << " unique targets.\n";

    squarem(ctx,
            par.reclassifyLambda,
            par.reclassifyMaxIterations,
            par.reclassifyTolerance,
            par.reclassifyAlpha,
            par.reclassifyAlphaTau,
            par.reclassifyCoveragePrior,
            par.threads);

    // reclassify only runs EM; low-abundance target dropping is handled by `mmseqs abundance --drop`.
    std::vector<TargetStats> allTargetStats = collectTargetStats(ctx);
    printAbundanceDistribution(allTargetStats);

    writeReclassifiedDb(ctx, reader.getDbtype(), par.db4, par.db4Index, par.threads, par.compressed);

    reader.close();
    return EXIT_SUCCESS;
}
