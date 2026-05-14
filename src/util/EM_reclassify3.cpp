// #include "Parameters.h"
// #include "DBReader.h"
// #include "DBWriter.h"
// #include "Debug.h"
// #include "Util.h"
// #include "Matcher.h"
// #include "FastSort.h"

// #include <algorithm>
// #include <cmath>
// #include <cctype>
// #include <limits>
// #include <numeric>
// #include <unordered_map>
// #include <unordered_set>
// #include <vector>
// #ifdef OPENMP
// #include <omp.h>
// #endif

// namespace {
// struct ReclassTaxEntry {
//     Matcher::result_t result;
//     double abundance;
//     double posterior;
//     double coverageConfidence;
// };

// typedef std::unordered_map<unsigned int, std::vector<ReclassTaxEntry> > MappingTable;

// struct Interval {
//     int start;
//     int end;
// };

// struct TargetStats {
//     unsigned int key;
//     unsigned int targetLength;
//     double abundance;
//     double coverageConfidence;
//     bool dropped;
//     std::vector<Interval> intervals;
// };

// struct ReclassTaxContext {
//     MappingTable mappingTable;
//     std::vector<unsigned int> queryOrder;
//     std::unordered_set<unsigned int> targetSet;
//     size_t queryCount;
//     bool hasBacktrace;
//     bool hasOrfPosition;

//     ReclassTaxContext() : queryCount(0), hasBacktrace(false), hasOrfPosition(false) {}
// };

// static const double STEP_MIN = -1.0;
// static const double STEP_MAX = 1.0;
// static const double EPS = 1e-12;
// static const size_t MIN_FILTER_TARGETS = 20;
// static const size_t MIN_TAIL_TARGETS = 2;

// static double clamp01(double value);

// static std::vector<unsigned int> targetListFromSet(const std::unordered_set<unsigned int> &targets) {
//     std::vector<unsigned int> out(targets.begin(), targets.end());
//     SORT_SERIAL(out.begin(), out.end());
//     return out;
// }

// static void loadAlignmentDb(DBReader<unsigned int> &reader, ReclassTaxContext &ctx) {
//     Debug::Progress progress(reader.getSize());
//     const char *entry[255];

//     for (size_t i = 0; i < reader.getSize(); ++i) {
//         progress.updateProgress();
//         const unsigned int queryKey = reader.getDbKey(i);
//         char *data = reader.getData(i, 0);

//         if (reader.getEntryLen(i) <= 1) {
//             continue;
//         }

//         std::vector<ReclassTaxEntry> &records = ctx.mappingTable[queryKey];
//         if (records.empty()) {
//             ctx.queryOrder.push_back(queryKey);
//         }
//         while (*data != '\0') {
//             const size_t columns = Util::getWordsOfLine(data, entry, 255);
//             if (columns < Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
//                 Debug(Debug::ERROR) << "Invalid alignment result record in query " << queryKey << ".\n";
//                 EXIT(EXIT_FAILURE);
//             }

//             if (columns == Matcher::ALN_RES_WITH_BT_COL_CNT || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT
//                 || columns == Matcher::ALN_RES_WITH_BT_COL_CNT + 1 || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT + 1) {
//                 ctx.hasBacktrace = true;
//             }
//             if (columns == Matcher::ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT
//                 || columns == Matcher::ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT + 1 || columns == Matcher::ALN_RES_WITH_ORF_AND_BT_COL_CNT + 1) {
//                 ctx.hasOrfPosition = true;
//             }

//             Matcher::result_t result = Matcher::parseAlignmentRecord(data, true);
//             records.push_back(ReclassTaxEntry{result, 0.0, 0.0, 0.0});
//             ctx.targetSet.insert(result.dbKey);
//             data = Util::skipLine(data);
//         }
//     }

//     ctx.queryCount = ctx.mappingTable.size();
// }

// static void initAbundance(MappingTable &mappingTable, const std::unordered_set<unsigned int> &targetSet, size_t queryCount) {
//     std::unordered_map<unsigned int, double> initAbundance;
//     initAbundance.reserve(targetSet.size());
//     for (std::unordered_set<unsigned int>::const_iterator it = targetSet.begin(); it != targetSet.end(); ++it) {
//         initAbundance[*it] = 0.0;
//     }

//     for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         double scoreSum = 0.0;
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             scoreSum += it->second[j].result.score;
//         }
//         if (scoreSum <= 0.0) {
//             continue;
//         }
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             initAbundance[it->second[j].result.dbKey] += it->second[j].result.score / scoreSum;
//         }
//     }

//     if (queryCount > 0) {
//         const double denom = static_cast<double>(queryCount);
//         for (std::unordered_map<unsigned int, double>::iterator it = initAbundance.begin(); it != initAbundance.end(); ++it) {
//             it->second /= denom;
//         }
//     }

//     for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             it->second[j].abundance = initAbundance[it->second[j].result.dbKey];
//         }
//     }
// }

// struct TargetHitRef {
//     const ReclassTaxEntry *entry;
//     double score;
//     double weight;
// };

// static void initCoverageConfidence(MappingTable &mappingTable,
//                                    const std::unordered_set<unsigned int> &targetSet,
//                                    int threads) {
//     (void)threads;
//     std::unordered_map<unsigned int, int> targetMin;
//     std::unordered_map<unsigned int, int> targetMax;
//     std::unordered_map<unsigned int, std::vector<TargetHitRef> > hitsByTarget;

//     for (std::unordered_set<unsigned int>::const_iterator it = targetSet.begin(); it != targetSet.end(); ++it) {
//         targetMin[*it] = std::numeric_limits<int>::max();
//         targetMax[*it] = std::numeric_limits<int>::min();
//         hitsByTarget.emplace(*it, std::vector<TargetHitRef>());
//     }

//     for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         double scoreSum = 0.0;
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             scoreSum += static_cast<double>(it->second[j].result.score);
//         }
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             const unsigned int target = it->second[j].result.dbKey;
//             if (it->second[j].result.dbStartPos < targetMin[target]) {
//                 targetMin[target] = it->second[j].result.dbStartPos;
//             }
//             if (it->second[j].result.dbEndPos > targetMax[target]) {
//                 targetMax[target] = it->second[j].result.dbEndPos;
//             }
//             const double score = static_cast<double>(it->second[j].result.score);
//             const double weight = (scoreSum > 0.0) ? (score / scoreSum) : 0.0;
//             hitsByTarget[target].push_back(TargetHitRef{&it->second[j], score, weight});
//         }
//     }

//     std::unordered_map<unsigned int, double> coverageFraction;
//     coverageFraction.reserve(targetSet.size());
//     const std::vector<unsigned int> targetList = targetListFromSet(targetSet);
//     std::vector<double> coverageFractionByIndex(targetList.size(), 0.0);

// #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
//     for (size_t i = 0; i < targetList.size(); ++i) {
//         const unsigned int target = targetList[i];
//         const int startPos = targetMin[target];
//         const int endPos = targetMax[target];
//         const int len = (endPos >= startPos) ? (endPos - startPos + 1) : 1;
//         std::vector<double> cov(static_cast<size_t>(len), 0.0);
//         std::vector<double> covConf(static_cast<size_t>(len), 0.0);

//         std::unordered_map<unsigned int, std::vector<TargetHitRef> >::const_iterator hitIt = hitsByTarget.find(target);
//         if (hitIt != hitsByTarget.end()) {
//             const std::vector<TargetHitRef> &hits = hitIt->second;
//             for (size_t h = 0; h < hits.size(); ++h) {
//                 const Matcher::result_t &result = hits[h].entry->result;
//                 const int targetLen = result.dbEndPos - result.dbStartPos + 1;
//                 if (targetLen <= 0) {
//                     continue;
//                 }

//                 const double mq = hits[h].score / static_cast<double>(targetLen);
//                 const int start = std::max(0, result.dbStartPos - startPos);
//                 const int end = std::min(len - 1, result.dbEndPos - startPos);
//                 for (int pos = start; pos <= end; ++pos) {
//                     cov[static_cast<size_t>(pos)] += mq;
//                     covConf[static_cast<size_t>(pos)] += hits[h].weight;
//                 }
//             }
//         }

//         double covered = 0.0;
//         double squaredCovered = 0.0;
//         for (size_t pos = 0; pos < covConf.size(); ++pos) {
//             const double clipped = std::min(1.0, covConf[pos]);
//             covered += clipped;
//             squaredCovered += clipped * clipped;
//         }
//         const double fraction = covered / static_cast<double>(len);
//         const double hhi = (covered > 0.0) ? (squaredCovered / (covered * covered)) : 1.0;
//         const double concentrationPenalty = 1.0 - hhi;
//         coverageFractionByIndex[i] = clamp01(fraction * concentrationPenalty);
//     }

//     for (size_t i = 0; i < targetList.size(); ++i) {
//         coverageFraction[targetList[i]] = coverageFractionByIndex[i];
//     }

//     for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             const unsigned int target = it->second[j].result.dbKey;
//             std::unordered_map<unsigned int, double>::const_iterator cf = coverageFraction.find(target);
//             it->second[j].coverageConfidence = (cf != coverageFraction.end()) ? cf->second : 0.0;
//         }
//     }
// }

// static double scoreTerm(const ReclassTaxEntry &entry, double lambda, double alpha, double gamma) {
//     const double bitScore = static_cast<double>(entry.result.score);
//     const double abundance = std::max(entry.abundance, EPS);
//     const double coverageConfidence = std::max(entry.coverageConfidence, EPS);
//     return std::exp(lambda * bitScore) * std::pow(abundance, alpha) * std::pow(coverageConfidence, gamma);
// }

// static void computePosterior(MappingTable &mappingTable, double lambda, double alpha, double gamma) {
//     for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         double denom = 0.0;
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             denom += scoreTerm(it->second[j], lambda, alpha, gamma);
//         }
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             const double value = scoreTerm(it->second[j], lambda, alpha, gamma);
//             it->second[j].posterior = (denom > 0.0) ? (value / denom) : 0.0;
//         }
//     }
// }

// static double logLikelihood(const MappingTable &mappingTable, double lambda, double alpha, double gamma, size_t queryCount) {
//     if (queryCount == 0) {
//         return 0.0;
//     }

//     double ll = 0.0;
//     for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         double denom = 0.0;
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             const double value = scoreTerm(it->second[j], lambda, alpha, gamma);
//             denom += value;
//             if (it->second[j].posterior > 0.0) {
//                 ll += it->second[j].posterior * std::log(value > 0.0 ? value : 1e-300);
//             }
//         }
//         ll -= std::log(denom > 0.0 ? denom : 1e-300);
//     }
//     return ll / static_cast<double>(queryCount);
// }

// static std::vector<double> abundanceVectorFromTable(const MappingTable &mappingTable, const std::vector<unsigned int> &targetList) {
//     std::unordered_map<unsigned int, double> abundance;
//     abundance.reserve(targetList.size());
//     for (size_t i = 0; i < targetList.size(); ++i) {
//         abundance[targetList[i]] = 0.0;
//     }

//     for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             abundance[it->second[j].result.dbKey] = it->second[j].abundance;
//         }
//     }

//     std::vector<double> out;
//     out.reserve(targetList.size());
//     for (size_t i = 0; i < targetList.size(); ++i) {
//         out.push_back(abundance[targetList[i]]);
//     }
//     return out;
// }

// static void setAbundance(MappingTable &mappingTable, const std::vector<unsigned int> &targetList, const std::vector<double> &abundanceVector) {
//     std::unordered_map<unsigned int, double> abundance;
//     abundance.reserve(targetList.size());
//     for (size_t i = 0; i < targetList.size(); ++i) {
//         abundance[targetList[i]] = abundanceVector[i];
//     }

//     for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             it->second[j].abundance = abundance[it->second[j].result.dbKey];
//         }
//     }
// }

// static std::vector<double> emUpdate(MappingTable &mappingTable,
//                                     double lambda,
//                                     const std::unordered_map<unsigned int, double> &fixedCoverageConfidence,
//                                     const std::vector<unsigned int> &targetList,
//                                     size_t queryCount,
//                                     double alpha,
//                                     double gamma) {
//     computePosterior(mappingTable, lambda, alpha, gamma);

//     std::unordered_map<unsigned int, double> nextAbundance;
//     nextAbundance.reserve(targetList.size());
//     for (size_t i = 0; i < targetList.size(); ++i) {
//         nextAbundance[targetList[i]] = 0.0;
//     }

//     for (MappingTable::const_iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             nextAbundance[it->second[j].result.dbKey] += it->second[j].posterior;
//         }
//     }

//     if (queryCount > 0) {
//         const double denom = static_cast<double>(queryCount);
//         for (std::unordered_map<unsigned int, double>::iterator it = nextAbundance.begin(); it != nextAbundance.end(); ++it) {
//             it->second /= denom;
//         }
//     }

//     for (MappingTable::iterator it = mappingTable.begin(); it != mappingTable.end(); ++it) {
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             const unsigned int target = it->second[j].result.dbKey;
//             it->second[j].abundance = nextAbundance[target];
//             std::unordered_map<unsigned int, double>::const_iterator fixed = fixedCoverageConfidence.find(target);
//             if (fixed != fixedCoverageConfidence.end()) {
//                 it->second[j].coverageConfidence = fixed->second;
//             } else {
//                 it->second[j].coverageConfidence = 0.0;
//             }
//         }
//     }

//     std::vector<double> out;
//     out.reserve(targetList.size());
//     for (size_t i = 0; i < targetList.size(); ++i) {
//         out.push_back(nextAbundance[targetList[i]]);
//     }
//     return out;
// }

// static std::vector<double> projectSimplex(const std::vector<double> &x) {
//     std::vector<double> projected = x;
//     for (size_t i = 0; i < projected.size(); ++i) {
//         if (projected[i] < 0.0) {
//             projected[i] = 0.0;
//         }
//     }

//     const double sum = std::accumulate(projected.begin(), projected.end(), 0.0);
//     if (sum > 0.0) {
//         for (size_t i = 0; i < projected.size(); ++i) {
//             projected[i] /= sum;
//         }
//     }
//     return projected;
// }

// static void squarem(ReclassTaxContext &ctx,
//                     double lambda,
//                     int maxIter,
//                     double tol,
//                     double alpha,
//                     double gamma,
//                     int threads) {
//     if (ctx.queryCount == 0 || ctx.targetSet.empty()) {
//         return;
//     }

//     initAbundance(ctx.mappingTable, ctx.targetSet, ctx.queryCount);
//     initCoverageConfidence(ctx.mappingTable, ctx.targetSet, threads);
//     Debug(Debug::INFO) << "Reclassify initialized coverage confidence." << "\n";

//     std::unordered_map<unsigned int, double> fixedCoverageConfidence;
//     fixedCoverageConfidence.reserve(ctx.targetSet.size());
//     for (MappingTable::const_iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end(); ++it) {
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             fixedCoverageConfidence[it->second[j].result.dbKey] = it->second[j].coverageConfidence;
//         }
//     }

//     const std::vector<unsigned int> targetList = targetListFromSet(ctx.targetSet);
//     std::vector<double> x0 = abundanceVectorFromTable(ctx.mappingTable, targetList);
//     std::vector<double> logLikelihoods;

//     for (int iter = 0; iter < maxIter; ++iter) {
//         const std::vector<double> x1 = emUpdate(ctx.mappingTable, lambda, fixedCoverageConfidence, targetList, ctx.queryCount, alpha, gamma);
//         const std::vector<double> x2 = emUpdate(ctx.mappingTable, lambda, fixedCoverageConfidence, targetList, ctx.queryCount, alpha, gamma);

//         std::vector<double> r(x0.size(), 0.0);
//         std::vector<double> v(x0.size(), 0.0);
//         for (size_t i = 0; i < x0.size(); ++i) {
//             r[i] = x1[i] - x0[i];
//             v[i] = x2[i] - x1[i] - r[i];
//         }

//         double normR = 0.0;
//         double normV = 0.0;
//         for (size_t i = 0; i < x0.size(); ++i) {
//             normR += r[i] * r[i];
//             normV += v[i] * v[i];
//         }
//         normR = std::sqrt(normR);
//         normV = std::sqrt(normV);

//         double accel = (normV == 0.0) ? -1.0 : -(normR / normV);
//         accel = std::max(STEP_MIN, std::min(STEP_MAX, accel));

//         std::vector<double> xNew(x0.size(), 0.0);
//         for (size_t i = 0; i < x0.size(); ++i) {
//             xNew[i] = x0[i] - 2.0 * accel * r[i] + accel * accel * v[i];
//         }
//         xNew = projectSimplex(xNew);

//         setAbundance(ctx.mappingTable, targetList, xNew);
//         computePosterior(ctx.mappingTable, lambda, alpha, gamma);
//         double currentLl = logLikelihood(ctx.mappingTable, lambda, alpha, gamma, ctx.queryCount);

//         if (!logLikelihoods.empty() && currentLl < logLikelihoods.back() - 1e-9) {
//             setAbundance(ctx.mappingTable, targetList, x2);
//             computePosterior(ctx.mappingTable, lambda, alpha, gamma);
//             currentLl = logLikelihood(ctx.mappingTable, lambda, alpha, gamma, ctx.queryCount);
//             xNew = x2;
//         }

//         logLikelihoods.push_back(currentLl);

//         double parameterChange = 0.0;
//         for (size_t i = 0; i < x0.size(); ++i) {
//             parameterChange = std::max(parameterChange, std::fabs(xNew[i] - x0[i]));
//         }

//         Debug(Debug::INFO) << "Reclassify iteration " << iter << ": LL=" << currentLl << " delta=" << parameterChange << "\n";
//         x0 = xNew;
//         if (parameterChange < tol && iter > 5) {
//             Debug(Debug::INFO) << "Reclassify converged after " << (iter + 1) << " iterations." << "\n";
//             break;
//         }
//     }
// }

// static void addInterval(std::vector<Interval> &intervals, int start, int end) {
//     Interval interval;
//     interval.start = std::min(start, end);
//     interval.end = std::max(start, end);
//     intervals.push_back(interval);
// }

// static std::vector<Interval> mergeIntervals(std::vector<Interval> intervals) {
//     if (intervals.empty()) {
//         return intervals;
//     }

//     std::sort(intervals.begin(), intervals.end(), [](const Interval &lhs, const Interval &rhs) {
//         if (lhs.start != rhs.start) {
//             return lhs.start < rhs.start;
//         }
//         return lhs.end < rhs.end;
//     });

//     std::vector<Interval> merged;
//     merged.push_back(intervals[0]);
//     for (size_t i = 1; i < intervals.size(); ++i) {
//         if (intervals[i].start <= merged.back().end + 1) {
//             merged.back().end = std::max(merged.back().end, intervals[i].end);
//         } else {
//             merged.push_back(intervals[i]);
//         }
//     }
//     return merged;
// }

// static std::vector<TargetStats> collectTargetStats(const ReclassTaxContext &ctx) {
//     std::unordered_map<unsigned int, TargetStats> statsByTarget;

//     for (MappingTable::const_iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end(); ++it) {
//         for (size_t j = 0; j < it->second.size(); ++j) {
//             const ReclassTaxEntry &entry = it->second[j];
//             TargetStats &stats = statsByTarget[entry.result.dbKey];
//             stats.key = entry.result.dbKey;
//             stats.targetLength = entry.result.dbLen;
//             stats.abundance = entry.abundance;
//             stats.coverageConfidence = entry.coverageConfidence;
//             stats.dropped = false;
//             addInterval(stats.intervals, entry.result.dbStartPos, entry.result.dbEndPos);
//         }
//     }

//     std::vector<TargetStats> out;
//     out.reserve(statsByTarget.size());
//     for (std::unordered_map<unsigned int, TargetStats>::iterator it = statsByTarget.begin(); it != statsByTarget.end(); ++it) {
//         it->second.intervals = mergeIntervals(it->second.intervals);
//         out.push_back(it->second);
//     }

//     std::sort(out.begin(), out.end(), [](const TargetStats &lhs, const TargetStats &rhs) {
//         if (lhs.abundance != rhs.abundance) {
//             return lhs.abundance > rhs.abundance;
//         }
//         return lhs.key < rhs.key;
//     });
//     return out;
// }

// static double clamp01(double value) {
//     return std::max(0.0, std::min(1.0, value));
// }

// static bool largestJumpCutoff(std::vector<double> values,
//                               bool useLowTail,
//                               double maxTailFraction,
//                               double &cutoff,
//                               size_t &tailCount) {
//     cutoff = 0.0;
//     tailCount = 0;
//     if (values.size() < MIN_FILTER_TARGETS) {
//         return false;
//     }

//     std::sort(values.begin(), values.end());
//     maxTailFraction = clamp01(maxTailFraction);
//     const double totalMass = std::accumulate(values.begin(), values.end(), 0.0);
//     if (totalMass <= EPS || maxTailFraction <= 0.0) {
//         return false;
//     }
//     const double maxTailMass = maxTailFraction * totalMass;

//     double bestGap = 0.0;
//     size_t bestIdx = 0;
//     double lowTailMass = 0.0;
//     for (size_t i = 0; i + 1 < values.size(); ++i) {
//         const size_t lowTailCount = i + 1;
//         const size_t highTailCount = values.size() - lowTailCount;
//         const size_t candidateTailCount = useLowTail ? lowTailCount : highTailCount;
//         lowTailMass += values[i];
//         const double highTailMass = totalMass - lowTailMass;
//         const double candidateTailMass = useLowTail ? lowTailMass : highTailMass;
//         if (candidateTailCount < MIN_TAIL_TARGETS || candidateTailMass > (maxTailMass + EPS)) {
//             continue;
//         }

//         const double gap = values[i + 1] - values[i];
//         if (gap > bestGap) {
//             bestGap = gap;
//             bestIdx = i;
//             tailCount = candidateTailCount;
//         }
//     }

//     if (bestGap <= EPS) {
//         return false;
//     }

//     cutoff = 0.5 * (values[bestIdx] + values[bestIdx + 1]);
//     return true;
// }

// static bool tailQuantileCutoff(std::vector<double> values,
//                                bool useLowTail,
//                                double maxTailFraction,
//                                double &cutoff,
//                                size_t &tailCount) {
//     cutoff = 0.0;
//     tailCount = 0;
//     if (values.size() < MIN_FILTER_TARGETS) {
//         return false;
//     }

//     std::sort(values.begin(), values.end());
//     maxTailFraction = clamp01(maxTailFraction);
//     const double totalMass = std::accumulate(values.begin(), values.end(), 0.0);
//     if (totalMass <= EPS || maxTailFraction <= 0.0) {
//         return false;
//     }
//     const double maxTailMass = maxTailFraction * totalMass;

//     double accumulatedMass = 0.0;
//     size_t maxTailCount = 0;
//     for (size_t i = 0; i < values.size(); ++i) {
//         const double candidate = accumulatedMass + values[i];
//         if (candidate > (maxTailMass + EPS)) {
//             break;
//         }
//         accumulatedMass = candidate;
//         ++maxTailCount;
//     }
//     if (maxTailCount < MIN_TAIL_TARGETS || maxTailCount >= values.size()) {
//         return false;
//     }

//     tailCount = maxTailCount;
//     if (useLowTail) {
//         cutoff = values[tailCount - 1];
//     } else {
//         cutoff = values[values.size() - tailCount];
//     }
//     return true;
// }

// static std::unordered_set<unsigned int> selectTailTargets(const std::vector<TargetStats> &stats,
//                                                           bool useLowTail,
//                                                           size_t tailCount,
//                                                           double maxTailFraction) {
//     std::vector<const TargetStats *> ordered;
//     ordered.reserve(stats.size());
//     for (size_t i = 0; i < stats.size(); ++i) {
//         ordered.push_back(&stats[i]);
//     }

//     std::sort(ordered.begin(), ordered.end(), [useLowTail](const TargetStats *lhs, const TargetStats *rhs) {
//         const double lhsValue = useLowTail ? lhs->abundance : lhs->coverageConfidence;
//         const double rhsValue = useLowTail ? rhs->abundance : rhs->coverageConfidence;
//         if (lhsValue != rhsValue) {
//             return useLowTail ? (lhsValue < rhsValue) : (lhsValue > rhsValue);
//         }
//         return lhs->key < rhs->key;
//     });

//     double totalMass = 0.0;
//     for (size_t i = 0; i < ordered.size(); ++i) {
//         const double value = useLowTail ? ordered[i]->abundance : ordered[i]->coverageConfidence;
//         totalMass += value;
//     }
//     const double maxTailMass = clamp01(maxTailFraction) * totalMass;

//     std::unordered_set<unsigned int> selected;
//     const size_t limit = std::min(tailCount, ordered.size());
//     double selectedMass = 0.0;
//     selected.reserve(limit);
//     for (size_t i = 0; i < limit; ++i) {
//         const double value = useLowTail ? ordered[i]->abundance : ordered[i]->coverageConfidence;
//         if (selected.size() >= MIN_TAIL_TARGETS && (selectedMass + value) > (maxTailMass + EPS)) {
//             break;
//         }
//         selectedMass += value;
//         selected.insert(ordered[i]->key);
//     }
//     return selected;
// }

// static std::unordered_set<unsigned int> selectDroppedTargets(const std::vector<TargetStats> &stats,
//                                                              double maxDropPercentage,
//                                                              double &abundanceCutoff) {
//     std::unordered_set<unsigned int> dropped;
//     if (stats.empty()) {
//         abundanceCutoff = 0.0;
//         return dropped;
//     }
//     if (stats.size() < MIN_FILTER_TARGETS) {
//         abundanceCutoff = 0.0;
//         return dropped;
//     }

//     std::vector<double> abundances;
//     abundances.reserve(stats.size());
//     for (size_t i = 0; i < stats.size(); ++i) {
//         abundances.push_back(stats[i].abundance);
//     }

//     const double maxTailFraction = clamp01(maxDropPercentage / 100.0);
//     size_t abundanceTailCount = 0;
//     bool hasAbundanceCutoff = largestJumpCutoff(abundances, true, maxTailFraction, abundanceCutoff, abundanceTailCount);
//     if (hasAbundanceCutoff == false) {
//         hasAbundanceCutoff = tailQuantileCutoff(abundances, true, maxTailFraction, abundanceCutoff, abundanceTailCount);
//     }
//     if (hasAbundanceCutoff == false) {
//         abundanceCutoff = 0.0;
//         return dropped;
//     }

//     const std::unordered_set<unsigned int> lowAbundanceTargets = selectTailTargets(stats, true, abundanceTailCount, maxTailFraction);
//     for (std::unordered_set<unsigned int>::const_iterator it = lowAbundanceTargets.begin(); it != lowAbundanceTargets.end(); ++it) {
//         dropped.insert(*it);
//     }
//     if (dropped.size() == stats.size()) {
//         dropped.clear();
//     }
//     return dropped;
// }

// static void applyDroppedTargets(ReclassTaxContext &ctx,
//                                 const std::unordered_set<unsigned int> &dropped,
//                                 size_t totalTargets,
//                                 double abundanceCutoff) {
//     if (dropped.empty()) {
//         Debug(Debug::INFO) << "Reclassify target filter kept all targets. abundance cutoff="
//                            << abundanceCutoff << "\n";
//         return;
//     }

//     for (MappingTable::iterator it = ctx.mappingTable.begin(); it != ctx.mappingTable.end();) {
//         std::vector<ReclassTaxEntry> &records = it->second;
//         records.erase(std::remove_if(records.begin(), records.end(), [&dropped](const ReclassTaxEntry &entry) {
//             return dropped.find(entry.result.dbKey) != dropped.end();
//         }), records.end());

//         if (records.empty()) {
//             it = ctx.mappingTable.erase(it);
//         } else {
//             ++it;
//         }
//     }

//     ctx.queryOrder.erase(std::remove_if(ctx.queryOrder.begin(), ctx.queryOrder.end(), [&ctx](unsigned int queryKey) {
//         return ctx.mappingTable.find(queryKey) == ctx.mappingTable.end();
//     }), ctx.queryOrder.end());

//     for (std::unordered_set<unsigned int>::const_iterator it = dropped.begin(); it != dropped.end(); ++it) {
//         ctx.targetSet.erase(*it);
//     }

//     const double removedPct = (totalTargets > 0)
//                                 ? (100.0 * static_cast<double>(dropped.size()) / static_cast<double>(totalTargets))
//                                 : 0.0;
//     Debug(Debug::INFO) << "Reclassify dropped " << dropped.size()
//                        << " of " << totalTargets
//                        << " targets (" << removedPct << "%)"
//                        << " using abundance <= " << abundanceCutoff << ".\n";
// }

// static bool compareByPosteriorThenBitScore(const ReclassTaxEntry &a, const ReclassTaxEntry &b) {
//     if (a.posterior != b.posterior) {
//         return a.posterior > b.posterior;
//     }
//     if (a.result.score != b.result.score) {
//         return a.result.score > b.result.score;
//     }
//     return Matcher::compareHits(a.result, b.result);
// }

// static void writeReclassifiedDb(const ReclassTaxContext &ctx,
//                                 int dbType,
//                                 const std::string &outDb,
//                                 const std::string &outIndex,
//                                 int threads,
//                                 bool compress) {
//     DBWriter writer(outDb.c_str(), outIndex.c_str(), threads, compress, dbType);
//     writer.open();

//     Debug::Progress progress(ctx.queryOrder.size());
// #pragma omp parallel
//     {
//         unsigned int thread_idx = 0;
// #ifdef OPENMP
//         thread_idx = static_cast<unsigned int>(omp_get_thread_num());
// #endif
//         char buffer[1024 + 32768 * 4];

// #pragma omp for schedule(dynamic, 5)
//         for (size_t i = 0; i < ctx.queryOrder.size(); ++i) {
//             progress.updateProgress();
//             const unsigned int queryKey = ctx.queryOrder[i];
//             MappingTable::const_iterator recordsIt = ctx.mappingTable.find(queryKey);
//             if (recordsIt == ctx.mappingTable.end()) {
//                 continue;
//             }

//             std::vector<ReclassTaxEntry> records = recordsIt->second;
//             SORT_SERIAL(records.begin(), records.end(), compareByPosteriorThenBitScore);

//             writer.writeStart(thread_idx);
//             for (size_t j = 0; j < records.size(); ++j) {
//                 Matcher::result_t res = records[j].result;
//                 res.seqId = static_cast<float>(records[j].posterior);
//                 size_t len = Matcher::resultToBuffer(buffer, res, ctx.hasBacktrace, ctx.hasOrfPosition);
//                 writer.writeAdd(buffer, len, thread_idx);
//             }
//             writer.writeEnd(queryKey, thread_idx);
//         }
//     }

//     writer.close();
// }
// }

// int emreclassify(int argc, const char **argv, const Command &command) {
//     Parameters &par = Parameters::getInstance();
//     par.parseParameters(argc, argv, command, true, 0, 0);

//     DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
//                                   DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
//     reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

//     ReclassTaxContext ctx;
//     loadAlignmentDb(reader, ctx);
//     Debug(Debug::INFO) << "Loaded " << ctx.queryCount << " queries with hits and " << ctx.targetSet.size() << " unique targets.\n";

//     squarem(ctx,
//             par.reclassifyLambda,
//             par.reclassifyMaxIterations,
//             par.reclassifyTolerance,
//             par.reclassifyAlpha,
//             par.reclassifyGamma,
//             par.threads);

//     std::vector<TargetStats> allTargetStats = collectTargetStats(ctx);
//     double abundanceCutoff = 0.0;
//     const std::unordered_set<unsigned int> dropped = selectDroppedTargets(allTargetStats,
//                                                                           par.reclassifyMaxDropPercentage,
//                                                                           abundanceCutoff);
//     applyDroppedTargets(ctx, dropped, allTargetStats.size(), abundanceCutoff);

//     writeReclassifiedDb(ctx, reader.getDbtype(), par.db4, par.db4Index, par.threads, par.compressed);

//     reader.close();
//     return EXIT_SUCCESS;
// }
