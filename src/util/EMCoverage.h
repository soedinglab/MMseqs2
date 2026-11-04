#ifndef MMSEQS_EMCOVERAGE_H
#define MMSEQS_EMCOVERAGE_H


#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace EMCoverage {

inline double clampUnit(double value) {
    return std::max(0.0, std::min(1.0, value));
}


inline int strandStep(int dbStartPos, int dbEndPos) {
    return (dbStartPos > dbEndPos) ? -1 : 1;
}

inline int targetStride(const std::string &compressedBacktrace, int dbStartPos, int dbEndPos) {
    if (dbStartPos < 0 || dbEndPos < 0) {
        return 1;
    }
    long targetColumns = 0;
    size_t count = 0;
    for (size_t i = 0; i < compressedBacktrace.size(); ++i) {
        const char c = compressedBacktrace[i];
        if (c >= '0' && c <= '9') {
            count = count * 10 + static_cast<size_t>(c - '0');
            continue;
        }
        const long n = (count == 0) ? 1 : static_cast<long>(count);
        count = 0;
        if (c == 'M' || c == 'D') {
            targetColumns += n;
        }
    }
    if (targetColumns <= 0) {
        return 1;
    }
    const long span = std::max(dbStartPos, dbEndPos) - std::min(dbStartPos, dbEndPos) + 1;
    if (span == targetColumns * 3) {
        return 3;
    }
    return 1;
}

inline void paintBacktrace(std::vector<double> &notCov,
                           std::vector<double> &depth,
                           const std::string &compressedBacktrace,
                           int dbStartPos,
                           int dbEndPos,
                           double gamma) {
    const int targetLength = static_cast<int>(notCov.size());
    const double keep = 1.0 - gamma;
    const int step = strandStep(dbStartPos, dbEndPos);
    const int stride = targetStride(compressedBacktrace, dbStartPos, dbEndPos);
    int targetPos = dbStartPos;
    size_t count = 0;
    for (size_t i = 0; i < compressedBacktrace.size(); ++i) {
        const char c = compressedBacktrace[i];
        if (c >= '0' && c <= '9') {
            count = count * 10 + static_cast<size_t>(c - '0');
            continue;
        }
        const int n = (count == 0) ? 1 : static_cast<int>(count);
        count = 0;
        if (c == 'M') {
            for (int k = 0; k < n * stride; ++k, targetPos += step) {
                if (targetPos >= 0 && targetPos < targetLength) {
                    notCov[static_cast<size_t>(targetPos)] *= keep;
                    depth[static_cast<size_t>(targetPos)] += gamma;
                }
            }
        } else if (c == 'D') {
            targetPos += n * stride * step;
        }
    }
}

inline void paintSpan(std::vector<double> &notCov,
                      std::vector<double> &depth,
                      int dbStartPos,
                      int dbEndPos,
                      double gamma) {
    const int targetLength = static_cast<int>(notCov.size());
    const double keep = 1.0 - gamma;
    const int start = std::max(0, std::min(dbStartPos, dbEndPos));
    const int end = std::min(targetLength - 1, std::max(dbStartPos, dbEndPos));
    for (int pos = start; pos <= end; ++pos) {
        notCov[static_cast<size_t>(pos)] *= keep;
        depth[static_cast<size_t>(pos)] += gamma;
    }
}

inline double evenness(const std::vector<double> &depth) {
    double sum = 0.0;
    double sumSquared = 0.0;
    for (size_t p = 0; p < depth.size(); ++p) {
        sum += depth[p];
        sumSquared += depth[p] * depth[p];
    }
    const double denominator = static_cast<double>(depth.size()) * sumSquared;
    return (denominator > 0.0) ? clampUnit((sum * sum) / denominator) : 0.0;
}

// Expected number of covered residues = sum_p P[p covered by >= 1 read] = sum_p (1 - notCov[p]).
inline double expectedCoveredLength(const std::vector<double> &notCov) {
    double covered = 0.0;
    for (size_t p = 0; p < notCov.size(); ++p) {
        covered += (1.0 - notCov[p]);
    }
    return covered;
}

// Mean expected read depth over the target = (1/L) * sum_p depth[p].
inline double meanDepth(const std::vector<double> &depth) {
    if (depth.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (size_t p = 0; p < depth.size(); ++p) {
        sum += depth[p];
    }
    return sum / static_cast<double>(depth.size());
}

}

#endif
