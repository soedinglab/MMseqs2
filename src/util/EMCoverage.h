#ifndef MMSEQS_EMCOVERAGE_H
#define MMSEQS_EMCOVERAGE_H

// Numeric coverage kernels shared by `mmseqs reclassify` and `mmseqs abundance`.
//
// directly and checks every formula below against hand-computed values. Keep it that way --
// these are the only places where a silent arithmetic change turns into a silently wrong
// abundance table.

#include <algorithm>
#include <string>
#include <vector>

namespace EMCoverage {

inline double clampUnit(double value) {
    return std::max(0.0, std::min(1.0, value));
}


inline int strandStep(int dbStartPos, int dbEndPos) {
    return (dbStartPos > dbEndPos) ? -1 : 1;
}

// How many target positions one backtrace column consumes.
//
// It is 1 whenever the backtrace and the target coordinates use the same alphabet (protein vs
// protein, nucleotide vs nucleotide), but 3 after `offsetalignment` has mapped a translated search
// back onto the parent nucleotide sequence: there the backtrace still counts amino-acid columns
// while dbStartPos/dbEndPos/dbLen have become nucleotide offsets on the contig. A real record from
// `search --search-type 2` against a nucleotide target DB looks like
//     dbStartPos=902  dbEndPos=303  dbLen=1203  backtrace="200M"
// -- 200 columns spanning 600 target positions. Walking one position per column would paint a third
// of the gene and leave the rest reading as uncovered.
//
// The stride is derived from the alignment rather than from a search-type flag, which the
// alignment DB does not carry: the target-consuming columns (M and D) have to span exactly
// dbEndPos - dbStartPos + 1 positions. Anything that is not a clean 1 or 3 falls back to 1, so a
// malformed or prefilter-only record degrades to the old behaviour instead of scattering coverage.
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

// Paint one hit's posterior-weighted coverage onto the two per-position accumulators:
//   notCov[p] *= (1 - gamma)   so 1 - notCov[p] is P[position covered by >= 1 read]
//   depth[p]  += gamma         so depth[p] is the expected read depth at p
// The backtrace is a run-length-encoded CIGAR ("153M2D10M"; a missing count means 1). Only 'M'
// columns place a read residue on a target position, so only they contribute coverage. 'M' and
// 'D' both consume a target residue and advance the target position; 'I' consumes the query only.
// Positions outside [0, L) are ignored, so a hit can never contribute coverage beyond the target.
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

// Fallback when the alignment DB carries no backtrace (search run without -a): treat the whole
// aligned span as covered, clipped to [0, L). Strand-safe by construction via min/max, and a
// prefilter-only record (dbStartPos == dbEndPos == -1) paints nothing.
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

// Span painting for the reclassify coverage prior, which weights by the within-query bit-score
// share rather than by an EM posterior and needs the raw read depth separately (see evenness()).
inline void paintSpanWeighted(std::vector<double> &weightedDepth,
                              std::vector<double> &readDepth,
                              int dbStartPos,
                              int dbEndPos,
                              double weight) {
    const int targetLength = static_cast<int>(weightedDepth.size());
    const int start = std::max(0, std::min(dbStartPos, dbEndPos));
    const int end = std::min(targetLength - 1, std::max(dbStartPos, dbEndPos));
    for (int pos = start; pos <= end; ++pos) {
        weightedDepth[static_cast<size_t>(pos)] += weight;
        readDepth[static_cast<size_t>(pos)] += 1.0;
    }
}

// Breadth F = (1/L) * sum_p min(1, d_p), in [0, 1].
// The denominator is the full target length, never the observed hit span: with the span, a target
// touched by a single read had span == that hit's length, so F ~= 1 and a one-read target scored
// as "broadly covered".
inline double breadth(const std::vector<double> &weightedDepth) {
    if (weightedDepth.empty()) {
        return 0.0;
    }
    double covered = 0.0;
    for (size_t p = 0; p < weightedDepth.size(); ++p) {
        covered += std::min(1.0, weightedDepth[p]);
    }
    return covered / static_cast<double>(weightedDepth.size());
}

// Evenness E = (sum_p d_p)^2 / (L * sum_p d_p^2), summed over the WHOLE target length.
// In [0, 1]: 1 iff the depth is uniform across the full length, and ~f when the depth is uniform
// over a fraction f of it and zero elsewhere, so a concentrated pileup scores low.
//
// The denominator is L, NOT the number of covered positions n_cov. With n_cov the term did the
// opposite of what it was written for: a pileup confined to a narrow window is "perfectly even
// over its own footprint" and scored E = 1, so it received no repeat penalty at all. Measured on
// a 500-residue target, C_t = F * E came out as
//     full depth-1 tiling                 F=1.00 E=1.000 -> C=1.000
//     full tiling + one 100x repeat spike F=1.00 E=0.068 -> C=0.068
//     100x pileup on a 50-residue repeat  F=0.10 E=1.000 -> C=0.100   <- artefact ranked highest
// i.e. the pure artefact outranked the genuinely covered target. Dividing by L instead gives
// 1.000 / 0.068 / 0.010 and restores the intended order.
//
// Known consequence: since E == f for a depth uniform over a fraction f, E is no longer
// independent of breadth and C = F * E is roughly quadratic in the covered fraction. That is
// unavoidable rather than a defect -- "concentrated" *means* "the covered region is small", so any
// term that can detect it has to see the uncovered positions too. Only the ratio C_t / sum(C)
// enters the prior, so this changes the prior's sharpness, not its scale.
//
// Feed this the RAW read depth, not the score-share-weighted depth. Weighted, E reacted to hit
// ambiguity rather than to depth raggedness: two reads at disjoint positions (depth is exactly one
// read everywhere either way) scored E = 1.000 at weights 1.0/1.0 but E = 0.599 at 1.0/0.1, while
// two equally vague reads at 0.1/0.1 scored E = 1.000 again. Ambiguity is already priced into F.
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
