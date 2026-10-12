#ifndef FWBW_SCALED_H
#define FWBW_SCALED_H

// Scaled-mantissa forward-backward for FwBwAligner::runFwBw<false, 0> (no profile), used instead of the
// log-space blocks when FwBwAligner::scaled is set.
//
// Same model as FwBwAligner:
//   ZM_ij = (ZM_{i-1,j-1} + ZE_{i-1,j-1} + ZF_{i-1,j-1} + 1) * exp(S_ij)      (S already divided by T)
//   ZE_ij = ZM_{i,j-1} * eGo + ZE_{i,j-1} * eGe
//   ZF_ij = ZM_{i-1,j} * eGo + ZF_{i-1,j} * eGe
//   Z = sum ZM,  P_ij = ZM_ij * ZMb_ij / (exp(S_ij) * Z) = ZM_ij * qb_ij / Z
// where ZMb is the recurrence on the reversed matrix and qb_ij = ZMb_ij / exp(S_ij) its
// "predecessor + 1" term, so the posterior needs no division and no second exp.
//
// Every value is a float mantissa times 2^s with one integer scale s >= 0 per (row, 16-column block).
// The local-start term 1 (= 2^-s in block units) bounds the values from below, so scales never go
// negative. A block is renormalised (exponent arithmetic only) when its values leave [2^-32, 2^32]; that
// test is a branch that is almost never taken, so it stays off the row-to-row dependency chain. One exp
// per cell (computed in the forward pass and kept for the backward pass), no log, no division.
//
// Non-finite results (exp saturation, NaN, Z = 0, P > 1 + 1e-3) are returned as flags; the caller then
// runs the log-space code, which keeps its overflow semantics (maxP = 0).

#include "simd.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

class FwBwScaled {
public:
    // return flags of run()
    static const unsigned int FLAG_NON_FINITE = 1;
    static const unsigned int FLAG_Z = 2;
    static const unsigned int FLAG_P = 4;

    // buffers reused across calls
    struct Workspace {
        std::vector<float> pM, pE, pF; // previous row, mantissas [colsPadded]
        std::vector<int> ps;           // previous row, scale per block [blocks]
        std::vector<int> fs;           // forward ZM scale per (row, block) [rows * blocks]
        std::vector<float> eS;         // exp(S) of the forward pass, reused by the backward pass [rows * colsPadded]
        std::vector<float> acc;        // partial Z per block [blocks * VECSIZE_FLOAT]
        std::vector<int> accS;         // its scale [blocks]
        double zv = 0;                 // Z = zv * 2^zs
        int zs = 0;
    };

    // score: rows row pointers of cols scores already divided by T (FwBwAligner::scoreForward)
    // prob: rows row pointers (>= cols floats each), receives the posteriors
    // returns 0 on success, else FLAG_* (prob and maxP are then invalid)
    static unsigned int run(float* const* score, int rows, int cols, float expGo, float expGe,
                            float* const* prob, float& maxP, Workspace& ws) {
        FlushDenormals ftz;
        const int blocks = (cols + BLOCK_COLS - 1) / BLOCK_COLS;
        const int colsPadded = blocks * BLOCK_COLS;
        if (ws.pM.size() < (size_t)colsPadded) {
            ws.pM.resize(colsPadded);
            ws.pE.resize(colsPadded);
            ws.pF.resize(colsPadded);
        }
        if (ws.ps.size() < (size_t)blocks) {
            ws.ps.resize(blocks);
            ws.accS.resize(blocks);
            ws.acc.resize((size_t)blocks * VECSIZE_FLOAT);
        }
        if (ws.fs.size() < (size_t)rows * blocks) {
            ws.fs.resize((size_t)rows * blocks);
        }
        if (ws.eS.size() < (size_t)rows * colsPadded) {
            ws.eS.resize((size_t)rows * colsPadded);
        }
        float unused = 0.0f;
        unsigned int flags = pass<false>(score, rows, cols, expGo, expGe, prob, ws, 0, 1.0f, unused);
        if (flags != 0) {
            return flags;
        }
        // Z = m * 2^(zs + e), m in [0.5, 1)
        int e = 0;
        const double m = std::frexp(ws.zv, &e);
        return pass<true>(score, rows, cols, expGo, expGe, prob, ws, ws.zs + e, (float)(1.0 / m), maxP);
    }

private:
    // columns per block, one scale per (row, block)
    static const int BLOCK_COLS = 16;
    static const int BLOCK_VECS = BLOCK_COLS / VECSIZE_FLOAT;

    // flush-to-zero / denormals-are-zero for the duration of a call (x86 only): values decaying below
    // 2^-126 of their scale (gap chains in low-scoring regions) are negligible next to the local-start
    // term, but denormal arithmetic costs 10-27% of the run time
    struct FlushDenormals {
#if defined(SIMDE_X86_SSE_NATIVE)
        unsigned int csr;
        FlushDenormals() : csr(_mm_getcsr()) { _mm_setcsr(csr | 0x8040u); }
        ~FlushDenormals() { _mm_setcsr(csr); }
#else
        // user-provided, so that an unused instance does not trigger -Wunused-variable
        FlushDenormals() {}
#endif
    };

#ifdef AVX2
    static inline float simdf32_first(simd_float a) {
        return _mm256_cvtss_f32(a);
    }

    // last lane in all lanes
    static inline simd_float simdf32_bcast_last(simd_float a) {
        return _mm256_permutevar8x32_ps(a, _mm256_set1_epi32(7));
    }
#else
    static inline float simdf32_first(simd_float a) {
        return _mm_cvtss_f32(a);
    }

    // last lane in all lanes
    static inline simd_float simdf32_bcast_last(simd_float a) {
        return _mm_shuffle_ps(a, a, 0xFF);
    }
#endif

    static inline simd_float simdf32_isnan(simd_float a) {
        return simdf32_cmp(a, a, 3);
    }

    // lanes shifted towards higher indices by one, zero shifted in
    static inline simd_float simdf32_shiftl1(simd_float a) {
        return simdi_i2fcast(simdi8_shiftl(simdf_f2icast(a), 4));
    }

    // y_t = sum_{k<t} x_k g^{t-1-k}: exclusive prefix sum decayed by g (cf. simdf32_prefixsum in Fwbw.cpp)
    static inline simd_float exclusiveScan(simd_float x, simd_float g1, simd_float g2, simd_float g4) {
        simd_float y = simdf32_shiftl1(x);
        y = simdf32_fmadd(simdi_i2fcast(simdi8_shiftl(simdf_f2icast(y), 4)), g1, y);
        y = simdf32_fmadd(simdi_i2fcast(simdi8_shiftl(simdf_f2icast(y), 8)), g2, y);
#ifdef AVX2
        y = simdf32_fmadd(simdi_i2fcast(simdi8_shiftl(simdf_f2icast(y), 16)), g4, y);
#else
        (void)g4;
#endif
        return y;
    }

    // 2^e for integer e, 0 below 2^-126, clamped at 2^127
    static inline float pow2f(int e) {
        e = std::max(-127, std::min(e, 127));
        const uint32_t b = (uint32_t)(e + 127) << 23;
        float f;
        std::memcpy(&f, &b, sizeof(f));
        return f;
    }

    // binary exponent of a normal float
    static inline int fexp(float x) {
        uint32_t b;
        std::memcpy(&b, &x, sizeof(b));
        return (int)((b >> 23) & 0xff) - 127;
    }

    static inline double hsum(simd_float a) {
        float t[VECSIZE_FLOAT];
        simdf32_storeu(t, a);
        double s = 0;
        for (int k = 0; k < VECSIZE_FLOAT; ++k) {
            s += t[k];
        }
        return s;
    }

    // row[p, p + VECSIZE_FLOAT), lanes at or after cols set to pad
    static inline simd_float loadCols(const float* row, int p, int cols, float pad) {
        if (p + VECSIZE_FLOAT <= cols) {
            return simdf32_loadu(row + p);
        }
        float t[VECSIZE_FLOAT];
        for (int k = 0; k < VECSIZE_FLOAT; ++k) {
            t[k] = (p + k < cols) ? row[p + k] : pad;
        }
        return simdf32_loadu(t);
    }

    // stores the lanes of a before cols to row[p, p + VECSIZE_FLOAT)
    static inline void storeCols(float* row, int p, int cols, simd_float a) {
        if (p + VECSIZE_FLOAT <= cols) {
            simdf32_storeu(row + p, a);
            return;
        }
        float t[VECSIZE_FLOAT];
        simdf32_storeu(t, a);
        for (int k = 0; k < VECSIZE_FLOAT && p + k < cols; ++k) {
            row[p + k] = t[k];
        }
    }

    // Z += v * 2^s
    static inline void zAdd(Workspace& ws, double v, int s) {
        if (v == 0) {
            return;
        }
        if (ws.zv == 0) {
            ws.zv = v;
            ws.zs = s;
        } else if (s > ws.zs) {
            ws.zv = std::ldexp(ws.zv, ws.zs - s) + v;
            ws.zs = s;
        } else {
            ws.zv += std::ldexp(v, s - ws.zs);
        }
    }

    // One pass. Forward: rows 0..rows-1, columns 0..colsPadded-1; stores the ZM mantissas in prob and their
    // scales in ws.fs, accumulates Z. Backward: the same recurrence on the reversed matrix (row rows-1-i,
    // column colsPadded-1-j, so backward block b covers forward block blocks-1-b); overwrites prob with P.
    template<bool backward>
    static unsigned int pass(float* const* score, int rows, int cols, float expGo, float expGe,
                             float* const* prob, Workspace& ws, int zExp, float zInvMant, float& maxP) {
        const int blocks = (cols + BLOCK_COLS - 1) / BLOCK_COLS;
        const int colsPadded = blocks * BLOCK_COLS;
        std::fill(ws.pM.begin(), ws.pM.begin() + colsPadded, 0.0f);
        std::fill(ws.pE.begin(), ws.pE.begin() + colsPadded, 0.0f);
        std::fill(ws.pF.begin(), ws.pF.begin() + colsPadded, 0.0f);
        std::fill(ws.ps.begin(), ws.ps.begin() + blocks, 0);
        if (!backward) {
            std::fill(ws.acc.begin(), ws.acc.begin() + (size_t)blocks * VECSIZE_FLOAT, 0.0f);
            std::fill(ws.accS.begin(), ws.accS.begin() + blocks, 0);
            ws.zv = 0;
            ws.zs = 0;
        }
        // gPow[k] = expGe^k
        float gPow[BLOCK_COLS + 1];
        gPow[0] = 1.0f;
        for (int k = 1; k <= BLOCK_COLS; ++k) {
            gPow[k] = gPow[k - 1] * expGe;
        }
        const simd_float vGo = simdf32_set(expGo);
        const simd_float vGe = simdf32_set(expGe);
        const simd_float vGe2 = simdf32_set(gPow[2]);
        const simd_float vGe4 = simdf32_set(gPow[4]);
        const simd_float vGeLast = simdf32_set(gPow[VECSIZE_FLOAT - 1]);
        // (1, g, ..., g^(VECSIZE_FLOAT-1))
        const simd_float vGeLanes = simdf32_loadu(gPow);
        // g^(v * VECSIZE_FLOAT): decay from the block start to vector v
        simd_float vGeVec[BLOCK_VECS];
        for (int v = 0; v < BLOCK_VECS; ++v) {
            vGeVec[v] = simdf32_set(gPow[v * VECSIZE_FLOAT]);
        }
        const simd_float hi = simdf32_set(4294967296.0f);
        const simd_float lo = simdf32_set(1.0f / 4294967296.0f);
        const int allLanes = (1 << VECSIZE_FLOAT) - 1;
        const simd_float lane0 = simdi_i2fcast(simdi32_insert(simdi_setzero(), -1, 0));
        simd_float bad = simdf32_setzero();
        simd_float vMaxP = simdf32_setzero();
        unsigned int flags = 0;
        float* pM = ws.pM.data();
        float* pE = ws.pE.data();
        float* pF = ws.pF.data();
        int* ps = ws.ps.data();

        simd_float e[BLOCK_VECS];
        simd_float q[BLOCK_VECS];
        simd_float M[BLOCK_VECS];
        simd_float E[BLOCK_VECS];
        simd_float F[BLOCK_VECS];
        for (int il = 0; il < rows; ++il) {
            const int i = backward ? rows - 1 - il : il;
            const float* scoreRow = score[i];
            float* probRow = prob[i];
            float* expRow = ws.eS.data() + (size_t)i * colsPadded;
            int* fsRow = ws.fs.data() + (size_t)i * blocks;
            // E entering the block (this row) and its scale
            float cv = 0.0f;
            int cs = 0;
            // D = M + E + F of the previous row left of the block and its scale
            float leftD = 0.0f;
            int leftS = 0;
            for (int b = 0; b < blocks; ++b) {
                const int j0 = b * BLOCK_COLS;
                const int sPrev = ps[b];
                int w = std::max(sPrev, leftS);
                const float fL = pow2f(sPrev - w);
                const float fN = pow2f(leftS - w);
                const simd_float one = simdf32_set(pow2f(-w));
                const simd_float vfL = simdf32_set(fL);
                // predecessor + 1: D of the previous row shifted right by one column, its first column from
                // the block to the left
                simd_float carryD = simdf32_set(leftD * fN);
                simd_float lastD = simdf32_setzero();
                for (int v = 0; v < BLOCK_VECS; ++v) {
                    const int j = j0 + v * VECSIZE_FLOAT;
                    // physical column of vector v
                    const int p = backward ? colsPadded - (v + 1) * VECSIZE_FLOAT - j0 : j;
                    if (!backward) {
                        // exp once; the backward pass reads it back (padding columns hold 0)
                        e[v] = simdf32_exp(loadCols(scoreRow, p, cols, -std::numeric_limits<float>::infinity()));
                        simdf32_storeu(expRow + p, e[v]);
                    } else {
                        e[v] = simdf32_reverse(simdf32_loadu(expRow + p));
                    }
                    const simd_float vPrevM = simdf32_loadu(pM + j);
                    const simd_float vPrevF = simdf32_loadu(pF + j);
                    const simd_float D = simdf32_add(simdf32_add(vPrevM, simdf32_loadu(pE + j)), vPrevF);
                    const simd_float DfL = simdf32_mul(D, vfL);
                    q[v] = simdf32_add(simdf32_blendv_ps(simdf32_shiftl1(DfL), carryD, lane0), one);
                    carryD = simdf32_bcast_last(DfL);
                    lastD = D;
                    M[v] = simdf32_mul(q[v], e[v]);
                    F[v] = simdf32_mul(simdf32_fmadd(vPrevM, vGo, simdf32_mul(vPrevF, vGe)), vfL);
                }
                const int wq = w;
                if (cv != 0.0f && cs > w) {
                    // carry on a larger scale: lift the block to it
                    const simd_float f = simdf32_set(pow2f(w - cs));
                    for (int v = 0; v < BLOCK_VECS; ++v) {
                        M[v] = simdf32_mul(M[v], f);
                        F[v] = simdf32_mul(F[v], f);
                    }
                    w = cs;
                }
                const float c = (cv == 0.0f) ? 0.0f : cv * pow2f(cs - w);
                // E along the row: decayed prefix sums of M * eGo, plus the carry c; t is E entering vector v
                // without c
                const simd_float vc = simdf32_set(c);
                simd_float t = simdf32_setzero();
                for (int v = 0; v < BLOCK_VECS; ++v) {
                    const simd_float x = simdf32_mul(M[v], vGo);
                    const simd_float y = exclusiveScan(x, vGe, vGe2, vGe4);
                    E[v] = simdf32_fmadd(simdf32_fmadd(vc, vGeVec[v], t), vGeLanes, y);
                    t = simdf32_fmadd(simdf32_fmadd(t, vGeLast, simdf32_bcast_last(y)), vGe, simdf32_bcast_last(x));
                }
                const float cOut = simdf32_first(t) + c * gPow[BLOCK_COLS];
                leftD = simdf32_first(simdf32_bcast_last(lastD));
                leftS = sPrev;

                if (backward) {
                    // posterior P = fwdM * qb / Z
                    const int k = fsRow[blocks - 1 - b] + wq - zExp;
                    const simd_float fct = simdf32_set(pow2f(k / 2) * pow2f(k - k / 2) * zInvMant);
                    for (int v = 0; v < BLOCK_VECS; ++v) {
                        const int p = colsPadded - (v + 1) * VECSIZE_FLOAT - j0;
                        const simd_float fwdM = simdf32_reverse(loadCols(probRow, p, cols, 0.0f));
                        const simd_float P = simdf32_mul(simdf32_mul(fwdM, q[v]), fct);
                        storeCols(probRow, p, cols, simdf32_reverse(P));
                        vMaxP = simdf32_max(vMaxP, P);
                        bad = simdf32_or(bad, simdf32_isnan(P));
                    }
                }

                // lazy renormalisation: only when some value leaves [2^-32, 2^32]
                simd_float mx = simdf32_setzero();
                for (int v = 0; v < BLOCK_VECS; ++v) {
                    const simd_float D = simdf32_add(simdf32_add(M[v], E[v]), F[v]);
                    mx = (v == 0) ? D : simdf32_max(mx, D);
                    bad = simdf32_or(bad, simdf32_isnan(D));
                }
                int sn = w;
                if (simdf32_movemask_ps(simdf32_gt(mx, hi)) != 0 ||
                    (w > 0 && simdf32_movemask_ps(simdf32_lt(mx, lo)) == allLanes)) {
                    const float m = simdf32_hmax(mx);
                    if (!(m < 1.0e38f)) {
                        flags |= FLAG_NON_FINITE;
                    } else {
                        sn = (m > 0.0f) ? std::max(0, w + fexp(m)) : 0;
                    }
                }
                if (sn != w) {
                    const simd_float f = simdf32_set(pow2f(w - sn));
                    for (int v = 0; v < BLOCK_VECS; ++v) {
                        M[v] = simdf32_mul(M[v], f);
                        E[v] = simdf32_mul(E[v], f);
                        F[v] = simdf32_mul(F[v], f);
                    }
                }
                for (int v = 0; v < BLOCK_VECS; ++v) {
                    const int j = j0 + v * VECSIZE_FLOAT;
                    simdf32_storeu(pM + j, M[v]);
                    simdf32_storeu(pE + j, E[v]);
                    simdf32_storeu(pF + j, F[v]);
                }
                ps[b] = sn;
                cv = cOut;
                cs = (cOut == 0.0f) ? 0 : w;

                if (!backward) {
                    simd_float sumM = M[0];
                    for (int v = 0; v < BLOCK_VECS; ++v) {
                        storeCols(probRow, j0 + v * VECSIZE_FLOAT, cols, M[v]);
                        if (v > 0) {
                            sumM = simdf32_add(sumM, M[v]);
                        }
                    }
                    fsRow[b] = sn;
                    float* acc = ws.acc.data() + (size_t)b * VECSIZE_FLOAT;
                    if (sn != ws.accS[b]) {
                        zAdd(ws, hsum(simdf32_loadu(acc)), ws.accS[b]);
                        simdf32_storeu(acc, simdf32_setzero());
                        ws.accS[b] = sn;
                    }
                    simdf32_storeu(acc, simdf32_add(simdf32_loadu(acc), sumM));
                }
            }
        }
        if (simdf32_movemask_ps(bad) != 0) {
            flags |= FLAG_NON_FINITE;
        }
        if (!backward) {
            for (int b = 0; b < blocks; ++b) {
                zAdd(ws, hsum(simdf32_loadu(ws.acc.data() + (size_t)b * VECSIZE_FLOAT)), ws.accS[b]);
            }
            if (!(ws.zv > 0) || !std::isfinite(ws.zv)) {
                flags |= FLAG_Z;
            }
        } else {
            maxP = simdf32_hmax(vMaxP);
            if (!(maxP <= 1.0f + 1e-3f)) {
                flags |= FLAG_P;
            }
        }
        return flags;
    }
};

#endif
