// Forward-backward kernel, textually included once per instruction set from FwbwDefault.cpp and
// FwbwAVX512.cpp. It is not a stand-alone translation unit and must not be added to the build
// directly.

struct States {
    const static uint8_t STOP=0;
    const static uint8_t M=1;
    const static uint8_t I=2;
    const static uint8_t D=3;
};

struct FWBWState {
    const static bool FORWARD = true;
    const static bool BACKWARD = false;
};

const char *simdArchitecture() {
#if defined(AVX512)
    return "AVX512 512-bit";
#elif defined(MMSEQS_AVX512VL_MASKS)
    return "AVX512VL 256-bit";
#elif defined(AVX2)
    return "AVX2 256-bit";
#else
    return "SSE4.1 128-bit";
#endif
}

inline void calculate_max4(float& max, float& term1, float& term2, float& term3, float& term4, uint8_t& state) {
    if (term1 > term2) { max = term1; state = States::STOP; }
    else { max = term2; state = States::M; }
    if (term3 > max) { max = term3; state = States::I; }
    if (term4 > max) { max = term4; state = States::D; }
}

inline simd_float simdf32_prefixsum(simd_float a) {
#if defined(AVX2)
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl_4(simdf_f2icast(a))));
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl_8(simdf_f2icast(a))));
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl_16(simdf_f2icast(a))));
    return a; 
#elif defined(AVX512)
    a = simdf32_add(a, simdi_i2fcast(_mm512_shift_left4(simdf_f2icast(a))));
    a = simdf32_add(a, simdi_i2fcast(_mm512_shift_left8(simdf_f2icast(a))));
    a = simdf32_add(a, simdi_i2fcast(_mm512_shift_left16(simdf_f2icast(a))));
    a = simdf32_add(a, simdi_i2fcast(_mm512_shift_left32(simdf_f2icast(a)))); 
    return a;

#else
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl(simdf_f2icast(a), 4)));
    a = simdf32_add(a, simdi_i2fcast(simdi8_shiftl(simdf_f2icast(a), 8)));
    return a;
#endif
}

// Helpers of ScaledFwBw and MacFill

// simdf32_movemask_ps with every lane set
const int ALL_LANES = (1 << VECSIZE_FLOAT) - 1;

// 2^e at compile time
constexpr float pow2c(int e) {
    return e == 0 ? 1.0f : (e < 0 ? 0.5f * pow2c(e + 1) : 2.0f * pow2c(e - 1));
}

// row[p, p + VECSIZE_FLOAT) with lanes at or after cols set to pad
inline simd_float loadCols(const float *row, int p, int cols, float pad) {
    if (p + VECSIZE_FLOAT <= cols) {
        return simdf32_loadu(row + p);
    }
    if (p >= cols) {
        return simdf32_set(pad);
    }
    float t[VECSIZE_FLOAT];
    for (int k = 0; k < VECSIZE_FLOAT; ++k) {
        t[k] = (p + k < cols) ? row[p + k] : pad;
    }
    return simdf32_loadu(t);
}

// Scaled-mantissa forward-backward, used by runFwBw (score matrix or profile, any backtrace mode) when
// FwBwAligner::setScaledMode is set. Same model as the log-space blocks below:
//   ZM_ij = (ZM_{i-1,j-1} + ZE_{i-1,j-1} + ZF_{i-1,j-1} + 1) * exp(S_ij)      (S already divided by T)
//   ZE_ij = ZM_{i,j-1} * eGo + ZE_{i,j-1} * eGe
//   ZF_ij = ZM_{i-1,j} * eGo + ZF_{i-1,j} * eGe
//   Z = sum ZM,  P_ij = ZM_ij * ZMb_ij / (exp(S_ij) * Z)
// ZMb is the same recurrence on the reversed matrix. ZMb_ij / exp(S_ij) is its "predecessor + 1" term,
// which the backward pass has at hand, so the posterior needs no division and no second exp.
//
// Values are stored as float mantissas with one power-of-two scale per (row, BLOCK_COLS-column block):
// 16 columns for a score matrix (LoL-align's steep scores lose paths far more often in wider blocks), 32 for
// a profile (substitution scores are flat enough; 12-17% faster per cell). A block is rescaled (exponent
// arithmetic only) when its values leave [2^-RENORM_BITS, 2^RENORM_BITS]; the local-start term 1 bounds
// the values from below, so scales never go negative.
//
// Lost paths: values sharing a block's scale can only lie within ~2^126 of the block's largest value.
// As a block is rescaled only when needed, it can keep a large scale after a strong path has left it, and a
// path that starts there underflows without any non-finite value. Z of the backward pass then differs from
// the forward pass's: run<true, ...> sums both and flags the difference (FLAG_ZCHECK), and runScaledKernel
// runs the call again in runScaledDouble().
//
// The kernel (lane helpers in simd.h, decayedScan) covers 4 (SSE, NEON), 8 (AVX2, AVX512VL at 256 bit)
// and 16 (AVX512 with MMSEQS_NO_SIMD256) lanes.
static_assert(VECSIZE_FLOAT == 4 || VECSIZE_FLOAT == 8 || VECSIZE_FLOAT == 16, "4, 8 or 16 float lanes");

template<int BLOCK_COLS>
class ScaledFwBw {
public:
    // flags of run() and runScaledDouble(); prob and maxP are invalid if any is set
    static const unsigned int FLAG_NON_FINITE = 1; // exp saturated or NaN
    static const unsigned int FLAG_Z = 2;          // Z is 0 or not finite
    static const unsigned int FLAG_P = 4;          // a posterior above 1
    static const unsigned int FLAG_ZCHECK = 8;     // backward Z differs from forward Z by more than 1e-4 (runScaledDouble: 1e-6)

    // a block is rescaled only when its values leave [2^-RENORM_BITS, 2^RENORM_BITS]
    static const int RENORM_BITS = 32;

    // score, prob: rows x cols, score already divided by T; with expScores, score holds exp(S) instead
    template<bool checkZ, bool expScores = false>
    unsigned int run(float **score, int rows, int cols, float expGo, float expGe, float **prob, float &maxP) {
        // flush-to-zero for the duration of the call (x86 only): values far below their block's scale are
        // negligible, but denormal arithmetic costs 10-27% of the run time
        SimdFlushDenormals ftz;
        const int blocks = (cols + BLOCK_COLS - 1) / BLOCK_COLS;
        const int colsPadded = blocks * BLOCK_COLS;
        if (prevScale.size() < (size_t)blocks) {
            prevM.resize(colsPadded);
            prevF.resize(colsPadded);
            diagRow[0].resize(colsPadded + 1);
            diagRow[1].resize(colsPadded + 1);
            prevScale.resize(blocks);
        }
        if (fwdScale.size() < (size_t)rows * blocks) {
            fwdScale.resize((size_t)rows * blocks);
        }
        // with expScores, the backward pass reads exp(S) from score itself
        if (!expScores && expScore.size() < (size_t)rows * colsPadded) {
            expScore.resize((size_t)rows * colsPadded);
        }
        float unused = 0.0f;
        unsigned int flags = pass<false, checkZ, expScores>(score, rows, cols, expGo, expGe, prob, 0, 1.0f, unused);
        if (flags != 0) {
            return flags;
        }
        // Z = m * 2^zExp with m in [0.5, 1)
        int e = 0;
        const double m = std::frexp(zFwd.value, &e);
        flags = pass<true, checkZ, expScores>(score, rows, cols, expGo, expGe, prob, zFwd.exp + e, (float)(1.0 / m), maxP);
        if (checkZ) {
            const double ratio = (zBwd.value / zFwd.value) * std::ldexp(1.0, zBwd.exp - zFwd.exp);
            if (!(std::fabs(ratio - 1.0) <= 1e-4)) {
                flags |= FLAG_ZCHECK;
            }
        }
        return flags;
    }

    // The same model with double mantissas, one power-of-two scale per (row, DOUBLE_BLOCK_COLS-column block),
    // rescaled in every block; scalar. Values sharing a scale can lie ~2^1022 (e^708) apart instead of
    // ~2^126, which keeps the paths run() loses next to much stronger ones. ~5x slower per cell than
    // run(), ~10x faster than runDouble(). score: S / T (never exp(S)); prob, maxP and the flags as for run().
    unsigned int runScaledDouble(float **score, int rows, int cols, float expGo, float expGe, float **prob, float &maxP) {
        const int blocks = (cols + DOUBLE_BLOCK_COLS - 1) / DOUBLE_BLOCK_COLS;
        if (dM.size() < (size_t)rows * cols) {
            dM.resize((size_t)rows * cols);
        }
        if (dExpScore.size() < (size_t)rows * cols) {
            dExpScore.resize((size_t)rows * cols);
        }
        if (dFwdScale.size() < (size_t)rows * blocks) {
            dFwdScale.resize((size_t)rows * blocks);
        }
        double unused = 0.0;
        unsigned int flags = passScaledDouble<false>(score, rows, cols, expGo, expGe, prob, 0, 1.0, unused);
        if (flags != 0) {
            return flags;
        }
        if (!(zFwd.value > 0) || !std::isfinite(zFwd.value)) {
            return FLAG_Z;
        }
        // Z = m * 2^zExp with m in [0.5, 1)
        int e = 0;
        const double m = std::frexp(zFwd.value, &e);
        double maxPd = 0.0;
        flags = passScaledDouble<true>(score, rows, cols, expGo, expGe, prob, zFwd.exp + e, 1.0 / m, maxPd);
        const double ratio = (zBwd.value / zFwd.value) * std::ldexp(1.0, zBwd.exp - zFwd.exp);
        if (!(std::fabs(ratio - 1.0) <= 1e-6)) {
            flags |= FLAG_ZCHECK;
        }
        if (!(maxPd <= 1.0 + 1e-3)) {
            flags |= FLAG_P;
        }
        maxP = (float)maxPd;
        return flags;
    }

    // frees the buffers of runScaledDouble() and runDouble() (rows x cols doubles): only the rare calls that
    // fall back need them
    void freeFallbackBuffers() {
        std::vector<double>().swap(dM);
        std::vector<double>().swap(dExpScore);
        std::vector<double>().swap(dPred);
        for (int k = 0; k < 6; ++k) {
            std::vector<double>().swap(dRow[k]);
        }
        std::vector<int>().swap(dFwdScale);
        std::vector<int>().swap(dPrevScale);
        std::vector<int>().swap(dCurrScale);
    }

    // The same model in double precision and log space with log-sum-exp per cell: nothing shares a
    // scale, so no path is lost. ~50x slower per cell than run(), the last resort for the calls
    // runScaledDouble() flags.
    // Z = 0 (all scores -inf) gives P = 0 and maxP = 0, like an overflow of the log-space blocks.
    void runDouble(float **score, int rows, int cols, float expGo, float expGe, float **prob, float &maxP) {
        const double NI = -std::numeric_limits<double>::infinity();
        const double logGo = std::log((double)expGo);
        const double logGe = std::log((double)expGe);
        if (dM.size() < (size_t)rows * cols) {
            dM.resize((size_t)rows * cols);
        }
        double *M = dM.data();
        std::vector<double> &Eprev = dRow[0], &Ecurr = dRow[1], &Fprev = dRow[2], &Fcurr = dRow[3];
        Eprev.assign(cols, NI);
        Ecurr.assign(cols, NI);
        Fprev.assign(cols, NI);
        Fcurr.assign(cols, NI);
        double logZ = NI;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                const double s = score[i][j];
                const double pred = (i > 0 && j > 0) ? lse3(M[(size_t)(i - 1) * cols + j - 1], Eprev[j - 1], Fprev[j - 1]) : NI;
                const double m = (s == NI) ? NI : s + lse2(pred, 0.0);
                M[(size_t)i * cols + j] = m;
                Fcurr[j] = (i > 0) ? lse2(M[(size_t)(i - 1) * cols + j] + logGo, Fprev[j] + logGe) : NI;
                Ecurr[j] = (j > 0) ? lse2(M[(size_t)i * cols + j - 1] + logGo, Ecurr[j - 1] + logGe) : NI;
                logZ = lse2(logZ, m);
            }
            std::swap(Eprev, Ecurr);
            std::swap(Fprev, Fcurr);
        }
        maxP = 0.0f;
        if (logZ == NI) {
            for (int i = 0; i < rows; ++i) {
                std::fill(prob[i], prob[i] + cols, 0.0f);
            }
            return;
        }
        // backward recurrence in forward coordinates, posterior on the fly
        std::vector<double> &BMnext = dRow[0], &BEnext = dRow[1], &BFnext = dRow[2];
        std::vector<double> &BMcurr = dRow[3], &BEcurr = dRow[4], &BFcurr = dRow[5];
        for (int k = 0; k < 6; ++k) {
            dRow[k].assign(cols, NI);
        }
        double maxPd = 0.0;
        for (int i = rows - 1; i >= 0; --i) {
            for (int j = cols - 1; j >= 0; --j) {
                const double s = score[i][j];
                const double pred = (i < rows - 1 && j < cols - 1) ? lse3(BMnext[j + 1], BEnext[j + 1], BFnext[j + 1]) : NI;
                const double bm = (s == NI) ? NI : s + lse2(pred, 0.0);
                BMcurr[j] = bm;
                BFcurr[j] = (i < rows - 1) ? lse2(BMnext[j] + logGo, BFnext[j] + logGe) : NI;
                BEcurr[j] = (j < cols - 1) ? lse2(BMcurr[j + 1] + logGo, BEcurr[j + 1] + logGe) : NI;
                const double p = (s == NI) ? 0.0 : std::exp(M[(size_t)i * cols + j] + bm - s - logZ);
                prob[i][j] = (float)p;
                maxPd = std::max(maxPd, p);
            }
            std::swap(BMnext, BMcurr);
            std::swap(BEnext, BEcurr);
            std::swap(BFnext, BFcurr);
        }
        maxP = (float)maxPd;
    }

private:
    static const int BLOCK_VECS = BLOCK_COLS / VECSIZE_FLOAT;
    static const int DOUBLE_BLOCK_COLS = 16; // runScaledDouble

    // Z = value * 2^exp. Each block sums its ZM in float at its current scale and hands the sum to the
    // double total whenever the scale changes.
    struct ZSum {
        std::vector<float> blockSum; // [blocks * VECSIZE_FLOAT]
        std::vector<int> blockScale; // [blocks]
        double value = 0;
        int exp = 0;

        void reset(int blocks) {
            if (blockScale.size() < (size_t)blocks) {
                blockScale.resize(blocks);
                blockSum.resize((size_t)blocks * VECSIZE_FLOAT);
            }
            std::fill(blockSum.begin(), blockSum.begin() + (size_t)blocks * VECSIZE_FLOAT, 0.0f);
            std::fill(blockScale.begin(), blockScale.begin() + blocks, 0);
            value = 0;
            exp = 0;
        }
        // add to this total
        void add(double v, int s) {
            addTo(v, s, value, exp);
        }
        // total * 2^totalExp += v * 2^s; a power-of-two product instead of ldexp, which would be a call inside
        // pass()'s block loop (GCC then saves and restores the vector registers around it in every block).
        // Both are exact for normal results; a shift below -1022 gives 0 instead of a denormal, which the
        // following addition of a normal value absorbs either way.
        static void addTo(double v, int s, double &total, int &totalExp) {
            if (v == 0) {
                return;
            }
            if (total == 0) {
                total = v;
                totalExp = s;
            } else if (s > totalExp) {
                total = total * pow2d(totalExp - s) + v;
                totalExp = s;
            } else {
                total += v * pow2d(s - totalExp);
            }
        }
        // sums, scales: blockSum.data(), blockScale.data(); total, totalExp: the running Z, all held by the caller
        // (see pass()) so that they stay in registers
        static void addBlock(float *sums, int *scales, int b, simd_float m, int scale, double &total, int &totalExp) {
            float *sum = sums + (size_t)b * VECSIZE_FLOAT;
            if (scale != scales[b]) {
                addTo(simdf32_hsum_double(simdf32_loadu(sum)), scales[b], total, totalExp);
                simdf32_storeu(sum, simdf32_setzero());
                scales[b] = scale;
            }
            simdf32_storeu(sum, simdf32_add(simdf32_loadu(sum), m));
        }
        void finish(int blocks) {
            for (int b = 0; b < blocks; ++b) {
                add(simdf32_hsum_double(simdf32_loadu(blockSum.data() + (size_t)b * VECSIZE_FLOAT)), blockScale[b]);
            }
        }
    };

    // One pass over the matrix. Forward: rows top to bottom, stores the ZM mantissas in prob and their
    // scales in fwdScale, sums Z. Backward: the same recurrence on the reversed matrix (row rows-1-i,
    // column colsPadded-1-j, so backward block b covers forward block blocks-1-b), overwrites prob with P.
    // zExp, zInvMant: Z of the forward pass (backward only)
    template<bool backward, bool checkZ, bool expScores>
    unsigned int pass(float **score, int rows, int cols, float expGo, float expGe, float **prob,
                      int zExp, float zInvMant, float &maxP) {
        const int blocks = (cols + BLOCK_COLS - 1) / BLOCK_COLS;
        const int colsPadded = blocks * BLOCK_COLS;
        std::fill(prevM.begin(), prevM.begin() + colsPadded, 0.0f);
        std::fill(prevF.begin(), prevF.begin() + colsPadded, 0.0f);
        std::fill(diagRow[0].begin(), diagRow[0].begin() + colsPadded + 1, 0.0f);
        std::fill(diagRow[1].begin(), diagRow[1].begin() + colsPadded + 1, 0.0f);
        std::fill(prevScale.begin(), prevScale.begin() + blocks, 0);
        if (!backward) {
            zFwd.reset(blocks);
        } else if (checkZ) {
            zBwd.reset(blocks);
        }
        // the buffers through local pointers: a vector store may alias any memory, the pointers inside the
        // std::vectors included, so GCC would load them again through this after every store
        float *const prevMPtr = prevM.data();
        float *const prevFPtr = prevF.data();
        int *const prevScalePtr = prevScale.data();
        float *const diagRowPtr[2] = {diagRow[0].data(), diagRow[1].data()};
        ZSum &z = backward ? zBwd : zFwd;
        float *const zBlockSum = z.blockSum.data();
        int *const zBlockScale = z.blockScale.data();
        // running Z = zTotal * 2^zTotalExp, in locals (registers) during the pass, written back to z before finish()
        double zTotal = 0.0;
        int zTotalExp = 0;
        // the block with padding columns: the last one forward, the first one backward
        const int edgeBlock = (colsPadded == cols) ? -1 : (backward ? 0 : blocks - 1);
        // gePow[k] = expGe^k
        float gePow[BLOCK_COLS + 1];
        gePow[0] = 1.0f;
        for (int k = 1; k <= BLOCK_COLS; ++k) {
            gePow[k] = gePow[k - 1] * expGe;
        }
        const simd_float vGo = simdf32_set(expGo);
        const simd_float vGe = simdf32_set(expGe);
        const simd_float geStep[4] = {simdf32_set(gePow[1]), simdf32_set(gePow[2]), simdf32_set(gePow[4]),
                                      simdf32_set(gePow[8])};
        const simd_float vGeLast = simdf32_set(gePow[VECSIZE_FLOAT - 1]);
        const simd_float vGeLanes = simdf32_loadu(gePow); // (1, ge, ..., ge^(VECSIZE_FLOAT-1))
        simd_float vGeVec[BLOCK_VECS];                     // ge^(v * VECSIZE_FLOAT): block start to vector v
        for (int v = 0; v < BLOCK_VECS; ++v) {
            vGeVec[v] = simdf32_set(gePow[v * VECSIZE_FLOAT]);
        }
        const simd_float vHi = simdf32_set(pow2c(RENORM_BITS));
        const simd_float vLo = simdf32_set(pow2c(-RENORM_BITS));
        const simd_float lane0 = simdi_i2fcast(simdi32_insert(simdi_setzero(), -1, 0));
        simd_float nan = simdf32_setzero();
        simd_float vMaxP = simdf32_setzero();
        unsigned int flags = 0;

        simd_float expS[BLOCK_VECS], pred[BLOCK_VECS], M[BLOCK_VECS], E[BLOCK_VECS], F[BLOCK_VECS], diag[BLOCK_VECS];
        for (int row = 0; row < rows; ++row) {
            const int i = backward ? rows - 1 - row : row;
            float *expRow = expScores ? nullptr : expScore.data() + (size_t)i * colsPadded;
            int *fwdScaleRow = fwdScale.data() + (size_t)i * blocks;
            const float *scoreRow = score[i];
            float *probRow = prob[i];
            // diag = M + E + F of the previous and of this row, column j at [j + 1] ([0] = 0 left of the matrix)
            const float *diagAbove = diagRowPtr[row & 1];
            float *diagCurr = diagRowPtr[(row + 1) & 1];
            // E entering the block from the left, and its scale
            float carryE = 0.0f;
            int carryScale = 0;
            // scale of the previous row in the block to the left
            int diagLeftScale = 0;
            // only the block at the matrix edge holds padding columns: the others load and store without
            // the column checks (two copies of the block body, edgeBlock < 0 when there are none).
            // processBlock updates the row state (carryE, carryScale, diagLeftScale) and the pass state (flags,
            // nan, vMaxP, the Z total); it is inlined into both copies so that this state stays in registers
            const auto processBlock = [&](const int b, auto edgeTag) __attribute__((__always_inline__)) {
                constexpr bool edge = decltype(edgeTag)::value;
                const int j0 = b * BLOCK_COLS;
                const int scaleAbove = prevScalePtr[b];
                int scale = std::max(scaleAbove, diagLeftScale);
                const simd_float fromAbove = simdf32_set(pow2fNonPos(scaleAbove - scale));
                const simd_float localStart = simdf32_set(pow2fNonPos(-scale));
                // pred = diag of the previous row one column to the left, plus the local start; the column left
                // of the block is on that block's scale
                const simd_float carryDiag = simdf32_set(diagAbove[j0] * pow2fNonPos(diagLeftScale - scale));
                for (int v = 0; v < BLOCK_VECS; ++v) {
                    const int j = j0 + v * VECSIZE_FLOAT;
                    if (!backward) {
                        // forward: exp(S) once, kept in expScore for the backward pass (padding columns get
                        // exp(-inf) = 0); with expScores, score already holds exp(S) and the backward pass reads it
                        if (expScores) {
                            expS[v] = loadColsIn<edge>(scoreRow, j, cols, 0.0f);
                        } else {
                            expS[v] = simdf32_exp(loadColsIn<edge>(scoreRow, j, cols, -std::numeric_limits<float>::infinity()));
                            simdf32_storeu(expRow + j, expS[v]);
                        }
                    } else {
                        expS[v] = simdf32_reverse(
                            expScores ? loadColsIn<edge>(scoreRow, colsPadded - j - VECSIZE_FLOAT, cols, 0.0f)
                                      : simdf32_loadu(expRow + colsPadded - j - VECSIZE_FLOAT));
                    }
                    const simd_float mAbove = simdf32_loadu(prevMPtr + j);
                    const simd_float fAbove = simdf32_loadu(prevFPtr + j);
                    const simd_float diagScaled =
                        simdf32_no_contract(simdf32_mul(simdf32_loadu(diagAbove + j), fromAbove));
                    pred[v] = simdf32_add((v == 0) ? simdf32_blendv_ps(diagScaled, carryDiag, lane0) : diagScaled, localStart);
                    M[v] = simdf32_mul(pred[v], expS[v]);
                    F[v] = simdf32_mul(simdf32_fmadd(mAbove, vGo, simdf32_mul(fAbove, vGe)), fromAbove);
                }
                const int predScale = scale;
                // without E from the left carryScale is 0, and scales are >= 0
                if (carryScale > scale) {
                    // E from the left is on a larger scale: lift the block to it
                    const simd_float f = simdf32_set(pow2fNonPos(scale - carryScale));
                    for (int v = 0; v < BLOCK_VECS; ++v) {
                        M[v] = simdf32_mul(M[v], f);
                        F[v] = simdf32_mul(F[v], f);
                    }
                    scale = carryScale;
                }
                // E along the row: ge-decayed prefix sums of M * go, plus the carry from the left
                const float eIn = carryE * pow2fNonPos(carryScale - scale); // 0 without E from the left
                const simd_float vEIn = simdf32_set(eIn);
                simd_float eRun = simdf32_setzero(); // E entering vector v, without eIn
                for (int v = 0; v < BLOCK_VECS; ++v) {
                    const simd_float x = simdf32_mul(M[v], vGo);
                    const simd_float y = decayedScan(x, geStep);
                    E[v] = simdf32_fmadd(simdf32_fmadd(vEIn, vGeVec[v], eRun), vGeLanes, y);
                    eRun = simdf32_fmadd(simdf32_fmadd(eRun, vGeLast, simdf32_broadcast_last(y)), vGe,
                                         simdf32_broadcast_last(x));
                }
                const float eOut = simdf32_first(eRun) + eIn * gePow[BLOCK_COLS];
                diagLeftScale = scaleAbove;

                // rescale the block if its values left [2^-RENORM_BITS, 2^RENORM_BITS]
                simd_float diagMax = simdf32_setzero();
                simd_float diagSum = simdf32_setzero(); // NaN iff a diag is NaN (all >= 0, so no inf - inf)
                for (int v = 0; v < BLOCK_VECS; ++v) {
                    diag[v] = simdf32_add(simdf32_add(M[v], E[v]), F[v]);
                    diagMax = (v == 0) ? diag[v] : simdf32_max(diagMax, diag[v]);
                    diagSum = (v == 0) ? diag[v] : simdf32_add(diagSum, diag[v]);
                }
                nan = simdf32_or(nan, simdf32_isnan(diagSum));
                int newScale = scale;
                if (simdf32_movemask_ps(simdf32_gt(diagMax, vHi)) != 0 ||
                    (scale > 0 && simdf32_movemask_ps(simdf32_lt(diagMax, vLo)) == ALL_LANES)) {
                    const float m = simdf32_hmax(diagMax);
                    if (!(m < 1.0e38f)) {
                        flags |= FLAG_NON_FINITE;
                    } else {
                        newScale = (m > 0.0f) ? std::max(0, scale + floatExp(m)) : 0;
                    }
                }
                if (newScale != scale) {
                    const simd_float f = simdf32_set(pow2f(scale - newScale));
                    for (int v = 0; v < BLOCK_VECS; ++v) {
                        M[v] = simdf32_no_contract(simdf32_mul(M[v], f));
                        E[v] = simdf32_no_contract(simdf32_mul(E[v], f));
                        F[v] = simdf32_no_contract(simdf32_mul(F[v], f));
                        // from the rescaled values, as the next row would sum them
                        diag[v] = simdf32_add(simdf32_add(M[v], E[v]), F[v]);
                    }
                }
                simd_float sumM = M[0];
                for (int v = 0; v < BLOCK_VECS; ++v) {
                    const int j = j0 + v * VECSIZE_FLOAT;
                    simdf32_storeu(prevMPtr + j, M[v]);
                    simdf32_storeu(prevFPtr + j, F[v]);
                    simdf32_storeu(diagCurr + j + 1, diag[v]);
                    if (v > 0) {
                        sumM = simdf32_add(sumM, M[v]);
                    }
                }
                prevScalePtr[b] = newScale;
                carryE = eOut;
                carryScale = (eOut == 0.0f) ? 0 : scale;

                if (!backward) {
                    // forward: keep ZM (mantissas in prob, scale in fwdScale) and sum Z
                    for (int v = 0; v < BLOCK_VECS; ++v) {
                        storeColsIn<edge>(probRow, j0 + v * VECSIZE_FLOAT, cols, M[v]);
                    }
                    fwdScaleRow[b] = newScale;
                    ZSum::addBlock(zBlockSum, zBlockScale, b, sumM, newScale, zTotal, zTotalExp);
                } else {
                    // backward: P = M_fwd * pred_bwd / Z overwrites ZM in prob; pred_bwd is at predScale,
                    // 2^pExp is split in two so neither factor overflows
                    const int pExp = fwdScaleRow[blocks - 1 - b] + predScale - zExp;
                    const simd_float pFactor = simdf32_set(pow2f(pExp / 2) * pow2f(pExp - pExp / 2) * zInvMant);
                    simd_float pSum = simdf32_setzero(); // NaN iff a P is NaN (all P >= 0)
                    for (int v = 0; v < BLOCK_VECS; ++v) {
                        const int p = colsPadded - j0 - (v + 1) * VECSIZE_FLOAT;
                        const simd_float mFwd = simdf32_reverse(loadColsIn<edge>(probRow, p, cols, 0.0f));
                        const simd_float P = simdf32_mul(simdf32_mul(mFwd, pred[v]), pFactor);
                        storeColsIn<edge>(probRow, p, cols, simdf32_reverse(P));
                        vMaxP = simdf32_max(vMaxP, P);
                        pSum = simdf32_add(pSum, P);
                    }
                    nan = simdf32_or(nan, simdf32_isnan(pSum));
                    if (checkZ) {
                        ZSum::addBlock(zBlockSum, zBlockScale, b, sumM, newScale, zTotal, zTotalExp);
                    }
                }
            };
            for (int b = 0; b < blocks; ++b) {
                if (b == edgeBlock) {
                    processBlock(b, std::true_type());
                } else {
                    processBlock(b, std::false_type());
                }
            }
        }
        if (simdf32_movemask_ps(nan) != 0) {
            flags |= FLAG_NON_FINITE;
        }
        z.value = zTotal;
        z.exp = zTotalExp;
        if (!backward) {
            z.finish(blocks);
            if (!(zFwd.value > 0) || !std::isfinite(zFwd.value)) {
                flags |= FLAG_Z;
            }
        } else {
            if (checkZ) {
                z.finish(blocks);
            }
            maxP = simdf32_hmax(vMaxP);
            if (!(maxP <= 1.0f + 1e-3f)) {
                flags |= FLAG_P;
            }
        }
        return flags;
    }

    // One pass of runScaledDouble, the recurrence of pass() cell by cell. Forward: stores exp(S) in dExpScore,
    // the ZM mantissas in dM and their scales in dFwdScale, sums Z. Backward: the same recurrence on the
    // reversed matrix (row rows-1-row, column cols-1-c), writes P to prob and its maximum to maxPd.
    // zExp, zInvMant: Z of the forward pass (backward only)
    template<bool backward>
    unsigned int passScaledDouble(float **score, int rows, int cols, double expGo, double expGe, float **prob,
                            int zExp, double zInvMant, double &maxPd) {
        const int blocks = (cols + DOUBLE_BLOCK_COLS - 1) / DOUBLE_BLOCK_COLS;
        // M, F of the previous and the current row at [c]; diag = M + E + F at [c + 1] ([0] = 0)
        std::vector<double> &mPrev = dRow[0], &fPrev = dRow[1], &diagPrev = dRow[2];
        std::vector<double> &mCurr = dRow[3], &fCurr = dRow[4], &diagCurr = dRow[5];
        for (int k = 0; k < 6; ++k) {
            dRow[k].assign(cols + 1, 0.0);
        }
        dPred.resize(cols);
        dPrevScale.assign(blocks, 0);
        dCurrScale.assign(blocks, 0);
        ZSum &z = backward ? zBwd : zFwd;
        z.reset(0); // value and exp only: passScaledDouble adds whole sums, no block sums
        unsigned int flags = 0;
        for (int row = 0; row < rows; ++row) {
            const int i = backward ? rows - 1 - row : row;
            double *expRow = dExpScore.data() + (size_t)i * cols;
            double *mFwdRow = dM.data() + (size_t)i * cols;
            int *fwdScaleRow = dFwdScale.data() + (size_t)i * blocks;
            // M and E of the column left of the block, and their scale
            double mLeft = 0.0, eLeft = 0.0;
            int leftScale = 0;
            for (int b = 0; b < blocks; ++b) {
                const int c0 = b * DOUBLE_BLOCK_COLS;
                const int c1 = std::min(cols, c0 + DOUBLE_BLOCK_COLS);
                // the largest scale among the inputs: the block above, the previous row's column left of the
                // block (diagPrev[c0]) and this row's column left of the block
                int scale = dPrevScale[b];
                if (b > 0 && diagPrev[c0] != 0.0) {
                    scale = std::max(scale, dPrevScale[b - 1]);
                }
                if (mLeft != 0.0 || eLeft != 0.0) {
                    scale = std::max(scale, leftScale);
                }
                const double fromAbove = pow2d(dPrevScale[b] - scale);
                const double fromDiagLeft = (b > 0) ? pow2d(dPrevScale[b - 1] - scale) : 0.0;
                const double fromLeft = pow2d(leftScale - scale);
                // 2^-scale split in two, so that a scale above 1022 gives a denormal instead of 0
                const int halfScale = scale / 2;
                const double localStart = pow2d(-halfScale) * pow2d(halfScale - scale);
                double m = mLeft * fromLeft;
                double e = eLeft * fromLeft;
                double diagMax = 0.0;
                // c: column in pass order, j: column in the matrix
                for (int c = c0; c < c1; ++c) {
                    const int j = backward ? cols - 1 - c : c;
                    double expS;
                    if (!backward) {
                        expS = std::exp((double)score[i][j]);
                        expRow[j] = expS;
                    } else {
                        expS = expRow[j];
                    }
                    const double pred = diagPrev[c] * ((c == c0) ? fromDiagLeft : fromAbove) + localStart;
                    const double f = (mPrev[c] * expGo + fPrev[c] * expGe) * fromAbove;
                    // E from M and E of column c - 1
                    e = m * expGo + e * expGe;
                    m = pred * expS;
                    const double diag = m + e + f;
                    mCurr[c] = m;
                    fCurr[c] = f;
                    diagCurr[c + 1] = diag;
                    dPred[c] = pred;
                    diagMax = std::max(diagMax, diag);
                }
                // rescale so that the largest diag is in [1, 2)
                const int newScale = (diagMax > 0.0) ? scale + std::ilogb(diagMax) : 0;
                const double rescale = pow2d(scale - newScale);
                double sumM = 0.0;
                for (int c = c0; c < c1; ++c) {
                    mCurr[c] *= rescale;
                    fCurr[c] *= rescale;
                    diagCurr[c + 1] *= rescale;
                    sumM += mCurr[c];
                }
                // max() drops a NaN, the sum keeps it
                if (!(diagMax < std::numeric_limits<double>::infinity()) || std::isnan(sumM)) {
                    flags |= FLAG_NON_FINITE;
                }
                mLeft = m * rescale;
                eLeft = e * rescale;
                leftScale = newScale;
                dCurrScale[b] = newScale;
                if (!backward) {
                    for (int c = c0; c < c1; ++c) {
                        mFwdRow[c] = mCurr[c];
                    }
                    fwdScaleRow[b] = newScale;
                } else {
                    // P = M_fwd * pred_bwd / Z; pred_bwd is at scale, 2^pExp is split in two so neither
                    // factor overflows
                    for (int c = c0; c < c1; ++c) {
                        const int j = cols - 1 - c;
                        const int pExp = fwdScaleRow[j / DOUBLE_BLOCK_COLS] + scale - zExp;
                        const double p = mFwdRow[j] * dPred[c] * zInvMant * pow2d(pExp / 2) * pow2d(pExp - pExp / 2);
                        prob[i][j] = (float)p;
                        if (std::isnan(p)) {
                            flags |= FLAG_NON_FINITE;
                        }
                        maxPd = std::max(maxPd, p);
                    }
                }
                z.add(sumM, newScale);
            }
            std::swap(mPrev, mCurr);
            std::swap(fPrev, fCurr);
            std::swap(diagPrev, diagCurr);
            std::swap(dPrevScale, dCurrScale);
        }
        return flags;
    }

    static double lse2(double a, double b) {
        const double m = std::max(a, b);
        if (m == -std::numeric_limits<double>::infinity()) {
            return m;
        }
        return m + std::log(std::exp(a - m) + std::exp(b - m));
    }

    static double lse3(double a, double b, double c) {
        const double m = std::max(a, std::max(b, c));
        if (m == -std::numeric_limits<double>::infinity()) {
            return m;
        }
        return m + std::log(std::exp(a - m) + std::exp(b - m) + std::exp(c - m));
    }

    // 2^e for integer e: 0 below 2^-126, clamped at 2^127
    static float pow2f(int e) {
        e = std::max(-127, std::min(e, 127));
        const uint32_t bits = (uint32_t)(e + 127) << 23;
        float f;
        memcpy(&f, &bits, sizeof(f));
        return f;
    }

    // 2^e for integer e: 0 below 2^-1022, clamped at 2^1023
    static double pow2d(int e) {
        if (e < -1022) {
            return 0.0;
        }
        e = std::min(e, 1023);
        const uint64_t bits = (uint64_t)(e + 1023) << 52;
        double d;
        memcpy(&d, &bits, sizeof(d));
        return d;
    }

    // pow2f for e <= 0 (one clamp)
    static float pow2fNonPos(int e) {
        e = std::max(-127, e);
        const uint32_t bits = (uint32_t)(e + 127) << 23;
        float f;
        memcpy(&f, &bits, sizeof(f));
        return f;
    }

    // binary exponent of a normal float
    static int floatExp(float x) {
        uint32_t bits;
        memcpy(&bits, &x, sizeof(bits));
        return (int)((bits >> 23) & 0xff) - 127;
    }

    // y_t = sum_{k<t} x_k * ge^(t-1-k): exclusive prefix sum decayed by ge; geStep[s] = ge^(2^s)
    static simd_float decayedScan(simd_float x, const simd_float *geStep) {
        simd_float y = simdf32_shiftl_lanes<1>(x);
        y = simdf32_fmadd(simdf32_shiftl_lanes<1>(y), geStep[0], y);
        y = simdf32_fmadd(simdf32_shiftl_lanes<2>(y), geStep[1], y);
#if VECSIZE_FLOAT >= 8
        y = simdf32_fmadd(simdf32_shiftl_lanes<4>(y), geStep[2], y);
#endif
#if VECSIZE_FLOAT == 16
        y = simdf32_fmadd(simdf32_shiftl_lanes<8>(y), geStep[3], y);
#endif
        return y;
    }

    // loadCols/storeCols, or a plain vector load/store inside the matrix (edge = false)
    template<bool edge>
    static simd_float loadColsIn(const float *row, int p, int cols, float pad) {
        return edge ? loadCols(row, p, cols, pad) : simdf32_loadu(row + p);
    }
    template<bool edge>
    static void storeColsIn(float *row, int p, int cols, simd_float a) {
        if (edge) {
            storeCols(row, p, cols, a);
        } else {
            simdf32_storeu(row + p, a);
        }
    }

    // stores the lanes of a before cols to row[p, ...)
    static void storeCols(float *row, int p, int cols, simd_float a) {
        if (p + VECSIZE_FLOAT <= cols) {
            simdf32_storeu(row + p, a);
            return;
        }
        if (p >= cols) {
            return;
        }
        float t[VECSIZE_FLOAT];
        simdf32_storeu(t, a);
        // one test per lane: GCC turns a copy loop into a memcpy call, and a call inside pass()'s block loop
        // makes it save and restore the vector registers around it in every block
        for (int k = 0; k < VECSIZE_FLOAT; ++k) {
            if (p + k < cols) {
                row[p + k] = t[k];
            }
        }
    }

    // buffers of run(), kept across calls and freed with the aligner
    std::vector<float> prevM, prevF;         // previous row's mantissas [colsPadded]
    std::vector<float> diagRow[2];           // M + E + F of the previous and the current row [colsPadded + 1]
    std::vector<int> prevScale;              // previous row's scale per block
    std::vector<int> fwdScale;               // forward ZM scale per (row, block)
    std::vector<float> expScore;             // exp(score) of the forward pass [rows * colsPadded], without expScores
    ZSum zFwd, zBwd;
    // buffers of runScaledDouble() and runDouble(), freed after each fallback (freeFallbackBuffers)
    std::vector<double> dM, dRow[6];         // runDouble: log-space M, E/F rows; runScaledDouble: ZM mantissas,
                                             // M/F/diag rows [cols + 1]
    std::vector<double> dExpScore, dPred;    // runScaledDouble: exp(score) [rows * cols], pred of the row [cols]
    std::vector<int> dFwdScale;              // runScaledDouble: forward ZM scale per (row, block)
    std::vector<int> dPrevScale, dCurrScale; // runScaledDouble: scale per block of the previous and the current row
};

// The matrix part of Impl::computeBacktrace<mode> (MAC), VECSIZE_FLOAT cells at a time, with the same floats,
// states and end cell as its cell-by-cell loop (rows, cols >= 1). Returns false, leaving the result to
// that loop, when P holds a NaN (log-space posteriors at a low temperature; max() would drop it) or when
// mact / 2 is not a float:
//   S_ij = max(P_ij - mact, S_i-1,j-1 + P_ij - mact, S_i,j-1 - mact/2, S_i-1,j - mact/2)
// Only the left term depends on the same row: S_j = max(T_j, S_j-1 - h) is a running maximum that loses
// h = mact/2 per column. Subtracting h commutes with max exactly (both are monotone), so a log-step scan
// and a carry that subtracts h once per lane give the floats of the loop. (The loop subtracts 0.5 * mact
// in double and rounds to float; double has more than 2 * 24 + 2 bits, so that equals the float subtraction.)
class MacFill {
public:
    template<int mode>
    bool run(float **P, int rows, int cols, float mact, uint8_t **bt, size_t &maxI, size_t &maxJ, float &score) {
        if (0.5f * mact * 2.0f != mact) {
            return false;
        }
        const float NI = -std::numeric_limits<float>::infinity();
        const float FMAX = std::numeric_limits<float>::max();
        const int colsVec = (cols + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT * VECSIZE_FLOAT;
        // S of the previous and current row at [j], j = 0..cols, -inf behind
        for (int r = 0; r < 2; ++r) {
            if (sRow[r].size() < (size_t)colsVec + VECSIZE_FLOAT + 1) {
                sRow[r].resize(colsVec + VECSIZE_FLOAT + 1);
            }
        }
        float *prev = sRow[0].data();
        float *curr = sRow[1].data();
        std::fill(prev, prev + colsVec + VECSIZE_FLOAT + 1, NI);
        std::fill(curr, curr + colsVec + VECSIZE_FLOAT + 1, NI);
        prev[0] = 0.0f;
        std::fill(prev + 1, prev + cols + 1, (mode == 3) ? -FMAX : 0.0f);

        float laneIdx[VECSIZE_FLOAT];
        for (int k = 0; k < VECSIZE_FLOAT; ++k) {
            laneIdx[k] = (float)k;
        }
        const simd_float vLane = simdf32_loadu(laneIdx);
        const simd_float vZero = simdf32_setzero();
        const simd_float vNI = simdf32_set(NI);
        const simd_float vMact = simdf32_set(mact);
        const simd_float vH = simdf32_set(0.5f * mact);
        const simd_float lt1 = simdf32_lt(vLane, simdf32_set(1.0f));
        const simd_float lt2 = simdf32_lt(vLane, simdf32_set(2.0f));
        const simd_float lt4 = simdf32_lt(vLane, simdf32_set(4.0f));
        const simd_float lt8 = simdf32_lt(vLane, simdf32_set(8.0f));
        // hFrom[k]: h in lanes >= k, 0 below (x - 0 is exact)
        simd_float hFrom[VECSIZE_FLOAT];
        for (int k = 0; k < VECSIZE_FLOAT; ++k) {
            hFrom[k] = simdf32_blendv_ps(vH, vZero, simdf32_lt(vLane, simdf32_set((float)k)));
        }
        // hSteps[k] = (k + 1) h: what the carry loses on the way to lane k; hSpan = VECSIZE_FLOAT h, the largest
        float hSteps[VECSIZE_FLOAT];
        for (int k = 0; k < VECSIZE_FLOAT; ++k) {
            hSteps[k] = (float)(k + 1) * (0.5f * mact);
        }
        const simd_float vHSteps = simdf32_loadu(hSteps);
        const float hSpan = (float)VECSIZE_FLOAT * 0.5f * mact;
        const simd_float stM = simdf32_set((float)States::M);
        const simd_float stStop = simdf32_set((float)States::STOP);
        const simd_float stI = simdf32_set((float)States::I);
        const simd_float stD = simdf32_set((float)States::D);

        float best = -FMAX;
        size_t bestI = 0, bestJ = 0;
        float lastColAbove = (mode == 3) ? -FMAX : 0.0f; // S[i-1][cols]
        float st[VECSIZE_FLOAT];
        for (int i = 1; i <= rows; ++i) {
            const float left = (mode == 3 && (i & 1)) ? -FMAX : 0.0f;
            curr[0] = left;
            float carry = left; // S[i][j0 - 1]
            simd_float rowMax = vNI;
            simd_float nan = vZero;
            const float *pRow = P[i - 1];
            uint8_t *btRow = bt[i];
            for (int j0 = 1; j0 <= cols; j0 += VECSIZE_FLOAT) {
                const simd_float p = loadCols(pRow, j0 - 1, cols, 0.0f);
                nan = simdf32_or(nan, simdf32_isnan(p));
                const simd_float t1 = simdf32_sub(p, vMact);                                               // STOP
                const simd_float t2 = simdf32_sub(simdf32_add(simdf32_loadu(prev + j0 - 1), p), vMact);  // M
                const simd_float t4 = simdf32_sub(simdf32_loadu(prev + j0), vH);                           // D
                // running maximum of T - k*h inside the vector, then the carry from the left
                simd_float y = simdf32_max(simdf32_max(t1, t2), t4);
                y = simdf32_max(y, subH(simdf32_blendv_ps(simdf32_shiftl_lanes<1>(y), vNI, lt1), vH, 1));
                y = simdf32_max(y, subH(simdf32_blendv_ps(simdf32_shiftl_lanes<2>(y), vNI, lt2), vH, 2));
#if VECSIZE_FLOAT >= 8
                y = simdf32_max(y, subH(simdf32_blendv_ps(simdf32_shiftl_lanes<4>(y), vNI, lt4), vH, 4));
#else
                (void)lt4;
#endif
#if VECSIZE_FLOAT == 16
                y = simdf32_max(y, subH(simdf32_blendv_ps(simdf32_shiftl_lanes<8>(y), vNI, lt8), vH, 8));
#else
                (void)lt8;
#endif
                // The carry wins only in lanes where carry - h, subtracted k + 1 times, is above y. That chain
                // rounds k + 1 times, by less than 2^-19 (|carry| + hSpan) in total; if y exceeds carry - (k + 1) h
                // by twice that (slack) in every lane, S = y exactly and the serial chain is skipped.
                simd_float S;
                const float slack = (std::fabs(carry) + hSpan) * pow2c(-18);
                const simd_float bound = simdf32_add(simdf32_sub(simdf32_set(carry), vHSteps), simdf32_set(slack));
                if (simdf32_movemask_ps(simdf32_gt(y, bound)) == ALL_LANES) {
                    S = y;
                } else {
                    simd_float c = simdf32_set(carry);
                    for (int k = 0; k < VECSIZE_FLOAT; ++k) {
                        c = simdf32_sub(c, hFrom[k]); // lane k: carry - h, k + 1 times
                    }
                    S = simdf32_max(y, c);
                }
                if (j0 + VECSIZE_FLOAT - 1 > cols) {
                    S = simdf32_blendv_ps(vNI, S, simdf32_lt(vLane, simdf32_set((float)(cols + 1 - j0))));
                }
                simdf32_storeu(curr + j0, S);
                // states in the order of calculate_max4
                const simd_float t3 = simdf32_sub(
                    simdf32_blendv_ps(simdf32_shiftl_lanes<1>(S), simdf32_set(carry), lt1), vH); // I
                // simdf32_max(a, b) is a > b ? a : b, the value the blends on c12 / c3 select (ties: b)
                simd_float c12 = simdf32_gt(t1, t2);
                simd_float mx = simdf32_max(t1, t2);
                simd_float state = simdf32_blendv_ps(stM, stStop, c12);
                const simd_float c3 = simdf32_gt(t3, mx);
                mx = simdf32_max(t3, mx);
                state = simdf32_blendv_ps(state, stI, c3);
                state = simdf32_blendv_ps(state, stD, simdf32_gt(t4, mx));
                if (j0 + VECSIZE_FLOAT - 1 <= cols) {
                    simdf32_storeu_u8(btRow + j0, state);
                } else {
                    simdf32_storeu(st, state);
                    for (int k = 0; j0 + k <= cols; ++k) {
                        btRow[j0 + k] = (uint8_t)st[k];
                    }
                }
                rowMax = simdf32_max(rowMax, S);
                carry = simdf32_first(simdf32_broadcast_last(S));
            }
            if (simdf32_movemask_ps(nan) != 0) {
                return false;
            }
            // end cell: the first strict improvement in row-major order
            if (mode == 1) {
                const float m = simdf32_hmax(rowMax);
                if (m > best) {
                    int j = 1;
                    while (curr[j] != m) {
                        ++j;
                    }
                    best = m;
                    bestI = i;
                    bestJ = j;
                }
            } else if (mode == 2) {
                for (int j = (i == rows) ? 1 : cols; j <= cols; ++j) {
                    if (curr[j] > best) {
                        best = curr[j];
                        bestI = i;
                        bestJ = j;
                    }
                }
            } else if (i < rows) {
                lastColAbove = curr[cols];
            }
            std::swap(prev, curr);
        }
        if (mode == 3) {
            // as computeBacktrace: the row buffer it reads after the last swap holds row rows - 1
            bestI = rows;
            bestJ = cols;
            best = lastColAbove;
        }
        maxI = bestI;
        maxJ = bestJ;
        score = best;
        return true;
    }

private:
    // a - h, n times one after the other (each step rounds, as in a cell-by-cell loop)
    static simd_float subH(simd_float a, simd_float h, int n) {
        for (int k = 0; k < n; ++k) {
            a = simdf32_sub(a, h);
        }
        return a;
    }

    std::vector<float> sRow[2]; // S of the previous and the current row
};

// FwBwAligner Constructor for general case: use profile scoring matrix

class Impl final : public FwBwAligner::Impl {
public:
    // Constructor for the general case: score from the substitution matrix profile
    Impl(SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float temperature, float mact,
         size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace);
    // Constructor for a user-defined scoring matrix
    Impl(float gapOpen, float gapExtend, float temperature, float mact,
         size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace);
    ~Impl();

    void reallocateProfile(size_t newColsCapacity);
    void resizeMatrix(bool profile, bool backtrace, size_t newRowsCapacity, size_t newColsCapacity);
    void resetParams(float newGapOpen, float newGapExtend, float newTemperature);
    void setTemperature(float newTemperature) { temperature = newTemperature; }
    void setScaledMode(FwBwAligner::ScaledMode mode) { scaledMode = mode; }
    void initProfile(unsigned char* colAANum, size_t colAALen);
    void initAlignment(unsigned char* targetNum, size_t targetLen, size_t queryLen);
    void initScoreMatrix(float** inputScoreMatrix, int* gaps);
    void runFwBw(bool profile, int backtraceMode);
    FwBwAligner::s_align getFwbwAlnResult() const;

    float getProbability(int i, int j) const { return P[i][j]; }
    float** getProbabiltyMatrix() { return P; }
    float getMaxP() const { return maxP; }
    float getTemperature() const { return temperature; }
    size_t getRowsCapacity() const { return rowsCapacity; }
    size_t getColsCapacity() const { return colsCapacity; }
    size_t getBlockLength() const { return length; }

    template<bool profile, int backtrace>
    void resizeMatrix(size_t newRowsCapacity, size_t newColsCapacity);

    template<bool profile, int backtrace>
    void runFwBw();

private:
    template<bool profile>
    void forward();

    template<bool profile>
    void backward();

    template<bool profile>
    void computeProbabilityMatrix();

    template<int backtrace>
    void computeBacktrace();

    bool runScaled(bool profile);

    unsigned char* colSeqAANum = nullptr;
    unsigned char* rowSeqAANum = nullptr;

    float maxP;
    float temperature;

    size_t length;
    float** zm;
    float** P;
    float* zmFirst;
    float* zeFirst;
    float* zfFirst;

    float** zInit;
    // Score Matrix
    float** scoreForwardProfile = nullptr; // profile true
    float** scoreForwardProfile_exp = nullptr; // profile true
    float** scoreBackwardProfile_exp = nullptr; // profile true

    float** scoreForward = nullptr; // profile false
    //backtrace
    uint8_t** btMatrix = nullptr; // backtrace true

    float* zmBlockPrev;
    float* zmBlockCurr;
    float* zeBlock;
    float* zfBlock;

    float* vj;
    float* wj;
    float* exp_ge_arr;

    float gapOpen;
    float gapExtend;
    float mact; // backtrace true

    size_t rowsCapacity;
    size_t colsCapacity;

    size_t blockCapacity;
    size_t blocks;
    size_t colSeqLen;
    size_t rowSeqLen;

    float** blosum = nullptr; // Profile true
    float* S_prev = nullptr; // backtrace true
    float* S_curr = nullptr; // backtrace true
    float exp_go;
    float exp_ge;
    float max_zm;
    float sum_exp;

    size_t colSeqLen_padding;

    FwBwAligner::s_align alignResult;

    FwBwAligner::ScaledMode scaledMode = FwBwAligner::SCALED_OFF;
    ScaledFwBw<16> scaled;        // score matrix
    ScaledFwBw<32> scaledProfile; // profile
    std::vector<float*> scaledRows; // profile mode: score row, then exp(score) row, of each target residue
    MacFill macFill;              // computeBacktrace: the MAC matrix VECSIZE_FLOAT cells at a time
};

Impl::Impl(SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float temperature, float mact, size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace)
                : temperature(temperature), length(length), gapOpen(gapOpen), gapExtend(gapExtend), mact(mact), rowsCapacity(rowsCapacity), colsCapacity(colsCapacity) {    
    blockCapacity = colsCapacity / length;
    // ZM
    zm = malloc_matrix<float>(rowsCapacity, colsCapacity);
    // Block
    zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zmBlockPrev = (float *) malloc_simd_float((length+1) * sizeof(float));
    zmBlockCurr = (float *) malloc_simd_float((length+1) * sizeof(float));
    zeBlock = (float *) malloc_simd_float((length+1) * sizeof(float));
    zfBlock = (float *) malloc_simd_float((length+1) * sizeof(float));
    // zInit forward & backward
    zInit = malloc_matrix<float>(3, rowsCapacity);
    // Score Matrix (targetProfile 21xcolsCapacity)
    scoreForwardProfile = malloc_matrix<float>(21, colsCapacity);
    scoreForwardProfile_exp = malloc_matrix<float>(21, colsCapacity);
    scoreBackwardProfile_exp = malloc_matrix<float>(21, colsCapacity);
    
    // V,J,exp_ge_arr for ZE
    vj = (float *) malloc_simd_float(length * sizeof(float));
    wj = (float *) malloc_simd_float(length * sizeof(float));
    exp_ge_arr = (float *) malloc_simd_float(length * sizeof(float));

    for (size_t i = 0; i < length; ++i) { 
        vj[i] = exp(((length - 1) * gapExtend + gapOpen - i * gapExtend) / temperature);
        wj[i] = exp(((length - 1) * gapExtend - i * gapExtend) / temperature);
    }
    for (size_t i = 0; i < length; ++i) {
        exp_ge_arr[i] = exp((i * gapExtend + gapExtend) / temperature);
    }
    // Gap open and extend
    exp_go = (static_cast<float>(exp(gapOpen / temperature))); 
    exp_ge = (static_cast<float>(exp(gapExtend / temperature)));
    // Blosum matrix
    blosum = malloc_matrix<float>(21, 21);
    for (int i = 0; i < subMat.alphabetSize; ++i) {
        for (int j = 0; j < subMat.alphabetSize; ++j) {
            blosum[i][j] = static_cast<float>(subMat.subMatrix[i][j]);
        }
    }
    if (backtrace) {
        btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
        S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
        S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    }
}

// FwBwAligner Constructor for user-defined scoring matrix
Impl::Impl(float gapOpen, float gapExtend, float temperature, float mact, size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace)
                    : temperature(temperature), length(length), gapOpen(gapOpen), gapExtend(gapExtend), mact(mact), rowsCapacity(rowsCapacity), colsCapacity(colsCapacity) {
    
    // scoreForward
    scoreForward = malloc_matrix<float>(rowsCapacity, colsCapacity);
    // ZM
    zm = malloc_matrix<float>(rowsCapacity, colsCapacity);
    // Block
    zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    zmBlockPrev = (float *) malloc_simd_float((length+1) * sizeof(float));
    zmBlockCurr = (float *) malloc_simd_float((length+1) * sizeof(float));
    zeBlock = (float *) malloc_simd_float((length+1) * sizeof(float));
    zfBlock = (float *) malloc_simd_float((length+1) * sizeof(float));

    // zInit forward & backward
    zInit = malloc_matrix<float>(3, rowsCapacity);
    
    // V,J,exp_ge_arr for ZE
    vj = (float *) malloc_simd_float(length * sizeof(float));
    wj = (float *) malloc_simd_float(length * sizeof(float));
    exp_ge_arr = (float *) malloc_simd_float(length * sizeof(float));

    for (size_t i = 0; i < length; ++i) { 
        vj[i] = exp(((length - 1) * gapExtend + gapOpen - i * gapExtend) / temperature);
        wj[i] = exp(((length - 1) * gapExtend - i * gapExtend) / temperature);
    }
    for (size_t i = 0; i < length; ++i) {
        exp_ge_arr[i] = exp((i * gapExtend + gapExtend) / temperature);
    }
    // Gap open and extend
    exp_go = (static_cast<float>(exp(gapOpen / temperature))); 
    exp_ge = (static_cast<float>(exp(gapExtend / temperature)));

    if (backtrace != 0) {
        blosum = nullptr;
        btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
        S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
        S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    }
}

Impl::~Impl(){
    // matrices used both in lolalign and general cases
    free(zm);
    free(zmBlockPrev);
    free(zmBlockCurr);
    free(zeBlock);
    free(zfBlock);
    free(zmFirst);
    free(zeFirst);
    free(zfFirst);
    free(vj);
    free(wj);
    free(zInit);
    free(exp_ge_arr);

    // matrices used in only one case
    if (scoreForward != nullptr) {
       free(scoreForward);
    }
    if (scoreForwardProfile != nullptr) {
        free(scoreForwardProfile);
    }
    if (scoreForwardProfile_exp != nullptr) {
        free(scoreForwardProfile_exp);
    }
    if (scoreBackwardProfile_exp != nullptr) {
        free(scoreBackwardProfile_exp);
    }
    if (btMatrix != nullptr) {
        free(btMatrix);
    }
    if (blosum != nullptr) {
        free(blosum);
    }
    if (S_prev != nullptr) {
        free(S_prev);
    }
    if (S_curr != nullptr) {
        free(S_curr);
    }
    
}

//Reallocatation or Resizing
void Impl::reallocateProfile(size_t newColsCapacity) { // reallocate profile when colSeqLen(queryLen in general) exceeds colsCapacity
    free(scoreForwardProfile); scoreForwardProfile = malloc_matrix<float>(21, newColsCapacity);
    free(scoreForwardProfile_exp); scoreForwardProfile_exp = malloc_matrix<float>(21, newColsCapacity);
    free(scoreBackwardProfile_exp); scoreBackwardProfile_exp = malloc_matrix<float>(21, newColsCapacity);
}

template<>
void Impl::resizeMatrix<true,true>(size_t newRowLen, size_t newColLen) {
    //profile = true(scoreForwardProfile, scoreForwardProfile_exp, scoreBackwardProfile_exp)
    //backtrace = true(btMatrix, S_prev, S_curr)

    //check need resizing
    bool resizeRows = newRowLen > rowsCapacity;
    bool resizeCols = newColLen > colsCapacity;
    size_t newRowsCapacity;
    size_t newColsCapacity;
    if (resizeRows || resizeCols) { 
        newRowsCapacity = resizeRows ? ((newRowLen + length - 1) / length) * length : rowsCapacity;
        newColsCapacity = resizeCols ? ((newColLen + length - 1) / length) * length : colsCapacity;
    } else {
        return; // no need to resize
    }

    rowsCapacity = newRowsCapacity;
    colsCapacity = newColsCapacity;
    blockCapacity = newColsCapacity / length;
    free(zm); zm = malloc_matrix<float>(rowsCapacity, colsCapacity);    
    free(zmFirst); zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zeFirst); zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zfFirst); zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zInit); zInit = malloc_matrix<float>(3, rowsCapacity);
    //backtrace
    free(S_prev); S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    free(S_curr); S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    free(btMatrix); btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
}

template<>
void Impl::resizeMatrix<true,false>(size_t newRowLen, size_t newColLen) {
    //profile = true(scoreForwardProfile, scoreForwardProfile_exp, scoreBackwardProfile_exp)
    
    //check need resizing
    bool resizeRows = newRowLen > rowsCapacity;
    bool resizeCols = newColLen > colsCapacity;
    size_t newRowsCapacity;
    size_t newColsCapacity;
    if (resizeRows || resizeCols) { 
        newRowsCapacity = resizeRows ? ((newRowLen + length - 1) / length) * length : rowsCapacity;
        newColsCapacity = resizeCols ? ((newColLen + length - 1) / length) * length : colsCapacity;
    } else {
        return; // no need to resize
    }

    rowsCapacity = newRowsCapacity;
    colsCapacity = newColsCapacity;
    blockCapacity = newColsCapacity / length;
    free(zm); zm = malloc_matrix<float>(rowsCapacity, colsCapacity);    
    free(zmFirst); zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zeFirst); zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zfFirst); zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zInit); zInit = malloc_matrix<float>(3, rowsCapacity);
}

template<>
void Impl::resizeMatrix<false,true>(size_t newRowLen, size_t newColLen) {
    //profile = false(scoreForward)
    //backtrace = true(btMatrix, S_prev, S_curr)

    //check need resizing
    bool resizeRows = newRowLen > rowsCapacity;
    bool resizeCols = newColLen > colsCapacity;
    size_t newRowsCapacity;
    size_t newColsCapacity;
    if (resizeRows || resizeCols) { 
        newRowsCapacity = resizeRows ? ((newRowLen + length - 1) / length) * length : rowsCapacity;
        newColsCapacity = resizeCols ? ((newColLen + length - 1) / length) * length : colsCapacity;
    } else {
        return; // no need to resize
    }

    rowsCapacity = newRowsCapacity;
    colsCapacity = newColsCapacity;
    blockCapacity = newColsCapacity / length;
    free(zm); zm = malloc_matrix<float>(rowsCapacity, colsCapacity);    
    free(zmFirst); zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zeFirst); zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zfFirst); zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zInit); zInit = malloc_matrix<float>(3, rowsCapacity);
    //profile
    free(scoreForward); scoreForward = malloc_matrix<float>(rowsCapacity, colsCapacity);
    //backtrace
    free(S_prev); S_prev = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    free(S_curr); S_curr = (float *) malloc_simd_float((colsCapacity+1) * sizeof(float));
    free(btMatrix); btMatrix = malloc_matrix<uint8_t>(rowsCapacity + 1, colsCapacity + 1);
}

template<>
void Impl::resizeMatrix<false,false>(size_t newRowLen, size_t newColLen) {
    //profile = false(scoreForward)
    //backtrace = false

    //check need resizing
    bool resizeRows = newRowLen > rowsCapacity;
    bool resizeCols = newColLen > colsCapacity;
    size_t newRowsCapacity;
    size_t newColsCapacity;
    if (resizeRows || resizeCols) { 
        newRowsCapacity = resizeRows ? ((newRowLen + length - 1) / length) * length : rowsCapacity;
        newColsCapacity = resizeCols ? ((newColLen + length - 1) / length) * length : colsCapacity;
    } else {
        return; // no need to resize
    }

    rowsCapacity = newRowsCapacity;
    colsCapacity = newColsCapacity;
    blockCapacity = newColsCapacity / length;
    free(zm); zm = malloc_matrix<float>(rowsCapacity, colsCapacity);    
    free(zmFirst); zmFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zeFirst); zeFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zfFirst); zfFirst = (float *) malloc_simd_float((rowsCapacity+1) * sizeof(float));
    free(zInit); zInit = malloc_matrix<float>(3, rowsCapacity);
    //profile
    free(scoreForward); scoreForward = malloc_matrix<float>(rowsCapacity, colsCapacity);
}

void Impl::resetParams(float newGapOpen, float newGapExtend, float newTemperature) {
    gapOpen = newGapOpen;
    gapExtend = newGapExtend;
    temperature = newTemperature;
    exp_go = (static_cast<float>(exp(gapOpen / temperature))); 
    exp_ge = (static_cast<float>(exp(gapExtend / temperature)));

    for (size_t i = 0; i < length; ++i) { 
        vj[i] = exp(((length - 1) * gapExtend + gapOpen - i * gapExtend) / temperature);
        wj[i] = exp(((length - 1) * gapExtend - i * gapExtend) / temperature);
    }
    for (size_t i = 0; i < length; ++i) {
        exp_ge_arr[i] = exp((i * gapExtend + gapExtend) / temperature);
    }
}


void Impl::initAlignment(unsigned char* targetAANum, size_t targetLen, size_t queryLen) {
    rowSeqAANum = targetAANum; rowSeqLen = targetLen;
    resizeMatrix<true, true>(targetLen, queryLen);
}

void Impl::initScoreMatrix(float** inputScoreMatrix, int* gaps) {
    //gaps: [rowStart, rowEnd, colStart, colEnd]
    rowSeqLen = gaps[1]-gaps[0]; colSeqLen = gaps[3]-gaps[2]; //row=tlen, col=qlen
    colSeqLen_padding = ((colSeqLen + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT; //colpadding
    blocks = (colSeqLen / length) + (colSeqLen % length != 0);
    simd_float vTemp = simdf32_set(temperature);
    if (colSeqLen > colsCapacity) {
        size_t newColsCapacity = ((colSeqLen + length-1)/length)* length;
        free(scoreForward); scoreForward = malloc_matrix<float>(rowsCapacity, newColsCapacity);
    }
    size_t colLoopCount = colSeqLen/VECSIZE_FLOAT;
    size_t colEndPos = colLoopCount*VECSIZE_FLOAT;
    for (size_t i = 0; i < rowSeqLen; ++i){
        for (size_t j = 0; j < colEndPos; j+=VECSIZE_FLOAT) {
            simd_float vScoreForward = simdf32_loadu(&inputScoreMatrix[i+gaps[0]][j+gaps[2]]);
            vScoreForward = simdf32_div(vScoreForward, vTemp);
            simdf32_store(&scoreForward[i][j], vScoreForward);
        }
        for(size_t j = colEndPos; j < colSeqLen; ++j){
            scoreForward[i][j] = inputScoreMatrix[i+gaps[0]][j+gaps[2]] / temperature;

        }
        for (size_t j = colSeqLen; j < colSeqLen_padding; ++j){
            scoreForward[i][j] = FLT_MIN_EXP;
        } 
    }
}

void Impl::initProfile(unsigned char* colAANum, size_t colAALen) {
    colSeqAANum = colAANum; colSeqLen = colAALen;
    colSeqLen_padding = ((colSeqLen + VECSIZE_FLOAT -1) / VECSIZE_FLOAT) * VECSIZE_FLOAT;
    blocks = (colSeqLen / length) + (colSeqLen % length != 0);
    if (colSeqLen > colsCapacity) {
        size_t newColsCapacity = ((colSeqLen + length-1)/length)* length;
        reallocateProfile(newColsCapacity);
    }
    
    //scoreForward : 21 * qlen
    for (size_t i=0; i<21; ++i){
        for (size_t j=0; j < colAALen; ++j) {
            float score = blosum[i][colSeqAANum[j]]/temperature;
            scoreForwardProfile[i][j] = score;
        }
        std::fill(&scoreForwardProfile[i][colSeqLen], &scoreForwardProfile[i][colSeqLen_padding], FLT_MIN_EXP);
    }
    for (size_t i=0; i<21; ++i){
        for (size_t j=0; j < colSeqLen_padding; j += VECSIZE_FLOAT) {
            simd_float vScoreForward = simdf32_load(&scoreForwardProfile[i][j]);
            vScoreForward = simdf32_exp(vScoreForward);
            simdf32_store(&scoreForwardProfile_exp[i][j], vScoreForward);
        }
        for (size_t j=0; j < colSeqLen; ++j){
            size_t reverse_j = colSeqLen - 1 - j;
            scoreBackwardProfile_exp[i][reverse_j] = scoreForwardProfile_exp[i][j];
        }
        //remainder 
        for (size_t j=colSeqLen; j < colSeqLen_padding; ++j){
            scoreBackwardProfile_exp[i][j] = 0;
        }
    }
}

template <bool profile>
void Impl::forward() {
    //Init zInit
    for (size_t i = 0 ; i < 3; ++i) {
        std::fill(zInit[i], zInit[i] + rowsCapacity, FLT_MIN_EXP); // rowsCapacity -> tlen
    }  
    max_zm = -std::numeric_limits<float>::max(); 
    simd_float vMax_zm = simdf32_set(max_zm);
    P = nullptr; // reset p. do we need this?
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, colSeqLen) - start;
        size_t cols = length;
        if (memcpy_cols != length) {
            cols = ((memcpy_cols + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT; //padding vecsize_float
        }
        //Init blocks
        memset(zmBlockPrev, 0, (length + 1) * sizeof(float));
        memset(zeBlock, 0, (length + 1) * sizeof(float));
        memset(zfBlock, 0, (length + 1) * sizeof(float));
            
        memcpy(zmFirst + 1, zInit[0], rowSeqLen * sizeof(float));
        memcpy(zeFirst + 1, zInit[1], rowSeqLen * sizeof(float));
        memcpy(zfFirst + 1, zInit[2], rowSeqLen * sizeof(float));

        //Init initial values
        zmBlockPrev[0] = 0; zfBlock[0] = 0; zeBlock[0] = 0;
        zmBlockCurr[0] = exp(zmFirst[1]);
        float ze_i0 = expf(zeFirst[1]);
        float current_max = 0;
        float zmMaxRowBlock = -std::numeric_limits<float>::max();
        float log_zmMax = 0;
        simd_float vZmMaxRowBlock;
        for (size_t i = 1; i <= rowSeqLen; ++i) {
            simd_float vExpMax = simdf32_set(exp(-current_max));
            simd_float vZeI0 = simdf32_set(ze_i0);
            simd_float vLastPrefixSum = simdf32_setzero(); 
            simd_float vZmax_tmp = simdf32_set(-std::numeric_limits<float>::max());
            // ZM calculation
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_load(&zmBlockPrev[j-1]);
                simd_float vZe = simdf32_load(&zeBlock[j-1]);
                simd_float vZf = simdf32_load(&zfBlock[j-1]);
                simd_float vScoreMatrix;
                if (profile) {
                    vScoreMatrix = simdf32_load(&scoreForwardProfile_exp[rowSeqAANum[i-1]][start + j - 1]);
                } else {
                    vScoreMatrix = simdf32_exp(simdf32_load(&scoreForward[i-1][start + j - 1]));
                } 
                simd_float vZmCurrUpdate = simdf32_add(simdf32_add(vZmPrev, vZe), simdf32_add(vZf, vExpMax));
                vZmCurrUpdate = simdf32_mul(vZmCurrUpdate, vScoreMatrix);
                vZmax_tmp = simdf32_max(vZmax_tmp, vZmCurrUpdate);
                simdf32_storeu(&zmBlockCurr[j], vZmCurrUpdate);
            }
            zmMaxRowBlock = simdf32_hmax(vZmax_tmp);
            vZmMaxRowBlock = simdf32_set(zmMaxRowBlock);

            // ZF calculation 
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_loadu(&zmBlockPrev[j]);
                simd_float vZf = simdf32_loadu(&zfBlock[j]);
                simd_float vZfUpdate = simdf32_add(
                                        simdf32_mul(vZmPrev, simdf32_set(exp_go)),
                                        simdf32_mul(vZf, simdf32_set(exp_ge))
                                        );
                vZfUpdate = simdf32_div(vZfUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zfBlock[j], vZfUpdate);
            }
            for (size_t j = 0; j < cols; j += VECSIZE_FLOAT) { 
                simd_float vZmCurr = simdf32_load(&zmBlockCurr[j]);
                simd_float vVj = simdf32_load(&vj[j]);
                simd_float vCumsumZm = simdf32_mul(vZmCurr, vVj);
                vCumsumZm = simdf32_prefixsum(vCumsumZm);
                vCumsumZm = simdf32_add(vCumsumZm, vLastPrefixSum);
                vLastPrefixSum = simdf32_set(vCumsumZm[(VECSIZE_FLOAT - 1)]);
                simd_float vWj = simdf32_load(&wj[j]);
                simd_float vExp_ge_arr = simdf32_load(&exp_ge_arr[j]);
                simd_float vZeUpdate = simdf32_add(
                                        simdf32_div(vCumsumZm, vWj),
                                        simdf32_mul(vZeI0, vExp_ge_arr)
                                        );
                // simd_float vZeUpdate = simdf32_fmadd(vZeI0, vExp_ge_arr, simdf32_div(vCumsumZm, vWj));
                vZeUpdate = simdf32_div(vZeUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zeBlock[j+1], vZeUpdate);
            }

            log_zmMax = log(zmMaxRowBlock);
            current_max += log_zmMax;
            simd_float vCurrMax = simdf32_set(current_max);
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT){
                simd_float vZmCurr = simdf32_loadu(&zmBlockCurr[j]);
                vZmCurr = simdf32_div(vZmCurr, vZmMaxRowBlock);
                simdf32_storeu(&zmBlockCurr[j], vZmCurr);
                vZmCurr = simdf32_add(simdf32_log(vZmCurr), vCurrMax);
                vMax_zm = simdf32_max(vMax_zm, vZmCurr);
                simdf32_store(&zm[i - 1][start + j-1], vZmCurr);
            }     

            #if defined(AVX512)
                simd_float vNextZinit = _mm512_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
            #elif defined(AVX2)
                simd_float vNextZinit = _mm256_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
            #else // Fallback to SSE
                simd_float vNextZinit = _mm_set_ps(
                    1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
            #endif


            vNextZinit = simdf32_log(vNextZinit);
            vNextZinit = simdf32_add(vNextZinit, vCurrMax);

            zInit[0][i-1] = vNextZinit[0]; zInit[1][i-1] = vNextZinit[1]; zInit[2][i-1] = vNextZinit[2];
            std::swap(zmBlockCurr, zmBlockPrev);
            
            if (i < rowSeqLen) {
                zmFirst[i+1] -= current_max;
                zeFirst[i+1] -= current_max;

#if defined(AVX512)
                simd_float vNextFirstExp = _mm512_set_ps(
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max,
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0];
                ze_i0 = vNextFirstExp[1];
                zmBlockPrev[0] = vNextFirstExp[2];
                zeBlock[0] = vNextFirstExp[3];
                zfBlock[0] = vNextFirstExp[4];
#elif defined(AVX2)
                simd_float vNextFirstExp = _mm256_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );        
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1]; 
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#else // Fallback to SSE
                simd_float vNextFirstExp1 = _mm_set_ps(
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );    
                simd_float vNextFirstExp2= _mm_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max    
                );    
                vNextFirstExp1 = simdf32_exp(vNextFirstExp1);
                vNextFirstExp2 = simdf32_exp(vNextFirstExp2);
                zmBlockCurr[0] = vNextFirstExp1[0]; ze_i0 = vNextFirstExp1[1]; 
                zmBlockPrev[0] = vNextFirstExp1[2]; zeBlock[0] = vNextFirstExp1[3];
                zfBlock[0] = vNextFirstExp2[0];
#endif 
            } else{
                zmBlockPrev[0] = exp(zmFirst[i] - log_zmMax);
                zeBlock[0] = exp(zeFirst[i] - log_zmMax); 
                zfBlock[0] = exp(zfFirst[i] - current_max); 
            }
        }  
    }

    sum_exp= 0.0; 
    simd_float vSum_exp = simdf32_setzero();

    //Calculate max_zm
    for (size_t i = 0; i < VECSIZE_FLOAT; ++i) {
        max_zm = std::max(max_zm, vMax_zm[i]);
    }
    vMax_zm = simdf32_set(max_zm);

    //Calculate sum_exp
    for (size_t i = 0; i < rowSeqLen; ++i) {
        //Removed remainder handling. Check needed
        for (size_t j = 0; j < colSeqLen_padding; j+= VECSIZE_FLOAT) {
            simd_float vZmForward = simdf32_load(&zm[i][j]);
            vZmForward = simdf32_exp(simdf32_sub(vZmForward, vMax_zm));
            vSum_exp = simdf32_add(vSum_exp, vZmForward);
        }
    }
    sum_exp += simdf32_hadd(vSum_exp);
}

template <bool profile>
void Impl::backward()  {
    //Init zInit
    size_t vecsize_float = static_cast<size_t>(VECSIZE_FLOAT);
    for (size_t i = 0 ; i < 3; ++i) {
        std::fill(zInit[i], zInit[i] + rowsCapacity, FLT_MIN_EXP); // rowsCapacity -> tlen
    }  
    for (size_t b = 0; b < blocks; ++b) {
        size_t start = b * length;
        size_t end = (b + 1) * length;
        size_t memcpy_cols = std::min(end, colSeqLen) - start;
        size_t cols = length;
        if (memcpy_cols != length) {
            cols = ((memcpy_cols + VECSIZE_FLOAT - 1) / VECSIZE_FLOAT) * VECSIZE_FLOAT; //padding vecsize_float
        }
        //Init blocks
        memset(zmBlockPrev, 0, (length + 1) * sizeof(float));
        memset(zeBlock, 0, (length + 1) * sizeof(float));
        memset(zfBlock, 0, (length + 1) * sizeof(float));

        memcpy(zmFirst + 1, zInit[0], rowSeqLen * sizeof(float));
        memcpy(zeFirst + 1, zInit[1], rowSeqLen * sizeof(float));
        memcpy(zfFirst + 1, zInit[2], rowSeqLen * sizeof(float));

        //Init initial values
        zmBlockPrev[0] = 0; zfBlock[0] = 0; zeBlock[0] = 0;
        zmBlockCurr[0] = exp(zmFirst[1]);
        float ze_i0 = expf(zeFirst[1]);
        float current_max = 0;
        float zmMaxRowBlock = -std::numeric_limits<float>::max();
        float log_zmMax = 0;
        simd_float vZmMaxRowBlock;
        for (size_t i = 1; i <= rowSeqLen; ++i) {
            simd_float vExpMax = simdf32_set(exp(-current_max));
            simd_float vZeI0 = simdf32_set(ze_i0);
            simd_float vLastPrefixSum = simdf32_setzero(); 
            simd_float vZmax_tmp = simdf32_set(-std::numeric_limits<float>::max());
            // ZM calculation
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_load(&zmBlockPrev[j-1]);
                simd_float vZe = simdf32_load(&zeBlock[j-1]);
                simd_float vZf = simdf32_load(&zfBlock[j-1]);
                simd_float vScoreMatrix;
                if (profile) {
                    vScoreMatrix = simdf32_load(&scoreBackwardProfile_exp[rowSeqAANum[rowSeqLen - i]][start + j - 1]);
                } else {
                    size_t reverse_i = rowSeqLen - i;
                    size_t reverse_j = colSeqLen - start - j + 1;
                    simd_float vScoreBackward;
                    if (reverse_j >= vecsize_float) {
                        vScoreBackward = simdf32_loadu(&scoreForward[reverse_i][reverse_j - vecsize_float]);
                    } else {
                        size_t elements_to_fill = reverse_j;
                        vScoreBackward = simdf32_set(FLT_MIN_EXP);
                        // fill from the back
                        for (size_t k = 0; k < elements_to_fill; ++k) {
                            vScoreBackward[VECSIZE_FLOAT - elements_to_fill + k] = scoreForward[reverse_i][k];
                        }
                    }
                    vScoreBackward = simdf32_reverse(vScoreBackward);
                    vScoreMatrix = simdf32_exp(vScoreBackward);
                }
                simd_float vZmCurrUpdate = simdf32_add(simdf32_add(vZmPrev, vZe), simdf32_add(vZf, vExpMax));
                vZmCurrUpdate = simdf32_mul(vZmCurrUpdate, vScoreMatrix);
                vZmax_tmp = simdf32_max(vZmax_tmp, vZmCurrUpdate);
                simdf32_storeu(&zmBlockCurr[j], vZmCurrUpdate);
            }
            zmMaxRowBlock = simdf32_hmax(vZmax_tmp);
            vZmMaxRowBlock = simdf32_set(zmMaxRowBlock);

            // ZF calculation 
            for (size_t j = 1; j <= cols; j += VECSIZE_FLOAT) {
                simd_float vZmPrev = simdf32_loadu(&zmBlockPrev[j]);
                simd_float vZf = simdf32_loadu(&zfBlock[j]);
                simd_float vZfUpdate = simdf32_add(
                                        simdf32_mul(vZmPrev, simdf32_set(exp_go)),
                                        simdf32_mul(vZf, simdf32_set(exp_ge))
                                        );
                vZfUpdate = simdf32_div(vZfUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zfBlock[j], vZfUpdate);
            }
            for (size_t j = 0; j < cols; j += VECSIZE_FLOAT) { 
                simd_float vZmCurr = simdf32_load(&zmBlockCurr[j]);
                simd_float vVj = simdf32_load(&vj[j]);
                simd_float vCumsumZm = simdf32_mul(vZmCurr, vVj);
                vCumsumZm = simdf32_prefixsum(vCumsumZm);
                vCumsumZm = simdf32_add(vCumsumZm, vLastPrefixSum);
                vLastPrefixSum = simdf32_set(vCumsumZm[(VECSIZE_FLOAT - 1)]);
                simd_float vWj = simdf32_load(&wj[j]);
                simd_float vExp_ge_arr = simdf32_load(&exp_ge_arr[j]);
                simd_float vZeUpdate = simdf32_add(
                                        simdf32_div(vCumsumZm, vWj),
                                        simdf32_mul(vZeI0, vExp_ge_arr)
                                        );
                vZeUpdate = simdf32_div(vZeUpdate, vZmMaxRowBlock);
                simdf32_storeu(&zeBlock[j+1], vZeUpdate);
            }

            log_zmMax = log(zmMaxRowBlock);
            current_max += log_zmMax;
            size_t adjusted_memcpycols = memcpy_cols - memcpy_cols % VECSIZE_FLOAT;
            size_t forwardBlockStart = colSeqLen - start;
            simd_float vCurrMax = simdf32_set(current_max);
            for (size_t j = 1; j <= adjusted_memcpycols; j += VECSIZE_FLOAT) {
                forwardBlockStart -= vecsize_float;
                size_t simd_index = forwardBlockStart;
                simd_float vZmCurr = simdf32_loadu(&zmBlockCurr[j]);
                vZmCurr = simdf32_div(vZmCurr, vZmMaxRowBlock);
                simdf32_storeu(&zmBlockCurr[j], vZmCurr);
                vZmCurr = simdf32_add(simdf32_log(vZmCurr), vCurrMax);

                simd_float vZmForward = simdf32_loadu(&zm[rowSeqLen - i][simd_index]);

                simd_float vZmCurr_reverse = simdf32_reverse(vZmCurr);
                simd_float vZmForward_Backward = simdf32_add(vZmForward, vZmCurr_reverse);
                simdf32_storeu(&zm[rowSeqLen - i][simd_index], vZmForward_Backward);
            }

            // Handle remainder
            if (memcpy_cols != length) {
                size_t remainder = memcpy_cols % VECSIZE_FLOAT;
                simd_float vZmCurr = simdf32_loadu(&zmBlockCurr[adjusted_memcpycols+1]);
                vZmCurr = simdf32_div(vZmCurr, vZmMaxRowBlock);
                simdf32_storeu(&zmBlockCurr[adjusted_memcpycols+1], vZmCurr);
                vZmCurr = simdf32_add(simdf32_log(vZmCurr), simdf32_set(current_max));
                for (size_t k = 0; k < remainder; ++k) {
                    size_t j_index = remainder - k - 1;
                    zm[rowSeqLen - i][j_index] += vZmCurr[k];
                }
            }

#if defined(AVX512)
                simd_float vNextZinit = _mm512_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
#elif defined(AVX2)
                simd_float vNextZinit = _mm256_set_ps(
                    1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
#else // Fallback to SSE
                simd_float vNextZinit = _mm_set_ps(
                    1.0f,
                    zfBlock[memcpy_cols],
                    zeBlock[memcpy_cols],
                    zmBlockCurr[memcpy_cols]
                );
#endif

            vNextZinit = simdf32_log(vNextZinit);
            vNextZinit = simdf32_add(vNextZinit, vCurrMax);

            zInit[0][i-1] = vNextZinit[0];
            zInit[1][i-1] = vNextZinit[1];
            zInit[2][i-1] = vNextZinit[2];
            std::swap(zmBlockCurr, zmBlockPrev);
            
            if (i < rowSeqLen) {
                zmFirst[i+1] -= current_max;
                zeFirst[i+1] -= current_max;

#if defined(AVX512)
                simd_float vNextFirstExp = _mm512_set_ps(
                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0];
                ze_i0 = vNextFirstExp[1];
                zmBlockPrev[0] = vNextFirstExp[2];
                zeBlock[0] = vNextFirstExp[3];
                zfBlock[0] = vNextFirstExp[4];
#elif defined(AVX2)
                simd_float vNextFirstExp = _mm256_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max, 
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );        
                vNextFirstExp = simdf32_exp(vNextFirstExp);
                zmBlockCurr[0] = vNextFirstExp[0]; ze_i0 = vNextFirstExp[1]; 
                zmBlockPrev[0] = vNextFirstExp[2]; zeBlock[0] = vNextFirstExp[3]; zfBlock[0] = vNextFirstExp[4];
#else // Fallback to SSE
                simd_float vNextFirstExp1 = _mm_set_ps(
                    zeFirst[i] - log_zmMax,
                    zmFirst[i] - log_zmMax,
                    zeFirst[i+1],
                    zmFirst[i+1]
                );    
                simd_float vNextFirstExp2= _mm_set_ps(
                    0.0f, 0.0f, 0.0f,
                    zfFirst[i] - current_max     
                );    
                vNextFirstExp1 = simdf32_exp(vNextFirstExp1);
                vNextFirstExp2 = simdf32_exp(vNextFirstExp2);
                zmBlockCurr[0] = vNextFirstExp1[0]; ze_i0 = vNextFirstExp1[1]; 
                zmBlockPrev[0] = vNextFirstExp1[2]; zeBlock[0] = vNextFirstExp1[3];
                zfBlock[0] = vNextFirstExp2[0];
#endif 
            } else{
                zmBlockPrev[0] = exp(zmFirst[i] - log_zmMax);
                zeBlock[0] = exp(zeFirst[i] - log_zmMax); 
                zfBlock[0] = exp(zfFirst[i] - current_max); 
            }
        }  
    }
}

template<>
// profile: true for profile, false for query*target matrix
// backtrace: 0 for no backtrace, 1 for local alignment, 2 for semi-global alignment, 3 for global alignment. Default: <true, 1>
void Impl::runFwBw<true, 0>() {
    forward<true>();
    backward<true>();
    computeProbabilityMatrix<true>();
}

template<>
void Impl::runFwBw<false, 0>() {
    forward<false>();
    backward<false>();
    computeProbabilityMatrix<false>();
}

template<>
void Impl::runFwBw<true, 1>() {
    forward<true>();
    backward<true>();
    computeProbabilityMatrix<true>();
    computeBacktrace<1>();
}

template<>
void Impl::runFwBw<false, 1>() {
    forward<false>();
    backward<false>();
    computeProbabilityMatrix<false>();
    computeBacktrace<1>();
}

template<>
void Impl::runFwBw<true, 2>() {
    forward<true>();
    backward<true>();
    computeProbabilityMatrix<true>();
    computeBacktrace<2>();
}

template<>
void Impl::runFwBw<false, 2>() {
    forward<false>();
    backward<false>();
    computeProbabilityMatrix<false>();
    computeBacktrace<2>();
}

template<>
void Impl::runFwBw<true, 3>() {
    forward<true>();
    backward<true>();
    computeProbabilityMatrix<true>();
    computeBacktrace<3>();
}

template<>
void Impl::runFwBw<false, 3>() {
    forward<false>();
    backward<false>();
    computeProbabilityMatrix<false>();
    computeBacktrace<3>();
}

template<bool profile>
void Impl::computeProbabilityMatrix() {
    float logsumexp_zm = max_zm + log(sum_exp);

    simd_float vLogsumexp_zm = simdf32_set(logsumexp_zm);
    size_t colLoopCount = colSeqLen / VECSIZE_FLOAT; 
    size_t colLoopEndPos = colLoopCount * VECSIZE_FLOAT;

    P = zm; //reuse zm matrix to store the probability matrix
    maxP = 0.0;
    simd_float vMaxP = simdf32_setzero();
    for (size_t i = 0; i < rowSeqLen; ++i) {
        // Fill the probability matrix
        //Removed remainder handling. Check needed
        for (size_t j = 0; j < colLoopEndPos; j += VECSIZE_FLOAT) {
            simd_float vZmForward_Backward = simdf32_load(&zm[i][j]);
            simd_float scoreForwardVal;
            if (profile) {
                scoreForwardVal = simdf32_load(&scoreForwardProfile[rowSeqAANum[i]][j]);
            } else {
                scoreForwardVal = simdf32_load(&scoreForward[i][j]);
            }
            // simd_float scoreForwardVal = simdf32_load(&scoreForward[targetNum[i]][j]);
            simd_float P_val = simdf32_exp(simdf32_sub(vZmForward_Backward, simdf32_add(scoreForwardVal, vLogsumexp_zm)));
            simdf32_store(&P[i][j], P_val);
            vMaxP = simdf32_max(vMaxP, P_val);
        }    
        for (size_t j = colLoopEndPos; j < colSeqLen; ++j) {
            if (profile) {
                P[i][j] = exp(zm[i][j] - scoreForwardProfile[rowSeqAANum[i]][j] - logsumexp_zm);
            } else {
                P[i][j] = exp(zm[i][j] - scoreForward[i][j] - logsumexp_zm);
            }
            // P[i][j] = exp(zm[i][j] - scoreForward[targetNum[i]][j] - logsumexp_zm);
            maxP = std::max(maxP, P[i][j]);
        }
    }

    // Calculate the maximum probability
    for (size_t k = 0; k < VECSIZE_FLOAT; ++k) {
        maxP = std::max(maxP, vMaxP[k]);
    }

}

template<int backtrace>
void Impl::computeBacktrace() {
    // MAC algorithm from HH-suite
    uint8_t val;
    size_t max_i = 0;
    size_t max_j = 0;
    float term1, term2, term3, term4 = 0.0f;
    float score_MAC = -std::numeric_limits<float>::max();

    memset(S_curr, 0, (colSeqLen + 1) * sizeof(float));

    switch (backtrace) {
        case 1: // local
            memset(S_prev, 0, (colSeqLen + 1) * sizeof(float));
            break;
        case 2: // semiglobal
            memset(S_prev, 0, (colSeqLen + 1) * sizeof(float));
            break;
        case 3: // global
            std::fill(S_prev, S_prev + colSeqLen + 1, -std::numeric_limits<float>::max());
            S_prev[0] = 0.0;
            S_curr[0] = -std::numeric_limits<float>::max();
            break;
    }


    for (size_t i = 0; i <= rowSeqLen; ++i) {
        btMatrix[i][0] = States::STOP;
    }
    for (size_t j = 0; j <= colSeqLen; ++j) {
        btMatrix[0][j] = States::STOP;
    }

    // SIMD version of the loop below: same S, states and end cell
    const bool simdFill = rowSeqLen > 0 && colSeqLen > 0 &&
                          macFill.run<backtrace>(P, (int)rowSeqLen, (int)colSeqLen, mact, btMatrix, max_i, max_j, score_MAC);
    for (size_t i = 1; i <= rowSeqLen && !simdFill; ++i) {
        for (size_t j = 1; j <= colSeqLen; ++j) {
            term1 = P[i - 1][j - 1] - mact; // STOP
            term2 = S_prev[j - 1] + P[i - 1][j - 1] - mact; // M
            term4 = S_prev[j] - 0.5 * mact; // D
            term3 = S_curr[j - 1] - 0.5 * mact; // I
            calculate_max4(S_curr[j], term1, term2, term3, term4, val);
            btMatrix[i][j] = val;
            
            switch (backtrace) {
                case 1: // local
                    if (S_curr[j] > score_MAC) {
                        max_i = i;
                        max_j = j;
                        score_MAC = S_curr[j];
                    }
                    break;
                case 2: // semiglobal
                    if ((i == rowSeqLen || j == colSeqLen) && S_curr[j] > score_MAC) { // only calculate for last column if j is last column
                        max_i = i;
                        max_j = j;
                        score_MAC = S_curr[j];
                    }
                    break;
                case 3: // global
                    break;
            }
        }
        std::swap(S_prev, S_curr);
        if (backtrace == 3 && i==1) {
            S_prev[0] = -std::numeric_limits<float>::max();
        }
    }
    // Set max_i and max_j for global alignment
    if (backtrace == 3 && !simdFill) {
        max_i = rowSeqLen;
        max_j = colSeqLen;
        score_MAC = S_curr[colSeqLen];
    }
    // traceback 
    alignResult = FwBwAligner::s_align();
    alignResult.cigar = "";
    alignResult.cigar.reserve(colSeqLen + rowSeqLen);
    alignResult.score1 = maxP;
    alignResult.score2 = score_MAC;

    alignResult.qEndPos1 = max_j - 1;
    alignResult.dbEndPos1 = max_i - 1;
    uint32_t aaIds = 0;
    bool exitLoop = false;
    while (max_i > 0 && max_j > 0 && !exitLoop) {
        uint8_t state = btMatrix[max_i][max_j];
        switch (state) {
            case States::M:
                --max_i;
                --max_j;
                alignResult.qStartPos1 = max_j;
                alignResult.dbStartPos1 = max_i;
                alignResult.cigar.push_back('M');
                aaIds += (rowSeqAANum[max_i] == colSeqAANum[max_j]);
                break;

            case States::I:
                --max_j;
                alignResult.cigar.push_back('I');
                break;

            case States::D:
                --max_i;
                alignResult.cigar.push_back('D');
                break;

            default:
                exitLoop = true; 
                break;
        }
    }
    while (!alignResult.cigar.empty() && alignResult.cigar.back() != 'M') {
        alignResult.cigar.pop_back();
    }
    alignResult.cigarLen = alignResult.cigar.length();
    std::reverse(alignResult.cigar.begin(), alignResult.cigar.end());
    alignResult.identicalAACnt = aaIds;
    alignResult.qCov = SmithWaterman::computeCov(alignResult.qStartPos1, alignResult.qEndPos1, colSeqLen);
    alignResult.dbCov = SmithWaterman::computeCov(alignResult.dbStartPos1, alignResult.dbEndPos1, rowSeqLen);
}

FwBwAligner::s_align Impl::getFwbwAlnResult() const {
    return alignResult;
}

// Runtime-parameter entry points, so that the template instantiation stays inside this
// translation unit and the dispatch interface carries no template arguments.
void Impl::resizeMatrix(bool profile, bool backtrace, size_t newRowsCapacity, size_t newColsCapacity) {
    if (profile) {
        if (backtrace) {
            resizeMatrix<true, true>(newRowsCapacity, newColsCapacity);
        } else {
            resizeMatrix<true, false>(newRowsCapacity, newColsCapacity);
        }
    } else {
        if (backtrace) {
            resizeMatrix<false, true>(newRowsCapacity, newColsCapacity);
        } else {
            resizeMatrix<false, false>(newRowsCapacity, newColsCapacity);
        }
    }
}

// run() first; a call it flags (any FLAG_*, in SCALED also a backward Z that disagrees) is run again with
// double mantissas (runScaledDouble), and one still flagged in log space (runDouble)
// expScore: exp(score) rows, or score itself when expScores is false
template<bool expScores, class Kernel>
static void runScaledKernel(Kernel &k, FwBwAligner::ScaledMode mode, float **score, float **expScore,
                            int rows, int cols, float expGo, float expGe, float **prob, float &maxP) {
    unsigned int flags = (mode == FwBwAligner::SCALED_FAST)
                         ? k.template run<false, expScores>(expScore, rows, cols, expGo, expGe, prob, maxP)
                         : k.template run<true, expScores>(expScore, rows, cols, expGo, expGe, prob, maxP);
    if (flags != 0) {
        flags = k.runScaledDouble(score, rows, cols, expGo, expGe, prob, maxP);
        if (flags != 0) {
            k.runDouble(score, rows, cols, expGo, expGe, prob, maxP);
        }
        k.freeFallbackBuffers();
    }
}

// Scaled-mantissa kernel in place of forward, backward and computeProbabilityMatrix: fills P = zm and maxP
bool Impl::runScaled(bool profile) {
    if (scaledMode == FwBwAligner::SCALED_OFF) {
        return false;
    }
    P = zm;
    const int rows = (int)rowSeqLen, cols = (int)colSeqLen;
    if (profile) {
        // the profile's score row (already divided by T) and its exp, computed once per query, for each
        // target residue
        scaledRows.resize(2 * rowSeqLen);
        for (size_t i = 0; i < rowSeqLen; ++i) {
            scaledRows[i] = scoreForwardProfile[rowSeqAANum[i]];
            scaledRows[rowSeqLen + i] = scoreForwardProfile_exp[rowSeqAANum[i]];
        }
        runScaledKernel<true>(scaledProfile, scaledMode, scaledRows.data(), scaledRows.data() + rowSeqLen, rows, cols,
                              exp_go, exp_ge, zm, maxP);
    } else {
        runScaledKernel<false>(scaled, scaledMode, scoreForward, scoreForward, rows, cols, exp_go, exp_ge, zm, maxP);
    }
    return true;
}

void Impl::runFwBw(bool profile, int backtraceMode) {
    if (runScaled(profile)) {
        // the backtrace only reads P and maxP
        switch (backtraceMode) {
            case 0: break;
            case 2: computeBacktrace<2>(); break;
            case 3: computeBacktrace<3>(); break;
            default: computeBacktrace<1>(); break;
        }
        return;
    }
    if (profile) {
        switch (backtraceMode) {
            case 0: runFwBw<true, 0>(); break;
            case 1: runFwBw<true, 1>(); break;
            case 2: runFwBw<true, 2>(); break;
            case 3: runFwBw<true, 3>(); break;
            default: runFwBw<true, 1>(); break;
        }
    } else {
        switch (backtraceMode) {
            case 0: runFwBw<false, 0>(); break;
            case 1: runFwBw<false, 1>(); break;
            case 2: runFwBw<false, 2>(); break;
            case 3: runFwBw<false, 3>(); break;
            default: runFwBw<false, 1>(); break;
        }
    }
}

// The initial capacity follows the vector width: the block kernels process VECSIZE_FLOAT
// columns at a time.
size_t suggestedCapacity() {
    return VECSIZE_FLOAT * sizeof(float) * 20;
}

FwBwAligner::Impl *createImpl(SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float temperature,
                              float mact, size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace) {
    return new Impl(subMat, gapOpen, gapExtend, temperature, mact, rowsCapacity, colsCapacity, length, backtrace);
}

FwBwAligner::Impl *createImpl(float gapOpen, float gapExtend, float temperature, float mact,
                              size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace) {
    return new Impl(gapOpen, gapExtend, temperature, mact, rowsCapacity, colsCapacity, length, backtrace);
}
