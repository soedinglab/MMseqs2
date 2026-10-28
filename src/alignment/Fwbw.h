#ifndef FWBW_H
#define FWBW_H

#include "CpuFeatures.h"

#include <cstdint>
#include <cstdlib>
#include <string>

class SubstitutionMatrix;

// Forward-backward aligner with runtime SIMD dispatch.
//
// The kernel exists once, in FwbwInternal.cpp, and is compiled into one translation unit per
// instruction set (FwbwDefault.cpp at the binary's baseline ISA, FwbwAVX512.cpp with AVX512
// flags). Each of those includes simd.h with its own defines, which only works as long as no
// SIMD type or vector-width-dependent layout escapes into this header -- hence the pimpl.
class FwBwAligner {
public:
    typedef struct {
        float score1;
        float score2;
        int32_t dbStartPos1;
        int32_t dbEndPos1;
        int32_t	qStartPos1;
        int32_t qEndPos1;
        int32_t ref_end2;
        float qCov;
        float dbCov;
        std::string cigar;
        double evalue;
        int identicalAACnt;
        int32_t cigarLen;
        int word;
    } s_align;

    // Forward-backward kernel of runFwBw (all modes; the backtrace runs on its posteriors).
    enum ScaledMode {
        SCALED_OFF = 0,  // log-space blocks
        SCALED = 1,      // scaled mantissas; a call whose backward Z disagrees with the forward Z is run again
                         // with double mantissas, and in log-space double if that is flagged as well
        SCALED_FAST = 2  // scaled mantissas without that check, a few percent faster. WARNING: a path that
                         // starts next to a much stronger one can then be lost without any flag, and the
                         // posteriors around it are wrong by up to 1, not just slightly off. Seen on score
                         // matrices whose per-cell scores (divided by T) have a standard deviation of about
                         // 4-5 nats or more, earlier with cheap gap extension. Use SCALED unless the scores
                         // are known to be flat.
    };

    class Impl {
    public:
        virtual ~Impl() {}
        virtual void reallocateProfile(size_t newColsCapacity) = 0;
        virtual void resizeMatrix(bool profile, bool backtrace, size_t newRowsCapacity, size_t newColsCapacity) = 0;
        virtual void resetParams(float newGapOpen, float newGapExtend, float newTemperature) = 0;
        virtual void setTemperature(float newTemperature) = 0;
        virtual void setScaledMode(ScaledMode mode) = 0;
        virtual void initProfile(unsigned char* colAANum, size_t colAALen) = 0;
        virtual void initAlignment(unsigned char* targetNum, size_t targetLen, size_t queryLen) = 0;
        virtual void initScoreMatrix(float** inputScoreMatrix, int* gaps) = 0;
        virtual void runFwBw(bool profile, int backtraceMode) = 0;
        virtual s_align getFwbwAlnResult() const = 0;
        virtual float getProbability(int i, int j) const = 0;
        virtual float** getProbabiltyMatrix() = 0;
        virtual float getMaxP() const = 0;
        virtual float getTemperature() const = 0;
        virtual size_t getRowsCapacity() const = 0;
        virtual size_t getColsCapacity() const = 0;
        virtual size_t getBlockLength() const = 0;
    };

    FwBwAligner(SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float temperature, float mact,
                size_t rowsCapacity = 320, size_t colsCapacity = 320, size_t length = 16, int backtrace = 0);
    FwBwAligner(float gapOpen, float gapExtend, float temperature, float mact,
                size_t rowsCapacity = 320, size_t colsCapacity = 320, size_t length = 16, int backtrace = 0);
    ~FwBwAligner() { delete impl; }

    FwBwAligner(const FwBwAligner&) = delete;
    FwBwAligner& operator=(const FwBwAligner&) = delete;

    // Row/column capacity to start from, sized to the vector width of the selected backend.
    static size_t suggestedCapacity();

    // Name of the SIMD backend the dispatcher selected, for debug output.
    static const char *simdArchitecture();

    void reallocateProfile(size_t newColsCapacity) { impl->reallocateProfile(newColsCapacity); }
    void resizeMatrix(bool profile, bool backtrace, size_t newRowsCapacity, size_t newColsCapacity) {
        impl->resizeMatrix(profile, backtrace, newRowsCapacity, newColsCapacity);
    }
    void resetParams(float newGapOpen, float newGapExtend, float newTemperature) {
        impl->resetParams(newGapOpen, newGapExtend, newTemperature);
    }
    void initProfile(unsigned char* colAANum, size_t colAALen) { impl->initProfile(colAANum, colAALen); }
    void initAlignment(unsigned char* targetNum, size_t targetLen, size_t queryLen) {
        impl->initAlignment(targetNum, targetLen, queryLen);
    }
    void initScoreMatrix(float** inputScoreMatrix, int* gaps) { impl->initScoreMatrix(inputScoreMatrix, gaps); }
    void setTemperature(float newTemperature) { impl->setTemperature(newTemperature); }
    // forward-backward kernel of runFwBw, see ScaledMode
    void setScaledMode(ScaledMode mode) { impl->setScaledMode(mode); }

    // profile: score from the query profile (true) or from a user-supplied score matrix (false).
    // backtraceMode: 0 no backtrace, 1 local, 2 semi-global, 3 global.
    void runFwBw(bool profile, int backtraceMode) { impl->runFwBw(profile, backtraceMode); }

    template<bool profile, int backtrace>
    void runFwBw() { impl->runFwBw(profile, backtrace); }

    s_align getFwbwAlnResult() const { return impl->getFwbwAlnResult(); }
    float getProbability(int i, int j) const { return impl->getProbability(i, j); }
    float** getProbabiltyMatrix() { return impl->getProbabiltyMatrix(); }
    float getMaxP() const { return impl->getMaxP(); }
    float getTemperature() const { return impl->getTemperature(); }
    size_t getRowsCapacity() const { return impl->getRowsCapacity(); }
    size_t getColsCapacity() const { return impl->getColsCapacity(); }
    size_t getBlockLength() const { return impl->getBlockLength(); }

private:
    Impl *impl;
};

namespace FwBwDefault {
FwBwAligner::Impl *createImpl(SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float temperature,
                              float mact, size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace);
FwBwAligner::Impl *createImpl(float gapOpen, float gapExtend, float temperature, float mact,
                              size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace);
size_t suggestedCapacity();
const char *simdArchitecture();
}

#ifdef MMSEQS_HAVE_AVX512_DISPATCH
namespace FwBwAVX512 {
FwBwAligner::Impl *createImpl(SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float temperature,
                              float mact, size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace);
FwBwAligner::Impl *createImpl(float gapOpen, float gapExtend, float temperature, float mact,
                              size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace);
size_t suggestedCapacity();
const char *simdArchitecture();
}
#endif

inline FwBwAligner::FwBwAligner(SubstitutionMatrix &subMat, float gapOpen, float gapExtend, float temperature,
                                float mact, size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace)
    : impl(NULL) {
#ifdef MMSEQS_HAVE_AVX512_DISPATCH
    if (cpuHasAVX512()) {
        impl = FwBwAVX512::createImpl(subMat, gapOpen, gapExtend, temperature, mact, rowsCapacity, colsCapacity, length, backtrace);
        return;
    }
#endif
    impl = FwBwDefault::createImpl(subMat, gapOpen, gapExtend, temperature, mact, rowsCapacity, colsCapacity, length, backtrace);
}

inline FwBwAligner::FwBwAligner(float gapOpen, float gapExtend, float temperature, float mact,
                                size_t rowsCapacity, size_t colsCapacity, size_t length, int backtrace)
    : impl(NULL) {
#ifdef MMSEQS_HAVE_AVX512_DISPATCH
    if (cpuHasAVX512()) {
        impl = FwBwAVX512::createImpl(gapOpen, gapExtend, temperature, mact, rowsCapacity, colsCapacity, length, backtrace);
        return;
    }
#endif
    impl = FwBwDefault::createImpl(gapOpen, gapExtend, temperature, mact, rowsCapacity, colsCapacity, length, backtrace);
}

inline size_t FwBwAligner::suggestedCapacity() {
#ifdef MMSEQS_HAVE_AVX512_DISPATCH
    if (cpuHasAVX512()) {
        return FwBwAVX512::suggestedCapacity();
    }
#endif
    return FwBwDefault::suggestedCapacity();
}

inline const char *FwBwAligner::simdArchitecture() {
#ifdef MMSEQS_HAVE_AVX512_DISPATCH
    if (cpuHasAVX512()) { return FwBwAVX512::simdArchitecture(); }
#endif
    return FwBwDefault::simdArchitecture();
}

#endif //FWBW_H
