#ifndef MMSEQS_UNGAPPEDALIGNER_H
#define MMSEQS_UNGAPPEDALIGNER_H

#include "CpuFeatures.h"

#include <cstddef>
#include <cstdint>

class BaseMatrix;
class Sequence;

// Gapless (ungapped) diagonal scoring with runtime SIMD dispatch.
//
// The kernel exists once, in UngappedAlignerInternal.cpp, and is compiled into one translation
// unit per instruction set (UngappedAlignerDefault.cpp at the binary's baseline ISA,
// UngappedAlignerAVX512.cpp with AVX512 flags). Each of those includes simd.h with its own
// defines, which only works as long as no SIMD type escapes into this header -- hence the pimpl.
class UngappedAligner {
public:
    class Impl {
    public:
        virtual ~Impl() {}
        virtual void initQuery(const Sequence *q, const int8_t *mat, const BaseMatrix *m) = 0;
        virtual int score(const unsigned char *db_sequence, int32_t db_length) const = 0;
        virtual int score(const unsigned char *db_sequence, int32_t db_length, int &bestDiagonal) const = 0;
    };

    UngappedAligner(size_t maxSequenceLength, int aaSize, bool aaBiasCorrection, float aaBiasCorrectionScale);
    ~UngappedAligner() { delete impl; }

    UngappedAligner(const UngappedAligner&) = delete;
    UngappedAligner& operator=(const UngappedAligner&) = delete;

    void initQuery(const Sequence *q, const int8_t *mat, const BaseMatrix *m) {
        impl->initQuery(q, mat, m);
    }

    // Max diagonal score of query vs. db_sequence.
    int score(const unsigned char *db_sequence, int32_t db_length) const {
        return impl->score(db_sequence, db_length);
    }

    // Same score, and additionally the diagonal (queryPos - dbPos) the maximum was found on.
    // Callers that want to score a second channel along the alignment this pass found would
    // otherwise have to rescan the whole matrix. The score saturates at 255 - bias (unsigned
    // 8-bit accumulator), so on strongly matching pairs several diagonals can reach the ceiling
    // and the reported one is the first to get there rather than the true argmax. Use the
    // 2-argument form whenever the diagonal is not needed, so the hot path stays untouched.
    int score(const unsigned char *db_sequence, int32_t db_length, int &bestDiagonal) const {
        return impl->score(db_sequence, db_length, bestDiagonal);
    }

    // Name of the SIMD backend the dispatcher selected, for debug output.
    static const char *simdArchitecture();

private:
    Impl *impl;
};

namespace UngappedDefault {
UngappedAligner::Impl *createImpl(size_t maxSequenceLength, int aaSize,
                                  bool aaBiasCorrection, float aaBiasCorrectionScale);
const char *simdArchitecture();
}

#ifdef MMSEQS_HAVE_AVX512_DISPATCH
namespace UngappedAVX512 {
UngappedAligner::Impl *createImpl(size_t maxSequenceLength, int aaSize,
                                  bool aaBiasCorrection, float aaBiasCorrectionScale);
const char *simdArchitecture();
}
#endif

inline UngappedAligner::UngappedAligner(size_t maxSequenceLength, int aaSize,
                                        bool aaBiasCorrection, float aaBiasCorrectionScale)
    : impl(NULL) {
#ifdef MMSEQS_HAVE_AVX512_DISPATCH
    if (cpuHasAVX512()) {
        impl = UngappedAVX512::createImpl(maxSequenceLength, aaSize, aaBiasCorrection, aaBiasCorrectionScale);
        return;
    }
#endif
    impl = UngappedDefault::createImpl(maxSequenceLength, aaSize, aaBiasCorrection, aaBiasCorrectionScale);
}

inline const char *UngappedAligner::simdArchitecture() {
#ifdef MMSEQS_HAVE_AVX512_DISPATCH
    if (cpuHasAVX512()) {
        return UngappedAVX512::simdArchitecture();
    }
#endif
    return UngappedDefault::simdArchitecture();
}

#endif //MMSEQS_UNGAPPEDALIGNER_H
