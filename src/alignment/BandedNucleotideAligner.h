//
// Written by Martin Steinegger
//
// Wrapper for KSW2 aligner.
// Local banded nucleotide aligner
//
#include <Parameters.h>
#include "StripedSmithWaterman.h"

#include "Util.h"
#include "SubstitutionMatrix.h"
#include "Debug.h"


class BandedNucleotideAligner {
public:


    BandedNucleotideAligner(BaseMatrix *subMat, size_t maxSequenceLength, int gapo, int gape);

    ~BandedNucleotideAligner();

    void initQuery(Sequence *q);

    s_align align(Sequence * targetSeqObj, short diagonal, EvalueComputation * evaluer);

private:
    SubstitutionMatrix::FastMatrix fastMatrix;
    uint8_t * targetSeq;
    uint8_t * targetSeqRev;
    uint8_t * querySeq;
    uint8_t * querySeqRev;
    Sequence * querySeqObj;
    int8_t * mat;
    uint32_t * cigar;
    int gapo;
    int gape;
};
