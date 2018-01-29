//
// Written by Martin Steinegger
//
// Wrapper for KSW2 aligner.
// Local banded nucleotide aligner
//
#include <Parameters.h>
#include <DistanceCalculator.h>
#include <ksw2/ksw2.h>
#include "BandedNucleotideAligner.h"

#include "Util.h"
#include "SubstitutionMatrix.h"
#include "Debug.h"
#include "StripedSmithWaterman.h"


BandedNucleotideAligner::BandedNucleotideAligner(BaseMatrix * subMat, size_t maxSequenceLength, int gapo, int gape) :
fastMatrix(SubstitutionMatrix::createAsciiSubMat(*subMat))
{

    targetSeq =  new uint8_t[maxSequenceLength];
    targetSeqRev =  new uint8_t[maxSequenceLength];
    querySeq =  new uint8_t[maxSequenceLength];
    querySeqRev =  new uint8_t[maxSequenceLength];
    mat = new int8_t[subMat->alphabetSize*subMat->alphabetSize];
    for (int i = 0; i < subMat->alphabetSize; i++) {
        for (int j = 0; j < subMat->alphabetSize; j++) {
            mat[i*subMat->alphabetSize + j] = subMat->subMatrix[i][j];
        }
    }
    this->gape = gape;
    this->gapo = gapo;
}

BandedNucleotideAligner::~BandedNucleotideAligner(){
    delete [] querySeq;
    delete [] targetSeq;
    delete [] targetSeqRev;
    delete [] querySeqRev;
    delete [] fastMatrix.matrixData;
    delete [] fastMatrix.matrix;
    delete [] mat;
}

void BandedNucleotideAligner::initQuery(Sequence * query){
    querySeqObj = query;
    for (int i = 0; i < query->L; ++i) {
        querySeq[i] = query->int_sequence[i];
    }
    SmithWaterman::seq_reverse((int8_t *)querySeqRev, (int8_t *)querySeq, query->L);
}



s_align BandedNucleotideAligner::align(Sequence * targetSeqObj, short diagonal,
                                       EvalueComputation * evaluer)
{
    for (int i = 0; i < targetSeqObj->L; ++i) {
        targetSeq[i] = targetSeqObj->int_sequence[i];
    }
    SmithWaterman::seq_reverse((int8_t *)targetSeqRev, (int8_t *)targetSeq, targetSeqObj->L);

    unsigned short distanceToDiagonal = abs(diagonal);
    DistanceCalculator::LocalAlignment alignment;
    int qUngappedStartPos, qUngappedEndPos, dbUngappedStartPos, dbUngappedEndPos;
    if (diagonal >= 0){
        int diagonalLen = std::min(targetSeqObj->L, querySeqObj->L - distanceToDiagonal);
        alignment = DistanceCalculator::computeSubstitutionStartEndDistance(
                querySeqObj->getSeqData() + distanceToDiagonal,
                targetSeqObj->getSeqData(),
                diagonalLen, fastMatrix.matrix);
        qUngappedStartPos = alignment.startPos + distanceToDiagonal;
        qUngappedEndPos = alignment.endPos + distanceToDiagonal;
        dbUngappedStartPos = alignment.startPos;
        dbUngappedEndPos = alignment.endPos;
    }else{
        int diagonalLen = std::min(targetSeqObj->L - distanceToDiagonal, querySeqObj->L);
        alignment = DistanceCalculator::computeSubstitutionStartEndDistance(querySeqObj->getSeqData(),
                                                                            targetSeqObj->getSeqData() +
                                                                            distanceToDiagonal,
                                                                            diagonalLen,
                                                                            fastMatrix.matrix);
        qUngappedStartPos = alignment.startPos;
        qUngappedEndPos = alignment.endPos;
        dbUngappedStartPos = alignment.startPos + distanceToDiagonal;
        dbUngappedEndPos = alignment.endPos + distanceToDiagonal;
    }

    if(qUngappedStartPos == 0 && qUngappedEndPos == querySeqObj->L -1
       && dbUngappedStartPos == 0 && dbUngappedEndPos == targetSeqObj->L - 1){
        s_align result;
        uint32_t * retCigar = new uint32_t[1];
        retCigar[0] = 0;
        retCigar[0] = querySeqObj->L << 4;
        result.cigar = retCigar;
        result.cigarLen = 1;
        result.score1 = alignment.score;
        result.qStartPos1 = qUngappedStartPos;
        result.qEndPos1 = qUngappedEndPos;
        result.dbEndPos1 = dbUngappedEndPos;
        result.dbStartPos1 = dbUngappedStartPos;
        result.qCov = SmithWaterman::computeCov(result.qStartPos1, result.qEndPos1, querySeqObj->L);
        result.tCov = SmithWaterman::computeCov(result.dbStartPos1, result.dbEndPos1, targetSeqObj->L);
        result.evalue = evaluer->computeEvalue(result.score1, querySeqObj->L);
        return result;
    }
//    printf("%d\t%d\t%d\n", alignment.score,  alignment.startPos, alignment.endPos);

    // get middle position of ungapped alignment
    int qStartRev = (querySeqObj->L  - qUngappedEndPos) - 1;
    int tStartRev = (targetSeqObj->L - dbUngappedEndPos) - 1;

    ksw_extz_t ez;
    int flag = 0;
    flag |= KSW_EZ_SCORE_ONLY;
    flag |= KSW_EZ_EXTZ_ONLY;
    ksw_extz2_sse(0, querySeqObj->L - qStartRev, querySeqRev + qStartRev, targetSeqObj->L - tStartRev, targetSeqRev + tStartRev, 5, mat, gapo, gape, 64, 40, flag, &ez);

    int qStartPos = querySeqObj->L  - ( qStartRev + ez.max_q ) -1 ;
    int tStartPos = targetSeqObj->L - ( tStartRev + ez.max_t ) -1;

    int alignFlag = 0;
    alignFlag |= KSW_EZ_EXTZ_ONLY;

    ksw_extz_t ezAlign;
//    ezAlign.cigar = cigar;
//    printf("%d %d\n", qStartPos, tStartPos);
    memset(&ezAlign, 0, sizeof(ksw_extz_t));
    ksw_extz2_sse(0, querySeqObj->L-qStartPos, querySeq+qStartPos, targetSeqObj->L-tStartPos, targetSeq+tStartPos, 5,
                  mat, gapo, gape, 64, 40, alignFlag, &ezAlign);

    std::string letterCode = "MID";
    uint32_t * retCigar = new uint32_t[ezAlign.n_cigar];
    for(int i = 0; i < ezAlign.n_cigar; i++){
        retCigar[i]=ezAlign.cigar[i];
    }
    s_align result;
    result.cigar = retCigar;
    result.cigarLen = ezAlign.n_cigar;
    result.score1 = ezAlign.max;
    result.qStartPos1 = qStartPos;
    result.qEndPos1 = ezAlign.max_q;
    result.dbEndPos1 = ezAlign.max_t;
    result.dbStartPos1 = tStartPos;
    result.qCov = SmithWaterman::computeCov(result.qStartPos1, result.qEndPos1, querySeqObj->L);
    result.tCov = SmithWaterman::computeCov(result.dbStartPos1, result.dbEndPos1, targetSeqObj->L);
    result.evalue = evaluer->computeEvalue(result.score1, querySeqObj->L);
    free(ezAlign.cigar);
    return result;
//        std::cout << static_cast<float>(aaIds)/ static_cast<float>(alignment.len) << std::endl;

}
