#include "BlockAligner.h"
#include "Sequence.h"
#include "Util.h"
#include "Parameters.h"
#include "SubstitutionMatrix.h"
#include "EvalueComputation.h"
#include "DistanceCalculator.h"
#include "simd.h"

#define MAX_SIZE 4096 //change
#define MIN_SIZE 32

#include <algorithm>

// bestDiagonal: one lane per diagonal, 8-bit first, a group reruns in 16 bits if a lane may have saturated
static const int DIAG_LANES = VECSIZE_INT * 2;
static const int DIAG_CHAINS = 6;
static const int DIAG_GROUP = DIAG_CHAINS * DIAG_LANES;
static const int DIAG_MAX_STEP = 8;
static const int DIAG_MIN_SCORE = -2;
static const int DIAG_BYTE_SATURATED = 127 + DIAG_MIN_SCORE * DIAG_MAX_STEP;

BlockAligner::BlockAligner(
    int dbtype,
    size_t maxSequenceLength,
    BaseMatrix *m, SubstitutionMatrix::FastMatrix* fastMatrix, EvalueComputation * evaluer,
    bool compBiasCorrection, float compBiasCorrectionScale,
    int8_t gapOpen, int8_t gapExtend
) : 
    maxSequenceLength(maxSequenceLength),
    gaps({gapOpen, gapExtend}),
    compBiasCorrection(compBiasCorrection),
    compBiasCorrectionScale(compBiasCorrectionScale),
    dbtype(dbtype),
    subMat((SubstitutionMatrix*) m),
    fastMatrix(fastMatrix),
    evaluer(evaluer) {
    range={MIN_SIZE, MAX_SIZE};
    query = block_new_padded_aa(maxSequenceLength, range.max);
    queryBias = block_new_pos_bias(maxSequenceLength, range.max);
    queryRevNumSeq = new int8_t[maxSequenceLength];
    memset(queryRevNumSeq, 0, maxSequenceLength * sizeof(int8_t));
    targetRevNumSeq = new int8_t[maxSequenceLength];
    
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE) == false) {
        querySeqType = Parameters::DBTYPE_AMINO_ACIDS;
        target = block_new_padded_aa(maxSequenceLength, range.max);
        targetBias = block_new_pos_bias(maxSequenceLength, range.max);
        matrix = block_new_simple_aamatrix(1, -1);
        for (int aa1 = 0; aa1 < subMat->alphabetSize; aa1++) {
			for (int aa2 = 0; aa2 < subMat->alphabetSize; aa2++) {
				// instead of num2aa, use aa directly
				block_set_aamatrix_num(matrix, aa1, aa2,
								subMat->subMatrix[aa1][aa2]);
			}
		}
    } else {
        Debug(Debug::ERROR) << "blockaligner for profile is not supported yet. \n";
        EXIT(EXIT_FAILURE);
    }
    blockTrace = block_new_aa_trace_xdrop(maxSequenceLength, maxSequenceLength, range.max);
    blockNoTrace = block_new_aa_xdrop(maxSequenceLength, maxSequenceLength, range.max);
    cigar = block_new_cigar(maxSequenceLength, maxSequenceLength);
    queryCompBias = new int16_t[maxSequenceLength];
    targetCompBias = new int16_t[maxSequenceLength];
    queryCompBiasRev = new int16_t[maxSequenceLength];
    // read as positional bias even when compBiasCorrection never fills them, so zero them
    memset(queryCompBias, 0, maxSequenceLength * sizeof(int16_t));
    memset(queryCompBiasRev, 0, maxSequenceLength * sizeof(int16_t));
    memset(targetCompBias, 0, maxSequenceLength * sizeof(int16_t));
    tmpCompBias   = new float[maxSequenceLength];
    diagProfileStride = maxSequenceLength + 4 * DIAG_GROUP;
    diagProfileWord = new int16_t[subMat->alphabetSize * diagProfileStride];
    diagProfileByte = new int8_t[subMat->alphabetSize * diagProfileStride];
    diagSubMatrix = new int8_t[subMat->alphabetSize * subMat->alphabetSize];
    for (int a = 0; a < subMat->alphabetSize; a++) {
        for (int b = 0; b < subMat->alphabetSize; b++) {
            diagSubMatrix[a * subMat->alphabetSize + b] = static_cast<int8_t>(std::max(subMat->subMatrix[a][b], (short) DIAG_MIN_SCORE));
        }
    }
    diagProfileStale = true;
    diagProfileWordStale = true;
}

BlockAligner::~BlockAligner() {
    block_free_cigar(cigar);
    block_free_aa_trace_xdrop(blockTrace);
    block_free_aa_xdrop(blockNoTrace);
    block_free_padded_aa(query);
    block_free_pos_bias(queryBias);
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE) == false) {
        block_free_padded_aa(target);
        block_free_pos_bias(targetBias);
        block_free_aamatrix(matrix);
    }
    delete[] queryCompBias;
    delete[] targetCompBias;
    delete[] tmpCompBias;
    delete[] queryRevNumSeq;
    delete[] targetRevNumSeq;
    delete[] queryCompBiasRev;
    delete[] diagProfileWord;
    delete[] diagProfileByte;
    delete[] diagSubMatrix;
}

void BlockAligner::initQuery(Sequence* query){
    currentQuery = query;
    querySeq = query->getSeqData();
    queryNumSeq = query->numSequence;
    queryLength = query->L;
    diagProfileStale = true;
    // the alignments only read the first queryLength positions of these buffers
    std::reverse_copy(queryNumSeq, queryNumSeq + queryLength, queryRevNumSeq);
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_AMINO_ACIDS)&& compBiasCorrection){
        SubstitutionMatrix::calcLocalAaBiasCorrection(subMat, queryNumSeq, queryLength, tmpCompBias, compBiasCorrectionScale);
        for (int i =0; i < queryLength; i++) { 
            queryCompBias[i] = (int16_t) (tmpCompBias[i]);
		}
    } else {
        memset(queryCompBias, 0, queryLength * sizeof(int16_t));
    }
    std::reverse_copy(queryCompBias, queryCompBias + queryLength, queryCompBiasRev);
}

// the ungapped alignment only reads the query, so it can skip the reverse and composition bias of initQuery
void BlockAligner::setQuery(Sequence* query){
    currentQuery = query;
    querySeq = query->getSeqData();
    queryNumSeq = query->numSequence;
    queryLength = query->L;
    diagProfileStale = true;
}

// note: traceback cigar string will be reversed, but UngappedAln_res will contain correct start and end positions
s_align BlockAligner::gappedLocalAlign(
    Sequence* currentTarget,
    int qIdx, int tIdx,
    Cigar* cigar, int32_t x_drop,
    float covThr,
    int covMode
) {
    s_align local_aln = {};

    AlignResult res;
    res.score = -1000000000;
    const unsigned char* qNum = currentQuery->numSequence;
    const unsigned char* tNum = currentTarget->numSequence;
    size_t qLen = currentQuery->L;
    size_t tLen = currentTarget->L;

    // forwards alignment starting at (qIdx, tIdx)
    int fqueryAlnLen = qLen - qIdx;
    int fqueryStartPos = qIdx;
    int ftargetAlnLen = tLen - tIdx;
    int ftargetStartPos = tIdx;
    block_set_bytes_padded_aa_numsequence(query, (uint8_t*)(qNum + fqueryStartPos), fqueryAlnLen, range.max);
    block_set_bytes_padded_aa_numsequence(target, (uint8_t*)(tNum + ftargetStartPos), ftargetAlnLen, range.max);
    
    // PosBias
    block_set_pos_bias(queryBias, queryCompBias + fqueryStartPos, fqueryAlnLen);
    block_set_pos_bias(targetBias, targetCompBias + ftargetStartPos, ftargetAlnLen);

    block_align_aa_xdrop_posbias(blockNoTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop); // forward with no trace
    // block_align_aa_xdrop_posbias(blockNoTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop); // forward with no trace(blockNoTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop); 
    res = block_res_aa_xdrop(blockNoTrace);

    // Validate the forward result before deriving positions/coverage (SIZE_MAX or 0 == failure).
    if (res.query_idx == SIZE_MAX || res.reference_idx == SIZE_MAX
        || res.query_idx == 0 || res.reference_idx == 0) {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // this should avoid that the hit is added
        return local_aln;
    }

    int qEnd = qIdx + res.query_idx - 1;
    int tEnd = tIdx + res.reference_idx - 1;
    float tmpqcov = SmithWaterman::computeCov(0, qEnd, currentQuery->L);
    float tmptcov = SmithWaterman::computeCov(0, tEnd, currentTarget->L);
    if (!Util::hasCoverage(covThr, covMode, tmpqcov, tmptcov)) {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f;
        return local_aln;
    }

    // Backward full alignment
    int bqueryAlnLen = qEnd +1;
    int bqueryStartPos = qLen - bqueryAlnLen;
    int btargetAlnLen = tEnd +1;
    int btargetStartPos = tLen - btargetAlnLen;

	std::reverse_copy(tNum, tNum + tLen, targetRevNumSeq);

    block_set_bytes_padded_aa_numsequence(query, (uint8_t*)(queryRevNumSeq + bqueryStartPos), bqueryAlnLen, range.max);
    block_set_bytes_padded_aa_numsequence(target, (uint8_t*)(targetRevNumSeq + btargetStartPos), btargetAlnLen, range.max);
    
    //PosBias
    block_set_pos_bias(queryBias, queryCompBiasRev + bqueryStartPos, bqueryAlnLen);
    block_set_pos_bias(targetBias, targetCompBias + btargetStartPos, btargetAlnLen); // 0 

    block_align_aa_trace_xdrop_posbias(blockTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(blockTrace);
    // Validate the backward result before tracing back (avoid stale CIGAR).
    if (res.query_idx == SIZE_MAX || res.reference_idx == SIZE_MAX
        || res.query_idx == 0 || res.reference_idx == 0) {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f;
        return local_aln;
    }
    block_cigar_aa_trace_xdrop(blockTrace, res.query_idx, res.reference_idx, cigar);

    int score = res.score;
    double evalue = evaluer->computeEvalue(score, qLen);
    int bitScore = static_cast<int>(evaluer->computeBitScore(score) + 0.5);
    local_aln.qStartPos1 = qEnd - res.query_idx + 1;
    local_aln.dbStartPos1 = tEnd - res.reference_idx + 1;
    local_aln.score1 = bitScore;
    local_aln.qEndPos1 = qEnd;
    local_aln.dbEndPos1 = tEnd;
    local_aln.evalue = evalue;
    return local_aln;
}

//unused
BlockAligner::UngappedAln_res BlockAligner::hammingDistance(
    Sequence* currentTarget, const unsigned short diagonal)
{
    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeUngappedAlignment(
                                                        querySeq, queryLength, currentTarget->getSeqData(), currentTarget->L,
                                                        diagonal, fastMatrix->matrix, Parameters::RESCORE_MODE_HAMMING
                                                    );
    int diagonalLen = alignment.diagonalLen;
    int alnLen = diagonalLen;
    double evalue = 0.0;
    int bitScore = 0;
    float qCov = static_cast<float>(diagonalLen) / static_cast<float>(currentQuery->L);
    float tCov = static_cast<float>(diagonalLen) / static_cast<float>(currentTarget->L);

    return UngappedAln_res(
        bitScore,
        qCov,
        tCov,
        evalue, 
        alnLen,
        0,
        0,
        0,
        0,
        alignment.diagonalLen,
        alignment.score,
        alignment.diagonal
    );
}


template <typename T>
static void diagonalMaxScores(const T* queryProfile, const unsigned char* dbSequence, size_t profileStride, int j, int dbEnd, T* maxScores) {
    const int vectorSize = VECSIZE_INT * 4 / sizeof(T);
    const int vectorCount = DIAG_GROUP / vectorSize;
    simd_int vScore[vectorCount], vMaxScore[vectorCount];
    for (int c = 0; c < vectorCount; c++) {
        vScore[c] = vMaxScore[c] = (sizeof(T) == 1) ? simdi8_set(-128) : simdi16_set(-32768);
    }
    for (; j + DIAG_MAX_STEP <= dbEnd; j += DIAG_MAX_STEP) {
        for (int k = 0; k < DIAG_MAX_STEP; k++) {
            const T* profileRow = queryProfile + dbSequence[j + k] * profileStride + j + k;
            for (int c = 0; c < vectorCount; c++) {
                const simd_int vP = simdi_loadu((const simd_int*) (profileRow + c * vectorSize));
                vScore[c] = (sizeof(T) == 1) ? simdi8_adds(vScore[c], vP) : simdi16_adds(vScore[c], vP);
            }
        }
        for (int c = 0; c < vectorCount; c++) {
            vMaxScore[c] = (sizeof(T) == 1) ? simdi8_max(vMaxScore[c], vScore[c]) : simdi16_max(vMaxScore[c], vScore[c]);
        }
    }
    for (; j < dbEnd; j++) {
        const T* profileRow = queryProfile + dbSequence[j] * profileStride + j;
        for (int c = 0; c < vectorCount; c++) {
            const simd_int vP = simdi_loadu((const simd_int*) (profileRow + c * vectorSize));
            vScore[c] = (sizeof(T) == 1) ? simdi8_adds(vScore[c], vP) : simdi16_adds(vScore[c], vP);
            vMaxScore[c] = (sizeof(T) == 1) ? simdi8_max(vMaxScore[c], vScore[c]) : simdi16_max(vMaxScore[c], vScore[c]);
        }
    }
    for (int c = 0; c < vectorCount; c++) {
        simdi_store((simd_int*) (maxScores + c * vectorSize), vMaxScore[c]);
    }
}

// diagonalMaxScores<int8_t> of two adjacent groups in one pass
static void diagonalMaxScoresTwoGroups(const int8_t* queryProfile, const unsigned char* dbSequence, size_t profileStride, int j, int dbEndSecond, int dbEnd, int8_t* maxScores) {
    const int vectorSize = VECSIZE_INT * 4;
    const int vectorCount = DIAG_GROUP / vectorSize;
    simd_int vScore[2 * vectorCount], vMaxScore[2 * vectorCount];
    for (int c = 0; c < 2 * vectorCount; c++) {
        vScore[c] = vMaxScore[c] = simdi8_set(-128);
    }
    const int lastBlockEndSecond = j + std::max(0, dbEndSecond - j) / DIAG_MAX_STEP * DIAG_MAX_STEP;
    for (; j + DIAG_MAX_STEP <= lastBlockEndSecond; j += DIAG_MAX_STEP) {
        for (int k = 0; k < DIAG_MAX_STEP; k++) {
            const int8_t* profileRow = queryProfile + dbSequence[j + k] * profileStride + j + k;
            for (int c = 0; c < 2 * vectorCount; c++) {
                vScore[c] = simdi8_adds(vScore[c], simdi_loadu((const simd_int*) (profileRow + c * vectorSize)));
            }
        }
        for (int c = 0; c < 2 * vectorCount; c++) {
            vMaxScore[c] = simdi8_max(vMaxScore[c], vScore[c]);
        }
    }
    for (; j + DIAG_MAX_STEP <= dbEnd; j += DIAG_MAX_STEP) {
        for (int k = 0; k < DIAG_MAX_STEP; k++) {
            const int8_t* profileRow = queryProfile + dbSequence[j + k] * profileStride + j + k;
            for (int c = 0; c < 2 * vectorCount; c++) {
                vScore[c] = simdi8_adds(vScore[c], simdi_loadu((const simd_int*) (profileRow + c * vectorSize)));
            }
            for (int c = vectorCount; c < 2 * vectorCount; c++) {
                vMaxScore[c] = simdi8_max(vMaxScore[c], vScore[c]);
            }
        }
        for (int c = 0; c < vectorCount; c++) {
            vMaxScore[c] = simdi8_max(vMaxScore[c], vScore[c]);
        }
    }
    for (; j < dbEnd; j++) {
        const int8_t* profileRow = queryProfile + dbSequence[j] * profileStride + j;
        for (int c = 0; c < 2 * vectorCount; c++) {
            vScore[c] = simdi8_adds(vScore[c], simdi_loadu((const simd_int*) (profileRow + c * vectorSize)));
            vMaxScore[c] = simdi8_max(vMaxScore[c], vScore[c]);
        }
    }
    for (int c = 0; c < 2 * vectorCount; c++) {
        simdi_store((simd_int*) (maxScores + c * vectorSize), vMaxScore[c]);
    }
}

unsigned short BlockAligner::bestDiagonal(Sequence* currentTarget, int minDiag, int maxDiag) {
    if (diagProfileStale) {
        // locals, since the int8_t stores may alias the members
        const unsigned char* querySeq = queryNumSeq;
        const int queryLen = queryLength;
        diagProfileStride = queryLen + 4 * DIAG_GROUP;
        for (int a = 0; a < subMat->alphabetSize; a++) {
            const int8_t* subScores = diagSubMatrix + a * subMat->alphabetSize;
            int8_t* profileRow = diagProfileByte + a * diagProfileStride;
            std::fill(profileRow, profileRow + 2 * DIAG_GROUP, DIAG_MIN_SCORE);
            for (int i = 0; i < queryLen; i++) {
                profileRow[2 * DIAG_GROUP + i] = subScores[querySeq[i]];
            }
            std::fill(profileRow + 2 * DIAG_GROUP + queryLen, profileRow + diagProfileStride, DIAG_MIN_SCORE);
        }
        diagProfileStale = false;
        diagProfileWordStale = true;
    }
    const unsigned char* dbSequence = currentTarget->numSequence;
    const int dbLength = currentTarget->L;
    const int endDiag = std::min(maxDiag + 1, queryLength);
    int maxScore = -1;
    int maxDiagonal = 0;
    int8_t maxScoresByte[2 * DIAG_GROUP] __attribute__((aligned(64)));
    int16_t maxScoresWord[DIAG_GROUP] __attribute__((aligned(64)));
    for (int d0 = std::max(minDiag, -(dbLength - 1)); d0 < endDiag;) {
        // profileRow[j + lane] subScores query position d0 + j + lane against dbSequence[j]
        const int dbStart = std::max(0, -(d0 + DIAG_GROUP - 1));
        const int dbEnd = std::min(dbLength, queryLength - d0);
        const int dbStartSecond = std::max(0, -(d0 + 2 * DIAG_GROUP - 1));
        const bool twoGroups = (d0 + DIAG_GROUP < endDiag) && ((dbStart - dbStartSecond) % DIAG_MAX_STEP == 0);
        if (twoGroups) {
            const int dbEndSecond = std::min(dbLength, queryLength - d0 - DIAG_GROUP);
            diagonalMaxScoresTwoGroups(diagProfileByte + 2 * DIAG_GROUP + d0, dbSequence, diagProfileStride, dbStartSecond, dbEndSecond, dbEnd, maxScoresByte);
        } else {
            diagonalMaxScores<int8_t>(diagProfileByte + 2 * DIAG_GROUP + d0, dbSequence, diagProfileStride, dbStart, dbEnd, maxScoresByte);
        }
        for (int group = 0; group < (twoGroups ? 2 : 1); group++, d0 += DIAG_GROUP) {
            const int8_t* groupMaxScores = maxScoresByte + group * DIAG_GROUP;
            const int diagonalCount = std::min(DIAG_GROUP, endDiag - d0);
            bool saturated = false;
            for (int k = 0; k < diagonalCount; k++) {
                maxScoresWord[k] = groupMaxScores[k] + 128 - 32768;
                saturated |= (groupMaxScores[k] >= DIAG_BYTE_SATURATED);
            }
            if (saturated) {
                if (diagProfileWordStale) {
                    std::copy(diagProfileByte, diagProfileByte + subMat->alphabetSize * diagProfileStride, diagProfileWord);
                    diagProfileWordStale = false;
                }
                diagonalMaxScores<int16_t>(diagProfileWord + 2 * DIAG_GROUP + d0, dbSequence, diagProfileStride,
                                       std::max(0, -(d0 + DIAG_GROUP - 1)), std::min(dbLength, queryLength - d0), maxScoresWord);
            }
            for (int k = 0; k < diagonalCount; k++) {
                if (maxScoresWord[k] + 32768 > maxScore) {
                    maxScore = maxScoresWord[k] + 32768;
                    maxDiagonal = d0 + k;
                }
            }
        }
    }
    return static_cast<unsigned short>(maxDiagonal);
}

BlockAligner::UngappedAln_res BlockAligner::ungappedAlign(
    Sequence* target, const unsigned short diagonal)
{
    const char* targetSeq = target->getSeqData();
    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeUngappedAlignment(
                                                        querySeq, queryLength, targetSeq, target->L,
                                                        diagonal, fastMatrix->matrix, Parameters::RESCORE_MODE_ALIGNMENT
                                                    );
    unsigned int distanceToDiagonal = alignment.distToDiagonal;
    int ungappedDiagonal = alignment.diagonal;
    int distance = alignment.score;

    int qUngappedStartPos, qUngappedEndPos, tUngappedStartPos, tUngappedEndPos;
    if (ungappedDiagonal >= 0) {
        qUngappedStartPos = alignment.startPos + distanceToDiagonal;
        qUngappedEndPos = alignment.endPos + distanceToDiagonal;
        tUngappedStartPos = alignment.startPos;
        tUngappedEndPos = alignment.endPos;
    } else {
        qUngappedStartPos = alignment.startPos;
        qUngappedEndPos = alignment.endPos;
        tUngappedStartPos = alignment.startPos + distanceToDiagonal;
        tUngappedEndPos = alignment.endPos + distanceToDiagonal;
    }

    //check ungapped criteria
    double evalue = evaluer->computeEvalue(alignment.score, currentQuery->L);
    int bitScore = static_cast<int>(evaluer->computeBitScore(distance) + 0.5);
    int alnLen = (qUngappedEndPos - qUngappedStartPos) + 1;

    float qCov = SmithWaterman::computeCov(qUngappedStartPos, qUngappedEndPos , queryLength);
    float tCov = SmithWaterman::computeCov(tUngappedStartPos, tUngappedEndPos , target->L);
    
    return UngappedAln_res(
        bitScore,
        qCov,
        tCov,
        evalue, 
        alnLen,
        qUngappedStartPos,
        tUngappedStartPos,
        qUngappedEndPos,
        tUngappedEndPos,
        alignment.diagonalLen,
        alignment.score,
        alignment.diagonal
    );
}


s_align
BlockAligner::align(
    Sequence* currentTarget,
    size_t qIdx,
    size_t tIdx,
    std::string& backtrace,
    int xdrop,
    float covThr,
    int covMode
) {
    s_align local_aln = gappedLocalAlign(currentTarget, qIdx, tIdx, cigar, xdrop, covThr, covMode);
    // gappedLocalAlign returns before writing the CIGAR when it fails, so do not walk a stale one
    if (local_aln.evalue < 0) {
        return local_aln;
    }

    int aaIds = 0;
    size_t cigarLength = block_len_cigar(cigar);
    int queryPos = 0;
    int targetPos = 0;
    int queryStartPos = local_aln.qEndPos1;
    int targetStartPos = local_aln.dbEndPos1;
    const char* qNum = currentQuery->getSeqData();
    const char* tNum = currentTarget->getSeqData();
    // std::cout << "qstart: " << local_aln.qStartPos1 << " qend: " << local_aln.qEndPos1 << " cigarLength: " << cigarLength << std::endl;
    for (size_t c = 0; c < cigarLength; c++) {
        OpLen o = block_get_cigar(cigar, c);
        if(o.op == 1){
            for(size_t j = 0; j < o.len; j++){
                // change traceback with int not char
                if(qNum[-queryPos - j + queryStartPos] == tNum[-targetPos - j + targetStartPos]){
                    aaIds++;
                }
            }
            queryPos += o.len;
            targetPos += o.len;
            backtrace.append(o.len,'M');
        }else if(o.op == 4){
            queryPos += o.len;
            backtrace.append(o.len,'I');
        }else if(o.op == 5){
            targetPos += o.len;
            backtrace.append(o.len,'D');
        }
    }

    std::reverse(backtrace.begin(), backtrace.end());
    local_aln.qStartPos1 = (local_aln.qEndPos1+1) - queryPos;
	local_aln.dbStartPos1 = (local_aln.dbEndPos1+1) - targetPos;
    local_aln.qCov = SmithWaterman::computeCov(local_aln.qStartPos1, local_aln.qEndPos1, currentQuery->L);
	local_aln.tCov = SmithWaterman::computeCov(local_aln.dbStartPos1, local_aln.dbEndPos1, currentTarget->L);
    local_aln.score2 = 0;
    local_aln.ref_end2 = -1;
    local_aln.identicalAACnt = aaIds;
    return local_aln;
}



s_align BlockAligner::gappedLocalAlignForward(
    Sequence* currentTarget,
    int qIdx, int tIdx,
    Cigar* cigar, int32_t x_drop
) {
    s_align local_aln = {};

    AlignResult res;
    res.score = -1000000000;
    const unsigned char* qNum = currentQuery->numSequence;
    const unsigned char* tNum = currentTarget->numSequence;
    size_t qLen = currentQuery->L;
    size_t tLen = currentTarget->L;

    // forwards alignment starting at (qIdx, tIdx)
    int fqueryAlnLen = qLen - qIdx;
    int fqueryStartPos = qIdx;
    int ftargetAlnLen = tLen - tIdx;
    int ftargetStartPos = tIdx;
    block_set_bytes_padded_aa_numsequence(query, (uint8_t*)(qNum + fqueryStartPos), fqueryAlnLen, range.max);
    block_set_bytes_padded_aa_numsequence(target, (uint8_t*)(tNum + ftargetStartPos), ftargetAlnLen, range.max);
    
    // PosBias
    block_set_pos_bias(queryBias, queryCompBias + fqueryStartPos, fqueryAlnLen);
    block_set_pos_bias(targetBias, targetCompBias + ftargetStartPos, ftargetAlnLen);

    block_align_aa_trace_xdrop_posbias(blockTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(blockTrace);

    // a failed x-drop yields SIZE_MAX indices, so validate before tracing back from an invalid cell
    if (res.query_idx == SIZE_MAX || res.reference_idx == SIZE_MAX
        || res.query_idx == 0 || res.reference_idx == 0) {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // failure marker; caller must not read the (unwritten) CIGAR
        return local_aln;
    }

    block_cigar_aa_trace_xdrop(blockTrace, res.query_idx, res.reference_idx, cigar);

    // res.query_idx / res.reference_idx are counts, so the last aligned index is start + count - 1
    local_aln.qEndPos1 = qIdx + res.query_idx - 1;
    local_aln.dbEndPos1 = tIdx + res.reference_idx - 1;
    local_aln.score1 = res.score;

    return local_aln;
}

// note: traceback cigar string will be reversed, but UngappedAln_res will contain correct start and end positions
s_align BlockAligner::gappedLocalAlignBackward(
    Sequence* currentTarget,
    int qIdx, int tIdx,
    Cigar* cigar, int32_t x_drop
) {
    s_align local_aln = {};

    AlignResult res;
    res.score = -1000000000;
    size_t qLen = currentQuery->L;
    size_t tLen = currentTarget->L;
    
    const unsigned char* tNum = currentTarget->numSequence;
    // reversed alignment starting at the max score location from forwards alignment
    int fqueryAlnLen = qIdx +1;
    int fqueryStartPos = qLen - fqueryAlnLen;
    int ftargetAlnLen = tIdx +1;
    int ftargetStartPos = tLen - ftargetAlnLen;
	std::reverse_copy(tNum, tNum + tLen, targetRevNumSeq);

    block_set_bytes_padded_aa_numsequence(query, (uint8_t*)(queryRevNumSeq + fqueryStartPos), fqueryAlnLen, range.max);
    block_set_bytes_padded_aa_numsequence(target, (uint8_t*)(targetRevNumSeq + ftargetStartPos), ftargetAlnLen, range.max);
    
    //PosBias
    block_set_pos_bias(queryBias, queryCompBiasRev + fqueryStartPos, fqueryAlnLen);
    block_set_pos_bias(targetBias, targetCompBias + ftargetStartPos, ftargetAlnLen);

    block_align_aa_trace_xdrop_posbias(blockTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(blockTrace);
    // Validate BEFORE tracing back.
    if (res.query_idx == SIZE_MAX || res.reference_idx == SIZE_MAX
        || res.query_idx == 0 || res.reference_idx == 0) {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // failure marker; caller must not read the (unwritten) CIGAR
        return local_aln;
    }

    block_cigar_aa_trace_xdrop(blockTrace, res.query_idx, res.reference_idx, cigar);
    local_aln.score1 = res.score;
    return local_aln;
}


s_align
BlockAligner::bandedalignForward(
    Sequence* currentTarget,
    size_t qIdx,
    size_t tIdx,
    std::string& backtrace,
    int xdrop
) {
    //forward
    s_align local_aln = gappedLocalAlignForward(currentTarget, qIdx, tIdx, cigar, xdrop);
    // Forward failed: the CIGAR was not (re)written, so do not walk it (it would be stale).
    if (local_aln.evalue < 0) {
        return local_aln;
    }

    int aaIds = 0;
    size_t cigarLength = block_len_cigar(cigar);
    int queryPos = 0;
    int targetPos = 0;
    const char* qNum = currentQuery->getSeqData();
    const char* tNum = currentTarget->getSeqData();
    int queryStartPos = qIdx;
    int targetStartPos = tIdx;
    // Note: 'M' signals either a match or mismatch
    for (size_t c = 0; c < cigarLength; c++) {
        OpLen o = block_get_cigar(cigar, c);
        if(o.op == 1){
            for(size_t j = 0; j < o.len; j++){
                // change traceback with int not char
                if(qNum[queryStartPos + j + queryPos] == tNum[targetStartPos + j + targetPos]){
                    aaIds++;
                }
            }
            queryPos += o.len;
            targetPos += o.len;
            backtrace.append(o.len,'M');
        }else if(o.op == 4){
            queryPos += o.len;
            backtrace.append(o.len,'I');
        }else if(o.op == 5){
            targetPos += o.len;
            backtrace.append(o.len,'D');
        }
    }

    local_aln.qEndPos1 = qIdx + queryPos -1;
    local_aln.dbEndPos1 = tIdx + targetPos -1;
    if (cigarLength ==0){
        local_aln.qEndPos1 = qIdx;
        local_aln.dbEndPos1 = tIdx;
    }
    local_aln.score2 = 0;
    local_aln.ref_end2 = -1;
    local_aln.identicalAACnt = aaIds;
    return local_aln;
}

s_align
BlockAligner::bandedalignBackward(
    Sequence* currentTarget,
    size_t qIdx,
    size_t tIdx,
    std::string& backtrace,
    int xdrop
) {
    //Backward
    s_align local_aln = gappedLocalAlignBackward(currentTarget, qIdx, tIdx, cigar, xdrop);
    // Backward failed: the CIGAR was not (re)written, so do not walk it (it would be stale).
    if (local_aln.evalue < 0) {
        return local_aln;
    }

    int aaIds = 0;
    size_t cigarLength = block_len_cigar(cigar);
    int queryPos = 0;
    int targetPos = 0;
    int queryStartPos = qIdx;
    int targetStartPos = tIdx;
    const char* qNum = currentQuery->getSeqData();
    const char* tNum = currentTarget->getSeqData();

    for (size_t c = 0; c < cigarLength; c++) {
        OpLen o = block_get_cigar(cigar, c);
        if(o.op == 1){
            for(size_t j = 0; j < o.len; j++){
                // change traceback with int not char
                if(qNum[-queryPos - j + queryStartPos] == tNum[-targetPos - j + targetStartPos]){
                    aaIds++;
                }
            }
            queryPos += o.len;
            targetPos += o.len;
            backtrace.append(o.len,'M');
        }else if(o.op == 4){
            queryPos += o.len;
            backtrace.append(o.len,'I');
        }else if(o.op == 5){
            targetPos += o.len;
            backtrace.append(o.len,'D');
        }
    }

    std::reverse(backtrace.begin(), backtrace.end());
    
    local_aln.qStartPos1 = (qIdx + 1) - queryPos;
	local_aln.dbStartPos1 = (tIdx + 1) - targetPos;
    if (cigarLength ==0){
        local_aln.qStartPos1 = qIdx;
        local_aln.dbStartPos1 = tIdx;
    }
    local_aln.score2 = 0;
    local_aln.ref_end2 = -1;
    local_aln.identicalAACnt = aaIds;
    return local_aln;
}


s_align
BlockAligner::bandedalign(
    Sequence* currentTarget,
    size_t qIdx,
    size_t tIdx,
    std::string& backtrace,
    int xdrop,
    float covThr,
    int covMode
) {
    //Forward
    s_align local_aln = {};   // zero every field so failure returns never expose garbage (e.g. identicalAACnt)
    std::string backtrace_forward;
    s_align local_aln_Forward = bandedalignForward(currentTarget, qIdx, tIdx, backtrace_forward, xdrop);
    if (local_aln_Forward.evalue < 0) {   // forward extension failed
        local_aln.evalue = -1.0f;
        return local_aln;
    }
    float tmpqcov = SmithWaterman::computeCov(0, local_aln_Forward.qEndPos1, currentQuery->L);
    float tmptcov = SmithWaterman::computeCov(0, local_aln_Forward.dbEndPos1, currentTarget->L);
    bool hasCov = Util::hasCoverage(covThr, covMode, tmpqcov, tmptcov);
    if (!hasCov) {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // this should avoid that the hit is added
        return local_aln;
    }

    //Backward
    s_align local_aln_Backward = {};
    std::string backtrace_backward;
    const bool backwardStitched = (qIdx > 0 && tIdx > 0);
    if (qIdx > 0 && tIdx > 0) {
        local_aln_Backward = bandedalignBackward(currentTarget, qIdx, tIdx, backtrace_backward, xdrop);

    } else {
        local_aln_Backward.qStartPos1 = qIdx;
        local_aln_Backward.dbStartPos1 = tIdx;
        local_aln_Backward.score1 = 0;
        local_aln_Backward.identicalAACnt = 0;
    }
    if (backwardStitched && local_aln_Backward.evalue < 0) {   // backward extension failed
        local_aln.evalue = -1.0f;
        return local_aln;
    }

    if (backtrace_backward.empty() == false && backtrace_backward.back() != 'M') {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // this should avoid that the hit is added
        return local_aln;
    }
    if (backtrace_forward.empty() == false && backtrace_forward.front() != 'M') {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // this should avoid that the hit is added
        return local_aln;
    }
    
    //combine backtrace
    if (backtrace_backward.empty() == false){
        backtrace = backtrace_backward;
        backtrace.pop_back(); // remove last character to avoid double counting the position at (qIdx, tIdx)
    }
    backtrace.append(backtrace_forward);
    
    local_aln.qStartPos1 = local_aln_Backward.qStartPos1;
    local_aln.dbStartPos1 = local_aln_Backward.dbStartPos1;
    local_aln.qEndPos1 = local_aln_Forward.qEndPos1;
    local_aln.dbEndPos1 = local_aln_Forward.dbEndPos1;
    // Keep the raw (pre-bitscore) score signed while stitching; seedScore can be negative.
    int rawScore = static_cast<int>(local_aln_Forward.score1) + static_cast<int>(local_aln_Backward.score1);
    local_aln.identicalAACnt = local_aln_Forward.identicalAACnt + local_aln_Backward.identicalAACnt;
    if (backwardStitched) {
        // Seed cell is counted by both forward and backward extensions; subtract one copy.
        const int seedScore = subMat->subMatrix[currentQuery->numSequence[qIdx]][currentTarget->numSequence[tIdx]]
                            + (compBiasCorrection ? queryCompBias[qIdx] : 0);
        rawScore -= seedScore;
        local_aln.identicalAACnt -= 1;
    }
    
    float qcov = SmithWaterman::computeCov(local_aln.qStartPos1, local_aln.qEndPos1, currentQuery->L);
    float tcov = SmithWaterman::computeCov(local_aln.dbStartPos1, local_aln.dbEndPos1, currentTarget->L);
    double evalue = evaluer->computeEvalue(rawScore, currentQuery->L);
    int bitScore = static_cast<int>(evaluer->computeBitScore(rawScore) + 0.5);
    
    local_aln.qCov = qcov;
    local_aln.tCov = tcov;
    local_aln.evalue = evalue;
    local_aln.score1 = bitScore;

    local_aln.score2 = 0;
    local_aln.ref_end2 = -1;
    return local_aln;
}
