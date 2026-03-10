#include "BlockAligner.h"
#include "Sequence.h"
#include "Util.h"
#include "Parameters.h"
#include "SubstitutionMatrix.h"
#include "EvalueComputation.h"
#include "DistanceCalculator.h"

#define MAX_SIZE 4096 //change
#define MIN_SIZE 32

#include <algorithm>

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
    memset(queryCompBiasRev, 0, maxSequenceLength * sizeof(int16_t));
    //set targetCompBias to 0
    memset(targetCompBias, 0, maxSequenceLength * sizeof(int16_t));
    tmpCompBias   = new float[maxSequenceLength];
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
    delete[] queryCompBiasRev;
    delete[] queryCompBiasRevArr;

}

void BlockAligner::initQuery(Sequence* query){
    currentQuery = query;
    querySeq = query->getSeqData();
    queryNumSeq = query->numSequence;
    queryLength = query->L;
    //memset queryRevNumSeq
    memset(queryRevNumSeq, 0, maxSequenceLength * sizeof(int8_t));
    
    std::reverse_copy(queryNumSeq, queryNumSeq + queryLength, queryRevNumSeq);
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_AMINO_ACIDS)&& compBiasCorrection){
        SubstitutionMatrix::calcLocalAaBiasCorrection(subMat, queryNumSeq, queryLength, tmpCompBias, compBiasCorrectionScale);
        for (int i =0; i < queryLength; i++) { 
            queryCompBias[i] = (int16_t) (tmpCompBias[i]);
		}
        memset(this->queryCompBiasRev, 0, this->maxSequenceLength * sizeof(int8_t));
        std::reverse_copy(queryCompBias, queryCompBias + this->queryLength, queryCompBiasRev);
    }
}

// note: traceback cigar string will be reversed, but UngappedAln_res will contain correct start and end positions
s_align BlockAligner::gappedLocalAlign(
    Sequence* currentTarget,
    int qIdx, int tIdx,
    Cigar* cigar, int32_t x_drop,
    float covThr,
    int covMode
) {
    s_align local_aln;

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

    int qEnd = qIdx + res.query_idx -1 ;
    int tEnd = tIdx + res.reference_idx -1;
    
    float tmpqcov = SmithWaterman::computeCov(0, qEnd, currentQuery->L);
    float tmptcov = SmithWaterman::computeCov(0, tEnd, currentTarget->L);
    bool hasCov = Util::hasCoverage(covThr, covMode, tmpqcov, tmptcov);
    
    if (res.query_idx == SIZE_MAX || res.reference_idx == SIZE_MAX || !hasCov) {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // this should avoid that the hit is added
        return local_aln;
    }
    
    // Backward full alignment
    int bqueryAlnLen = qEnd +1;
    int bqueryStartPos = qLen - bqueryAlnLen;
    int btargetAlnLen = tEnd +1;
    int btargetStartPos = tLen - btargetAlnLen;

    int8_t* targetRevNumSeq = new int8_t[tLen];
	std::reverse_copy(tNum, tNum + tLen, targetRevNumSeq);

    block_set_bytes_padded_aa_numsequence(query, (uint8_t*)(queryRevNumSeq + bqueryStartPos), bqueryAlnLen, range.max);
    block_set_bytes_padded_aa_numsequence(target, (uint8_t*)(targetRevNumSeq + btargetStartPos), btargetAlnLen, range.max);
    
    //PosBias
    block_set_pos_bias(queryBias, queryCompBiasRev + bqueryStartPos, bqueryAlnLen);
    block_set_pos_bias(targetBias, targetCompBias + btargetStartPos, btargetAlnLen); // 0 

    block_align_aa_trace_xdrop_posbias(blockTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(blockTrace);
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
    delete[] targetRevNumSeq;
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
    s_align local_aln;

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

    block_align_aa_trace_xdrop_posbias(blockTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop); // forward with no trace
    res = block_res_aa_xdrop(blockTrace);
    block_cigar_aa_trace_xdrop(blockTrace, res.query_idx, res.reference_idx, cigar);

    int qEnd = qIdx + res.query_idx;
    int tEnd = tIdx + res.reference_idx;
    
    if (res.query_idx == SIZE_MAX || res.reference_idx == SIZE_MAX) {
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // this should avoid that the hit is added
        return local_aln;
    }

    
    local_aln.qEndPos1 = qEnd;
    local_aln.dbEndPos1 = tEnd;
    local_aln.score1 = res.score;

    return local_aln;
}

// note: traceback cigar string will be reversed, but UngappedAln_res will contain correct start and end positions
s_align BlockAligner::gappedLocalAlignBackward(
    Sequence* currentTarget,
    int qIdx, int tIdx,
    Cigar* cigar, int32_t x_drop
) {
    s_align local_aln;

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
    int8_t* targetRevNumSeq = new int8_t[tLen];
	std::reverse_copy(tNum, tNum + tLen, targetRevNumSeq);

    block_set_bytes_padded_aa_numsequence(query, (uint8_t*)(queryRevNumSeq + fqueryStartPos), fqueryAlnLen, range.max);
    block_set_bytes_padded_aa_numsequence(target, (uint8_t*)(targetRevNumSeq + ftargetStartPos), ftargetAlnLen, range.max);
    
    //PosBias
    block_set_pos_bias(queryBias, queryCompBiasRev + fqueryStartPos, fqueryAlnLen);
    block_set_pos_bias(targetBias, targetCompBias + ftargetStartPos, ftargetAlnLen);

    block_align_aa_trace_xdrop_posbias(blockTrace, query, queryBias, target, targetBias, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(blockTrace);
    block_cigar_aa_trace_xdrop(blockTrace, res.query_idx, res.reference_idx, cigar);
    if (res.query_idx == SIZE_MAX || res.reference_idx == SIZE_MAX) {
        // Debug(Debug::ERROR) << "wrong end position: " << qEnd << "\t" << tEnd << "\n";
        local_aln.score1 = 0.0f;
        local_aln.qCov = 0.0f;
        local_aln.tCov = 0.0f;
        local_aln.evalue = -1.0f; // this should avoid that the hit is added
        return local_aln;
    }

    local_aln.score1 = res.score;
    delete[] targetRevNumSeq;
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
    s_align local_aln;
    std::string backtrace_forward;
    s_align local_aln_Forward = bandedalignForward(currentTarget, qIdx, tIdx, backtrace_forward, xdrop);
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
    s_align local_aln_Backward;
    std::string backtrace_backward;
    if (qIdx > 0 && tIdx > 0) {
        local_aln_Backward = bandedalignBackward(currentTarget, qIdx, tIdx, backtrace_backward, xdrop);

    } else {
        local_aln_Backward.qStartPos1 = qIdx;
        local_aln_Backward.dbEndPos1 = tIdx;
        local_aln_Backward.score1 = 0;
        local_aln_Backward.identicalAACnt = 0;
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
    local_aln.score1 = local_aln_Forward.score1 + local_aln_Backward.score1; // temporary
    local_aln.identicalAACnt = local_aln_Forward.identicalAACnt + local_aln_Backward.identicalAACnt -1 ;
    
    float qcov = SmithWaterman::computeCov(local_aln.qStartPos1, local_aln.qEndPos1, currentQuery->L);
    float tcov = SmithWaterman::computeCov(local_aln.dbStartPos1, local_aln.dbEndPos1, currentTarget->L);
    double evalue = evaluer->computeEvalue(local_aln.score1, currentQuery->L);
    int bitScore = static_cast<int>(evaluer->computeBitScore(local_aln.score1) + 0.5);
    
    local_aln.qCov = qcov;
    local_aln.tCov = tcov;
    local_aln.evalue = evalue;
    local_aln.score1 = bitScore;

    local_aln.score2 = 0;
    local_aln.ref_end2 = -1;
    return local_aln;
}
