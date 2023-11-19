#include "BlockAligner.h"
#include "Sequence.h"
#include "Util.h"
#include "Parameters.h"
#include "SubstitutionMatrix.h"
#include "EvalueComputation.h"
#include "DistanceCalculator.h"

BlockAligner::BlockAligner(
    size_t maxSequenceLength,
    uintptr_t min,
    uintptr_t max,
    int8_t gapOpen,
    int8_t gapExtend,
    BaseMatrix& subMat,
    int dbtype
) : range({min, max}),
    gaps({gapOpen, gapExtend}),
    subMat(subMat),
    dbtype(dbtype),
    fastMatrix(SubstitutionMatrix::createAsciiSubMat(subMat)) {
    a = block_new_padded_aa(maxSequenceLength, max);
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        b = block_new_padded_aa(maxSequenceLength, max);
        matrix = block_new_simple_aamatrix(1, -1);
        for (int i = 0; i < subMat.alphabetSize; i++) {
            for (int j = 0; j < subMat.alphabetSize; j++) {
                block_set_aamatrix(
                    matrix,
                    subMat.num2aa[i],
                    subMat.num2aa[j],
                    subMat.subMatrix[i][j]
                );
            }
        }
    } else {
        bProfile = block_new_aaprofile(maxSequenceLength, max, gaps.extend);
    }
    blockTrace = block_new_aa_trace_xdrop(maxSequenceLength, maxSequenceLength, max);
    blockNoTrace = block_new_aa_xdrop(maxSequenceLength, maxSequenceLength, max);
    cigar = block_new_cigar(maxSequenceLength, maxSequenceLength);
    sAlnCigar = new uint32_t[maxSequenceLength];
}

BlockAligner::~BlockAligner() {
    block_free_cigar(cigar);
    block_free_aa_trace_xdrop(blockTrace);
    block_free_aa_xdrop(blockNoTrace);
    block_free_padded_aa(a);
    if (dbtype == Parameters::DBTYPE_AMINO_ACIDS) {
        block_free_padded_aa(b);
        block_free_aamatrix(matrix);
    } else {
        block_free_aaprofile(bProfile);
    }
    delete[] sAlnCigar;
}

void BlockAligner::initQuery(
    const char* querySeq,
    unsigned int queryLength,
    int querySeqType
) {
    this->querySeq = querySeq;
    this->queryLength = queryLength;
    this->querySeqType = querySeqType;

    queryConsensus = "";
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_HMM_PROFILE)) {
        Sequence::extractProfileConsensus(querySeq, queryLength * Sequence::PROFILE_READIN_SIZE, subMat, queryConsensus);
    }
}

// note: traceback cigar string will be reversed, but LocalAln will contain correct start and end positions
BlockAligner::LocalAln align_local(
    BlockHandle block_trace, BlockHandle block_no_trace,
    const char* a_str, size_t a_len, PaddedBytes* a,
    const char* b_str, size_t b_len, PaddedBytes* b,
    const AAMatrix* matrix, Gaps gaps,
    int a_idx, int b_idx,
    Cigar* cigar, SizeRange range, int32_t x_drop
) {
    BlockAligner::LocalAln res_aln;
    AlignResult res;

    // forwards alignment starting at (a_idx, b_idx)
    block_set_bytes_padded_aa(a, (uint8_t*)(a_str + a_idx), a_len - a_idx, range.max);
    block_set_bytes_padded_aa(b, (uint8_t*)(b_str + b_idx), b_len - b_idx, range.max);

    block_align_aa_xdrop(block_no_trace, a, b, matrix, gaps, range, x_drop);
    res = block_res_aa_xdrop(block_no_trace);

    res_aln.a_end = a_idx + res.query_idx;
    res_aln.b_end = b_idx + res.reference_idx;

    // reversed alignment starting at the max score location from forwards alignment
    block_set_bytes_rev_padded_aa(a, (uint8_t*)a_str, res_aln.a_end, range.max);
    block_set_bytes_rev_padded_aa(b, (uint8_t*)b_str, res_aln.b_end, range.max);

    block_align_aa_trace_xdrop(block_trace, a, b, matrix, gaps, range, x_drop);
    res = block_res_aa_trace_xdrop(block_trace);
    block_cigar_eq_aa_trace_xdrop(block_trace, a, b, res.query_idx, res.reference_idx, cigar);

    res_aln.a_start = res_aln.a_end - res.query_idx;
    res_aln.b_start = res_aln.b_end - res.reference_idx;
    res_aln.score = res.score;
    return res_aln;
}

BlockAligner::LocalAln align_local_profile(
    BlockHandle block_trace, BlockHandle block_no_trace,
    const char* a_str, const size_t a_len, PaddedBytes* a,
    const char* b_str, const size_t b_len, AAProfile* bProfile,
    Gaps gaps, BaseMatrix& subMat,
    int a_idx, int b_idx,
    Cigar* cigar, SizeRange range, int32_t x_drop
) {
    BlockAligner::LocalAln res_aln;
    AlignResult res;

    // forwards alignment starting at (a_idx, b_idx)
    block_set_bytes_padded_aa(a, (uint8_t*)(a_str + a_idx), a_len - a_idx, range.max);

    // assign extra profile columns to 'U', which is unused
    int aa = Sequence::PROFILE_READIN_SIZE;
    uint8_t order[Sequence::PROFILE_READIN_SIZE];
    memset(order, 'U', Sequence::PROFILE_READIN_SIZE);
    memcpy(order, (uint8_t*)subMat.num2aa, Sequence::PROFILE_AA_SIZE);

    block_clear_aaprofile(bProfile, b_len - b_idx, range.max);
    // note: scores are divided by 4 by shifting right by 2
    block_set_all_aaprofile(bProfile, order, aa, (int8_t*)(b_str + b_idx * aa), (b_len - b_idx) * aa, 0, 2);
    block_set_all_gap_open_C_aaprofile(bProfile, gaps.open);
    block_set_all_gap_close_C_aaprofile(bProfile, 0);
    block_set_all_gap_open_R_aaprofile(bProfile, gaps.open);

    block_align_profile_aa_xdrop(block_no_trace, a, bProfile, range, x_drop);
    res = block_res_aa_xdrop(block_no_trace);

    res_aln.a_end = a_idx + res.query_idx;
    res_aln.b_end = b_idx + res.reference_idx;

    // reversed alignment starting at the max score location from forwards alignment
    block_set_bytes_rev_padded_aa(a, (uint8_t*)a_str, res_aln.a_end, range.max);

    block_clear_aaprofile(bProfile, res_aln.b_end, range.max);
    block_set_all_rev_aaprofile(bProfile, order, aa, (int8_t*)b_str, res_aln.b_end * aa, 0, 2);
    block_set_all_gap_open_C_aaprofile(bProfile, gaps.open);
    block_set_all_gap_close_C_aaprofile(bProfile, 0);
    block_set_all_gap_open_R_aaprofile(bProfile, gaps.open);

    block_align_profile_aa_trace_xdrop(block_trace, a, bProfile, range, x_drop);
    res = block_res_aa_trace_xdrop(block_trace);
    block_cigar_aa_trace_xdrop(block_trace, res.query_idx, res.reference_idx, cigar);

    res_aln.a_start = res_aln.a_end - res.query_idx;
    res_aln.b_start = res_aln.b_end - res.reference_idx;
    res_aln.score = res.score;
    return res_aln;
}


BlockAligner::LocalAln BlockAligner::ungappedAlign(
    const char* targetData, size_t targetLen,
    int diagonal
) {
    const char* consSeq = queryConsensus.c_str();
    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeUngappedAlignment(
        queryConsensus.length() > 0 ? consSeq : querySeq, queryLength, targetData, targetLen,
        diagonal, fastMatrix.matrix, Parameters::RESCORE_MODE_ALIGNMENT
    );

    unsigned int distanceToDiagonal = alignment.distToDiagonal;
    diagonal = alignment.diagonal;

    int qUngappedStartPos, qUngappedEndPos, dbUngappedStartPos, dbUngappedEndPos;
    if (diagonal >= 0) {
        qUngappedStartPos = alignment.startPos + distanceToDiagonal;
        qUngappedEndPos = alignment.endPos + distanceToDiagonal;
        dbUngappedStartPos = alignment.startPos;
        dbUngappedEndPos = alignment.endPos;
    } else {
        qUngappedStartPos = alignment.startPos;
        qUngappedEndPos = alignment.endPos;
        dbUngappedStartPos = alignment.startPos + distanceToDiagonal;
        dbUngappedEndPos = alignment.endPos + distanceToDiagonal;
    }

    return LocalAln(
        qUngappedStartPos,
        dbUngappedStartPos,
        qUngappedEndPos,
        dbUngappedEndPos,
        alignment.score
    );
}

s_align
BlockAligner::align(
    const char* targetSeq,
    unsigned int targetLength,
    int diagonal,
    std::string& backtrace,
    EvalueComputation *evaluer,
    int xdrop
) {
    BlockAligner::LocalAln pos = ungappedAlign(
        targetSeq, targetLength,
        diagonal
    );

    LocalAln local_aln;
    if (queryConsensus.length() == 0) {
        local_aln = align_local(blockTrace, blockNoTrace, querySeq, queryLength, a, targetSeq, targetLength, b, matrix, gaps, pos.a_end, pos.b_end, cigar, range, xdrop);
    } else {
        // local_aln = align_local_profile(blockTrace, blockNoTrace, querySeq, queryLength, a, targetSeq, targetLength, bProfile, gaps, subMat, pos.a_end, pos.b_end, cigar, range, xdrop);
        local_aln = align_local_profile(blockTrace, blockNoTrace, targetSeq, targetLength, a, querySeq, queryLength, bProfile, gaps, subMat, pos.b_end, pos.a_end, cigar, range, xdrop);
        std::swap(local_aln.a_start, local_aln.b_start);
        std::swap(local_aln.a_end, local_aln.b_end);
    }

    float qcov = SmithWaterman::computeCov(local_aln.a_start, local_aln.a_end, queryLength);
    float dbcov = SmithWaterman::computeCov(local_aln.b_start, local_aln.b_end, targetLength);
    
    int bitScore = static_cast<int>(evaluer->computeBitScore(local_aln.score) + 0.5);
    double evalue = evaluer->computeEvalue(local_aln.score, queryLength);

    // Note: 'M' signals either a match or mismatch
    char ops_char[] = {' ', 'M', '=', 'X', 'I', 'D'};

    // int alnLength = Matcher::computeAlnLength(local_aln.a_start, local_aln.a_end, local_aln.b_start, local_aln.b_end);
    int alnLength = 0;

    size_t cigarLength = block_len_cigar(cigar);
    size_t aaIds = 0;
    if (cigarLength > 0) {
        int32_t targetPos = 0, queryPos = 0;
        for (unsigned long c = 0; c < cigarLength; ++c) {
            OpLen o = block_get_cigar(cigar, cigarLength - 1 - c);
            char letter = ops_char[o.op];
            uint32_t length = o.len;

            switch (letter) {
                case '=':
                    aaIds += length;
                    // FALLTHROUGH
                case 'X':
                    // FALLTHROUGH
                case 'M':
                    queryPos += length;
                    targetPos += length;
                    backtrace.append(length, 'M');
                    sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'M');
                    break;
                case 'I':
                    queryPos += length;
                    backtrace.append(length, 'I');
                    sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'I');
                    break;
                case 'D':
                    targetPos += length;
                    backtrace.append(length, 'D');
                    sAlnCigar[c] = SmithWaterman::to_cigar_int(length, 'D');
                    break;
            }
            alnLength += length;
        }
    }

    s_align r;
    r.score1 = bitScore;
    r.qStartPos1 = local_aln.a_start;
    r.qEndPos1 = local_aln.a_end;
    r.dbStartPos1 = local_aln.b_start;
    r.dbEndPos1 = local_aln.b_end;
    r.score2 = 0;
    r.ref_end2 = -1;
    r.qCov = qcov;
    r.tCov = dbcov;
    r.cigar = sAlnCigar;
    r.cigarLen = cigarLength;
    r.evalue = evalue;
    r.identicalAACnt = aaIds;
    return r;
}
