/* The MIT License
   Copyright (c) 2012-1015 Boston College.
   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
   Written by Michael Farrar, 2006 (alignment), Mengyao Zhao (SSW Library) and Martin Steinegger (change structure add aa composition, profile and AVX2 support).
   Please send bug reports and/or suggestions to martin.steinegger@snu.ac.kr.
*/
#include "Parameters.h"
#include "StripedSmithWaterman.h"
#include "UngappedAlignment.h"

#include "Util.h"
#include "SubstitutionMatrix.h"
#include "Debug.h"
#include "block_aligner.h"
#include <iostream>

#define MAX_SIZE 4096

struct s_block{
	PaddedBytes* query;
	PosBias* query_bias;
	AAMatrix* mat_aa;
	BlockHandle block_trace;
	int16_t* query_bias_arr;
};

SmithWaterman::SmithWaterman(size_t maxSequenceLength, int aaSize, bool aaBiasCorrection,
                             float aaBiasCorrectionScale, SubstitutionMatrix * subMat) {
	maxSequenceLength += 1;
	this->subMat = subMat;
    this->aaBiasCorrectionScale = aaBiasCorrectionScale;
	this->aaBiasCorrection = aaBiasCorrection;

	// int32_t alignment needs larger seqSize, was +7/8 for word before
    segSize = (maxSequenceLength+3)/4;
	vHStore = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	vHLoad  = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	vE      = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	vHmax   = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));

    isQueryProfile = false;

    // setting up query
	profile = new s_profile();
	// query profile
	profile->profile_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
	profile->profile_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_int = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));

	profile->profile_rev_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_rev_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
	profile->profile_rev_int = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));

	// query sequence
    profile->query_sequence     = new int8_t[maxSequenceLength];
	profile->query_rev_sequence = new int8_t[maxSequenceLength];
	profile->composition_bias   = new int8_t[maxSequenceLength];
	profile->composition_bias_rev   = new int8_t[maxSequenceLength];
	profile->profile_word_linear = new short*[aaSize];
	profile_word_linear_data = new short[aaSize*maxSequenceLength];
	profile->profile_int_linear = new int32_t*[aaSize];
	profile_int_linear_data = new int32_t[aaSize*maxSequenceLength];
	profile->mat_rev            = new int8_t[std::max(maxSequenceLength, (size_t)aaSize) * aaSize * 2]; // why multiply 2?
	profile->mat                = new int8_t[std::max(maxSequenceLength, (size_t)aaSize) * aaSize * 2];
	tmp_composition_bias   = new float[maxSequenceLength];
    scorePerCol = new int8_t[maxSequenceLength];
    /* array to record the largest score of each reference position */
	maxColumn = new uint8_t[maxSequenceLength*sizeof(uint32_t)];
	memset(maxColumn, 0, maxSequenceLength*sizeof(uint32_t));
	memset(profile->query_sequence, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->query_rev_sequence, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->mat_rev, 0, maxSequenceLength * aaSize);
	memset(profile->composition_bias, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->composition_bias_rev, 0, maxSequenceLength * sizeof(int8_t));

	// blockaligner
	block = new s_block();
	block->query = block_new_padded_aa(maxSequenceLength, MAX_SIZE);
	block->query_bias = block_new_pos_bias(maxSequenceLength, MAX_SIZE);
	block->mat_aa = block_new_simple_aamatrix(1, -1);
	block->block_trace = block_new_aa_trace_xdrop(maxSequenceLength, maxSequenceLength, MAX_SIZE);
	block->query_bias_arr = new int16_t[maxSequenceLength];

	profile->pos_aa_rev = new int8_t[maxSequenceLength * 32];
}

SmithWaterman::~SmithWaterman(){
	free(vHStore);
	free(vHLoad);
	free(vE);
	free(vHmax);
	free(profile->profile_byte);
	free(profile->profile_word);
	free(profile->profile_int);
	free(profile->profile_rev_byte);
	free(profile->profile_rev_word);
	free(profile->profile_rev_int);
	delete [] profile->query_rev_sequence;
	delete [] profile->query_sequence;
	delete [] profile->composition_bias;
	delete [] profile->composition_bias_rev;
	delete [] profile->profile_word_linear;
	delete [] profile_word_linear_data;
	delete [] profile->profile_int_linear;
	delete [] profile_int_linear_data;
	delete [] profile->mat_rev;
	delete [] profile->mat;
	delete [] tmp_composition_bias;
	delete [] scorePerCol;
	delete [] maxColumn;
	delete [] profile->pos_aa_rev;
	delete profile;

	block_free_padded_aa(block->query);
	block_free_pos_bias(block->query_bias);
	block_free_aamatrix(block->mat_aa);
	block_free_aa_trace_xdrop(block->block_trace);
	delete [] block->query_bias_arr;
	delete block;
}


/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch. */
template <typename T, size_t Elements, unsigned int type>
void SmithWaterman::createQueryProfile(simd_int *profile, const int8_t *query_sequence, const int8_t * composition_bias, const int8_t *mat,
        const int32_t query_length, const int32_t aaSize, uint8_t bias, const int32_t offset, const int32_t entryLength) {
	const int32_t segLen = (query_length + Elements - 1) / Elements;
	T* t = (T*) profile;
    for (int32_t nt = 0; LIKELY(nt < aaSize); nt++) {
		for (int32_t i = 0; i < segLen; i++) {
			int32_t  j = i;
			for (size_t segNum = 0; LIKELY(segNum < Elements) ; segNum++) {
				// if will be optmized out by compiler
				if(type == SUBSTITUTIONMATRIX) {    // substitution score for query_seq constrained by nt
					// query_sequence starts from 1 to n
					// *t++ = ( j >= query_length) ? bias : mat[nt * aaSize + query_sequence[j + offset ]] + composition_bias[j + offset] + bias; // mat[nt][q[j]] mat eq 20*20
					if (j >= query_length) {
						*t++ = bias;
					} else {
						const int q = query_sequence[j + offset];           
						const float cb = composition_bias[j + offset];      
					
						*t++ = mat[nt * aaSize + q] + cb + bias;        
					}
				} if(type == PROFILE) {
                    // profile starts by 0
//                    *t++ = (j >= query_length) ? bias : (mat[nt * entryLength + (j + (offset - 1))] + bias); //mat eq L*20  // mat[nt][j]
                    *t++ = (j >= query_length) ? bias : mat[nt * entryLength + j + offset] + bias;
//					// profile starts by 0 // TODO: offset?
//					*t++ = (j >= query_length) ? bias : mat[nt * entryLength + j + offset] + bias; //mat eq L*20  // mat[nt][j]
//					printf("(%1d, %1d) ", j , *(t-1));
				}
				j += segLen;
			}
			// std::cout << std::endl;
		}
	}
}

template <typename T, size_t Elements>
void SmithWaterman::updateQueryProfile(simd_int *profile, const int32_t query_length, const int32_t aaSize,
        uint8_t shift) {
    const int32_t segLen = (query_length + Elements - 1) / Elements;
    T* t = (T*) profile;
    for (uint32_t i = 0; i < segLen * Elements * aaSize; i++) {
        t[i] += shift;
    }
//    for (int32_t nt = 0; LIKELY(nt < aaSize); nt++) {
//        for (int32_t i = 0; i < segLen; i++) {
//            int32_t j = i;
//            for (size_t segNum = 0; LIKELY(segNum < Elements); segNum++) {
//                t* = t* + shift;
//                t++;
//                j += segLen;
//            }
//        }
//    }
}


uint8_t SmithWaterman::computeBias(const int32_t target_length, const int8_t *mat, const int32_t aaSize) {
    int8_t db_bias = 0;
    int32_t matSize = target_length * aaSize;
    for (int32_t i = 0; i < matSize; i++) {
        db_bias = std::min(mat[i], db_bias);
    }
    db_bias = abs(db_bias);
    return (db_bias);
}

// TODO: check the output (sanity check)
void SmithWaterman::reverseMat(int8_t *mat_rev, const int8_t *mat, const int32_t aaSize, const int32_t target_length) {
    for (int32_t i = 0; i < aaSize; i++) {
        const int8_t *startToRead = mat + (i * target_length);
        int8_t *startToWrite = mat_rev + (i * target_length);
        std::reverse_copy(startToRead, startToRead + target_length, startToWrite);
    }
}

s_align SmithWaterman::ssw_align (
        const unsigned char *db_num_sequence,
        int32_t db_length,
        std::string & backtrace,
        const uint8_t gap_open,
        const uint8_t gap_extend,
        const uint8_t alignmentMode,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
        const double  evalueThr,
        EvalueComputation * evaluer,
        const int covMode, const float covThr, const float correlationScoreWeight,
        const int32_t maskLen) {
    s_align alignment;
    // check if both query and target are profiles
	if (profile->isProfile) {
        alignment = ssw_align_private<SmithWaterman::PROFILE_SEQ>(db_num_sequence, db_length, backtrace, gap_open,
                                                                  gap_extend, alignmentMode, evalueThr, evaluer, covMode, covThr, correlationScoreWeight, maskLen);
    } else {
        alignment = ssw_align_private<SmithWaterman::SEQ_SEQ>(db_num_sequence, db_length, backtrace, gap_open,
                                                              gap_extend, alignmentMode, evalueThr, evaluer, covMode, covThr, correlationScoreWeight, maskLen);
    }
    return alignment;
}


template <unsigned int type>
s_align SmithWaterman::ssw_align_private (
	const unsigned char *db_sequence,
	int32_t db_length,
	std::string & backtrace,
	const uint8_t gap_open,
	const uint8_t gap_extend,
	const uint8_t alignmentMode,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
	const double  evalueThr,
	EvalueComputation * evaluer,
	const int covMode, const float covThr, const float correlationScoreWeight,
	const int32_t maskLen) {

	int32_t query_length = profile->query_length;

	// find the alignment position
	s_align align = alignScoreEndPos<type>(db_sequence, db_length, gap_open, gap_extend, maskLen);

	// no residue could be aligned
	if (align.dbEndPos1 == -1) {
		return align;
	}

	align.qCov = computeCov(0, align.qEndPos1, query_length);
	align.tCov = computeCov(0, align.dbEndPos1, db_length);

	bool hasLowerCoverage = !(Util::hasCoverage(covThr, covMode, align.qCov, align.tCov));
	align.evalue = evaluer->computeEvalue(align.score1, query_length);
	bool hasLowerEvalue = align.evalue > evalueThr;

	if (alignmentMode == 0 || ((alignmentMode == 2 || alignmentMode == 1) && (hasLowerEvalue || hasLowerCoverage))) {
		return align;
	}

	// run very shot and long overflowing alignments with SW instead of block aligner
	// short alignments are very fast with byte SW, long alignments produce slightly different scores FIXME
	if (align.word != 1) {
		return alignStartPosBacktrace<type>(db_sequence, db_length, gap_open, gap_extend, alignmentMode, backtrace, align, evaluer, covMode, covThr, correlationScoreWeight, maskLen);
	}

	bool blockAlignFailed = false;
	s_align alignTmp = alignStartPosBacktraceBlock<type>(db_sequence, db_length, gap_open, gap_extend, backtrace, align);
	if (align.score1 == UINT32_MAX) {
		blockAlignFailed = true;
	} else {
		align = alignTmp;
	}

	if (blockAlignFailed) {
		Debug(Debug::WARNING) << "Block alignment failed, falling back to Smith-Waterman\n";
		align = alignStartPosBacktrace<type>(db_sequence, db_length, gap_open, gap_extend, alignmentMode, backtrace, align, evaluer, covMode, covThr, correlationScoreWeight, maskLen);
	}

	// Check is needed (Below is for alignStartPosBacktraceBlock not for alignStartPosBacktrace since as it's already handled internally.)
	// align.qCov = computeCov(align.qStartPos1, align.qEndPos1, query_length);
    // align.tCov = computeCov(align.dbStartPos1, align.dbEndPos1, db_length);
    // hasLowerCoverage = !(Util::hasCoverage(covThr, covMode, align.qCov, align.tCov));

	return align;
}

template <unsigned int type>
s_align SmithWaterman::alignScoreEndPos (
		const unsigned char *db_sequence,
		int32_t db_length,
		const uint8_t gap_open,
		const uint8_t gap_extend,
		const int32_t maskLen) {
	int32_t query_length = profile->query_length;

	s_align r;
	r.dbStartPos1 = -1;
	r.qStartPos1 = -1;
	r.cigar = 0;
	r.cigarLen = 0;

	std::pair<alignment_end, alignment_end> bests;
	if (!profile->profile_byte) {
		Debug(Debug::ERROR) << "Not initialized profile\n";
	}
	// 1. byte
	bests = sw_sse2_byte<type>(db_sequence, 0, db_length, query_length, gap_open, gap_extend,
				profile->profile_byte, UCHAR_MAX, profile->bias, maskLen);
	r.word = 0;
	// 2. word
	if (bests.first.score == 255) {
		bests = sw_sse2_word<type>(db_sequence, 0, db_length, query_length, gap_open, gap_extend,
                    profile->profile_word, USHRT_MAX, maskLen);
        r.word = 1;
	}
	// 3. int
	// Comment out int32_t now for benchmark
	if (bests.first.score == INT16_MAX) {
		bests = sw_sse2_int<type>(db_sequence, 0, db_length, query_length, gap_open, gap_extend,
					profile->profile_int, USHRT_MAX, maskLen);
		r.word = 2;
	}

	r.score1 = bests.first.score;
    r.dbEndPos1 = bests.first.ref;
    r.qEndPos1 = bests.first.read;

	if (maskLen >= 15) {
		r.score2 = bests.second.score;
		r.ref_end2 = bests.second.ref;
	} else {
		r.score2 = 0;
		r.ref_end2 = -1;
	}
	return r;
}

template <unsigned int type>
s_align SmithWaterman::alignStartPosBacktraceBlock(
	const unsigned char *db_sequence,
	int32_t db_length,
	const uint8_t gap_open,
	const uint8_t gap_extend,
	std::string & backtrace,
	s_align r) {
	size_t query_len = profile->query_length;
	size_t target_len = db_length;
	Gaps gaps;
	gaps.open   = -gap_open;
	gaps.extend = -gap_extend;

	int32_t target_score = r.score1;
	AAProfile* queryProfile = nullptr;
	PosBias* target_bias = nullptr;

	// set query
	int32_t queryAlnLen = r.qEndPos1 + 1;
	int32_t queryStartPos = query_len - queryAlnLen;
	if (type == PROFILE_SEQ) {
		queryProfile = block_new_aaprofile(queryAlnLen, MAX_SIZE, gaps.extend);
		// Fill pos_aa_block and aa_pos_block with the relevant range(0-qEndPos, queryAlnLen)
		// extracted from the matrix initialized in ssw_init(0-qLen) for blockaligner
		// Replaced block_set_aaprofile, which set every position independently
		int8_t* pos_aa_block = aaprofile_pos_aa(queryProfile);
		int16_t* aa_pos_block = aaprofile_aa_pos(queryProfile);
		size_t curr_len_block = block_get_curr_len_aaprofile(queryProfile);

		for (int i = 0; i < queryAlnLen; i++) {
			// source: &profile->pos_aa_rev[(queryStartPos + i) * 32]
			// dest:   &pos_aa_block[(i + 1) * 32]
			memcpy(
				&pos_aa_block[(i + 1) * 32],
				&profile->pos_aa_rev[(queryStartPos + i) * 32],
				subMat->alphabetSize
			);
		}
		// transpose pos_aa_block to aa_pos_block
		for (int i = 0; i <= queryAlnLen; i++) {
			for (int b = 0; b < subMat->alphabetSize; b++) { // or 32 'A'-'Z'
				int8_t val = pos_aa_block[i * 32 + b];
				aa_pos_block[b * curr_len_block + i] = static_cast<int16_t>(val);
			}
		}
		// set all gap open and close values, including the costs of padding
		block_set_all_gap_open_C_aaprofile(queryProfile, gaps.open);
		block_set_all_gap_close_C_aaprofile(queryProfile, 0);
		block_set_all_gap_open_R_aaprofile(queryProfile, gaps.open);
	} else if (type == SEQ_SEQ) {
		// Since we use num in traceback, we don't need to convert num to aa.
		// block_set_bytes_padded_aa(block->query,  (const uint8_t*) (block->query_sequence_str.data() + queryStartPos), queryAlnLen, MAX_SIZE);
		block_set_bytes_padded_aa_numsequence(block->query, (const uint8_t*) (profile->query_rev_sequence + queryStartPos), queryAlnLen, MAX_SIZE);
		block_set_pos_bias(block->query_bias, block->query_bias_arr + queryStartPos, queryAlnLen);
	}
	// set target
	int32_t targetAlnLen = r.dbEndPos1 + 1;
	int32_t targetStartPos = target_len - targetAlnLen;
	PaddedBytes* target = block_new_padded_aa(target_len, MAX_SIZE);

	// Since we use num in traceback, we don't need to convert num to aa.
	// std::string db_sequence_str;
	// // copy this db_sequence,db_sequence + r.dbEndPos1 + 1 in reverse order to db_sequence_str and mappping to ascii using subMat->num2aa
	// for(int i = targetAlnLen - 1; i >= 0; i--){
	// 	db_sequence_str.push_back(subMat->num2aa[db_sequence[i]]);
	// }
	// block_set_bytes_padded_aa(target, (const uint8_t*) db_sequence_str.data(), targetAlnLen, MAX_SIZE);
	int8_t* db_rev_sequence = new int8_t[target_len];
	std::reverse_copy(db_sequence, db_sequence + db_length, db_rev_sequence);
	block_set_bytes_padded_aa_numsequence(target, (const uint8_t*) (db_rev_sequence+targetStartPos), targetAlnLen, MAX_SIZE);

	if (type == SEQ_SEQ){
		target_bias = block_new_pos_bias(targetAlnLen, MAX_SIZE); // fill 0
	}

	// substitute profile values to queryProfile
	AlignResult res;
	size_t min_size = 32;
	res.score = -1000000000;
	res.query_idx = -1;
	res.reference_idx = -1;

	if (type == SEQ_SEQ){
		while (min_size <= MAX_SIZE && res.score < target_score) {
			// allow max block size to grow
			SizeRange range;
			range.min = min_size;
			range.max = MAX_SIZE;
			// estimated x-drop threshold
			int32_t x_drop = -(min_size * gaps.extend + gaps.open);
			block_align_aa_trace_xdrop_posbias(block->block_trace, block->query, block->query_bias, target, target_bias,
										block->mat_aa, gaps, range, x_drop);
			res = block_res_aa_trace_xdrop(block->block_trace);
			min_size *= 2;
		}
	} else if (type == PROFILE_SEQ) {
		while (min_size <= MAX_SIZE && res.score < target_score) {
			// allow max block size to grow
			SizeRange range;
			range.min = min_size;
			range.max = MAX_SIZE;
			// estimated x-drop threshold
			int32_t x_drop = -(min_size * gaps.extend + gaps.open);
			block_align_profile_aa_trace_xdrop(block->block_trace, target, queryProfile, range, x_drop);
			res = block_res_aa_trace_xdrop(block->block_trace);
			min_size *= 2;
		}
	}

	size_t cigar_len, queryPos, targetPos;
	uint32_t aaIds;

	Cigar* cigar = block_new_cigar(res.query_idx, res.reference_idx);
	// char ops_char[] = {' ', 'M', '=', 'X', 'I', 'D'};
	if (res.score != target_score && !(target_score == INT16_MAX && res.score >= target_score)) {
		r.score1 = UINT32_MAX;
		goto cleanup;
	}

	block_cigar_aa_trace_xdrop(block->block_trace, res.query_idx, res.reference_idx, cigar);
	cigar_len = block_len_cigar(cigar);

	// Note: 'M' signals either query match or mismatch
	aaIds = 0;
	queryPos = 0;
	targetPos = 0;

	for (size_t i = 0; i < cigar_len; i++) {
		OpLen o = block_get_cigar(cigar, i);
		if(o.op == 1){
			for(size_t j = 0; j < o.len; j++){
				// change traceback with int not char
				if(profile->query_rev_sequence[queryPos + j + queryStartPos] == db_rev_sequence[targetPos + j + targetStartPos]){
					aaIds++;
				}
			}
			queryPos += o.len;
			targetPos += o.len;
			backtrace.append(o.len,'M');
		}else if(o.op == 4){
			switch (type) {
				case SEQ_SEQ:
					queryPos += o.len;
					backtrace.append(o.len,'I');
					break;
				case PROFILE_SEQ:
					targetPos += o.len;
					backtrace.append(o.len,'D');
					break;
			}
		}else if(o.op == 5){
			switch (type) {
				case SEQ_SEQ:
					targetPos += o.len;
					backtrace.append(o.len,'D');
					break;
				case PROFILE_SEQ:
					queryPos += o.len;
					backtrace.append(o.len,'I');
					break;
			}
		}
	}
	r.identicalAACnt = aaIds;

	//reverse backtrace
	std::reverse(backtrace.begin(), backtrace.end());
	r.qStartPos1 = (r.qEndPos1 + 1) - queryPos;
	r.dbStartPos1 = (r.dbEndPos1 + 1) - targetPos;

	r.qCov = computeCov(r.qStartPos1, r.qEndPos1, query_len);
	r.tCov = computeCov(r.dbStartPos1, r.dbEndPos1, db_length);

cleanup:
	block_free_padded_aa(target);
	block_free_cigar(cigar);
	if (type == PROFILE_SEQ) {
		block_free_aaprofile(queryProfile);
	} else if (type == SEQ_SEQ) {
		block_free_pos_bias(target_bias);
	}
	delete [] db_rev_sequence;
	return r;
}

template <unsigned int type>
s_align SmithWaterman::alignStartPosBacktrace (
        const unsigned char *db_sequence,
        int32_t db_length,
        const uint8_t gap_open,
        const uint8_t gap_extend,
        const uint8_t alignmentMode,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
        std::string & backtrace,
        s_align r,
		EvalueComputation * evaluer,
        const int covMode, const float covThr,
		const float correlationScoreWeight,
        const int32_t maskLen) {
    int32_t query_length = profile->query_length;
    int32_t queryOffset = query_length - r.qEndPos1 - 1;

	std::pair<alignment_end, alignment_end> bests_reverse;

    // Find the beginning position of the best alignment.
    if (r.word == 0) {
        if (type == PROFILE_SEQ) { // or type == PROFILE_SEQ
            createQueryProfile<int8_t, VECSIZE_INT * 4, PROFILE>(profile->profile_rev_byte, profile->query_rev_sequence, NULL, profile->mat_rev,
																 r.qEndPos1 + 1, profile->alphabetSize, profile->bias, queryOffset, profile->query_length);
		} else if (type == SEQ_SEQ) { // or type == SEQ_SEQ
			createQueryProfile<int8_t, VECSIZE_INT * 4, SUBSTITUTIONMATRIX>(profile->profile_rev_byte, profile->query_rev_sequence, profile->composition_bias_rev, profile->mat,
																 r.qEndPos1 + 1, profile->alphabetSize, profile->bias, queryOffset, 0);
		} else {
			fprintf(stderr, "Unknown type in alignStartPosBacktrace: %d\n", type);
			EXIT(EXIT_FAILURE);
		}
		bests_reverse = sw_sse2_byte<type>(db_sequence, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1, gap_open,
										   gap_extend, profile->profile_rev_byte,
										   r.score1, profile->bias, maskLen);
    } else if (r.word == 1) {
        if (type == PROFILE_SEQ) {
            createQueryProfile<int16_t, VECSIZE_INT * 2, PROFILE>(profile->profile_rev_word, profile->query_rev_sequence, NULL, profile->mat_rev,
																  r.qEndPos1 + 1, profile->alphabetSize, 0, queryOffset, profile->query_length);
		} else if (type == SEQ_SEQ) {
			createQueryProfile<int16_t, VECSIZE_INT * 2, SUBSTITUTIONMATRIX>(profile->profile_rev_word, profile->query_rev_sequence, profile->composition_bias_rev, profile->mat,
																  r.qEndPos1 + 1, profile->alphabetSize, 0, queryOffset, 0);
		} else {
			fprintf(stderr, "Unknown type in alignStartPosBacktrace: %d\n", type);
			EXIT(EXIT_FAILURE);
		}
		bests_reverse = sw_sse2_word<type>(db_sequence, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1, gap_open,
										   gap_extend, profile->profile_rev_word,
										   r.score1, maskLen);
	}
	// Comment out int32_t now for benchmark
	else if (r.word == 2) {
        if ((type == PROFILE_SEQ)) {
			createQueryProfile<int32_t, VECSIZE_INT * 1, PROFILE>(profile->profile_rev_int, profile->query_rev_sequence, NULL, profile->mat_rev,
																	r.qEndPos1 + 1, profile->alphabetSize, 0, queryOffset, profile->query_length);
		}  else if (type == SEQ_SEQ) {
			createQueryProfile<int32_t, VECSIZE_INT * 1, SUBSTITUTIONMATRIX>(profile->profile_rev_int, profile->query_rev_sequence, profile->composition_bias_rev, profile->mat,
																	r.qEndPos1 + 1, profile->alphabetSize, 0, queryOffset, 0);
		}
		bests_reverse = sw_sse2_int<type>(db_sequence, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1, gap_open,
											gap_extend, profile->profile_rev_int,
											r.score1, maskLen);
	}

    if(bests_reverse.first.score != r.score1){
		Debug(Debug::ERROR) << "r.word: " << r.word << "\n";
        Debug(Debug::ERROR) << "bests_reverse.first.score: " << bests_reverse.first.score << "\n";
		Debug(Debug::ERROR) << "r.score1: " << r.score1 << "\n";
		Debug(Debug::ERROR) << "Score of forward/backward SW differ. This should not happen.\n";
        Debug(Debug::ERROR) << "Start: Q: " << (r.qEndPos1 - bests_reverse.first.read) << ", T: " << bests_reverse.first.ref << ". End: Q: " << r.qEndPos1 << ", T " << r.dbEndPos1 << "\n";
		//  if qry is not a profile, just exit
        if (!(type == PROFILE_SEQ)) {
            EXIT(EXIT_FAILURE);
        }
    }

    r.dbStartPos1 = bests_reverse.first.ref;
    r.qStartPos1 = r.qEndPos1 - bests_reverse.first.read;

    if (r.dbStartPos1 == -1) {
        fprintf(stderr, "Target start position is -1. This should not happen.\n");
        EXIT(EXIT_FAILURE);
    }

    r.qCov = computeCov(r.qStartPos1, r.qEndPos1, query_length);
    r.tCov = computeCov(r.dbStartPos1, r.dbEndPos1, db_length);
    bool hasLowerCoverage = !(Util::hasCoverage(covThr, covMode, r.qCov, r.tCov));

    // only start and end point are needed
    if (alignmentMode == 1 || hasLowerCoverage) {
        return r;
    }

    // Generate cigar.
    db_length = r.dbEndPos1 - r.dbStartPos1 + 1;
    query_length = r.qEndPos1 - r.qStartPos1 + 1;
    int32_t band_width = abs(db_length - query_length) + 1;

    cigar* path;
    if (type == PROFILE_SEQ) {
        path = banded_sw<type>(db_sequence + r.dbStartPos1, profile->query_sequence + r.qStartPos1, profile->composition_bias + r.qStartPos1, db_length,
							   query_length, r.qStartPos1, r.score1, gap_open, gap_extend,
							   band_width, profile->mat, profile->query_length);
	} else {
        path = banded_sw<type>(db_sequence + r.dbStartPos1, profile->query_sequence + r.qStartPos1, profile->composition_bias + r.qStartPos1, db_length,
							   query_length, r.qStartPos1, r.score1, gap_open, gap_extend,
							   band_width, profile->mat, profile->alphabetSize);
		db_length = r.dbEndPos1 - r.dbStartPos1 + 1;
		query_length = r.qEndPos1 - r.qStartPos1 + 1;
		band_width = abs(db_length - query_length) + 1;
	}

    if (path != NULL) {
        r.cigar = path->seq;
        r.cigarLen = path->length;
    }

    uint32_t aaIds = 0;
    size_t mStateCnt = 0;
	// Need check below for Profile_seq
	computerBacktrace(profile, db_sequence, r, backtrace, aaIds, scorePerCol, mStateCnt);
    r.identicalAACnt = aaIds;
	if(correlationScoreWeight > 0.0){
        int correlationScore = computeCorrelationScore(scorePerCol, mStateCnt);
        r.score1 += static_cast<float>(correlationScore) * correlationScoreWeight;
        r.evalue = evaluer->computeEvalue(r.score1, query_length);
    }
	if(path != NULL) {
        delete path;
    }
    return r;
}

template
s_align SmithWaterman::ssw_align_private<SmithWaterman::SEQ_SEQ>(const unsigned char*, int32_t, std::string&, const uint8_t, const uint8_t, const uint8_t, const double, EvalueComputation*, const int, const float, const float, const int32_t);
template
s_align SmithWaterman::ssw_align_private<SmithWaterman::PROFILE_SEQ>(const unsigned char*, int32_t, std::string&, const uint8_t, const uint8_t, const uint8_t, const double, EvalueComputation*, const int, const float, const float, const int32_t);

template
s_align SmithWaterman::alignScoreEndPos<SmithWaterman::SEQ_SEQ>(const unsigned char*, int32_t, const uint8_t, const uint8_t, const int32_t);
template
s_align SmithWaterman::alignScoreEndPos<SmithWaterman::PROFILE_SEQ>(const unsigned char*, int32_t, const uint8_t, const uint8_t, const int32_t);

template
s_align SmithWaterman::alignStartPosBacktraceBlock<SmithWaterman::SEQ_SEQ>(const unsigned char*, int32_t, const uint8_t, const uint8_t, std::string&, s_align);
template
s_align SmithWaterman::alignStartPosBacktraceBlock<SmithWaterman::PROFILE_SEQ>(const unsigned char*, int32_t, const uint8_t, const uint8_t, std::string&, s_align);

template
s_align SmithWaterman::alignStartPosBacktrace<SmithWaterman::SEQ_SEQ>(const unsigned char*, int32_t, const uint8_t, const uint8_t, const uint8_t, std::string&, s_align, EvalueComputation*, const int, const float, const float, const int32_t);
template
s_align SmithWaterman::alignStartPosBacktrace<SmithWaterman::PROFILE_SEQ>(const unsigned char*, int32_t, const uint8_t, const uint8_t, const uint8_t, std::string&, s_align, EvalueComputation*, const int, const float, const float, const int32_t);

void SmithWaterman::computerBacktrace(s_profile * query, const unsigned char * db_sequence,
                                      s_align & alignment, std::string & backtrace,
                                      uint32_t & aaIds, int8_t * scorePerCol, size_t & mStatesCnt){
    int32_t targetPos = alignment.dbStartPos1, queryPos = alignment.qStartPos1;
    for (int32_t c = 0; c < alignment.cigarLen; ++c) {
        char letter = SmithWaterman::cigar_int_to_op(alignment.cigar[c]);
        uint32_t length = SmithWaterman::cigar_int_to_len(alignment.cigar[c]);
        backtrace.reserve(length);
        for (uint32_t i = 0; i < length; ++i){
            if (letter == 'M') {
                aaIds += (db_sequence[targetPos] == query->query_sequence[queryPos]);
				scorePerCol[mStatesCnt] = query->mat[query->query_sequence[queryPos] * query->alphabetSize + db_sequence[targetPos]] + query->composition_bias[queryPos];
                ++mStatesCnt;
                ++queryPos;
                ++targetPos;
                backtrace.append("M");
            } else {
                if (letter == 'I') {
                    ++queryPos;
                    backtrace.append("I");
                }
                else{
                    ++targetPos;
                    backtrace.append("D");
                }
            }
        }
    }
}

char SmithWaterman::cigar_int_to_op(uint32_t cigar_int) {
	uint8_t letter_code = cigar_int & 0xfU;
	static const char map[] = {
			'M',
			'I',
			'D',
			'N',
			'S',
			'H',
			'P',
			'=',
			'X',
	};

	if (letter_code >= (sizeof(map)/sizeof(map[0]))) {
		return 'M';
	}

	return map[letter_code];
}

uint32_t SmithWaterman::cigar_int_to_len (uint32_t cigar_int)
{
	uint32_t res = cigar_int >> 4;
	return res;
}


int SmithWaterman::computeCorrelationScore(int8_t * scorePreCol, size_t length){
    int corrScore1=0;
    int corrScore2=0;
    int corrScore3=0;
    int corrScore4=0;
    if(length >= 2){
        corrScore1 += scorePreCol[1] * scorePreCol[0];
    }
    if(length >= 3){
        corrScore1 += scorePreCol[2] * scorePreCol[1];
        corrScore2 += scorePreCol[2] * scorePreCol[0];
    }
    if(length >= 4){
        corrScore1 += scorePreCol[3] * scorePreCol[2];
        corrScore2 += scorePreCol[3] * scorePreCol[1];
        corrScore3 += scorePreCol[3] * scorePreCol[0];
    }
    for (size_t step = 4; step < length; step++){
        corrScore1 += scorePreCol[step] * scorePreCol[step - 1];
        corrScore2 += scorePreCol[step] * scorePreCol[step - 2];
        corrScore3 += scorePreCol[step] * scorePreCol[step - 3];
        corrScore4 += scorePreCol[step] * scorePreCol[step - 4];
    }
    return (corrScore1+corrScore2+corrScore3+corrScore4);
}


template <unsigned int type>
std::pair<SmithWaterman::alignment_end, SmithWaterman::alignment_end> SmithWaterman::sw_sse2_byte (
														   const unsigned char *db_sequence,
														   int8_t ref_dir,	// 0: forward ref; 1: reverse ref
														   int32_t db_length,
														   int32_t query_length,
														   const uint8_t gap_open, /* will be used as - */
														   const uint8_t gap_extend, /* will be used as - */
														   const simd_int* query_profile_byte, /* profile_byte loaded in ssw_init */
														   uint8_t terminate,	/* the best alignment score: used to terminate
                                                         the matrix calculation when locating the
                                                         alignment beginning point. If this score
                                                         is set to 0, it will not be used */
														   uint8_t bias,  /* Shift 0 point to a positive value. */
														   int32_t maskLen) {
	uint8_t max = 0;		                     /* the max alignment score */
	int32_t end_query = query_length - 1;
	int32_t end_db = -1; /* 0_based best alignment ending point; Initialized as isn't aligned -1. */
	const int SIMD_SIZE = VECSIZE_INT * 4;
	int32_t segLen = (query_length + SIMD_SIZE-1) / SIMD_SIZE; /* number of segment */
	/* array to record the largest score of each reference position */
	memset(this->maxColumn, 0, db_length * sizeof(uint8_t));
	uint8_t * maxColumn = (uint8_t *) this->maxColumn;

	/* Define 16 byte 0 vector. */
	simd_int vZero = simdi32_set(0);
	simd_int* pvHStore = vHStore;
	simd_int* pvHLoad = vHLoad;
	simd_int* pvE = vE;
	simd_int* pvHmax = vHmax;

    memset(pvHStore,0,segLen*sizeof(simd_int));
	memset(pvHLoad,0,segLen*sizeof(simd_int));
	memset(pvE,0,segLen*sizeof(simd_int));
	memset(pvHmax,0,segLen*sizeof(simd_int));

	int32_t i, j;
    /* 16 byte insertion begin vector */
	simd_int vGapO = simdi8_set(gap_open);

    /* 16 byte insertion extension vector */
	simd_int vGapE = simdi8_set(gap_extend);

	/* 16 byte bias vector */
	simd_int vBias = simdi8_set(bias);

	simd_int vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
	simd_int vMaxMark = vZero; /* Trace the highest score till the previous column. */
	simd_int vTemp;
	int32_t edge, begin = 0, end = db_length, step = 1;

    //fprintf(stderr, "start alignment of length %d [%u]\n", query_length, segLen * SIMD_SIZE);

	/* outer loop to process the reference sequence */
	if (ref_dir == 1) {
		begin = db_length - 1;
		end = -1;
		step = -1;
	}


	for (i = begin; LIKELY(i != end); i += step) {
//	    cnt = i;
		simd_int e, vF = vZero, vMaxColumn = vZero; /* Initialize F value to 0.
                                                    Any errors to vH values will be corrected in the Lazy_F loop.
                                                    */

		simd_int vH = pvHStore[segLen - 1];
		vH = simdi8_shiftl (vH, 1); /* Shift the 128-bit value in vH left by 1 byte. */
		const simd_int* vP = query_profile_byte + db_sequence[i] * segLen; /* Right part of the query_profile_byte */

		/* Swap the 2 H buffers. */
		simd_int* pv = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pv;

		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); ++j) {
		    simd_int score = simdi_load(vP + j);

            vH = simdui8_adds(vH, score);
            vH = simdui8_subs(vH, vBias);   /* vH will be always > 0 */

            /* Get max from vH, vE and vF. */
			e = simdi_load(pvE + j);
            vH = simdui8_max(vH, e);
            vH = simdui8_max(vH, vF);
			vMaxColumn = simdui8_max(vMaxColumn, vH);

			/* Save vH values. */
			simdi_store(pvHStore + j, vH);

			/* Update vE value. */
			vH = simdui8_subs(vH, vGapO); /* saturation arithmetic, result >= 0 */
			e = simdui8_subs(e, vGapE);
			e = simdui8_max(e, vH);
			simdi_store(pvE + j, e);

			/* Update vF value. */
			vF = simdui8_subs(vF, vGapE);
			vF = simdui8_max(vF, vH);
			/* Load the next vH. */
			vH = simdi_load(pvHLoad + j);
        }

		/* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
		/* reset pointers to the start of the saved data */
		j = 0;
		vH = simdi_load (pvHStore + j);
		/*  the computed vF value is for the given column.  since */
		/*  we are at the end, we need to shift the vF value over */
		/*  to the next column. */
		vF = simdi8_shiftl (vF, 1);
		vTemp = simdui8_subs(vH, vGapO);
		vTemp = simdui8_subs (vF, vTemp);
		while (simd_any(vTemp)) {
			vH = simdui8_max (vH, vF);
			vMaxColumn = simdui8_max(vMaxColumn, vH);
			simdi_store (pvHStore + j, vH);

			vF = simdui8_subs (vF, vGapE);
			j++;
			if (j >= segLen)
			{
				j = 0;
				vF = simdi8_shiftl (vF, 1);
			}
			vH = simdi_load (pvHStore + j);
			vTemp = simdui8_subs(vH, vGapO);
			vTemp = simdui8_subs (vF, vTemp);
		}

		vMaxScore = simdui8_max(vMaxScore, vMaxColumn);
		if (!simd_eq_all(vMaxMark, vMaxScore)) {
			uint8_t temp;
			vMaxMark = vMaxScore;
			temp = simdi8_hmax(vMaxScore);

			if (LIKELY(temp > max)) {
			    max = temp;
				if (max + bias >= 255) {
				    break;	//overflow
				}
				end_db = i;

				/* Store the column with the highest alignment score in order to trace the alignment ending position on read. */
				for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
			}
		}

        //uint8_t *t = (uint8_t *)pvHStore;
        //for (int ti = 0; ti < segLen * SIMD_SIZE; ++ti) {
        //    fprintf(stderr, "%d ", t[ti / segLen + ti % segLen * SIMD_SIZE]);
        //}
        //fprintf(stderr, "\n");

		/* Record the max score of current column. */
		maxColumn[i] = simdi8_hmax(vMaxColumn);
		//		fprintf(stderr, "maxColumn[%d]: %d\n", i, maxColumn[i]);
		if (maxColumn[i] == terminate) break;
	}

	/* Trace the alignment ending position on read. */
	uint8_t *t = (uint8_t*)pvHmax;
	int32_t column_len = segLen * SIMD_SIZE;
	for (i = 0; LIKELY(i < column_len); ++i, ++t) {
		int32_t temp;
		if (*t == max) {
			temp = i / SIMD_SIZE + i % SIMD_SIZE * segLen;
			if (temp < end_query) end_query = temp;
		}
	}

	/* Find the most possible 2nd best alignment. */
	alignment_end best0;
    best0.score = max + bias >= 255 ? 255 : max;
    best0.ref = end_db;
    best0.read = end_query;

    alignment_end best1;
    best1.score = 0;
    best1.ref = 0;
    best1.read = 0;

	edge = (end_db - maskLen) > 0 ? (end_db - maskLen) : 0;
	for (i = 0; i < edge; i ++) {
		if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
		}
	}
	edge = (end_db + maskLen) > db_length ? db_length : (end_db + maskLen);
	for (i = edge + 1; i < db_length; i ++) {
		if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
		}
	}
	return std::make_pair(best0, best1);
}

template <unsigned int type>
std::pair<SmithWaterman::alignment_end, SmithWaterman::alignment_end> SmithWaterman::sw_sse2_word (const unsigned char* db_sequence,
														   int8_t ref_dir,	// 0: forward ref; 1: reverse ref
														   int32_t db_length,
														   int32_t query_length,
														   const uint8_t gap_open, /* will be used as - */
														   const uint8_t gap_extend, /* will be used as - */
														   const simd_int* query_profile_word,
														   uint16_t terminate,
                                                           int32_t maskLen) {

	uint16_t max = 0;		                     /* the max alignment score */
	int32_t end_read = query_length - 1;
	int32_t end_ref = 0; /* 1_based best alignment ending point; Initialized as isn't aligned - 0. */
	const unsigned int SIMD_SIZE = VECSIZE_INT * 2;
	int32_t segLen = (query_length + SIMD_SIZE-1) / SIMD_SIZE; /* number of segment */
	/* array to record the alignment read ending position of the largest score of each reference position */
	memset(this->maxColumn, 0, db_length * sizeof(uint16_t));
	uint16_t * maxColumn = (uint16_t *) this->maxColumn;

	/* Define 16 byte 0 vector. */
	simd_int vZero = simdi32_set(0);
	simd_int* pvHStore = vHStore;
	simd_int* pvHLoad = vHLoad;
	simd_int* pvE = vE;
	simd_int* pvHmax = vHmax;
	memset(pvHStore,0,segLen*sizeof(simd_int));
	memset(pvHLoad,0, segLen*sizeof(simd_int));
	memset(pvE,0,     segLen*sizeof(simd_int));
	memset(pvHmax,0,  segLen*sizeof(simd_int));

	int32_t i, j, k;

    /* 16 byte insertion begin vector */
    simd_int vGapO = simdi16_set(gap_open);

	/* 16 byte insertion extension vector */
	simd_int vGapE = simdi16_set(gap_extend);

	simd_int vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
	simd_int vMaxMark = vZero; /* Trace the highest score till the previous column. */
	int32_t edge, begin = 0, end = db_length, step = 1;

    //fprintf(stderr, "start alignment of length %d [%d]\n", query_length, segLen * SIMD_SIZE);

	/* outer loop to process the reference sequence */
	if (ref_dir == 1) {
		begin = db_length - 1;
		end = -1;
		step = -1;
	}


	for (i = begin; LIKELY(i != end); i += step) {
		simd_int e, vF = vZero, vMaxColumn = vZero; /* Initialize F value to 0.
                                Any errors to vH values will be corrected in the Lazy_F loop.
                                */

		simd_int vH = pvHStore[segLen - 1];
		vH = simdi8_shiftl (vH, 2); /* Shift the 128-bit value in vH left by 2 byte. */
		const simd_int* vP = query_profile_word + db_sequence[i] * segLen; /* Right part of the query_profile_byte */

        /* Swap the 2 H buffers. */
        simd_int* pv = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pv;

		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); ++j) {
		    simd_int score = simdi_load(vP + j);
			vH = simdi16_adds(vH, score);

			/* Get max from vH, vE and vF. */
			e = simdi_load(pvE + j);
            vH = simdi16_max(vH, e);
			vH = simdi16_max(vH, vF);
			vMaxColumn = simdi16_max(vMaxColumn, vH);

			/* Save vH values. */
			simdi_store(pvHStore + j, vH);

			/* Update vE value. */
			vH = simdui16_subs(vH, vGapO); /* saturation arithmetic, result >= 0 */
			e = simdui16_subs(e, vGapE);
			e = simdi16_max(e, vH);
			simdi_store(pvE + j, e);

			/* Update vF value. */
			vF = simdui16_subs(vF, vGapE);
			vF = simdi16_max(vF, vH);

			/* Load the next vH. */
			vH = simdi_load(pvHLoad + j);
		}

		/* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
		for (k = 0; LIKELY(k < (int32_t) SIMD_SIZE); ++k) {
			vF = simdi8_shiftl (vF, 2);
			for (j = 0; LIKELY(j < segLen); ++j) {
				vH = simdi_load(pvHStore + j);
				vH = simdi16_max(vH, vF);
				vMaxColumn = simdi16_max(vMaxColumn, vH); //newly added line
				simdi_store(pvHStore + j, vH);
				vH = simdui16_subs(vH, vGapO);
				vF = simdui16_subs(vF, vGapE);
				if (UNLIKELY(!simd_any(simdi16_gt(vF, vH)))) goto end;
			}
		}

		end:
		vMaxScore = simdi16_max(vMaxScore, vMaxColumn);
		if (!simd_eq_all(vMaxMark, vMaxScore)) {
			uint16_t temp;
			vMaxMark = vMaxScore;
			temp = simdi16_hmax(vMaxScore);

			if (LIKELY(temp > max)) {
				max = temp;
				end_ref = i;
				for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
			}
		}

        //uint16_t *t = (uint16_t *)pvHStore;
        //for (size_t ti = 0; ti < segLen * SIMD_SIZE; ++ti) {
        //    fprintf(stderr, "%d ", t[ti / segLen + ti % segLen * SIMD_SIZE]);
        //}
        //fprintf(stderr, "\n");

		/* Record the max score of current column. */
		maxColumn[i] = simdi16_hmax(vMaxColumn);
		if (maxColumn[i] == terminate) break;
	}

	/* Trace the alignment ending position on read. */
	uint16_t *t = (uint16_t*)pvHmax;
	int32_t column_len = segLen * SIMD_SIZE;
	for (i = 0; LIKELY(i < column_len); ++i, ++t) {
		int32_t temp;
		if (*t == max) {
			temp = i / SIMD_SIZE + i % SIMD_SIZE * segLen;
			if (temp < end_read) end_read = temp;
		}
	}

	/* Find the most possible 2nd best alignment. */
	alignment_end best0;
    best0.score = max;
    best0.ref = end_ref;
    best0.read = end_read;

    alignment_end best1;
    best1.score = 0;
    best1.ref = 0;
    best1.read = 0;

	edge = (end_ref - maskLen) > 0 ? (end_ref - maskLen) : 0;
	for (i = 0; i < edge; i ++) {
		if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
		}
	}
	edge = (end_ref + maskLen) > db_length ? db_length : (end_ref + maskLen);
	for (i = edge; i < db_length; i ++) {
		if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
		}
	}

	return std::make_pair(best0, best1);
}

template <unsigned int type>
std::pair<SmithWaterman::alignment_end, SmithWaterman::alignment_end> SmithWaterman::sw_sse2_int (const unsigned char* db_sequence,
														   int8_t ref_dir,	// 0: forward ref; 1: reverse ref
														   int32_t db_length,
														   int32_t query_length,
														   const uint8_t gap_open, /* will be used as - */
														   const uint8_t gap_extend, /* will be used as - */
														   const simd_int* query_profile_int,
														   uint32_t terminate,
														   int32_t maskLen) {
#define max4(m, vm) ((m) = simdi32_hmax((vm)));

 	uint32_t max = 0; /* the max alignment score */
    int32_t end_read = query_length - 1;
    int32_t end_ref = 0; /* 1_based best alignment ending point; Initialized as isn't aligned - 0. */
    const unsigned int SIMD_SIZE = VECSIZE_INT;
    int32_t segLen = (query_length + SIMD_SIZE-1) / SIMD_SIZE; /* number of segment */
    /* array to record the alignment read ending position of the largest score of each reference position */
    memset(this->maxColumn, 0, db_length * sizeof(uint32_t));
    uint32_t * maxColumn = (uint32_t *) this->maxColumn;

    /* Define 16 byte 0 vector. */
    simd_int vZero = simdi32_set(0);
    simd_int* pvHStore = vHStore;
    simd_int* pvHLoad = vHLoad;
    simd_int* pvE = vE;
    simd_int* pvHmax = vHmax;
    memset(pvHStore, 0, segLen*sizeof(simd_int));
    memset(pvHLoad,  0, segLen*sizeof(simd_int));
    memset(pvE,      0, segLen*sizeof(simd_int));
    memset(pvHmax,   0, segLen*sizeof(simd_int));

    int32_t i, j, k;
    /* 16 byte insertion begin vector */
	simd_int vGapO = simdi32_set(gap_open);

	/* 16 byte insertion extension vector */
	simd_int vGapE = simdi32_set(gap_extend);

  	simd_int vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
    simd_int vMaxMark = vZero; /* Trace the highest score till the previous column. */
    simd_int vTemp;
    int32_t edge, begin = 0, end = db_length, step = 1;

    /* outer loop to process the reference sequence */
    if (ref_dir == 1) {
        begin = db_length - 1;
        end = -1;
        step = -1;
    }


	for (i = begin; LIKELY(i != end); i += step) {
		simd_int e, vF = vZero;
		simd_int vMaxColumn = vZero; /* Initialize F value to 0.
                                Any errors to vH values will be corrected in the Lazy_F loop.
                                */
		simd_int vH = pvHStore[segLen - 1];
        vH = simdi8_shiftl(vH, 4); /* Shift the 128-bit value in vH left by 4 byte. */

		/* Swap the 2 H buffers. */
        simd_int* pv = pvHLoad;

		const simd_int* vP = query_profile_int + db_sequence[i] * segLen; /* Right part of the query_profile_byte */

		pvHLoad = pvHStore;
		pvHStore = pv;

		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); ++j) {
		    simd_int score = simdi_load(vP + j);
			// vH = simdi32_adds(vH, score);
			vH = simdi32_add(vH, score);
			/* Get max from vH, vE and vF. */
			e = simdi_load(pvE + j);
            vH = simdi32_max(vH, e);
			vH = simdi32_max(vH, vF);

			vMaxColumn = simdi32_max(vMaxColumn, vH);

			/* Save vH values. */
			simdi_store(pvHStore + j, vH);

			/* Update vE value. */
			vH = simdui32_subs(vH, vGapO); /* saturation arithmetic, result >= 0 */
			e = simdui32_subs(e, vGapE);
			e = simdi32_max(e, vH);
			simdi_store(pvE + j, e);

			/* Update vF value. */
			vF = simdui32_subs(vF, vGapE);
			vF = simdi32_max(vF, vH);

			/* Load the next vH. */
			vH = simdi_load(pvHLoad + j);
		}

		/* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
        for (k = 0; LIKELY(k < (int32_t) SIMD_SIZE); ++k) {
            vF = simdi8_shiftl(vF, 4);
            for (j = 0; LIKELY(j < segLen); ++j) {
                vH = simdi_load(pvHStore + j);
				vH = simdi32_max(vH, vF);
                vMaxColumn = simdi32_max(vMaxColumn, vH); //newly added line
                simdi_store(pvHStore + j, vH);
				vH = simdui32_subs(vH, vGapO);
				vF = simdui32_subs(vF, vGapE);
				if (UNLIKELY(! simdi8_movemask(simdi32_gt(vF, vH)))) goto end;
			}
		}
		end:
        vMaxScore = simdi32_max(vMaxScore, vMaxColumn);
        vTemp = simdi32_eq(vMaxMark, vMaxScore);
        uint32_t cmp = simdi8_movemask(vTemp);
        if (cmp != SIMD_MOVEMASK_MAX) {
            uint32_t temp;
            vMaxMark = vMaxScore;
            max4(temp, vMaxScore);
            vMaxScore = vMaxMark;

            if (LIKELY(temp > max)) {
                max = temp;
                end_ref = i;
                for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
            }
        }

        /* Record the max score of current column. */
        max4(maxColumn[i], vMaxColumn);
        if (maxColumn[i] == terminate) break;
	}

    /* Trace the alignment ending position on read. */
    uint32_t *t = (uint32_t*)pvHmax;
    int32_t column_len = segLen * SIMD_SIZE;
    for (i = 0; LIKELY(i < column_len); ++i, ++t) {
        int32_t temp;
        if (*t == max) {
            temp = i / SIMD_SIZE + i % SIMD_SIZE * segLen;
            if (temp < end_read) end_read = temp;
        }
    }

    /* Find the most possible 2nd best alignment. */
    alignment_end best0;
    best0.score = max;
    best0.ref = end_ref;
    best0.read = end_read;

    alignment_end best1;
    best1.score = 0;
    best1.ref = 0;
    best1.read = 0;

    edge = (end_ref - maskLen) > 0 ? (end_ref - maskLen) : 0;
    for (i = 0; i < edge; i ++) {
        if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
        }
    }
    edge = (end_ref + maskLen) > db_length ? db_length : (end_ref + maskLen);
    for (i = edge; i < db_length; i ++) {
        if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
        }
    }

    return std::make_pair(best0, best1);
#undef max4
}

void SmithWaterman::ssw_init(const Sequence* q,
							 const int8_t* mat,
							 const BaseMatrix *m) {
	//init profile
    profile->bias = 0;
	profile->query_length = q->L;
	const int32_t alphabetSize = m->alphabetSize;
    profile->alphabetSize = m->alphabetSize;

	bool isProfile = Parameters::isEqualDbtype(q->getSequenceType(), Parameters::DBTYPE_HMM_PROFILE);
	profile->isProfile = isProfile;
	int32_t compositionBias = 0;
	if (!isProfile && aaBiasCorrection) {
		SubstitutionMatrix::calcLocalAaBiasCorrection(m, q->numSequence, q->L, tmp_composition_bias, aaBiasCorrectionScale);
		for (int i =0; i < q->L; i++) {
			profile->composition_bias[i] = (int8_t) (tmp_composition_bias[i] < 0.0)? tmp_composition_bias[i] - 0.5: tmp_composition_bias[i] + 0.5;
			compositionBias = (compositionBias < profile->composition_bias[i]) ? compositionBias : profile->composition_bias[i];
		}
		compositionBias = std::min(compositionBias, 0);
	} else {
		memset(profile->composition_bias, 0, q->L* sizeof(int8_t));
	}

	// copy memory to local memory
	if (isProfile) {
		memcpy(profile->mat, mat, q->L * Sequence::PROFILE_AA_SIZE * sizeof(int8_t));
		// set neutral state 'X' (score=0)
		memset(profile->mat + ((alphabetSize - 1) * q->L), 0, q->L * sizeof(int8_t));
	} else {
		memcpy(profile->mat, mat, alphabetSize * alphabetSize * sizeof(int8_t));
	}
	memcpy(profile->query_sequence, q->numSequence, q->L);

	int32_t bias = 0;
	int32_t matSize = isProfile ? q->L * Sequence::PROFILE_AA_SIZE : alphabetSize * alphabetSize;

    for (int32_t i = 0; i < matSize; i++) {
        if (mat[i] < bias) {
            bias = mat[i];
        }
    }
    bias = abs(bias) + abs(compositionBias);
    profile->bias = bias;

    if (isProfile) {
        // offset = 1 when createQueryProfile
        // create byte version of profiles
        createQueryProfile<int8_t, VECSIZE_INT * 4, PROFILE>(profile->profile_byte, profile->query_sequence, NULL,
                                                            profile->mat, q->L, alphabetSize, profile->bias, 0, q->L);
		// create word version of profiles
		createQueryProfile<int16_t, VECSIZE_INT * 2, PROFILE>(profile->profile_word, profile->query_sequence, NULL,
															profile->mat, q->L, alphabetSize, 0, 0, q->L);
		// create int version of profiles
		createQueryProfile<int32_t, VECSIZE_INT * 1, PROFILE>(profile->profile_int, profile->query_sequence, NULL,
															profile->mat, q->L, alphabetSize, 0, 0, q->L);
 	  	// create linear version of word profile
    	for (int32_t i = 0; i< alphabetSize; i++) {
            profile->profile_word_linear[i] = &profile_word_linear_data[i*q->L];
            for (int j = 0; j < q->L; j++) {
                profile->profile_word_linear[i][j] = mat[i * q->L + j];
            }
        }
    } else {
        // create byte version of query profile
        createQueryProfile<int8_t, VECSIZE_INT * 4, SUBSTITUTIONMATRIX>(profile->profile_byte, profile->query_sequence, profile->composition_bias, profile->mat, q->L, alphabetSize, bias, 0, 0);
        // create word version of query profile
        createQueryProfile<int16_t, VECSIZE_INT * 2, SUBSTITUTIONMATRIX>(profile->profile_word, profile->query_sequence, profile->composition_bias, profile->mat, q->L, alphabetSize, 0, 0, 0);
        // create int version of query profile
		createQueryProfile<int32_t, VECSIZE_INT * 1, SUBSTITUTIONMATRIX>(profile->profile_int, profile->query_sequence, profile->composition_bias, profile->mat, q->L, alphabetSize, 0, 0, 0);
		// create linear version of word profile
        for (int32_t i = 0; i< alphabetSize; i++) {
            profile->profile_word_linear[i] = &profile_word_linear_data[i*q->L];
            for (int j = 0; j < q->L; j++) {
                profile->profile_word_linear[i][j] = mat[i * alphabetSize + q->numSequence[j]] + profile->composition_bias[j];
            }
        }
    }

	// create reverse structures
	std::reverse_copy(profile->query_sequence, profile->query_sequence + q->L, profile->query_rev_sequence);
	std::reverse_copy(profile->composition_bias, profile->composition_bias + q->L, profile->composition_bias_rev);

	if (isProfile) {
		for (int32_t i = 0; i < alphabetSize; i++) {
            const int8_t *startToRead = profile->mat + (i * q->L);
            int8_t *startToWrite = profile->mat_rev + (i * q->L);
            std::reverse_copy(startToRead, startToRead + q->L, startToWrite);
        }
		memset(profile->pos_aa_rev, 0x80, q->L*32); // do we need it?
		int rowIdx = 0;
		for (int i = 0; i < profile->query_length; i++) {
			for (int aa = 0; aa < m->alphabetSize; aa++) {
				int score = profile->mat_rev[aa * profile->query_length + i];
				// int idx = rowIdx + (m->num2aa[aa] - 'A'); //orig
				int idx = rowIdx + aa; //new
				profile->pos_aa_rev[idx] = static_cast<int8_t>(score);
			}
			rowIdx += 32;
		}
	}
	else {
		for (int i = 0; i < q->L; i++) {
			block->query_bias_arr[i] = profile->composition_bias_rev[i];
		}

		for (int aa1 = 0; aa1 < m->alphabetSize; aa1++) {
			for (int aa2 = 0; aa2 < m->alphabetSize; aa2++) {
				// instead of num2aa, use aa directly
				block_set_aamatrix_num(block->mat_aa, aa1, aa2, m->subMatrix[aa1][aa2]);
			}
		}
	}
}

template <unsigned int type>
SmithWaterman::cigar * SmithWaterman::banded_sw(const unsigned char *db_sequence, const int8_t *query_sequence,
                                                const int8_t * compositionBias,
												int32_t db_length, int32_t query_length, int32_t queryStart,
												int32_t score, const uint32_t gap_open,
												const uint32_t gap_extend,
												int32_t band_width, const int8_t *mat,
												const int32_t qry_n) {
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

    /* Convert the coordinate in the scoring matrix into the coordinate in one line of the band. */
#define set_u(u, w, i, j) { int x=(i)-(w); x=x>0?x:0; (u)=(j)-x+1; }

    /* Convert the coordinate in the direction matrix into the coordinate in one line of the band. */
#define set_d(u, w, i, j, p) { int x=(i)-(w); x=x>0?x:0; x=(j)-x; (u)=x*3+p; }

    uint32_t *c = (uint32_t*)malloc(16 * sizeof(uint32_t)), *c1;
    int32_t i, j, e, f, temp1, temp2, s = 16, s1 = 8, l, max = 0;
    int64_t s2 = 1024;
    char op, prev_op;
    int64_t width, width_d;
    int32_t *h_b, *e_b, *h_c;
    int8_t *direction, *direction_line;
    cigar* result = new cigar();
    h_b = (int32_t*)malloc(s1 * sizeof(int32_t));
    e_b = (int32_t*)malloc(s1 * sizeof(int32_t));
    h_c = (int32_t*)malloc(s1 * sizeof(int32_t));
    direction = (int8_t*)malloc(s2 * sizeof(int8_t));

	do {
		width = band_width * 2 + 3, width_d = band_width * 2 + 1;
		while (width >= s1) {
			++s1;
			kroundup32(s1);
			h_b = (int32_t*)realloc(h_b, s1 * sizeof(int32_t));
			e_b = (int32_t*)realloc(e_b, s1 * sizeof(int32_t));
			h_c = (int32_t*)realloc(h_c, s1 * sizeof(int32_t));
		}
		int64_t targetSize = width_d * query_length * 3;
		while (targetSize >= s2) {
			++s2;
			kroundup32(s2);
			if (s2 < 0) {
				fprintf(stderr, "Alignment score and position are not consensus.\n");
				EXIT(1);
			}
			direction = (int8_t*)realloc(direction, s2 * sizeof(int8_t));
		}
		direction_line = direction;
		for (j = 1; LIKELY(j < width - 1); j ++) {
		    h_b[j] = 0;
		}
		for (i = 0; LIKELY(i < query_length); i ++) {
			int32_t beg = 0, end = db_length - 1, u = 0, edge;
			j = i - band_width;	beg = beg > j ? beg : j; // band start
			j = i + band_width; end = end < j ? end : j; // band end
			edge = end + 1 < width - 1 ? end + 1 : width - 1;
			f = h_b[0] = e_b[0] = h_b[edge] = e_b[edge] = h_c[0] = 0;
			int64_t directionOffset = width_d * i * 3;
			direction_line = direction + directionOffset;

			for (j = beg; LIKELY(j <= end); j ++) {
				int32_t b, e1, f1, d, de, df, dh;
				set_u(u, band_width, i, j);
				set_u(e, band_width, i - 1, j);
				set_u(b, band_width, i, j - 1);
				set_u(d, band_width, i - 1, j - 1);
				set_d(de, band_width, i, j, 0);
				set_d(df, band_width, i, j, 1);
				set_d(dh, band_width, i, j, 2);

				temp1 = i == 0 ? -gap_open : h_b[e] - gap_open;
                temp2 = i == 0 ? -gap_extend : e_b[e] - gap_extend;
				e_b[u] = temp1 > temp2 ? temp1 : temp2;
				direction_line[de] = temp1 > temp2 ? 3 : 2;

				temp1 = h_c[b] - gap_open;
                temp2 = f - gap_extend;
				f = temp1 > temp2 ? temp1 : temp2;
				direction_line[df] = temp1 > temp2 ? 5 : 4;

                f1 = f > 0 ? f : 0;
				e1 = e_b[u] > 0 ? e_b[u] : 0;

				temp1 = e1 > f1 ? e1 : f1;

				//TODO: careful with the variable names
				if (type == PROFILE_SEQ) {
				    // db_sequence is a numerical sequence
                    temp2 = h_b[d] + mat[db_sequence[j] * qry_n + (queryStart + i)];
                } else {
                    temp2 = h_b[d] + mat[query_sequence[i] * qry_n + db_sequence[j]] + compositionBias[i];
				}

				h_c[u] = temp1 > temp2 ? temp1 : temp2;
				if (h_c[u] > max) max = h_c[u];

				if (temp1 <= temp2) direction_line[dh] = 1;
				else direction_line[dh] = e1 > f1 ? direction_line[de] : direction_line[df];
			}
			for (j = 1; j <= u; j ++) {
			    h_b[j] = h_c[j];
			    //fprintf(stderr, "%d ", h_b[j]);
			}
			//fprintf(stderr, "\n");
		}
        //TODO make band_width dependet on how far the score was off
		band_width *= 2;
	} while (LIKELY(max < score));
	band_width /= 2;


	// trace back
	i = query_length - 1;
	j = db_length - 1;
	e = 0;	// Count the number of M, D or I.
	l = 0;	// record length of current cigar
	op = prev_op = 'M';
	temp2 = 2;	// h
	while (LIKELY(i > 0) || LIKELY(j > 0)) {
		set_d(temp1, band_width, i, j, temp2);
		switch (direction_line[temp1]) {
			case 1:
				--i;
				--j;
				temp2 = 2;
				direction_line -= width_d * 3;
				op = 'M';
				break;
			case 2:
				--i;
				temp2 = 0;	// e
				direction_line -= width_d * 3;
				op = 'I';
				break;
			case 3:
				--i;
				temp2 = 2;
				direction_line -= width_d * 3;
				op = 'I';
				break;
			case 4:
				--j;
				temp2 = 1;
				op = 'D';
				break;
			case 5:
				--j;
				temp2 = 2;
				op = 'D';
				break;
			default:
				fprintf(stderr, "Trace back error: %d.\n", direction_line[temp1 - 1]);
				free(direction);
				free(h_c);
				free(e_b);
				free(h_b);
				free(c);
				delete result;
				return 0;
		}
		if (op == prev_op) ++e;
		else {
			++l;
			while (l >= s) {
				++s;
				kroundup32(s);
				c = (uint32_t*)realloc(c, s * sizeof(uint32_t));
			}
			c[l - 1] = to_cigar_int(e, prev_op);
			prev_op = op;
			e = 1;
		}
	}
	if (op == 'M') {
		++l;
		while (l >= s) {
			++s;
			kroundup32(s);
			c = (uint32_t*)realloc(c, s * sizeof(uint32_t));
		}
		c[l - 1] = to_cigar_int(e + 1, op);
	}else {
		l += 2;
		while (l >= s) {
			++s;
			kroundup32(s);
			c = (uint32_t*)realloc(c, s * sizeof(uint32_t));
		}
		c[l - 2] = to_cigar_int(e, op);
		c[l - 1] = to_cigar_int(1, 'M');
	}

	// reverse cigar
	c1 = (uint32_t*)new uint32_t[l * sizeof(uint32_t)];
	s = 0;
	e = l - 1;
	while (LIKELY(s <= e)) {
		c1[s] = c[e];
		c1[e] = c[s];
		++ s;
		-- e;
	}
	result->seq = c1;
	result->length = l;

	free(direction);
	free(h_c);
	free(e_b);
	free(h_b);
	free(c);
	return result;
#undef kroundup32
#undef set_u
#undef set_d
}

uint32_t SmithWaterman::to_cigar_int (uint32_t length, char op_letter)
{
	uint32_t res;
	uint8_t op_code;

	switch (op_letter) {
		case 'M': /* alignment match (can be a sequence match or mismatch */
		default:
			op_code = 0;
			break;
		case 'I': /* insertion to the reference */
			op_code = 1;
			break;
		case 'D': /* deletion from the reference */
			op_code = 2;
			break;
		case 'N': /* skipped region from the reference */
			op_code = 3;
			break;
		case 'S': /* soft clipping (clipped sequences present in SEQ) */
			op_code = 4;
			break;
		case 'H': /* hard clipping (clipped sequences NOT present in SEQ) */
			op_code = 5;
			break;
		case 'P': /* padding (silent deletion from padded reference) */
			op_code = 6;
			break;
		case '=': /* sequence match */
			op_code = 7;
			break;
		case 'X': /* sequence mismatch */
			op_code = 8;
			break;
	}

	res = (length << 4) | op_code;
	return res;
}

void SmithWaterman::printVector(__m128i v){
	for (int i = 0; i < 8; i++)
		printf("%d ", ((short) (sse2_extract_epi16(v, i)) + 32768));
	std::cout << "\n";
}

void SmithWaterman::printVectorUS(__m128i v){
	for (int i = 0; i < 8; i++)
		printf("%d ", (unsigned short) sse2_extract_epi16(v, i));
	std::cout << "\n";
}

unsigned short SmithWaterman::sse2_extract_epi16(__m128i v, int pos) {
	switch(pos){
		case 0: return _mm_extract_epi16(v, 0);
		case 1: return _mm_extract_epi16(v, 1);
		case 2: return _mm_extract_epi16(v, 2);
		case 3: return _mm_extract_epi16(v, 3);
		case 4: return _mm_extract_epi16(v, 4);
		case 5: return _mm_extract_epi16(v, 5);
		case 6: return _mm_extract_epi16(v, 6);
		case 7: return _mm_extract_epi16(v, 7);
	}
	std::cerr << "Fatal error in QueryScore: position in the vector is not in the legal range (pos = " << pos << ")\n";
	EXIT(1);
	// never executed
	return 0;
}

float SmithWaterman::computeCov(unsigned int startPos, unsigned int endPos, unsigned int len) {
	return (std::min(len, std::max(startPos, endPos)) - std::min(startPos, endPos) + 1) / (float) len;
}

s_align SmithWaterman::scoreIdentical(unsigned char *dbSeq, int L, EvalueComputation * evaluer,
                                      int alignmentMode, std::string &backtrace) {
	if(profile->query_length != L){
		std::cerr << "scoreIdentical has different length L: "
				  << L << " query_length: " << profile->query_length
				  << "\n";
		EXIT(1);
	}

	s_align r;
	// to be compatible with --alignment-mode 1 (score only)
	if(alignmentMode == 0){
		r.dbStartPos1 = -1;
		r.qStartPos1 = -1;
	}else{
		r.qStartPos1 = 0;
		r.dbStartPos1 = 0;
	}

	r.qEndPos1 = L -1;
	r.dbEndPos1 = L -1;
	r.cigarLen = L;
	r.qCov =  1.0;
	r.tCov = 1.0;
    r.cigar = NULL;
	short score = 0;
	for(int pos = 0; pos < L; pos++){
		int currScore = profile->profile_word_linear[dbSeq[pos]][pos];
		score += currScore;
        backtrace.push_back('M');
	}
	r.score1=score;
	r.evalue = evaluer->computeEvalue(r.score1, profile->query_length);
    r.identicalAACnt = L;
	return r;
}

template <typename F>
inline F simd_hmax(const F * in, unsigned int n) {
    F current = std::numeric_limits<F>::min();
    do {
        current = std::max(current, *in++);
    } while(--n);

    return current;
}

int SmithWaterman::ungapped_alignment(const unsigned char *db_sequence, int32_t db_length) {
#define SWAP(tmp, arg1, arg2) tmp = arg1; arg1 = arg2; arg2 = tmp;

    int i; // position in query bands (0,..,W-1)
    int j; // position in db sequence (0,..,dbseq_length-1)
    int element_count = (VECSIZE_INT * 4);
    const int W = (profile->query_length + (element_count - 1)) /
                  element_count; // width of bands in query and score matrix = hochgerundetes LQ/16

    simd_int *p;
    simd_int S;              // 16 unsigned bytes holding S(b*W+i,j) (b=0,..,15)
    simd_int Smax = simdi_setzero();
    simd_int Soffset; // all scores in query profile are shifted up by Soffset to obtain pos values
    simd_int *s_prev, *s_curr; // pointers to Score(i-1,j-1) and Score(i,j), resp.
    simd_int *qji;             // query profile score in row j (for residue x_j)
    simd_int *s_prev_it, *s_curr_it;
    simd_int *query_profile_it = (simd_int *) profile->profile_byte;

    // Load the score offset to all 16 unsigned byte elements of Soffset
    Soffset = simdi8_set(profile->bias);
    s_curr = vHStore;
    s_prev = vHLoad;

    memset(vHStore, 0, W * sizeof(simd_int));
    memset(vHLoad, 0, W * sizeof(simd_int));

    for (j = 0; j < db_length; ++j) // loop over db sequence positions
    {

        // Get address of query scores for row j
        qji = query_profile_it + db_sequence[j] * W;

        // Load the next S value
        S = simdi_load(s_curr + W - 1);
        S = simdi8_shiftl(S, 1);

        // Swap s_prev and s_curr, smax_prev and smax_curr
        SWAP(p, s_prev, s_curr);

        s_curr_it = s_curr;
        s_prev_it = s_prev;

        for (i = 0; i < W; ++i) // loop over query band positions
        {
            // Saturated addition and subtraction to score S(i,j)
            S = simdui8_adds(S, *(qji++)); // S(i,j) = S(i-1,j-1) + (q(i,x_j) + Soffset)
            S = simdui8_subs(S, Soffset);       // S(i,j) = max(0, S(i,j) - Soffset)
            simdi_store(s_curr_it++, S);       // store S to s_curr[i]
            Smax = simdui8_max(Smax, S);       // Smax(i,j) = max(Smax(i,j), S(i,j))

            // Load the next S and Smax values
            S = simdi_load(s_prev_it++);
        }
    }
    int score = simd_hmax((unsigned char *) &Smax, element_count);

    /* return largest score */
    return score;
#undef SWAP
}
