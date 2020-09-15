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

#include "Util.h"
#include "SubstitutionMatrix.h"
#include "Debug.h"

#include <iostream>

SmithWaterman::SmithWaterman(size_t maxSequenceLength, int aaSize, bool aaBiasCorrection) {
	maxSequenceLength += 1;
	this->aaBiasCorrection = aaBiasCorrection;
	const int segSize = (maxSequenceLength+7)/8;
	vHStore = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	vHLoad  = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	vE      = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	vHmax   = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	profile = new s_profile();
	profile->profile_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
	profile->profile_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
	profile->profile_rev_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
	profile->profile_rev_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
	profile->query_rev_sequence = new int8_t[maxSequenceLength];
	profile->query_sequence     = new int8_t[maxSequenceLength];
	profile->composition_bias   = new int8_t[maxSequenceLength];
	profile->composition_bias_rev   = new int8_t[maxSequenceLength];
	profile->profile_word_linear = new short*[aaSize];
	profile_word_linear_data = new short[aaSize*maxSequenceLength];
	profile->mat_rev            = new int8_t[maxSequenceLength * aaSize * 2];
	profile->mat                = new int8_t[maxSequenceLength * aaSize * 2];
	tmp_composition_bias   = new float[maxSequenceLength];
	/* array to record the largest score of each reference position */
	maxColumn = new uint8_t[maxSequenceLength*sizeof(uint16_t)];
	memset(maxColumn, 0, maxSequenceLength*sizeof(uint16_t));

	memset(profile->query_sequence, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->query_rev_sequence, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->mat_rev, 0, maxSequenceLength * aaSize);
	memset(profile->composition_bias, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->composition_bias_rev, 0, maxSequenceLength * sizeof(int8_t));
}

SmithWaterman::~SmithWaterman(){
	free(vHStore);
	free(vHLoad);
	free(vE);
	free(vHmax);
	free(profile->profile_byte);
	free(profile->profile_word);
	free(profile->profile_rev_byte);
	free(profile->profile_rev_word);
	delete [] profile->query_rev_sequence;
	delete [] profile->query_sequence;
	delete [] profile->composition_bias;
	delete [] profile->composition_bias_rev;
	delete [] profile->profile_word_linear;
	delete [] profile_word_linear_data;
	delete [] profile->mat_rev;
	delete [] profile->mat;
	delete [] tmp_composition_bias;
	delete [] maxColumn;
	delete profile;
}


/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch. */
template <typename T, size_t Elements, const unsigned int type>
void SmithWaterman::createQueryProfile(simd_int *profile, const int8_t *query_sequence, const int8_t * composition_bias, const int8_t *mat,
									   const int32_t query_length, const int32_t aaSize, uint8_t bias,
									   const int32_t offset, const int32_t entryLength) {

	const int32_t segLen = (query_length+Elements-1)/Elements;
	T* t = (T*)profile;

	/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch */
	for (int32_t nt = 0; LIKELY(nt < aaSize); nt++) {
//		printf("{");
		for (int32_t i = 0; i < segLen; i ++) {
			int32_t  j = i;
//			printf("(");
			for (size_t segNum = 0; LIKELY(segNum < Elements) ; segNum ++) {
				// if will be optmized out by compiler
				if(type == SUBSTITUTIONMATRIX) {     // substitution score for query_seq constrained by nt
					// query_sequence starts from 1 to n
					*t++ = ( j >= query_length) ? bias : mat[nt * aaSize + query_sequence[j + offset ]] + composition_bias[j + offset] + bias; // mat[nt][q[j]] mat eq 20*20
//					printf("(%1d, %1d) ", query_sequence[j ], *(t-1));

				} if(type == PROFILE) {
					// profile starts by 0
					*t++ = ( j >= query_length) ? bias : mat[nt * entryLength  + (j + (offset - 1) )] + bias; //mat eq L*20  // mat[nt][j]
//					printf("(%1d, %1d) ", j , *(t-1));
				}
				j += segLen;
			}
//			printf(")");
		}
//		printf("}\n");
	}
//	printf("\n");
//	std::flush(std::cout);

}


s_align SmithWaterman::ssw_align (
		const unsigned char *db_sequence,
		int32_t db_length,
		const uint8_t gap_open,
		const uint8_t gap_extend,
		const uint8_t alignmentMode,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
		const double  evalueThr,
		EvalueComputation * evaluer,
		const int covMode, const float covThr,
		const int32_t maskLen) {

	int32_t word = 0, query_length = profile->query_length;
	int32_t band_width = 0;
	cigar* path;
	s_align r;
	r.dbStartPos1 = -1;
	r.qStartPos1 = -1;
	r.cigar = 0;
	r.cigarLen = 0;
	//if (maskLen < 15) {
	//	fprintf(stderr, "When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.\n");
	//}

    std::pair<alignment_end, alignment_end> bests;
    std::pair<alignment_end, alignment_end> bests_reverse;
    // Find the alignment scores and ending positions
	if (profile->profile_byte) {
		bests = sw_sse2_byte(db_sequence, 0, db_length, query_length, gap_open, gap_extend, profile->profile_byte, UCHAR_MAX, profile->bias, maskLen);

		if (profile->profile_word && bests.first.score == 255) {
			bests = sw_sse2_word(db_sequence, 0, db_length, query_length, gap_open, gap_extend, profile->profile_word, USHRT_MAX, maskLen);
			word = 1;
		} else if (bests.first.score == 255) {
			fprintf(stderr, "Please set 2 to the score_size parameter of the function ssw_init, otherwise the alignment results will be incorrect.\n");
			EXIT(EXIT_FAILURE);
		}
	}else if (profile->profile_word) {
		bests = sw_sse2_word(db_sequence, 0, db_length, query_length, gap_open, gap_extend, profile->profile_word, USHRT_MAX, maskLen);
		word = 1;
	}else {
		fprintf(stderr, "Please call the function ssw_init before ssw_align.\n");
		EXIT(EXIT_FAILURE);
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

    const bool isProfile = Parameters::isEqualDbtype(profile->sequence_type, Parameters::DBTYPE_HMM_PROFILE)
                         || Parameters::isEqualDbtype(profile->sequence_type, Parameters::DBTYPE_PROFILE_STATE_PROFILE);
    // no residue could be aligned
    if (r.dbEndPos1 == -1) {
        return r;
    }
    int32_t queryOffset = query_length - r.qEndPos1;
	r.evalue = evaluer->computeEvalue(r.score1, query_length);
    bool hasLowerEvalue = r.evalue > evalueThr;
	r.qCov = computeCov(0, r.qEndPos1, query_length);
	r.tCov = computeCov(0, r.dbEndPos1, db_length);
    bool hasLowerCoverage = !(Util::hasCoverage(covThr, covMode, r.qCov, r.tCov));

	if (alignmentMode == 0 || ((alignmentMode == 2 || alignmentMode == 1) && (hasLowerEvalue || hasLowerCoverage))) {
        return r;
	}

	// Find the beginning position of the best alignment.
	if (word == 0) {
		if (isProfile) {
			createQueryProfile<int8_t, VECSIZE_INT * 4, PROFILE>(profile->profile_rev_byte, profile->query_rev_sequence, NULL, profile->mat_rev,
																 r.qEndPos1 + 1, profile->alphabetSize, profile->bias, queryOffset, profile->query_length);
		} else {
			createQueryProfile<int8_t, VECSIZE_INT * 4, SUBSTITUTIONMATRIX>(profile->profile_rev_byte, profile->query_rev_sequence, profile->composition_bias_rev, profile->mat,
																			r.qEndPos1 + 1, profile->alphabetSize, profile->bias, queryOffset, 0);
		}
		bests_reverse = sw_sse2_byte(db_sequence, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1, gap_open, gap_extend, profile->profile_rev_byte,
									 r.score1, profile->bias, maskLen);
	} else {
		if (isProfile) {
			createQueryProfile<int16_t, VECSIZE_INT * 2, PROFILE>(profile->profile_rev_word, profile->query_rev_sequence, NULL, profile->mat_rev,
																  r.qEndPos1 + 1, profile->alphabetSize, 0, queryOffset, profile->query_length);

		} else {
			createQueryProfile<int16_t, VECSIZE_INT * 2, SUBSTITUTIONMATRIX>(profile->profile_rev_word, profile->query_rev_sequence, profile->composition_bias_rev, profile->mat,
																			 r.qEndPos1 + 1, profile->alphabetSize, 0, queryOffset, 0);
		}
		bests_reverse = sw_sse2_word(db_sequence, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1, gap_open, gap_extend, profile->profile_rev_word,
									 r.score1, maskLen);
	}
	if(bests_reverse.first.score != r.score1){
		fprintf(stderr, "Score of forward/backward SW differ. This should not happen.\n");
		EXIT(EXIT_FAILURE);
	}

	r.dbStartPos1 = bests_reverse.first.ref;
	r.qStartPos1 = r.qEndPos1 - bests_reverse.first.read;

    if (r.dbStartPos1 == -1) {
        fprintf(stderr, "Target start position is -1. This should not happen.\n");
        EXIT(EXIT_FAILURE);
    }

	r.qCov = computeCov(r.qStartPos1, r.qEndPos1, query_length);
	r.tCov = computeCov(r.dbStartPos1, r.dbEndPos1, db_length);
	hasLowerCoverage = !(Util::hasCoverage(covThr, covMode, r.qCov, r.tCov));
    // only start and end point are needed
    if (alignmentMode == 1 || hasLowerCoverage) {
        return r;
    }

	// Generate cigar.
	db_length = r.dbEndPos1 - r.dbStartPos1 + 1;
	query_length = r.qEndPos1 - r.qStartPos1 + 1;
	band_width = abs(db_length - query_length) + 1;

	if (isProfile) {
		path = banded_sw<PROFILE>(db_sequence + r.dbStartPos1, profile->query_sequence + r.qStartPos1,
								  NULL, db_length, query_length,
								  r.qStartPos1, r.score1, gap_open, gap_extend, band_width,
								  profile->mat, profile->query_length);
	} else {
		path = banded_sw<SUBSTITUTIONMATRIX>(db_sequence + r.dbStartPos1,
											 profile->query_sequence + r.qStartPos1,
											 profile->composition_bias + r.qStartPos1,
											 db_length, query_length, r.qStartPos1, r.score1,
											 gap_open, gap_extend, band_width,
											 profile->mat, profile->alphabetSize);
	}
	if (path != NULL) {
		r.cigar = path->seq;
		r.cigarLen = path->length;
	}	delete path;

	return r;
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

std::pair<SmithWaterman::alignment_end, SmithWaterman::alignment_end> SmithWaterman::sw_sse2_byte (const unsigned char* db_sequence,
														   int8_t ref_dir,	// 0: forward ref; 1: reverse ref
														   int32_t db_length,
														   int32_t query_length,
														   const uint8_t gap_open, /* will be used as - */
														   const uint8_t gap_extend, /* will be used as - */
														   const simd_int* query_profile_byte,
														   uint8_t terminate,	/* the best alignment score: used to terminate
                                                         the matrix calculation when locating the
                                                         alignment beginning point. If this score
                                                         is set to 0, it will not be used */
														   uint8_t bias,  /* Shift 0 point to a positive value. */
														   int32_t maskLen) {
#define max16(m, vm) ((m) = simdi8_hmax((vm)));

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
	//	int32_t distance = query_length * 2 / 3;
	//	int32_t distance = query_length / 2;
	//	int32_t distance = query_length;

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
		//		max16(maxColumn[i], vMaxColumn);
		//		fprintf(stderr, "middle[%d]: %d\n", i, maxColumn[i]);

		simd_int vH = pvHStore[segLen - 1];
		vH = simdi8_shiftl (vH, 1); /* Shift the 128-bit value in vH left by 1 byte. */
		const simd_int* vP = query_profile_byte + db_sequence[i] * segLen; /* Right part of the query_profile_byte */
		//	int8_t* t;
		//	int32_t ti;
		//        fprintf(stderr, "i: %d of %d:\t ", i,segLen);
		//for (t = (int8_t*)vP, ti = 0; ti < segLen; ++ti) fprintf(stderr, "%d\t", *t++);
		//fprintf(stderr, "\n");

		/* Swap the 2 H buffers. */
		simd_int* pv = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pv;

		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); ++j) {
			vH = simdui8_adds(vH, simdi_load(vP + j));
			vH = simdui8_subs(vH, vBias); /* vH will be always > 0 */
			//	max16(maxColumn[i], vH);
			//	fprintf(stderr, "H[%d]: %d\n", i, maxColumn[i]);
			//	int8_t* t;
			//	int32_t ti;
			//for (t = (int8_t*)&vH, ti = 0; ti < 16; ++ti) fprintf(stderr, "%d\t", *t++);

			/* Get max from vH, vE and vF. */
			e = simdi_load(pvE + j);
			vH = simdui8_max(vH, e);
			vH = simdui8_max(vH, vF);
			vMaxColumn = simdui8_max(vMaxColumn, vH);

			//	max16(maxColumn[i], vMaxColumn);
			//	fprintf(stderr, "middle[%d]: %d\n", i, maxColumn[i]);
			//	for (t = (int8_t*)&vMaxColumn, ti = 0; ti < 16; ++ti) fprintf(stderr, "%d\t", *t++);

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
		vTemp = simdui8_subs (vH, vGapO);
		vTemp = simdui8_subs (vF, vTemp);
		vTemp = simdi8_eq (vTemp, vZero);
		uint32_t cmp = simdi8_movemask (vTemp);
		while (cmp != SIMD_MOVEMASK_MAX) {
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

			vTemp = simdui8_subs (vH, vGapO);
			vTemp = simdui8_subs (vF, vTemp);
			vTemp = simdi8_eq (vTemp, vZero);
			cmp  = simdi8_movemask (vTemp);
		}

		vMaxScore = simdui8_max(vMaxScore, vMaxColumn);
		vTemp = simdi8_eq(vMaxMark, vMaxScore);
		cmp = simdi8_movemask(vTemp);
		if (cmp != SIMD_MOVEMASK_MAX) {
			uint8_t temp;
			vMaxMark = vMaxScore;
			max16(temp, vMaxScore);
			vMaxScore = vMaxMark;

			if (LIKELY(temp > max)) {
				max = temp;
				if (max + bias >= 255) break;	//overflow
				end_db = i;

				/* Store the column with the highest alignment score in order to trace the alignment ending position on read. */
				for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
			}
		}

		/* Record the max score of current column. */
		max16(maxColumn[i], vMaxColumn);
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
		//			fprintf (stderr, "maxColumn[%d]: %d\n", i, maxColumn[i]);
		if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
		}
	}
	edge = (end_db + maskLen) > db_length ? db_length : (end_db + maskLen);
	for (i = edge + 1; i < db_length; i ++) {
		//			fprintf (stderr, "db_length: %d\tmaxColumn[%d]: %d\n", db_length, i, maxColumn[i]);
		if (maxColumn[i] > best1.score) {
            best1.score = maxColumn[i];
            best1.ref = i;
		}
	}

	return std::make_pair(best0, best1);
#undef max16
}


std::pair<SmithWaterman::alignment_end, SmithWaterman::alignment_end> SmithWaterman::sw_sse2_word (const unsigned char* db_sequence,
														   int8_t ref_dir,	// 0: forward ref; 1: reverse ref
														   int32_t db_length,
														   int32_t query_length,
														   const uint8_t gap_open, /* will be used as - */
														   const uint8_t gap_extend, /* will be used as - */
														   const simd_int*query_profile_word,
														   uint16_t terminate,
														   int32_t maskLen) {

#define max8(m, vm) ((m) = simdi16_hmax((vm)));

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
	simd_int vTemp;
	int32_t edge, begin = 0, end = db_length, step = 1;

	/* outer loop to process the reference sequence */
	if (ref_dir == 1) {
		begin = db_length - 1;
		end = -1;
		step = -1;
	}
	for (i = begin; LIKELY(i != end); i += step) {
		simd_int e, vF = vZero; /* Initialize F value to 0.
                                Any errors to vH values will be corrected in the Lazy_F loop.
                                */
		simd_int vH = pvHStore[segLen - 1];
		vH = simdi8_shiftl (vH, 2); /* Shift the 128-bit value in vH left by 2 byte. */

		/* Swap the 2 H buffers. */
		simd_int* pv = pvHLoad;

		simd_int vMaxColumn = vZero; /* vMaxColumn is used to record the max values of column i. */

		const simd_int* vP = query_profile_word + db_sequence[i] * segLen; /* Right part of the query_profile_byte */
		pvHLoad = pvHStore;
		pvHStore = pv;

		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); j ++) {
			vH = simdi16_adds(vH, simdi_load(vP + j));

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
				if (UNLIKELY(! simdi8_movemask(simdi16_gt(vF, vH)))) goto end;
			}
		}

		end:
		vMaxScore = simdi16_max(vMaxScore, vMaxColumn);
		vTemp = simdi16_eq(vMaxMark, vMaxScore);
		uint32_t cmp = simdi8_movemask(vTemp);
		if (cmp != SIMD_MOVEMASK_MAX) {
			uint16_t temp;
			vMaxMark = vMaxScore;
			max8(temp, vMaxScore);
			vMaxScore = vMaxMark;

			if (LIKELY(temp > max)) {
				max = temp;
				end_ref = i;
				for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
			}
		}

		/* Record the max score of current column. */
		max8(maxColumn[i], vMaxColumn);
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
#undef max8
}

void SmithWaterman::ssw_init(const Sequence* q,
							 const int8_t* mat,
							 const BaseMatrix *m,
							 const int8_t score_size) {

	profile->bias = 0;
	profile->sequence_type = q->getSequenceType();
    const int32_t alphabetSize = m->alphabetSize;
	int32_t compositionBias = 0;
	bool isProfile = Parameters::isEqualDbtype(q->getSequenceType(), Parameters::DBTYPE_HMM_PROFILE)
	               || Parameters::isEqualDbtype(q->getSequenceType(), Parameters::DBTYPE_PROFILE_STATE_PROFILE);
	if (!isProfile && aaBiasCorrection) {
		SubstitutionMatrix::calcLocalAaBiasCorrection(m, q->numSequence, q->L, tmp_composition_bias);
		for (int i =0; i < q->L; i++) {
			profile->composition_bias[i] = (int8_t) (tmp_composition_bias[i] < 0.0)? tmp_composition_bias[i] - 0.5: tmp_composition_bias[i] + 0.5;
			compositionBias = (compositionBias < profile->composition_bias[i]) ? compositionBias : profile->composition_bias[i];
		}
		compositionBias = std::min(compositionBias, 0);
	} else {
		memset(profile->composition_bias, 0, q->L* sizeof(int8_t));
	}
	// copy memory to local memory
	if (Parameters::isEqualDbtype(profile->sequence_type, Parameters::DBTYPE_HMM_PROFILE)) {
		memcpy(profile->mat, mat, q->L * Sequence::PROFILE_AA_SIZE * sizeof(int8_t));
		// set neutral state 'X' (score=0)
		memset(profile->mat + ((alphabetSize - 1) * q->L), 0, q->L * sizeof(int8_t ));
	}else if(Parameters::isEqualDbtype(profile->sequence_type, Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
		memcpy(profile->mat, mat, q->L * alphabetSize * sizeof(int8_t));
	} else {
		memcpy(profile->mat, mat, alphabetSize * alphabetSize * sizeof(int8_t));
	}
	memcpy(profile->query_sequence, q->numSequence, q->L);
	if (score_size == 0 || score_size == 2) {
		/* Find the bias to use in the substitution matrix */
		int32_t bias = 0;
		int32_t matSize =  alphabetSize * alphabetSize;
		if (Parameters::isEqualDbtype(q->getSequenceType(), Parameters::DBTYPE_HMM_PROFILE)) {
			matSize = q->L * Sequence::PROFILE_AA_SIZE;
		}else if(Parameters::isEqualDbtype(q->getSequenceType(), Parameters::DBTYPE_PROFILE_STATE_PROFILE)) {
			matSize = q->L * alphabetSize;
		}

		for (int32_t i = 0; i < matSize; i++){
			if (mat[i] < bias){
				bias = mat[i];
			}
		}
		bias = abs(bias) + abs(compositionBias);
		profile->bias = bias;
		if (isProfile) {
			createQueryProfile<int8_t, VECSIZE_INT * 4, PROFILE>(profile->profile_byte, profile->query_sequence, NULL, profile->mat, q->L, alphabetSize, bias, 1, q->L);
		} else {
			createQueryProfile<int8_t, VECSIZE_INT * 4, SUBSTITUTIONMATRIX>(profile->profile_byte, profile->query_sequence, profile->composition_bias, profile->mat, q->L, alphabetSize, bias, 0, 0);
		}
	}
	if (score_size == 1 || score_size == 2) {
		if (isProfile) {
			createQueryProfile<int16_t, VECSIZE_INT * 2, PROFILE>(profile->profile_word, profile->query_sequence, NULL, profile->mat, q->L, alphabetSize, 0, 1, q->L);
			for (int32_t i = 0; i< alphabetSize; i++) {
				profile->profile_word_linear[i] = &profile_word_linear_data[i*q->L];
				for (int j = 0; j < q->L; j++) {
					//TODO is this right? :O
					profile->profile_word_linear[i][j] = mat[i * q->L + j];
				}
			}
		}else{
			createQueryProfile<int16_t, VECSIZE_INT * 2, SUBSTITUTIONMATRIX>(profile->profile_word, profile->query_sequence, profile->composition_bias, profile->mat, q->L, alphabetSize, 0, 0, 0);
			for(int32_t i = 0; i< alphabetSize; i++) {
				profile->profile_word_linear[i] = &profile_word_linear_data[i*q->L];
				for (int j = 0; j < q->L; j++) {
					profile->profile_word_linear[i][j] = mat[i * alphabetSize + q->numSequence[j]] + profile->composition_bias[j];
				}
			}
		}


	}
	// create reverse structures
	seq_reverse( profile->query_rev_sequence, profile->query_sequence, q->L);
	seq_reverse( profile->composition_bias_rev, profile->composition_bias, q->L);

	if (isProfile) {
		for (int32_t i = 0; i < alphabetSize; i++) {
			const int8_t *startToRead = profile->mat + (i * q->L);
			int8_t *startToWrite      = profile->mat_rev + (i * q->L);
			std::reverse_copy(startToRead, startToRead + q->L, startToWrite);
		}
	}
	profile->query_length = q->L;
	profile->alphabetSize = alphabetSize;
}
template <const unsigned int type>
SmithWaterman::cigar * SmithWaterman::banded_sw(const unsigned char *db_sequence, const int8_t *query_sequence, const int8_t * compositionBias,
												int32_t db_length, int32_t query_length, int32_t queryStart,
												int32_t score, const uint32_t gap_open,
												const uint32_t gap_extend, int32_t band_width, const int8_t *mat, int32_t n) {
	/*! @function
     @abstract  Round an integer to the next closest power-2 integer.
     @param  x  integer to be rounded (in place)
     @discussion x will be modified.
     */
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
		for (j = 1; LIKELY(j < width - 1); j ++) h_b[j] = 0;
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
				set_u(u, band_width, i, j);	set_u(e, band_width, i - 1, j);
				set_u(b, band_width, i, j - 1); set_u(d, band_width, i - 1, j - 1);
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

				e1 = e_b[u] > 0 ? e_b[u] : 0;
				f1 = f > 0 ? f : 0;
				temp1 = e1 > f1 ? e1 : f1;
				if(type == SUBSTITUTIONMATRIX){
					temp2 = h_b[d] + mat[query_sequence[i] * n + db_sequence[j]] + compositionBias[i];
				}
				if(type == PROFILE) {
					temp2 = h_b[d] + mat[db_sequence[j] * n + (queryStart + i)];
				}
				h_c[u] = temp1 > temp2 ? temp1 : temp2;

				if (h_c[u] > max) max = h_c[u];

				if (temp1 <= temp2) direction_line[dh] = 1;
				else direction_line[dh] = e1 > f1 ? direction_line[de] : direction_line[df];
			}
			for (j = 1; j <= u; j ++) h_b[j] = h_c[j];
		}
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

s_align SmithWaterman::scoreIdentical(unsigned char *dbSeq, int L, EvalueComputation * evaluer, int alignmentMode) {
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
	r.cigar = new uint32_t[L];
	short score = 0;
	for(int pos = 0; pos < L; pos++){
		int currScore = profile->profile_word_linear[dbSeq[pos]][pos];
		score += currScore;
		r.cigar[pos] = 'M';
	}
	r.score1=score;
	r.evalue = evaluer->computeEvalue(r.score1, profile->query_length);

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
	const int W = (profile->query_length + (element_count - 1)) / element_count; // width of bands in query and score matrix = hochgerundetes LQ/16

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

	memset(vHStore,0,W*sizeof(simd_int));
	memset(vHLoad,0,W*sizeof(simd_int));

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

