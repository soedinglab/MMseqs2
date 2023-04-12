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
#include <iostream>

SmithWaterman::SmithWaterman(size_t maxSequenceLength, int aaSize, bool aaBiasCorrection,
                             float aaBiasCorrectionScale, int targetSeqType) {
	maxSequenceLength += 1;
    this->aaBiasCorrectionScale = aaBiasCorrectionScale;
	this->aaBiasCorrection = aaBiasCorrection;

	int segmentSize = (maxSequenceLength+7)/8;
    segSize = segmentSize;
	vHStore = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	vHLoad  = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	vE      = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	vHmax   = (simd_int*) mem_align(ALIGN_INT, segSize * sizeof(simd_int));

	// setting up target
	target_profile_byte = (simd_int*) mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));


    isTargetProfile = Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_HMM_PROFILE);
    isQueryProfile = false;

    // setting up query
	profile = new s_profile();
	// query profile
	profile->profile_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
	profile->profile_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_rev_byte = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
    profile->profile_rev_word = (simd_int*)mem_align(ALIGN_INT, aaSize * segSize * sizeof(simd_int));
#ifdef GAP_POS_SCORING
    profile->profile_gDelOpen_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelOpen_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelClose_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelClose_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gIns_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gIns_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelOpen_rev_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelOpen_rev_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelClose_rev_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gDelClose_rev_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gIns_rev_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->profile_gIns_rev_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->gDelOpen = new uint8_t[maxSequenceLength];
    profile->gDelClose = new uint8_t[maxSequenceLength];
    profile->gDelOpen_rev = new uint8_t[maxSequenceLength];
    profile->gDelClose_rev = new uint8_t[maxSequenceLength];
    profile->gIns_rev = new uint8_t[maxSequenceLength];
#endif
    // query consensus profile
	profile->consens_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	profile->consens_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
	profile->consens_rev_byte = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    profile->consens_rev_word = (simd_int*)mem_align(ALIGN_INT, segSize * sizeof(simd_int));
    // query sequence
    profile->query_sequence     = new int8_t[maxSequenceLength];
	profile->query_rev_sequence = new int8_t[maxSequenceLength];
    profile->query_consens_sequence = new int8_t[maxSequenceLength];
    profile->query_consens_rev_sequence = new int8_t[maxSequenceLength];
	profile->composition_bias   = new int8_t[maxSequenceLength];
	profile->composition_bias_rev   = new int8_t[maxSequenceLength];
	profile->profile_word_linear = new short*[aaSize];
	profile_word_linear_data = new short[aaSize*maxSequenceLength];
	profile->mat_rev            = new int8_t[std::max(maxSequenceLength, (size_t)aaSize) * aaSize * 2];
	profile->mat                = new int8_t[std::max(maxSequenceLength, (size_t)aaSize) * aaSize * 2];
	tmp_composition_bias   = new float[maxSequenceLength];
    scorePerCol = new int8_t[maxSequenceLength];
    /* array to record the largest score of each reference position */
	maxColumn = new uint8_t[maxSequenceLength*sizeof(uint16_t)];
	memset(maxColumn, 0, maxSequenceLength*sizeof(uint16_t));
	memset(profile->query_sequence, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->query_rev_sequence, 0, maxSequenceLength * sizeof(int8_t));
    memset(profile->query_consens_sequence, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->query_consens_rev_sequence, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->mat_rev, 0, maxSequenceLength * aaSize);
	memset(profile->composition_bias, 0, maxSequenceLength * sizeof(int8_t));
	memset(profile->composition_bias_rev, 0, maxSequenceLength * sizeof(int8_t));
}

SmithWaterman::~SmithWaterman(){
	free(vHStore);
	free(vHLoad);
	free(vE);
	free(vHmax);
	free(target_profile_byte);
	free(profile->profile_byte);
	free(profile->profile_word);
	free(profile->profile_rev_byte);
	free(profile->profile_rev_word);
	free(profile->consens_byte);
	free(profile->consens_word);
	free(profile->consens_rev_byte);
	free(profile->consens_rev_word);
#ifdef GAP_POS_SCORING
    free(profile->profile_gDelOpen_byte);
    free(profile->profile_gDelOpen_word);
    free(profile->profile_gDelClose_byte);
    free(profile->profile_gDelClose_word);
    free(profile->profile_gIns_byte);
    free(profile->profile_gIns_word);
    free(profile->profile_gDelOpen_rev_byte);
    free(profile->profile_gDelOpen_rev_word);
    free(profile->profile_gDelClose_rev_byte);
    free(profile->profile_gDelClose_rev_word);
    free(profile->profile_gIns_rev_byte);
    free(profile->profile_gIns_rev_word);
    delete[] profile->gDelOpen;
    delete[] profile->gDelClose;
    delete[] profile->gDelOpen_rev;
    delete[] profile->gDelClose_rev;
    delete[] profile->gIns_rev;
#endif
	delete [] profile->query_rev_sequence;
	delete [] profile->query_sequence;
	delete [] profile->query_consens_sequence;
	delete [] profile->query_consens_rev_sequence;
	delete [] profile->composition_bias;
	delete [] profile->composition_bias_rev;
	delete [] profile->profile_word_linear;
	delete [] profile_word_linear_data;
	delete [] profile->mat_rev;
	delete [] profile->mat;
	delete [] tmp_composition_bias;
	delete [] scorePerCol;
	delete [] maxColumn;
	delete profile;
}


/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch. */
template <typename T, size_t Elements, const unsigned int type>
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
                    *t++ = ( j >= query_length) ? bias : mat[nt * aaSize + query_sequence[j + offset ]] + composition_bias[j + offset] + bias; // mat[nt][q[j]] mat eq 20*20
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
		}
	}
}

template <typename T, size_t Elements>
void SmithWaterman::createConsensProfile(simd_int *profile, const int8_t *consens_sequence, const int32_t query_length, const int32_t offset) {
    const int32_t segLen = (query_length + Elements - 1) / Elements;
    T* t = (T*) profile;
    for (int32_t i = 0; i < segLen; i++) {
        int32_t j = i;
        for (size_t segNum = 0; LIKELY(segNum < Elements); segNum++) {
            // beyond the length of query so pad with neutral consensus values
            *t++ = (j >= query_length) ? 20 :  consens_sequence[j + offset];
//            *t++ = (j >= query_length) ? 20 :  consens_sequence[j + (offset - 1)];
            j += segLen;
        }
    }

}


void SmithWaterman::createTargetProfile(simd_int *profile, const int8_t *mat, const int target_length,
        const int32_t aaSize, uint8_t bias) {
    const int32_t segSize = 32;
    int8_t* t = (int8_t*) profile;
    for (int i = 0; i < target_length; i++) {
        for (int j = 0; j < segSize; j++) {
            // beyond the length of amino acids so pad with neutral weights
            *t++ = (j >= aaSize) ? bias : mat[i + j * target_length] + bias;
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

#ifdef GAP_POS_SCORING
template <typename T, size_t Elements>
void createGapProfile(simd_int* profile_gDelOpen, simd_int* profile_gDelClose, simd_int* profile_gIns,
                      const uint8_t* gDelOpen, const uint8_t* gDelClose, const uint8_t* gIns,
                      const int32_t query_length, const int32_t offset) {
    const int32_t segLen = (query_length - offset + Elements - 1) / Elements;
    T* delOpen = (T*) profile_gDelOpen;
    T* delClose = (T*) profile_gDelClose;
    T* ins = (T*) profile_gIns;
    for (int32_t i = 0; LIKELY(i < segLen); ++i) {
        int32_t j = i;
        for (size_t segNum = 0; LIKELY(segNum < Elements); ++segNum) {
            *delOpen++ = (j < query_length) ? gDelOpen[j + offset + 1] : 0; // offset + 1 because it calculates F for the next column
            *delClose++ = (j < query_length) ? gDelClose[j + offset + 1] : 0;
            *ins++ = (j < query_length) ? gIns[j + offset] : 0;
            j += segLen;
        }
    }
}
#endif

s_align SmithWaterman::ssw_align (
        const unsigned char *db_num_sequence,
        const unsigned char *db_consens_sequence,
        const int8_t *db_mat,
        int32_t db_length,
        std::string & backtrace,
        const uint8_t gap_open,
        const uint8_t gap_extend,
        const uint8_t alignmentMode,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
        const double  evalueThr,
        EvalueComputation * evaluer,
        const int covMode, const float covThr, const float correlationScoreWeight,
        const int32_t maskLen, const size_t id) {
    s_align alignment;
    // check if both query and target are profiles
    if (isQueryProfile && isTargetProfile) {
        alignment = ssw_align_private<SmithWaterman::PROFILE_PROFILE, true>(db_consens_sequence, db_mat, db_length, backtrace, gap_open,
                                                                       gap_extend, alignmentMode, evalueThr, evaluer, covMode, covThr, correlationScoreWeight, maskLen, id);
    } else if (isQueryProfile && !isTargetProfile) {
        alignment = ssw_align_private<SmithWaterman::PROFILE_SEQ, true>(db_num_sequence, db_mat, db_length, backtrace, gap_open,
                                                                  gap_extend, alignmentMode, evalueThr, evaluer, covMode, covThr, correlationScoreWeight, maskLen, id);
    } else if (!isQueryProfile && isTargetProfile) {
        alignment = ssw_align_private<SmithWaterman::SEQ_PROFILE, false>(db_num_sequence, db_mat, db_length, backtrace, gap_open,
                                                                  gap_extend, alignmentMode, evalueThr, evaluer, covMode, covThr, correlationScoreWeight, maskLen, id);
    } else {
        alignment = ssw_align_private<SmithWaterman::SEQ_SEQ, false>(db_num_sequence, db_mat, db_length, backtrace, gap_open,
                                                              gap_extend, alignmentMode, evalueThr, evaluer, covMode, covThr, correlationScoreWeight, maskLen, id);
    }
    return alignment;
}

template <const unsigned int type, const bool posSpecificGaps>
s_align SmithWaterman::ssw_align_private (
		const unsigned char *db_sequence,
		const int8_t *db_mat,
		int32_t db_length,
        std::string & backtrace,
        const uint8_t gap_open,
		const uint8_t gap_extend,
		const uint8_t alignmentMode,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
		const double  evalueThr,
		EvalueComputation * evaluer,
		const int covMode, const float covThr, const float correlationScoreWeight,
		const int32_t maskLen, const size_t id) {

    target_id = id;
	int32_t word = 0, query_length = profile->query_length;
	int32_t band_width = 0;
	cigar* path;
	s_align r;
	r.dbStartPos1 = -1;
	r.qStartPos1 = -1;
	r.cigar = 0;
	r.cigarLen = 0;

    std::pair<alignment_end, alignment_end> bests;
    std::pair<alignment_end, alignment_end> bests_reverse;

    simd_int* db_profile_byte = target_profile_byte;

    const int32_t qry_n = profile->query_length;
    const int32_t db_n = db_length;
    const unsigned char * db_consens_seq = db_sequence;
    const int8_t *db_matrix = db_mat;

    // find the alignment position
    if (profile->profile_byte) {
        if (type == PROFILE_PROFILE) {
            uint8_t db_bias = computeBias(db_length, db_mat, profile->alphabetSize);
            if (db_bias > profile->bias) {
                uint8_t shift = abs(profile->bias - db_bias);
                updateQueryProfile<int8_t, VECSIZE_INT * 4>(profile->profile_byte, profile->query_length, profile->alphabetSize, shift);
            }
            profile->bias = std::max(db_bias, profile->bias);
            createTargetProfile(db_profile_byte, db_mat, db_length, profile->alphabetSize - 1, profile->bias);
        }
        bests = sw_sse2_byte<type,posSpecificGaps>(db_sequence, db_profile_byte, 0, db_length, query_length, gap_open, gap_extend,
                profile->profile_byte, profile->consens_byte,
#ifdef GAP_POS_SCORING
				profile->profile_gDelOpen_byte, profile->profile_gDelClose_byte, profile->profile_gIns_byte,
#endif
				UCHAR_MAX, profile->bias, maskLen);
        if (bests.first.score == 255) {
            bests = sw_sse2_word<type,posSpecificGaps>(db_sequence, db_profile_byte, 0, db_length, query_length, gap_open, gap_extend,
                    profile->profile_word, profile->consens_word,
#ifdef GAP_POS_SCORING
					profile->profile_gDelOpen_word, profile->profile_gDelClose_word, profile->profile_gIns_word,
#endif
					USHRT_MAX, profile->bias, maskLen);
            word = 1;
        }
    } else {
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

    // no residue could be aligned
    if (r.dbEndPos1 == -1) {
        return r;
    }
    int32_t queryOffset = query_length - r.qEndPos1 -1;
	r.evalue = evaluer->computeEvalue(r.score1, query_length);
    bool hasLowerEvalue = r.evalue > evalueThr;
	r.qCov = computeCov(0, r.qEndPos1, query_length);
	r.tCov = computeCov(0, r.dbEndPos1, db_length);
    bool hasLowerCoverage = !(Util::hasCoverage(covThr, covMode, r.qCov, r.tCov));

	if (alignmentMode == 0 || ((alignmentMode == 2 || alignmentMode == 1) && (hasLowerEvalue || hasLowerCoverage))) {
        return r;
	}

	if (word == 0) {
	    if (type == PROFILE_SEQ || type == PROFILE_PROFILE) {
	        createQueryProfile<int8_t, VECSIZE_INT * 4, PROFILE>(profile->profile_rev_byte, profile->query_rev_sequence, NULL, profile->mat_rev,
                                                                     r.qEndPos1 + 1, profile->alphabetSize, profile->bias, queryOffset, profile->query_length);
#ifdef GAP_POS_SCORING
			if (posSpecificGaps) {
				createGapProfile<int8_t, VECSIZE_INT * 4>(profile->profile_gDelOpen_rev_byte, profile->profile_gDelClose_rev_byte, profile->profile_gIns_rev_byte,
														profile->gDelOpen_rev, profile->gDelClose_rev, profile->gIns_rev, profile->query_length, queryOffset);
			}
#endif
	        if (type == PROFILE_PROFILE) {
                createConsensProfile<int8_t, VECSIZE_INT * 4>(profile->consens_rev_byte, profile->query_consens_rev_sequence,
                                                              r.qEndPos1 + 1, queryOffset);
	        }
	    } else {
            createQueryProfile<int8_t, VECSIZE_INT * 4, SUBSTITUTIONMATRIX>(profile->profile_rev_byte,
                                                                            profile->query_rev_sequence,
                                                                            profile->composition_bias_rev,
                                                                            profile->mat,r.qEndPos1 + 1,
                                                                            profile->alphabetSize, profile->bias,
                                                                            queryOffset, 0);
	    }
        bests_reverse = sw_sse2_byte<type,posSpecificGaps>(db_sequence, db_profile_byte, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1, gap_open,
                                           gap_extend, profile->profile_rev_byte, profile->consens_rev_byte,
#ifdef GAP_POS_SCORING
                                           profile->profile_gDelOpen_rev_byte, profile->profile_gDelClose_rev_byte, profile->profile_gIns_rev_byte,
#endif
										   r.score1, profile->bias, maskLen);
	} else {
        if (type == PROFILE_SEQ || type == PROFILE_PROFILE) {
            createQueryProfile<int16_t, VECSIZE_INT * 2, PROFILE>(profile->profile_rev_word,
                                                                  profile->query_rev_sequence, NULL, profile->mat_rev,
                                                                  r.qEndPos1 + 1, profile->alphabetSize, 0, queryOffset,
                                                                  profile->query_length);
#ifdef GAP_POS_SCORING
			if (posSpecificGaps) {
				createGapProfile<int16_t, VECSIZE_INT * 2>(profile->profile_gDelOpen_rev_word,
														profile->profile_gDelClose_rev_word,
														profile->profile_gIns_rev_word, profile->gDelOpen_rev,
														profile->gDelClose_rev, profile->gIns_rev,
														profile->query_length, queryOffset);
			}
#endif
            if (type == PROFILE_PROFILE) {
                createConsensProfile<int16_t, VECSIZE_INT * 2>(profile->consens_rev_word,
                                                               profile->query_consens_rev_sequence,
                                                               r.qEndPos1 + 1, queryOffset);
            }
        } else {
            createQueryProfile<int16_t, VECSIZE_INT * 2, SUBSTITUTIONMATRIX>(profile->profile_rev_word,
                                                                             profile->query_rev_sequence,
                                                                             profile->composition_bias_rev,
                                                                             profile->mat,
                                                                             r.qEndPos1 + 1, profile->alphabetSize, 0,
                                                                             queryOffset, 0);
        }
        bests_reverse = sw_sse2_word<type,posSpecificGaps>(db_sequence, db_profile_byte, 1, r.dbEndPos1 + 1, r.qEndPos1 + 1, gap_open,
                                           gap_extend,
                                           profile->profile_rev_word, profile->consens_rev_word,
#ifdef GAP_POS_SCORING
                                           profile->profile_gDelOpen_rev_word, profile->profile_gDelClose_rev_word, profile->profile_gIns_rev_word,
#endif
										   r.score1, profile->bias, maskLen);
    }


	if(bests_reverse.first.score != r.score1){
        fprintf(stderr, "Score of forward/backward SW differ: %d %d. Q: %lu T: %lu.\n", r.score1, bests_reverse.first.score, query_id, target_id);
        fprintf(stderr, "Start: Q: %d, T: %d. End: Q: %d, T %d\n", r.qEndPos1 - bests_reverse.first.read, bests_reverse.first.ref, r.qEndPos1, r.dbEndPos1);
        //  if qry is not a profile, just exit
        if (!(type == PROFILE_SEQ) || !(type == PROFILE_PROFILE)) {
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
	hasLowerCoverage = !(Util::hasCoverage(covThr, covMode, r.qCov, r.tCov));
    // only start and end point are needed
    if (alignmentMode == 1 || hasLowerCoverage) {
        return r;
    }

	// Generate cigar.
// db_length and query_length updated
    db_length = r.dbEndPos1 - r.dbStartPos1 + 1;
    query_length = r.qEndPos1 - r.qStartPos1 + 1;
    band_width = abs(db_length - query_length) + 1;


	// TODO: fix banded_sw
	if (type == PROFILE_PROFILE) {
	    path = banded_sw<type,posSpecificGaps>(db_consens_seq + r.dbStartPos1, profile->query_sequence + r.qStartPos1, profile->query_consens_sequence + r.qStartPos1, NULL, db_length,
	            query_length, r.qStartPos1, r.dbStartPos1, r.score1, gap_open, gap_extend,
#ifdef GAP_POS_SCORING
				profile->gDelOpen + r.qStartPos1, profile->gDelClose + r.qStartPos1, profile->gIns + r.qStartPos1,
#endif
				band_width, profile->mat, db_matrix, qry_n, db_n);
	} else if (type == PROFILE_SEQ) {
        path = banded_sw<type,posSpecificGaps>(db_sequence + r.dbStartPos1, profile->query_sequence + r.qStartPos1, NULL, profile->composition_bias + r.qStartPos1, db_length,
                query_length, r.qStartPos1, r.dbStartPos1, r.score1, gap_open, gap_extend,
#ifdef GAP_POS_SCORING
				profile->gDelOpen + r.qStartPos1, profile->gDelClose + r.qStartPos1, profile->gIns + r.qStartPos1,
#endif
				band_width, profile->mat, NULL, profile->query_length, 0);
	} else {
        path = banded_sw<type,posSpecificGaps>(db_sequence + r.dbStartPos1, profile->query_sequence + r.qStartPos1, NULL, profile->composition_bias + r.qStartPos1, db_length,
                query_length, r.qStartPos1, r.dbStartPos1, r.score1, gap_open, gap_extend,
#ifdef GAP_POS_SCORING
				nullptr, nullptr, nullptr,
#endif
		band_width, profile->mat, NULL, profile->alphabetSize, 0);
		db_length = r.dbEndPos1 - r.dbStartPos1 + 1;
		query_length = r.qEndPos1 - r.qStartPos1 + 1;
		band_width = abs(db_length - query_length) + 1;
	}

    if (path != NULL) {
        r.cigar = path->seq;
        r.cigarLen = path->length;
    }
    delete path;

    uint32_t aaIds = 0;
    size_t mStateCnt = 0;
    if (type == PROFILE_PROFILE) {
        computerBacktrace<PROFILE>(profile, db_sequence, r, backtrace, aaIds, scorePerCol, mStateCnt);
    }else{
        computerBacktrace<SUBSTITUTIONMATRIX>(profile, db_sequence, r, backtrace, aaIds, scorePerCol, mStateCnt);
    }
    r.identicalAACnt = aaIds;
    if(correlationScoreWeight > 0.0){
        int correlationScore = computeCorrelationScore(scorePerCol, mStateCnt);
        r.score1 += static_cast<float>(correlationScore) * correlationScoreWeight;
        r.evalue = evaluer->computeEvalue(r.score1, query_length);
    }
    return r;
}

template <const unsigned int type>
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
                if(type == PROFILE){
                    scorePerCol[mStatesCnt] = query->mat[db_sequence[targetPos] * query->query_length + queryPos];
                }
                if(type == SUBSTITUTIONMATRIX){
                    scorePerCol[mStatesCnt] = query->mat[query->query_sequence[queryPos] * query->alphabetSize + db_sequence[targetPos]] + query->composition_bias[queryPos];
                }
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


template <const unsigned int type, const bool posSpecificGaps>
std::pair<SmithWaterman::alignment_end, SmithWaterman::alignment_end> SmithWaterman::sw_sse2_byte (
														   const unsigned char *db_sequence,
														   const simd_int* db_profile_byte,
														   int8_t ref_dir,	// 0: forward ref; 1: reverse ref
														   int32_t db_length,
														   int32_t query_length,
														   const uint8_t gap_open, /* will be used as - */
														   const uint8_t gap_extend, /* will be used as - */
														   const simd_int* query_profile_byte, /* profile_byte loaded in ssw_init */
														   const simd_int* query_consens_byte, /* profile_consens_byte loaded in ssw_init */
#ifdef GAP_POS_SCORING
                                                           const simd_int *gap_open_del,
                                                           const simd_int *gap_close_del,
                                                           const simd_int *gap_open_ins,
#endif
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

    //fprintf(stderr, "start alignment of length %d [%u]\n", query_length, segLen * SIMD_SIZE);

	/* outer loop to process the reference sequence */
	if (ref_dir == 1) {
		begin = db_length - 1;
		end = -1;
		step = -1;
	}

#ifndef AVX2
    const simd_int sixten  = simdi8_set(16);
    const simd_int fiveten = simdi8_set(15);
#endif

    // store the query consensus profile
    const simd_int* vQueryCons = query_consens_byte;
	for (i = begin; LIKELY(i != end); i += step) {
//	    cnt = i;
		simd_int e, vF = vZero, vMaxColumn = vZero; /* Initialize F value to 0.
                                                    Any errors to vH values will be corrected in the Lazy_F loop.
                                                    */

		simd_int vH = pvHStore[segLen - 1];
		vH = simdi8_shiftl (vH, 1); /* Shift the 128-bit value in vH left by 1 byte. */
		const simd_int* vP = query_profile_byte + db_sequence[i] * segLen; /* Right part of the query_profile_byte */

#ifdef AVX2
        simd_int target_scores1 = simdi8_set(0);
        if (type == PROFILE_PROFILE) {
            target_scores1 = simdi_load(&db_profile_byte[i]);
        }
#else
        simd_int target_scores1 = simdi8_set(0);
        simd_int target_scores2 = simdi8_set(0);
        if (type == PROFILE_PROFILE) {
            target_scores1 = simdi_load(&db_profile_byte[i]);
            target_scores2 = simdi_load(&db_profile_byte[i + 16]);
        }
#endif

		/* Swap the 2 H buffers. */
		simd_int* pv = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pv;

		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); ++j) {
		    simd_int score = simdi8_set(0);
		    if (type == PROFILE_PROFILE) {
#ifdef AVX2
                __m256i scoreLookup = UngappedAlignment::Shuffle(target_scores1, simdi_load(vQueryCons + j));
#else
		        const __m128i vQueryConsJ = _mm_load_si128(vQueryCons + j);
                __m128i score01 = _mm_shuffle_epi8(target_scores1, vQueryConsJ);
                __m128i score16 = _mm_shuffle_epi8(target_scores2, vQueryConsJ);
                __m128i lookup_mask01 = _mm_cmplt_epi8(vQueryConsJ, sixten);
                __m128i lookup_mask16 = _mm_cmplt_epi8(fiveten, vQueryConsJ);
                score01 = _mm_and_si128(lookup_mask01, score01);
                score16 = _mm_and_si128(lookup_mask16, score16);
                __m128i scoreLookup = _mm_add_epi8(score01, score16);
#endif
                //score = simdui8_max(scoreLookup, simdi_load(vP + j));
                score = simdui8_avg(scoreLookup, simdi_load(vP + j));
		    } else {
		        score = simdi_load(vP + j);
		    }


            vH = simdui8_adds(vH, score);
            vH = simdui8_subs(vH, vBias);   /* vH will be always > 0 */

            /* Get max from vH, vE and vF. */
			e = simdi_load(pvE + j);
            vH = simdui8_max(vH, e);
#ifdef GAP_POS_SCORING
            if (posSpecificGaps) {
                vH = simdui8_max(vH, simdui8_subs(vF, simdi_load(gap_close_del + j)));
            } else {
#endif
                vH = simdui8_max(vH, vF);
#ifdef GAP_POS_SCORING
            }
#endif
			vMaxColumn = simdui8_max(vMaxColumn, vH);

			/* Save vH values. */
			simdi_store(pvHStore + j, vH);

			/* Update vE value. */
#ifdef GAP_POS_SCORING
            if (posSpecificGaps) {
                // copy vH for update of vF
                vTemp = vH;
                vH = simdui8_subs(vH, simdi_load(gap_open_ins + j)); /* saturation arithmetic, result >= 0 */
            } else {
#endif
                vH = simdui8_subs(vH, vGapO); /* saturation arithmetic, result >= 0 */
#ifdef GAP_POS_SCORING
            }
#endif

			e = simdui8_subs(e, vGapE);
			e = simdui8_max(e, vH);
			simdi_store(pvE + j, e);

			/* Update vF value. */
			vF = simdui8_subs(vF, vGapE);
#ifdef GAP_POS_SCORING
            if (posSpecificGaps) {
                vF = simdui8_max(vF, simdui8_subs(vTemp, simdi_load(gap_open_del + j)));
            } else {
#endif
                vF = simdui8_max(vF, vH);
#ifdef GAP_POS_SCORING
            }
#endif

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
#ifdef GAP_POS_SCORING
        if (posSpecificGaps) {
            vTemp = simdui8_subs(vH, simdi_load(gap_open_del + j));
        } else {
#endif
            vTemp = simdui8_subs(vH, vGapO);
#ifdef GAP_POS_SCORING
        }
#endif
		vTemp = simdui8_subs (vF, vTemp);
		vTemp = simdi8_eq (vTemp, vZero);
		uint32_t cmp = simdi8_movemask (vTemp);
		while (cmp != SIMD_MOVEMASK_MAX) {
#ifdef GAP_POS_SCORING
            if (posSpecificGaps) {
                vH = simdui8_max (vH, simdui8_subs(vF, simdi_load(gap_close_del + j)));
                simdi_store(pvE + j, simdui8_max(simdi_load(pvE + j), simdui8_subs(vH, simdi_load(gap_open_ins + j))));
            } else {
#endif
                vH = simdui8_max (vH, vF);
#ifdef GAP_POS_SCORING
            }
#endif

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

#ifdef GAP_POS_SCORING
            if (posSpecificGaps) {
                vTemp = simdui8_subs(vH, simdi_load(gap_open_del + j));
            } else {
#endif
                vTemp = simdui8_subs(vH, vGapO);
#ifdef GAP_POS_SCORING
            }
#endif
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
#undef max16
}

template <const unsigned int type, const bool posSpecificGaps>
std::pair<SmithWaterman::alignment_end, SmithWaterman::alignment_end> SmithWaterman::sw_sse2_word (const unsigned char* db_sequence,
														   const simd_int* db_profile_byte,
														   int8_t ref_dir,	// 0: forward ref; 1: reverse ref
														   int32_t db_length,
														   int32_t query_length,
														   const uint8_t gap_open, /* will be used as - */
														   const uint8_t gap_extend, /* will be used as - */
														   const simd_int* query_profile_word,
														   const simd_int* query_consens_word,
#ifdef GAP_POS_SCORING
                                                           const simd_int *gap_open_del,
                                                           const simd_int *gap_close_del,
                                                           const simd_int *gap_open_ins,
#endif
														   uint16_t terminate,
                                                           const uint16_t bias,
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

	/* 16 byte bias vector */
	simd_int vBias = simdi16_set(bias);
	//simd_int vBias = simdi16_set(-bias);    // set as a negative value for simd use
	simd_int vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
	simd_int vMaxMark = vZero; /* Trace the highest score till the previous column. */
	simd_int vTemp;
	int32_t edge, begin = 0, end = db_length, step = 1;

    //fprintf(stderr, "start alignment of length %d [%d]\n", query_length, segLen * SIMD_SIZE);

	/* outer loop to process the reference sequence */
	if (ref_dir == 1) {
		begin = db_length - 1;
		end = -1;
		step = -1;
	}

    // store the query consensus profile
	const simd_int* vQueryCons = query_consens_word;

#ifndef AVX2
    const simd_int sixten  = simdi8_set(16);
    const simd_int fiveten = simdi8_set(15);
#endif

	for (i = begin; LIKELY(i != end); i += step) {
		simd_int e, vF = vZero, vMaxColumn = vZero; /* Initialize F value to 0.
                                Any errors to vH values will be corrected in the Lazy_F loop.
                                */

		simd_int vH = pvHStore[segLen - 1];
		vH = simdi8_shiftl (vH, 2); /* Shift the 128-bit value in vH left by 2 byte. */
		const simd_int* vP = query_profile_word + db_sequence[i] * segLen; /* Right part of the query_profile_byte */

#ifdef AVX2
        simd_int target_scores1 = simdi16_set(0);
        if (type == PROFILE_PROFILE) {
            target_scores1 = simdi_load(&db_profile_byte[i]);
        }
#else
        simd_int target_scores1 = simdi16_set(0);
        simd_int target_scores2 = simdi16_set(0);
        if (type == PROFILE_PROFILE) {
            target_scores1 = simdi_load(&db_profile_byte[i]);
            target_scores2 = simdi_load(&db_profile_byte[i + 16]);
        }
#endif

        /* Swap the 2 H buffers. */
        simd_int* pv = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pv;

		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); j ++) {
		    simd_int score = simdi16_set(0);
		    if (type == PROFILE_PROFILE) {
#ifdef AVX2
                __m256i scoreLookup = UngappedAlignment::Shuffle(target_scores1, simdi_load(vQueryCons + j));
#else
                const __m128i vQueryConsJ = _mm_load_si128(vQueryCons + j);
                __m128i score01 = _mm_shuffle_epi8(target_scores1, vQueryConsJ);
                __m128i score16 = _mm_shuffle_epi8(target_scores2, vQueryConsJ);
                __m128i lookup_mask01 = _mm_cmplt_epi8(vQueryConsJ, sixten);
                __m128i lookup_mask16 = _mm_cmplt_epi8(fiveten, vQueryConsJ);
                score01 = _mm_and_si128(lookup_mask01, score01);
                score16 = _mm_and_si128(lookup_mask16, score16);
                __m128i scoreLookup = _mm_add_epi8(score01, score16);
#endif
                scoreLookup = simdi_and(scoreLookup, simdi16_set(0x00FF));
                score = simdui16_avg(scoreLookup, simdi16_add(simdi_load(vP + j), vBias));
				score = simdi16_sub(score, vBias);
				//scoreLookup = simdi16_add(scoreLookup, vBias);
                //score = simdi16_max(scoreLookup, simdi_load(vP + j));
		    } else {
		        score = simdi_load(vP + j);
		    }

			vH = simdi16_adds(vH, score);

			/* Get max from vH, vE and vF. */
			e = simdi_load(pvE + j);
            vH = simdi16_max(vH, e);
#ifdef GAP_POS_SCORING
            if (posSpecificGaps) {
                vH = simdi16_max(vH, simdui16_subs(vF, simdi_load(gap_close_del + j)));
            } else {
#endif
                vH = simdi16_max(vH, vF);
#ifdef GAP_POS_SCORING
            }
#endif

			vMaxColumn = simdi16_max(vMaxColumn, vH);

			/* Save vH values. */
			simdi_store(pvHStore + j, vH);

			/* Update vE value. */
#ifdef GAP_POS_SCORING
            if (posSpecificGaps) {
                // copy vH for update of vF
                vTemp = vH;
                vH = simdui16_subs(vH, simdi_load(gap_open_ins + j)); /* saturation arithmetic, result >= 0 */
            } else {
#endif
                vH = simdui16_subs(vH, vGapO); /* saturation arithmetic, result >= 0 */
#ifdef GAP_POS_SCORING
            }
#endif
			e = simdui16_subs(e, vGapE);
			e = simdi16_max(e, vH);
			simdi_store(pvE + j, e);

			/* Update vF value. */
			vF = simdui16_subs(vF, vGapE);
#ifdef GAP_POS_SCORING
            if (posSpecificGaps) {
                vF = simdi16_max(vF, simdui16_subs(vTemp, simdi_load(gap_open_del + j)));
            } else {
#endif
                vF = simdi16_max(vF, vH);
#ifdef GAP_POS_SCORING
            }
#endif

			/* Load the next vH. */
			vH = simdi_load(pvHLoad + j);
		}

		/* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
		for (k = 0; LIKELY(k < (int32_t) SIMD_SIZE); ++k) {
			vF = simdi8_shiftl (vF, 2);
			for (j = 0; LIKELY(j < segLen); ++j) {
				vH = simdi_load(pvHStore + j);
#ifdef GAP_POS_SCORING
                if (posSpecificGaps) {
                    vH = simdi16_max(vH, simdui16_subs(vF, simdi_load(gap_close_del + j)));
                    simdi_store(pvE + j, simdi16_max(simdi_load(pvE + j), simdui16_subs(vH, simdi_load(gap_open_ins + j))));
                } else {
#endif
                    vH = simdi16_max(vH, vF);
#ifdef GAP_POS_SCORING
                }
#endif

				vMaxColumn = simdi16_max(vMaxColumn, vH); //newly added line
				simdi_store(pvHStore + j, vH);
#ifdef GAP_POS_SCORING
                if (posSpecificGaps) {
                    vH = simdui16_subs(vH, simdi_load(gap_open_del + j));
                } else {
#endif
                    vH = simdui16_subs(vH, vGapO);
#ifdef GAP_POS_SCORING
                }
#endif
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

        //uint16_t *t = (uint16_t *)pvHStore;
        //for (size_t ti = 0; ti < segLen * SIMD_SIZE; ++ti) {
        //    fprintf(stderr, "%d ", t[ti / segLen + ti % segLen * SIMD_SIZE]);
        //}
        //fprintf(stderr, "\n");

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
							 const BaseMatrix *m) {

    query_id = q->getId();
	profile->bias = 0;
    profile->query_length = q->L;
	profile->sequence_type = q->getSequenceType();
	isQueryProfile = (Parameters::isEqualDbtype(profile->sequence_type, Parameters::DBTYPE_HMM_PROFILE));
    const int32_t alphabetSize = m->alphabetSize;
    profile->alphabetSize = m->alphabetSize;

	int32_t compositionBias = 0;
	bool isProfile = Parameters::isEqualDbtype(q->getSequenceType(), Parameters::DBTYPE_HMM_PROFILE);
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
	if (Parameters::isEqualDbtype(profile->sequence_type, Parameters::DBTYPE_HMM_PROFILE)) {
		memcpy(profile->mat, mat, q->L * Sequence::PROFILE_AA_SIZE * sizeof(int8_t));
		// set neutral state 'X' (score=0)
		memset(profile->mat + ((alphabetSize - 1) * q->L), 0, q->L * sizeof(int8_t ));
	} else {
		memcpy(profile->mat, mat, alphabetSize * alphabetSize * sizeof(int8_t));
	}
	memcpy(profile->query_sequence, q->numSequence, q->L);
	// numConsensusSequence points to NULL if not profile
	if (isQueryProfile) {
	    memcpy(profile->query_consens_sequence, q->numConsensusSequence, q->L);
	}
	// create gap-penalties profile
    if (isProfile) {
#ifdef GAP_POS_SCORING
        profile->gIns = q->gIns;
        // insertion penalties are shifted by one position for the reverse direction (2nd to last becomes first)
        std::reverse_copy(q->gIns, q->gIns + q->L - 1, profile->gIns_rev);

        for (int32_t i = 0; i < q->L; ++i) {
            profile->gDelOpen[i] = q->gDel[i] & 0xF;
            profile->gDelClose[i] = q->gDel[i] >> 4;
        }
        profile->gDelClose_rev[0] = 0;
        profile->gDelOpen_rev[0] = 0;
        std::reverse_copy(profile->gDelOpen + 1, profile->gDelOpen + q->L, profile->gDelClose_rev + 1);
        std::reverse_copy(profile->gDelClose + 1, profile->gDelClose + q->L, profile->gDelOpen_rev + 1);
#endif
        for (int32_t i = 0; i < alphabetSize; i++) {
            const int8_t *startToRead = profile->mat + (i * q->L);
            int8_t *startToWrite      = profile->mat_rev + (i * q->L);
            std::reverse_copy(startToRead, startToRead + q->L, startToWrite);
        }
    }

	int32_t bias = 0;
	int32_t matSize = alphabetSize * alphabetSize;
    if (Parameters::isEqualDbtype(q->getSequenceType(), Parameters::DBTYPE_HMM_PROFILE)) {
        matSize = q->L * Sequence::PROFILE_AA_SIZE;
    }
    for (int32_t i = 0; i < matSize; i++) {
        if (mat[i] < bias) {
            bias = mat[i];
        }
    }
    bias = abs(bias) + abs(compositionBias);
    profile->bias = bias;

    if (isProfile) {
        // offset = 1 when createQueryProfile + createConsensProfile is old versions
        // create byte version of profiles
        createQueryProfile<int8_t, VECSIZE_INT * 4, PROFILE>(profile->profile_byte, profile->query_sequence, NULL,
                                                             profile->mat, q->L, alphabetSize, profile->bias, 0, q->L);
        createConsensProfile<int8_t, VECSIZE_INT * 4>(profile->consens_byte, profile->query_consens_sequence, q->L, 0);
#ifdef GAP_POS_SCORING
        createGapProfile<int8_t, VECSIZE_INT * 4>(profile->profile_gDelOpen_byte, profile->profile_gDelClose_byte,
                                                  profile->profile_gIns_byte, profile->gDelOpen, profile->gDelClose, q->gIns, q->L, 0);
#endif
        // create word version of profiles
        createQueryProfile<int16_t, VECSIZE_INT * 2, PROFILE>(profile->profile_word, profile->query_sequence, NULL,
                                                              profile->mat, q->L, alphabetSize, 0, 0, q->L);
        createConsensProfile<int16_t, VECSIZE_INT * 2>(profile->consens_word, profile->query_consens_sequence, q->L, 0);
#ifdef GAP_POS_SCORING
        createGapProfile<int16_t, VECSIZE_INT * 2>(profile->profile_gDelOpen_word, profile->profile_gDelClose_word, profile->profile_gIns_word,
                                                   profile->gDelOpen, profile->gDelClose, q->gIns, q->L, 0);
#endif
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
        // create linear version of word profile
        for (int32_t i = 0; i< alphabetSize; i++) {
            profile->profile_word_linear[i] = &profile_word_linear_data[i*q->L];
            for (int j = 0; j < q->L; j++) {
                profile->profile_word_linear[i][j] = mat[i * alphabetSize + q->numSequence[j]] + profile->composition_bias[j];
            }
        }
    }

	// create reverse structures
	if (isProfile) {
        std::reverse_copy(profile->query_sequence, profile->query_sequence + q->L, profile->query_rev_sequence);
        std::reverse_copy(profile->query_consens_sequence, profile->query_consens_sequence + q->L, profile->query_consens_rev_sequence);
	} else {
        std::reverse_copy(profile->query_sequence, profile->query_sequence + q->L, profile->query_rev_sequence);
        std::reverse_copy(profile->composition_bias, profile->composition_bias + q->L, profile->composition_bias_rev);
	}


	if (isProfile) {
        for (int32_t i = 0; i < alphabetSize; i++) {
            const int8_t *startToRead = profile->mat + (i * q->L);
            int8_t *startToWrite = profile->mat_rev + (i * q->L);
            std::reverse_copy(startToRead, startToRead + q->L, startToWrite);
        }
    }
}

template <const unsigned int type, const bool posSpecificGaps>
SmithWaterman::cigar * SmithWaterman::banded_sw(const unsigned char *db_sequence, const int8_t *query_sequence,
                                                const int8_t *query_consens_sequence, const int8_t * compositionBias,
												int32_t db_length, int32_t query_length, int32_t queryStart,
												int32_t targetStart, int32_t score, const uint32_t gap_open,
												const uint32_t gap_extend,
#ifdef GAP_POS_SCORING
												uint8_t *gDelOpen, uint8_t *gDelClose, uint8_t *gIns,
#endif
												int32_t band_width, const int8_t *mat, const int8_t *db_mat,
												const int32_t qry_n, const int32_t tgt_n) {
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
    int32_t subScore = 0;

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

#ifdef GAP_POS_SCORING
                if (posSpecificGaps) {
                    temp1 = i == 0 ? -gap_open : h_b[e] - gDelOpen[i];
                } else {
#endif
                    temp1 = i == 0 ? -gap_open : h_b[e] - gap_open;
#ifdef GAP_POS_SCORING
                }
#endif

                temp2 = i == 0 ? -gap_extend : e_b[e] - gap_extend;
				e_b[u] = temp1 > temp2 ? temp1 : temp2;
				direction_line[de] = temp1 > temp2 ? 3 : 2;

#ifdef GAP_POS_SCORING
                if (posSpecificGaps) {
                    temp1 = h_c[b] - gIns[i];
                } else {
#endif
                    temp1 = h_c[b] - gap_open;
#ifdef GAP_POS_SCORING
                }
#endif

                temp2 = f - gap_extend;
				f = temp1 > temp2 ? temp1 : temp2;
				direction_line[df] = temp1 > temp2 ? 5 : 4;

                f1 = f > 0 ? f : 0;
#ifdef GAP_POS_SCORING
                if (posSpecificGaps) {
                    e1 = std::max(0, e_b[u] - gDelClose[i + 1]);
                } else {
#endif
                    e1 = e_b[u] > 0 ? e_b[u] : 0;
#ifdef GAP_POS_SCORING
                }
#endif

				temp1 = e1 > f1 ? e1 : f1;

				//TODO: careful with the variable names
				if (type == PROFILE_PROFILE) {
				    // both db_sequence and query_sequence must be the consensus sequence
                    int32_t minScore = std::min(mat[db_sequence[j] * qry_n + (queryStart + i)], db_mat[query_consens_sequence[i] * tgt_n + (targetStart + j)]);
                    int32_t absMinScore = abs(minScore);
                    int32_t maxScore = std::max(mat[db_sequence[j] * qry_n + (queryStart + i)], db_mat[query_consens_sequence[i] * tgt_n + (targetStart + j)]);
                    subScore = (absMinScore + minScore) + (absMinScore + maxScore) ;
                    subScore = ((subScore + 1) / 2) - absMinScore;
					temp2 = h_b[d] + subScore;
				} else if (type == PROFILE_SEQ) {
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
