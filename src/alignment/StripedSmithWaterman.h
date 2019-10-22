
/* $Id: smith_waterman_sse2.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

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
  Written by Michael Farrar, 2006.
  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
*/

#ifndef SMITH_WATERMAN_SSE2_H
#define SMITH_WATERMAN_SSE2_H

#include <climits>

#include <cstdio>
#include <cstdlib>

#if !defined(__APPLE__) && !defined(__llvm__)
#include <malloc.h>
#endif

#include "simd.h"
#include "BaseMatrix.h"

#include "Sequence.h"
#include "EvalueComputation.h"
typedef struct {
    short qStartPos;
    short dbStartPos;
    short qEndPos;
    short dbEndPos;
} aln_t;


typedef struct {
    uint32_t score1;
    uint32_t score2;
    int32_t dbStartPos1;
    int32_t dbEndPos1;
    int32_t	qStartPos1;
    int32_t qEndPos1;
    int32_t ref_end2;
    float qCov;
    float tCov;
    uint32_t* cigar;
    int32_t cigarLen;
    double evalue;
} s_align;

class SmithWaterman{
public:

    SmithWaterman(size_t maxSequenceLength, int aaSize, bool aaBiasCorrection);
    ~SmithWaterman();

    // prints a __m128 vector containing 8 signed shorts
    static void printVector (__m128i v);

    // prints a __m128 vector containing 8 unsigned shorts, added 32768
    static void printVectorUS (__m128i v);

    static unsigned short sse2_extract_epi16(__m128i v, int pos);

    // The dynamic programming matrix entries for the query and database sequences are stored sequentially (the order see the Farrar paper).
    // This function calculates the index within the dynamic programming matrices for the given query and database sequence position.
    static inline int midx (int qpos, int dbpos, int iter){
        return dbpos * (8 * iter) + (qpos % iter) * 8 + (qpos / iter);
    }

    // @function	ssw alignment.
    /*!	@function	Do Striped Smith-Waterman alignment.

     @param	db_sequence	pointer to the target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of
     function ssw_init

     @param	db_length	length of the target sequence

     @param	gap_open	the absolute value of gap open penalty

     @param	gap_extend	the absolute value of gap extension penalty

     @param	alignmentMode	bitwise FLAG; (from high to low) bit 5: when setted as 1, function ssw_align will return the best alignment
     beginning position; bit 6: when setted as 1, if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1
     < filterd), (whatever bit 5 is setted) the function will return the best alignment beginning position and
     cigar; bit 7: when setted as 1, if the best alignment score >= filters, (whatever bit 5 is setted) the function
     will return the best alignment beginning position and cigar; bit 8: when setted as 1, (whatever bit 5, 6 or 7 is
     setted) the function will always return the best alignment beginning position and cigar. When flag == 0, only
     the optimal and sub-optimal scores and the optimal alignment ending position will be returned.

     @param	filters	score filter: when bit 7 of flag is setted as 1 and bit 8 is setted as 0, filters will be used (Please check the
     decription of the flag parameter for detailed usage.)

     @param	filterd	distance filter: when bit 6 of flag is setted as 1 and bit 8 is setted as 0, filterd will be used (Please check
     the decription of the flag parameter for detailed usage.)

     @param	maskLen	The distance between the optimal and suboptimal alignment ending position >= maskLen. We suggest to use
     readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function will NOT
     return the suboptimal alignment information. Detailed description of maskLen: After locating the optimal
     alignment ending position, the suboptimal alignment score can be heuristically found by checking the second
     largest score in the array that contains the maximal score of each column of the SW matrix. In order to avoid
     picking the scores that belong to the alignments sharing the partial best alignment, SSW C library masks the
     reference loci nearby (mask length = maskLen) the best alignment ending position and locates the second largest
     score from the unmasked elements.

     @return	pointer to the alignment result structure

     @note	Whatever the parameter flag is setted, this function will at least return the optimal and sub-optimal alignment score,
     and the optimal alignment ending positions on target and query sequences. If both bit 6 and 7 of the flag are setted
     while bit 8 is not, the function will return cigar only when both criteria are fulfilled. All returned positions are
     0-based coordinate.
     */
    s_align  ssw_align (const unsigned char*db_sequence,
                        int32_t db_length,
                        const uint8_t gap_open,
                        const uint8_t gap_extend,
                        const uint8_t alignmentMode,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
                        const double filters,
                        EvalueComputation * filterd,
                        const int covMode, const float covThr,
                        const int32_t maskLen);


    /*!	@function computed ungapped alignment score

   @param	db_sequence	pointer to the target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of
   function ssw_init

   @param	db_length	length of the target sequence
   @return	max diagonal score
   */
   int ungapped_alignment(const unsigned char *db_sequence,
                          int32_t db_length);

  /*!	@function	Create the query profile using the query sequence.
   @param	read	pointer to the query sequence; the query sequence needs to be numbers
   @param	readLen	length of the query sequence
   @param	mat	pointer to the substitution matrix; mat needs to be corresponding to the read sequence
   @param	n	the square root of the number of elements in mat (mat has n*n elements)
   @param	score_size	estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if
   your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
   @return	pointer to the query profile structure
   @note	example for parameter read and mat:
   If the query sequence is: ACGTATC, the sequence that read points to can be: 1234142
   Then if the penalty for match is 2 and for mismatch is -2, the substitution matrix of parameter mat will be:
   //A  C  G  T
   2 -2 -2 -2 //A
   -2  2 -2 -2 //C
   -2 -2  2 -2 //G
   -2 -2 -2  2 //T
   mat is the pointer to the array {2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2}
   */
    void ssw_init(const Sequence *q, const int8_t *mat, const BaseMatrix *m, const int8_t score_size);


    static char cigar_int_to_op (uint32_t cigar_int);

    static uint32_t cigar_int_to_len (uint32_t cigar_int);


    static float computeCov(unsigned int startPos, unsigned int endPos, unsigned int len);

    s_align scoreIdentical(unsigned char *dbSeq, int L, EvalueComputation * evaluer, int alignmentMode);

    static void seq_reverse(int8_t * reverse, const int8_t* seq, int32_t end)	/* end is 0-based alignment ending position */
    {
        int32_t start = 0;
        while (LIKELY(start <= end)) {
            reverse[start] = seq[end];
            reverse[end] = seq[start];
            ++start;
            --end;
        }
    }


private:

    struct s_profile{
        simd_int* profile_byte;	// 0: none
        simd_int* profile_word;	// 0: none
        simd_int* profile_rev_byte;	// 0: none
        simd_int* profile_rev_word;	// 0: none
        int8_t* query_sequence;
        int8_t* query_rev_sequence;
        int8_t* composition_bias;
        int8_t* composition_bias_rev;
        int8_t* mat;
        // Memory layout of if mat + queryProfile is qL * AA
        //    Query length
        // A  -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
        // C  -1  -4   2   5  -3  -2   0  -3   1  -3  -2   0  -1   2   0   0  -1  -3  -4  -2
        // ...
        // Y -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
        // Memory layout of if mat + sub is AA * AA
        //     A   C    ...                                                                Y
        // A  -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
        // C  -1  -4   2   5  -3  -2   0  -3   1  -3  -2   0  -1   2   0   0  -1  -3  -4  -2
        // ...
        // Y -1  -3  -2  -1  -4  -2  -2  -3  -1  -3  -2  -2   7  -1  -2  -1  -1  -2  -5  -3
        int8_t* mat_rev; // needed for queryProfile
        int32_t query_length;
        int32_t sequence_type;
        int32_t alphabetSize;
        uint8_t bias;
        short ** profile_word_linear;
    };
    simd_int* vHStore;
    simd_int* vHLoad;
    simd_int* vE;
    simd_int* vHmax;
    uint8_t * maxColumn;

    typedef struct {
        uint16_t score;
        int32_t ref;	 //0-based position
        int32_t read;    //alignment ending position on read, 0-based
    } alignment_end;


    typedef struct {
        uint32_t* seq;
        int32_t length;
    } cigar;

    /* Striped Smith-Waterman
     Record the highest score of each reference position.
     Return the alignment score and ending position of the best alignment, 2nd best alignment, etc.
     Gap begin and gap extension are different.
     wight_match > 0, all other weights < 0.
     The returned positions are 0-based.
     */
    std::pair<alignment_end, alignment_end> sw_sse2_byte (const unsigned char*db_sequence,
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
                                 int32_t maskLen);

    std::pair<alignment_end, alignment_end> sw_sse2_word (const unsigned char* db_sequence,
                                 int8_t ref_dir,	// 0: forward ref; 1: reverse ref
                                 int32_t db_length,
                                 int32_t query_length,
                                 const uint8_t gap_open, /* will be used as - */
                                 const uint8_t gap_extend, /* will be used as - */
                                 const simd_int*query_profile_byte,
                                 uint16_t terminate,
                                 int32_t maskLen);

    template <const unsigned int type>
    SmithWaterman::cigar *banded_sw(const unsigned char *db_sequence, const int8_t *query_sequence, const int8_t * compositionBias, int32_t db_length, int32_t query_length, int32_t queryStart, int32_t score, const uint32_t gap_open, const uint32_t gap_extend, int32_t band_width, const int8_t *mat, int32_t n);

    /*!	@function		Produce CIGAR 32-bit unsigned integer from CIGAR operation and CIGAR length
     @param	length		length of CIGAR
     @param	op_letter	CIGAR operation character ('M', 'I', etc)
     @return			32-bit unsigned integer, representing encoded CIGAR operation and length
     */
    inline uint32_t to_cigar_int (uint32_t length, char op_letter);

    s_profile* profile;


    const static unsigned int SUBSTITUTIONMATRIX = 1;
    const static unsigned int PROFILE = 2;

    template <typename T, size_t Elements, const unsigned int type>
    void createQueryProfile(simd_int *profile, const int8_t *query_sequence, const int8_t * composition_bias, const int8_t *mat, const int32_t query_length, const int32_t aaSize, uint8_t bias, const int32_t offset, const int32_t entryLength);

    float *tmp_composition_bias;
    short * profile_word_linear_data;
    bool aaBiasCorrection;
};
#endif /* SMITH_WATERMAN_SSE2_H */
