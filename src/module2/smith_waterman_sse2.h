
/* $Id: smith_waterman_sse2.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

/******************************************************************
  Copyright 2006 by Michael Farrar.  All rights reserved.
  This program may not be sold or incorporated into a commercial product,
  in whole or in part, without written consent of Michael Farrar.  For 
  further information regarding permission for use or reproduction, please 
  contact: Michael Farrar at farrar.michael@gmail.com.
*******************************************************************/

/*
  Written by Michael Farrar, 2006.
  Please send bug reports and/or suggestions to farrar.michael@gmail.com.
*/

#ifndef SMITH_WATERMAN_SSE2_H
#define SMITH_WATERMAN_SSE2_H

#ifdef __SUNPRO_C
#include <sunmedia_intrin.h>
#else
#include <emmintrin.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>


typedef struct {
    short qStartPos;
    short dbStartPos;
    short qEndPos;
    short dbEndPos;
} aln_t;

// calculates 8 cells of the SW dynamic programming matrix in parallel
int smith_waterman_sse2_word(const unsigned char *     query_sequence,
                         unsigned short *    query_profile_word,
                         const int                 query_length,
                         const unsigned char *     db_sequence,
                         const int                 db_length,
                         unsigned short      gap_open,
                         unsigned short      gap_extend,
                         void *              workspace_void,
                         void *              Hmatrix,
                         void *              Ematrix,
                         void *              Fmatrix,
                         short *             qmaxpos,
                         short *            dbmaxpos
                         );


// calculates 16 cells of the SW dynamic programming matrix in parallel
// deprecated (no traceback possible)
int smith_waterman_sse2_byte(const unsigned char *     query_sequence,
                         unsigned char *     query_profile_byte,
                         const int                 query_length,
                         const unsigned char *     db_sequence,
                         const int                 db_length,
                         unsigned char       bias,
                         unsigned char       gap_open,
                         unsigned char       gap_extend,
                         void *              workspace_void);

// calculates the start positions of the alignment in the query and in the database sequences given the dynamic programming matrices and the end positions of the alignment
void traceback_word(short* Hmatrix, 
        short* Ematrix,
        short* Fmatrix,
        unsigned char* query_sequence,
        unsigned short * query_profile_word,
        short qLen,
        const unsigned char* db_sequence,
        short dbLen,
        short qmaxpos, 
        short dbmaxpos, 
        unsigned short gap_open, 
        unsigned short gap_extend,
        short* qstartpos,
        short* dbstartpos);

// prints a __m128 vector containing 8 signed shorts
void printvector (__m128i v);

// The dynamic programming matrix entries for the query and database sequences are stored sequentially (the order see the Farrar paper).
// This function calculates the index within the dynamic programming matrices for the given query and database sequence position.
inline int midx (int qpos, int dbpos, int iter){
    return dbpos * (8 * iter) + (qpos % iter) * 8 + (qpos / iter);
}
#endif /* SMITH_WATERMAN_SSE2_H */
