
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
extern "C" {
#include <sunmedia_intrin.h>
}
#else
extern "C" {
#include <emmintrin.h>
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <limits.h> 

#include "../commons/Sequence.h"

typedef struct {
    short qStartPos;
    short dbStartPos;
    short qEndPos;
    short dbEndPos;
} aln_t;

// calculates 8 cells of the SW dynamic programming matrix in parallel
int smith_waterman_sse2_word(char* query_id,
                         int *     query_sequence,
                         unsigned short *    query_profile_word,
                         const int                 query_length,
                         char* db_id,
                         int *     db_sequence,
                         const int                 db_length,
                         unsigned short      gap_open,
                         unsigned short      gap_extend,
                         void *              workspace_void,
                         void *              Hmatrix,
                         void *              Ematrix,
                         void *              Fmatrix,
                         unsigned short *             qmaxpos,
                         unsigned short *            dbmaxpos
                         );


// calculates the start positions of the alignment in the query and in the database sequences given the dynamic programming matrices and the end positions of the alignment
void traceback_word(short* Hmatrix, 
        short* Ematrix,
        short* Fmatrix,
        Sequence* querySeq,
        Sequence* dbSeq,
        unsigned short * query_profile_word,
        unsigned short qmaxpos, 
        unsigned short dbmaxpos, 
        unsigned short gap_open, 
        unsigned short gap_extend,
        unsigned short* qstartpos,
        unsigned short* dbstartpos, 
        int* aaIds,
        int* overflow_warning);

// prints a __m128 vector containing 8 signed shorts
void printVector (__m128i v);

// prints a __m128 vector containing 8 unsigned shorts, added 32768
void printVectorUS (__m128i v);

unsigned short sse2_extract_epi16(__m128i v, int pos);

// The dynamic programming matrix entries for the query and database sequences are stored sequentially (the order see the Farrar paper).
// This function calculates the index within the dynamic programming matrices for the given query and database sequence position.
inline int midx (int qpos, int dbpos, int iter){
    return dbpos * (8 * iter) + (qpos % iter) * 8 + (qpos / iter);
}
#endif /* SMITH_WATERMAN_SSE2_H */
