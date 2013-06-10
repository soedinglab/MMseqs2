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

#include "smith_waterman_sse2.h"

// amino acids in the sequences are numerically coded
int
smith_waterman_sse2_word(const unsigned char *     query_sequence,
                         unsigned short *    query_profile_word,
                         const int                 query_length,
                         const unsigned char *     db_sequence,
                         const int                 db_length,
                         unsigned short      gap_open,
                         unsigned short      gap_extend,
                         void *             workspace_void,
                         void *             Hmatrix,
                         void *             Ematrix,
                         void *             Fmatrix,
                         short *            qmaxpos,
                         short *            dbmaxpos)
{
    int     i, j, k;
    short   score;

    int     cmp;
    int     iter = (query_length + 7) / 8;
    
    __m128i *p;
    __m128i *workspace = (__m128i *) workspace_void;

    __m128i E, F, H;

    __m128i v_maxscore;
    __m128i v_gapopen;
    __m128i v_gapextend;

    __m128i v_min;
    __m128i v_minimums;
    __m128i v_temp;

    __m128i *pHLoad, *pHStore;
    __m128i *pE;

    __m128i *pScore;

    __m128i v_maxtemp;
    __m128i v_maxmask;
    __m128i v_qpos = _mm_setzero_si128();
    __m128i v_dbpos = _mm_setzero_si128();
    __m128i v_qmaxpos = _mm_setzero_si128();
    __m128i v_dbmaxpos = _mm_setzero_si128();

    __m128i* HStore = (__m128i*) Hmatrix;
    __m128i* EStore = (__m128i*) Ematrix;
    __m128i* FStore = (__m128i*) Fmatrix;
    __m128i f_temp;

    /* Load gap opening penalty to all elements of a constant */
    v_gapopen = _mm_insert_epi16 (v_gapopen, gap_open, 0);
    v_gapopen = _mm_shufflelo_epi16 (v_gapopen, 0);
    v_gapopen = _mm_shuffle_epi32 (v_gapopen, 0);

    /* Load gap extension penalty to all elements of a constant */
    v_gapextend = _mm_insert_epi16 (v_gapextend, gap_extend, 0);
    v_gapextend = _mm_shufflelo_epi16 (v_gapextend, 0);
    v_gapextend = _mm_shuffle_epi32 (v_gapextend, 0);

    /* load v_maxscore with the zeros.  since we are using signed */
    /*  math, we will bias the maxscore to -32768 so we have the */
    /*  full range of the short. */
    v_maxscore = _mm_cmpeq_epi16 (v_maxscore, v_maxscore);
    v_maxscore = _mm_slli_epi16 (v_maxscore, 15);

    v_minimums = _mm_shuffle_epi32 (v_maxscore, 0);

    v_min = _mm_shuffle_epi32 (v_maxscore, 0);
    v_min = _mm_srli_si128 (v_min, 14);

    /* Zero out the storage vector */
    k = 2 * iter;

    p = workspace;
    for (i = 0; i < k; i++)
    {
        _mm_store_si128 (p++, v_maxscore);
    }

    pE = workspace;
    pHStore = pE + iter;
    pHLoad = pHStore + iter;

    for (i = 0; i < db_length; ++i)
    {

        /* fetch first data asap. */
        // get the score row for the amino acid i (the row contains the scores for the amino acid i at each position in the query sequence)
        pScore = ((__m128i *) query_profile_word) + (db_sequence[i] * iter);

        /* bias all elements in F to -32768 */
        F = _mm_cmpeq_epi16 (F, F);
        F = _mm_slli_epi16 (F, 15);

        /* load the next h value */
        H = _mm_load_si128 (pHStore + iter - 1);
        H = _mm_slli_si128 (H, 2);
        H = _mm_or_si128 (H, v_min);

        p = pHLoad;
        pHLoad = pHStore;
        pHStore = p;

        for (j = 0; j < iter; j++)
        {
            // save F values for the current position
            _mm_store_si128 (FStore + j, F);

            // set current query sequence positions
            v_qpos = _mm_setr_epi16(j, iter+j, 2*iter+j, 3*iter+j, 4*iter+j, 5*iter+j, 6*iter+j, 7*iter+j);

            // set current db sequence position
            v_dbpos = _mm_set1_epi16(i);

            /* load E values */
            E = _mm_load_si128 (pE + j);

            // save E values for the current position
            _mm_store_si128 (EStore + j, E);

            /* add score to H */
            H = _mm_adds_epi16 (H, *pScore++);

            /* Update highest score encountered this far */
            v_maxscore = _mm_max_epi16(v_maxscore, H);

            // save the positions of the highest score
            v_maxmask = _mm_cmpeq_epi16(H, v_maxscore);
            
            v_dbpos = _mm_and_si128(v_dbpos, v_maxmask);
            v_qpos = _mm_and_si128(v_qpos, v_maxmask);
            
            v_dbmaxpos = _mm_andnot_si128(v_maxmask, v_dbmaxpos);
            v_qmaxpos = _mm_andnot_si128(v_maxmask, v_qmaxpos);

            v_dbmaxpos = _mm_max_epi16(v_dbmaxpos, v_dbpos);
            v_qmaxpos = _mm_max_epi16(v_qmaxpos, v_qpos);

            /* get max from H, E and F */
            H = _mm_max_epi16 (H, E);
            H = _mm_max_epi16 (H, F);

            /* save H values */
            _mm_store_si128 (pHStore + j, H);
            _mm_store_si128 (HStore + j, H);

            /* subtract the gap open penalty from H */
            H = _mm_subs_epi16 (H, v_gapopen);

            /* update E value */
            E = _mm_subs_epi16 (E, v_gapextend);
            E = _mm_max_epi16 (E, H);

            /* update F value */
            F = _mm_subs_epi16 (F, v_gapextend);
            F = _mm_max_epi16 (F, H);

            /* save E values */
            _mm_store_si128 (pE + j, E);

            /* load the next h value */
            H = _mm_load_si128 (pHLoad + j);
        }

        /* reset pointers to the start of the saved data */
        j = 0;
        H = _mm_load_si128 (pHStore + j);

        /*  the computed F value is for the given column.  since */
        /*  we are at the end, we need to shift the F value over */
        /*  to the next column. */
        F = _mm_slli_si128 (F, 2);
        F = _mm_or_si128 (F, v_min);

        v_temp = _mm_subs_epi16 (H, v_gapopen);
        v_temp = _mm_cmpgt_epi16 (F, v_temp);
        cmp  = _mm_movemask_epi8 (v_temp);

        while (cmp != 0x0000) 
        {
            E = _mm_load_si128 (pE + j);

            // save the new maximum for F in the matrix
            f_temp = _mm_load_si128 (FStore + j);
            f_temp = _mm_max_epi16 (f_temp, F);
            _mm_store_si128 (FStore + j, f_temp);

            H = _mm_max_epi16 (H, F);

            /* save H values */
            _mm_store_si128 (pHStore + j, H);
            _mm_store_si128 (HStore + j, H);

            /* update E in case the new H value would change it */
            // for the next column
            H = _mm_subs_epi16 (H, v_gapopen);
            E = _mm_max_epi16 (E, H);
            _mm_store_si128 (pE + j, E);

            /* update F value */
            F = _mm_subs_epi16 (F, v_gapextend);

            j++;
            if (j >= iter)
            {
                j = 0;
                F = _mm_slli_si128 (F, 2);
                F = _mm_or_si128 (F, v_min);
            }
            H = _mm_load_si128 (pHStore + j);

            v_temp = _mm_subs_epi16 (H, v_gapopen);
            v_temp = _mm_cmpgt_epi16 (F, v_temp);
            cmp  = _mm_movemask_epi8 (v_temp);
        }

        // move to the next matrix column
        HStore = HStore + iter;
        EStore = EStore + iter;
        FStore = FStore + iter;
    }

    v_maxtemp = _mm_setzero_si128();
    v_maxtemp = _mm_add_epi16(v_maxtemp, v_maxscore);

    /* find largest score in the v_maxscore vector */
    v_temp = _mm_srli_si128 (v_maxscore, 8);
    v_maxscore = _mm_max_epi16 (v_maxscore, v_temp);
    v_temp = _mm_srli_si128 (v_maxscore, 4);
    v_maxscore = _mm_max_epi16 (v_maxscore, v_temp);
    v_temp = _mm_srli_si128 (v_maxscore, 2);
    v_maxscore = _mm_max_epi16 (v_maxscore, v_temp);

    /* extract the largest score */
    score = _mm_extract_epi16 (v_maxscore, 0);

    // find the sequence positions of the largest score
    v_maxmask = _mm_set1_epi16(score);
    v_maxmask = _mm_cmpeq_epi16(v_maxtemp, v_maxmask);

    v_qmaxpos = _mm_and_si128(v_qmaxpos, v_maxmask);
    v_dbmaxpos = _mm_and_si128(v_dbmaxpos, v_maxmask);

    // find the single non-zero entry in the query and database position vectors
    v_temp = _mm_srli_si128 (v_qmaxpos, 8);
    v_qmaxpos = _mm_max_epi16 (v_qmaxpos, v_temp);
    v_temp = _mm_srli_si128 (v_qmaxpos, 4);
    v_qmaxpos = _mm_max_epi16 (v_qmaxpos, v_temp);
    v_temp = _mm_srli_si128 (v_qmaxpos, 2);
    v_qmaxpos = _mm_max_epi16 (v_qmaxpos, v_temp);

    v_temp = _mm_srli_si128 (v_dbmaxpos, 8);
    v_dbmaxpos = _mm_max_epi16 (v_dbmaxpos, v_temp);
    v_temp = _mm_srli_si128 (v_dbmaxpos, 4);
    v_dbmaxpos = _mm_max_epi16 (v_dbmaxpos, v_temp);
    v_temp = _mm_srli_si128 (v_dbmaxpos, 2);
    v_dbmaxpos = _mm_max_epi16 (v_dbmaxpos, v_temp);

    // extract the query and db sequence positions
    *qmaxpos = _mm_extract_epi16 (v_qmaxpos, 0);
    *dbmaxpos = _mm_extract_epi16 (v_dbmaxpos, 0);

    /* return largest score biased by 32768 */
    return score + 32768;
}

void traceback_word(short* H, 
        short* E, 
        short* F,
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
        short* dbstartpos){

    int iter = (qLen + 7) / 8;

    int qpos, dbpos, idx;

    __m128i *pScore;

    qpos = qmaxpos;
    dbpos = dbmaxpos;
/*
    int imin = 0;
    int imax = qLen;

    int jmin = 0;
    int jmax = 14;

    if (imax <= qLen && jmax <= dbLen){
        printf("           \t");
        for (int j = jmin; j< jmax; j++){
            printf("%11d\t", j);
        }
        printf("\n");
        for (int i = imin; i < imax; i++){
            printf("%11d\t", i);
            for (int j = jmin; j< jmax; j++){
                idx = j * (8 * iter) + (i % iter) * 8 + (i / iter); 

                printf("%2d|%2d|%2d|%2d\t", H[idx] + 32768, E[idx] + 32768, F[idx] + 32768, *((short*)query_profile_word + (db_sequence[j] * iter * 8 + i%iter * 8 + i/iter)));
            }
            printf("\n");
        }
        printf("           \t");
        for (int j = jmin; j< jmax; j++){
            printf("%11d\t", j);
        }
        printf("\n");

    }
*/
    idx = midx(qpos, dbpos, iter); 
//    printf("qpos:%d, dbpos:%d, index: %d\n", qpos, dbpos, idx);
    while (qpos != 0 && dbpos != 0 && H[idx] + 32768 != 0){
        printf("%d %d %d\n", qpos, dbpos, H[idx] + 32768);
        // match between q[i] and db[j]
        if (H[idx] == (H[midx(qpos-1, dbpos-1, iter)] + *((short*)query_profile_word + db_sequence[dbpos] * iter * 8 + qpos%iter * 8 + qpos/iter)) ){ // H[i][j] == H[i-1][j-1] + score(q[i], db[j])
            qpos--;
            dbpos--;
            idx = midx(qpos, dbpos, iter);
        }
        // continue with the E matrix = gap in the db sequence
        else if (H[idx] == E[idx]){
//            printf("Entering E\n");
            while (qpos != 0 && dbpos != 0 && H[idx] + 32768 != 0){
                if (E[idx] == E[midx(qpos, dbpos-1, iter)] - gap_extend){
                    dbpos--;
                    idx = midx(qpos, dbpos, iter);
//                    printf("%d %d %d\n", qpos, dbpos, E[idx] + 32768);
                }
                else if (E[idx] == H[midx(qpos, dbpos-1, iter)] - gap_open){
                    // leave E matrix
                    dbpos--;
                    idx = midx(qpos, dbpos, iter);
//                    printf("%d %d %d\n", qpos, dbpos, E[idx] + 32768);
//                    printf("Leaving E\n");
                    break;
                }
                else{
                    printf("ERROR 1\n");
                    exit(1);
                }
            }
        }
        // continue with the F matrix = gap in the query sequence
        else if (H[idx] == F[idx]){
//            printf("Entering F\n");
            while (qpos != 0 && dbpos != 0 && H[idx] + 32768 != 0){
                if (F[idx] == F[midx(qpos-1, dbpos, iter)] - gap_extend){
                    qpos--;
                    idx = midx(qpos, dbpos, iter);
//                    printf("%d %d %d\n", qpos, dbpos, F[idx] + 32768);
                }
                else if (F[idx] == H[midx(qpos-1, dbpos, iter)] - gap_open){
                    // leave F matrix
                    qpos--;
                    idx = midx(qpos, dbpos, iter);
//                    printf("%d %d %d\n", qpos, dbpos, F[idx] + 32768);
//                    printf("Leaving F\n");
                    break;
                }
                else{
                    printf("ERROR 2\n");
                    exit(1);
                }
            }
        }
        else{
            printf("ERROR 3\n");
            exit(1);
        }
    }

    *qstartpos = qpos;
    *dbstartpos = dbpos;
}

void printvector(__m128i v){

    short val = _mm_extract_epi16 (v, 0);
    printf("%d %d \n", 0, val);
   
    val = _mm_extract_epi16 (v, 1);
    printf("%d %d\n", 1, val);

    val = _mm_extract_epi16 (v, 2);
    printf("%d %d\n", 2, val);

    val = _mm_extract_epi16 (v, 3);
    printf("%d %d\n", 3, val);

    val = _mm_extract_epi16 (v, 4);
    printf("%d %d\n", 4, val);

    val = _mm_extract_epi16 (v, 5);
    printf("%d %d\n", 5, val);

    val = _mm_extract_epi16 (v, 6);
    printf("%d %d\n", 6, val);

    val = _mm_extract_epi16 (v, 7);
    printf("%d %d\n", 7, val);

}



int
smith_waterman_sse2_byte(const unsigned char *     query_sequence,
                         unsigned char *     query_profile_byte,
                         const int                 query_length,
                         const unsigned char *     db_sequence,
                         const int                 db_length,
                         unsigned char       bias,
                         unsigned char       gap_open,
                         unsigned char       gap_extend,
                         void *             workspace_void)
{
    int     i, j, k;
    int     score;

    int     dup;
    int     cmp;
    int     iter = (query_length + 15) / 16;
    
    __m128i *p;
    __m128i *workspace = (__m128i *) workspace_void;

    __m128i E, F, H;

    __m128i v_maxscore;
    __m128i v_bias;
    __m128i v_gapopen;
    __m128i v_gapextend;

    __m128i v_temp;
    __m128i v_zero;

    __m128i *pHLoad, *pHStore;
    __m128i *pE;

    __m128i *pScore;

    /* Load the bias to all elements of a constant */
    dup    = ((short) bias << 8) | bias;
    v_bias = _mm_insert_epi16 (v_bias, dup, 0);
    v_bias = _mm_shufflelo_epi16 (v_bias, 0);
    v_bias = _mm_shuffle_epi32 (v_bias, 0);

    /* Load gap opening penalty to all elements of a constant */
    dup  = ((short) gap_open << 8) | gap_open;
    v_gapopen = _mm_insert_epi16 (v_gapopen, dup, 0);
    v_gapopen = _mm_shufflelo_epi16 (v_gapopen, 0);
    v_gapopen = _mm_shuffle_epi32 (v_gapopen, 0);

    /* Load gap extension penalty to all elements of a constant */
    dup  = ((short) gap_extend << 8) | gap_extend;
    v_gapextend = _mm_insert_epi16 (v_gapextend, dup, 0);
    v_gapextend = _mm_shufflelo_epi16 (v_gapextend, 0);
    v_gapextend = _mm_shuffle_epi32 (v_gapextend, 0);

    /* initialize the max score */
    v_maxscore = _mm_xor_si128 (v_maxscore, v_maxscore);

    /* create a constant of all zeros for comparison */
    v_zero = _mm_xor_si128 (v_zero, v_zero);

    /* Zero out the storage vector */
    k = iter * 2;

    p = workspace;
    for (i = 0; i < k; i++)
    {
        _mm_store_si128 (p++, v_maxscore);
    }

    pE = workspace;
    pHStore = pE + iter;
    pHLoad = pHStore + iter;

    for (i = 0; i < db_length; ++i)
    {
        /* fetch first data asap. */
        pScore = (__m128i *) query_profile_byte + db_sequence[i] * iter;

        /* zero out F value. */
        F = _mm_xor_si128 (F, F);

        /* load the next h value */
        H = _mm_load_si128 (pHStore + iter - 1);
        H = _mm_slli_si128 (H, 1);

        p = pHLoad;
        pHLoad = pHStore;
        pHStore = p;

        for (j = 0; j < iter; j++)
        {
            /* load values E. */
            E = _mm_load_si128 (pE + j);

            /* add score to H */
            H = _mm_adds_epu8 (H, *pScore++);
            H = _mm_subs_epu8 (H, v_bias);

            /* Update highest score encountered this far */
            v_maxscore = _mm_max_epu8 (v_maxscore, H);

            /* get max from H, E and F */
            H = _mm_max_epu8 (H, E);
            H = _mm_max_epu8 (H, F);

            /* save H values */
            _mm_store_si128 (pHStore + j, H);

            /* subtract the gap open penalty from H */
            H = _mm_subs_epu8 (H, v_gapopen);

            /* update E value */
            E = _mm_subs_epu8 (E, v_gapextend);
            E = _mm_max_epu8 (E, H);

            /* update F value */
            F = _mm_subs_epu8 (F, v_gapextend);
            F = _mm_max_epu8 (F, H);

            /* save E values */
            _mm_store_si128 (pE + j, E);

            /* load the next h value */
            H = _mm_load_si128 (pHLoad + j);
        }

        /* reset pointers to the start of the saved data */
        j = 0;
        H = _mm_load_si128 (pHStore + j);

        /*  the computed F value is for the given column.  since */
        /*  we are at the end, we need to shift the F value over */
        /*  to the next column. */
        F = _mm_slli_si128 (F, 1);
        v_temp = _mm_subs_epu8 (H, v_gapopen);
        v_temp = _mm_subs_epu8 (F, v_temp);
        v_temp = _mm_cmpeq_epi8 (v_temp, v_zero);
        cmp  = _mm_movemask_epi8 (v_temp);

        while (cmp != 0xffff) 
        {
            E = _mm_load_si128 (pE + j);

            H = _mm_max_epu8 (H, F);

            /* save H values */
            _mm_store_si128 (pHStore + j, H);

            /* update E in case the new H value would change it */
            H = _mm_subs_epu8 (H, v_gapopen);
            E = _mm_max_epu8 (E, H);
            _mm_store_si128 (pE + j, E);

            /* update F value */
            F = _mm_subs_epu8 (F, v_gapextend);

            j++;
            if (j >= iter)
            {
                j = 0;
                F = _mm_slli_si128 (F, 1);
            }
            H = _mm_load_si128 (pHStore + j);

            v_temp = _mm_subs_epu8 (H, v_gapopen);
            v_temp = _mm_subs_epu8 (F, v_temp);
            v_temp = _mm_cmpeq_epi8 (v_temp, v_zero);
            cmp  = _mm_movemask_epi8 (v_temp);
        }
    }

    /* find largest score in the v_maxscore vector */
    v_temp = _mm_srli_si128 (v_maxscore, 8);
    v_maxscore = _mm_max_epu8 (v_maxscore, v_temp);
    v_temp = _mm_srli_si128 (v_maxscore, 4);
    v_maxscore = _mm_max_epu8 (v_maxscore, v_temp);
    v_temp = _mm_srli_si128 (v_maxscore, 2);
    v_maxscore = _mm_max_epu8 (v_maxscore, v_temp);
    v_temp = _mm_srli_si128 (v_maxscore, 1);
    v_maxscore = _mm_max_epu8 (v_maxscore, v_temp);

    /* store in temporary variable */
    score = _mm_extract_epi16 (v_maxscore, 0);
    score = score & 0x00ff;

    /*  check if we might have overflowed */
    if (score + bias >= 255)
    {
        score = 255;
    }

    /* return largest score */
    return score;
}
