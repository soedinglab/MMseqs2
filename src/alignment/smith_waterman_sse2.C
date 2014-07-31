/******************************************************************
  Copyright 2006 by Michael Farrar.  All rights reserved.
  This program may not be sold or incorporated into a commercial product,
  in whole or in part, without written consent of Michael Farrar. 
 *******************************************************************/

/*
   Written by Michael Farrar, 2006 (alignment), and Maria Hauser, 2012 (traceback, modified alignment).
   Please send bug reports and/or suggestions to mhauser@genzentrum.lmu.de.
   */

#include "smith_waterman_sse2.h"

// amino acids in the sequences are numerically coded
    int
smith_waterman_sse2_word(char* query_id, 
        int *     query_sequence,
        unsigned short *    query_profile_word,
        const int                 query_length,
        char* db_id, 
        int *     db_sequence,
        const int                 db_length,
        unsigned short      gap_open,
        unsigned short      gap_extend,
        void *             workspace_void,
        void *             Hmatrix,
        void *             Ematrix,
        void *             Fmatrix,
        unsigned short *            qmaxpos,
        unsigned short *            dbmaxpos)
{
    unsigned short     i, j, k;
    short   score;

    int     cmp;
    // iteration length
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

    __m128i v1 = _mm_setzero_si128();
    __m128i v2 = _mm_setzero_si128();
    __m128i v3 = _mm_setzero_si128();

    v1 = v_maxscore;
    v2 = v_qmaxpos;
    v3 = v_dbmaxpos;

    v_maxtemp = _mm_setzero_si128();
    v_maxtemp = _mm_adds_epi16(v_maxtemp, v_maxscore);

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

    /* return the largest score biased by 32768 */
    return score + 32768;
}

void traceback_word(short* H, 
        short* E, 
        short* F,
        Sequence* query,
        Sequence* dbSeq,
        unsigned short * query_profile_word,
        unsigned short qmaxpos, 
        unsigned short dbmaxpos, 
        unsigned short gap_open, 
        unsigned short gap_extend,
        unsigned short* qstartpos,
        unsigned short* dbstartpos,
        int* aaIds,
        int* overflow_warning){

    int* query_sequence = query->int_sequence;
    char* queryDbKey = query->getDbKey();
    int qLen = query->L;
    int* db_sequence = dbSeq->int_sequence;
    char* dbDbKey = dbSeq->getDbKey();
    int dbLen = dbSeq->L;

    // number of iterations
    int iter = (qLen + 7) / 8;

    int qpos, dbpos, idx;

    qpos = qmaxpos;
    dbpos = dbmaxpos;

    /*
    for (int i = 0; i < query->L; i++){
        for (int j = 0; j < dbSeq->L; j++){
            int idx = midx(i, j, iter);
            printf("%d  ", H[idx] + 32768);
        }
        printf("\n");
    }
    printf("\n");
    */

    *aaIds = (query_sequence[qpos] == db_sequence[dbpos]);
    // dynamic programming matrix index, depending on positions in the sequences and iteration length
    idx = midx(qpos, dbpos, iter); 
    while (qpos > 0 && dbpos > 0){
        // match between q[i] and db[j]
        if (H[idx] == (H[midx(qpos-1, dbpos-1, iter)] + *((short*)query_profile_word + db_sequence[dbpos] * iter * 8 + qpos%iter * 8 + qpos/iter)) ){ // H[i][j] == H[i-1][j-1] + score(q[i], db[j])
            qpos--;
            dbpos--;
            idx = midx(qpos, dbpos, iter);
            *aaIds += (query_sequence[qpos] == db_sequence[dbpos]);
        }
        // if H[idx] = 0 then only a match is checked (this is meant for allowing X-X matches at the beginning of the alignment, they have score 0)
        // else, identical sequences with x at the beginning of the sequence will never produce alignment coverage 1.0
        // if there is no match with score 0 then traceback is aborted
        else if (H[idx] + 32768 == 0)
            break;
        // continue with the E matrix = gap in the db sequence
        else if (H[idx] == E[idx]){
            while (qpos != 0 && dbpos != 0 && H[idx] + 32768 != 0){
                if (E[idx] == E[midx(qpos, dbpos-1, iter)] - gap_extend){
                    dbpos--;
                    idx = midx(qpos, dbpos, iter);
                }
                else if (E[idx] == H[midx(qpos, dbpos-1, iter)] - gap_open){
                    // leave E matrix
                    dbpos--;
                    idx = midx(qpos, dbpos, iter);
                    break;
                }
                else{
#pragma omp critical
                    {
                        printf("ERROR 1\n");
                        printf("query: %s, db seq: %s\n", queryDbKey, dbDbKey);
                        printf("qLen: %d, dbLen: %d\nqmaxpos: %d, dbmaxpos: %d\nqpos: %d, dbpos: %d\n", qLen, dbLen, qmaxpos, dbmaxpos, qpos, dbpos);
                        printf("score: %d\n", H[midx(qpos, dbpos, iter)] + 32768);
                        printf("short overflow warning set: %d\n", *overflow_warning);
                        exit(1);
                    }
                }
            }
        }
        // continue with the F matrix = gap in the query sequence
        else if (H[idx] == F[idx]){
            while (qpos != 0 && dbpos != 0 && H[idx] + 32768 != 0){
                if (F[idx] == F[midx(qpos-1, dbpos, iter)] - gap_extend){
                    qpos--;
                    idx = midx(qpos, dbpos, iter);
                }
                else if (F[idx] == H[midx(qpos-1, dbpos, iter)] - gap_open){
                    // leave F matrix
                    qpos--;
                    idx = midx(qpos, dbpos, iter);
                    break;
                }
                else{
#pragma omp critical
                    {
                        printf("ERROR 2\n");
                        printf("query: %s, db seq: %s\n", queryDbKey, dbDbKey);
                        printf("qLen: %d, dbLen: %d\nqmaxpos: %d, dbmaxpos: %d\nqpos: %d, dbpos: %d\n", qLen, dbLen, qmaxpos, dbmaxpos, qpos, dbpos);
                        printf("score: %d\n", H[midx(qpos, dbpos, iter)] + 32768);
                        printf("short overflow warning set: %d\n", *overflow_warning);
                        exit(1);
                    }
                }
            }
        }
        else{
            if (H[idx] == SHRT_MAX){
                // if we did not find any trace, then here is a short overflow
                qpos--;
                dbpos--;
                *overflow_warning = 1;
            }
            else{
#pragma omp critical
                {
                    printf("ERROR 3\n");
                    printf("query: %s, db seq: %s\n", queryDbKey, dbDbKey);
                    printf("qLen: %d, dbLen: %d\nqmaxpos: %d, dbmaxpos: %d\nqpos: %d, dbpos: %d\n", qLen, dbLen, qmaxpos, dbmaxpos, qpos, dbpos);
                    printf("score: %d\n", H[midx(qpos, dbpos, iter)] + 32768);
                    printf("short overflow warning set: %d\n", *overflow_warning);
                    exit(1);
                }
            }
        }
    }

    // don't terrify the user ;-)
//    if (*overflow_warning == 1)
//        printf("\nWARNING: short range overflow (query: %s, db seq: %s)\nThe alignment might be inaccurate!\n", queryDbKey, dbDbKey);
    
    *qstartpos = (unsigned short) qpos;
    *dbstartpos = (unsigned short) dbpos;
}

void printVector(__m128i v){
    for (int i = 0; i < 8; i++)
       printf("%d ", ((short) (sse2_extract_epi16(v, i)) + 32768));
    std::cout << "\n";
}

void printVectorUS(__m128i v){
    for (int i = 0; i < 8; i++)
        printf("%d ", (unsigned short) sse2_extract_epi16(v, i));
    std::cout << "\n";
}

unsigned short sse2_extract_epi16(__m128i v, int pos) {
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
    exit(1);
    // never executed
    return 0;
}
