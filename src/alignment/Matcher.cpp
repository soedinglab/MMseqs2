#include "Matcher.h"

Matcher::Matcher(SubstitutionMatrix* m){

    this->m = m;

    // memory for the sequence profile is allocated only once, 40000 is longer than the length of the currently longest sequence (36805 for Q3ASY8_CHLCH)
    // = 1.6 MB
    this->queryProfileWord = (unsigned short*) memalign (16, m->alphabetSize * 40000 * sizeof(unsigned short));

    // workspace memory for the alignment calculation (see Farrar code)
    void * workspace_memory  = (void *)memalign(16, 2 * 40000 * sizeof(__m128i) + 256);
    workspace = (void *) ((((size_t) workspace_memory) + 255) & (~0xff));
}

Matcher::~Matcher(){
    free(this->queryProfileWord);
    free(this->workspace);
}

Matcher::result_t Matcher::getSWResult(Sequence* query, Sequence* dbSeq, int seqDbSize){

    unsigned short gap_open = 11;
    unsigned short gap_extend = 1;

    calcQueryProfileWord(query);
    
    // allocate memory for the three dynamic programming matrices
    void* Hmatrix = memalign(16,(query->L + 7)/8 * dbSeq->L * sizeof(__m128i));   // 2GB für 36805*36805 (Q3ASY8_CHLCH)
    void* Ematrix = memalign(16,(query->L + 7)/8 * dbSeq->L * sizeof(__m128i));
    void* Fmatrix = memalign(16,(query->L + 7)/8 * dbSeq->L * sizeof(__m128i));
    
    short qStartPos = 0;
    short qEndPos = 0;
    short dbStartPos = 0;
    short dbEndPos = 0;
    // calculation of the score and traceback of the alignment
    int s = smith_waterman_sse2_word(query->int_sequence, this->queryProfileWord, query->L, 
            dbSeq->int_sequence, dbSeq->L, 
            gap_open, gap_extend, 
            workspace, 
            Hmatrix, Ematrix, Fmatrix, 
            &qEndPos, &dbEndPos);
    if (s > 0)
        traceback_word((short*) Hmatrix, (short*) Ematrix, (short*) Fmatrix, 
                query->int_sequence, this->queryProfileWord, query->L, 
                dbSeq->int_sequence, dbSeq->L, 
                qEndPos, dbEndPos, 
                gap_open, gap_extend, 
                &qStartPos, &dbStartPos);

    // calculation of the coverage and e-value
    float qcov = (qEndPos - qStartPos)/(float)query->L;
    float dbcov = (dbEndPos - dbStartPos)/(float)dbSeq->L;
    // bit factor is 2.0 (also used in the substitution matrix calculation)
    double evalue = (double)(seqDbSize * query->L * dbSeq->L) * fpow2((double)-s/2.0);

    free(Hmatrix);
    free(Ematrix);
    free(Fmatrix);

    result_t result = {std::string(dbSeq->dbKey), s, qcov, dbcov, evalue};
    return result;
}

void Matcher::calcQueryProfileWord(Sequence* query){

    int segLen = (query->L + 7) / 8;

    int a,h,i,j,k;
    for (a = 0; a < m->alphabetSize; ++a)
    {   
        h = a * segLen * 8;
        for (i = 0; i < segLen; ++i)
        {   
            j = i;
            for (k = 0; k < 8; ++k)
            {   
                if (j >= query->L)
                    queryProfileWord[h] = 0;
                else {
                    queryProfileWord[h] = (unsigned short) (m->subMatrix[query->int_sequence[j]][a]);
                }
                ++h;
                j += segLen;
            }
        }
    }
}
