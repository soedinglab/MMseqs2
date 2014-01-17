#include "Matcher.h"

Matcher::Matcher(BaseMatrix* m, int maxSeqLen){

    this->m = m;
    this->maxSeqLen = maxSeqLen;

    // memory for the sequence profile is allocated only once
    // = 1.6 MB
    this->queryProfileWord = (unsigned short*) memalign (16, m->alphabetSize * maxSeqLen * sizeof(unsigned short));

    // workspace memory for the alignment calculation (see Farrar code)
    this->workspace_memory  = (void *)memalign(16, 2 * maxSeqLen * sizeof(__m128i) + 256);
    this->workspace = (void *) ((((size_t) workspace_memory) + 255) & (~0xff));
}

Matcher::~Matcher(){
    free(this->queryProfileWord);
    free(this->workspace_memory);
}

Matcher::result_t Matcher::getSWResult(Sequence* query, Sequence* dbSeq, int seqDbSize){

    unsigned short gap_open = 10;
    unsigned short gap_extend = 1;

    calcQueryProfileWord(query);

    // allocate memory for the three dynamic programming matrices
    void* Hmatrix = memalign(16,(query->L + 7)/8 * dbSeq->L * sizeof(__m128i));   // 2.5GB für 36805*36805 (Q3ASY8_CHLCH)
    void* Ematrix = memalign(16,(query->L + 7)/8 * dbSeq->L * sizeof(__m128i));
    void* Fmatrix = memalign(16,(query->L + 7)/8 * dbSeq->L * sizeof(__m128i));

    unsigned short qStartPos = 0;
    unsigned short qEndPos = 0;
    unsigned short dbStartPos = 0;
    unsigned short dbEndPos = 0;
    int aaIds = 0;

    // calculation of the score and traceback of the alignment
    int s = smith_waterman_sse2_word(query->getDbKey(), query->int_sequence, this->queryProfileWord, query->L, 
            dbSeq->getDbKey(), dbSeq->int_sequence, dbSeq->L, 
            gap_open, gap_extend, 
            workspace, 
            Hmatrix, Ematrix, Fmatrix, 
            &qEndPos, &dbEndPos);

    traceback_word((short*) Hmatrix, (short*) Ematrix, (short*) Fmatrix, 
            query, dbSeq,
            this->queryProfileWord,
            qEndPos, dbEndPos, 
            gap_open, gap_extend, 
            &qStartPos, &dbStartPos, &aaIds);

    // calculation of the coverage and e-value
    float qcov = (std::min(query->L, (int) qEndPos) - qStartPos + 1)/ (float)query->L;
    float dbcov = (std::min(dbSeq->L, (int) dbEndPos) - dbStartPos + 1)/(float)dbSeq->L;
    double evalue = (double)(seqDbSize * query->L * dbSeq->L) * fpow2((double)-s/m->getBitFactor());
    float seqId = (float)aaIds/(float)(std::min(query->L, dbSeq->L));

    free(Hmatrix);
    free(Ematrix);
    free(Fmatrix);

    result_t result = {std::string(dbSeq->getDbKey()), s, qcov, dbcov, seqId, evalue};
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
