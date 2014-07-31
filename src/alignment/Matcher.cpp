#include "Matcher.h"
#include "../commons/Util.h"

Matcher::Matcher(BaseMatrix* m, int maxSeqLen){

    this->m = m;
    this->maxSeqLen = maxSeqLen;

    // memory for the sequence profile is allocated only once
    // = 1.6 MB
    this->queryProfileWord = (unsigned short*) Util::mem_align (16, m->alphabetSize * maxSeqLen * sizeof(unsigned short));

    // workspace memory for the alignment calculation (see Farrar code)
    this->workspace_memory  = (void *) Util::mem_align(16, 2 * maxSeqLen * sizeof(__m128i) + 256);
    this->workspace = (void *) ((((size_t) workspace_memory) + 255) & (~0xff));

    maxAllocatedLen = 1000;
    
    H_workspace = Util::mem_align(16,(maxAllocatedLen + 7)/8 * maxAllocatedLen * sizeof(__m128i));
    E_workspace = Util::mem_align(16,(maxAllocatedLen + 7)/8 * maxAllocatedLen * sizeof(__m128i));
    F_workspace = Util::mem_align(16,(maxAllocatedLen + 7)/8 * maxAllocatedLen * sizeof(__m128i));
}

Matcher::~Matcher(){
    free(this->queryProfileWord);
    free(this->workspace_memory);
    free(this->H_workspace);
    free(this->E_workspace);
    free(this->F_workspace);
}

Matcher::result_t Matcher::getSWResult(Sequence* query, Sequence* dbSeq, int seqDbSize){

    unsigned short gap_open = 10;
    unsigned short gap_extend = 1;

    calcQueryProfileWord(query);

    void* Hmatrix;
    void* Ematrix;
    void* Fmatrix;

    if ((query->L > maxAllocatedLen) || (dbSeq->L > maxAllocatedLen)){
        // allocate new memory
        Hmatrix = Util::mem_align(16,(query->L + 7)/8 * dbSeq->L * sizeof(__m128i));   // 2.5GB für 36805*36805 (Q3ASY8_CHLCH)
        Ematrix = Util::mem_align(16,(query->L + 7)/8 * dbSeq->L * sizeof(__m128i));
        Fmatrix = Util::mem_align(16,(query->L + 7)/8 * dbSeq->L * sizeof(__m128i));
    }
    else {
        // use the already allocated memory
        Hmatrix = H_workspace;
        Ematrix = E_workspace;
        Fmatrix = F_workspace;
    }

    unsigned short qStartPos = 0;
    unsigned short qEndPos = 0;
    unsigned short dbStartPos = 0;
    unsigned short dbEndPos = 0;
    int aaIds = 0;
    int overflow_warning = 0;

    // calculation of the score
    int s = smith_waterman_sse2_word(query->getDbKey(), query->int_sequence, this->queryProfileWord, query->L, 
            dbSeq->getDbKey(), dbSeq->int_sequence, dbSeq->L, 
            gap_open, gap_extend, 
            workspace, 
            Hmatrix, Ematrix, Fmatrix, 
            &qEndPos, &dbEndPos);

    // traceback of the alignment
    traceback_word((short*) Hmatrix, (short*) Ematrix, (short*) Fmatrix, 
            query, dbSeq,
            this->queryProfileWord,
            qEndPos, dbEndPos, 
            gap_open, gap_extend, 
            &qStartPos, &dbStartPos, &aaIds, &overflow_warning);

    // calculation of the coverage and e-value
    float qcov;
    float dbcov;
    float seqId;

    if (overflow_warning == 0){
        qcov = (std::min(query->L, (int) qEndPos) - qStartPos + 1)/ (float)query->L;
        dbcov = (std::min(dbSeq->L, (int) dbEndPos) - dbStartPos + 1)/(float)dbSeq->L;
        seqId = (float)aaIds/(float)(std::min(query->L, dbSeq->L));
    }
    else{
        // there was an unsigned short range overflow during the alignment
        // we cannot say anything reliable about sequence coverage and identity
        qcov = 1.0;
        dbcov = 1.0;
        seqId = 1.0;
    }

    double evalue = ((double) (query->L * dbSeq->L)) * pow (2.0, ((double)(-s)/(double)m->getBitFactor())); // fpow2((double)-s/m->getBitFactor());
    evalue = evalue * (double)(seqDbSize);

    if ((query->L > maxAllocatedLen) || (dbSeq->L > maxAllocatedLen)){
        free(Hmatrix);
        free(Ematrix);
        free(Fmatrix);
    }

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
