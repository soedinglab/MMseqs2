#include "Matcher.h"
#include "Util.h"

Matcher::Matcher(BaseMatrix* m, int maxSeqLen){
    
    this->m = m;
    this->maxSeqLen = maxSeqLen;
    this->tinySubMat = new int8_t[m->alphabetSize*m->alphabetSize];
    for (size_t i = 0; i < m->alphabetSize; i++) {
        for (size_t j = 0; j < m->alphabetSize; j++) {
            tinySubMat[i*m->alphabetSize + j] = m->subMatrix[i][j];
        }
    }
    aligner = new SmithWaterman(maxSeqLen, m->alphabetSize);
}

Matcher::~Matcher(){
    delete aligner;
    delete [] tinySubMat;
}

void Matcher::initQuery(Sequence* query){
    currentQuery = query;
    aligner->ssw_init(query, this->tinySubMat, this->m->alphabetSize, 2);
}

Matcher::result_t Matcher::getSWResult(Sequence* dbSeq,const size_t seqDbSize,const double evalThr){
    
    
    unsigned short qStartPos = 0;
    unsigned short qEndPos = 0;
    unsigned short dbStartPos = 0;
    unsigned short dbEndPos = 0;
    int aaIds = 0;
    int overflow_warning = 0;
    
    // calculation of the score and traceback of the alignment
    int32_t maskLen = currentQuery->L / 2;
    
    // calcuate stop score
    const double qL = static_cast<double>(currentQuery->L);
    const double dbL = static_cast<double>(dbSeq->L);
    // avoid nummerical issues -log(evalThr/(qL*dbL*seqDbSize))
    double datapoints = -log(static_cast<double>(seqDbSize)) - log(qL) - log(dbL) + log(evalThr);
    uint16_t scoreThr = (uint16_t) (m->getBitFactor() * -(datapoints));
    //std::cout <<datapoints << " " << m->getBitFactor() <<" "<< evalThr << " " << seqDbSize << " " << currentQuery->L << " " << dbSeq->L<< " " << scoreThr << " " << std::endl;
    s_align * alignment = aligner->ssw_align(dbSeq->int_sequence, dbSeq->L, GAP_OPEN, GAP_EXTEND, 2, scoreThr, 0, maskLen);
    // calculation of the coverage and e-value
    float qcov;
    float dbcov;
    float seqId;
    // compute sequence identity
    if(alignment->cigar){
        int32_t targetPos = alignment->dbStartPos1, queryPos = alignment->qStartPos1;
        for (int32_t c = 0; c < alignment->cigarLen; ++c) {
            char letter = SmithWaterman::cigar_int_to_op(alignment->cigar[c]);
            uint32_t length = SmithWaterman::cigar_int_to_len(alignment->cigar[c]);
            for (int i = 0; i < length; ++i){
                if (letter == 'M') {
                    if (dbSeq->int_sequence[targetPos] == currentQuery->int_sequence[queryPos]){
                        aaIds++;
                    }
                    ++queryPos;
                    ++targetPos;
                } else {
                    if (letter == 'I') ++queryPos;
                    else ++targetPos;
                }
            }
        }
    }
    
    qStartPos = alignment->qStartPos1;
    dbStartPos = alignment->dbStartPos1;
    qEndPos = alignment->qEndPos1;
    dbEndPos = alignment->dbEndPos1;
    
    if (overflow_warning == 0){
        qcov = (std::min(currentQuery->L, (int) qEndPos) - qStartPos + 1)/ (float)currentQuery->L;
        dbcov = (std::min(dbSeq->L, (int) dbEndPos) - dbStartPos + 1)/(float)dbSeq->L;
        seqId = (float)aaIds/(float)(std::min(currentQuery->L, dbSeq->L)); //TODO
    }
    else{
        // there was an unsigned short range overflow during the alignment
        // we cannot say anything reliable about sequence coverage and identity
        qcov = 1.0;
        dbcov = 1.0;
        seqId = 1.0;
    }
    double evalue = ( static_cast<double>(qL * dbL)) * pow (2.71828, ((double)(-alignment->score1)/(double)m->getBitFactor())); // fpow2((double)-s/m->getBitFactor());
    evalue = evalue * (double)(seqDbSize);
    result_t result = {std::string(dbSeq->getDbKey()), alignment->score1, qcov, dbcov, seqId, evalue};
    delete [] alignment->cigar;
    delete alignment;
    return result;
}


