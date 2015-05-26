#include <BaseMatrix.h>
#include "Matcher.h"
#include "Util.h"
#include "smith_waterman_sse2.h"

Matcher::Matcher(int maxSeqLen, BaseMatrix *m) {
    this->m = m;
    this->tinySubMat = NULL;
    setSubstitutionMatrix(m);
    this->maxSeqLen = maxSeqLen;
    aligner = new SmithWaterman(maxSeqLen, m->alphabetSize);
}


void Matcher::setSubstitutionMatrix(BaseMatrix *m){
    this->tinySubMat = new int8_t[m->alphabetSize*m->alphabetSize];
    for (int i = 0; i < m->alphabetSize; i++) {
        for (int j = 0; j < m->alphabetSize; j++) {
            tinySubMat[i*m->alphabetSize + j] = m->subMatrix[i][j];
        }
    }
}

Matcher::~Matcher(){
    delete aligner;
    if(tinySubMat != NULL){
        delete [] tinySubMat;
        tinySubMat = NULL;
    }
}

void Matcher::initQuery(Sequence* query){
    currentQuery = query;
    if(query->getSeqType() == Sequence::HMM_PROFILE)
        aligner->ssw_init(query, query->getAlignmentProfile(), this->m->alphabetSize, 2);
    else
        aligner->ssw_init(query, this->tinySubMat, this->m->alphabetSize, 2);
}


Matcher::result_t Matcher::getSWResult(Sequence* dbSeq, const size_t seqDbSize,
                                        const double evalThr, const unsigned int mode){
    unsigned int qStartPos = 0;
    unsigned int qEndPos = 0;
    unsigned int dbStartPos = 0;
    unsigned int dbEndPos = 0;
    int aaIds = 0;

    // calculation of the score and traceback of the alignment
    int32_t maskLen = currentQuery->L / 2;
    
    // calcuate stop score
    const double qL = static_cast<double>(currentQuery->L);
    const double dbL = static_cast<double>(dbSeq->L);

    // avoid nummerical issues -log(evalThr/(qL*dbL*seqDbSize))
    double datapoints = -log(static_cast<double>(seqDbSize)) - log(qL) - log(dbL) + log(evalThr);
    uint16_t scoreThr = (uint16_t) (m->getBitFactor() * -(datapoints));
    if(evalThr == 0.0)
        scoreThr = 0;
    //std::cout << seqDbSize << " " << 100 << " " << scoreThr << std::endl;
    //std::cout <<datapoints << " " << m->getBitFactor() <<" "<< evalThr << " " << seqDbSize << " " << currentQuery->L << " " << dbSeq->L<< " " << scoreThr << " " << std::endl;
    s_align * alignment = aligner->ssw_align(dbSeq->int_sequence, dbSeq->L, GAP_OPEN, GAP_EXTEND, mode, scoreThr, 0, maskLen);
    // calculation of the coverage and e-value
    float qcov = 0.0;
    float dbcov = 0.0;
    float seqId = 0.0;
    // compute sequence identity
    std::string backtrace;
    if(mode == SCORE_COV_SEQID){
        if(alignment->cigar){
            backtrace.reserve(alignment->cigarLen);
            int32_t targetPos = alignment->dbStartPos1, queryPos = alignment->qStartPos1;
            for (int32_t c = 0; c < alignment->cigarLen; ++c) {
                char letter = SmithWaterman::cigar_int_to_op(alignment->cigar[c]);
                uint32_t length = SmithWaterman::cigar_int_to_len(alignment->cigar[c]);
                for (uint32_t i = 0; i < length; ++i){
                    if (letter == 'M') {
                        if (dbSeq->int_sequence[targetPos] == currentQuery->int_sequence[queryPos]){
                            aaIds++;
                        }
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
            // compute sequence id
            seqId = (float)aaIds/(float)(std::min(currentQuery->L, dbSeq->L)); //TODO
        }
    }
    
    qStartPos = alignment->qStartPos1;
    dbStartPos = alignment->dbStartPos1;
    qEndPos = alignment->qEndPos1;
    dbEndPos = alignment->dbEndPos1;
    // normalize score
    alignment->score1 = alignment->score1 - log2(dbSeq->L);
    if(mode == SCORE_COV || mode == SCORE_COV_SEQID) {
        qcov = (std::min(currentQuery->L, (int) qEndPos) - qStartPos + 1) / (float) currentQuery->L;
        dbcov = (std::min(dbSeq->L, (int) dbEndPos) - dbStartPos + 1) / (float) dbSeq->L;
    }
    // 100 because the score is normalized

    // Karlin-Altschul statistics
    //  E =  qL dL * exp(-S)
    double evalue =  pow (exp(1), ((double)(-(alignment->score1)/(double)m->getBitFactor())));
    evalue *= (qL * seqDbSize * dbSeq->L);
    result_t result(std::string(dbSeq->getDbKey()), alignment->score1, qcov, dbcov, seqId, evalue, qStartPos, qEndPos, dbStartPos, dbEndPos, backtrace);
    delete [] alignment->cigar;
    delete alignment;
    return result;
}


