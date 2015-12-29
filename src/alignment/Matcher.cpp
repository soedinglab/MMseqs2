#include "Matcher.h"
#include <BaseMatrix.h>
#include "Util.h"
#include "smith_waterman_sse2.h"

Matcher::Matcher(int maxSeqLen, BaseMatrix *m, size_t dbLen, size_t dbSize, bool aaBiasCorrection){
    this->m = m;
    this->tinySubMat = NULL;
    setSubstitutionMatrix(m);
    this->maxSeqLen = maxSeqLen;
    aligner = new SmithWaterman(maxSeqLen, m->alphabetSize, aaBiasCorrection);
    kmnByLen = new double[maxSeqLen];
    BlastScoreUtils::BlastStat stats = BlastScoreUtils::getAltschulStatsForMatrix(m->getMatrixName(), GAP_OPEN, GAP_EXTEND);
    for(size_t len = 0; len < maxSeqLen; len++){
        kmnByLen[len] = BlastScoreUtils::computeKmn(len, stats.K, stats.lambda, stats.alpha, stats.beta, dbLen, dbSize);
    }
    double logK = log( stats.K);
    this->lambda = stats.lambda;
    this->lambdaLog2 =  stats.lambda / log(2.0);
    this->logKLog2 = logK / log(2.0);
    //std::cout << "lambda=" << lambdaLog2 << " logKLog2=" << logKLog2 << std::endl;
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
    delete[] kmnByLen;
}

void Matcher::initQuery(Sequence* query){
    currentQuery = query;
    if(query->getSeqType() == Sequence::HMM_PROFILE){
        aligner->ssw_init(query, query->getAlignmentProfile(), this->m, this->m->alphabetSize, 2);
    }else{
        aligner->ssw_init(query, this->tinySubMat, this->m, this->m->alphabetSize, 2);
    }
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
            seqId = (float)aaIds/(float)(std::min(currentQuery->L, dbSeq->L));
        }
    }

    qStartPos = alignment->qStartPos1;
    dbStartPos = alignment->dbStartPos1;
    qEndPos = alignment->qEndPos1;
    dbEndPos = alignment->dbEndPos1;
    // normalize score
//    alignment->score1 = alignment->score1 - log2(dbSeq->L);
    if(mode == SCORE_COV || mode == SCORE_COV_SEQID) {
        qcov  = computeCov(qStartPos, qEndPos, currentQuery->L);
        dbcov = computeCov(dbStartPos, dbEndPos, dbSeq->L);
    }

    // statistics
    //  E =  qL dL * 2^(-S)
    double evalue = BlastScoreUtils::computeEvalue(alignment->score1, kmnByLen[currentQuery->L], this->lambda);
    int bitScore =(short) (BlastScoreUtils::computeBitScore(alignment->score1, lambdaLog2, logKLog2)+0.5);
    size_t alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, dbStartPos, dbEndPos);

    //blast stat
//    double lambda= 0.267;
//    double K= 0.041;
//    double Kmn=(qL * seqDbSize * dbSeq->L);
//    double evalue = Kmn * exp(-(alignment->score1 * lambda));
    result_t result(dbSeq->getDbKey(), bitScore, qcov, dbcov, seqId, evalue, alnLength, qStartPos, qEndPos, currentQuery->L, dbStartPos, dbEndPos, dbSeq->L, backtrace);
    delete [] alignment->cigar;
    delete alignment;
    return result;
}


float Matcher::computeCov(unsigned int startPos, unsigned int endPos, unsigned int len) {
    return (std::min(len, endPos) - startPos + 1) / (float) len;
}

std::vector<Matcher::result_t> Matcher::readAlignmentResults(char *data) {
    std::vector<Matcher::result_t> ret;
    while(*data != '\0'){
        char * entry[255];
        char key[255];
        size_t elmCount = Util::getWordsOfLine(data, entry, 255 );
        ptrdiff_t keySize =  (entry[1] - data);
        strncpy(key, data, keySize);
        key[keySize] = '\0';

        unsigned int targetId = (unsigned int) strtoul(key, NULL, 10);
        double score = strtod(entry[1],NULL);
        double seqId = strtod(entry[2],NULL);
        double eval = strtod(entry[3],NULL);


        size_t qStart = strtoull(entry[4],NULL,0);
        size_t qEnd = strtoull(entry[5],NULL,0);
        size_t qLen = strtoull(entry[6],NULL,0);
        size_t dbStart = strtoull(entry[7],NULL,0);
        size_t dbEnd = strtoull(entry[8],NULL,0);
        size_t dbLen = strtoull(entry[9],NULL,0);
        double qCov = Matcher::computeCov(qStart, qEnd, qLen);
        double dbCov = Matcher::computeCov(dbStart, dbEnd, dbLen);
        size_t alnLength = Matcher::computeAlnLength(qStart, qEnd, dbStart, dbEnd);

        Matcher::result_t result(targetId, score, qCov, dbCov, seqId, eval, alnLength, qStart, qEnd, qLen, dbStart, dbEnd,
                                 dbLen, "");
        ret.push_back(result);

        data = Util::skipLine(data);
    }
    return ret;
}

size_t Matcher::computeAlnLength(size_t qStart, size_t qEnd, size_t dbStart, size_t dbEnd) {
    return std::max(qEnd - qStart, dbEnd - dbStart);
}
