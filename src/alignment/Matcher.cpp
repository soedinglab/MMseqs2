#include <iomanip>
#include <itoa.h>
#include "Matcher.h"
#include "Util.h"
#include "Parameters.h"
#include "StripedSmithWaterman.h"


Matcher::Matcher(int querySeqType, int targetSeqType, int maxSeqLen, BaseMatrix *m, EvalueComputation * evaluer,
                 bool aaBiasCorrection, float aaBiasCorrectionScale, int gapOpen, int gapExtend, float correlationScoreWeight, int zdrop)
                 : gapOpen(gapOpen), gapExtend(gapExtend), correlationScoreWeight(correlationScoreWeight), m(m), evaluer(evaluer), tinySubMat(NULL)  {
    setSubstitutionMatrix(m);

    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        nuclaligner = new BandedNucleotideAligner(m, maxSeqLen, gapOpen, gapExtend, zdrop);
        aligner = NULL;
    } else {
        nuclaligner = NULL;
        aligner = new SmithWaterman(maxSeqLen, m->alphabetSize, aaBiasCorrection,
                                    aaBiasCorrectionScale, targetSeqType);
    }
    //std::cout << "lambda=" << lambdaLog2 << " logKLog2=" << logKLog2 << std::endl;
}


void Matcher::setSubstitutionMatrix(BaseMatrix *m){
    tinySubMat = new int8_t[m->alphabetSize * m->alphabetSize];
    for (int i = 0; i < m->alphabetSize; i++) {
        for (int j = 0; j < m->alphabetSize; j++) {
            tinySubMat[i*m->alphabetSize + j] = m->subMatrix[i][j];
        }
    }
}

Matcher::~Matcher(){
    if(aligner != NULL){
        delete aligner;
    }
    if(nuclaligner != NULL){
        delete nuclaligner;
    }
    if(tinySubMat != NULL){
        delete [] tinySubMat;
        tinySubMat = NULL;
    }
}

void Matcher::initQuery(Sequence* query){
    currentQuery = query;
    if(Parameters::isEqualDbtype(query->getSequenceType(), Parameters::DBTYPE_NUCLEOTIDES)){
        nuclaligner->initQuery(query);
    }else if(Parameters::isEqualDbtype(query->getSeqType(), Parameters::DBTYPE_HMM_PROFILE)){
        aligner->ssw_init(query, query->getAlignmentProfile(), this->m);
    }else{
        aligner->ssw_init(query, this->tinySubMat, this->m);
    }
}

Matcher::result_t Matcher::getSWResult(Sequence* dbSeq, const int diagonal, bool isReverse, const int covMode, const float covThr,
                                       const double evalThr, unsigned int alignmentMode, unsigned int seqIdMode, bool isIdentity,
                                       bool wrappedScoring){

    // calculation of the score and traceback of the alignment
    int32_t maskLen = currentQuery->L / 2;
    int origQueryLen = wrappedScoring? currentQuery->L / 2 : currentQuery->L ;

    s_align alignment;
    // compute sequence identity
    std::string backtrace;
    if(Parameters::isEqualDbtype(dbSeq->getSequenceType(), Parameters::DBTYPE_NUCLEOTIDES)){
        if(diagonal==INT_MAX){
            Debug(Debug::ERROR) << "Query sequence " << currentQuery->getDbKey() << " has a result with no diagonal information. Please check your database.\n";
            EXIT(EXIT_FAILURE);
        }
        alignment = nuclaligner->align(dbSeq, diagonal, isReverse, backtrace, evaluer, wrappedScoring);
        alignmentMode = Matcher::SCORE_COV_SEQID;
    } else {
        if (isIdentity == false) {
            alignment = aligner->ssw_align(dbSeq->numSequence, dbSeq->numConsensusSequence,
                                           dbSeq->getAlignmentProfile(), dbSeq->L, backtrace,
                                           gapOpen, gapExtend, alignmentMode, evalThr, evaluer, covMode,
                                           covThr, correlationScoreWeight, maskLen, dbSeq->getId());
        } else {
            alignment = aligner->scoreIdentical(dbSeq->numSequence, dbSeq->L, evaluer, alignmentMode, backtrace);
        }
    }

    // calculation of the coverage and e-value
    float qcov = 0.0;
    float dbcov = 0.0;
    float seqId = 0.0;

    const unsigned int qStartPos = alignment.qStartPos1;
    const unsigned int dbStartPos = alignment.dbStartPos1;
    const unsigned int qEndPos = alignment.qEndPos1;
    const unsigned int dbEndPos = alignment.dbEndPos1;
    // normalize score
//    alignment->score1 = alignment->score1 - log2(dbSeq->L);
    if(alignmentMode == Matcher::SCORE_COV || alignmentMode == Matcher::SCORE_COV_SEQID) {
        qcov  = alignment.qCov;
        dbcov = alignment.tCov;
    }

    unsigned int alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, dbStartPos, dbEndPos);
    // try to estimate sequence id
    if(alignmentMode == Matcher::SCORE_COV_SEQID){
        // compute sequence id
        if(alignment.cigar){
            // OVERWRITE alnLength with gapped value
            alnLength = backtrace.size();
        }
        seqId = Util::computeSeqId(seqIdMode, alignment.identicalAACnt, origQueryLen, dbSeq->L, alnLength);
    }else if( alignmentMode == Matcher::SCORE_COV){
        // "20%   30%   40%   50%   60%   70%   80%   90%   99%"
        // "0.52  1.12  1.73  2.33  2.93  3.53  4.14  4.74  5.28"
        unsigned int qAlnLen = std::max(qEndPos - qStartPos, static_cast<unsigned int>(1));
        unsigned int dbAlnLen = std::max(dbEndPos - dbStartPos, static_cast<unsigned int>(1));
        //seqId = (alignment.score1 / static_cast<float>(std::max(qAlnLength, dbAlnLength)))  * 0.1656 + 0.1141;
        seqId = estimateSeqIdByScorePerCol(alignment.score1, qAlnLen, dbAlnLen);
    }else if ( alignmentMode == Matcher::SCORE_ONLY){
        unsigned int qAlnLen = std::max(qEndPos, static_cast<unsigned int>(1));
        unsigned int dbAlnLen = std::max(dbEndPos, static_cast<unsigned int>(1));
        //seqId = (alignment.score1 / static_cast<float>(std::max(dbAlnLen, qAlnLen)))  * 0.1656 + 0.1141;
        seqId = estimateSeqIdByScorePerCol(alignment.score1, qAlnLen, dbAlnLen);
    }

    //  E =  qL dL * exp^(-S/lambda)
    double evalue = alignment.evalue;
    int bitScore = static_cast<int>(evaluer->computeBitScore(alignment.score1)+0.5);

    result_t result;
    if(isReverse){
        result = result_t(dbSeq->getDbKey(), bitScore, qcov, dbcov, seqId, evalue, alnLength, qStartPos, qEndPos, origQueryLen, dbEndPos, dbStartPos, dbSeq->L, backtrace);
    }else{
        result = result_t(dbSeq->getDbKey(), bitScore, qcov, dbcov, seqId, evalue, alnLength, qStartPos, qEndPos, origQueryLen, dbStartPos, dbEndPos, dbSeq->L, backtrace);
    }


    delete [] alignment.cigar;
    return result;
}


void Matcher::readAlignmentResults(std::vector<result_t> &result, char *data, bool readCompressed) {
    if(data == NULL) {
        return;
    }

    while(*data != '\0'){
        result.emplace_back(parseAlignmentRecord(data, readCompressed));
        data = Util::skipLine(data);
    }
}

int Matcher::computeAlnLength(int qStart, int qEnd, int dbStart, int dbEnd) {
    return std::max(abs(qEnd - qStart), abs(dbEnd - dbStart)) + 1;
}

float Matcher::estimateSeqIdByScorePerCol(uint16_t score, unsigned int qLen, unsigned int tLen) {
    float estimatedSeqId = (score / static_cast<float>(std::max(qLen, tLen))) * 0.1656 + 0.1141;
    estimatedSeqId = std::min(estimatedSeqId, 1.0f);
    return std::max(0.0f, estimatedSeqId);
}

std::string Matcher::compressAlignment(const std::string& bt) {
    std::string ret;
    char state = 'M';
    size_t counter = 0;
    for (size_t i = 0; i < bt.size(); ++i) {
        if (bt[i] != state) {
            // we could leave this out if counter == 1
            // to save a few byte (~5% of total cigar strings)
            ret.append(SSTR(counter));
            ret.push_back(state);
            state = bt[i];
            counter = 1;
        } else {
            counter++;
        }
    }
    ret.append(SSTR(counter));
    ret.push_back(state);
    return ret;
}

std::string Matcher::uncompressAlignment(const std::string &cbt) {
    std::string bt;
    bt.reserve(cbt.size());
    size_t count = 0;
    for (size_t i = 0; i < cbt.size(); ++i) {
        char c = cbt[i];
        if (c >= '0' && c <= '9') {
            count = count * 10 + c - '0';
        } else {
            bt.append(count == 0 ? 1 : count, c);
            count = 0;
        }
    }
    return bt;
}

Matcher::result_t Matcher::parseAlignmentRecord(const char *data, bool readCompressed) {
    const char *entry[255];
    size_t columns = Util::getWordsOfLine(data, entry, 255);
    if (columns < ALN_RES_WITHOUT_BT_COL_CNT) {
        Debug(Debug::ERROR) << "Invalid alignment result record.\n";
        EXIT(EXIT_FAILURE);
    }

    char key[255];
    ptrdiff_t keySize =  (entry[1] - data);
    strncpy(key, data, keySize);
    key[keySize] = '\0';

    unsigned int targetId = Util::fast_atoi<unsigned int>(key);
    int score = Util::fast_atoi<int>(entry[1]);
    double seqId = strtod(entry[2],NULL);
    double eval = strtod(entry[3],NULL);

    int qStart =  Util::fast_atoi<int>(entry[4]);
    int qEnd = Util::fast_atoi<int>(entry[5]);
    int qLen = Util::fast_atoi<int>(entry[6]);
    int dbStart = Util::fast_atoi<int>(entry[7]);
    int dbEnd = Util::fast_atoi<int>(entry[8]);
    int dbLen = Util::fast_atoi<int>(entry[9]);
    int adjustQstart = (qStart==-1)? 0 : qStart;
    int adjustDBstart = (dbStart==-1)? 0 : dbStart;
    double qCov = SmithWaterman::computeCov(adjustQstart, qEnd, qLen);
    double dbCov = SmithWaterman::computeCov(adjustDBstart, dbEnd, dbLen);
    size_t alnLength = Matcher::computeAlnLength(adjustQstart, qEnd, adjustDBstart, dbEnd);

    switch(columns) {
        // 10 no backtrace
        case ALN_RES_WITHOUT_BT_COL_CNT:
            return Matcher::result_t(targetId, score, qCov, dbCov, seqId, eval,
                                 alnLength, qStart, qEnd, qLen, dbStart, dbEnd, dbLen, -1, -1, -1, -1, "");
        // 11 with backtrace
        case ALN_RES_WITH_BT_COL_CNT:
            if (readCompressed) {
                return Matcher::result_t(targetId, score, qCov, dbCov, seqId, eval,
                                         alnLength, qStart, qEnd, qLen, dbStart, dbEnd,
                                         dbLen, -1, -1, -1, -1, std::string(entry[10], entry[11] - entry[10]));
            } else {
                return Matcher::result_t(targetId, score, qCov, dbCov, seqId, eval,
                                         alnLength, qStart, qEnd, qLen, dbStart, dbEnd,
                                         dbLen, -1, -1, -1, -1,
                                         uncompressAlignment(std::string(entry[10], entry[11] - entry[10])));
            }
        // 12 without backtrace but qOrfStart dbOrfStart
        case ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT:
            return Matcher::result_t(targetId, score, qCov, dbCov, seqId, eval,
                                     alnLength, qStart, qEnd, qLen, dbStart, dbEnd,
                                     dbLen, Util::fast_atoi<int>(entry[10]), Util::fast_atoi<int>(entry[11]),
                                     Util::fast_atoi<int>(entry[12]), Util::fast_atoi<int>(entry[13]), "");
        // 13 without backtrace but qOrfStart dbOrfStart
        case ALN_RES_WITH_ORF_AND_BT_COL_CNT:
            if (readCompressed) {
                return Matcher::result_t(targetId, score, qCov, dbCov, seqId, eval,
                                         alnLength, qStart, qEnd, qLen, dbStart, dbEnd,
                                         dbLen, Util::fast_atoi<int>(entry[10]), Util::fast_atoi<int>(entry[11]),
                                         Util::fast_atoi<int>(entry[12]), Util::fast_atoi<int>(entry[13]),
                                         std::string(entry[14], entry[15] - entry[14]));
            } else {
                return Matcher::result_t(targetId, score, qCov, dbCov, seqId, eval,
                                         alnLength, qStart, qEnd, qLen, dbStart, dbEnd,
                                         dbLen, Util::fast_atoi<int>(entry[10]), Util::fast_atoi<int>(entry[11]),
                                         Util::fast_atoi<int>(entry[12]), Util::fast_atoi<int>(entry[13]),
                                         uncompressAlignment(std::string(entry[14], entry[15] - entry[14])));
            }
        default:
            Debug(Debug::ERROR) << "Invalid column count in alignment.\n";
            EXIT(EXIT_FAILURE);
    }
}


size_t Matcher::resultToBuffer(char * buff1, const result_t &result, bool addBacktrace, bool compress, bool addOrfPosition) {
    char * basePos = buff1;
    char * tmpBuff = Itoa::u32toa_sse2((uint32_t) result.dbKey, buff1);
    *(tmpBuff-1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.score, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = Util::fastSeqIdToBuffer(result.seqId, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff += snprintf(tmpBuff, 32, "%.3E", result.eval);
    tmpBuff++;
    *(tmpBuff-1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.qStartPos, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.qEndPos, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.qLen, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.dbStartPos, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.dbEndPos, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.dbLen, tmpBuff);
    if(addOrfPosition){
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(result.queryOrfStartPos, tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(result.queryOrfEndPos, tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(result.dbOrfStartPos, tmpBuff);
        *(tmpBuff-1) = '\t';
        tmpBuff = Itoa::i32toa_sse2(result.dbOrfEndPos, tmpBuff);
    }
    if(addBacktrace == true){
        if(compress){
            *(tmpBuff-1) = '\t';
            std::string compressedCigar = Matcher::compressAlignment(result.backtrace);
            tmpBuff = strncpy(tmpBuff, compressedCigar.c_str(), compressedCigar.length());
            tmpBuff += compressedCigar.length()+1;
        }else{
            *(tmpBuff-1) = '\t';
            tmpBuff = strncpy(tmpBuff, result.backtrace.c_str(), result.backtrace.length());
            tmpBuff += result.backtrace.length()+1;
        }
    }
    *(tmpBuff-1) = '\n';
    *(tmpBuff) = '\0';
    return tmpBuff - basePos;
}

void Matcher::updateResultByRescoringBacktrace(const char *querySeq, const char *targetSeq, const char **subMat, EvalueComputation &evaluer,
                                                int gapOpen, int gapExtend, result_t &result) {
    int maxScore = 0;
    int maxBtEndPos = 0;
    int maxBtStartPos = 0;
    int maxQueryEndPos = 0;
    int maxQueryStartPos = 0;
    int maxTargetStartPos = 0;
    int maxTargetEndPos = 0;
    int minPos = -1;
    int minQueryPos = result.qStartPos-1;
    int minTargetPos = result.dbStartPos-1;
    int score = 0;
    int identicalAAs = 0;
    int maxIdAaCnt = 0;
    int queryPos = result.qStartPos;
    int targetPos = result.dbStartPos;
    bool isGapOpen = false;

    for(unsigned int pos = 0; pos < result.backtrace.size(); pos++){
        char letter = result.backtrace[pos];
        int curr;
        if (letter == 'M') {
            curr = subMat[static_cast<int>(querySeq[queryPos])][static_cast<int>(targetSeq[targetPos])];
            identicalAAs += (querySeq[queryPos] == targetSeq[targetPos]);
            isGapOpen = false;
        } else {
            curr = (isGapOpen) ? -gapExtend : -gapOpen;
            isGapOpen = (isGapOpen == false) ? true : isGapOpen;
        }
        score = curr + score;
        // minimum
        const bool isMinScore = (score <= 0);
        score = (isMinScore) ? 0 : score;
        identicalAAs = (isMinScore) ? 0 : identicalAAs;
        if(isMinScore){
            minPos = pos;
            minQueryPos = (letter == 'D') ? queryPos - 1 : queryPos;
            minTargetPos = (letter == 'I') ? targetPos - 1 : targetPos;
        }
        // new max
        const bool isNewMaxScore = (score > maxScore);
        if(isNewMaxScore){
            maxBtEndPos = pos;
            maxQueryEndPos = queryPos;
            maxTargetEndPos = targetPos;
            maxBtStartPos = minPos + 1;
            maxQueryStartPos = minQueryPos + 1;
            maxTargetStartPos = minTargetPos + 1;
            maxScore = score;
            maxIdAaCnt = identicalAAs;
        }
        queryPos += (letter == 'M' || letter == 'I') ? 1 : 0;
        targetPos += (letter == 'M' || letter == 'D') ? 1 : 0;
    }

    result.qStartPos = maxQueryStartPos;
    result.qEndPos = maxQueryEndPos;
    result.dbStartPos = maxTargetStartPos;
    result.dbEndPos = maxTargetEndPos;
    double bitScore = evaluer.computeBitScore(maxScore);
    double evalue = evaluer.computeEvalue(maxScore, result.qLen);
    result.score = bitScore;
    result.eval = evalue;
    result.alnLength = (maxBtEndPos - maxBtStartPos) + 1;
    result.seqId = static_cast<float>(maxIdAaCnt) / static_cast<float>(result.alnLength);
    result.backtrace = result.backtrace.substr(maxBtStartPos, result.alnLength);

}


