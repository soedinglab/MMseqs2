#ifndef MATCHER_H
#define MATCHER_H

//
// Written by Martin Steinegger & Maria Hauser, mhauser@genzentrum.lmu.de
//
// Calls SSE2 parallelized calculation of Smith-Waterman alignment and non-parallelized traceback afterwards.
//

#include <cfloat>
#include <algorithm>
#include <vector>
#include "itoa.h"

#include "Sequence.h"
#include "BaseMatrix.h"
#include "StripedSmithWaterman.h"
#include "EvalueComputation.h"
#include "BandedNucleotideAligner.h"

class Matcher{

public:
    static const unsigned int SCORE_ONLY = 0;
    static const unsigned int SCORE_COV = 1;
    static const unsigned int SCORE_COV_SEQID = 2;
    const static int ALN_RES_WITHOUT_BT_COL_CNT = 10;
    const static int ALN_RES_WITH_BT_COL_CNT = 11;
    const static int ALN_RES_WITH_ORF_POS_WITHOUT_BT_COL_CNT = 14;
    const static int ALN_RES_WITH_ORF_AND_BT_COL_CNT = 15;

    struct result_t {
        unsigned int dbKey;
        int score;
        float qcov;
        float dbcov;
        float seqId;
        double eval;
        unsigned int alnLength;
        int qStartPos;
        int qEndPos;
        unsigned int qLen;
        int dbStartPos;
        int dbEndPos;
        unsigned int dbLen;
        int queryOrfStartPos;
        int queryOrfEndPos;
        int dbOrfStartPos;
        int dbOrfEndPos;
        std::string backtrace;
        result_t(unsigned int dbkey,int score,
                 float qcov, float dbcov,
                 float seqId, double eval,
                 unsigned int alnLength,
                 int qStartPos,
                 int qEndPos,
                 unsigned int qLen,
                 int dbStartPos,
                 int dbEndPos,
                 unsigned int dbLen,
                 int queryOrfStartPos,
                 int queryOrfEndPos,
                 int dbOrfStartPos,
                 int dbOrfEndPos,
                 std::string backtrace) : dbKey(dbkey), score(score), qcov(qcov),
                                          dbcov(dbcov), seqId(seqId), eval(eval), alnLength(alnLength),
                                          qStartPos(qStartPos), qEndPos(qEndPos), qLen(qLen),
                                          dbStartPos(dbStartPos), dbEndPos(dbEndPos), dbLen(dbLen),
                                          queryOrfStartPos(queryOrfStartPos), queryOrfEndPos(queryOrfEndPos),
                                          dbOrfStartPos(dbOrfStartPos), dbOrfEndPos(dbOrfEndPos),
                                          backtrace(backtrace) {};

        result_t(unsigned int dbkey,int score,
                 float qcov, float dbcov,
                 float seqId, double eval,
                 unsigned int alnLength,
                 int qStartPos,
                 int qEndPos,
                 unsigned int qLen,
                 int dbStartPos,
                 int dbEndPos,
                 unsigned int dbLen,
                 std::string backtrace) : dbKey(dbkey), score(score), qcov(qcov),
                                          dbcov(dbcov), seqId(seqId), eval(eval), alnLength(alnLength),
                                          qStartPos(qStartPos), qEndPos(qEndPos), qLen(qLen),
                                          dbStartPos(dbStartPos), dbEndPos(dbEndPos), dbLen(dbLen),
                                          queryOrfStartPos(-1), queryOrfEndPos(-1),
                                          dbOrfStartPos(-1), dbOrfEndPos(-1),
                                          backtrace(backtrace) {};

        result_t(){};

        static void swapResult(result_t & res, EvalueComputation &evaluer, bool hasBacktrace){
            double rawScore = evaluer.computeRawScoreFromBitScore(res.score);
            res.eval = evaluer.computeEvalue(rawScore, res.dbLen);

            unsigned int qstart = res.qStartPos;
            unsigned int qend = res.qEndPos;
            unsigned int qLen = res.qLen;
            res.qStartPos = res.dbStartPos;
            res.qEndPos = res.dbEndPos;
            res.qLen = res.dbLen;
            res.dbStartPos = qstart;
            res.dbEndPos = qend;
            res.dbLen = qLen;
            if (hasBacktrace) {
                for (size_t j = 0; j < res.backtrace.size(); j++) {
                    if (res.backtrace.at(j) == 'I') {
                        res.backtrace.at(j) = 'D';
                    } else if (res.backtrace.at(j) == 'D') {
                        res.backtrace.at(j) = 'I';
                    }
                }
            }
        }

        static void protein2nucl(std::string & backtrace, std::string &newBacktrace) {
            char buffer[256];
            for (size_t pos = 0; pos < backtrace.size(); pos++) {
                int cnt =0;
                if (isdigit(backtrace[pos])){
                    cnt += Util::fast_atoi<int>(backtrace.c_str()+pos);
                    while (isdigit(backtrace[pos])){
                        pos++;
                    }
                }
                bool update = false;
                switch (backtrace[pos]) {
                    case 'M':
                    case 'D':
                    case 'I':
                        update = true;
                        break;
                }
                if (update) {
                    char *buffNext = Itoa::i32toa_sse2(cnt*3, buffer);
                    size_t len = buffNext - buffer;
                    newBacktrace.append(buffer, len - 1);
                    newBacktrace.push_back(backtrace[pos]);
                }
            }
        }
    };

    Matcher(int querySeqType, int maxSeqLen, BaseMatrix *m,
            EvalueComputation * evaluer, bool aaBiasCorrection,
            int gapOpen, int gapExtend, int zdrop = 40);

    ~Matcher();

    // run SSE2 parallelized Smith-Waterman alignment calculation and traceback
    result_t getSWResult(Sequence* dbSeq, const int diagonal, bool isReverse, const int covMode, const float covThr, const double evalThr,
                         unsigned int alignmentMode, unsigned int seqIdMode, bool isIdentical, bool wrappedScoring=false);

    // need for sorting the results
    static bool compareHits(const result_t &first, const result_t &second) {
        if (first.eval != second.eval) {
            return first.eval < second.eval;
        }
        if (first.score != second.score) {
            return first.score > second.score;
        }
        if (first.dbLen != second.dbLen) {
            return first.dbLen < second.dbLen;
        }
        return first.dbKey < second.dbKey;
    }

    static bool compareHitByPos(const result_t &first, const result_t &second) {
        int firstQStartPos  = std::min(first.qStartPos, first.qEndPos);
        int secondQStartPos = std::min(second.qStartPos, second.qEndPos);
        return firstQStartPos < secondQStartPos;
    }

    // need for sorting the results
    static bool compareHitsByPosAndStrand(const result_t &first, const result_t &second) {
        if (first.dbKey != second.dbKey) {
            return first.dbKey < second.dbKey;
        }
        bool qFirstRev = (first.qStartPos > first.qEndPos);
        bool qSecondRev = (second.qStartPos > second.qEndPos);
        if (qSecondRev < qFirstRev)
            return false;
        if (qFirstRev < qSecondRev)
            return true;
        bool dbFirstRev = (first.dbStartPos > first.dbEndPos);
        bool dbSecondRev = (second.dbStartPos > second.dbEndPos);
        if (dbSecondRev < dbFirstRev)
            return false;
        if (dbFirstRev < dbSecondRev)
            return true;
        int firstQStartPos  = std::min(first.qStartPos, first.qEndPos);
        int secondQStartPos = std::min(second.qStartPos, second.qEndPos);
        int firstDbStart    = std::min(first.dbStartPos, first.dbEndPos);
        int secondDbStart   = std::min(second.dbStartPos, second.dbEndPos);
        int firstDiagonal  = firstQStartPos - firstDbStart;
        int secondDiagonal = secondQStartPos - secondDbStart;
        if (firstDiagonal != secondDiagonal) {
            return firstDiagonal < secondDiagonal;
        }
        return firstDbStart < secondDbStart;
    }

    // map new query into memory (create queryProfile, ...)
    void initQuery(Sequence* query);

    static result_t parseAlignmentRecord(const char *data, bool readCompressed=false);

    static void readAlignmentResults(std::vector<result_t> &result, char *data, bool readCompressed = false);

    static float estimateSeqIdByScorePerCol(uint16_t score, unsigned int qLen, unsigned int tLen);

    static std::string compressAlignment(const std::string &bt);

    static std::string uncompressAlignment(const std::string &cbt);


    static size_t resultToBuffer(char * buffer, const result_t &result, bool addBacktrace, bool compress  = true, bool addOrfPosition = false);

    static int computeAlnLength(int anEnd, int start, int dbEnd, int dbStart);


private:

    // costs to open a gap
    int gapOpen;
    // costs to extend a gap
    int gapExtend;

    // holds values of the current active query
    Sequence * currentQuery;

    // aligner Class
    SmithWaterman * aligner;
    // aligner for nucl
    BandedNucleotideAligner * nuclaligner;
    // substitution matrix
    BaseMatrix* m;
    // evalue
    EvalueComputation * evaluer;
    // byte version of substitution matrix
    int8_t * tinySubMat;
    // set substituion matrix
    void setSubstitutionMatrix(BaseMatrix *m);

};

#endif
