// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>

#include "DistanceCalculator.h"
#include "Matcher.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"
#include <set>
#include <limits>
#include <cstdint>
#include <queue>

#ifdef OPENMP
#include <omp.h>
#endif

class CompareResultBySeqId {
public:
    bool operator() (Matcher::result_t & r1, Matcher::result_t & r2) {
        return r1.seqId < r2.seqId;
    }
};
typedef std::priority_queue<Matcher::result_t, std::vector<Matcher::result_t> , CompareResultBySeqId> SeqIdQueue;
Matcher::result_t selectBestExtentionFragment(SeqIdQueue & alignments,
                                              unsigned int queryKey) {
    // results are ordered by score
    while(alignments.empty() == false){
        Matcher::result_t res = alignments.top();
        alignments.pop();
        size_t dbKey = res.dbKey;
        const bool notRightStartAndLeftStart = !(res.dbStartPos == 0 && res.qStartPos == 0);
        const bool rightStart = res.dbStartPos == 0 && (res.dbEndPos != res.dbLen-1);
        const bool leftStart = res.qStartPos == 0   && (res.qEndPos != res.qLen-1);
        const bool isNotIdentity = (dbKey != queryKey);
        if((rightStart|| leftStart)  && notRightStartAndLeftStart && isNotIdentity){
            return res;
        }
    }
    return Matcher::result_t(UINT_MAX,0,0,0,0,0,0,0,0,0,0,0,0,"");
}


int doassembly(Parameters &par) {

    DBReader<unsigned int> *sequenceDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str());
    sequenceDbr->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> * alnReader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
    alnReader->open(DBReader<unsigned int>::NOSORT);

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    resultWriter.open();
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, 0.0f);
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(subMat);

    unsigned char * wasExtended = new unsigned char[sequenceDbr->getSize()];
    std::fill(wasExtended, wasExtended+sequenceDbr->getSize(), 0);

#pragma omp parallel
    {
#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
            Debug::printProgress(id);
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            unsigned int queryId = sequenceDbr->getDbKey(id);
            char *querySeq = sequenceDbr->getData(id);
            unsigned int querySeqLen = sequenceDbr->getSeqLens(id) - 2;
            unsigned int leftQueryOffset = 0;
            unsigned int rightQueryOffset = 0;
            std::string query(querySeq, querySeqLen); // no /n/0
            char *alnData = alnReader->getDataByDBKey(queryId);
            std::vector<Matcher::result_t> alignments = Matcher::readAlignmentResults(alnData);
            SeqIdQueue alnQueue;
            for (size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++) {
                alnQueue.push(alignments[alnIdx]);
            }

            bool queryCouldBeExtended = false;
            Matcher::result_t besttHitToExtend;
            while( (besttHitToExtend = selectBestExtentionFragment(alnQueue, queryId)).dbKey != UINT_MAX) {
                querySeqLen = query.size();
                querySeq = (char *) query.c_str();

//                querySeq.mapSequence(id, queryKey, query.c_str());
                unsigned int targetId = sequenceDbr->getId(besttHitToExtend.dbKey);
                if(targetId==UINT_MAX){
                    Debug(Debug::ERROR) << "Could not find targetId  " << besttHitToExtend.dbKey
                                        << " in database "<< sequenceDbr->getDataFileName() <<  "\n";
                    EXIT(EXIT_FAILURE);
                }
                char *targetSeq = sequenceDbr->getData(targetId);
                unsigned int targetSeqLen = sequenceDbr->getSeqLens(targetId) - 2;
                // check if alignment still make sense (can extend the query)
                if (besttHitToExtend.dbStartPos == 0   ) {
                    if((targetSeqLen - (besttHitToExtend.dbEndPos + 1)) <= rightQueryOffset){
                        continue;
                    }
                }else if ( besttHitToExtend.qStartPos == 0 ){
                    if(besttHitToExtend.dbStartPos <= leftQueryOffset){
                        continue;
                    }
                }
                int qStartPos, qEndPos, dbStartPos, dbEndPos;
                int diagonal = (leftQueryOffset + besttHitToExtend.qStartPos) - besttHitToExtend.dbStartPos;
                size_t dist = std::max(abs(diagonal), 0);
                if (diagonal >= 0) {
//                    targetSeq.mapSequence(targetId, besttHitToExtend.dbKey, dbSeq);
                    size_t diagonalLen = std::min(targetSeqLen, querySeqLen - abs(diagonal));
                    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeSubstituionStartEndDistance(
                            querySeq + abs(diagonal),
                            targetSeq, diagonalLen, fastMatrix.matrix);
                    qStartPos = alignment.startPos + dist;
                    qEndPos = alignment.endPos + dist;
                    dbStartPos = alignment.startPos;
                    dbEndPos = alignment.endPos;
                } else {
                    size_t diagonalLen = std::min(targetSeqLen - abs(diagonal), querySeqLen);
                    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeSubstituionStartEndDistance(
                            querySeq,
                            targetSeq + abs(diagonal),
                            diagonalLen, fastMatrix.matrix);
                    qStartPos = alignment.startPos;
                    qEndPos = alignment.endPos;
                    dbStartPos = alignment.startPos + dist;
                    dbEndPos = alignment.endPos + dist;
                }

//                std::cout << "\t" << besttHitToExtend.dbKey << std::endl;
                //std::cout << "Query : " << query << std::endl;
                //std::cout << "Target: " << std::string(dbSeq, targetSeq.L)  << std::endl;

                if (dbStartPos == 0 && qEndPos == (querySeqLen -1)  ) {
                    size_t dbFragLen = (targetSeqLen - dbEndPos) - 1; // -1 get not aligned element
                    std::string fragment = std::string(targetSeq + dbEndPos + 1, dbFragLen);
                    if (fragment.size() + query.size() >= par.maxSeqLen) {
                        Debug(Debug::WARNING) << "Sequence too long in query id: " << queryId << ". "
                                "Max length allowed would is " << par.maxSeqLen << "\n";
                        break;
                    }
                    //update that dbKey was used in assembly
                    __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));
                    queryCouldBeExtended = true;
                    query += fragment;
                    rightQueryOffset += dbFragLen;

                } else if (qStartPos == 0 && dbEndPos == (targetSeqLen - 1)) {
                    std::string fragment = std::string(targetSeq, dbStartPos); // +1 get not aligned element
                    if (fragment.size() + query.size() >= par.maxSeqLen) {
                        Debug(Debug::WARNING) << "Sequence too long in query id: " << queryId << ". "
                                "Max length allowed would is " << par.maxSeqLen << "\n";
                        break;
                    }
                    // update that dbKey was used in assembly
                    __sync_or_and_fetch(&wasExtended[targetId], static_cast<unsigned char>(0x80));
                    queryCouldBeExtended = true;
                    query = fragment + query;
                    leftQueryOffset += dbStartPos;
                }

            }
            if (queryCouldBeExtended == true) {
                query.push_back('\n');
                __sync_or_and_fetch(&wasExtended[queryId], static_cast<unsigned char>(0x80));
                resultWriter.writeData(query.c_str(), query.size(), queryId, thread_idx);
            }
        }
    } // end parallel
// add sequences that are not yet assembled
#pragma omp parallel for schedule(dynamic, 10000)
    for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        if(wasExtended[id] == 0){
            char *querySeqData = sequenceDbr->getData(id);
            unsigned int queryLen = sequenceDbr->getSeqLens(id) - 1; //skip null byte
            resultWriter.writeData(querySeqData, queryLen, sequenceDbr->getDbKey(id), thread_idx);
        }
    }

// cleanup
    resultWriter.close();
    alnReader->close();
    delete [] wasExtended;
    delete alnReader;
    delete [] fastMatrix.matrix;
    delete [] fastMatrix.matrixData;
    sequenceDbr->close();
    delete sequenceDbr;
    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}

int assembleresult(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::WARNING) << "Compute assembly.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

    int retCode = doassembly(par);

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for processing: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    return retCode;
}
