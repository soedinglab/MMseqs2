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
#ifdef OPENMP
#include <omp.h>
#endif


Matcher::result_t  selectBestExtentionFragment(DBReader<unsigned int>  * sequenceDbr,
                                              std::vector<Matcher::result_t> alignments,
                                              std::set<unsigned int> &prevFound,
                                              unsigned int queryKey) {
    // results are ordered by score
    double evalue = DBL_MAX;
    size_t bestIndex = SIZE_MAX;
    for(size_t i = 0; i < alignments.size(); i++){
        const bool foundPrev = (prevFound.find(alignments[i].dbKey) != prevFound.end());
        if(foundPrev == false){
            size_t dbKey = alignments[i].dbKey;
            const bool notRightStartAndLeftStart = !(alignments[i].dbStartPos == 0 && alignments[i].qStartPos == 0);
            const bool betterEvalue = (alignments[i].eval < evalue);
            const bool rightStart = alignments[i].dbStartPos == 0 ;
            const bool leftStart = alignments[i].qStartPos == 0;
            const bool isNotIdentity = (dbKey != queryKey);
            if((rightStart|| leftStart) && betterEvalue && notRightStartAndLeftStart && isNotIdentity){
                // mark label as visited
                evalue=alignments[i].eval;
                bestIndex = i;
            }
        }
    }
    if(bestIndex == SIZE_MAX) {
        return Matcher::result_t(UINT_MAX,0,0,0,0,0,0,0,0,0,0,0,0,"");
    } else {
        return alignments[bestIndex];
    }
}


int doassembly(Parameters &par) {

    DBReader<unsigned int> *sequenceDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str());
    sequenceDbr->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> * alnReader = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
    alnReader->open(DBReader<unsigned int>::NOSORT);

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    resultWriter.open();
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, 0.0f);

//    unsigned char * readLabel = new unsigned char[sequenceDbr->getSize()];
//    std::fill(readLabel, readLabel+sequenceDbr->getSize(), 0);

    // do some fancy quality control
    //labaleSequences(par, sequenceDbr, alnReader, subMat, readLabel);


#pragma omp parallel
    {
        Sequence querySeq(par.maxSeqLen, subMat.aa2int, subMat.int2aa, par.querySeqType,
                          0, false, false);
        Sequence targetSeq(par.maxSeqLen, subMat.aa2int, subMat.int2aa, par.querySeqType,
                          0, false, false);
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < sequenceDbr->getSize(); id++) {
            Debug::printProgress(id);
//            std::cout << id << std::endl;
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
//            if(readLabel[id]!=FULLYCOVERED){
//                continue;
//            }
            // set label as visited
//            __sync_or_and_fetch(&readLabel[id], static_cast<unsigned char>(0x80));

            unsigned int queryKey = sequenceDbr->getDbKey(id);
            char * querySeqData = sequenceDbr->getData(id);
            unsigned int queryLen = sequenceDbr->getSeqLens(id);

            unsigned int queryOffset = 0;
            std::string query(querySeqData, queryLen-2); // no /n/0
            std::vector<Matcher::result_t> alignments = Matcher::readAlignmentResults (  alnReader->getDataByDBKey(queryKey) );
            std::set<unsigned int> prevFound;
            for(size_t alnIdx = 0; alnIdx < alignments.size(); alnIdx++){
                querySeq.mapSequence(id, queryKey, query.c_str());
                Matcher::result_t besttHitToExtend = selectBestExtentionFragment(sequenceDbr,
                                                                                 alignments,
                                                                                 prevFound,
                                                                                 queryKey);
                if(besttHitToExtend.dbKey == UINT_MAX)
                    continue;
                prevFound.emplace(besttHitToExtend.dbKey);
                char * dbSeq = sequenceDbr->getDataByDBKey(besttHitToExtend.dbKey);
                targetSeq.mapSequence(besttHitToExtend.dbKey, besttHitToExtend.dbKey, dbSeq);
                int qStartPos, qEndPos, dbStartPos, dbEndPos;
                int diagonal = (queryOffset+besttHitToExtend.qStartPos) - besttHitToExtend.dbStartPos;
                size_t dist = std::max(abs(diagonal) - 1, 0);
                if (diagonal >= 0){
                    size_t diagonalLen = std::min(targetSeq.L, querySeq.L - abs(diagonal));
                    DistanceCalculator::LocalAlignment alignment = DistanceCalculator::computeSubstituionStartEndDistance(querySeq.int_sequence + abs(diagonal),
                                                                                            targetSeq.int_sequence, diagonalLen, subMat.subMatrix);
                    qStartPos = alignment.startPos + dist;
                    dbStartPos = alignment.startPos;
                    dbEndPos = alignment.endPos;
                }else{
                    size_t diagonalLen = std::min(targetSeq.L - abs(diagonal), querySeq.L);
                    DistanceCalculator::LocalAlignment alignment  = DistanceCalculator::computeSubstituionStartEndDistance(querySeq.int_sequence,
                                                                                            targetSeq.int_sequence + abs(diagonal),
                                                                                       diagonalLen, subMat.subMatrix);
                    qStartPos = alignment.startPos;
                    dbStartPos = alignment.startPos + dist;
                    dbEndPos = alignment.endPos + dist;
//                    __sync_fetch_and_add(&readLabel[id], 1);

                }


//                std::cout << "\t" << besttHitToExtend.dbKey << std::endl;

//                std::cout << "Query : " << query << std::endl;
//                std::cout << "Target: " << std::string(dbSeq, targetSeq.L)  << std::endl;

                if(dbStartPos==0){
                    size_t dbFragLen = (targetSeq.L-dbEndPos) - 1; // -1 get not aligned element
                    std::string fragment = std::string(dbSeq+dbEndPos+1, dbFragLen);
                    if(fragment.size() + query.size() >= par.maxSeqLen){
                        Debug(Debug::WARNING) << "Sequence too long in query id: "  << queryKey << ". "
                                                 "Max length allowed would is " << par.maxSeqLen << "\n";
                        break;
                    }
//                    std::cout << "Fargm1: "  << fragment << std::endl;
                    query += fragment;
                }else if(qStartPos==0){
                    std::string fragment = std::string(dbSeq, dbStartPos + 1); // +1 get not aligned element
                    if(fragment.size() + query.size() >= par.maxSeqLen){
                        Debug(Debug::WARNING) << "Sequence too long in query id: "  << queryKey << ". "
                                                 "Max length allowed would is " << par.maxSeqLen << "\n";
                        break;
                    }
//                    std::cout << "Fargm2: "  << fragment << std::endl;
                    query = fragment + query;
                    queryOffset += dbStartPos + 1;
                }

            }
            query.push_back('\n');
            resultWriter.writeData(query.c_str(),query.size(), queryKey, thread_idx);
        }

    }

// cleanup
    resultWriter.close();
    alnReader->close();
//    delete [] readLabel;
    delete alnReader;
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
