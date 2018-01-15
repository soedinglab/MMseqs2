#include "DistanceCalculator.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "QueryMatcher.h"
#include "CovSeqidQscPercMinDiag.out.h"
#include "CovSeqidQscPercMinDiagTargetCov.out.h"
#include "QueryMatcher.h"

#include <string>
#include <vector>
#include <sys/time.h>
#include <NucleotideMatrix.h>

#ifdef OPENMP
#include <omp.h>
#endif

float parsePrecisionLib(std::string scoreFile, double targetSeqid, double targetCov, double targetPrecision){
    std::stringstream in(scoreFile);
    std::string line;
    // find closest lower seq. id in a grid of size 5
    int intTargetSeqid = static_cast<int>((targetSeqid+0.0001)*100);
    int seqIdRest = (intTargetSeqid % 5);
    targetSeqid = static_cast<float>(intTargetSeqid- seqIdRest)/100;
    // find closest lower cov. id in a grid of size 10
    targetCov = static_cast<float>(static_cast<int>((targetCov + 0.0001)* 10))/10;
    while (std::getline(in, line)) {
        std::vector<std::string> values = Util::split(line, " ");
        float cov = strtod(values[0].c_str(),NULL);
        float seqid = strtod(values[1].c_str(),NULL);
        float scorePerCol = strtod(values[2].c_str(),NULL);
        float precision = strtod(values[3].c_str(),NULL);
        if(MathUtil::AreSame(cov, targetCov)
           && MathUtil::AreSame(seqid, targetSeqid)
           && precision >= targetPrecision){
            return scorePerCol;
        }
    }
    Debug(Debug::WARNING) << "Could not find any score per column for "
                             "cov="<< targetCov << " seq.id=" << targetSeqid << "\n"
                             "No hit will be filtered.\n";

    return 0;
}

int rescorediagonal(int argc, const char **argv, const Command &command) {
    Debug(Debug::INFO) << "Rescore diagonals.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> *qdbr = NULL;

    Debug(Debug::INFO) << "Query  file: " << par.db1 << "\n";
    qdbr = new DBReader<unsigned int>(par.db1.c_str(), (par.db1 + ".index").c_str());
    qdbr->open(DBReader<unsigned int>::NOSORT);
    int querySeqType  =  qdbr->getDbtype();

    DBReader<unsigned int> *tdbr = NULL;
    BaseMatrix *subMat;
    if (querySeqType == Sequence::NUCLEOTIDES) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
    } else {
        // keep score bias at 0.0 (improved ROC)
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    }


    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);

//    int sequenceType = Sequence::AMINO_ACIDS;
//    if (par.queryProfile == true) {
//        sequenceType = Sequence::HMM_PROFILE;
//    }


    qdbr->readMmapedDataInMemory();
    Debug(Debug::INFO) << "Target  file: " << par.db2 << "\n";
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), (par.db2 + ".index").c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        tdbr->readMmapedDataInMemory();
    }
    float scorePerColThr = 0.0;
    if(par.filterHits){
        if(par.rescoreMode == Parameters::RESCORE_MODE_HAMMING){
            Debug(Debug::WARNING) << "HAMMING distance can not be used to filter hits. Using --rescore-mode 1\n";
            par.rescoreMode = Parameters::RESCORE_MODE_SUBSTITUTION;
        }

        std::string libraryString = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
           ? std::string((const char*)CovSeqidQscPercMinDiag_out, CovSeqidQscPercMinDiag_out_len)
           : std::string((const char*)CovSeqidQscPercMinDiagTargetCov_out, CovSeqidQscPercMinDiagTargetCov_out_len);
        scorePerColThr = parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.99);
    }
    double * kmnByLen = new double[par.maxSeqLen];
    EvalueComputation evaluer(tdbr->getAminoAcidDBSize(), subMat, Matcher::GAP_OPEN, Matcher::GAP_EXTEND, false);
    DistanceCalculator globalAliStat;
    if (par.globalAlignment)
    {
        globalAliStat.prepareGlobalAliParam(*subMat);
    }
    Debug(Debug::WARNING) << "Prefilter database: " << par.db3 << "\n";
    DBReader<unsigned int> dbr_res(par.db3.c_str(), std::string(par.db3 + ".index").c_str());
    dbr_res.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    Debug(Debug::INFO) << "Result database: " << par.db4 << "\n";
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads);
    resultWriter.open();
    const size_t flushSize = 100000000;
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbr_res.getSize()) / static_cast<double>(flushSize)));
    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(dbr_res.getSize() - (i * flushSize), flushSize);
#pragma omp parallel
        {
#pragma omp for schedule(dynamic, 1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                Debug::printProgress(id);
                char buffer[1024+32768];
                std::string prefResultsOutString;
                prefResultsOutString.reserve(1000000);
                unsigned int thread_idx = 0;
#ifdef OPENMP
                thread_idx = (unsigned int) omp_get_thread_num();
#endif
                char *data = dbr_res.getData(id);
                unsigned int queryId = qdbr->getId(dbr_res.getDbKey(id));
                char *querySeq = qdbr->getData(queryId);
                int queryLen = std::max(0, static_cast<int>(qdbr->getSeqLens(queryId)) - 2);


//                if(par.rescoreMode != Parameters::RESCORE_MODE_HAMMING){
//                    query.mapSequence(id, queryId, querySeq);
//                    queryLen = query.L;
//                }else{
                    // -2 because of \n\0 in sequenceDB
//                }
                std::vector<hit_t> results = QueryMatcher::parsePrefilterHits(data);
                for (size_t entryIdx = 0; entryIdx < results.size(); entryIdx++) {
                    unsigned int targetId = tdbr->getId(results[entryIdx].seqId);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;
                    char * targetSeq = tdbr->getData(targetId);
                    int targetLen = std::max(0, static_cast<int>(tdbr->getSeqLens(targetId)) - 2);
                    short diagonal = results[entryIdx].diagonal;
                    unsigned short distanceToDiagonal = abs(diagonal);
                    unsigned int diagonalLen = 0;
                    unsigned int distance = 0;
                    DistanceCalculator::LocalAlignment alignment;

                    if (diagonal >= 0 && distanceToDiagonal < queryLen) {
                        diagonalLen = std::min(targetLen, queryLen - distanceToDiagonal);
                        if(par.rescoreMode == Parameters::RESCORE_MODE_HAMMING){
                            distance = DistanceCalculator::computeHammingDistance(querySeq + distanceToDiagonal,
                                                                                  targetSeq, diagonalLen);
                        }else if(par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION) {
                            distance = DistanceCalculator::computeSubstitutionDistance(querySeq + distanceToDiagonal,
                                                                                       targetSeq,
                                                                                       diagonalLen, fastMatrix.matrix,
                                                                                       par.globalAlignment);
                        }else if(par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT){
                            alignment = DistanceCalculator::computeSubstitutionStartEndDistance(
                                    querySeq + distanceToDiagonal,
                                    targetSeq,
                                    diagonalLen, fastMatrix.matrix);
                            distance = alignment.score;
                        }
                    } else if (diagonal < 0 && distanceToDiagonal < targetLen) {
                        diagonalLen = std::min(targetLen - distanceToDiagonal, queryLen);
                        if(par.rescoreMode == Parameters::RESCORE_MODE_HAMMING){
                            distance = DistanceCalculator::computeHammingDistance(querySeq,
                                                                                  targetSeq + distanceToDiagonal,
                                                                                  diagonalLen);
                        }else if(par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION){
                            distance = DistanceCalculator::computeSubstitutionDistance(querySeq,
                                                                                       targetSeq + distanceToDiagonal,
                                                                                       diagonalLen, fastMatrix.matrix,
                                                                                       par.globalAlignment);
                        }else if(par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT){
                            alignment = DistanceCalculator::computeSubstitutionStartEndDistance(querySeq,
                                                                                                targetSeq +
                                                                                                distanceToDiagonal,
                                                                                                diagonalLen,
                                                                                                fastMatrix.matrix);
                            distance = alignment.score;
                        }
                    }
                    double seqId = 0;
                    double evalue = 0.0;
                    float targetCov = static_cast<float>(diagonalLen) / static_cast<float>(targetLen);
                    float queryCov = static_cast<float>(diagonalLen) / static_cast<float>(queryLen);
                    Matcher::result_t result;
                    if(par.rescoreMode == Parameters::RESCORE_MODE_HAMMING){
                        seqId = (static_cast<float>(diagonalLen) - static_cast<float>(distance)) /
                                static_cast<float>(diagonalLen);
                    }else if(par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION || par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT){
                        //seqId = exp(static_cast<float>(distance) / static_cast<float>(diagonalLen));
                        if (par.globalAlignment)
                            seqId = globalAliStat.getPvalGlobalAli((float)distance,diagonalLen);
                        else{
                            evalue = evaluer.computeEvalue(distance, queryLen);
                            int bitScore =static_cast<short>(evaluer.computeBitScore(distance)+0.5);

                            if(par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT) {
                                int alnLen = alignment.endPos - alignment.startPos;
                                int qStartPos, qEndPos, dbStartPos, dbEndPos;
                                // -1 since diagonal is computed from sequence Len which starts by 1
                                if (diagonal >= 0) {
                                    qStartPos = alignment.startPos + distanceToDiagonal;
                                    qEndPos = alignment.endPos + distanceToDiagonal;
                                    dbStartPos = alignment.startPos;
                                    dbEndPos = alignment.endPos;
                                } else {
                                    qStartPos = alignment.startPos;
                                    qEndPos = alignment.endPos;
                                    dbStartPos = alignment.startPos + distanceToDiagonal;
                                    dbEndPos = alignment.endPos + distanceToDiagonal;
                                }
                                int qAlnLen = std::max(qEndPos - qStartPos, static_cast<int>(1));
                                int dbAlnLen = std::max(dbEndPos - dbStartPos, static_cast<int>(1));
                                //seqId = (alignment.score1 / static_cast<float>(std::max(qAlnLength, dbAlnLength)))  * 0.1656 + 0.1141;
                                seqId = Matcher::estimateSeqIdByScorePerCol(distance, qAlnLen, dbAlnLen);
                                // compute seq.id if hit fulfills e-value but not by seqId criteria
                                if(evalue <= par.evalThr && seqId > par.seqIdThr - 0.15){
                                    int idCnt = 0;
                                    for(int i = qStartPos; i < qEndPos; i++){
                                        idCnt += (querySeq[i] == targetSeq[dbStartPos+(i-qStartPos)]) ? 1 : 0;
                                    }
                                    seqId = static_cast<double>(idCnt) / (static_cast<double>(qEndPos) - static_cast<double>(qStartPos));
                                }
                                result = Matcher::result_t(results[entryIdx].seqId, bitScore, queryCov, targetCov, seqId, evalue, alnLen,
                                                           qStartPos, qEndPos, queryLen, dbStartPos, dbEndPos, targetLen, std::string());
                            }
                        }
                    }

                    //float maxSeqLen = std::max(static_cast<float>(targetLen), static_cast<float>(queryLen));
                    float currScorePerCol = static_cast<float>(distance)/static_cast<float>(diagonalLen);
                    // query/target cov mode
                    bool hasCov = Util::hasCoverage(par.covThr, par.covMode, queryCov, targetCov);
                    // --min-seq-id
                    bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    bool hasEvalue = (evalue <= par.evalThr);
                    // --filter-hits
                    bool hasToFilter = (par.filterHits == true  && currScorePerCol >= scorePerColThr);
                    if (isIdentity || hasToFilter || (hasCov && hasSeqId && hasEvalue))
                    {
                        int len  = 0;
                        if(par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT) {
                            size_t len = Matcher::resultToBuffer(buffer, result, false);
//                            len = snprintf(buffer, 100, "%s", alnString.c_str());
                            prefResultsOutString.append(buffer, len);
                        } else if(par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION){
                            len = snprintf(buffer, 100, "%u\t%.3e\t%d\n", results[entryIdx].seqId, evalue, diagonal);
                            prefResultsOutString.append(buffer, len);
                        } else {
                            len = snprintf(buffer, 100, "%u\t%.2f\t%d\n", results[entryIdx].seqId, seqId, diagonal);
                            prefResultsOutString.append(buffer, len);
                        }
                    }
                }
                // write prefiltering results string to ffindex database
                const size_t prefResultsLength = prefResultsOutString.length();
                char *prefResultsOutData = (char *) prefResultsOutString.c_str();
                resultWriter.writeData(prefResultsOutData, prefResultsLength, qdbr->getDbKey(queryId), thread_idx);
            }
        }
        dbr_res.remapData();
    }
    Debug(Debug::INFO) << "Done." << "\n";
    dbr_res.close();
    resultWriter.close();
    qdbr->close();
    delete qdbr;
    delete subMat;
    delete [] kmnByLen;
    delete [] fastMatrix.matrix;
    delete [] fastMatrix.matrixData;
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }
    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for diagonal calculation: " << (sec / 3600) << " h "
                       << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    return EXIT_SUCCESS;
}


