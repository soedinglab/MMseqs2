//
// Created by mad on 10/21/15.
//
#include <string>
#include <vector>
#include <Prefiltering.h>
#include <DistanceCalculator.h>
#include <sys/time.h>
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "CovSeqidQscPercMinDiag.out.h"
#include "CovSeqidQscPercMinDiagTargetCov.out.h"


float parsePrecisionLib(std::string scoreFile, double targetSeqid, double targetCov, double targetPrecision){
    std::stringstream in(scoreFile);
    std::string line;
    double qsc = 0;
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
        float percision = strtod(values[3].c_str(),NULL);
        if(MathUtil::AreSame(cov, targetCov)
           && MathUtil::AreSame(seqid, targetSeqid)
           && percision >= targetPrecision){
            return scorePerCol;
        }
    }
    Debug(Debug::WARNING) << "Could not find any score per column for "
                             "cov="<< targetCov << " seq.id=" << targetSeqid << "\n"
                             "No hit will be filtered.\n";

    return qsc;
}

int rescorediagonal(int argc, const char **argv, const Command &command) {
    Debug(Debug::WARNING) << "Rescore diagonals.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    DBReader<unsigned int> *qdbr = NULL;
    DBReader<unsigned int> *tdbr = NULL;
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, 0.0f);
    int sequenceType = Sequence::AMINO_ACIDS;
    if (par.queryProfile == true) {
        sequenceType = Sequence::HMM_PROFILE;
    }

    Debug(Debug::WARNING) << "Query  file: " << par.db1 << "\n";
    qdbr = new DBReader<unsigned int>(par.db1.c_str(), (par.db1 + ".index").c_str());
    qdbr->open(DBReader<unsigned int>::NOSORT);
    qdbr->readMmapedDataInMemory();
    Debug(Debug::WARNING) << "Target  file: " << par.db2 << "\n";
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
        if(par.rescoreMode==Parameters::RESCORE_MODE_HAMMING){
            Debug(Debug::ERROR) << "HAMMING distance can not be used to filter hits. "
                                   "Please use --rescore-mode 1\n";
            EXIT(EXIT_FAILURE);
        }
        if(par.targetCovThr > 0.0){
            std::string libraryString((const char*)CovSeqidQscPercMinDiagTargetCov_out,CovSeqidQscPercMinDiagTargetCov_out_len);
            scorePerColThr = parsePrecisionLib(libraryString, par.seqIdThr, par.targetCovThr, 0.99);
        }else{
            std::string libraryString((const char*)CovSeqidQscPercMinDiag_out,CovSeqidQscPercMinDiag_out_len);
            scorePerColThr = parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.99);
        }
    }
    double * kmnByLen = new double[par.maxSeqLen];
    BlastScoreUtils::BlastStat stats = BlastScoreUtils::getAltschulStatsForMatrix(subMat.getMatrixName(), 11, 1, false);
    for(int len = 0; len < par.maxSeqLen; len++){
        kmnByLen[len] = BlastScoreUtils::computeKmn(len, stats.K, stats.lambda, stats.alpha, stats.beta,
                                                    tdbr->getAminoAcidDBSize(), tdbr->getSize());
    }

    Debug(Debug::WARNING) << "Prefilter database: " << par.db3 << "\n";
    DBReader<unsigned int> dbr_res(par.db3.c_str(), std::string(par.db3 + ".index").c_str());
    dbr_res.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    Debug(Debug::WARNING) << "Result database: " << par.db4 << "\n";
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBWriter::BINARY_MODE);
    resultWriter.open();
    const size_t flushSize = 100000000;
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbr_res.getSize()) / static_cast<double>(flushSize)));
    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(dbr_res.getSize() - (i * flushSize), flushSize);
#pragma omp parallel
        {
            Sequence query(par.maxSeqLen, subMat.aa2int, subMat.int2aa, sequenceType, 0, false, par.compBiasCorrection);
            Sequence target(par.maxSeqLen, subMat.aa2int, subMat.int2aa, sequenceType, 0, false, par.compBiasCorrection);
#pragma omp for schedule(dynamic, 1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                Debug::printProgress(id);
                char buffer[100];
                std::string prefResultsOutString;
                prefResultsOutString.reserve(1000000);
                unsigned int thread_idx = 0;
#ifdef OPENMP
                thread_idx = (unsigned int) omp_get_thread_num();
#endif
                char *data = dbr_res.getData(id);
                unsigned int queryId = qdbr->getId(dbr_res.getDbKey(id));
                char *querySeq = qdbr->getData(queryId);
                query.mapSequence(id, queryId, querySeq);
                unsigned int queryLen = query.L;
                std::vector<hit_t> results = Prefiltering::readPrefilterResults(data);
                for (size_t entryIdx = 0; entryIdx < results.size(); entryIdx++) {
                    unsigned int targetId = tdbr->getId(results[entryIdx].seqId);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB))? true : false;

                    target.mapSequence(0, targetId, tdbr->getData(targetId));
                    unsigned int targetLen = target.L;
                    short diagonal = results[entryIdx].diagonal;
                    unsigned short distanceToDiagonal = abs(diagonal);
                    unsigned int diagonalLen = 0;
                    unsigned int distance = 0;
                    DistanceCalculator::LocalAlignment alignment;

                    if (diagonal >= 0 && distanceToDiagonal < queryLen) {
                        diagonalLen = std::min(targetLen, queryLen - distanceToDiagonal);
                        if(par.rescoreMode == Parameters::RESCORE_MODE_HAMMING){
                            distance = DistanceCalculator::computeHammingDistance(querySeq + distanceToDiagonal,
                                                                                  tdbr->getData(targetId), diagonalLen);
                        }else if(par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION) {
                            distance = DistanceCalculator::computeSubstituionDistance(query.int_sequence + distanceToDiagonal,
                                                                                      target.int_sequence,
                                                                                      diagonalLen, subMat.subMatrix,par.globalAlignment);
                        }else if(par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT){
                            alignment = DistanceCalculator::computeSubstituionStartEndDistance(query.int_sequence + distanceToDiagonal,
                                                                                               target.int_sequence,
                                                                                               diagonalLen, subMat.subMatrix);
                            distance = alignment.score;
                        }
                    } else if (diagonal < 0 && distanceToDiagonal < targetLen) {
                        diagonalLen = std::min(targetLen - distanceToDiagonal, queryLen);
                        if(par.rescoreMode == Parameters::RESCORE_MODE_HAMMING){
                            distance = DistanceCalculator::computeHammingDistance(querySeq,
                                                                                  tdbr->getData(targetId) + distanceToDiagonal,
                                                                                  diagonalLen);
                        }else if(par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION){
                            distance = DistanceCalculator::computeSubstituionDistance(query.int_sequence,
                                                                                      target.int_sequence + distanceToDiagonal,
                                                                                      diagonalLen, subMat.subMatrix,par.globalAlignment);
                        }else if(par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT){
                            alignment = DistanceCalculator::computeSubstituionStartEndDistance(query.int_sequence,
                                                                                  target.int_sequence + distanceToDiagonal,
                                                                                  diagonalLen, subMat.subMatrix);
                            distance = alignment.score;
                        }
                    }
                    double seqId = 0;
                    if(par.rescoreMode == Parameters::RESCORE_MODE_HAMMING){
                        seqId = (static_cast<float>(diagonalLen) - static_cast<float>(distance)) /
                                static_cast<float>(diagonalLen);
                    }else if(par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION || par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT){
                        //seqId = exp(static_cast<float>(distance) / static_cast<float>(diagonalLen));
                        seqId = BlastScoreUtils::computeEvalue(distance, kmnByLen[queryLen], stats.lambda);
                    }
                    float targetCov = static_cast<float>(diagonalLen) / static_cast<float>(targetLen);
                    float queryCov = static_cast<float>(diagonalLen) / static_cast<float>(queryLen);
                    //float maxSeqLen = std::max(static_cast<float>(targetLen), static_cast<float>(queryLen));
                    float currScorePerCol = static_cast<float>(distance)/static_cast<float>(diagonalLen);
                    // --target-cov
                    bool hasTargetCov =  targetCov >= (par.targetCovThr - std::numeric_limits<float>::epsilon());
                    // -c
                    bool hasCov = queryCov >= par.covThr && targetCov >= par.covThr;
                    // --min-seq-id
                    bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    // --filter-hits
                    bool hasToFilter = (par.filterHits == true  && currScorePerCol >= scorePerColThr);
                    if (isIdentity || hasToFilter || (hasTargetCov && hasCov && hasSeqId))
                    {
                        int len  = 0;
                        if(par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT) {
                            int alnLen = alignment.endPos - alignment.startPos;
                            int qStartPos = 0;
                            int qEndPos = 0;
                            int dbStartPos = 0;
                            int dbEndPos = 0;
                            // -1 since diagonal is computed from sequence Len which starts by 1
                            size_t dist = std::max(distanceToDiagonal - 1, 0);
                            if(diagonal >= 0){
                                qStartPos = alignment.startPos + dist;
                                qEndPos = alignment.endPos + dist ;
                                dbStartPos = alignment.startPos;
                                dbEndPos = alignment.endPos;
                            }else{
                                qStartPos = alignment.startPos;
                                qEndPos = alignment.endPos;
                                dbStartPos = alignment.startPos + dist;
                                dbEndPos = alignment.endPos + dist;
                            }

                            Matcher::result_t result(results[entryIdx].seqId, distance, queryCov, targetCov, seqId, seqId, alnLen, qStartPos, qEndPos, queryLen, dbStartPos, dbEndPos, targetLen, std::string());
                            std::string alnString = Matcher::resultToString(result, false);
                            len = snprintf(buffer, 100, "%s", alnString.c_str());
                        } else   if(par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION){
                            if (par.globalAlignment) // in case of global alignemnts, eval does not make sense->write the score
                                len = snprintf(buffer, 100, "%u\t%u\t%d\n", results[entryIdx].seqId, distance, diagonal);
                            else
                                len = snprintf(buffer, 100, "%u\t%.3e\t%d\n", results[entryIdx].seqId, seqId, diagonal);
                            
                        } else {
                            len = snprintf(buffer, 100, "%u\t%.2f\t%d\n", results[entryIdx].seqId, seqId, diagonal);
                        }
                        prefResultsOutString.append(buffer, len);
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
    Debug(Debug::WARNING) << "Done." << "\n";
    dbr_res.close();
    resultWriter.close();
    qdbr->close();
    delete qdbr;
    delete [] kmnByLen;
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for diagonal calculation: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    return 0;
}


