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
#include "NucleotideMatrix.h"

#ifdef OPENMP
#include <omp.h>
#endif

float parsePrecisionLib(const std::string &scoreFile, double targetSeqid, double targetCov, double targetPrecision) {
    std::stringstream in(scoreFile);
    std::string line;
    // find closest lower seq. id in a grid of size 5
    int intTargetSeqid = static_cast<int>((targetSeqid + 0.0001) * 100);
    int seqIdRest = (intTargetSeqid % 5);
    targetSeqid = static_cast<float>(intTargetSeqid - seqIdRest) / 100;
    // find closest lower cov. id in a grid of size 10
    targetCov = static_cast<float>(static_cast<int>((targetCov + 0.0001) * 10)) / 10;
    while (std::getline(in, line)) {
        std::vector<std::string> values = Util::split(line, " ");
        float cov = strtod(values[0].c_str(), NULL);
        float seqid = strtod(values[1].c_str(), NULL);
        float scorePerCol = strtod(values[2].c_str(), NULL);
        float precision = strtod(values[3].c_str(), NULL);
        if (MathUtil::AreSame(cov, targetCov) && MathUtil::AreSame(seqid, targetSeqid) && precision >= targetPrecision) {
            return scorePerCol;
        }
    }
    Debug(Debug::WARNING) << "Could not find any score per column for cov "
                          << targetCov << " seq.id. " << targetSeqid << ". No hit will be filtered.\n";

    return 0;
}

int doRescorediagonal(Parameters &par,
                      DBWriter &resultWriter,
                      DBReader<unsigned int> &resultReader,
              const size_t dbFrom, const size_t dbSize) {
    Debug(Debug::INFO) << "Query database: " << par.db1 << "\n";
    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str());
    qdbr.open(DBReader<unsigned int>::NOSORT);
    const int querySeqType = qdbr.getDbtype();
    if (par.noPreload == false) {
        qdbr.readMmapedDataInMemory();
    }

    BaseMatrix *subMat;
    if (querySeqType == Sequence::NUCLEOTIDES) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
    } else {
        // keep score bias at 0.0 (improved ROC)
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    }

    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(*subMat);

    Debug(Debug::INFO) << "Target database: " << par.db2 << "\n";
    DBReader<unsigned int> *tdbr = NULL;
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.noPreload == false) {
            tdbr->readMmapedDataInMemory();
        }
    }

    float scorePerColThr = 0.0;
    if (par.filterHits) {
        if (par.rescoreMode == Parameters::RESCORE_MODE_HAMMING) {
            Debug(Debug::WARNING) << "HAMMING distance can not be used to filter hits. Using --rescore-mode 1\n";
            par.rescoreMode = Parameters::RESCORE_MODE_SUBSTITUTION;
        }

        std::string libraryString = (par.covMode == Parameters::COV_MODE_BIDIRECTIONAL)
                                    ? std::string((const char*)CovSeqidQscPercMinDiag_out, CovSeqidQscPercMinDiag_out_len)
                                    : std::string((const char*)CovSeqidQscPercMinDiagTargetCov_out, CovSeqidQscPercMinDiagTargetCov_out_len);
        scorePerColThr = parsePrecisionLib(libraryString, par.seqIdThr, par.covThr, 0.99);
    }

    EvalueComputation evaluer(tdbr->getAminoAcidDBSize(), subMat, par.gapOpen, par.gapExtend, false);
    DistanceCalculator globalAliStat;
    if (par.globalAlignment) {
        globalAliStat.prepareGlobalAliParam(*subMat);
    }


    Debug(Debug::INFO) << "Result database: " << par.db4 << "\n";


    size_t totalMemory = Util::getTotalSystemMemory();
    size_t flushSize = 100000000;
    if (totalMemory > resultReader.getDataSize()) {
        flushSize = resultReader.getSize();
    }
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbSize) / static_cast<double>(flushSize)));
    for (size_t i = 0; i < iterations; i++) {
        size_t start = dbFrom + (i * flushSize);
        size_t bucketSize = std::min(dbSize - (i * flushSize), flushSize);
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            char buffer[1024+32768];

            std::string resultBuffer;
            resultBuffer.reserve(1000000);

            std::vector<Matcher::result_t> alnResults;
            alnResults.reserve(300);
            std::vector<hit_t> shortResults;
            shortResults.reserve(300);

#pragma omp for schedule(dynamic, 1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                Debug::printProgress(id);

                char *data = resultReader.getData(id);
                size_t queryKey = resultReader.getDbKey(id);
                unsigned int queryId = qdbr.getId(queryKey);
                char *querySeq = qdbr.getData(queryId);
                int queryLen = std::max(0, static_cast<int>(qdbr.getSeqLens(queryId)) - 2);

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
                    int dbLen = std::max(0, static_cast<int>(tdbr->getSeqLens(targetId)) - 2);

                    float queryLength = static_cast<float>(queryLen);
                    float targetLength = static_cast<float>(dbLen);
                    if(Util::canBeCovered(par.covThr, par.covMode, queryLength, targetLength)==false){
                        continue;
                    }
                    short diagonal = results[entryIdx].diagonal;
                    unsigned short distanceToDiagonal = abs(diagonal);
                    unsigned int diagonalLen = 0;
                    unsigned int distance = 0;
                    DistanceCalculator::LocalAlignment alignment;
                    if (diagonal >= 0 && distanceToDiagonal < queryLen) {
                        diagonalLen = std::min(dbLen, queryLen - distanceToDiagonal);
                        if (par.rescoreMode == Parameters::RESCORE_MODE_HAMMING) {
                            distance = DistanceCalculator::computeHammingDistance(
                                    querySeq + distanceToDiagonal, targetSeq, diagonalLen);
                        } else if (par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION) {
                            distance = DistanceCalculator::computeSubstitutionDistance(
                                    querySeq + distanceToDiagonal, targetSeq, diagonalLen, fastMatrix.matrix, par.globalAlignment);
                        } else if (par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT) {
                            alignment = DistanceCalculator::computeSubstitutionStartEndDistance(
                                    querySeq + distanceToDiagonal, targetSeq, diagonalLen, fastMatrix.matrix);
                            distance = alignment.score;
                        }
                    } else if (diagonal < 0 && distanceToDiagonal < dbLen) {
                        diagonalLen = std::min(dbLen - distanceToDiagonal, queryLen);
                        if (par.rescoreMode == Parameters::RESCORE_MODE_HAMMING) {
                            distance = DistanceCalculator::computeHammingDistance(
                                    querySeq, targetSeq + distanceToDiagonal, diagonalLen);
                        } else if (par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION) {
                            distance = DistanceCalculator::computeSubstitutionDistance(
                                    querySeq, targetSeq + distanceToDiagonal, diagonalLen, fastMatrix.matrix, par.globalAlignment);
                        } else if (par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT) {
                            alignment = DistanceCalculator::computeSubstitutionStartEndDistance(
                                    querySeq, targetSeq + distanceToDiagonal, diagonalLen, fastMatrix.matrix);
                            distance = alignment.score;
                        }
                    }

                    double seqId = 0;
                    double evalue = 0.0;
                    float targetCov = static_cast<float>(diagonalLen) / static_cast<float>(dbLen);
                    float queryCov = static_cast<float>(diagonalLen) / static_cast<float>(queryLen);

                    Matcher::result_t result;
                    if (par.rescoreMode == Parameters::RESCORE_MODE_HAMMING) {
                        int idCnt = (static_cast<float>(diagonalLen) - static_cast<float>(distance));
                        seqId = Util::computeSeqId(par.seqIdMode, idCnt, queryLen, dbLen, diagonalLen);
                    } else if (par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION || par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT) {
                        //seqId = exp(static_cast<float>(distance) / static_cast<float>(diagonalLen));
                        if (par.globalAlignment) {
                            // FIXME: value is never written to file
                            seqId = globalAliStat.getPvalGlobalAli((float)distance, diagonalLen);
                        } else {
                            evalue = evaluer.computeEvalue(distance, queryLen);
                            int bitScore = static_cast<short>(evaluer.computeBitScore(distance)+0.5);

                            if (par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT) {
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
//                                int qAlnLen = std::max(qEndPos - qStartPos, static_cast<int>(1));
//                                int dbAlnLen = std::max(dbEndPos - dbStartPos, static_cast<int>(1));
//                                seqId = (alignment.score1 / static_cast<float>(std::max(qAlnLength, dbAlnLength)))  * 0.1656 + 0.1141;

                                // compute seq.id if hit fulfills e-value but not by seqId criteria
                                if (evalue <= par.evalThr) {
                                    int idCnt = 0;
                                    for (int i = qStartPos; i <= qEndPos; i++) {
                                        idCnt += (querySeq[i] == targetSeq[dbStartPos+(i-qStartPos)]) ? 1 : 0;
                                    }
                                    unsigned int alnLength = Matcher::computeAlnLength(qStartPos, qEndPos, dbStartPos, dbEndPos);
                                    seqId = Util::computeSeqId(par.seqIdMode, idCnt, queryLen, dbLen, alnLength);
                                }

                                char *end = Itoa::i32toa_sse2(qEndPos-qStartPos, buffer);
                                size_t len = end - buffer;
                                std::string backtrace(buffer, len - 1);
                                backtrace.push_back('M');
                                queryCov = SmithWaterman::computeCov(qStartPos, qEndPos, queryLen);
                                targetCov = SmithWaterman::computeCov(dbStartPos, dbEndPos, dbLen);

                                result = Matcher::result_t(results[entryIdx].seqId, bitScore, queryCov, targetCov, seqId, evalue, alnLen,
                                                           qStartPos, qEndPos, queryLen, dbStartPos, dbEndPos, dbLen, backtrace);
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
                    bool hasToFilter = (par.filterHits == true && currScorePerCol >= scorePerColThr);
                    if (isIdentity || hasToFilter || (hasCov && hasSeqId && hasEvalue)) {
                        if (par.rescoreMode == Parameters::RESCORE_MODE_ALIGNMENT) {
                            alnResults.emplace_back(result);
                        } else if(par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION) {
                            hit_t hit;
                            hit.seqId = results[entryIdx].seqId;
                            hit.pScore = evalue;
                            hit.diagonal = diagonal;
                            shortResults.emplace_back(hit);
                        } else {
                            hit_t hit;
                            hit.seqId = results[entryIdx].seqId;
                            hit.pScore = seqId;
                            hit.diagonal = diagonal;
                            shortResults.emplace_back(hit);
                        }
                    }
                }

                if (par.sortResults > 0 && alnResults.size() > 1) {
                    std::sort(alnResults.begin(), alnResults.end(), Matcher::compareHits);
                }
                for (size_t i = 0; i < alnResults.size(); ++i) {
                    size_t len = Matcher::resultToBuffer(buffer, alnResults[i], true, false);
                    resultBuffer.append(buffer, len);
                }

                if (par.sortResults > 0 && shortResults.size() > 1) {
                    std::sort(shortResults.begin(), shortResults.end(), hit_t::compareHitsByPValueAndId);
                }
                for (size_t i = 0; i < shortResults.size(); ++i) {
                    if (par.rescoreMode == Parameters::RESCORE_MODE_SUBSTITUTION) {
                        size_t len = snprintf(buffer, 100, "%u\t%.3e\t%d\n", shortResults[i].seqId, shortResults[i].pScore, shortResults[i].diagonal);
                        resultBuffer.append(buffer, len);
                    } else {
                        size_t len = snprintf(buffer, 100, "%u\t%.2f\t%d\n", shortResults[i].seqId, shortResults[i].pScore, shortResults[i].diagonal);
                        resultBuffer.append(buffer, len);
                    }
                }

                resultWriter.writeData(resultBuffer.c_str(), resultBuffer.length(), queryKey, thread_idx);
                resultBuffer.clear();
                shortResults.clear();
                alnResults.clear();
            }
        }
        resultReader.remapData();
    }
    Debug(Debug::INFO) << "\nDone.\n";
    qdbr.close();
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }

    delete[] fastMatrix.matrix;
    delete[] fastMatrix.matrixData;
    delete subMat;
    return 0;
}

int rescorediagonal(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);


    Debug(Debug::INFO) << "Prefilter database: " << par.db3 << "\n";
    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str());
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
#ifdef HAVE_MPI
    size_t dbFrom = 0;
    size_t dbSize = 0;

    Util::decomposeDomainByAminoAcid(resultReader.getAminoAcidDBSize(), resultReader.getSeqLens(), resultReader.getSize(),
                                     MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db4, par.db4Index, MMseqsMPI::rank);

    DBWriter resultWriter(tmpOutput.first.c_str(), tmpOutput.second.c_str(), par.threads);
    resultWriter.open();
    int status = doRescorediagonal(par, resultWriter, resultReader, dbFrom, dbSize);
    resultWriter.close();

    MPI_Barrier(MPI_COMM_WORLD);
    if(MMseqsMPI::rank == 0) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for(unsigned int proc = 0; proc < MMseqsMPI::numProc; ++proc){
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(par.db4, par.db4Index, proc);
            splitFiles.push_back(std::make_pair(tmpFile.first,  tmpFile.second));
        }
        DBWriter::mergeResults(par.db4, par.db4Index, splitFiles);
    }
#else
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads);
    resultWriter.open();
    int status = doRescorediagonal(par, resultWriter, resultReader, 0, resultReader.getSize());
    resultWriter.close();

#endif




    resultReader.close();



    return status;
}


