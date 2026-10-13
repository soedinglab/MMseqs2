// The "fwbw" command. The aligner itself lives in FwbwInternal.cpp, compiled once per
// instruction set; this file stays free of simd.h so it can be built with the baseline flags.
#include "Fwbw.h"
#include "Alignment.h"
#include "Command.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "FastSort.h"
#include "Matcher.h"
#include "Parameters.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "Util.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

int fwbw(int argc, const char **argv, const Command &command) {
    //Prepare the parameters & DB
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    DBReader<DBKeyType> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX);
    qdbr.open(DBReader<DBKeyType>::NOSORT);
    DBReader<DBKeyType> tdbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX);
    tdbr.open(DBReader<DBKeyType>::NOSORT);
    DBReader<DBKeyType> alnRes (par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX);
    alnRes.open(DBReader<DBKeyType>::LINEAR_ACCCESS);

    DBWriter fwbwAlnWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    fwbwAlnWriter.open();
    const int querySeqType = qdbr.getDbtype();
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        Debug(Debug::ERROR) << "Invalid datatype. Nucleotide.\n";
        EXIT(EXIT_FAILURE);
    }
    SubstitutionMatrix subMat = SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias); // Check : par.scoreBias = 0.0
    
    const size_t flushSize = 100000000;
    size_t iterations = static_cast<int>(ceil(static_cast<double>(alnRes.getSize()) / static_cast<double>(flushSize)));
    Debug(Debug::INFO) << "Processing " << iterations << " iterations\n";
    Debug(Debug::INFO) << "SIMD architecture: " << FwBwAligner::simdArchitecture() << "\n";
    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(alnRes.getSize() - (i * flushSize), flushSize);
        Debug::Progress progress(bucketSize);
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            size_t length = par.blocklen;
            Sequence qSeq(par.maxSeqLen, qdbr.getDbtype(), &subMat, 0, false, false);
            Sequence dbSeq(par.maxSeqLen, tdbr.getDbtype(), &subMat, 0, false, false);
            
            const size_t assignSeqLen = FwBwAligner::suggestedCapacity();
            FwBwAligner fwbwaligner(subMat, -par.fwbwGapopen, -par.fwbwGapextend, par.temperature, par.mact, assignSeqLen, assignSeqLen, length, true);
            char entrybuffer[1024 + 32768*4];
            std::string alnResultsOutString;
            char buffer[1024 + 32768*4];
            std::vector<Matcher::result_t> localFwbwResults;
            localFwbwResults.reserve(300);

#pragma omp for schedule(dynamic,1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();
                DBKeyType key = alnRes.getDbKey(id);
                const size_t queryId = qdbr.getId(key);
                char *alnData = alnRes.getData(id, thread_idx);
                localFwbwResults.clear();

                const char* querySeq = qdbr.getData(queryId, thread_idx);
                size_t queryLen = qdbr.getSeqLen(queryId);

                qSeq.mapSequence(queryId, key, querySeq, queryLen);
                fwbwaligner.initProfile(qSeq.numSequence, queryLen);
                fwbwAlnWriter.writeStart(thread_idx);

                while (*alnData != '\0'){
                    Util::parseKey(alnData, entrybuffer);
                    DBKeyType targetKey = Util::fast_atoi<DBKeyType>(entrybuffer);
                    const size_t targetId = tdbr.getId(targetKey);
                    const char* targetSeq = tdbr.getData(targetId, thread_idx);
                    size_t targetLen = tdbr.getSeqLen(targetId);

                    dbSeq.mapSequence(targetId, targetKey, targetSeq, targetLen);
                    //Init target & Resizing memory
                    fwbwaligner.initAlignment(dbSeq.numSequence, targetLen, queryLen); 
                    // backtrace modes 2 (semi-global) and 3 (global) exist but are not exposed
                    fwbwaligner.runFwBw(true, par.fwbwBacktraceMode);

                    // Map s_align values to result_t 
                    FwBwAligner::s_align fwbwAlignment = fwbwaligner.getFwbwAlnResult();
                    
                    float qcov = fwbwAlignment.qCov;
                    float dbcov = fwbwAlignment.dbCov;
                    float evalue = 0;
                    const int score = fwbwAlignment.score2;
                    const unsigned int qStartPos = fwbwAlignment.qStartPos1;
                    const unsigned int dbStartPos = fwbwAlignment.dbStartPos1;
                    const unsigned int qEndPos = fwbwAlignment.qEndPos1;
                    const unsigned int dbEndPos = fwbwAlignment.dbEndPos1;
                    std::string backtrace = fwbwAlignment.cigar;
                    unsigned int alnLength = backtrace.size();
                    float seqId = Util::computeSeqId(par.seqIdMode, fwbwAlignment.identicalAACnt, queryLen, targetLen, alnLength);
                    Matcher::result_t res = Matcher::result_t(targetKey, score, qcov, dbcov, seqId, evalue, alnLength, qStartPos, qEndPos, queryLen, dbStartPos, dbEndPos, targetLen, backtrace);
                    if (Alignment::checkCriteria(res, 0, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                        localFwbwResults.emplace_back(res);
                    }
                    alnData = Util::skipLine(alnData);
                }

                // sort local results. They will currently be sorted by first fwbwscore, then targetlen, then by targetkey.
                SORT_SERIAL(localFwbwResults.begin(), localFwbwResults.end(), Matcher::compareHits);
                for (size_t result = 0; result < localFwbwResults.size(); result++) {
                    size_t len = Matcher::resultToBuffer(buffer, localFwbwResults[result], true, true);
                    alnResultsOutString.append(buffer, len);
                }
                fwbwAlnWriter.writeData(alnResultsOutString.c_str(), alnResultsOutString.length(), alnRes.getDbKey(id), thread_idx);
                alnResultsOutString.clear();
                localFwbwResults.clear();            
            }
        }
        alnRes.remapData();
        
    }
    fwbwAlnWriter.close();
    alnRes.close();
    qdbr.close();
    tdbr.close();

    return EXIT_SUCCESS;
}