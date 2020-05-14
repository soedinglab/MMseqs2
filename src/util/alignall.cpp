#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "Alignment.h"
#include "itoa.h"

#ifdef OPENMP
#include <omp.h>
#endif

int alignall(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.overrideParameterDescription(par.PARAM_ALIGNMENT_MODE, "How to compute the alignment:\n0: automatic\n1: only score and end_pos\n2: also start_pos and cov\n3: also seq.id", NULL, 0);
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED) {
        Debug(Debug::ERROR) << "Use rescorediagonal for ungapped alignment mode.\n";
        EXIT(EXIT_FAILURE);
    }
    if (par.addBacktrace == true) {
        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }
    unsigned int swMode = Alignment::initSWMode(par.alignmentMode, par.covThr, par.seqIdThr);

    DBReader<unsigned int> tdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    tdbr.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        tdbr.readMmapedDataInMemory();
    }
    const int targetSeqType = tdbr.getDbtype();

    int gapOpen, gapExtend;
    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(targetSeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, par.scoreBias);
        gapOpen = par.gapOpen.nucleotides;
        gapExtend = par.gapExtend.nucleotides;
    } else {
        // keep score bias at 0.0 (improved ROC)
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, par.scoreBias);
        gapOpen = par.gapOpen.aminoacids;
        gapExtend = par.gapExtend.aminoacids;
    }

    DBReader<unsigned int> dbr_res(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    dbr_res.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    resultWriter.open();

    EvalueComputation evaluer(tdbr.getAminoAcidDBSize(), subMat, gapOpen, gapExtend);
    const size_t flushSize = 100000000;
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbr_res.getSize()) / static_cast<double>(flushSize)));

    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(dbr_res.getSize() - (i * flushSize), flushSize);
        Debug::Progress progress(bucketSize);
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif

            Matcher matcher(targetSeqType, par.maxSeqLen, subMat, &evaluer, par.compBiasCorrection, gapOpen, gapExtend, par.zdrop);

            Sequence query(par.maxSeqLen, targetSeqType, subMat, 0, false, par.compBiasCorrection);
            Sequence target(par.maxSeqLen, targetSeqType, subMat, 0, false, par.compBiasCorrection);

            char buffer[1024 + 32768];

            std::vector<unsigned int> results;
            results.reserve(300);

#pragma omp for schedule(dynamic, 1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                progress.updateProgress();

                const unsigned int key = dbr_res.getDbKey(id);
                char *data = dbr_res.getData(id, thread_idx);

                results.clear();
                while (*data != '\0') {
                    Util::parseKey(data, buffer);
                    const unsigned int key = (unsigned int) strtoul(buffer, NULL, 10);
                    results.push_back(key);
                    data = Util::skipLine(data);
                }

                resultWriter.writeStart(thread_idx);
                for (size_t entryIdx1 = 0; entryIdx1 < results.size(); entryIdx1++) {
                    const unsigned int queryId = tdbr.getId(results[entryIdx1]);
                    const unsigned int queryKey = tdbr.getDbKey(queryId);
                    char *querySeq = tdbr.getData(queryId, thread_idx);
                    query.mapSequence(queryId, queryKey, querySeq, tdbr.getSeqLen(queryId));
                    matcher.initQuery(&query);

                    char * tmpBuff = Itoa::u32toa_sse2((uint32_t) queryKey, buffer);
                    *(tmpBuff-1) = '\t';
                    const unsigned int queryIdLen = tmpBuff - buffer;

                    for (size_t entryIdx = 0; entryIdx < results.size(); entryIdx++) {
                        const unsigned int targetId = tdbr.getId(results[entryIdx]);
                        const unsigned int targetKey = tdbr.getDbKey(targetId);
                        char *targetSeq = tdbr.getData(targetId, thread_idx);
                        target.mapSequence(id, targetKey, targetSeq, tdbr.getSeqLen(targetId));

                        if (Util::canBeCovered(par.covThr, par.covMode, query.L, target.L) == false) {
                            continue;
                        }

                        const bool isIdentity = (queryId == targetId && par.includeIdentity) ? true : false;
                        Matcher::result_t result = matcher.getSWResult(&target, INT_MAX, false, par.covMode, par.covThr, par.evalThr,
                                                                       swMode, par.seqIdMode, isIdentity);
                        // checkCriteria and Util::canBeCovered always work together
                        if (Alignment::checkCriteria(result, isIdentity, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode, par.covThr)) {
                            size_t len = Matcher::resultToBuffer(tmpBuff, result, par.addBacktrace);
                            resultWriter.writeAdd(buffer, queryIdLen + len, thread_idx);
                        }
                    }
                }
                resultWriter.writeEnd(key, thread_idx);
            }
        }
        dbr_res.remapData();
    }
    resultWriter.close();
    dbr_res.close();
    delete subMat;
    tdbr.close();

    return EXIT_SUCCESS;
}


