#include "DistanceCalculator.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "QueryMatcher.h"
#include "QueryMatcher.h"
#include "NucleotideMatrix.h"

#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif



int alignall(int argc, const char **argv, const Command &command) {
    Debug(Debug::INFO) << "Rescore diagonals.\n";
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> *qdbr = NULL;

    Debug(Debug::INFO) << "Query  file: " << par.db1 << "\n";
    qdbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str());
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

    qdbr->readMmapedDataInMemory();
    Debug(Debug::INFO) << "Target  file: " << par.db2 << "\n";
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        tdbr->readMmapedDataInMemory();
    }

    EvalueComputation evaluer(tdbr->getAminoAcidDBSize(), subMat, par.gapOpen, par.gapExtend, false);

    Debug(Debug::INFO) << "Prefilter database: " << par.db3 << "\n";
    DBReader<unsigned int> dbr_res(par.db3.c_str(), par.db3Index.c_str());
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
            Matcher matcher(querySeqType, par.maxSeqLen, subMat,
                            &evaluer, par.compBiasCorrection,
                            par.gapOpen, par.gapExtend);
            Sequence query(par.maxSeqLen, querySeqType, subMat,
                           par.kmerSize, par.spacedKmer, par.compBiasCorrection);
            Sequence target(par.maxSeqLen, querySeqType, subMat,
                           par.kmerSize, par.spacedKmer, par.compBiasCorrection);
#pragma omp for schedule(dynamic, 1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                Debug::printProgress(id);
                char buffer[1024 + 32768];

                unsigned int thread_idx = 0;
#ifdef OPENMP
                thread_idx = (unsigned int) omp_get_thread_num();
#endif
                char *data = dbr_res.getData(id);
                unsigned int qId = qdbr->getId(dbr_res.getDbKey(id));
                unsigned int qKey = qdbr->getDbKey(qId);
                std::vector<unsigned int> results;
                while (*data != '\0') {
                    Util::parseKey(data, buffer);
                    const unsigned int key = (unsigned int) strtoul(buffer, NULL, 10);
                    results.push_back(key);
                    data = Util::skipLine(data);
                }
                    resultWriter.writeStart(thread_idx);
                for (size_t entryIdx1 = 0; entryIdx1 < results.size(); entryIdx1++) {
                    unsigned int queryId = tdbr->getId(results[entryIdx1]);
                    unsigned int queryKey = tdbr->getDbKey(queryId);
                    char *querySeq = tdbr->getData(queryId);
                    query.mapSequence(id, queryKey, querySeq);
                    matcher.initQuery(&query);
                    char * tmpBuff = Itoa::u32toa_sse2((uint32_t) queryKey, buffer);
                    *(tmpBuff-1) = '\t';
                    const unsigned int queryIdLen = tmpBuff - buffer;

                    for (size_t entryIdx = 0; entryIdx < results.size(); entryIdx++) {
                        unsigned int targetId = tdbr->getId(results[entryIdx]);
                        unsigned int targetKey = tdbr->getDbKey(targetId);
                        char *targetSeq = tdbr->getData(targetId);
                        target.mapSequence(id, targetKey, targetSeq);

                        const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB)) ? true : false;
                        if(Util::canBeCovered(par.covThr, par.covMode, query.L, target.L)==false){
                            continue;
                        }
                        double seqId = 0;
                        double evalue = 0.0;
                        Matcher::result_t result = matcher.getSWResult(&target, INT_MAX, par.covMode, par.covThr, FLT_MAX,
                                                                       par.alignmentMode, par.seqIdMode, isIdentity);
                        // query/target cov mode
                        bool hasCov = Util::hasCoverage(par.covThr, par.covMode, result.qcov, result.dbcov);
                        // --min-seq-id
                        bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                        bool hasEvalue = (evalue <= par.evalThr);
                        // --filter-hits
                        if (isIdentity || (hasCov && hasSeqId && hasEvalue)) {
                            size_t len = Matcher::resultToBuffer(tmpBuff, result, true, false);
                            resultWriter.writeAdd(buffer, queryIdLen + len, thread_idx);
                        }
                    }
                }
                resultWriter.writeEnd(qKey, thread_idx, true);
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

    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }

    return EXIT_SUCCESS;
}


