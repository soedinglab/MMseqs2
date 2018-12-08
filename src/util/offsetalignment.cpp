#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Orf.h"
#include "AlignmentSymmetry.h"
#include "Timer.h"
#include "IndexReader.h"

#ifdef OPENMP
#include <omp.h>
#endif


void chainAlignmentHits(std::vector<Matcher::result_t> &results, std::vector<Matcher::result_t> &tmp) {
    if(results.size() > 1){
        std::stable_sort(results.begin(), results.end(), Matcher::compareHitsByPos);
        int prevDiagonal = INT_MAX;
        Matcher::result_t currResults = results[0];
        for (size_t resIdx = 0; resIdx < results.size(); resIdx++) {
//            bool queryRev = (results[resIdx].qStartPos > results[resIdx].qEndPos);
//            bool targetRev = (results[resIdx].qStartPos > results[resIdx].qEndPos);

            int currDiagonal = results[resIdx].qStartPos  - results[resIdx].dbStartPos;
            int nextDiagonal = UINT_MAX;
            const bool isDifferentKey = (currResults.dbKey != results[resIdx].dbKey);
            const bool isLastElement  = (resIdx == results.size() - 1);
            if (isLastElement == false) {
                nextDiagonal = results[resIdx+1].qStartPos  - results[resIdx+1].dbStartPos;
            }

            const bool queryIsOverlapping  = currResults.qEndPos >= results[resIdx].qStartPos   && currResults.qEndPos <= results[resIdx].qEndPos;
            const bool targetIsOverlapping = currResults.dbEndPos >= results[resIdx].dbStartPos && currResults.dbEndPos <= results[resIdx].dbEndPos;
            const bool sameNextDiagonal = (currDiagonal == nextDiagonal);
            const bool samePrevDiagonal = (currDiagonal == prevDiagonal);
            if ( (sameNextDiagonal || samePrevDiagonal ) && queryIsOverlapping && targetIsOverlapping) {
                currResults.qStartPos = std::min(currResults.qStartPos, results[resIdx].qStartPos);
                currResults.qEndPos = std::max(currResults.qEndPos,  results[resIdx].qEndPos);
                currResults.dbStartPos  = std::min(currResults.dbStartPos, results[resIdx].dbStartPos);
                currResults.dbEndPos  = std::max(currResults.dbEndPos, results[resIdx].dbEndPos);
            }

            prevDiagonal = currDiagonal;
            bool isDifferentNextDiagonal = (nextDiagonal != currDiagonal);
            if(isDifferentKey || isLastElement || isDifferentNextDiagonal){
                tmp.push_back(currResults);
                if(isLastElement == false) {
                    currResults = results[resIdx + 1];
                }
            }
        }

    }
}

void updateOffset(char* data, std::vector<Matcher::result_t> &results, const Orf::SequenceLocation *qloc, IndexReader& tOrfDBr, bool isNuclNucl, int thread_idx) {
    size_t startPos = results.size();
    Matcher::readAlignmentResults(results, data, true);
    size_t endPos = results.size();
    for (size_t i = startPos; i < endPos; i++) {
        Matcher::result_t &res = results[i];
        if (isNuclNucl == true || qloc == NULL) {
            size_t targetId = tOrfDBr.headerReader->getId(res.dbKey);
            char *header = tOrfDBr.headerReader->getData(targetId, thread_idx);

            Orf::SequenceLocation tloc = Orf::parseOrfHeader(header);
            res.dbKey   = (tloc.id != UINT_MAX) ? tloc.id : res.dbKey;
            size_t from = (tloc.id != UINT_MAX) ? tloc.from : (tloc.strand == Orf::STRAND_MINUS) ? 0 : res.dbLen - 1;

            int dbStartPos = isNuclNucl ? res.dbStartPos : res.dbStartPos * 3;
            int dbEndPos   = isNuclNucl ? res.dbEndPos : res.dbEndPos * 3;

            if (tloc.strand == Orf::STRAND_MINUS && tloc.id != UINT_MAX) {
                res.dbStartPos = from - dbStartPos;
                res.dbEndPos   = from - dbEndPos;
                if(isNuclNucl == false){
                    res.dbEndPos = res.dbEndPos - 2;
                }
            } else {
                res.dbStartPos = from + dbStartPos;
                res.dbEndPos   = from + dbEndPos;
                if(isNuclNucl == false){
                    res.dbEndPos = res.dbEndPos + 2;
                }
            }
        }
        if (qloc != NULL) {
            int qStartPos = isNuclNucl ? res.qStartPos : res.qStartPos * 3;
            int qEndPos   = isNuclNucl ? res.qEndPos : res.qEndPos * 3;

            size_t from = (qloc->id != UINT_MAX) ? qloc->from : (qloc->strand == Orf::STRAND_MINUS) ? 0 : res.qLen - 1;

            if (qloc->strand == Orf::STRAND_MINUS && qloc->id != UINT_MAX) {
                res.qStartPos  = from - qStartPos;
                res.qEndPos    = from - qEndPos;
            } else {
                res.qStartPos  = from + qStartPos;
                res.qEndPos    = from + qEndPos;
            }
        }
    }
}

void updateLengths(std::vector<Matcher::result_t> &results, unsigned int qSourceLen, IndexReader* tSourceDbr) {
    for (size_t i = 0; i < results.size(); ++i) {
        Matcher::result_t &res = results[i];
        if (qSourceLen != UINT_MAX) {
            res.qLen = qSourceLen;
        }
        if (tSourceDbr != NULL) {
            size_t targetId = tSourceDbr->sequenceReader->getId(res.dbKey);
            res.dbLen = std::max(tSourceDbr->sequenceReader->getSeqLens(targetId), static_cast<size_t>(2)) - 2;
        }
    }
}

int offsetalignment(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 6);

    const bool touch = par.preloadMode != Parameters::PRELOAD_MODE_MMAP;


    int queryDbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    int targetDbType = DBReader<unsigned int>::parseDbType(par.db3.c_str());

    Debug(Debug::INFO) << "Query database: " << par.db2 << "\n";
    IndexReader qOrfDbr(par.db2.c_str(), par.threads, IndexReader::NEED_ALT_HEADERS, touch);
    if (queryDbType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        return EXIT_FAILURE;
    }
    const bool queryNucl = Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES);
    IndexReader *qSourceDbr = NULL;
    if (queryNucl) {
        Debug(Debug::INFO) << "Source Query database: " << par.db1 << "\n";
        qSourceDbr = new IndexReader(par.db1.c_str(), par.threads, IndexReader::NEED_SEQ_INDEX, touch);
    }

    Debug(Debug::INFO) << "Target database: " << par.db4 << "\n";
    IndexReader tOrfDbr(par.db4.c_str(), par.threads, IndexReader::NEED_ALT_HEADERS, touch);
    if (targetDbType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        return EXIT_FAILURE;
    }
    const bool targetNucl = Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES);
    IndexReader *tSourceDbr = NULL;
    if (targetNucl) {
        Debug(Debug::INFO) << "Source Target database: " << par.db3 << "\n";
        tSourceDbr = new IndexReader(par.db3.c_str(), par.threads, IndexReader::NEED_SEQ_INDEX, touch);
    }

    Debug(Debug::INFO) << "Result database: " << par.db5 << "\n";
    DBReader<unsigned int> alnDbr(par.db5.c_str(), par.db5Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

#ifdef OPENMP
    unsigned int totalThreads = par.threads;
#else
    unsigned int totalThreads = 1;
#endif

    unsigned int localThreads = totalThreads;
    if (alnDbr.getSize() <= totalThreads) {
        localThreads = alnDbr.getSize();
    }

    // Compute mapping from contig -> orf[] from orf[]->contig in headers
    unsigned int *contigLookup = NULL;
    unsigned int *contigOffsets = NULL;
    char *contigExists = NULL;
    unsigned int maxContigKey = 0;
    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
        Timer timer;
        Debug(Debug::INFO) << "Computing ORF lookup...\n";
        unsigned int maxOrfKey = alnDbr.getLastKey();
        unsigned int *orfLookup = new unsigned int[maxOrfKey + 2]();
#pragma omp parallel num_threads(localThreads)
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 10)
            for (size_t i = 0; i <= maxOrfKey; ++i) {
                size_t queryId = qOrfDbr.headerReader->getId(i);
                if (queryId == UINT_MAX) {
                    orfLookup[i] = UINT_MAX;
                    continue;
                }
                unsigned int queryKey = qOrfDbr.headerReader->getDbKey(queryId);
                char *header = qOrfDbr.headerReader->getData(queryId, thread_idx);
                Orf::SequenceLocation qloc = Orf::parseOrfHeader(header);
                unsigned int id = (qloc.id != UINT_MAX) ? qloc.id : queryKey;
                orfLookup[i] = id;
            }
        }
        Debug(Debug::INFO) << "Computing contig offsets...\n";
        maxContigKey = qSourceDbr->sequenceReader->getLastKey();
        unsigned int *contigSizes = new unsigned int[maxContigKey + 2]();
#pragma omp parallel for schedule(static) num_threads(localThreads)
        for (size_t i = 0; i <= maxOrfKey ; ++i) {
            if(orfLookup[i] == UINT_MAX){
                continue;
            }
            __sync_fetch_and_add(&(contigSizes[orfLookup[i]]), 1);
        }
        contigOffsets = contigSizes;

        AlignmentSymmetry::computeOffsetFromCounts(contigOffsets, maxContigKey + 1);

        contigExists = new char[maxContigKey + 1]();
#pragma omp parallel for schedule(static) num_threads(localThreads)
        for (size_t i = 0; i < qSourceDbr->sequenceReader->getSize(); ++i) {
            contigExists[qSourceDbr->sequenceReader->getDbKey(i)] = 1;
        }

        Debug(Debug::INFO) << "Computing contig lookup...\n";
        contigLookup = new unsigned int[maxOrfKey + 2]();
#pragma omp parallel for schedule(static) num_threads(localThreads)
        for (size_t i = 0; i <= maxOrfKey; ++i) {
            if(orfLookup[i] == UINT_MAX){
                continue;
            }
            size_t offset = __sync_fetch_and_add(&(contigOffsets[orfLookup[i]]), 1);
            contigLookup[offset] = i;
        }
        delete[] orfLookup;

        for (unsigned int i = maxContigKey + 1; i > 0; --i) {
            contigOffsets[i] = contigOffsets[i - 1];
        }
        contigOffsets[0] = 0;
        Debug(Debug::INFO) << "Time for contig lookup: " << timer.lap() << "\n";
    }

    Debug(Debug::INFO) << "Writing results to: " << par.db6 << "\n";
    DBWriter resultWriter(par.db6.c_str(), par.db6Index.c_str(), localThreads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();

    const bool isNuclNucl = queryNucl && targetNucl;

#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char * buffer = new char[65536];

        std::string ss;
        ss.reserve(1024);

        std::vector<Matcher::result_t> results;
        std::vector<Matcher::result_t> tmp;
        results.reserve(300);
        tmp.reserve(300);
        size_t entryCount = alnDbr.getSize();
        if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
            entryCount = maxContigKey + 1;
        }

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < entryCount; ++i) {
            Debug::printProgress(i);
            unsigned int queryKey=UINT_MAX;
            if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
                queryKey = i;
                if (contigExists[i] == 0) {
                    continue;
                }
                unsigned int *orfKeys = &contigLookup[contigOffsets[i]];
                size_t orfCount = contigOffsets[i + 1] - contigOffsets[i];
                for (unsigned int j = 0; j < orfCount; ++j) {
                    unsigned int orfKey = orfKeys[j];
                    size_t orfId = alnDbr.getId(orfKey);
                    char *data = alnDbr.getData(orfId, thread_idx);
                    size_t queryId = qOrfDbr.headerReader->getId(orfKey);
                    char *header = qOrfDbr.headerReader->getData(queryId, thread_idx);
                    Orf::SequenceLocation qloc = Orf::parseOrfHeader(header);
                    updateOffset(data, results, &qloc, tOrfDbr, isNuclNucl, thread_idx);
                }
            } else if (Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
                queryKey = alnDbr.getDbKey(i);
                char *data = alnDbr.getData(i, thread_idx);
                updateOffset(data, results, NULL, tOrfDbr, isNuclNucl, thread_idx);
            }
            unsigned int qLen = UINT_MAX;
            if (qSourceDbr != NULL) {
                size_t queryId = qSourceDbr->sequenceReader->getId(queryKey);
                qLen = std::max(qSourceDbr->sequenceReader->getSeqLens(queryId), static_cast<size_t>(2)) - 2;
            }
            updateLengths(results, qLen, tSourceDbr);
            if(par.chainAlignment == false){
                std::stable_sort(results.begin(), results.end(), Matcher::compareHits);
                for(size_t i = 0; i < results.size(); i++){
                    Matcher::result_t &res = results[i];
                    bool hasBacktrace = (res.backtrace.size() > 0);
                    size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
                    ss.append(buffer, len);
                }
                resultWriter.writeData(ss.c_str(), ss.length(), queryKey, thread_idx);
            } else if(par.chainAlignment == true){
                chainAlignmentHits(results, tmp);
                for(size_t i = 0; i < tmp.size(); i++){
                    Matcher::result_t &res = tmp[i];
                    bool hasBacktrace = (res.backtrace.size() > 0);
                    size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
                    ss.append(buffer, len);
                }
                resultWriter.writeData(ss.c_str(), ss.length(), queryKey, thread_idx);
            }
            ss.clear();
            results.clear();
            tmp.clear();
        }
        delete[] buffer;
    }
    Debug(Debug::INFO) << "\n";
    resultWriter.close();

    if (contigLookup != NULL) {
        delete[] contigLookup;
    }

    if (contigOffsets != NULL) {
        delete[] contigOffsets;
    }

    if (contigExists != NULL) {
        delete[] contigExists;
    }

    if(tSourceDbr != NULL){
        delete tSourceDbr;
    }

    if(qSourceDbr != NULL){
        delete qSourceDbr;
    }
    alnDbr.close();

    return EXIT_SUCCESS;
}
