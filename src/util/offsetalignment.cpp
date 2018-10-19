#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Orf.h"
#include "AlignmentSymmetry.h"
#include "Timer.h"
#include "HeaderIdReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

void updateOffset(char* data, std::vector<Matcher::result_t> &results,
                  const Orf::SequenceLocation *qloc,
                  HeaderIdReader& tHeaderDbr,
                  bool nucleotide) {
    size_t startPos = results.size();
    Matcher::readAlignmentResults(results, data, true);
    size_t endPos = results.size();
    for (size_t i = startPos; i < endPos; i++) {
        Matcher::result_t &res = results[i];
        if (qloc == NULL) {
            size_t targetId = tHeaderDbr.getReader()->getId(res.dbKey);
            char *header = tHeaderDbr.getReader()->getData(targetId);
            Orf::SequenceLocation tloc = Orf::parseOrfHeader(header);
            res.dbKey = tloc.id;
            int dbStartPos = (nucleotide) ? res.dbStartPos : res.dbStartPos * 3;
            res.dbStartPos = tloc.from + dbStartPos;
            int dbEndPos = (nucleotide) ? res.dbEndPos + 1 : (res.dbEndPos + 1) * 3;
            res.dbEndPos   = tloc.from + dbEndPos;

            if (tloc.strand == Orf::STRAND_MINUS) {
                int start = res.dbStartPos;
                res.dbStartPos = res.dbEndPos;
                res.dbEndPos = start;
            }
            res.dbLen = (nucleotide) ? res.dbLen :  res.dbLen * 3;
        } else {
            int qStartPos = (nucleotide) ? res.qStartPos  : res.qStartPos * 3;
            int qEndPos = (nucleotide) ? (res.qEndPos+1)   : (res.qEndPos+1) * 3;

            res.qStartPos = qloc->from + qStartPos;
            res.qEndPos = qloc->from + qEndPos;

            if (qloc->strand == Orf::STRAND_MINUS) {
                int start = res.qStartPos;
                res.qStartPos = res.qEndPos;
                res.qEndPos = start;
            }
            res.qLen = (nucleotide) ? res.qLen : res.qLen * 3;
        }
    }
}

int offsetalignment(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    const int queryDbType = DBReader<unsigned int>::parseDbType(par.db1.c_str());
    const int targetDbType = DBReader<unsigned int>::parseDbType(par.db2.c_str());
    if (queryDbType == -1 || targetDbType == -1) {
        Debug(Debug::ERROR)
                << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        return EXIT_FAILURE;
    }

    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    Debug(Debug::INFO) << "Query database: " << par.db1 << "\n";
    HeaderIdReader qHeaderDbr(par.db1.c_str(), touch);

    Debug(Debug::INFO) << "Target database: " << par.db2 << "\n";
    HeaderIdReader tHeaderDbr(par.db2.c_str(), touch);



//    DBReader<unsigned int> qHeaderDbr(par.hdr1.c_str(), par.hdr1Index.c_str());
//    qHeaderDbr.open(DBReader<unsigned int>::NOSORT);

//    DBReader<unsigned int> tHeaderDbr(par.hdr2.c_str(), par.hdr2Index.c_str());
//    tHeaderDbr.open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << "Result database: " << par.db3 << "\n";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str());
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
    if (queryDbType == Sequence::NUCLEOTIDES) {
        Timer timer;
        Debug(Debug::INFO) << "Computing ORF lookup...\n";
        unsigned int maxOrfKey = alnDbr.getLastKey();
        unsigned int *orfLookup = new unsigned int[maxOrfKey + 2]();
#pragma omp parallel for schedule(dynamic, 10) num_threads(localThreads) reduction(max:maxContigKey)
        for (size_t i = 0; i <= maxOrfKey; ++i) {
            size_t queryId = qHeaderDbr.getReader()->getId(i);
            char *header = qHeaderDbr.getReader()->getData(queryId);
            Orf::SequenceLocation qloc = Orf::parseOrfHeader(header);
            orfLookup[i] = qloc.id;
            maxContigKey = std::max(maxContigKey, qloc.id);
        }

        Debug(Debug::INFO) << "Computing contig offsets...\n";
        unsigned int *contigSizes = new unsigned int[maxContigKey + 2]();
        contigExists = new char[maxContigKey + 1]();
#pragma omp parallel for schedule(static) num_threads(localThreads)
        for (size_t i = 0; i <= maxOrfKey; ++i) {
            __sync_fetch_and_add(&(contigSizes[orfLookup[i]]), 1);
            contigExists[orfLookup[i]] = 1;
        }
        contigOffsets = contigSizes;
        AlignmentSymmetry::computeOffsetFromCounts(contigOffsets, maxContigKey + 1);

        Debug(Debug::INFO) << "Computing contig lookup...\n";
        contigLookup = new unsigned int[maxOrfKey + 2]();
#pragma omp parallel for schedule(static) num_threads(localThreads)
        for (size_t i = 0; i <= maxOrfKey; ++i) {
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

    Debug(Debug::INFO) << "Writing results to: " << par.db4 << "\n";
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads);
    resultWriter.open();
    bool isNucl = queryDbType == Sequence::NUCLEOTIDES && targetDbType == Sequence::NUCLEOTIDES;
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[1024];

        std::string ss;
        ss.reserve(1024);

        std::vector<Matcher::result_t> results;
        results.reserve(300);

        size_t entryCount = alnDbr.getSize();
        if (queryDbType == Sequence::NUCLEOTIDES) {
            entryCount = maxContigKey + 1;
        }

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < entryCount; ++i) {
            Debug::printProgress(i);

            unsigned int queryKey;
            if (queryDbType == Sequence::NUCLEOTIDES) {
                queryKey = i;
                if (contigExists[i] == 0) {
                    continue;
                }
                unsigned int *orfKeys = &contigLookup[contigOffsets[i]];
                size_t orfCount = contigOffsets[i + 1] - contigOffsets[i];
                for (unsigned int j = 0; j < orfCount; ++j) {
                    unsigned int orfKey = orfKeys[j];
                    size_t orfId = alnDbr.getId(orfKey);
                    char *data = alnDbr.getData(orfId);
                    size_t queryId = qHeaderDbr.getReader()->getId(i);
                    char *header = qHeaderDbr.getReader()->getData(queryId);
                    Orf::SequenceLocation qloc = Orf::parseOrfHeader(header);
                    updateOffset(data, results, &qloc, tHeaderDbr, isNucl);
                }
            }
            if (targetDbType == Sequence::NUCLEOTIDES) {
                queryKey = alnDbr.getDbKey(i);
                char *data = alnDbr.getData(i);
                updateOffset(data, results, NULL, tHeaderDbr, isNucl);
            }
            std::stable_sort(results.begin(), results.end(), Matcher::compareHits);
            for(size_t i = 0; i < results.size(); i++){
                Matcher::result_t &res = results[i];
                bool hasBacktrace = (res.backtrace.size() > 0);
                size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false);
                ss.append(buffer, len);
            }
            resultWriter.writeData(ss.c_str(), ss.length(), queryKey, thread_idx);
            ss.clear();
            results.clear();
        }
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

    alnDbr.close();

    return EXIT_SUCCESS;
}

