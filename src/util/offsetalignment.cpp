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
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif


void chainAlignmentHits(std::vector<Matcher::result_t> &results, std::vector<Matcher::result_t> &tmp) {
    if(results.size() > 1){
        std::stable_sort(results.begin(), results.end(), Matcher::compareHitsByPosAndStrand);
        int prevDiagonal = INT_MAX;
        Matcher::result_t  currRegion;
        currRegion.dbKey = UINT_MAX;
        for (size_t resIdx = 0; resIdx < results.size(); resIdx++) {
            bool currQueryStrang = (results[resIdx].qStartPos > results[resIdx].qEndPos);
            int qStartPos = std::min(results[resIdx].qStartPos,  results[resIdx].qEndPos);
            int qEndPos = std::max(results[resIdx].qStartPos,  results[resIdx].qEndPos);
            bool currTargetStrand = (results[resIdx].dbStartPos > results[resIdx].dbEndPos);
            int dbStartPos = std::min(results[resIdx].dbStartPos, results[resIdx].dbEndPos);
            int dbEndPos = std::max(results[resIdx].dbStartPos, results[resIdx].dbEndPos);
            std::cout << results[resIdx].dbKey<< "\t" << qStartPos<< "\t" << qEndPos<< "\t" << dbStartPos<< "\t" << dbEndPos << "\t" << std::endl;
            if(currRegion.dbKey == UINT_MAX){
                currRegion = results[resIdx];
                currRegion.qStartPos = qStartPos;
                currRegion.qEndPos = qEndPos;
                currRegion.dbStartPos = dbStartPos;
                currRegion.dbEndPos = dbEndPos;
            }
            int currDiagonal = qStartPos - dbStartPos;
            int nextDiagonal = INT_MAX;
            bool nextQueryStrand = true;
            bool nextTargetStrand = true;
            const bool isDifferentKey = (currRegion.dbKey != results[resIdx].dbKey);
            const bool isLastElement  = (resIdx == results.size() - 1);
            if (isLastElement == false) {
                int nextqStartPos = std::min(results[resIdx+1].qStartPos,  results[resIdx+1].qEndPos);
                int nextdbStartPos = std::min(results[resIdx+1].dbStartPos, results[resIdx+1].dbEndPos);
                nextDiagonal = nextqStartPos - nextdbStartPos;
                nextQueryStrand =  (results[resIdx+1].qStartPos > results[resIdx+1].qEndPos);
                nextTargetStrand =  (results[resIdx+1].dbStartPos > results[resIdx+1].dbEndPos);
            }
            const bool queryIsOverlapping  = currRegion.qEndPos >= qStartPos   && currRegion.qEndPos <= qEndPos;
            const bool targetIsOverlapping = currRegion.dbEndPos >= dbStartPos && currRegion.dbEndPos <= dbEndPos;
            const bool sameNextDiagonal = (currDiagonal == nextDiagonal);
            const bool samePrevDiagonal = (currDiagonal == prevDiagonal);
            if ( (sameNextDiagonal || samePrevDiagonal ) && queryIsOverlapping && targetIsOverlapping) {
                currRegion.qStartPos = std::min(currRegion.qStartPos, qStartPos);
                currRegion.qEndPos = std::max(currRegion.qEndPos,  qEndPos);
                currRegion.dbStartPos  = std::min(currRegion.dbStartPos, dbStartPos);
                currRegion.dbEndPos  = std::max(currRegion.dbEndPos, dbEndPos);
            }

            prevDiagonal = currDiagonal;
            bool isDifferentNextDiagonal = (nextDiagonal != currDiagonal);
            bool isDifferentStrand = (nextQueryStrand != currQueryStrang ) || (nextTargetStrand != currTargetStrand );
            if(isDifferentKey || isLastElement || isDifferentNextDiagonal || isDifferentStrand){
                if(currQueryStrang){
                    std::swap(currRegion.qStartPos, currRegion.qEndPos);
                }
                if(currTargetStrand) {
                    std::swap(currRegion.dbStartPos, currRegion.dbEndPos);
                }
                tmp.push_back(currRegion);
                currRegion.dbKey = UINT_MAX;
            }
        }
    }
}

//

// We have serval options to consider
// Update query and target
//   Nucl/Nucl
//   TransNucl/TransNucl
// target update
//   Prot/Nucl
// query update
//   Nucl/Prot
void updateOffset(char* data, std::vector<Matcher::result_t> &results, const Orf::SequenceLocation *qloc,
                  IndexReader& tOrfDBr, bool targetNeedsUpdate, bool isNucleotideSearch, int thread_idx) {
    size_t startPos = results.size();
    Matcher::readAlignmentResults(results, data, true);
    size_t endPos = results.size();
    for (size_t i = startPos; i < endPos; i++) {
        Matcher::result_t &res = results[i];
        res.queryOrfStartPos = -1;
        res.queryOrfEndPos = -1;
        res.dbOrfStartPos = -1;
        res.dbOrfEndPos = -1;
        if (targetNeedsUpdate == true || qloc == NULL) {
            size_t targetId = tOrfDBr.sequenceReader->getId(res.dbKey);
            char *header = tOrfDBr.sequenceReader->getData(targetId, thread_idx);

            Orf::SequenceLocation tloc = Orf::parseOrfHeader(header);
            res.dbKey   = (tloc.id != UINT_MAX) ? tloc.id : res.dbKey;
            size_t from = (tloc.id != UINT_MAX) ? tloc.from : (tloc.strand == Orf::STRAND_MINUS) ? res.dbLen - 1 : 0;

            int dbStartPos = isNucleotideSearch ? res.dbStartPos : res.dbStartPos * 3;
            int dbEndPos   = isNucleotideSearch ? res.dbEndPos : res.dbEndPos * 3;
            res.dbOrfStartPos = from;
            res.dbOrfEndPos = tloc.to;
            if (tloc.strand == Orf::STRAND_MINUS) {
                res.dbStartPos = from - dbStartPos;
                res.dbEndPos   = from - dbEndPos;
                // account for last orf
                //  GGCACC
                //  GGCA
                //     ^
                //     last codon position
                //  GGCACC
                //    GGCA
                //    ^
                //     last codon position
                if(isNucleotideSearch == false){
                    res.dbEndPos = res.dbEndPos - 2;
                }
            } else {
                res.dbStartPos = from + dbStartPos;
                res.dbEndPos   = from + dbEndPos;
                if(isNucleotideSearch == false){
                    res.dbEndPos = res.dbEndPos + 2;
                }
            }
        }
        if (qloc != NULL) {
            int qStartPos = isNucleotideSearch ? res.qStartPos : res.qStartPos * 3;
            int qEndPos   = isNucleotideSearch ? res.qEndPos : res.qEndPos * 3;

            size_t from = (qloc->id != UINT_MAX) ? qloc->from : (qloc->strand == Orf::STRAND_MINUS) ? 0 : res.qLen - 1;
            res.queryOrfStartPos = from;
            res.queryOrfEndPos =  qloc->to;

            if (qloc->strand == Orf::STRAND_MINUS && qloc->id != UINT_MAX) {
                res.qStartPos  = from - qStartPos;
                res.qEndPos    = from - qEndPos;
                if(isNucleotideSearch == false){
                    res.qEndPos = res.qEndPos - 2;
                }
            } else {
                res.qStartPos  = from + qStartPos;
                res.qEndPos    = from + qEndPos;
                if(isNucleotideSearch == false){
                    res.qEndPos = res.qEndPos + 2;
                }
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
            res.dbLen = tSourceDbr->sequenceReader->getSeqLen(targetId);
        }
    }
}

int offsetalignment(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const bool touch = par.preloadMode != Parameters::PRELOAD_MODE_MMAP;
    int queryDbType = FileUtil::parseDbType(par.db1.c_str());
    if(Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_INDEX_DB)){
        DBReader<unsigned int> idxdbr(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        idxdbr.open(DBReader<unsigned int>::NOSORT);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&idxdbr);
        queryDbType=data.srcSeqType;
        idxdbr.close();
    }
    int targetDbType = FileUtil::parseDbType(par.db3.c_str());
    if(Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_INDEX_DB)){
        DBReader<unsigned int> idxdbr(par.db3.c_str(), par.db3Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        idxdbr.open(DBReader<unsigned int>::NOSORT);
        PrefilteringIndexData data = PrefilteringIndexReader::getMetadata(&idxdbr);
        targetDbType=data.srcSeqType;
        idxdbr.close();
    }

    IndexReader qOrfDbr(par.db2.c_str(), par.threads, IndexReader::HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    if (queryDbType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        return EXIT_FAILURE;
    }
    const bool queryNucl = Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES);
    IndexReader *qSourceDbr = NULL;
    if (queryNucl) {
        qSourceDbr = new IndexReader(par.db1.c_str(), par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX) : 0, DBReader<unsigned int>::USE_INDEX);
    }

    IndexReader * tOrfDbr;
    bool isSameOrfDB = (par.db2.compare(par.db4) == 0);
    if(isSameOrfDB){
        tOrfDbr = &qOrfDbr;
    }else{
        tOrfDbr = new IndexReader(par.db4.c_str(), par.threads, IndexReader::HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    if (targetDbType == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence/profile database.\n";
        return EXIT_FAILURE;
    }
    const bool targetNucl = Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES);
    IndexReader *tSourceDbr = NULL;
    bool isSameSrcDB = (par.db3.compare(par.db1) == 0);
    bool isNuclNuclSearch = false;
    bool isTransNucTransNucSearch = false;
    bool isTransNuclAln = false;
    if (targetNucl) {
        bool seqtargetNuc = true;
        if(isSameSrcDB){
            tSourceDbr = qSourceDbr;
        }else{
            tSourceDbr = new IndexReader(par.db3.c_str(), par.threads, IndexReader::SRC_SEQUENCES, (touch) ? IndexReader::PRELOAD_INDEX : 0, DBReader<unsigned int>::USE_INDEX );
        }

        if(Parameters::isEqualDbtype(tSourceDbr->getDbtype(), Parameters::DBTYPE_INDEX_DB)){
            IndexReader tseqDbr(par.db3, par.threads, IndexReader::SEQUENCES, 0, IndexReader::PRELOAD_INDEX);
            seqtargetNuc = Parameters::isEqualDbtype(tseqDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
            isTransNucTransNucSearch = Parameters::isEqualDbtype(tseqDbr.sequenceReader->getDbtype(), Parameters::DBTYPE_AMINO_ACIDS);
        }else{
            if(par.searchType == Parameters::SEARCH_TYPE_AUTO && (targetNucl == true && queryNucl == true )){
                Debug(Debug::WARNING) << "Assume that nucleotide search was performed\n";
                Debug(Debug::WARNING) << "If this is not correct than please provide the parameter --search-type 2 (translated) or 3 (nucleotide)\n";
            } else if(par.searchType == Parameters::SEARCH_TYPE_TRANSLATED){
                seqtargetNuc = false;
                isTransNucTransNucSearch = true;
            } else if(par.searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES){
                seqtargetNuc = true;
                isTransNucTransNucSearch = false;
            } else if(par.searchType == Parameters::SEARCH_TYPE_TRANS_NUCL_ALN){
                isTransNuclAln = true;
                seqtargetNuc = false;
                isTransNucTransNucSearch = true;
            }
        }

        isNuclNuclSearch = (queryNucl && targetNucl && seqtargetNuc);
    }

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
        Debug(Debug::INFO) << "Computing ORF lookup\n";
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
                size_t queryId = qOrfDbr.sequenceReader->getId(i);
                if (queryId == UINT_MAX) {
                    orfLookup[i] = UINT_MAX;
                    continue;
                }
                unsigned int queryKey = qOrfDbr.sequenceReader->getDbKey(queryId);
                char *header = qOrfDbr.sequenceReader->getData(queryId, thread_idx);
                Orf::SequenceLocation qloc = Orf::parseOrfHeader(header);
                unsigned int id = (qloc.id != UINT_MAX) ? qloc.id : queryKey;
                orfLookup[i] = id;
            }
        }
        Debug(Debug::INFO) << "Computing contig offsets\n";
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

        Debug(Debug::INFO) << "Computing contig lookup\n";
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

    size_t entryCount = alnDbr.getSize();
    if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
        entryCount = maxContigKey + 1;
    }
    Debug::Progress progress(entryCount);

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

        std::string newBacktrace;
        newBacktrace.reserve(300);

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < entryCount; ++i) {
            progress.updateProgress();
            unsigned int queryKey=UINT_MAX;
            unsigned int qLen = UINT_MAX;

            if (Parameters::isEqualDbtype(queryDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
                queryKey = i;
                if (contigExists[i] == 0) {
                    continue;
                }
                if (qSourceDbr != NULL) {
                    size_t queryId = qSourceDbr->sequenceReader->getId(queryKey);
                    qLen = qSourceDbr->sequenceReader->getSeqLen(queryId);
                }
                unsigned int *orfKeys = &contigLookup[contigOffsets[i]];
                size_t orfCount = contigOffsets[i + 1] - contigOffsets[i];
                for (unsigned int j = 0; j < orfCount; ++j) {
                    unsigned int orfKey = orfKeys[j];
                    size_t orfId = alnDbr.getId(orfKey);
                    // this is needed when alnDbr does not contain all identifier of the queryDB
                    if(orfId==UINT_MAX){
                        continue;
                    }
                    char *data = alnDbr.getData(orfId, thread_idx);
                    size_t queryId = qOrfDbr.sequenceReader->getId(orfKey);
                    char *header = qOrfDbr.sequenceReader->getData(queryId, thread_idx);
                    Orf::SequenceLocation qloc = Orf::parseOrfHeader(header);
                    if(qloc.id == UINT_MAX){
                        updateOffset(data, results, NULL, *tOrfDbr, (isNuclNuclSearch||isTransNucTransNucSearch), isNuclNuclSearch, thread_idx);
                    }else{
                        updateOffset(data, results, &qloc, *tOrfDbr, (isNuclNuclSearch||isTransNucTransNucSearch), isNuclNuclSearch, thread_idx);
                    }
                    // do not merge entries
                    if(par.mergeQuery == false){
                        for(size_t i = 0; i < results.size(); i++) {
                            Matcher::result_t &res = results[i];
                            bool hasBacktrace = (res.backtrace.size() > 0);
                            if (isTransNuclAln == true && isNuclNuclSearch == false && isTransNucTransNucSearch == true && hasBacktrace) {
                                newBacktrace.reserve(res.backtrace.length() * 3);
                                Matcher::result_t::protein2nucl(res.backtrace, newBacktrace);
                                res.backtrace = newBacktrace;
                                newBacktrace.clear();
                            }
                            size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false, true);
                            ss.append(buffer, len);
                        }
                        resultWriter.writeData(ss.c_str(), ss.length(), queryKey, thread_idx);
                        results.clear();
                        ss.clear();
                        tmp.clear();
                    }
                }
            } else if (Parameters::isEqualDbtype(targetDbType, Parameters::DBTYPE_NUCLEOTIDES)) {
                queryKey = alnDbr.getDbKey(i);
                if (qSourceDbr != NULL) {
                    size_t queryId = qSourceDbr->sequenceReader->getId(queryKey);
                    qLen = qSourceDbr->sequenceReader->getSeqLen(queryId);
                }
                char *data = alnDbr.getData(i, thread_idx);
                updateOffset(data, results, NULL, *tOrfDbr, true, isNuclNuclSearch, thread_idx);
            }
            if(par.mergeQuery == true){
                updateLengths(results, qLen, tSourceDbr);
                if(par.chainAlignment == false){
                    std::stable_sort(results.begin(), results.end(), Matcher::compareHits);
                    for(size_t i = 0; i < results.size(); i++){
                        Matcher::result_t &res = results[i];
                        bool hasBacktrace = (res.backtrace.size() > 0);
                        if (isTransNuclAln == true && isNuclNuclSearch == false && isTransNucTransNucSearch == true && hasBacktrace) {
                            newBacktrace.reserve(res.backtrace.length() * 3);
                            Matcher::result_t::protein2nucl(res.backtrace, newBacktrace);
                            res.backtrace = newBacktrace;
                            newBacktrace.clear();
                        }
                        size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false, true);
                        ss.append(buffer, len);
                    }
                    resultWriter.writeData(ss.c_str(), ss.length(), queryKey, thread_idx);
                } else if(par.chainAlignment == true){
                    chainAlignmentHits(results, tmp);
                    for(size_t i = 0; i < tmp.size(); i++){
                        Matcher::result_t &res = tmp[i];
                        bool hasBacktrace = (res.backtrace.size() > 0);
                        if (isTransNuclAln == true && isNuclNuclSearch == false && isTransNucTransNucSearch == true && hasBacktrace) {
                            newBacktrace.reserve(res.backtrace.length() * 3);
                            Matcher::result_t::protein2nucl(res.backtrace, newBacktrace);
                            res.backtrace = newBacktrace;
                            newBacktrace.clear();
                        }
                        size_t len = Matcher::resultToBuffer(buffer, res, hasBacktrace, false, true);
                        ss.append(buffer, len);
                    }
                    resultWriter.writeData(ss.c_str(), ss.length(), queryKey, thread_idx);
                }
                ss.clear();
                results.clear();
                tmp.clear();
            }
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

    if(isSameOrfDB == false){
        delete tOrfDbr;
    }
    if(tSourceDbr != NULL){
        if(isSameSrcDB==false){
            delete tSourceDbr;
        }
    }

    if(qSourceDbr != NULL){
        delete qSourceDbr;
    }
    alnDbr.close();

    return EXIT_SUCCESS;
}
