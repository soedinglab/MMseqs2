#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "MemoryMapped.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"

#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>

#include <map>

#ifdef OPENMP
#include <omp.h>
#endif

int pairaln(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);

    std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
    MappingReader* mapping = new MappingReader(db2NoIndexName);
    const int dbaccessMode = (DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader *tDbr;
    if (sameDB) {
        tDbr = &qDbr;
    } else {
        tDbr = new IndexReader(par.db2, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::NOSORT);
    if(alnDbr.getSize() % 2 != 0){
        Debug(Debug::ERROR) << "Alignment database does not seem to be paired\n";
        EXIT(EXIT_FAILURE);
    }
    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), localThreads,  par.compressed, alnDbr.getDbtype());
    resultWriter.open();

    Debug::Progress progress(alnDbr.getSize());
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[1024 + 32768*4];
        std::string result;
        result.reserve(1024 * 1024);
        std::vector<Matcher::result_t> resultA;
        resultA.reserve(100000);
        std::vector<Matcher::result_t> resultB;
        resultB.reserve(100000);
        std::unordered_map<unsigned int, Matcher::result_t *> findPair;
        std::string outputA;
        outputA.reserve(100000);
        std::string outputB;
        outputB.reserve(100000);
#pragma omp  for schedule(static, 2)
        for (size_t i = 0; i < alnDbr.getSize(); i+=2) {
            progress.updateProgress();
            progress.updateProgress();
            resultA.clear();
            resultB.clear();
            outputA.clear();
            outputB.clear();
            Matcher::readAlignmentResults(resultA, alnDbr.getData(i, thread_idx), true);
            Matcher::readAlignmentResults(resultB, alnDbr.getData(i+1, thread_idx), true);

            for (size_t resIdx = 0; resIdx < resultA.size(); ++resIdx) {
                unsigned int taxon = mapping->lookup(resultA[resIdx].dbKey);
                if(findPair.find(taxon) == findPair.end()){
                    findPair.insert({taxon, &resultA[resIdx]});
                }
            }
            // find pairs
            for (size_t resIdx = 0; resIdx < resultB.size(); ++resIdx) {
                unsigned int taxon = mapping->lookup(resultB[resIdx].dbKey);
                // found pair!
                if(findPair.find(taxon) != findPair.end()) {
                    bool hasBacktrace = (resultB[resIdx].backtrace.size() > 0);
                    size_t len = Matcher::resultToBuffer(buffer, *findPair[taxon], hasBacktrace, false, false);
                    outputA.append(buffer, len);
                    len = Matcher::resultToBuffer(buffer, resultB[resIdx], hasBacktrace, false, false);
                    outputB.append(buffer, len);
                    findPair.erase (taxon);
                }
            }
            resultWriter.writeData(outputA.c_str(), outputA.length(), alnDbr.getDbKey(i), thread_idx);
            resultWriter.writeData(outputB.c_str(), outputB.length(), alnDbr.getDbKey(i+1), thread_idx);
            findPair.clear();
        }
    }

    resultWriter.close();
    // clean up
    delete mapping;
    alnDbr.close();
    if (sameDB == false) {
        delete tDbr;
    }

    return 0;
}

