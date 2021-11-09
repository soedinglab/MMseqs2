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

#include <map>

#ifdef OPENMP
#include <omp.h>
#endif


// need for sorting the results
static bool compareByTaxId(const Matcher::result_t &first, const Matcher::result_t &second) {
    return (first.dbOrfStartPos < second.dbOrfStartPos);
}

int pairaln(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
    MappingReader* mapping = new MappingReader(db2NoIndexName);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::NOSORT);

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, par.compressed, alnDbr.getDbtype());
    resultWriter.open();

    Debug::Progress progress(alnDbr.getSize());
    do {
        int thread_idx = 0;
        char buffer[1024 + 32768*4];
        std::vector<Matcher::result_t> result;
        result.reserve(100000);
        std::unordered_map<unsigned int, size_t> findPair;
        std::string output;
        output.reserve(100000);

        // find intersection between all proteins
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            result.clear();
            Matcher::readAlignmentResults(result, alnDbr.getData(i, thread_idx), true);
            for (size_t resIdx = 0; resIdx < result.size(); ++resIdx) {
                unsigned int taxon = mapping->lookup(result[resIdx].dbKey);
                // we don't want to introduce a new field, reuse existing unused field here
                result[resIdx].dbOrfStartPos = taxon;
            }
            std::stable_sort(result.begin(), result.end(), compareByTaxId);
            unsigned int prevTaxon = UINT_MAX;
            // find pairs
            for (size_t resIdx = 0; resIdx < result.size(); ++resIdx) {
                unsigned int taxon = mapping->lookup(result[resIdx].dbKey);
                if (taxon == prevTaxon){
                    continue;
                }
                if (findPair.find(taxon) != findPair.end()) {
                    findPair[taxon]++;
                }else{
                    findPair.emplace(taxon, 1);
                }
                prevTaxon = taxon;
            }
        }

        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            progress.updateProgress();
            result.clear();
            output.clear();
            Matcher::readAlignmentResults(result, alnDbr.getData(i, thread_idx), true);
            // find pairs
            for (size_t resIdx = 0; resIdx < result.size(); ++resIdx) {
                unsigned int taxon = mapping->lookup(result[resIdx].dbKey);
                // we don't want to introduce a new field, reuse existing unused field here
                result[resIdx].dbOrfStartPos = taxon;
            }
            // stable sort is required to assure that best hit is first per taxon
            std::stable_sort(result.begin(), result.end(), compareByTaxId);
            unsigned int prevTaxon = UINT_MAX;
            for (size_t resIdx = 0; resIdx < result.size(); ++resIdx) {
                unsigned int taxon = result[resIdx].dbOrfStartPos;
                // found pair!
                if (taxon != prevTaxon
                    && findPair.find(taxon) != findPair.end()
                    && findPair[taxon] == alnDbr.getSize()) {
                    bool hasBacktrace = result[resIdx].backtrace.size() > 0;
                    size_t len = Matcher::resultToBuffer(buffer, result[resIdx], hasBacktrace, false, false);
                    output.append(buffer, len);
                }
                prevTaxon = taxon;
            }
            resultWriter.writeData(output.c_str(), output.length(), alnDbr.getDbKey(i), thread_idx);
        }
    } while(false);

    resultWriter.close();
    // clean up
    delete mapping;
    alnDbr.close();
    return EXIT_SUCCESS;
}

