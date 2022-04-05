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

    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_LOOKUP_REV);
    qdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int>::LookupEntry* lookup = qdbr.getLookup();
    unsigned int maxFileNumber = 0;
    for (unsigned int i = 0; i < qdbr.getLookupSize(); i++) {
        maxFileNumber = std::max(maxFileNumber, lookup[i].fileNumber);
    }
    //build a mapping from file number to ids from lookup
    std::vector<std::vector<unsigned int>> fileToIds(maxFileNumber + 1, std::vector<unsigned int>());
    for (size_t i = 0; i < qdbr.getLookupSize(); ++i) {
        fileToIds[lookup[i].fileNumber].push_back(lookup[i].id);
    }

    std::string db2NoIndexName = PrefilteringIndexReader::dbPathWithoutIndex(par.db2);
    MappingReader* mapping = new MappingReader(db2NoIndexName);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::NOSORT);

    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, alnDbr.getDbtype());
    resultWriter.open();

    Debug::Progress progress(fileToIds.size());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        std::vector<Matcher::result_t> result;
        result.reserve(100000);
        std::unordered_map<unsigned int, size_t> findPair;
        std::string output;
        output.reserve(100000);
#pragma omp for schedule(dynamic, 1)
	for (size_t fileNumber = 0; fileNumber < fileToIds.size(); fileNumber++) {
            char buffer[1024 + 32768 * 4];
            findPair.clear();
            progress.updateProgress();

            // find intersection between all proteins
            for (size_t i = 0; i < fileToIds[fileNumber].size(); i++) {
                result.clear();
                size_t id = fileToIds[fileNumber][i];
                Matcher::readAlignmentResults(result, alnDbr.getData(id, thread_idx), true);
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
                    if (taxon == prevTaxon) {
                        continue;
                    }
                    if (findPair.find(taxon) != findPair.end()) {
                        findPair[taxon]++;
                    } else {
                        findPair.emplace(taxon, 1);
                    }
                    prevTaxon = taxon;
                }
            }

            for (size_t i = 0; i < fileToIds[fileNumber].size(); i++) {
                result.clear();
                output.clear();
                size_t id = fileToIds[fileNumber][i];
                Matcher::readAlignmentResults(result, alnDbr.getData(id, thread_idx), true);
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
                        && findPair[taxon] == fileToIds[fileNumber].size()) {
                        bool hasBacktrace = result[resIdx].backtrace.size() > 0;
                        size_t len = Matcher::resultToBuffer(buffer, result[resIdx], hasBacktrace, false, false);
                        output.append(buffer, len);
                    }
                    prevTaxon = taxon;
                }
                resultWriter.writeData(output.c_str(), output.length(), alnDbr.getDbKey(id), thread_idx);
            }
        }
    }
    resultWriter.close();
    qdbr.close();
    // clean up
    delete mapping;
    alnDbr.close();
    return EXIT_SUCCESS;
}

