#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "TaxonomyExpression.h"
#include <map>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

static bool compareToFirstInt(const std::pair<unsigned int, unsigned int>& lhs, const std::pair<unsigned int, unsigned int>& rhs){
    return (lhs.first <= rhs.first);
}

int filtertaxseqdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    
    // open mapping (dbKey to taxid)
    std::vector<std::pair<unsigned int, unsigned int>> mapping;
    if (FileUtil::fileExists(std::string(par.db1 + "_mapping").c_str()) == false) {
        Debug(Debug::ERROR) << par.db1 + "_mapping" << " does not exist. Please create the taxonomy mapping!\n";
        EXIT(EXIT_FAILURE);
    }
    bool isSorted = Util::readMapping(par.db1 + "_mapping", mapping);
    if (isSorted == false) {
        std::stable_sort(mapping.begin(), mapping.end(), compareToFirstInt);
    }

    // open taxonomy - evolutionary relationships amongst taxa
    NcbiTaxonomy * t = NcbiTaxonomy::openTaxonomy(par.db1);
    
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    const bool isCompressed = reader.isCompressed();

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    writer.open();

    // a few NCBI taxa are blacklisted by default, they contain unclassified sequences (e.g. metagenomes) or other sequences (e.g. plasmids)
    // if we do not remove those, a lot of sequences would be classified as Root, even though they have a sensible LCA

    Debug::Progress progress(reader.getSize());

    Debug(Debug::INFO) << "Computing LCA\n";
    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        TaxonomyExpression taxonomyExpression(par.taxonList, *t);
        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            size_t offset = reader.getOffset(i);
            size_t length = reader.getEntryLen(i);
            unsigned int taxon = 0;

            // match dbKey to its taxon based on mapping
            std::pair<unsigned int, unsigned int> val;
            val.first = key;
            std::vector<std::pair<unsigned int, unsigned int>>::iterator mappingIt;
            mappingIt = std::upper_bound(mapping.begin(), mapping.end(), val, compareToFirstInt);
            if (mappingIt == mapping.end() || mappingIt->first != val.first) {
                taxon = 0;
            } else {
                taxon = mappingIt->second;
            }

            // if taxon is a descendent of the requested taxid, it will be retained.
            // e.g. if in taxonomyExpression taxid=2 (bacteria) and taxon=562 (E.coli) 
            // then the check will return "true" cause taxon is descendent of taxid
            if (taxonomyExpression.isAncestor(taxon)) {
                if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
                    writer.writeIndexEntry(key, offset, length, thread_idx);
                } else {
                    char* data = reader.getDataUncompressed(i);
                    size_t originalLength = reader.getEntryLen(i);
                    size_t entryLength = std::max(originalLength, static_cast<size_t>(1)) - 1;

                    if (isCompressed) {
                        // copy also the null byte since it contains the information if compressed or not
                        entryLength = *(reinterpret_cast<unsigned int *>(data)) + sizeof(unsigned int) + 1;
                        writer.writeData(data, entryLength, key, thread_idx, false, false);
                    } else {
                        writer.writeData(data, entryLength, key, thread_idx, true, false);
                    }
                    writer.writeIndexEntry(key, writer.getStart(thread_idx), originalLength, thread_idx);
                }
            }
        }
    };
    Debug(Debug::INFO) << "\n";

    writer.close(true);
    if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
        DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_NO_DATA_INDEX);
    } else {
        DBWriter::writeDbtypeFile(par.db2.c_str(), reader.getDbtype(), isCompressed);
        DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_ANCILLARY);
    }

    reader.close();
    delete t;

    return EXIT_SUCCESS;
}
