#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif


static bool compareToFirstInt(const std::pair<unsigned int, unsigned int> &lhs, const std::pair<unsigned int, unsigned int> &rhs) {
    return (lhs.first <= rhs.first);
}

int addtaxonomy(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::vector<std::pair<unsigned int, unsigned int>> mapping;
    if (FileUtil::fileExists((par.db1 + "_mapping").c_str()) == false) {
        Debug(Debug::ERROR) << par.db1 << "_mapping does not exist. Run createtaxdb to create taxonomy mapping.\n";
        EXIT(EXIT_FAILURE);
    }
    const bool isSorted = Util::readMapping(par.db1 + "_mapping", mapping);
    if (isSorted == false) {
        std::stable_sort(mapping.begin(), mapping.end(), compareToFirstInt);
    }
    if (mapping.size() == 0) {
        Debug(Debug::ERROR) << par.db1 << "_mapping is empty. Rerun createtaxdb to recreate taxonomy mapping.\n";
        EXIT(EXIT_FAILURE);
    }
    NcbiTaxonomy *t = NcbiTaxonomy::openTaxonomy(par.db1);
    std::vector<std::string> ranks = NcbiTaxonomy::parseRanks(par.lcaRanks);

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    size_t taxonNotFound = 0;
    size_t deletedNodes = 0;
    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        const char *entry[255];
        std::string result;
        result.reserve(4096);

#pragma omp for schedule(dynamic, 10) reduction (+: deletedNodes, taxonNotFound)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getEntryLen(i);

            if (length == 1) {
                continue;
            }
            std::pair<unsigned int, unsigned int> val;
            std::vector<std::pair<unsigned int, unsigned int> >::iterator mappingIt;
            if (par.pickIdFrom == Parameters::EXTRACT_QUERY) {
                val.first = key;
                mappingIt = std::upper_bound(mapping.begin(), mapping.end(), val, compareToFirstInt);
                if (mappingIt == mapping.end() || mappingIt->first != val.first) {
                    taxonNotFound++;
                    continue;
                }
            }

            while (*data != '\0') {
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns == 0) {
                    Debug(Debug::WARNING) << "Empty entry: " << i << "\n";
                    data = Util::skipLine(data);
                    continue;
                }
                if (par.pickIdFrom == Parameters::EXTRACT_TARGET) {
                    unsigned int id = Util::fast_atoi<unsigned int>(entry[0]);
                    val.first = id;
                    mappingIt = std::upper_bound(mapping.begin(), mapping.end(), val, compareToFirstInt);
                    if (mappingIt == mapping.end() || mappingIt->first != val.first) {
                        taxonNotFound++;
                        data = Util::skipLine(data);
                        continue;
                    }
                }
                unsigned int taxon = mappingIt->second;
                TaxonNode const *node = t->taxonNode(taxon, false);
                if (node == NULL) {
                    deletedNodes++;
                    data = Util::skipLine(data);
                    continue;
                }
                char *nextData = Util::skipLine(data);
                size_t dataSize = nextData - data;
                result.append(data, dataSize - 1);
                result += '\t' + SSTR(node->taxId) + '\t' + node->rank + '\t' + node->name;
                if (!ranks.empty()) {
                    std::string lcaRanks = Util::implode(t->AtRanks(node, ranks), ';');
                    result += '\t' + lcaRanks;
                }
                if (par.showTaxLineage == 1) {
                    result += '\t' + t->taxLineage(node, true);
                }
                if (par.showTaxLineage == 2) {
                    result += '\t' + t->taxLineage(node, false);
                }
                result += '\n';
                data = Util::skipLine(data);
            }
            writer.writeData(result.c_str(), result.size(), key, thread_idx);
            result.clear();
        }
    }
    Debug(Debug::INFO) << "Taxonomy for " << taxonNotFound << " entries not found and " << deletedNodes << " are deleted\n";
    writer.close();
    reader.close();
    delete t;
    return EXIT_SUCCESS;
}

