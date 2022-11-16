#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "MappingReader.h"

#ifdef OPENMP
#include <omp.h>
#endif

int addtaxonomy(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    NcbiTaxonomy *t = NcbiTaxonomy::openTaxonomy(par.db1);
    MappingReader mapping(par.db1);
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
            unsigned int taxon = 0;
            if (par.pickIdFrom == Parameters::EXTRACT_QUERY) {
                taxon = mapping.lookup(key);
                if (taxon == 0) {
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
                    taxon = mapping.lookup(id);
                    if (taxon == 0) {
                        taxonNotFound++;
                        data = Util::skipLine(data);
                        continue;
                    }
                }
                TaxonNode const *node = t->taxonNode(taxon, false);
                if (node == NULL) {
                    deletedNodes++;
                    data = Util::skipLine(data);
                    continue;
                }
                char *nextData = Util::skipLine(data);
                size_t dataSize = nextData - data;
                result.append(data, dataSize - 1);
                result.append(1, '\t');
                result.append(SSTR(node->taxId));
                result.append(1, '\t');
                result.append(t->getString(node->rankIdx));
                result.append(1, '\t');
                result.append(t->getString(node->nameIdx));
                if (!ranks.empty()) {
                    result.append(1, '\t');
                    result.append(Util::implode(t->AtRanks(node, ranks), ';'));
                }
                if (par.showTaxLineage == 1) {
                    result.append(1, '\t');
                    result.append(t->taxLineage(node, true));
                }
                if (par.showTaxLineage == 2) {
                    result.append(1, '\t');
                    result.append(t->taxLineage(node, false));
                }
                result.append(1, '\n');
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

