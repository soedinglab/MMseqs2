#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "filterdb.h"
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

int lca(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

    std::string nodesFile = par.db1 + "_nodes.dmp";
    std::string namesFile = par.db1 + "_names.dmp";
    std::string mergedFile = par.db1 + "_merged.dmp";
    std::string delnodesFile = par.db1 + "_delnodes.dmp";
    if (FileUtil::fileExists(nodesFile.c_str())
        && FileUtil::fileExists(namesFile.c_str())
        && FileUtil::fileExists(mergedFile.c_str())
        && FileUtil::fileExists(delnodesFile.c_str())) {
    } else if (FileUtil::fileExists("nodes.dmp")
               && FileUtil::fileExists("names.dmp")
               && FileUtil::fileExists("merged.dmp")
               && FileUtil::fileExists("delnodes.dmp")) {
        nodesFile = "nodes.dmp";
        namesFile = "names.dmp";
        mergedFile = "merged.dmp";
        delnodesFile = "delnodes.dmp";
    } else {
        Debug(Debug::ERROR) << "names.dmp, nodes.dmp, merged.dmp or delnodes.dmp from NCBI taxdump could not be found!\n";
        EXIT(EXIT_FAILURE);
    }
    std::vector<std::pair<unsigned int, unsigned int>> mapping;
    if(FileUtil::fileExists(std::string(par.db1 + "_mapping").c_str()) == false){
        Debug(Debug::ERROR) << par.db1 + "_mapping" << " does not exist. Please create the taxonomy mapping!\n";
        EXIT(EXIT_FAILURE);
    }
    bool isSorted = Util::readMapping( par.db1 + "_mapping", mapping);
    if(isSorted == false){
        std::stable_sort(mapping.begin(), mapping.end(), ffindexFilter::compareFirstInt());
    }

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_TAXONOMICAL_RESULT);
    writer.open();

    std::vector<std::string> ranks = Util::split(par.lcaRanks, ":");

    // a few NCBI taxa are blacklisted by default, they contain unclassified sequences (e.g. metagenomes) or other sequences (e.g. plasmids)
    // if we do not remove those, a lot of sequences would be classified as Root, even though they have a sensible LCA
    std::vector<std::string> blacklist = Util::split(par.blacklist, ",");
    const size_t taxaBlacklistSize = blacklist.size();
    int* taxaBlacklist = new int[taxaBlacklistSize];
    for (size_t i = 0; i < taxaBlacklistSize; ++i) {
        taxaBlacklist[i] = Util::fast_atoi<int>(blacklist[i].c_str());
    }
    Debug::Progress progress(reader.getSize());
    Debug(Debug::INFO) << "Loading NCBI taxonomy\n";
    NcbiTaxonomy t(namesFile, nodesFile, mergedFile, delnodesFile);
    size_t taxonNotFound=0;
    Debug(Debug::INFO) << "Computing LCA\n";
    #pragma omp parallel
    {
        const char *entry[255];
        char buffer[1024];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getSeqLens(i);

            std::vector<int> taxa;
            while (*data != '\0') {
                int taxon;
                unsigned int id;
                std::pair<unsigned int, unsigned int> val;
                std::vector<std::pair<unsigned int, unsigned int>>::iterator mappingIt;
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns == 0) {
                    Debug(Debug::WARNING) << "Empty entry: " << i << "!";
                    goto next;
                }

                id = Util::fast_atoi<unsigned int>(entry[0]);

                val.first = id;
                mappingIt = std::upper_bound(
                        mapping.begin(), mapping.end(), val,  ffindexFilter::compareToFirstInt);

                if (mappingIt->first != val.first) {
                    __sync_fetch_and_add(&taxonNotFound, 1);
                    data = Util::skipLine(data);
                    continue;
                }
                taxon = mappingIt->second;

                // remove blacklisted taxa
                for (size_t j = 0; j < taxaBlacklistSize; ++j) {
                    if (t.IsAncestor(taxaBlacklist[j], taxon)) {
                        goto next;
                    }
                }

                taxa.emplace_back(taxon);

                next:
                data = Util::skipLine(data);
            }

            if(length == 1){
                snprintf(buffer, 1024, "0\tno rank\tunclassified\n");
                writer.writeData(buffer, strlen(buffer), key, thread_idx);
                continue;
            }

            TaxonNode* node = t.LCA(taxa);
            if (node == NULL) {
                snprintf(buffer, 1024, "0\tno rank\tunclassified\n");
                writer.writeData(buffer, strlen(buffer), key, thread_idx);
                continue;
            }

            if (ranks.empty() == false) {
                std::string lcaRanks = Util::implode(t.AtRanks(node, ranks), ':');
                snprintf(buffer, 1024, "%d\t%s\t%s\t%s\n",
                         node->taxon, node->rank.c_str(), node->name.c_str(), lcaRanks.c_str());
                writer.writeData(buffer, strlen(buffer), key, thread_idx);
            } else {
                snprintf(buffer, 1024, "%d\t%s\t%s\n",
                         node->taxon, node->rank.c_str(), node->name.c_str());
                writer.writeData(buffer, strlen(buffer), key, thread_idx);
            }
        }
    };
    Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Taxonomy for " << taxonNotFound << " entries not found.\n";

    writer.close();
    reader.close();

    delete[] taxaBlacklist;

    return EXIT_SUCCESS;
}
