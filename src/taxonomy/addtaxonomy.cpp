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


int addtaxonomy(int argc, const char **argv, const Command& command) {
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
    std::vector<std::string> ranks = Util::split(par.lcaRanks, ":");

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    Debug(Debug::INFO) << "Loading NCBI taxonomy\n";
    NcbiTaxonomy t(namesFile, nodesFile, mergedFile, delnodesFile);

    Debug(Debug::INFO) << "Add taxonomy information \n";
    size_t taxonNotFound=0;
    Debug::Progress progress(reader.getSize());

    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        const char *entry[255];
        char buffer[10000];
        std::string resultData;
        resultData.reserve(4096);

        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getSeqLens(i);

            if (length == 1) {
                continue;
            }

            std::vector<int> taxa;
            while (*data != '\0') {
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns == 0) {
                    Debug(Debug::WARNING) << "Empty entry: " << i << "!";
                    data = Util::skipLine(data);
                    continue;
                }
                unsigned int id = Util::fast_atoi<unsigned int>(entry[0]);
                std::pair<unsigned int, unsigned int> val;
                val.first = id;
                std::vector<std::pair<unsigned int, unsigned int>>::iterator mappingIt = std::upper_bound(
                        mapping.begin(), mapping.end(), val,  ffindexFilter::compareToFirstInt);

                if (mappingIt->first != val.first) {
                     __sync_fetch_and_add(&taxonNotFound, 1);
//                    Debug(Debug::WARNING) << "No taxon mapping provided for id " << id << "\n";
                    data = Util::skipLine(data);
                    continue;
                }
                unsigned int taxon = mappingIt->second;
                TaxonNode* node = t.findNode(taxon);
                if(node == NULL){
                    Debug(Debug::WARNING) << "Deleted node " << taxon << "!\n";
                    data = Util::skipLine(data);
                    continue;
                }
                char * nextData = Util::skipLine(data);
                size_t dataSize = nextData - data;
                resultData.append(data, dataSize-1);
                resultData.push_back('\t');
                std::string lcaRanks = Util::implode(t.AtRanks(node, ranks), ':');
                int len;
                if (ranks.empty() == false) {
                    len = snprintf(buffer, 10000, "%d\t%s\t%s\n",
                                   node->taxon, node->rank.c_str(), node->name.c_str());
                } else {
                    len = snprintf(buffer, 10000, "%d\t%s\t%s\t%s\n",
                         node->taxon, node->rank.c_str(), node->name.c_str(), lcaRanks.c_str());
                }
                if(len < 0){
                    Debug(Debug::WARNING) << "Taxon record could not be written. Entry: " << i << "\t" << columns << "!\n";
                    data = Util::skipLine(data);
                    continue;
                }
                resultData.append(buffer, len),
                data = Util::skipLine(data);
            }
            writer.writeData(resultData.c_str(), resultData.size(), key, thread_idx);
            resultData.clear();
        }
    }
    Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Taxonomy for " << taxonNotFound << " entries  not found.\n";

    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}
