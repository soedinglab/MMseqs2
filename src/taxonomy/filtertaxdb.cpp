#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int filtertaxdb(int argc, const char **argv, const Command& command) {
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


    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    std::vector<std::string> ranks = Util::split(par.lcaRanks, ":");

    // a few NCBI taxa are blacklisted by default, they contain unclassified sequences (e.g. metagenomes) or other sequences (e.g. plasmids)
    // if we do not remove those, a lot of sequences would be classified as Root, even though they have a sensible LCA
    std::vector<std::string> list = Util::split(par.taxonList, ",");
    const size_t taxListSize = list.size();
    int* taxalist = new int[taxListSize];
    for (size_t i = 0; i < taxListSize; ++i) {
        taxalist[i] = Util::fast_atoi<int>(list[i].c_str());
    }


    bool invertSelection = par.invertSelection;
    Debug(Debug::INFO) << "Loading NCBI taxonomy\n";
    NcbiTaxonomy t(namesFile, nodesFile, mergedFile, delnodesFile);
    Debug::Progress progress(reader.getSize());

    Debug(Debug::INFO) << "Computing LCA\n";
    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        const char *entry[255];

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
                int taxon;
                bool isAncestor = false;
                size_t j;
                std::vector<std::pair<unsigned int, unsigned int>>::iterator mappingIt;
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns == 0) {
                    Debug(Debug::WARNING) << "Empty entry: " << i << "!";
                    goto next;
                }

                taxon = Util::fast_atoi<unsigned int>(entry[0]);
                writer.writeStart(thread_idx);


                // remove blacklisted taxa
                for (j = 0; j < taxListSize && isAncestor == false; ++j) {
                    if(taxalist[j] == 0){
                        isAncestor = std::max(isAncestor, (taxon == 0) ? true : false);
                        continue;
                    }
                    if(taxon == 0){
                        continue;
                    }
                    isAncestor = std::max(isAncestor, t.IsAncestor(taxalist[j], taxon));
                }

                if(invertSelection == true) {
                    if (isAncestor == true) {
                        goto next;
                    } else if (isAncestor == false) {
                        char * nextData = Util::skipLine(data);
                        size_t dataSize = nextData - data;
                        writer.writeAdd(data, dataSize, thread_idx);
                    }
                } else {
                    if (isAncestor == true) {
                        char * nextData = Util::skipLine(data);
                        size_t dataSize = nextData - data;
                        writer.writeAdd(data, dataSize, thread_idx);
                    } else if (isAncestor == false) {
                        goto next;
                    }
                }
                next:
                data = Util::skipLine(data);
            }
            writer.writeEnd(key, thread_idx);
        }
    };

    Debug(Debug::INFO) << "\n";

    writer.close();
    reader.close();

    delete[] taxalist;

    return EXIT_SUCCESS;
}
