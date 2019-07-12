#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include <algorithm>
#include <unordered_map>

#ifdef OPENMP
#include <omp.h>
#endif

static bool compareToFirstInt(const std::pair<unsigned int, unsigned int>& lhs, const std::pair<unsigned int, unsigned int>&  rhs){
    return (lhs.first <= rhs.first);
}

template<typename K, typename V>
V at(const std::unordered_map<K, V>& map, K key, V default_value = V()) {
    typename std::unordered_map<K, V>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return default_value;
    } else {
        return it->second;
    }
}

unsigned int cladeCountVal(const std::unordered_map<TaxID, TaxonCounts>& map, TaxID key) {
    typename std::unordered_map<TaxID, TaxonCounts>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return 0;
    } else {
        return it->second.cladeCount;
    }
}


void taxReport(FILE* FP,
        const NcbiTaxonomy& taxDB,
        const std::unordered_map<TaxID, TaxonCounts> & cladeCounts,
        unsigned long totalReads,
        TaxID taxID = 0, int depth = 0) {

    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end()? 0 : it->second.cladeCount;
    unsigned int taxCount = it == cladeCounts.end()? 0 : it->second.taxCount;
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "%.4f\t%i\t%i\tno rank\t0\tunclassified\n",
                    100 * cladeCount / double(totalReads),
                    cladeCount, taxCount);
        }
        taxReport(FP, taxDB, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode* taxon = taxDB.taxonNode(taxID);
        fprintf(FP, "%.4f\t%i\t%i\t%s\t%i\t%s%s\n",
                100*cladeCount/double(totalReads), cladeCount, taxCount,
                taxon->rank.c_str(), taxID, std::string(2*depth, ' ').c_str(), taxon->name.c_str());

        std::vector<TaxID> children = it->second.children;
        std::sort(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts,b); });
        for (TaxID childTaxId : children) {
            if (cladeCounts.count(childTaxId)) {
                taxReport(FP, taxDB, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}

int taxonomyreport(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

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
        std::stable_sort(mapping.begin(), mapping.end(), compareToFirstInt);
    }

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // TODO: Better way to get file specified by param3?
    FILE *resultFP = fopen(par.db3.c_str(), "w");

    // 1. Read taxonomy
    Debug(Debug::INFO) << "Loading NCBI taxonomy\n";
    NcbiTaxonomy taxDB(namesFile, nodesFile, mergedFile);

    // 2. Read LCA file
    Debug::Progress progress(reader.getSize());
    Debug(Debug::INFO) << "Reading LCA results\n";


    std::unordered_map<TaxID, unsigned int> taxCounts;

//  Currentlly not parallel
//    #pragma omp parallel
    {
        const char *entry[255];
        //char buffer[1024];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

//        #pragma omp for schedule(dynamic, 10) reduction (+:taxonNotFound, found)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            char *data = reader.getData(i, thread_idx);

            const size_t columns = Util::getWordsOfLine(data, entry, 255);
            if (columns == 0) {
                Debug(Debug::WARNING) << "Empty entry: " << i << "!";
            } else {
                int taxon = Util::fast_atoi<int>(entry[0]);
                ++taxCounts[taxon];
                //__sync_fetch_and_add(&(offsets[kmerIdx]), 1);
            }
        }
    };
    Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Found " << taxCounts.size() << " different taxa for " << reader.getSize() << " different reads.\n";
    Debug(Debug::INFO) << taxCounts.at(0) << " reads are unclassified.\n";

    std::unordered_map<TaxID, TaxonCounts> cladeCounts = taxDB.getCladeCounts(taxCounts);
    taxReport(resultFP, taxDB, cladeCounts, reader.getSize());

    reader.close();
    return EXIT_SUCCESS;
}
