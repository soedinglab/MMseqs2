#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "krona_prelude.html.h"

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


void taxReport(FILE* FP, const NcbiTaxonomy& taxDB, const std::unordered_map<TaxID, TaxonCounts> & cladeCounts,unsigned long totalReads,TaxID taxID = 0, int depth = 0) {
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

std::string escapeAttribute(const std::string& data) {
    std::string buffer;
    buffer.reserve(data.size() * 1.1);
    for (size_t i = 0; i < data.size(); ++i) {
        switch (data[i]) {
            case '&':
                buffer.append("&amp;");
                break;
            case '\"':
                buffer.append("&quot;");
                break;
            case '\'':
                buffer.append("&apos;");
                break;
            case '<':
                buffer.append("&lt;");
                break;
            case '>':
                buffer.append("&gt;");
                break;
            default:
                buffer.append(1, data[i]);
                break;
        }
    }
    return buffer;
}

void kronaReport(FILE* FP, const NcbiTaxonomy& taxDB, const std::unordered_map<TaxID, TaxonCounts> & cladeCounts,unsigned long totalReads,TaxID taxID = 0, int depth = 0) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end()? 0 : it->second.cladeCount;
//    unsigned int taxCount = it == cladeCounts.end()? 0 : it->second.taxCount;
    if (cladeCount == 0) {
        return;
    }
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "<node name=\"unclassified\"><magnitude><val>%d</val></magnitude></node>", cladeCount);
        }
        kronaReport(FP, taxDB, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode* taxon = taxDB.taxonNode(taxID);
        std::string escapedName = escapeAttribute(taxon->name);
        fprintf(FP, "<node name=\"%s\"><magnitude><val>%d</val></magnitude>", escapedName.c_str(), cladeCount);
        std::vector<TaxID> children = it->second.children;
        std::sort(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts,b); });
        for (TaxID childTaxId : children) {
            if (cladeCounts.count(childTaxId)) {
                kronaReport(FP, taxDB, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
        fprintf(FP, "</node>");
    }
}

int taxonomyreport(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // 1. Read taxonomy
    NcbiTaxonomy * taxDB = NcbiTaxonomy::openTaxonomy(par.db1);

    std::vector<std::pair<unsigned int, unsigned int>> mapping;
    if(FileUtil::fileExists(std::string(par.db1 + "_mapping").c_str()) == false){
        Debug(Debug::ERROR) << par.db1 + "_mapping" << " does not exist. Please create the taxonomy mapping!\n";
        EXIT(EXIT_FAILURE);
    }
    bool isSorted = Util::readMapping( par.db1 + "_mapping", mapping);
    if(isSorted == false){
        std::stable_sort(mapping.begin(), mapping.end(), compareToFirstInt);
    }

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // TODO: Better way to get file specified by param3?
    FILE *resultFP = fopen(par.db3.c_str(), "w");


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
    unsigned int unknownCnt = (taxCounts.find(0) != taxCounts.end()) ? taxCounts.at(0) : 0;
    Debug(Debug::INFO) << unknownCnt << " reads are unclassified.\n";

    std::unordered_map<TaxID, TaxonCounts> cladeCounts = taxDB->getCladeCounts(taxCounts);
    if (par.reportMode == 0) {
        taxReport(resultFP, *taxDB, cladeCounts, reader.getSize());
    } else {
        fwrite(krona_prelude_html, krona_prelude_html_len, sizeof(char), resultFP);
        fprintf(resultFP, "<node name=\"all\"><magnitude><val>%zu</val></magnitude>", reader.getSize());
        kronaReport(resultFP, *taxDB, cladeCounts, reader.getSize());
        fprintf(resultFP, "</node></krona></div></body></html>");
    }
    delete taxDB;
    reader.close();
    return EXIT_SUCCESS;
}

