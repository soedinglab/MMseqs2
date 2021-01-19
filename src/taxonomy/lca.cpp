#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

static bool compareToFirstInt(const std::pair<unsigned int, unsigned int>& lhs, const std::pair<unsigned int, unsigned int>&  rhs){
    return (lhs.first <= rhs.first);
}

int dolca(int argc, const char **argv, const Command& command, bool majority) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    NcbiTaxonomy * t = NcbiTaxonomy::openTaxonomy(par.db1);

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

    if (majority) {
        if (par.voteMode != Parameters::AGG_TAX_UNIFORM && Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_CLUSTER_RES)) {
            Debug(Debug::WARNING) << "Cluster input can only be used with --vote-mode 0\nContinuing with --vote-mode 0\n";
            par.voteMode = Parameters::AGG_TAX_UNIFORM;
        } else if (par.voteMode == Parameters::AGG_TAX_MINUS_LOG_EVAL && (Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_PREFILTER_RES) || Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_PREFILTER_REV_RES))) {
            Debug(Debug::WARNING) << "Prefilter input can only be used with --vote-mode 0 or 2\nContinuing with --vote-mode 0\n";
            par.voteMode = Parameters::AGG_TAX_UNIFORM;
        }
    }

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_TAXONOMICAL_RESULT);
    writer.open();

    std::vector<std::string> ranks = NcbiTaxonomy::parseRanks(par.lcaRanks);

    // a few NCBI taxa are blacklisted by default, they contain unclassified sequences (e.g. metagenomes) or other sequences (e.g. plasmids)
    // if we do not remove those, a lot of sequences would be classified as Root, even though they have a sensible LCA
    std::vector<TaxID> blacklist;
    std::vector<std::string> splits = Util::split(par.blacklist, ",");
    for (size_t i = 0; i < splits.size(); ++i) {
        TaxID taxon = Util::fast_atoi<int>(splits[i].c_str());
        if (taxon == 0) {
            Debug(Debug::WARNING) << "Cannot block root taxon 0\n";
            continue;
        }
        if (t->nodeExists(taxon) == false) {
            Debug(Debug::WARNING) << "Ignoring missing blocked taxon " << taxon << "\n";
            continue;
        }

        const char *split;
        if ((split = strchr(splits[i].c_str(), ':')) != NULL) {
            const char* name = split + 1;
            const TaxonNode* node = t->taxonNode(taxon, false);
            if (node == NULL) {
                Debug(Debug::WARNING) << "Ignoring missing blocked taxon " << taxon << "\n";
                continue;
            }
            const char* nodeName = t->getString(node->nameIdx);
            if (strcmp(nodeName, name) != 0) {
                Debug(Debug::WARNING) << "Node name '" << name << "' does not match to be blocked name '" << nodeName << "'\n";
                continue;
            }
        }
        blacklist.emplace_back(taxon);
    }

    // will be used when no hits
    std::string noTaxResult = "0\tno rank\tunclassified";
    if (!ranks.empty()) {
        noTaxResult += '\t';
    }
    if (par.showTaxLineage > 0) {
        noTaxResult += '\t';
    }
    noTaxResult += '\n';


    size_t taxonNotFound = 0;
    size_t found = 0;
    Debug::Progress progress(reader.getSize());
    #pragma omp parallel
    {
        const char *entry[255];
        std::string result;
        result.reserve(4096);
        unsigned int thread_idx = 0;

#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        #pragma omp for schedule(dynamic, 10) reduction (+:taxonNotFound, found)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getEntryLen(i);

            std::vector<int> taxa;
            std::vector<WeightedTaxHit> weightedTaxa;
            while (*data != '\0') {
                TaxID taxon;
                unsigned int id;
                std::pair<unsigned int, unsigned int> val;
                std::vector<std::pair<unsigned int, unsigned int>>::iterator mappingIt;
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                data = Util::skipLine(data);
                if (columns == 0) {
                    Debug(Debug::WARNING) << "Empty entry: " << i << "!";
                    continue;
                }

                id = Util::fast_atoi<unsigned int>(entry[0]);
                val.first = id;
                mappingIt = std::upper_bound(mapping.begin(), mapping.end(), val, compareToFirstInt);

                if (mappingIt == mapping.end() || mappingIt->first != val.first) {
                    // TODO: Check which taxa were not found
                    taxonNotFound += 1;
                    continue;
                }
                found++;
                taxon = mappingIt->second;

                // remove blacklisted taxa
                bool isBlacklisted = false;
                for (size_t j = 0; j < blacklist.size(); ++j) {
                    if (blacklist[j] == 0) {
                        continue;
                    }
                    if (t->IsAncestor(blacklist[j], taxon)) {
                        isBlacklisted = true;
                        break;
                    }
                }

                if (isBlacklisted == false) {
                    if (majority) {
                        float weight = FLT_MAX;
                        if (par.voteMode == Parameters::AGG_TAX_MINUS_LOG_EVAL) {
                            if (columns <= 3) {
                                Debug(Debug::ERROR) << "No alignment result for taxon " << taxon << " found\n";
                                EXIT(EXIT_FAILURE);
                            }
                            weight = strtod(entry[3], NULL);
                        } else if (par.voteMode == Parameters::AGG_TAX_SCORE) {
                            if (columns <= 1) {
                                Debug(Debug::ERROR) << "No alignment result for taxon " << taxon << " found\n";
                                EXIT(EXIT_FAILURE);
                            }
                            weight = strtod(entry[1], NULL);
                        }
                        weightedTaxa.emplace_back(taxon, weight, par.voteMode);
                    } else {
                        taxa.emplace_back(taxon);
                    }
                }
            }

            if (length == 1) {
                writer.writeData(noTaxResult.c_str(), noTaxResult.size(), key, thread_idx);
                continue;
            }

            TaxonNode const * node = NULL;
            if (majority) {
                WeightedTaxResult result = t->weightedMajorityLCA(weightedTaxa, par.majorityThr);
                node = t->taxonNode(result.taxon, false);
            } else {
                node = t->LCA(taxa);
            }
            if (node == NULL) {
                writer.writeData(noTaxResult.c_str(), noTaxResult.size(), key, thread_idx);
                continue;
            }

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
            writer.writeData(result.c_str(), result.size(), key, thread_idx);
            result.clear();
        }
    }
    Debug(Debug::INFO) << "Taxonomy for " << taxonNotFound << " out of " << taxonNotFound+found << " entries not found\n";
    writer.close();
    reader.close();
    delete t;

    return EXIT_SUCCESS;
}

int lca(int argc, const char **argv, const Command& command) {
    return dolca(argc, argv, command, false);
}

int majoritylca(int argc, const char **argv, const Command& command) {
    return dolca(argc, argv, command, true);
}
