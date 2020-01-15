#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include <map>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

TaxID selectTaxForSet (const std::vector<TaxID> &setTaxa, NcbiTaxonomy const *taxonomy, const float majorityCutoff) {
    // count num occurences of each ancestor 
    std::map<TaxID,unsigned int> ancTaxIdsCounts;
    size_t totalAssignedSeqs = 0;

    for (size_t i = 0; i < setTaxa.size(); ++i) {
        TaxID currTaxId = setTaxa[i];
        // ignore unassigned sequences
        if (currTaxId == 0) {
            continue;
        }
        TaxonNode const * node = taxonomy->taxonNode(currTaxId, false);
        if (node == NULL) {
            Debug(Debug::ERROR) << "taxonid: " << currTaxId << " does not match a legal taxonomy node.\n";
            EXIT(EXIT_FAILURE);
        }
        totalAssignedSeqs++;

        // add count
        if (ancTaxIdsCounts.find(currTaxId) != ancTaxIdsCounts.end()) {
            ancTaxIdsCounts[currTaxId]++;
        } else {
            ancTaxIdsCounts.insert(std::pair<TaxID,unsigned int>(currTaxId,1));
        }

        // iterate all ancestors up to the root
        TaxID currParentTaxId = node->parentTaxId;
        while (currParentTaxId != currTaxId) {
            TaxonNode const * node = taxonomy->taxonNode(currParentTaxId, false);
            currTaxId = currParentTaxId;
            if (ancTaxIdsCounts.find(currTaxId) != ancTaxIdsCounts.end()) {
                ancTaxIdsCounts[currTaxId]++;
            } else {
                ancTaxIdsCounts.insert(std::pair<TaxID,unsigned int>(currTaxId,1));
            }
            currParentTaxId = node->parentTaxId;
        }
    }

    // select the lowest ancestor that meets the cutoff
    int minRank = INT_MAX;
    TaxID selctedTaxon = 0;
    float selectedPercent = 0;

    for (std::map<TaxID,unsigned int>::iterator it = ancTaxIdsCounts.begin(); it != ancTaxIdsCounts.end(); it++) {
        float currPercent = float(it->second) / totalAssignedSeqs;
        if (currPercent >= majorityCutoff) {
            TaxID currTaxId = it->first;
            TaxonNode const * node = taxonomy->taxonNode(currTaxId, false);
            int currRankInd = NcbiTaxonomy::findRankIndex(node->rank);
            if (currRankInd > 0) {
                if ((currRankInd < minRank) || ((currRankInd == minRank) && (currPercent > selectedPercent))) {
                    selctedTaxon = currTaxId;
                    minRank = currRankInd;
                    selectedPercent = currPercent;
                }
            }
        }
    }

    return (selctedTaxon);
}

int aggregatetax(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // open taxonomy - evolutionary relationships amongst taxa
    NcbiTaxonomy * t = NcbiTaxonomy::openTaxonomy(par.db1);
    
    // open mapping of set to sequence
    DBReader<unsigned int> setToSeqReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    setToSeqReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    // open tax assignments per sequence
    DBReader<unsigned int> taxSeqReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    taxSeqReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_TAXONOMICAL_RESULT);
    writer.open();

    std::vector<std::string> ranks = NcbiTaxonomy::parseRanks(par.lcaRanks);

    Debug::Progress progress(taxSeqReader.getSize());

    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        // per thread variables
        const char *entry[2048];
        std::vector<TaxID> setTaxa;
        std::string setTaxStr;
        setTaxStr.reserve(4096);

        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < setToSeqReader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int setKey = setToSeqReader.getDbKey(i);

            char *results = setToSeqReader.getData(i, thread_idx);

            // process a specific set
            while (*results != '\0') {
                Util::getWordsOfLine(results, entry, 255);
                unsigned int seqKey = Util::fast_atoi<unsigned int>(entry[0]);

                char *seqToTaxData = taxSeqReader.getDataByDBKey(seqKey, thread_idx);
                Util::getWordsOfLine(seqToTaxData, entry, 255);
                TaxID taxon = Util::fast_atoi<int>(entry[0]);

                setTaxa.emplace_back(taxon);
                results = Util::skipLine(results);
            }

            // aggregate
            TaxID setSelectedTaxon = selectTaxForSet(setTaxa, t, par.majorityThr);
            TaxonNode const * node = t->taxonNode(setSelectedTaxon, false);
            
            // prepare write
            if ((setSelectedTaxon == 0) || (node == NULL)) {
                setTaxStr = "0\tno rank\tunclassified";
                if (!ranks.empty()) {
                    setTaxStr += '\t';
                }
                if (par.showTaxLineage) {
                    setTaxStr += '\t';
                }
            } else {
                setTaxStr = SSTR(node->taxId) + '\t' + node->rank + '\t' + node->name;
                if (!ranks.empty()) {
                    std::string lcaRanks = Util::implode(t->AtRanks(node, ranks), ';');
                    setTaxStr += '\t' + lcaRanks;
                }
                if (par.showTaxLineage) {
                    setTaxStr += '\t' + t->taxLineage(node);
                }
            }
            setTaxStr += '\n';

            writer.writeData(setTaxStr.c_str(), setTaxStr.size(), setKey, thread_idx);
            setTaxStr.clear();

            // ready to move to the next set
            setTaxa.clear();
        }
    };
    Debug(Debug::INFO) << "\n";

    writer.close();
    taxSeqReader.close();
    setToSeqReader.close();
    delete t;

    return EXIT_SUCCESS;
}
