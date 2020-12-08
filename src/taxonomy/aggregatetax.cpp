#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include <map>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

int aggregate(const bool useAln, int argc, const char **argv, const Command& command) {
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

    // open alignment per sequence - will be used only if useAln
    DBReader<unsigned int>* alnSeqReader = NULL;
    if (useAln == true) {
        alnSeqReader = new DBReader<unsigned int>(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        alnSeqReader->open(DBReader<unsigned int>::NOSORT);
    }

    // output is either db4 or db5
    std::string outDbStr = par.db4;
    std::string outDbIndexStr = par.db4Index;
    if (useAln == true) {
        outDbStr = par.db5;
        outDbIndexStr = par.db5Index;
    }

    DBWriter writer(outDbStr.c_str(), outDbIndexStr.c_str(), par.threads, par.compressed, Parameters::DBTYPE_TAXONOMICAL_RESULT);
    writer.open();

    std::vector<std::string> ranks = NcbiTaxonomy::parseRanks(par.lcaRanks);

    Debug::Progress progress(setToSeqReader.getSize());

    #pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        // per thread variables
        const char *entry[255];
        std::vector<WeightedTaxHit> setTaxa;

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

                size_t seqId = taxSeqReader.getId(seqKey);
                if (seqId == UINT_MAX) {
                    Debug(Debug::ERROR) << "Missing key " << seqKey << " in tax result\n";
                    EXIT(EXIT_FAILURE);
                }
                char *seqToTaxData = taxSeqReader.getData(seqId, thread_idx);
                TaxID taxon = Util::fast_atoi<int>(seqToTaxData);

                float evalue = FLT_MAX;
                if (useAln == true && taxon != 0) {
                    size_t alnId = alnSeqReader->getId(seqKey);
                    if (alnId == UINT_MAX) {
                        Debug(Debug::ERROR) << "Missing key " << alnId << " in alignment result\n";
                        EXIT(EXIT_FAILURE);
                    }
                    char *seqToAlnData = alnSeqReader->getData(alnId, thread_idx);
                    size_t numCols = Util::getWordsOfLine(seqToAlnData, entry, 255);
                    if (numCols < Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                        Debug(Debug::ERROR) << "No alignment result for taxon " << taxon << " found\n";
                        EXIT(EXIT_FAILURE);
                    }
                    evalue = strtod(entry[3], NULL);
                }

                setTaxa.emplace_back(taxon, evalue);
                results = Util::skipLine(results);
            }

            // aggregate - the counters will be filled by the selection function:
            WeightedTaxResult result = t->weightedMajorityLCA(setTaxa, par.majorityThr);
            TaxonNode const * node = t->taxonNode(result.taxon, false);

            size_t totalNumSeqs = result.assignedSeqs + result.unassignedSeqs;
            
            // prepare write
            if ((result.taxon == 0) || (node == NULL)) {
                setTaxStr.append(SSTR(0));
                setTaxStr.append(1, '\t');
                setTaxStr.append("no rank");
                setTaxStr.append(1, '\t');
                setTaxStr.append("unclassified");
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(totalNumSeqs));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(result.assignedSeqs));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(result.seqsAgreeWithSelectedTaxon));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(roundf(result.selectedPercent * 100) / 100));

                if (!ranks.empty()) {
                    setTaxStr += '\t';
                }
                if (par.showTaxLineage > 0) {
                    setTaxStr += '\t';
                }
            } else {
                setTaxStr.append(SSTR(node->taxId));
                setTaxStr.append(1, '\t');
                setTaxStr.append(node->rank);
                setTaxStr.append(1, '\t');
                setTaxStr.append(node->name);
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(totalNumSeqs));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(result.assignedSeqs));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(result.seqsAgreeWithSelectedTaxon));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(roundf(result.selectedPercent * 100) / 100));

                if (!ranks.empty()) {
                    std::string lcaRanks = Util::implode(t->AtRanks(node, ranks), ';');
                    setTaxStr += '\t' + lcaRanks;
                }
                if (par.showTaxLineage == 1) {
                    setTaxStr += '\t' + t->taxLineage(node, true);
                }
                if (par.showTaxLineage == 2) {
                    setTaxStr += '\t' + t->taxLineage(node, false);
                }
            }
            setTaxStr += '\n';

            writer.writeData(setTaxStr.c_str(), setTaxStr.size(), setKey, thread_idx);
            setTaxStr.clear();

            // ready to move to the next set
            setTaxa.clear();
        }
    }

    writer.close();
    taxSeqReader.close();
    setToSeqReader.close();
    if (alnSeqReader != NULL) {
        alnSeqReader->close();
        delete alnSeqReader;
    }
    delete t;

    return EXIT_SUCCESS;

}

int aggregatetaxweights(int argc, const char **argv, const Command& command) {
    return aggregate(true, argc, argv, command);
}

int aggregatetax(int argc, const char **argv, const Command& command) {
    return aggregate(false, argc, argv, command);
}
