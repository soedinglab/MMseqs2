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

const double MAX_WEIGHT = 1000;
const TaxID ROOT_TAXID = 1;
const int ROOT_RANK = INT_MAX;

struct taxHit {
    void setByEntry(const TaxID & taxonInput, const bool useAln, const char ** taxHitData, const size_t numCols, const int voteMode) {
        taxon = taxonInput;
        evalue = 1.0;
        weight = 0.0;

        // if voteMode is evalue-based, all tax-assigned sequences should have alignment info...
        if ((taxon != 0) && (numCols < Matcher::ALN_RES_WITHOUT_BT_COL_CNT) && (useAln == true)) {
            Debug(Debug::ERROR) << "voteMode is evalue-based but taxonid: " << taxon << " does not have alignment info.\n";
            EXIT(EXIT_FAILURE);
        }

        // extract from alignment info
        if (useAln == true) {
            evalue = strtod(taxHitData[3],NULL);
        }

        // update weight according to mode
        if (voteMode == Parameters::AGG_TAX_UNIFORM) {
            weight = 1.0;
        } else if (voteMode == Parameters::AGG_TAX_MINUS_LOG_EVAL) {
            if (evalue > 0) {
                weight = -log(evalue);
            } else {
                weight = MAX_WEIGHT;
            }
        }
    }

    TaxID taxon;
    double evalue;
    double weight;
};

struct taxNode {
    void set(const double weightInput, const bool isCandidateInput, const TaxID & childTaxonInput) {
        weight = weightInput;
        isCandidate = isCandidateInput;
        childTaxon = childTaxonInput;
    }

    void update(const double weightToAdd, const TaxID & childTaxonInput) {
        if (childTaxon != childTaxonInput) {
            isCandidate = true;
            childTaxon = childTaxonInput;
        }
        weight += weightToAdd;
    }

    // these will be filled when iterating over all contributing lineages
    double weight;
    bool isCandidate;
    TaxID childTaxon;
};

TaxID selectTaxForSet (const std::vector<taxHit> &setTaxa, NcbiTaxonomy const *taxonomy, const float majorityCutoff, 
                        size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent) {
    // count num occurences of each ancestor, possibly weighted 
    std::map<TaxID,taxNode> ancTaxIdsCounts;

    // initialize counters and weights
    numAssignedSeqs = 0;
    numUnassignedSeqs = 0;
    numSeqsAgreeWithSelectedTaxon = 0;
    selectedPercent = 0;
    double totalAssignedSeqsWeights = 0.0;

    for (size_t i = 0; i < setTaxa.size(); ++i) {
        TaxID currTaxId = setTaxa[i].taxon;
        double currWeight = setTaxa[i].weight;
        // ignore unassigned sequences
        if (currTaxId == 0) {
            numUnassignedSeqs++;
            continue;
        }
        TaxonNode const * node = taxonomy->taxonNode(currTaxId, false);
        if (node == NULL) {
            Debug(Debug::ERROR) << "taxonid: " << currTaxId << " does not match a legal taxonomy node.\n";
            EXIT(EXIT_FAILURE);
        }
        totalAssignedSeqsWeights += currWeight;
        numAssignedSeqs++;

        // each start of a path due to an orf is a candidate
        if (ancTaxIdsCounts.find(currTaxId) != ancTaxIdsCounts.end()) {
            ancTaxIdsCounts[currTaxId].update(currWeight, 0);
        } else {
            taxNode currNode;
            currNode.set(currWeight, true, 0);
            ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currTaxId, currNode));
        }

        // iterate all ancestors up to root (including). add currWeight and candidate status to each
        TaxID currParentTaxId = node->parentTaxId;
        while (currParentTaxId != currTaxId) {
            if (ancTaxIdsCounts.find(currParentTaxId) != ancTaxIdsCounts.end()) {
                ancTaxIdsCounts[currParentTaxId].update(currWeight, currTaxId);
            } else {
                taxNode currParentNode;
                currParentNode.set(currWeight, false, currTaxId);
                ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currParentTaxId, currParentNode));
            }
            // move up:
            currTaxId = currParentTaxId;
            node = taxonomy->taxonNode(currParentTaxId, false);
            currParentTaxId = node->parentTaxId;
        }     
    }

    // select the lowest ancestor that meets the cutoff
    int minRank = INT_MAX;
    TaxID selctedTaxon = 0;

    for (std::map<TaxID,taxNode>::iterator it = ancTaxIdsCounts.begin(); it != ancTaxIdsCounts.end(); it++) {
        // consider only candidates:
        if (it->second.isCandidate == false) {
            continue;
        }

        double currPercent = float(it->second.weight) / totalAssignedSeqsWeights;
        if (currPercent >= majorityCutoff) {
            // iterate all ancestors to find lineage min rank (the candidate is a descendant of a node with this rank)
            TaxID currTaxId = it->first;
            TaxonNode const * node = taxonomy->taxonNode(currTaxId, false);
            int currMinRank = ROOT_RANK;
            TaxID currParentTaxId = node->parentTaxId;
            while (currParentTaxId != currTaxId) {
                int currRankInd = NcbiTaxonomy::findRankIndex(node->rank);
                if ((currRankInd > 0) && (currRankInd < currMinRank)) {
                    currMinRank = currRankInd;
                    // the rank can only go up on the way to the root, so we can break
                    break;
                }
                // move up:
                currTaxId = currParentTaxId;
                node = taxonomy->taxonNode(currParentTaxId, false);
                currParentTaxId = node->parentTaxId;
            }

            if ((currMinRank < minRank) || ((currMinRank == minRank) && (currPercent > selectedPercent))) {
                selctedTaxon = it->first;
                minRank = currMinRank;
                selectedPercent = currPercent;
            }
        }
    }

    // count the number of seqs who have selectedTaxon in their ancestors (agree with selection):
    if (selctedTaxon == ROOT_TAXID) {
        // all agree with "root"
        numSeqsAgreeWithSelectedTaxon = numAssignedSeqs;
        return (selctedTaxon);
    }
    if (selctedTaxon == 0) {
        // nothing informative
        return (selctedTaxon);
    }
    // otherwise, iterate over all seqs
    for (size_t i = 0; i < setTaxa.size(); ++i) {
        TaxID currTaxId = setTaxa[i].taxon;
        // ignore unassigned sequences
        if (currTaxId == 0) {
            continue;
        }
        TaxonNode const * node = taxonomy->taxonNode(currTaxId, false);

        // iterate all ancestors up to the root
        TaxID currParentTaxId = node->parentTaxId;
        while (currParentTaxId != currTaxId) {
            if (currTaxId == selctedTaxon) {
                numSeqsAgreeWithSelectedTaxon++;
                break;
            }
            currTaxId = currParentTaxId;
            node = taxonomy->taxonNode(currParentTaxId, false);
            currParentTaxId = node->parentTaxId;
        }
    }

    return (selctedTaxon);
}

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
    } else if (par.voteMode == Parameters::AGG_TAX_MINUS_LOG_EVAL) {
        Debug(Debug::ERROR) << "voteMode is evalue-based but no alignment database was provided. consider calling aggregatetaxweights\n";
        EXIT(EXIT_FAILURE);
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
        const char *entry[2048];
        std::vector<taxHit> setTaxa;
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
                Util::getWordsOfLine(seqToTaxData, entry, 255);
                TaxID taxon = Util::fast_atoi<int>(entry[0]);
                size_t numCols = 0;

                if (useAln == true) {
                    size_t alnId = alnSeqReader->getId(seqKey);
                    if (alnId == UINT_MAX) {
                        Debug(Debug::ERROR) << "Missing key " << alnId << " in alignment result\n";
                        EXIT(EXIT_FAILURE);
                    }
                    char *seqToAlnData = alnSeqReader->getData(alnId, thread_idx);
                    numCols = Util::getWordsOfLine(seqToAlnData, entry, 255);
                }

                taxHit currTaxHit;
                currTaxHit.setByEntry(taxon, useAln, entry, numCols, par.voteMode);
                setTaxa.emplace_back(currTaxHit);
                results = Util::skipLine(results);
            }

            // aggregate - the counters will be filled by the selection function:
            size_t numAssignedSeqs = 0;
            size_t numUnassignedSeqs = 0;
            size_t numSeqsAgreeWithSelectedTaxon = 0;
            double selectedPercent = 0;

            TaxID setSelectedTaxon = selectTaxForSet(setTaxa, t, par.majorityThr, numAssignedSeqs, numUnassignedSeqs, numSeqsAgreeWithSelectedTaxon, selectedPercent);
            TaxonNode const * node = t->taxonNode(setSelectedTaxon, false);

            size_t totalNumSeqs = numAssignedSeqs + numUnassignedSeqs;
            
            // prepare write
            if ((setSelectedTaxon == 0) || (node == NULL)) {
                setTaxStr.append(SSTR(0));
                setTaxStr.append(1, '\t');
                setTaxStr.append("no rank");
                setTaxStr.append(1, '\t');
                setTaxStr.append("unclassified");
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(totalNumSeqs));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(numAssignedSeqs));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(numSeqsAgreeWithSelectedTaxon));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(roundf(selectedPercent * 100) / 100));

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
                setTaxStr.append(SSTR(numAssignedSeqs));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(numSeqsAgreeWithSelectedTaxon));
                setTaxStr.append(1, '\t');
                setTaxStr.append(SSTR(roundf(selectedPercent * 100) / 100));

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
