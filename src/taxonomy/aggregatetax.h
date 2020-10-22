//
// Created by Jaebeom Kim on 22/10/2020.
//

#ifndef ADCLASSIFIER2_AGGREGATETAX_H
#define ADCLASSIFIER2_AGGREGATETAX_H
#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "Matcher.h"
#include <map>
#include <algorithm>

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
        if (childTaxon != childTaxonInput) { //isCandidate가 뭐야??
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
                       size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent);
TaxID selectLcaFromTaxIdList(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
                             size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent);
int aggregate(const bool useAln, int argc, const char **argv, const Command& command);
int aggregatetaxweights(int argc, const char **argv, const Command& command);
int aggregatetax(int argc, const char **argv, const Command& command);

#endif //ADCLASSIFIER2_AGGREGATETAX_H
