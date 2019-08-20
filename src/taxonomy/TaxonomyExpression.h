//
// Created by Martin Steinegger on 2019-08-20.
//
#ifndef MMSEQS_TAXONOMYEXPRESSION_H
#define MMSEQS_TAXONOMYEXPRESSION_H

#include "NcbiTaxonomy.h"
#include "Debug.h"
#include <vector>

class TaxonomyExpression{
    struct TaxonomyTerm{
        unsigned int taxId;
        bool shouldBeAncestor;
        TaxonomyTerm( unsigned int taxId, bool negaitve)
                :taxId(taxId), shouldBeAncestor(negaitve){}
    };
    std::vector<std::vector<TaxonomyTerm>> taxTerms;

public:
    TaxonomyExpression(std::string expression){
        bool inBracket = false;
        bool shouldBeAncestor = true;
        unsigned int taxId;
        std::vector<TaxonomyTerm> term;
        for(size_t pos = 0; pos < expression.size(); pos++){
            switch(expression[pos]) {
                case '(':
                    if(inBracket == true){
                        Debug(Debug::ERROR) << "Error in expression " << expression << ". It is not allowed to open another bracket within a bracket terms\n";
                        EXIT(EXIT_FAILURE);
                    }
                    inBracket = true;
                    break;
                case ')':
                    inBracket = false;
                    break;
                case '!':
                    shouldBeAncestor = false;
                    break;
                case '&':
                    shouldBeAncestor = true;
                    break;
                case ',':
                    shouldBeAncestor = true;
                    if(inBracket == false){
                        taxTerms.push_back(term);
                        term.clear();
                    }else{
                        Debug(Debug::ERROR) << "Error in expression " << expression << ". It is not allowed to use , within bracket terms only use &\n";
                        EXIT(EXIT_FAILURE);
                    }
                    break;
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                    taxId = Util::fast_atoi<unsigned int>(&expression[pos]);
                    while(pos < expression.size() && expression[pos] >= '0' && expression[pos] <= '9'){
                        pos++;
                    }
                    pos--;
                    term.emplace_back(taxId, shouldBeAncestor);
                    shouldBeAncestor = true;
                    break;
                default:
                    Debug(Debug::ERROR) << "Wrong character in expression: " << expression[pos] << "\n";
                    EXIT(EXIT_FAILURE);
                    break;
            }
        }
        taxTerms.push_back(term);

    }

    bool isAncestor(NcbiTaxonomy &taxonomy, std::vector<TaxonomyTerm> &termToCheck, unsigned int taxId){
        size_t ancestorCnt = 0;
        for (size_t j = 0; j < termToCheck.size(); ++j) {
            ancestorCnt += (taxonomy.IsAncestor(termToCheck[j].taxId, taxId) == termToCheck[j].shouldBeAncestor);
        }
        return (ancestorCnt == termToCheck.size());
    }
    // this function returns the index of the term that fulfils the criteria
    // -1 means no term fulfils the criteria
    int isAncestorOf(NcbiTaxonomy &taxonomy, unsigned int taxId){
        int index = -1;
        bool ancestor = false;
        for (size_t j = 0; j < taxTerms.size() && !ancestor; ++j) {
            ancestor |= isAncestor(taxonomy, taxTerms[j], taxId);
            index = (ancestor) ? j : index;
        }
        return index;
    }


    std::vector<std::vector<TaxonomyTerm>> getTaxTerms(){
        return taxTerms;
    }

};
#endif //MMSEQS_TAXONOMYEXPRESSION_H
