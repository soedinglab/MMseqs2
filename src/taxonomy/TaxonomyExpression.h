//
// Created by Martin Steinegger on 2019-08-20.
//
#ifndef MMSEQS_TAXONOMYEXPRESSION_H
#define MMSEQS_TAXONOMYEXPRESSION_H

#include "NcbiTaxonomy.h"
#include "Debug.h"
#include <vector>
#include <ctype.h>
#include "ExpressionParser.h"

// class need one instance per thread
class TaxonomyExpression{

private:
    struct TaxContext {
        NcbiTaxonomy* t;
        TaxID taxId;
    };
    TaxContext tc;
    ExpressionParser * parser;
    std::vector<te_variable> vars;

public:

    static double acst(void* context, double a) {
        TaxContext* o = (TaxContext*)context;
        bool retVal = o->t->IsAncestor((TaxID)a, o->taxId);
        return (retVal) ? 1.0 : 0.0;
    }

    TaxonomyExpression(std::string expression, NcbiTaxonomy &taxonomy){
        std::string bracketExpression;
        bool inNumber = false;
        // make brackets around numbers for tinyexpr
        for(size_t i = 0; i< expression.size(); i++){
            if(isdigit(expression[i]) && inNumber == true){
                bracketExpression.push_back(expression[i]);
            }else if(isdigit(expression[i]) && inNumber == false){
                bracketExpression.append("a(");
                bracketExpression.push_back(expression[i]);
                inNumber=true;
            } else if(inNumber == true) {
                bracketExpression.append(")");
                bracketExpression.push_back(expression[i]);
                inNumber=false;
            }else{
                bracketExpression.push_back(expression[i]);
            }
        }
        if(inNumber == true){
            bracketExpression.append(")");
        }
        tc.t = &taxonomy;
        te_variable var;
        var.name = "a";
        // GCC 4.8 does not like casting functions to void*
        // GCC > 4.8 are fine with this
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
        var.address = (const void *) &acst;
#pragma GCC diagnostic pop
        var.type = TE_CLOSURE1;
        var.context = (void *) &tc;
        vars.push_back(var);
        parser=new ExpressionParser(bracketExpression.c_str(), vars);
    }

    ~TaxonomyExpression(){
        delete parser;
    }

    bool isAncestor(TaxID taxId){
        tc.taxId = taxId;
        const double result = parser->evaluate();
        return (result != 0);
    }
};
#endif //MMSEQS_TAXONOMYEXPRESSION_H
