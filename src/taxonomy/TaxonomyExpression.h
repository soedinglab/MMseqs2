//
// Created by Martin Steinegger on 2019-08-20.
//
#ifndef MMSEQS_TAXONOMYEXPRESSION_H
#define MMSEQS_TAXONOMYEXPRESSION_H

#include "NcbiTaxonomy.h"
#include "ExpressionParser.h"

#include <vector>
#include <ctype.h>

// class need one instance per thread
class TaxonomyExpression {
public:
    enum CommaMeaning {
        COMMA_IS_COMMA,
        COMMA_IS_OR,
        COMMA_IS_AND
    };

    TaxonomyExpression(const std::string &expression, NcbiTaxonomy &taxonomy, CommaMeaning commaMeaning = COMMA_IS_OR) {
        std::string bracketExpression;
        bool inNumber = false;
        for (size_t i = 0; i < expression.size(); i++) {
            // make brackets around numbers for tinyexpr
            const bool isDigit = isdigit(expression[i]);
            if (isDigit && inNumber == true) {
                bracketExpression.push_back(expression[i]);
            } else if (isDigit && inNumber == false) {
                bracketExpression.append("a(");
                bracketExpression.push_back(expression[i]);
                inNumber = true;
            } else {
                if (inNumber == true) {
                    bracketExpression.append(")");
                    inNumber = false;
                }
                if (commaMeaning != COMMA_IS_COMMA && expression[i] == ',') {
                    if (commaMeaning == COMMA_IS_OR) {
                        bracketExpression.append("||");
                    } else if (commaMeaning == COMMA_IS_AND) {
                        bracketExpression.append("&&");
                    }
                } else {
                    bracketExpression.push_back(expression[i]);
                }
            }
        }
        if (inNumber == true) {
            bracketExpression.append(")");
        }
        tc.t = &taxonomy;
        te_variable var;
        var.name = "a";
        // GCC 4.8 does not like casting functions to void*
        // GCC > 4.8 is fine with this
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
        var.address = (const void *) &acst;
#pragma GCC diagnostic pop
        var.type = TE_CLOSURE1;
        var.context = (void *) &tc;
        vars.push_back(var);
        parser = new ExpressionParser(bracketExpression.c_str(), vars);
    }

    ~TaxonomyExpression() {
        delete parser;
    }

    bool isAncestor(TaxID taxId) {
        tc.taxId = taxId;
        const double result = parser->evaluate();
        return (result != 0);
    }

private:
    struct TaxContext {
        NcbiTaxonomy *t;
        TaxID taxId;
    };
    TaxContext tc;
    ExpressionParser *parser;
    std::vector<te_variable> vars;

    static double acst(void *context, double a) {
        TaxContext *o = (TaxContext *) context;
        bool retVal = o->t->IsAncestor((TaxID) a, o->taxId);
        return (retVal) ? 1.0 : 0.0;
    }
};
#endif //MMSEQS_TAXONOMYEXPRESSION_H
