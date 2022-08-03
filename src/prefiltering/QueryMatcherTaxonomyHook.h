#ifndef QUERY_MATCHER_TAXONOMY_HOOK_H
#define QUERY_MATCHER_TAXONOMY_HOOK_H

#include "QueryMatcher.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"
#include "DBReader.h"
#include "PrefilteringIndexReader.h"
#include "TaxonomyExpression.h"

class QueryMatcherTaxonomyHook : public QueryMatcherHook {
public:
    QueryMatcherTaxonomyHook(std::string targetPath, DBReader<unsigned int>* targetReader, const std::string& expressionString)
        : targetReader(targetReader), dbFrom(0) {
        std::string targetName = PrefilteringIndexReader::dbPathWithoutIndex(targetPath);
        taxonomy = NcbiTaxonomy::openTaxonomy(targetName);
        taxonomyMapping = new MappingReader(targetName);
        expression = new TaxonomyExpression(expressionString, *taxonomy);
    }

    ~QueryMatcherTaxonomyHook() {
        delete taxonomy;
        delete taxonomyMapping;
        delete expression;
    }

    void setDbFrom(unsigned int from) {
        dbFrom = from;
    }

    size_t afterDiagonalMatchingHook(QueryMatcher& matcher, size_t resultSize) {
        size_t writePos = 0;
        for (size_t i = 0; i < resultSize; i++) {
            unsigned int currId = matcher.foundDiagonals[i].id;
            unsigned int key = targetReader->getDbKey(dbFrom + currId);
            TaxID currTax = taxonomyMapping->lookup(key);
            if (expression->isAncestor(currTax)) {
                if (i != writePos) {
                    matcher.foundDiagonals[writePos] = matcher.foundDiagonals[i];
                }
                writePos++;
            }
        }
        return writePos;
    }

    NcbiTaxonomy* taxonomy;
    MappingReader* taxonomyMapping;
    DBReader<unsigned int>* targetReader;
    TaxonomyExpression* expression;

    unsigned int dbFrom;
};

#endif