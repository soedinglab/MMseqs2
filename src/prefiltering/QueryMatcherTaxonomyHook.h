#ifndef QUERY_MATCHER_TAXONOMY_HOOK_H
#define QUERY_MATCHER_TAXONOMY_HOOK_H

#include "QueryMatcher.h"
#include "NcbiTaxonomy.h"
#include "MappingReader.h"
#include "DBReader.h"
#include "TaxonomyExpression.h"

#ifdef OPENMP
#include <omp.h>
#endif

class QueryMatcherTaxonomyHook : public QueryMatcherHook {
public:
    QueryMatcherTaxonomyHook(std::string targetPath, DBReader<unsigned int>* targetReader, const std::string& expressionString, unsigned int threads)
        : targetReader(targetReader), dbFrom(0), threads(threads) {
        std::string targetName = dbPathWithoutIndex(targetPath);
        taxonomy = NcbiTaxonomy::openTaxonomy(targetName);
        taxonomyMapping = new MappingReader(targetName);
        expression = new TaxonomyExpression*[threads];
        for (unsigned int i = 0; i < threads; i++) {
            expression[i] = new TaxonomyExpression(expressionString, *taxonomy);
        }
    }

    ~QueryMatcherTaxonomyHook() {
        delete taxonomy;
        delete taxonomyMapping;
        for (unsigned int i = 0; i < threads; i++) {
            delete expression[i];
        }
        delete[] expression;
    }

    void setDbFrom(unsigned int from) {
        dbFrom = from;
    }

    size_t afterDiagonalMatchingHook(QueryMatcher& matcher, size_t resultSize) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        size_t writePos = 0;
        for (size_t i = 0; i < resultSize; i++) {
            unsigned int currId = matcher.foundDiagonals[i].id;
            unsigned int key = targetReader->getDbKey(dbFrom + currId);
            TaxID currTax = taxonomyMapping->lookup(key);
            if (expression[thread_idx]->isAncestor(currTax)) {
                if (i != writePos) {
                    matcher.foundDiagonals[writePos] = matcher.foundDiagonals[i];
                }
                writePos++;
            }
        }
        return writePos;
    }

    static std::string dbPathWithoutIndex(const std::string& dbname) {
        static const std::vector<std::string> suffices = {
            "_ss.idx",
            "_ss.linidx",
            "_ss",
            ".idx",
            ".linidx"
        };
        for (size_t i = 0; i < suffices.size(); ++i) {
            size_t lastpos = dbname.rfind(suffices[i]);
            if (lastpos != std::string::npos && dbname.size() - lastpos == suffices[i].length()){
                return dbname.substr(0, lastpos);
            }
        }
        return dbname;
    }

    NcbiTaxonomy* taxonomy;
    MappingReader* taxonomyMapping;
    DBReader<unsigned int>* targetReader;
    TaxonomyExpression** expression;

    unsigned int dbFrom;
    unsigned int threads;
};

#endif