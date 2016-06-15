#ifndef DOMAIN_H
#define DOMAIN_H

#include <string>
#include <ostream>

struct Domain {
    std::string query;

    unsigned int qStart;
    unsigned int qEnd;
    unsigned int qLength;

    std::string target;

    unsigned int tStart;
    unsigned int tEnd;
    unsigned int tLength;

    double eValue;

    Domain(const std::string &query, unsigned int qStart, unsigned int qEnd, unsigned int qLength,
           const std::string &target, unsigned int tStart, unsigned int tEnd, unsigned int tLength, double eValue) :
            query(query), qStart(qStart), qEnd(qEnd), qLength(qLength),
            target(target), tStart(tStart), tEnd(tEnd), tLength(tLength), eValue(eValue) { }

    friend bool operator<(const Domain &h1, const Domain &h2) {
        return h1.eValue < h2.eValue;
    }

    void writeResult(std::ostream &out) const {
        const char sep = '\t';
        out << query << sep << target << sep << qStart << sep << qEnd << sep << qLength;
        out << sep << tStart << sep << tEnd << sep << tLength << sep << eValue;
    }
};

#endif
