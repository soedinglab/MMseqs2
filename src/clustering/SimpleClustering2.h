//
// Created by lars on 28.07.15.
//

#ifndef MMSEQS_SIMPLECLUSTERING2_H
#define MMSEQS_SIMPLECLUSTERING2_H

#include <DBReader.h>
#include <DBWriter.h>
#include <set>
#include <list>
#include "SetElement.h"

class SimpleClustering2 {
public:
    SimpleClustering2(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage);

    std::list<set *>  execute();
    ~SimpleClustering2();
private:
    DBReader* seqDbr;

    DBReader* alnDbr;

    float seqIdThr;

    float coverage;
//datastructures



//methods


};

#endif //MMSEQS_SIMPLECLUSTERING2_H
