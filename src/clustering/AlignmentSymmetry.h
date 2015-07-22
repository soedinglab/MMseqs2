//
// Created by lars on 10.06.15.
//

#ifndef MMSEQS_ALIGNMENTSYMMETRY_H
#define MMSEQS_ALIGNMENTSYMMETRY_H
#include <DBReader.h>
#include <DBWriter.h>
#include <bits/stl_set.h>
#include <bits/stl_list.h>
#include "SetElement.h"

class AlignmentSymmetry {
public:
    AlignmentSymmetry(DBReader * seqDbr, DBReader * alnDbr, DBWriter* alnWr, float seqIdThr, float coverage);

    void  execute();
private:
    DBReader* seqDbr;

    DBReader* alnDbr;

    DBWriter* alnWr;

    float seqIdThr;

    float coverage;
//datastructures

    unsigned int dbSize;




};
#endif //MMSEQS_ALIGNMENTSYMMETRY_H
