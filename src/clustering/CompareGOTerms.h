//
// Created by lars on 12.04.15.
//


#ifndef MMSEQS_COMPAREGOTERMS_H
#define MMSEQS_COMPAREGOTERMS_H

#include <string>

#include "DBReader.h"
#include "Debug.h"

class CompareGOTerms{

    ~CompareGOTerms();

public:

    void init();
    CompareGOTerms(std::string go_ffindex,std::string go_ffindex_indexfile,std::string protid_go_ffindex,std::string protid_go_ffindex_indexfile);

private:

    DBReader* go_ffindex_reader;

    DBReader* protid_go_ffindex_reader;


};

#endif //MMSEQS_COMPAREGOTERMS_H
