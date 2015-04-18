//
// Created by lars on 12.04.15.
//

#include "CompareGOTerms.h"


CompareGOTerms::CompareGOTerms(std::string go_ffindex,std::string go_ffindex_indexfile,std::string protid_go_ffindex,std::string protid_go_ffindex_indexfile){

    Debug(Debug::INFO) << "Opening sequence database...\n";
    go_ffindex_reader = new DBReader(go_ffindex.c_str(), go_ffindex_indexfile.c_str());
    go_ffindex_reader->open(DBReader::NOSORT);
    Debug(Debug::INFO) << "Opening sequence database...\n";
    protid_go_ffindex_reader = new DBReader(protid_go_ffindex.c_str(), protid_go_ffindex_indexfile.c_str());
    protid_go_ffindex_reader->open(DBReader::NOSORT);

}

CompareGOTerms::~CompareGOTerms() {
    go_ffindex_reader->~DBReader();
    protid_go_ffindex_reader->~DBReader();
}

void CompareGOTerms::init() {

}