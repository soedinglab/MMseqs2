//
// Created by mad on 5/9/15.
//

#ifndef MMSEQS_MERGECLUSTER_H
#define MMSEQS_MERGECLUSTER_H

#include <string>
#include <list>

void mergeClusteringResults(std::string seqDB, std::string outDB, std::list<std::string> cluSteps);


int mergecluster(int argn,const char **argv);

#endif //MMSEQS_MERGECLUSTER_H
