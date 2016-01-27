#ifndef WORKFLOW_FUNCTIONS_H
#define WORKFLOW_FUNCTIONS_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <list>

void copy(std::string inFile, std::string outFile);

void deleteTmpFiles(std::list<std::string>* tmpFiles);

#endif
