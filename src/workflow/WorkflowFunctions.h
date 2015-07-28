#ifndef WORKFLOW_FUNCTIONS_H
#define WORKFLOW_FUNCTIONS_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <list>

#include "Prefiltering.h"
#include "Alignment.h"
#include "Clustering.h"
#include "Parameters.h"


extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}


std::string runPrefilter(std::string inDBData, std::string inDBWorkingIndex, std::string targetDBData,
                         std::string targetDBIndex, std::string tmpDir,
                         Parameters par, int step_num, std::list<std::string>* tmpFiles);

std::string runAlignment(std::string inDBData, std::string inDBWorkingIndex, std::string targetDBData,
                         std::string perfResult, std::string targetDBIndex, std::string tmpDir,
                         Parameters par, int step_num, std::list<std::string>* tmpFiles);

std::string runStep(std::string inDBData, std::string inDBWorkingIndex, std::string targetDBData, std::string targetDBIndex, std::string tmpDir,
        Parameters par,
        int step_num, int restart, bool search, std::list<std::string>* tmpFiles);

void copy(std::string inFile, std::string outFile);

float getZscoreForSensitivity (float sensitivity);

void deleteTmpFiles(std::list<std::string>* tmpFiles);

#endif
