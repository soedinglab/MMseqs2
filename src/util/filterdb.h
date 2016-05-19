#ifndef FILTERDB_H
#define FILTERDB_H

// Written by Martin Steinegger & Clovis Galiez
//
// Filter a ffindex based on a RegEx or a filtering file.
//

#include <cstddef>
#include <utility>
#include <string>
#include <fstream>
#include <iostream>
#include <regex.h>
#include <vector>
#include <algorithm> 
#include <algorithm> 
#include <sstream>


#define REGEX_FILTERING 0
#define FILE_FILTERING 1
#define FILE_MAPPING 2


class ffindexFilter {

public:
	// Constructor for RegEx Filtering
	ffindexFilter(std::string inDB,
					std::string outDB,
                  int threads,
					size_t column,
                  std::string regexStr);
	// Constructor for File based Filtering			  
	ffindexFilter(std::string inDB,
					std::string outDB,
					std::string filterFile,
                  int threads,
					size_t column,
					bool positiveFiltering);
	// Constructor for ID mapping
	ffindexFilter(std::string inDB,
					std::string outDB,
					std::string filterFile,
                  int threads,
					size_t column);
				  
	~ffindexFilter();
	
	int runFilter();
private:
	int mode;
	int threads;
	size_t column;
	
	std::string inDB,outDB;
	DBWriter* dbw;
	DBReader<unsigned int>* dataDb;
	
	std::string regexStr;
	regex_t regex;
	std::string filterFile;
	std::vector<std::string> filter;
	bool positiveFiltering; // positiveFilter = true => outDB = inDB \intersect filter ; othw : outDB = inDB - filter
	
	
	std::vector<std::pair<std::string,std::string>> mapping;
	
	int initFiles();
	
	

	struct compareString {
		bool operator() (const std::string& lhs, const std::string& rhs) const{
			return (lhs.compare(rhs)<=0);
		}
	};


	struct compareFirstString {
		bool operator() (const std::pair<std::string, std::string>& lhs, const std::pair<std::string,std::string>& rhs) const{
			return (lhs.first.compare(rhs.first)<=0);
		}
	};

	struct compareToFirstString {
		bool operator() (const std::string& lhs, const std::pair<std::string,std::string>& rhs) const{
			return (lhs.compare(rhs.first)<=0);
		}
	};

};

#endif
