#ifndef FILTERDB_H
#define FILTERDB_H

// Written by Martin Steinegger & Clovis Galiez
//
// Filter a ffindex based on a RegEx or a filtering file.
//

#include <cstddef>
#include <utility>
#include <string>
#include <vector>
#include <regex.h>

#define REGEX_FILTERING 0
#define FILE_FILTERING 1
#define FILE_MAPPING 2
#define GET_FIRST_LINES 3
#define NUMERIC_COMPARISON 4

#define GREATER_OR_EQUAL "ge"
#define LOWER_OR_EQUAL "le"
#define EQUAL "e"

class ffindexFilter {
public:
	// Constructor for RegEx Filtering
	ffindexFilter(std::string inDB,
				  std::string outDB,
				  int threads,
				  size_t column,
				  std::string regexStr,
				  bool trimToOneColumn = false);

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

	// Constructor for extract first line
	ffindexFilter(std::string inDB,
				  std::string outDB,
				  int threads,
				  size_t column,
                 int numberOfLines);
                 
	// Constructor for numeric comparison
	ffindexFilter(std::string inDB,
				  std::string outDB,
				  int threads,
				  size_t column,
                 float compValue,
                 std::string compOperator);
				  
	~ffindexFilter();
	
	int runFilter();

private:
	std::string inDB;
	std::string outDB;
    std::string filterFile;

	int threads;
	size_t column;
    std::string regexStr;
    bool trimToOneColumn;
    // positiveFilter = true => outDB = inDB \intersect filter ; othw : outDB = inDB - filter
    bool positiveFiltering;
    int numberOfLines;
    int mode;
    float compValue;
	std::string compOperator;

    DBWriter* dbw;
	DBReader<unsigned int>* dataDb;
	
	regex_t regex;
	std::vector<std::string> filter;

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
		bool operator() (const std::pair<std::string,std::string>& lhs,const std::string& rhs) const{
			return (lhs.first.compare(rhs)<0);
		}
	};

};

#endif
