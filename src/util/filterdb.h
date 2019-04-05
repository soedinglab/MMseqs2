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
#include "ExpressionParser.h"

#define REGEX_FILTERING 0
#define FILE_FILTERING 1
#define FILE_MAPPING 2
#define GET_FIRST_LINES 3
#define NUMERIC_COMPARISON 4
#define SORT_ENTRIES 5
#define BEATS_FIRST 6
#define JOIN_DB 7
#define COMPUTE_POSITIONS 8
#define TRANSITIVE_REPLACE 9
#define EXPRESSION_FILTERING 10

#define GREATER_OR_EQUAL "ge"
#define LOWER_OR_EQUAL "le"
#define EQUAL "e"

#define INCREASING  1
#define DECREASING  2
#define SHUFFLE     3

class ffindexFilter {
public:


    ffindexFilter(Parameters &par); 
	~ffindexFilter();
	
	int runFilter();

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

    struct compareFirstInt {
        bool operator() (const std::pair<unsigned int, unsigned int>& lhs, const std::pair<unsigned int, unsigned int>& rhs) const{
            return (lhs.first < rhs.first);
        }
    };

	struct compareToFirstString {
		bool operator() (const std::pair<std::string,std::string>& lhs,const std::string& rhs) const{
			return (lhs.first.compare(rhs)<0);
		}
	};


    struct compareToFirstInt {
        bool operator() (const std::pair<unsigned int, unsigned int>& lhs,const unsigned int rhs) const{
            return (lhs.first <= rhs);
        }
    };

	static bool compareToFirstInt(const std::pair<unsigned int, unsigned int>& lhs, const std::pair<unsigned int, unsigned int>&  rhs){
		return (lhs.first <= rhs.first);
	}


private:
	std::string inDB;
	std::string outDB;
    std::string filterFile;

    int sortingMode;
    
	int threads;
	int compressed;

	size_t column;
	int columnToTake;
    std::string regexStr;
    bool trimToOneColumn;
    // positiveFilter = true => outDB = inDB \intersect filter ; othw : outDB = inDB - filter
    bool positiveFiltering;
    int numberOfLines;
    int mode;
	double compValue;
	std::string compOperator;
    bool shouldAddSelfMatch;
    ExpressionParser* parser;
    std::vector<int> bindableParserColumns;

    DBWriter* dbw;
	DBReader<unsigned int>* dataDb;
	DBReader<unsigned int>* joinDB;
    DBReader<unsigned int>* swapDB;
    DBReader<unsigned int>* clusterDB;
	
	regex_t regex;
	std::vector<std::string> filter;

	std::vector<std::pair<std::string,std::string>> mapping;
	
	int initFiles();
	

};

#endif
