#ifndef RESULT2PROFILE_H
#define RESULT2PROFILE_H


#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"
#include "DBConcat.h"


#define STAT_LINECOUNT_STR "linecount"
#define STAT_MEAN_STR "mean"
#define STAT_DOOLITTLE_STR "doolittle"
#define STAT_CHARGES_STR "charges"
#define STAT_SEQLEN_STR "seqlen"

class statsComputer {
public:
	statsComputer(Parameters &par);
	~statsComputer();
	int run();
private:
    const size_t LINE_BUFFER_SIZE = 1000000;
    
    enum {  STAT_LINECOUNT,
            STAT_MEAN,
            STAT_DOOLITTLE,
            STAT_CHARGES,
            STAT_SEQLEN,
            STAT_UNKNOWN
    };
	
    Parameters *par;
    int stat;

    DBReader<unsigned int> *queryReader;
    DBReader<unsigned int> *queryHeaderReader;
    DBReader<unsigned int> *targetReader;
    DBReader<unsigned int> *targetHeaderReader;
    DBReader<unsigned int> *resultReader;

    DBWriter *statWriter;
    
    float pH;
    std::unordered_map<char,float> doolittleValues;
    std::unordered_map<char,float> chargeValues;
    
    float averageValueOnAminoAcids(std::unordered_map<char,float> values,char* seq);
    int sequenceWise(float (statsComputer::*statFunction)(char*));
    
    int countNumberOfLines();
    int meanValue();
    float doolittle(char *);
    float charges(char *);
    float strlen(char *);
	
};


#endif