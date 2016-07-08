#ifndef RESULT2PROFILE_H
#define RESULT2PROFILE_H


#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Log.h"
#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"
#include "DBConcat.h"


#define STAT_LINECOUNT_STR "linecount"


class statsComputer {
public:
	statsComputer(Parameters &par);
	~statsComputer();
	int run();
private:

    enum {  STAT_LINECOUNT,
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
    
    int countNumberOfLines();
	
};


#endif