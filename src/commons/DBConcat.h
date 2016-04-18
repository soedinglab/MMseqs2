#ifndef DBCONCAT_H
#define DBCONCAT_H

#include "DBWriter.h"
#include "DBReader.h"

#include "Log.h"
#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"

#include <algorithm> 

#ifdef OPENMP
#include <omp.h>
#endif

class DBConcat : public DBReader<unsigned int> {

public:
	DBConcat(const char* dataFileNameA, const char* indexFileNameA,
				const char* dataFileNameB, const char* indexFileNameB,
				const char* dataFileNameC, const char* indexFileNameC, int threads, int dataMode = USE_DATA|USE_INDEX);
	~DBConcat();
	
	void concat();
	
	unsigned int dbAKeyMap(unsigned int);
	unsigned int dbBKeyMap(unsigned int);
	
private:
	DBWriter *concatWriter;
	
	DBReader<unsigned int> *dbA;
	DBReader<unsigned int> *dbB;
			
	size_t indexSizeA;
	size_t indexSizeB;
	
	std::pair<unsigned int, unsigned int> *keysA,*keysB;
	
	int threadsNb;
	
	bool sameDatabase;

    struct compareFirstEntry {
        bool operator() (const std::pair<unsigned  int, unsigned  int>& lhs, const std::pair<unsigned  int, unsigned  int>& rhs) const{
            return (lhs.first < rhs.first);
        }
	};
	
    struct compareKeyToFirstEntry {
        bool operator() (const unsigned  int& lhs, const std::pair<unsigned  int, unsigned  int>& rhs) const{
            return (lhs <= rhs.first);
        }
	};

};


#endif
