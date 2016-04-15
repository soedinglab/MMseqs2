#ifndef DBCONCAT_H
#define DBCONCAT_H

#include "DBWriter.h"
#include "DBReader.h"

#include "Log.h"
#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

class DBConcat : public DBReader<unsigned int> {

public:
	DBConcat(const char* dataFileNameA, const char* indexFileNameA,
				const char* dataFileNameB, const char* indexFileNameB,
				const char* dataFileNameC, const char* indexFileNameC, int threads, int dataMode = USE_DATA|USE_INDEX);
	void concat();
	
	unsigned int dbAKeyMap(unsigned int);
	unsigned int dbBKeyMap(unsigned int);
	
private:
	DBWriter *concatWriter;
	
	DBReader<unsigned int> *dbA;
	DBReader<unsigned int> *dbB;
			
	size_t indexSizeA;
	size_t indexSizeB;
	
	unsigned int *keysA,*keysB;
	
	int threadsNb;
	
	bool sameDatabase;

};


#endif
