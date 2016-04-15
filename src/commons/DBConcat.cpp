
#include "DBConcat.h"
#include "DBReader.h"


DBConcat::DBConcat(const char* dataFileNameA, const char* indexFileNameA,const char* dataFileNameB, const char* indexFileNameB,const char* dataFileNameC, const char* indexFileNameC, int threads, int dataMode)
			: DBReader(std::strcmp(dataFileNameA,dataFileNameB) == 0 ? dataFileNameA : dataFileNameC,std::strcmp(indexFileNameA,indexFileNameB) == 0 ? indexFileNameA : indexFileNameC,dataMode)
{
	if (std::strcmp(dataFileNameA,dataFileNameB) == 0)
		sameDatabase = true;
	else
		sameDatabase = false;
		
		
	dbA = new DBReader<unsigned int>(dataFileNameA, indexFileNameA);
	dbB = new DBReader<unsigned int>(dataFileNameB, indexFileNameB);
	
#ifdef OPENMP
    omp_set_num_threads(threads);
#endif

	concatWriter = new DBWriter(dataFileNameC, indexFileNameC, threads, DBWriter::BINARY_MODE);

}

 // If dbA != dbB, then Concatenate dbA and dbB in concatWriter ("dataFileNameC")
 // and "this" will be a reader on "dataFileNameC" after calling open()
 // otherwise, do nothing and "this"  will be a reader on "dataFileNameA"
void DBConcat::concat(){
	
	
	unsigned int currentKey = 0;
	
	if (!sameDatabase)
	{

		dbA->open(DBReader<unsigned int>::NOSORT);
		dbB->open(DBReader<unsigned int>::NOSORT);
		
		indexSizeA = dbA->getSize();
		indexSizeB = dbB->getSize();
		
		keysA = new unsigned int[indexSizeA];
		keysB = new unsigned int[indexSizeB];
		
		
		concatWriter->open();
		
	
	#pragma omp parallel
	{
		
		
		
	#pragma omp for schedule(static)
		for (size_t id = 0; id < indexSizeA;id++)
		{
			char *data;
			
			Log::printProgress(id);
			int thread_idx = 0;
	#ifdef OPENMP
			thread_idx = omp_get_thread_num();
	
	#endif
			data = dbA->getData(id);
			concatWriter->write(data, dbA->getSeqLens(id) -1 , SSTR(id).c_str(), thread_idx);
			keysA[id] = dbA->getDbKey(id);
			
	
			currentKey++;
		}


	#pragma omp for schedule(static)
		for (size_t id = 0; id < indexSizeB;id++)
		{
			char *data;
			
			Log::printProgress(id);
			int thread_idx = 0;
	#ifdef OPENMP
			thread_idx = omp_get_thread_num();
	#endif
			data = dbB->getData(id);
			concatWriter->write(data, dbB->getSeqLens(id) -1, SSTR(id+indexSizeA).c_str(), thread_idx);
			keysB[id] = dbB->getDbKey(id);
			
			currentKey++;
		}
	}
		concatWriter->close();
		dbA->close();
		dbB->close();
		
	}
}


unsigned int DBConcat::dbAKeyMap(unsigned int key) {
	if (sameDatabase)
		return key;
	
	for (unsigned int i = 0;i < indexSizeA; i++) {
		if (keysA[i] == key)
			return i;
	}
	Debug(Debug::ERROR) << "Cannot find in dbA the key " << key << " !\n";
	EXIT(EXIT_FAILURE);
}

unsigned int DBConcat::dbBKeyMap(unsigned int key) {
	if (sameDatabase)
		return key;
		
	for (unsigned int i = 0;i < indexSizeB; i++) {
		if (keysB[i] == key)
			return i + indexSizeA;
	}
	Debug(Debug::ERROR) << "Cannot find in dbB the key " << key << " !\n";
	EXIT(EXIT_FAILURE);
}


