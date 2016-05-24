
#include "DBConcat.h"
#include "DBReader.h"


DBConcat::DBConcat(const char* dataFileNameA, const char* indexFileNameA,const char* dataFileNameB, const char* indexFileNameB,const char* dataFileNameC, const char* indexFileNameC, int threads, int dataMode, bool preserveKeysA)
			: preserveKeysA(preserveKeysA),DBReader(std::strcmp(dataFileNameA,dataFileNameB) == 0 ? dataFileNameA : dataFileNameC,std::strcmp(indexFileNameA,indexFileNameB) == 0 ? indexFileNameA : indexFileNameC,dataMode)
{
	lastKeyA = 0;
	
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
	
	if (!sameDatabase)
	{
	
		dbA->open(DBReader<unsigned int>::NOSORT);
		dbB->open(DBReader<unsigned int>::NOSORT);
		
		indexSizeA = dbA->getSize();
		indexSizeB = dbB->getSize();
		
		// keys paris are like : (key,i) where key is the ith key in the ffindex
		keysA = new std::pair<unsigned int, unsigned int> [indexSizeA];
		keysB = new std::pair<unsigned int, unsigned int> [indexSizeB];
		
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
			unsigned int newKey;
			
			if (preserveKeysA)
				newKey = dbA->getDbKey(id);
			else
				newKey = id;
				
			concatWriter->write(data, dbA->getSeqLens(id) -1 , SSTR(newKey).c_str(), thread_idx);
			keysA[id] = std::make_pair(dbA->getDbKey(id),newKey);// need to store the index, because it'll be sorted out by keys later
			
			lastKeyA = std::max(lastKeyA,newKey);
	
		}

	lastKeyA++;
	
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
			concatWriter->write(data, dbB->getSeqLens(id) -1, SSTR(id+lastKeyA).c_str(), thread_idx);
			keysB[id] = std::make_pair(dbB->getDbKey(id),id+lastKeyA);// need to store the index, because it'll be sorted out by keys later
			
		}
	}
	
		//sort by key
		std::stable_sort(keysA, keysA  + indexSizeA, compareFirstEntry());
		std::stable_sort(keysB, keysB  + indexSizeB, compareFirstEntry());
		
		concatWriter->close();
		delete concatWriter;
		dbA->close();
		delete dbA;
		dbB->close();
		delete dbB;
	}
}


unsigned int DBConcat::dbAKeyMap(unsigned int key) {
	if (sameDatabase)
		return key;
	
	std::pair<unsigned  int, unsigned  int>* originalMap = std::upper_bound(keysA, keysA + indexSizeA, key, compareKeyToFirstEntry());
	return (*originalMap).second;
}

unsigned int DBConcat::dbBKeyMap(unsigned int key) {
	if (sameDatabase)
		return key;
		
	std::pair<unsigned  int, unsigned  int>*originalMap = std::upper_bound(keysB, keysB + indexSizeB, key, compareKeyToFirstEntry());
	return (*originalMap).second;
}



DBConcat::~DBConcat() {
	if (!sameDatabase)
	{
		delete keysA;
		delete keysB;
	}
}



int dbconcat(int argc, const char **argv) {
    std::string usage("Concatenates two ffindex DB.\n");
    usage.append("USAGE: <DB1> <DB2> <outDB>\n");
    usage.append("\nDesigned and implemented by Clovis Galiez <clovis.galiez@mpibpc.mpg.de>\n");


    Parameters par;
    par.parseParameters(argc, argv, usage, par.dbconcat, 3);


    
    struct timeval start, end;
    gettimeofday(&start, NULL);
    
	//if (par.compressMSA)
	
	DBConcat outDB(par.db1.c_str(), par.db1Index.c_str(),par.db2.c_str(), par.db2Index.c_str(),par.db3.c_str(),par.db3Index.c_str(),1,DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX,true);
	outDB.concat();
	
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for concatenating DBs: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    return 0;

}
