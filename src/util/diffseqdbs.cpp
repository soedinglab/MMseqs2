// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>
#include <algorithm>

#include "Parameters.h"

#include "DBReader.h"
#include "DBConcat.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif



struct compareFirstEntry {
	bool operator() (const std::pair<std::string, unsigned  int>& lhs, const std::pair<std::string, unsigned  int>& rhs) const{
		return (lhs.first.compare(rhs.first)<=0);
	}
};

struct compareKeyToFirstEntry {
	bool operator() (const std::string& lhs, const std::pair<std::string, unsigned  int>& rhs) const{
		return (lhs.compare(rhs.first)<=0);
	}
};

int diffseqdbs(int argc, const char **argv, const Command& command) {
	Parameters& par = Parameters::getInstance();
	par.parseParameters(argc, argv, command, 5);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    std::string headerDBold(par.db1);
    headerDBold.append("_h");
    std::string headerDBoldIndex(par.db1);
    headerDBoldIndex.append("_h.index");
	
    std::string headerDBnew(par.db2);
    headerDBnew.append("_h");
	std::string headerDBnewIndex(par.db2);
    headerDBnewIndex.append("_h.index");
	
	
    DBReader<unsigned int> *DBoldReader = new DBReader<unsigned int>(headerDBold.c_str(), headerDBoldIndex.c_str());
    DBoldReader->open(DBReader<unsigned int>::NOSORT);

	
    DBReader<unsigned int> *DBnewReader = new DBReader<unsigned int>(headerDBnew.c_str(), headerDBnewIndex.c_str());
    DBnewReader->open(DBReader<unsigned int>::NOSORT);
	
	std::ofstream removedSeqDBWriter,keptSeqDBWriter,newSeqDBWriter;	
	removedSeqDBWriter.open(par.db3);
	keptSeqDBWriter.open(par.db4);
	newSeqDBWriter.open(par.db5);

	size_t indexSizeOld = DBoldReader->getSize();
	size_t indexSizeNew = DBnewReader->getSize();

	// keys pairs are like : (headerID,key) where key is the ffindex key corresponding to the header
	std::pair<std::string, unsigned int> *keysOld = new std::pair<std::string, unsigned int> [indexSizeOld];
	std::pair<std::string, unsigned int> *keysNew = new std::pair<std::string, unsigned int> [indexSizeNew];
	bool *checkedNew = new bool[indexSizeNew]; // store if the sequence has been seen in the old DB

	// Fill up the hash tables for the old and new DB
#pragma omp for schedule(static)
	for (size_t id = 0; id < indexSizeOld; id++) {
		keysOld[id] = std::make_pair(std::string(DBoldReader->getData(id)),DBoldReader->getDbKey(id));
	}
#pragma omp for schedule(static)
	for (size_t id = 0; id < indexSizeNew; id++) {
		keysNew[id] = std::make_pair(std::string(DBnewReader->getData(id)),DBnewReader->getDbKey(id));
		checkedNew[id] = false;
	}

	//sort by header
	std::stable_sort(keysOld, keysOld  + indexSizeOld, compareFirstEntry());
	std::stable_sort(keysNew, keysNew  + indexSizeNew, compareFirstEntry());
	
	#pragma omp for schedule(static)
	for (size_t id = 0; id < indexSizeOld; id++) {
		
		std::string keyToSearch = std::string(keysOld[id].first);
		std::pair<std::string, unsigned  int>* mappedKey = std::upper_bound(keysNew, keysNew + indexSizeNew, keyToSearch, compareKeyToFirstEntry());


		if (mappedKey != keysNew + indexSizeNew && keyToSearch == mappedKey->first)
		{
			// Found
			size_t indexInNewDB = (mappedKey - keysNew);// / sizeof(std::pair<std::string, unsigned  int>);	
			//std::cout << indexInNewDB <<std::endl;
			checkedNew[indexInNewDB] = true;
			keptSeqDBWriter << keysOld[id].second << "\t" << keysNew[indexInNewDB].second << std::endl;
			
		} else {
			// not found
			removedSeqDBWriter << keysOld[id].second << std::endl;
		}
	}
	
	#pragma omp for schedule(static)
	for (size_t id = 0; id < indexSizeNew; id++) {
		if (!checkedNew[id])
			newSeqDBWriter << keysNew[id].second << std::endl;
	}
	


	DBoldReader->close();
	DBnewReader->close();
	
	delete DBoldReader;
	delete DBnewReader;
	
	delete[] keysOld;
	delete[] keysNew;
	
	delete[] checkedNew;
	
	removedSeqDBWriter.close();
	keptSeqDBWriter.close();
	newSeqDBWriter.close();
	
	return 0;
}

