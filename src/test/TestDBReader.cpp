#include <iostream>
#include <list>
#include <algorithm>
#include <cmath>

#include "Clustering.h"
#include "DBReader.h"
#include "DBWriter.h"

const char* binary_name = "test_dbreader";

int main (int, const char**) {
    // DBReader test
    DBReader<unsigned int> reader("dataLinear", "dataLinear.index", 1, 0);
    reader.open(0);
    reader.readMmapedDataInMemory();
    reader.printMagicNumber();
    std::cout << reader.getSize() << std::endl;
    for(size_t i = 0; i < reader.getSize(); i++){
        std::cout << reader.getSeqLen(i) << std::endl;
        std::cout << reader.getData(i, 0) << std::endl;
    }
    reader.close();
    DBReader<unsigned int> reader2("dataGap", "dataGap.index", 1, 0);
    reader2.open(0);
    std::cout << reader2.getSize() << std::endl;
    for(size_t i = 0; i < reader2.getSize(); i++){
        std::cout << reader2.getSeqLen(i) << std::endl;
        std::cout << reader2.getData(i, 0) << std::endl;
    }
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey(2,0) << std::endl;
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey(6,0) << std::endl;
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey(1,0) << std::endl;
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey(111,0) << std::endl;
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey(12,0) << std::endl;

    std::cout << "Check getId: " << reader2.getId(2) << std::endl;
    std::cout << "Check getId: " << reader2.getId(6) << std::endl;
    std::cout << "Check getId: " << reader2.getId(1) << std::endl;
    std::cout << "Check not found getId: " << reader2.getId(8) << std::endl;
    std::cout << "Check length: " << (reader2.getSeqLen(reader2.getId(111)) == 13) << std::endl;
    std::cout << "Check length: " << (reader2.getSeqLen(reader2.getId(12)) == 10) << std::endl;
    reader2.close();
    // test sort mode
    DBReader<unsigned int> reader3("dataGap", "dataGap.index", 1, 0);
    reader3.open(DBReader<unsigned int>::SORT_BY_LENGTH);
    for(size_t i = 0; i < reader3.getSize(); i++){
        size_t id =  reader3.getDbKey(i);
        std::cout << id <<  "\t" << reader3.getSeqLen(i) << "\t" << reader3.getData(i, 0) ;
    }
    std::cout << reader3.getId(111) << "\t" << reader3.getDataByDBKey(111,0);
    std::cout << reader3.getId(12) << "\t" <<  reader3.getDataByDBKey(12,0);
    std::cout << reader3.getId(1) << "\t" <<  reader3.getDataByDBKey(1,0);
    std::cout << reader3.getId(2) << "\t" <<  reader3.getDataByDBKey(2,0);
    std::cout << reader3.getId(3) << "\t" <<  reader3.getDataByDBKey(3,0);
    std::cout << reader3.getId(4) << "\t" <<  reader3.getDataByDBKey(4,0);
    std::cout << reader3.getId(5) << "\t" <<   reader3.getDataByDBKey(5,0);
}
