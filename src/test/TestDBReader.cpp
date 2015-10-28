#include <iostream>
#include <list>
#include <algorithm>
#include <math.h>


#include "Clustering.h"
#include "SetElement.h"

#include "DBReader.h"
#include "DBWriter.h"


int main(int argc, char **argv)
{
    // DBReader test
    DBReader reader("dataLinear", "dataLinear.index");
    reader.open(0);
    std::cout << reader.getSize() << std::endl;
    for(size_t i = 0; i < reader.getSize(); i++){
        std::cout << reader.getSeqLens(i) << std::endl;
        std::cout << reader.getData(i) << std::endl;
    }
    reader.close();
    DBReader reader2("dataGap", "dataGap.index");
    reader2.open(0);
    std::cout << reader2.getSize() << std::endl;
    for(size_t i = 0; i < reader2.getSize(); i++){
        std::cout << reader2.getSeqLens(i) << std::endl;
        std::cout << reader2.getData(i) << std::endl;
    }
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey("2") << std::endl;
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey("6") << std::endl;
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey("1") << std::endl;
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey("111") << std::endl;
    std::cout << "Check getDataByDBKey: " << reader2.getDataByDBKey("12") << std::endl;

    std::cout << "Check getId: " << reader2.getId("2") << std::endl;
    std::cout << "Check getId: " << reader2.getId("6") << std::endl;
    std::cout << "Check getId: " << reader2.getId("1") << std::endl;
    std::cout << "Check not found getId: " << reader2.getId("8") << std::endl;
    std::cout << "Check length: " << (reader2.getSeqLens(reader2.getId("111")) == 13) << std::endl;
    std::cout << "Check length: " << (reader2.getSeqLens(reader2.getId("12")) == 10) << std::endl;



    reader2.close();

}
