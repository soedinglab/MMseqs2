#include "WorkflowFunctions.h"
#include "Util.h"

void copy(std::string inFile, std::string outFile){
    // buffer
    std::ifstream inf(inFile.c_str(), std::ios::binary);
    std::ofstream outf(outFile.c_str(), std::ios::binary);
    if (inf.is_open()) {
        outf << inf.rdbuf();
        inf.close();
        outf.close();
    }
    else{
        std::cerr << "Could not open file " << inFile << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void deleteTmpFiles(std::list<std::string>* tmpFiles){
    for (std::list<std::string>::iterator it = tmpFiles->begin(); it != tmpFiles->end(); it++){
        std::cout << "Deleting " << *it << "\n";
        if( remove((*it).c_str()) != 0 )
            std::cerr << "Error deleting file " << *it << "\n";
    }
}

