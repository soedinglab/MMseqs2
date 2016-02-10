#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

#include <climits>
#include <sstream>
#include <fstream>

int clusteredges(int argc, const char * argv[])
{
    std::string usage("Extracts all child to parent edges from clustering N to clustering N+1.\n");
    usage.append("Written by Milot Mirdita <milot@mirdita.de>.\n\n");
    usage.append("USAGE: <edgeListOut> <clusteringSteps...> \n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.onlyverbosity, 1, true, true);
    Debug::setDebugLevel(par.verbosity);

    size_t clusterings = static_cast<size_t >(argc - 1);

    if(clusterings < 2) {
        Debug(Debug::ERROR) << "Need at least two cluster results to extract edges\n";
        EXIT(EXIT_FAILURE);
    }

    typedef DBReader<unsigned int> Reader;

    Reader* readers[clusterings];
    for(int i = 0; i < clusterings; ++i) {
        std::pair<std::string, std::string> name = Util::databaseNames(argv[i+1]);
        readers[i] = new Reader(name.first.c_str(), name.second.c_str());
    }

    if(FileUtil::fileExists(par.db1.c_str())) {
        errno = EEXIST;
        perror(par.db1.c_str());
        EXIT(EXIT_FAILURE);
    }

    std::ofstream edgeList;
    edgeList.open(par.db1);

    if(edgeList.fail()) {
        Debug(Debug::ERROR) << "Could not open " << par.db1 << " for writing.";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Start writing file to " << par.db1 << "\n";
    for(size_t i = 1; i < clusterings; ++i){
        Reader* parentReader = readers[i - 1];
        parentReader->open(Reader::NOSORT);
        Reader* childReader  = readers[i];
        childReader->open(Reader::NOSORT);

        for(size_t j = 0; j < parentReader->getSize(); ++j) {
            unsigned int parentKey = parentReader->getDbKey(j);

            if(i == 1) {
                // add the implicit root explicitly
                edgeList << parentKey << "\tNULL\n";
            }

            std::istringstream ss(parentReader->getData(j));

            std::string line;
            std::getline(ss, line);
            while(std::getline(ss, line)) {
                if(line.length() < 1)
                    continue;
                unsigned int childKey = static_cast<unsigned int>(strtoul(line.c_str(), NULL, 10));
                //char* childData = childReader.getDataByDBKey(childKey);

                edgeList << childKey << "\t" << parentKey << "\n";
            }
        }

        childReader->close();
        parentReader->close();
    }

    for(int i = 0; i < clusterings; ++i) {
        delete readers[i];
    }

    edgeList.close();

    return EXIT_SUCCESS;
}
