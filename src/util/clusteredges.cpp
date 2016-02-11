#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

#include <fstream>

int clusteredges(int argc, const char *argv[]) {
    std::string usage("Extracts all child to parent edges from clustering N to clustering N+1.\n");
    usage.append("Written by Milot Mirdita <milot@mirdita.de>.\n\n");
    usage.append("USAGE: <edgeListOut> <clusteringSteps...> \n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.onlyverbosity, 1, true, true);
    Debug::setDebugLevel(par.verbosity);

    size_t clusterings = static_cast<size_t >(argc - 1);

    if (clusterings == 0) {
        Debug(Debug::ERROR) << "Need at least one cluster result to extract edges\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> *readers[clusterings];
    for (int i = 0; i < clusterings; ++i) {
        std::pair<std::string, std::string> name = Util::databaseNames(argv[i + 1]);
        readers[i] = new DBReader<unsigned int>(name.first.c_str(), name.second.c_str());
    }

    if (FileUtil::fileExists(par.db1.c_str())) {
        errno = EEXIST;
        perror(par.db1.c_str());
        EXIT(EXIT_FAILURE);
    }

    std::ofstream edgeList;
    edgeList.open(par.db1);

    if (edgeList.fail()) {
        Debug(Debug::ERROR) << "Could not open " << par.db1 << " for writing.";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Start writing file to " << par.db1 << "\n";
    for (size_t i = 0; i < clusterings; ++i) {
        DBReader<unsigned int> *reader = readers[i];
        reader->open(DBReader<unsigned int>::NOSORT);

        for (size_t j = 0; j < reader->getSize(); ++j) {
            unsigned int parentKey = reader->getDbKey(j);

            // add the implicit root explicitly
            if (i == 0) {
                edgeList << parentKey << "\tNULL\n";
            }

            std::istringstream ss(reader->getData(j));

            std::string line;
            std::getline(ss, line);
            while (std::getline(ss, line)) {
                if (line.length() < 1)
                    continue;
                unsigned int childKey = static_cast<unsigned int>(strtoul(line.c_str(), NULL, 10));

                edgeList << childKey << "\t" << parentKey << "\n";
            }
        }

        reader->close();
    }

    for (int i = 0; i < clusterings; ++i) {
        delete readers[i];
    }

    edgeList.close();

    return EXIT_SUCCESS;
}
