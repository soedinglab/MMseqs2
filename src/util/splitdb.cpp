#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include <sstream>

int splitdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    if(par.split < 1) {
        Debug(Debug::ERROR) << "Cannot split databases into 0 or negative chunks.";
        EXIT(EXIT_FAILURE);
    }

    DBReader<std::string> dbr(par.db1.c_str(), par.db1Index.c_str());
    dbr.open(DBReader<std::string>::NOSORT);

    unsigned int* sizes = dbr.getSeqLens();
    size_t size = dbr.getSize();

    if((size_t)par.split > size) {
        Debug(Debug::ERROR) << "Cannot split databases into more chunks than database contains.";
        EXIT(EXIT_FAILURE);
    }

    for (int split = 0; split < par.split; split++) {
        std::ostringstream outName;
        outName << par.db2 << "_" << split << "_" << par.split;
        std::string outString = outName.str();
        DBWriter writer(outString.c_str(), std::string(outName.str() + ".index").c_str());
        writer.open();

        size_t startIndex = 0;
        size_t domainSize = 0;

        if(par.splitAA) {
            Util::decomposeDomainByAminoAcid(dbr.getAminoAcidDBSize(), sizes, size, split, par.split, &startIndex,
                                             &domainSize);
        } else {
            Util::decomposeDomain(size, split, par.split, &startIndex, &domainSize);
        }

        for(size_t i = startIndex; i < (startIndex + domainSize); i++){
            std::string outerKey = dbr.getDbKey(i);
            char * data = dbr.getData(i);
            writer.writeData(data, sizes[i], (char *) outerKey.c_str());
        }
        writer.close();
    }

    dbr.close();
    return EXIT_SUCCESS;
}
