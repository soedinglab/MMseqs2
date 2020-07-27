#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

int splitdb(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    if (par.split < 1) {
        Debug(Debug::ERROR) << "Cannot split databases into 0 or negative chunks.";
        EXIT(EXIT_FAILURE);
    }

    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::NOSORT);


    if ((size_t) par.split > dbr.getSize()) {
        Debug(Debug::ERROR) << "Cannot split databases into more chunks than database contains.";
        EXIT(EXIT_FAILURE);
    }

    for (int split = 0; split < par.split; split++) {
        std::string outDb = par.db2 + "_" + SSTR(split) + "_" + SSTR(par.split);
        DBWriter writer(outDb.c_str(), std::string(outDb + ".index").c_str(), 1, par.compressed, dbr.getDbtype());
        writer.open();

        size_t startIndex = 0;
        size_t domainSize = 0;
        if (par.splitAA) {
            dbr.decomposeDomainByAminoAcid(split, par.split, &startIndex, &domainSize);
        } else {
            Util::decomposeDomain(dbr.getSize(), split, par.split, &startIndex, &domainSize);
        }

        for (size_t i = startIndex; i < (startIndex + domainSize); i++) {
            unsigned int outerKey = dbr.getDbKey(i);
            char *data = dbr.getData(i, 0);
            writer.writeData(data, dbr.getEntryLen(i) - 1, outerKey);
        }
        writer.close();
    }

    dbr.close();
    return EXIT_SUCCESS;
}
