#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include <vector>
#include <sstream>

int splitffindex (int argc, const char * argv[])
{
    int err = EXIT_SUCCESS;

    std::string usage;
    usage.append("\nSplits a ffindex database into multiple ffindex databases.\n");
    usage.append("Written by Milot Mirdita (milot@mirdita.de).\n\n");
    usage.append("USAGE: <ffindexInDB> <ffindexOutDB>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.splitffindex, 2);

    if(par.split < 1) {
        Debug(Debug::ERROR) << "Cannot split databases into 0 or negative chunks.";
        EXIT(EXIT_FAILURE);
    }

    DBReader dbr(par.db1.c_str(), par.db1Index.c_str());
    dbr.open(DBReader::NOSORT);

    unsigned int* sizes = dbr.getSeqLens();
    size_t size = dbr.getSize();

    if(par.split > size) {
        Debug(Debug::ERROR) << "Cannot split databases into more chunks than database contains.";
        EXIT(EXIT_FAILURE);
    }

    for (size_t split = 0; split < par.split; split++) {
        std::ostringstream outName;
        outName << par.db2 << "_" << split << "_" << par.split;
        DBWriter writer(outName.str().c_str(), std::string(outName.str() + ".index").c_str());
        writer.open();

        size_t startIndex = 0;
        size_t domainSize = 0;

        if(par.splitAA) {
            Util::decomposeDomainByAminoaAcid(dbr.getAminoAcidDBSize(), sizes, size, split, par.split, &startIndex, &domainSize);
        } else {
            Util::decompose_domain(size, split, par.split, &startIndex, &domainSize);
        }

        for(size_t i = startIndex; i < (startIndex + domainSize); i++){
            char * outerKey = dbr.getDbKey(i);
            char * data = dbr.getData(i);
            writer.write(data, sizes[i], outerKey);
        }
        writer.close();
    }

    dbr.close();
    return err;
}
