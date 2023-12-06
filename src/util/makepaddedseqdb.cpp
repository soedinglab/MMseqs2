#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

int makepaddedseqdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str(), 1,
                               DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::SORT_BY_LENGTH);
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), 1, false, dbr.getDbtype());
    writer.open();
    std::string result;
    const int ALIGN = 4;
    for (long id = dbr.getSize() - 1; id >= 0; id--) {
        unsigned int key = dbr.getDbKey(id);
        char *data = dbr.getData(id, 0);
        size_t seqLen = dbr.getSeqLen(id);
        const size_t sequencepadding = (seqLen % ALIGN == 0) ? 0 : ALIGN - seqLen % ALIGN;
        result.append(data, seqLen);
        result.append(sequencepadding, ' ');
        writer.writeData(data, seqLen + sequencepadding, key, 0, false);
    }
    writer.close(true);
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_ANCILLARY);

    dbr.close();
    return EXIT_SUCCESS;
}