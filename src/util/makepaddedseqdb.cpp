#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

int makepaddedseqdb(int argc, const char **argv, const Command &command) {
    const char AA_TO_20[256] = {
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20,  0, 20,  4,  3,  6, 13,  7,  8,  9, 20, 11, 10, 12,  2, 20,
        14,  5,  1, 15, 16, 20, 19, 17, 20, 18, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20
    };

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
        for (size_t i = 0; i < seqLen; i++) {
            result.append(1, AA_TO_20[(unsigned char)data[i]]);
        }
        result.append(sequencepadding, static_cast<char>(20));
        writer.writeData(result.c_str(), result.size(), key, 0, false, false);
        writer.writeIndexEntry(key, writer.getStart(0), seqLen, 0);
        result.clear();
    }
    writer.close(true, false);
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_ANCILLARY);

    dbr.close();
    return EXIT_SUCCESS;
}