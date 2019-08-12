#include "Parameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MemoryMapped.h"

#include <climits>

int createsubdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string& file = par.db1Index;
    if (FileUtil::fileExists(file.c_str()) == false) {
        file = par.db1;
        if (FileUtil::fileExists(file.c_str()) == false) {
            Debug(Debug::ERROR) << "File " << file << " does not exist.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    MemoryMapped order(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char* data = (char *) order.getData();

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);
    const bool isCompressed = reader.isCompressed();

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, 0, Parameters::DBTYPE_OMIT_FILE);
    writer.open();

    char dbKey[256];
    while (*data != '\0') {
        Util::parseKey(data, dbKey);
        data = Util::skipLine(data);

        const unsigned int key = Util::fast_atoi<unsigned int>(dbKey);
        const size_t id = reader.getId(key);
        if (id >= UINT_MAX) {
            Debug(Debug::WARNING) << "Key " << dbKey << " not found in database\n";
            continue;
        }
        if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
            writer.writeIndexEntry(key, reader.getOffset(id), reader.getSeqLens(id), 0);
        } else {
            char* data = reader.getDataUncompressed(id);
            size_t originalLength = reader.getSeqLens(id);
            size_t entryLength = std::max(originalLength, static_cast<size_t>(1)) - 1;

            if (isCompressed) {
                // copy also the null byte since it contains the information if compressed or not
                entryLength = *(reinterpret_cast<unsigned int *>(data)) + sizeof(unsigned int) + 1;
                writer.writeData(data, entryLength, key, 0, false, false);
            } else {
                writer.writeData(data, entryLength, key, 0, true, false);
            }
            // do not write null byte since
            writer.writeIndexEntry(key, writer.getStart(0), originalLength, 0);
        }
    }
    writer.close();
    if (par.subDbMode == Parameters::SUBDB_MODE_SOFT) {
        DBReader<unsigned int>::softLink(reader, par.db3);
    }
    DBWriter::writeDbtypeFile(par.db3.c_str(), reader.getDbtype(), isCompressed);

    reader.close();
    order.close();

    return EXIT_SUCCESS;
}
