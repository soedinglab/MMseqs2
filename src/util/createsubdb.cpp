#include "Parameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#include <climits>

int createsubdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 3);

    FILE *orderFile = NULL;
    if (FileUtil::fileExists(par.db1Index.c_str())) {
        orderFile = fopen(par.db1Index.c_str(), "r");
    } else {
        if(FileUtil::fileExists(par.db1.c_str())){
            orderFile = fopen(par.db1.c_str(), "r");
        }else{
            Debug(Debug::ERROR) << "File " << par.db1 << " does not exist.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);
    const bool isCompressed = reader.isCompressed();

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, 0, reader.getDbtype() & ~(1 << 31));
    writer.open();

    char *line = (char*)malloc(1024);
    size_t len = 0;
    char dbKey[256];
    while (getline(&line, &len, orderFile) != -1) {
        Util::parseKey(line, dbKey);
        const unsigned int key = Util::fast_atoi<unsigned int>(dbKey);
        const size_t id = reader.getId(key);
        if (id >= UINT_MAX) {
            Debug(Debug::WARNING) << "Key " << line << " not found in database\n";
            continue;
        }
        char* data = reader.getDataUncompressed(id);
        size_t originalLength = reader.getSeqLens(id);
        size_t entryLength = std::max(originalLength, static_cast<size_t>(1)) - 1;
        if (isCompressed) {
            // copy also the null byte since it contains the information if compressed or not
            entryLength = *(reinterpret_cast<unsigned int*>(data)) + sizeof(unsigned int) + 1;
            writer.writeData(data, entryLength, key, 0, false, false);
        }else{
            writer.writeData(data, entryLength, key, 0, true, false);
        }
        // do not write null byte since
        writer.writeIndexEntry(key, writer.getStart(0), originalLength, 0);
    }
    writer.close();
    DBWriter::writeDbtypeFile(par.db3.c_str(), reader.getDbtype(), isCompressed);
    free(line);
    reader.close();
    fclose(orderFile);

    return EXIT_SUCCESS;
}
