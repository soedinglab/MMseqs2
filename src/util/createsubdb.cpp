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

    FILE *orderFile =  fopen(par.db1.c_str(), "r");
    if (FileUtil::fileExists(par.db1Index.c_str())) {
        orderFile = fopen(par.db1Index.c_str(), "r");
    }

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str());
    reader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str());
    writer.open();

    Debug(Debug::INFO) << "Start writing to file " << par.db3 << "\n";
    char * line = new char[65536];
    char dbKey[255 + 1];
    size_t len = 0;
    while (getline(&line, &len, orderFile) != -1) {
        Util::parseKey(line, dbKey);
        const unsigned int key = Util::fast_atoi<unsigned int>(dbKey);
        size_t id = reader.getId(key);
        if(id >= UINT_MAX) {
            Debug(Debug::WARNING) << "Key " << line << " not found in database\n";
            continue;
        }

        const char* data = reader.getData(id);
        // discard null byte
        size_t length = reader.getSeqLens(id) - 1;
        writer.writeData(data, length, key);
    }

    if(FileUtil::fileExists((par.db2 + ".dbtype").c_str())){
        FileUtil::copyFile((par.db2 + ".dbtype").c_str(), (par.db3 + ".dbtype").c_str());
    }
    writer.close();

    delete[] line;
    reader.close();
    fclose(orderFile);

    return EXIT_SUCCESS;
}
