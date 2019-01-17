//
// Created by Martin Steinegger on 2019-01-17.
//

#include "Parameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#include <climits>

int view(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 1, false);
    std::vector<std::string> ids = Util::split(par.idList, ",");
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);
    char dbKey[256];
    for (size_t i = 0; i< ids.size(); i++) {
        strncpy(dbKey, ids[i].c_str(), ids[i].size());
        dbKey[ids[i].size()]='\0';
        const unsigned int key = Util::fast_atoi<unsigned int>(dbKey);
        const size_t id = reader.getId(key);
        if (id >= UINT_MAX) {
            Debug(Debug::WARNING) << "Key " << ids[i] << " not found in database\n";
            continue;
        }
        char* data = reader.getData(id, 0);
        Debug(Debug::INFO) << data;
    }
    reader.close();
    EXIT(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}
