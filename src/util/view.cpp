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
#include <IndexReader.h>

int view(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.verbosity = 1;
    par.parseParameters(argc, argv, command, true, 0, 0);
    std::vector<std::string> ids = Util::split(par.idList, ",");
    int indexSrcType = IndexReader::SEQUENCES;
    switch(par.idxEntryType){
        case 0:
            indexSrcType = IndexReader::SEQUENCES;
            break;
        case 1:
            indexSrcType = IndexReader::SRC_SEQUENCES;
            break;
        case 2:
            indexSrcType = IndexReader::HEADERS;
            break;
        case 3:
            indexSrcType = IndexReader::SRC_HEADERS;
            break;
    }
    IndexReader reader(par.db1, par.threads, indexSrcType, 0);
    char dbKey[256];
    for (size_t i = 0; i< ids.size(); i++) {
        strncpy(dbKey, ids[i].c_str(), ids[i].size());
        dbKey[ids[i].size()]='\0';
        const unsigned int key = Util::fast_atoi<unsigned int>(dbKey);
        const size_t id = reader.sequenceReader->getId(key);
        if (id >= UINT_MAX) {
            Debug(Debug::WARNING) << "Key " << ids[i] << " not found in database\n";
            continue;
        }
        char* data = reader.sequenceReader->getData(id, 0);
        std::cout << data;
    }
    EXIT(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}
