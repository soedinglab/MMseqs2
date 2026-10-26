#include "Parameters.h"
#include "IndexReader.h"
#include "Debug.h"
#include "Util.h"

int view(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, 0, 0);
    std::vector<std::string> ids = Util::split(par.idList, ",");
    int indexSrcType = IndexReader::SEQUENCES;
    switch (par.idxEntryType) {
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
    const bool lookupMode = par.dbIdMode == Parameters::ID_MODE_LOOKUP;
    int dbMode = DBReader<DBKeyType>::USE_INDEX|DBReader<DBKeyType>::USE_DATA;
    if (lookupMode) {
        dbMode |= DBReader<DBKeyType>::USE_LOOKUP_REV;
    }
    IndexReader reader(par.db1, par.threads, indexSrcType, false, dbMode);
    for (size_t i = 0; i < ids.size(); ++i) {
        DBKeyType key;
        std::string& ref = ids[i];
        if (lookupMode) {
            size_t lookupId = reader.sequenceReader->getLookupIdByAccession(ref);
            if (lookupId == SIZE_MAX) {
                Debug(Debug::WARNING) << "Could not find " << ref << " in lookup\n";
                continue;
            }
            key = reader.sequenceReader->getLookupKey(lookupId);
        } else {
            key = Util::fast_atoi<DBKeyType>(ref.c_str());
        }

        const size_t id = reader.sequenceReader->getId(key);
        if (id == DB_ENTRY_NOT_FOUND) {
            Debug(Debug::ERROR) << "Key " << ids[i] << " not found in database\n";
            continue;
        }
        char* data = reader.sequenceReader->getData(id, 0);
        size_t size = reader.sequenceReader->getEntryLen(id) - 1;
        fwrite(data, sizeof(char), size, stdout);
    }
    EXIT(EXIT_SUCCESS);
}
