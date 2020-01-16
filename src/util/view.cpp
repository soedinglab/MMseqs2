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
    IndexReader reader(par.db1, par.threads, indexSrcType, 0);
    for (size_t i = 0; i < ids.size(); ++i) {
        const unsigned int key = Util::fast_atoi<unsigned int>(ids[i].c_str());
        const size_t id = reader.sequenceReader->getId(key);
        if (id >= UINT_MAX) {
            Debug(Debug::ERROR) << "Key " << ids[i] << " not found in database\n";
            continue;
        }
        char* data = reader.sequenceReader->getData(id, 0);
        size_t size = reader.sequenceReader->getEntryLen(id) - 1;
        fwrite(data, sizeof(char), size, stdout);
    }
    EXIT(EXIT_SUCCESS);
}
