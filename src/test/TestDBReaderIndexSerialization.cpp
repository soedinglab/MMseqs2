#include "Debug.h"
#include "DBReader.h"

const char* binary_name = "test_dbreaderindexserialization";

int main (int, const char**) {
    DBReader<IdType> reader("", "/Users/mirdita/tmp/db.index", 1, DBReader<IdType>::USE_INDEX);
    reader.open(DBReader<IdType>::NOSORT);

    Debug(Debug::INFO) << reader.getSize() << " " << reader.getAminoAcidDBSize() << "\n";
    Debug(Debug::INFO) << reader.getIndex()[0].id  << " " << reader.getIndex()[0].offset  << " " << reader.getIndex()[0].length  << "\n";

    char* data = DBReader<IdType>::serialize(reader);
    DBReader<IdType>* newdbr = DBReader<IdType>::unserialize(data, 1);
    newdbr->open(DBReader<IdType>::NOSORT);

    Debug(Debug::INFO) << newdbr->getSize() << " " << newdbr->getAminoAcidDBSize() << "\n";
    Debug(Debug::INFO) << newdbr->getIndex()[0].id  << " " << newdbr->getIndex()[0].offset  << " " << newdbr->getIndex()[0].length << "\n";
    free(data);

    newdbr->close();
    delete newdbr;
    reader.close();
    return EXIT_SUCCESS;
}
