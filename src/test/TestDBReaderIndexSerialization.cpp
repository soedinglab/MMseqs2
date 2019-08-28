#include "Debug.h"
#include "DBReader.h"

const char* binary_name = "test_dbreaderindexserialization";

int main (int, const char**) {
    DBReader<unsigned int> reader("", "/Users/mirdita/tmp/db.index", 1, DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << reader.getSize() << " " << reader.getAminoAcidDBSize() << "\n";
    Debug(Debug::INFO) << reader.getIndex()[0].id  << " " << reader.getIndex()[0].offset  << " " << reader.getIndex()[0].length  << "\n";

    char* data = DBReader<unsigned int>::serialize(reader);
    DBReader<unsigned int>* newdbr = DBReader<unsigned int>::unserialize(data, 1);
    newdbr->open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << newdbr->getSize() << " " << newdbr->getAminoAcidDBSize() << "\n";
    Debug(Debug::INFO) << newdbr->getIndex()[0].id  << " " << newdbr->getIndex()[0].offset  << " " << newdbr->getIndex()[0].length << "\n";
    free(data);

    newdbr->close();
    delete newdbr;
    reader.close();
    return EXIT_SUCCESS;
}
