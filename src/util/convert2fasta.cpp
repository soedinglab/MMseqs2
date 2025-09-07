/*
 * convert2fasta
 * written by Milot Mirdita <milot@mirdita.de>
 */

#include <cstring>
#include <cstdio>

#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"

const char headerStart[] = {'>'};
const char newline[] = {'\n'};

int convert2fasta(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<KeyType> db(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<KeyType>::USE_DATA | DBReader<KeyType>::USE_INDEX);
    db.open(DBReader<KeyType>::NOSORT);

    DBReader<KeyType> db_header(par.hdr1.c_str(), par.hdr1Index.c_str(), 1, DBReader<KeyType>::USE_DATA | DBReader<KeyType>::USE_INDEX);
    db_header.open(DBReader<KeyType>::NOSORT);

    FILE* fastaFP = fopen(par.db2.c_str(), "w");
    if(fastaFP == NULL) {
        perror(par.db2.c_str());
        EXIT(EXIT_FAILURE);
    }


    DBReader<KeyType>* from = &db;
    if(par.useHeaderFile) {
        from = &db_header;
    }

    Debug(Debug::INFO) << "Start writing file to " << par.db2 << "\n";
    for(size_t i = 0; i < from->getSize(); i++){
        KeyType key = from->getDbKey(i);
        KeyType headerKey = db_header.getId(key);
        const char* headerData = db_header.getData(headerKey, 0);
        const size_t headerLen = db_header.getEntryLen(headerKey);

        fwrite(headerStart, sizeof(char), 1, fastaFP);
        fwrite(headerData, sizeof(char), headerLen - 2, fastaFP);
        fwrite(newline, sizeof(char), 1, fastaFP);

        KeyType bodyKey = db.getId(key);
        const char* bodyData = db.getData(bodyKey, 0);
        const size_t bodyLen = db.getEntryLen(bodyKey);
        fwrite(bodyData, sizeof(char), bodyLen - 2, fastaFP);
        fwrite(newline, sizeof(char), 1, fastaFP);
    }
    if (fclose(fastaFP) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << par.db2 << "\n";
        EXIT(EXIT_FAILURE);
    }
    db_header.close();
    db.close();

    return EXIT_SUCCESS;
}
