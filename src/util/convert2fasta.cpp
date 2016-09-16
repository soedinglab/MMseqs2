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

const char header_start[] = {'>'};
const char newline[] = {'\n'};

int convert2fasta(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    std::string data_filename = par.db1;
    std::string index_filename = par.db1Index;

    std::string data_filename_hdr(data_filename);
    data_filename_hdr.append("_h");

    std::string index_filename_hdr(data_filename);
    index_filename_hdr.append("_h.index");

    DBReader<std::string> db(data_filename.c_str(), index_filename.c_str());
    db.open(DBReader<std::string>::NOSORT);

    DBReader<std::string> db_header(data_filename_hdr.c_str(), index_filename_hdr.c_str());
    db_header.open(DBReader<std::string>::NOSORT);

    FILE *fastaFP =  FileUtil::openFileOrDie(par.db2.c_str(), "w", false);

    DBReader<std::string>& from = db;
    if(par.useHeaderFile) {
        from = db_header;
    }

    Debug(Debug::INFO) << "Start writing file to " << par.db2 << "\n";
    for(size_t i = 0; i < from.getSize(); i++){
        std::string key = from.getDbKey(i);

        const char* header_data = db_header.getDataByDBKey(key.c_str());

        fwrite(header_start, sizeof(char), 1, fastaFP);
        fwrite(header_data, sizeof(char), strlen(header_data) - 1, fastaFP);
        fwrite(newline, sizeof(char), 1, fastaFP);

        const char* body_data = db.getDataByDBKey(key.c_str());
        fwrite(body_data, sizeof(char), strlen(body_data) - 1, fastaFP);
        fwrite(newline, sizeof(char), 1, fastaFP);

    }

    fclose(fastaFP);
    db_header.close();
    db.close();

    return EXIT_SUCCESS;
}
