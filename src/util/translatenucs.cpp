#include <unistd.h>
#include <string>
#include <limits.h>
#include <stdlib.h>

#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#include "TranslateNucl.h"

int translatenucs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    std::string in_header_filename = std::string(par.db1 + "_h");
    std::string in_header_index_filename = std::string(par.db1 + "_h.index");
    std::string out_header_filename  = std::string(par.db2 + "_h");
    std::string out_header_index_filename = std::string(par.db2 + "_h.index");

    // set links to header
    Debug(Debug::INFO) << "Set sym link from " << in_header_filename << " to " << out_header_filename << "\n";
    char *abs_in_header_filename = realpath(in_header_filename.c_str(), NULL);
    symlink(abs_in_header_filename, out_header_filename.c_str());
    free(abs_in_header_filename);
    char *abs_in_header_index_filename = realpath(in_header_index_filename.c_str(), NULL);
    Debug(Debug::INFO) << "Set sym link from " << in_header_index_filename << " to " << out_header_index_filename << "\n";
    symlink(abs_in_header_index_filename, out_header_index_filename.c_str());
    free(abs_in_header_index_filename);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str());
    writer.open();

    size_t entries = reader.getSize();
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
    char* aa = new char[par.maxSeqLen/3 + 1];
    for (size_t i = 0; i < entries; ++i) {
        unsigned int key = reader.getDbKey(i);
        char* data = reader.getData(i);

        // ignore null char at the end
        // needs to be int in order to be able to check
        int length = reader.getSeqLens(i) - 1;

        if((data[length] != '\n' && length % 3 != 0) && (data[length - 1] == '\n' && (length - 1) % 3 != 0)) {
            Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is not divisible by three. Adjust length to (lenght=" <<  length - (length % 3) << ").\n";
            length = length - (length % 3);
        }

        if(length < 3)  {
            Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is too short. Skipping entry.\n";
            continue;
        }
//        std::cout << data << std::endl;
        translateNucl.translate(aa, data, length);
        aa[length/3] = '\n';
//        std::cout << aa << std::endl;
        writer.writeData(aa, (length / 3) + 1, (char *) SSTR(key).c_str());
    }
    delete[] aa;

    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}

