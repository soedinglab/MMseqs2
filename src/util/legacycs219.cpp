#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>
#include <cstdio>

#include <string>

#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"


// from hh-suite cslib
const int kIntToChar[] = {
        33, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 47, 48, 49, 50, 51,
        52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68,
        69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
        85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
        101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
        117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132,
        133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
        149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
        165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
        181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196,
        197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
        213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228,
        229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244,
        245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 42, 45, 45
};

int legacycs219(int argn, const char **argv)
{
    std::string usage;
    usage.append("Turns a cs219 FFindex database into the lagacy cs219 format.\n");
    usage.append("USAGE: <ffindexCS219InDB> <ffindexA3MInDB> <lagacyCS219Out>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");
    std::vector<MMseqsParameter> orf_par = {
        Parameters::PARAM_V
    };

    Parameters par;
    par.parseParameters(argn, argv, usage, orf_par, 3);

    Debug::setDebugLevel(par.verbosity);

    std::string in_filename = std::string(par.db1 + ".ffdata");
    std::string in_index_filename = std::string(par.db1 + ".ffindex");

    DBReader reader(in_filename.c_str(), in_index_filename.c_str());
    reader.open(DBReader::NOSORT);

    std::string a3m_filename = std::string(par.db2 + ".ffdata");
    std::string a3m_index_filename = std::string(par.db2 + ".ffindex");

    DBReader a3m_reader(a3m_filename.c_str(), a3m_index_filename.c_str());
    a3m_reader.open(DBReader::NOSORT);

    std::string out_filename  = std::string(par.db3 + ".cs219");
    std::string out_sizes_filename = std::string(par.db3 + ".cs219.sizes");

    FILE* out_file = Util::openFileOrDie(out_filename.c_str(), "wb");
    FILE* out_size = Util::openFileOrDie(out_sizes_filename.c_str(), "w");

    size_t entries = reader.getSize();
    unsigned int* entries_length = reader.getSeqLens();

    size_t num_cs = 0;
    for (size_t i = 0; i < entries; ++i) {
        // get the fasta header from the a3m file, skip the first char (its always a #)
        char* a3m = a3m_reader.getData(i) + 1;
        // find the end of the header but include the newline
        char* header_end = strchr(a3m, '\n') + 1;
        ptrdiff_t header_length = header_end - a3m;

        char* data = reader.getData(i);

        // ignore the null byte at the end
        unsigned int length = entries_length[i] - 1;

        char cs[length];
        for(size_t i = 0; i < length; ++i) {
            cs[i] = kIntToChar[(uint8_t) data[i]];
        }

        fputc('>', out_file);
        fwrite(a3m, sizeof(char), header_length, out_file);
        fwrite(cs, sizeof(char), length, out_file);
        fputc('\n', out_file);
        num_cs += length;
    }


    fprintf(out_size, "%zu %zu", entries, num_cs);

    fclose(out_size);
    fclose(out_file);
    reader.close();

    return EXIT_SUCCESS;
}

