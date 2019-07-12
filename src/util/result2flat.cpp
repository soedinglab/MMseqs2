#include <cstdio>
#include <Parameters.h>

#include "DBReader.h"
#include "Debug.h"
#include "Util.h"


int result2flat(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> querydb_header(par.hdr1.c_str(), par.hdr1Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    querydb_header.open(DBReader<unsigned int>::NOSORT);
    querydb_header.readMmapedDataInMemory();

    DBReader<unsigned int> targetdb_header(par.hdr2.c_str(), par.hdr2Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    targetdb_header.open(DBReader<unsigned int>::NOSORT);
    targetdb_header.readMmapedDataInMemory();

    DBReader<unsigned int> dbr_data(par.db3.c_str(), par.db3Index.c_str(),  1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    dbr_data.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    FILE *fastaFP = fopen(par.db4.c_str(), "w");

    char header_start[] = {'>'};
    char newline[] = {'\n'};

    char *dbKeyBuffer = new char[par.maxSeqLen * 20];
    for (size_t i = 0; i < dbr_data.getSize(); i++) {

        // Write the header, taken from the original queryDB
        fwrite(header_start, sizeof(char), 1, fastaFP);
        unsigned int key = dbr_data.getDbKey(i);
        char *header_data = querydb_header.getDataByDBKey(key, 0);

        std::string headerStr;
        if (par.useHeader == true){
            headerStr = header_data;
            if (headerStr.length() > 0) {
                if (headerStr[headerStr.length() - 1] == '\n') {
                    headerStr[headerStr.length() - 1] = ' ';
                }
            }
        }else{
            headerStr=Util::parseFastaHeader(header_data);
        }
        fwrite(headerStr.c_str(), sizeof(char), headerStr.length(), fastaFP);
        fwrite(newline, sizeof(char), 1, fastaFP);

        // write data
        char *data = dbr_data.getData(i, 0);
        while (*data != '\0') {
            // dbKeyBuffer can contain sequence
            Util::parseKey(data, dbKeyBuffer);
            const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
            char *header_data = targetdb_header.getDataByDBKey(dbKey, 0);
            std::string dataStr;
            if (par.useHeader == true && header_data != NULL && dbr_data.getDbtype() == -1)
            {
                dataStr = Util::parseFastaHeader(header_data);
                char *endLenData = Util::skipLine(data);
                size_t keyLen = strlen(dbKeyBuffer);
                char *dataWithoutKey = data + keyLen;
                size_t dataToCopySize = endLenData - dataWithoutKey;
                std::string data(dataWithoutKey, dataToCopySize);
                dataStr.append(data);
            } else {
                char *startLine = data;
                char *endLine = Util::skipLine(data);
                size_t n = endLine - startLine;
                dataStr = std::string(startLine, n);
            }

            // newline at the end
            if (dataStr.length() > 0) {
                if (dataStr[dataStr.length() - 1] != '\n') {
                    dataStr.push_back('\n');
                }
            }

//            std::cout << dataStr << std::endl;
            fwrite(dataStr.c_str(), sizeof(char), dataStr.length(), fastaFP);
            dataStr.clear();
            data = Util::skipLine(data);
        }
    }
    delete[] dbKeyBuffer;


    fclose(fastaFP);
    targetdb_header.close();
    querydb_header.close();
    dbr_data.close();


    return 0;
}
