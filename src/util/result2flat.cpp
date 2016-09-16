#include <cstdio>
#include <Parameters.h>

#include "DBReader.h"
#include "Debug.h"
#include "Util.h"


int result2flat(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    Debug(Debug::WARNING) << "Query file is " <<  par.db1 << "\n";
    std::string queryHeaderDB =  par.db1 + "_h";
    DBReader<std::string> querydb_header(queryHeaderDB.c_str(), std::string(queryHeaderDB+".index").c_str());
    querydb_header.open(DBReader<std::string>::NOSORT);
    querydb_header.readMmapedDataInMemory();

    Debug(Debug::WARNING) << "Target file is " << par.db2 << "\n";
    std::string targetHeaderDB =  par.db2 + "_h";
    DBReader<std::string> targetdb_header(targetHeaderDB.c_str(), std::string(targetHeaderDB+".index").c_str());
    targetdb_header.open(DBReader<std::string>::NOSORT);
    targetdb_header.readMmapedDataInMemory();

    Debug(Debug::WARNING) << "Data file is " << par.db3 << "\n";
    DBReader<std::string> dbr_data( par.db3.c_str(), std::string( par.db3+".index").c_str());
    dbr_data.open(DBReader<std::string>::NOSORT);

    FILE *fastaFP =  fopen(par.db4.c_str(), "w");
    char header_start[] = {'>'};
    char newline[] = {'\n'};
    Debug(Debug::WARNING) << "Start writing file to " << par.db4 << "\n";
    char * dbKey = new char[par.maxSeqLen * 20];
    for(size_t i = 0; i < dbr_data.getSize(); i++){

        // Write the header, taken from the originial queryDB
        fwrite(header_start, sizeof(char), 1, fastaFP);
        std::string key = dbr_data.getDbKey(i);
        char * header_data = querydb_header.getDataByDBKey(key);
        std::string headerStr = Util::parseFastaHeader(header_data);
        fwrite(headerStr.c_str(), sizeof(char), headerStr.length(), fastaFP);
        //fwrite(header_data, sizeof(char), strlen(header_data) - 1, fastaFP);
        fwrite(newline, sizeof(char), 1, fastaFP);

        // write data
        char * data = dbr_data.getData(i);
        while(*data != '\0') {
            Util::parseKey(data, dbKey);
            char * header_data = targetdb_header.getDataByDBKey(dbKey);
            std::string dataStr;
            if(header_data != NULL){
                dataStr = Util::parseFastaHeader(header_data);
            }
            if( par.useHeader == true ) {
                dataStr = Util::parseFastaHeader(header_data);
                char * endLenData = Util::skipLine(data);
                size_t keyLen = strlen(dbKey);
                char * dataWithoutKey = data + keyLen;
                size_t dataToCopySize = endLenData-dataWithoutKey;
                std::string data(dataWithoutKey, dataToCopySize);
                dataStr.append(data);
            }else{
                char * startLine = data;
                char * endLine = Util::skipLine(data);
                size_t n = endLine - startLine;
                dataStr = std::string(startLine, n);
            }
			
			 // newline at the end
			 if(dataStr.length() > 0){
					if(dataStr[dataStr.length()-1] != '\n'){
						dataStr.push_back('\n');
					}
			 }
			
//            std::cout << dataStr << std::endl;
            fwrite(dataStr.c_str(), sizeof(char), dataStr.length(), fastaFP);
            dataStr.clear();
            data = Util::skipLine(data);
        }
    }
    delete [] dbKey;
    Debug(Debug::WARNING) << "Done." << "\n";

    fclose(fastaFP);
    targetdb_header.close();
    querydb_header.close();
    dbr_data.close();


    return 0;
}
