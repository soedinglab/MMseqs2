#include "DBReader.h"
#include "Debug.h"
#include <stdio.h>
void printUsage(){
    std::string usage("\nMerge multiple ffindex files based on simular id into one file. \n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE: ffindex_database_merge ffindexDB fastaOutDB [ffindexHeaderDB]\n");
    Debug(Debug::ERROR) << usage;
}

void parseArgs(int argc, const char** argv,
               std::string* ffindexSeqDB,
               std::string* fastaOutDB,
               std::string* ffindexHeaderDB){
    if (argc < 2){
        printUsage();
        exit(EXIT_FAILURE);
    }
    ffindexSeqDB->assign(argv[1]);
    fastaOutDB->assign(argv[2]);

    if ( argc > 3){
        ffindexHeaderDB->assign(argv[3]);
    }
}



int main (int argc, const char * argv[])
{
    
    std::string ffindexSeqDB = "";
    std::string fastaOutDB = "";
    std::string ffindexHeaderDB = "";

    parseArgs(argc, argv, &ffindexSeqDB, &fastaOutDB, &ffindexHeaderDB);
    Debug(Debug::WARNING) << "Data file is " << ffindexSeqDB << "\n";
    DBReader dbr_data(ffindexSeqDB.c_str(), std::string(ffindexSeqDB+".index").c_str());
    dbr_data.open(DBReader::NOSORT);
    DBReader * dbr_header = NULL;
    if(ffindexHeaderDB.length() > 0) {
        Debug(Debug::WARNING) << "Header file is " << ffindexHeaderDB << "\n";
        dbr_header = new DBReader(ffindexHeaderDB.c_str(), std::string(ffindexHeaderDB+".index").c_str());
        dbr_header->open(DBReader::NOSORT);
    }
    FILE *fastaFP =  fopen(fastaOutDB.c_str(), "w");
    char header_start[] = {'>'};
    char newline[] = {'\n'};
    Debug(Debug::WARNING) << "Start writing file to " << fastaOutDB << "\n";

    for(size_t i = 0; i < dbr_data.getSize(); i++){
        
        fwrite(header_start, sizeof(char), 1, fastaFP);
        char * key = dbr_data.getDbKey(i);
        if(dbr_header != NULL){
            char * header_data =dbr_header->getDataByDBKey(key);
            fwrite(header_data, sizeof(char), strlen(header_data) - 1, fastaFP);
        }else{
            fwrite(key, sizeof(char), strlen(key), fastaFP);
        }
        fwrite(newline, sizeof(char), 1, fastaFP);
        // write data
        char * data = dbr_data.getData(i);
        fwrite(data, sizeof(char), strlen(data), fastaFP);
    }
    Debug(Debug::WARNING) << "Done." << "\n";

    fclose(fastaFP);
    
    dbr_data.close();
    if (dbr_header != NULL) {
        dbr_header->close();
        delete dbr_header;
    }

    return 0;
}
