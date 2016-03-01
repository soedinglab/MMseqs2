#include <cstddef>
#include <stdio.h>
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Prefiltering.h"
#include "Util.h"
typedef std::map<unsigned int, std::string *>::iterator SwapIt;

void readAllKeysIntoMap(std::string datafile,
                        std::map<unsigned int, std::string *> &map) {
    char dbKey[255 + 1];
    std::ifstream resultFile(datafile);
    std::string line;
    while(std::getline(resultFile, line)){
        size_t i;
        for(i = 0; i < line.size() - 1; i++ ){
            if(line[i] != '\0') {
                break;
            }
        }
        Util::parseKey((char*)line.c_str() + i, dbKey);
        unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);

        if(dbKey[0]=='\0')
            continue;
        SwapIt it = map.find(key);
        if (it == map.end()) {
            map[key] = new std::string();
        }
    }
    resultFile.close();
}

void processSplit(DBReader<unsigned int> &dbr,
                  FILE * dataFile, std::map<unsigned int, std::string *> &map,
                  size_t startIndex, size_t domainSize) {
    std::string resultData;
    for (size_t i = startIndex; i < (startIndex + domainSize); i++) {
        char dbKey[255 + 1];
        std::string outerKey = SSTR(dbr.getDbKey(i));
        int c1;
        while(( c1=fgetc(dataFile)) != EOF && c1 != (int) '\0'){
            char file1Char = (char) c1;
            resultData.push_back(file1Char);
        }
        resultData.push_back('\0');
        char * data = (char *) resultData.c_str();
        while (*data != '\0') {
            // extract key from results (ids must be always at the first position)
            Util::parseKey(data, dbKey);
            unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
            std::string *entry = NULL;
            entry = map[key];
            // write data to map
            entry->append(outerKey);
            // next db key
            char *endPosOfId = data + Util::skipNoneWhitespace(data);
            data = Util::skipLine(data);
            entry->append(endPosOfId, data);
        }
        resultData.clear();
    }
}

int swapresults (int argc, const char * argv[]){
    std::string usage("Swaps results of ffindex database. A -> A, B, C to A->A, B->A, C->A \n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de).\n\n");
    usage.append("USAGE: <ffindexDB> <fastaDB> [ffindexHeaderDB]\n");
    Parameters par;
    par.parseParameters(argc, argv, usage, par.swapresults, 2);
    size_t splitSize = par.split;
    Debug(Debug::INFO) << "FFindex input file is " << par.db1 << "\n";
    std::pair<std::string, std::string> name = Util::databaseNames(par.db1);
    Debug(Debug::INFO) << "Start to swap results. Write to " << par.db2 << ".\n";
    size_t entries_num = 0;
    std::vector<std::pair<std::string, std::string> > filesToDelete;
    std::map<unsigned int, std::string *> swapMap;

    // read all keys
    readAllKeysIntoMap(name.first, swapMap);
    DBReader<unsigned int> dbr(name.first.c_str(), name.second.c_str());
    dbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    FILE * dataFile = fopen(name.first.c_str(), "r");
    for (size_t split = 0; split < splitSize; split++) {
        // create and open db write
        // create splite file name
        std::string splitName(par.db2);
        splitName.append("_").append(SSTR(split));
        std::pair<std::string, std::string> splitNames = Util::databaseNames(splitName);
        filesToDelete.push_back(std::pair<std::string, std::string>(splitNames.first, splitNames.second));

        DBWriter splitWrite(splitNames.first.c_str(), splitNames.second.c_str(), 1);
        splitWrite.open();

        size_t startIndex = 0;
        size_t domainSize = 0;
        Util::decomposeDomain(dbr.getSize(), split, splitSize, &startIndex, &domainSize);
        Debug(Debug::INFO) << "Process split " << split << " from " << startIndex << " to " << (startIndex + domainSize) << " ... ";

        processSplit(dbr, dataFile, swapMap, startIndex, domainSize);
        // Write sorted results and delete memory
        for (SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
            splitWrite.write((char *) iterator->second->c_str(), iterator->second->size(),
                             (char *) SSTR(iterator->first).c_str(), 0);
            entries_num++;
            // remove just the value (string *) not the keys
            // the keys are needed for the merging step later
            delete iterator->second;
            iterator->second = new std::string();
        }
        Debug(Debug::INFO) << "Done.\n";
        splitWrite.close();
    }
    dbr.close();
    fclose(dataFile);
    for (SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
        delete iterator->second;
    }
    // merge output of all swap splits
    Prefiltering::mergeOutput(par.db2, std::string(par.db2 + ".index"), filesToDelete);
    for (size_t i = 0; i < filesToDelete.size(); i++) {
        remove(filesToDelete[i].first.c_str());
        remove(filesToDelete[i].second.c_str());
    }
    return 0;
}
