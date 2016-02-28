#include "DBReader.h"
#include "DBWriter.h"

#include "Debug.h"
#include <cstddef>
#include <stdio.h>
#include <Prefiltering.h>
#include "Util.h"


int swapresults (int argc, const char * argv[]){

    std::string usage("Swaps results of ffindex database. A -> A, B, C to A->A, B->A, C->A \n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de).\n\n");
    usage.append("USAGE: <ffindexDB> <fastaDB> [ffindexHeaderDB]\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.swapresults, 2);

    std::string ffindexResDB = par.db1;
    std::string outDB = par.db2;
    size_t splitSize = par.split;
    Debug(Debug::INFO) << "FFindex input file is " << ffindexResDB << "\n";
    std::pair<std::string, std::string> name = Util::databaseNames(ffindexResDB);
    char dbKey[255 + 1];
    Debug(Debug::INFO) << "Start to swap results. Write to " << outDB << ".\n";
    size_t entries_num = 0;
    std::vector<std::pair<std::string, std::string> > filenames_to_delete;
    std::map<std::string, std::string *> swapMap;
    typedef std::map<std::string, std::string *>::iterator SwapIt;
    std::ifstream dataFile(name.first);
    std::string line;
    while(std::getline(dataFile, line)){
        if(line[0] == '\0'){
            size_t i;
            for(i = 0; i < line.size(); i++ ){
                if(line[i] != '\0') {
                    break;
                }
            }
            Util::parseKey((char*)line.c_str() + i, dbKey);
        }else
            Util::parseKey((char*)line.c_str(), dbKey);
        if(dbKey[0]=='\0')
            continue;
        SwapIt it = swapMap.find(dbKey);
//        std::cout << dbKey << std::endl;
        if (it == swapMap.end()) {
            swapMap[dbKey] = new std::string();
        }
    }
    dataFile.close();
    DBReader<unsigned int> dbr(name.first.c_str(), name.second.c_str());
    dbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    for (size_t split = 0; split < splitSize; split++) {
        Debug(Debug::INFO) << "Process split " << split << " ... ";
        // create and open db write
        // create splite file name
        std::string splitName(outDB);
        splitName.append("_");
        splitName.append(SSTR(split));
        std::pair<std::string, std::string> splitNames = Util::databaseNames(splitName);
        DBWriter splitWrite(splitNames.first.c_str(), splitNames.second.c_str(), 1);
        splitWrite.open();
        filenames_to_delete.push_back(std::pair<std::string, std::string>(splitNames.first, splitNames.second));

        size_t startIndex = 0;
        size_t domainSize = 0;
        Util::decomposeDomain(dbr.getSize(), split, splitSize, &startIndex, &domainSize);
        for (size_t i = startIndex; i < (startIndex + domainSize); i++) {
            std::string outerKey = SSTR(dbr.getDbKey(i));
            char *data = dbr.getData(i);
            if (*data == '\0') { // check if file contains entry
                Debug(Debug::INFO) << "\nSequence " << outerKey
                << " does not contain any sequence!\n";
                continue;
            }

            while (*data != '\0') {
                // extract key from results (ids must be always at the first position)
                Util::parseKey(data, dbKey);
                std::string *entry = NULL;
                entry = swapMap[dbKey];
                // write data to map
                entry->append(outerKey);
                // next db key
                char *endPosOfId = data + Util::skipNoneWhitespace(data);
                data = Util::skipLine(data);
                entry->append(endPosOfId, data);
            }
        }
        // Write sorted results and delete memory
        for (SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
            splitWrite.write((char *) iterator->second->c_str(), iterator->second->size(),
                             (char *) iterator->first.c_str(), 0);
            entries_num++;
            // remove just the value (string *) not the keys
            // the keys are needed for the merging step later
            delete iterator->second;
            iterator->second = new std::string();
        }
        splitWrite.close();
        Debug(Debug::INFO) << "Done.\n";
    }
    dbr.close();

    // merge output of all swap splits
    Prefiltering::mergeOutput(outDB, std::string(outDB + ".index"), filenames_to_delete);

    for (size_t i = 0; i < filenames_to_delete.size(); i++) {
        remove(filenames_to_delete[i].first.c_str());
        remove(filenames_to_delete[i].second.c_str());
    }
    for (SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
        delete iterator->second;

    }

    return 0;
}
