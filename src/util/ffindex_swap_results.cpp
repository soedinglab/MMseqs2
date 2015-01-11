#include "DBReader.h"
#include "DBWriter.h"

#include "Debug.h"
#include <cstddef>
#include <stdio.h>
#include "Util.h"

void printUsageFFindexSwapResults(){
    std::string usage("Swaps results of ffindex database. A -> A, B, C to A->A, B->A, C->A \n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de).\n\n");
    usage.append("USAGE: <ffindexDB> <fastaDB> [ffindexHeaderDB]\n");
    Debug(Debug::ERROR) << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexSeqDB, std::string* fastaOutDB) {
    if (argc < 2) {
        printUsageFFindexSwapResults();
        EXIT(EXIT_FAILURE);
    }
    ffindexSeqDB->assign(argv[1]);
    fastaOutDB->assign(argv[2]);
}

int swapresults (int argc, const char * argv[])
{

    std::string ffindexSeqDB = "";
    std::string outDB = "";
    size_t splitSize = 1;
    parseArgs(argc, argv, &ffindexSeqDB, &outDB);
    Debug(Debug::WARNING) << "FFindex input file is " << ffindexSeqDB << "\n";
    DBReader dbr(ffindexSeqDB.c_str(), std::string(ffindexSeqDB+".index").c_str());
    dbr.open(DBReader::NOSORT);

    char dbKey[255+1];
    Debug(Debug::WARNING) << "Start to swap results. Write to " << outDB << ".\n";
    size_t entries_num = 0;
    std::vector<std::pair<std::string,std::string > > filenames_to_delete;
    std::map<std::string, std::string * > swapMap;
    typedef std::map<std::string, std::string * >::iterator SwapIt;
    for(size_t split = 0; split < splitSize; split++){
        // create splite file name
        std::string out_name       = outDB + "_" + SSTR(split);
        std::string out_name_index = (out_name + ".index");
        std::cout << "Process split " << split  << " ... ";
        // create and open db write
        DBWriter splitWrite(out_name.c_str(), out_name_index.c_str(), 1);
        splitWrite.open();
        filenames_to_delete.push_back(std::pair<std::string, std::string>(out_name,out_name_index));

        size_t startIndex = 0;
        size_t endIndex = 0;
        Util::decompose_domain(dbr.getSize(), split, splitSize, &startIndex, &endIndex);
        for(size_t i = startIndex; i < endIndex; i++){
            char * outerKey = dbr.getDbKey(i);
            char * data = dbr.getData(i);
            if(*data == '\0'){ // check if file contains entry
                Debug(Debug::ERROR) << "ERROR: Sequence " << outerKey
                        << " does not containe any sequence!\n";
                continue;
            }

            while (*data != '\0')
            {
                // extract key from results (ids must be always at the first position)
                Util::parseKey(data, dbKey);
                std::string * entry = NULL;
                SwapIt it = swapMap.find(dbKey);
                if(it == swapMap.end()|| it->second == NULL){
                    entry = new std::string();
                    entry->reserve(1620);
                    swapMap[dbKey] = entry;
                }else{
                    entry = swapMap[dbKey];
                }
                // write data to map
                entry->append(outerKey);
                // next db key
                char * endPosOfId    = data + Util::skipNoneWhitespace(data);
                data = Util::skipLine(data);
                entry->append(endPosOfId, data);
            }
        }
        // write results and delete memory
        for(SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
            splitWrite.write((char *) iterator->second->c_str(),  iterator->second->size(), (char *) iterator->first.c_str() ,0);
            entries_num++;
            // remove just the value (string *) not the keys
            // the keys are needed for the merging step later
            delete iterator->second;
            iterator->second = NULL;
        }
        splitWrite.close();
        std::cout << "Done." << std::endl;
    }
    dbr.close();

    DBWriter writer(outDB.c_str(), std::string( outDB +".index").c_str());
    writer.open();
    //write new index with all ids (A -> B) of site B
    std::string tmp_name       = outDB + "_all_ids_index";
    std::string tmp_name_index = (tmp_name + ".index");
    FILE* all_index = fopen(tmp_name_index.c_str(), "w");
    for(SwapIt it = swapMap.begin(); it != swapMap.end(); it++) {
        fprintf(all_index, "%s\t%zd\t%zd\n", it->first.c_str(), 0, 0);
    }
    fclose(all_index);
    swapMap.clear();
    // make temp. DBReader with all ids
    DBReader all_ids(ffindexSeqDB.c_str(), tmp_name_index.c_str());
    all_ids.open(DBReader::NOSORT);
    writer.mergeFiles(&all_ids, filenames_to_delete, 1000000);
    all_ids.close();
    remove(tmp_name_index.c_str());
    for (size_t i = 0; i < filenames_to_delete.size(); i++) {
        remove(filenames_to_delete[i].first.c_str());
        remove(filenames_to_delete[i].second.c_str());
    }
    writer.close();

    return 0;
}
