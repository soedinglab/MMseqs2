#include "DBReader.h"
#include "DBWriter.h"

#include "Debug.h"
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
    DBReader qdbr(ffindexSeqDB.c_str(), std::string(ffindexSeqDB+".index").c_str());
    qdbr.open(DBReader::NOSORT);

    char newline[] = {'\n'};
    char dbKey[255+1];
    Debug(Debug::WARNING) << "Start to swap results. Write to " << outDB << ".\n";
    size_t entries_num = 0;
    const unsigned int ELEMENTS_IN_RECORD = 2;
    char * words[ELEMENTS_IN_RECORD];
    std::vector<std::pair<std::string,std::string > > all_filenames;
    for(size_t split = 0; split < splitSize; split++){
        // create splite file name
        std::string out_name       = outDB + "_" + SSTR(split);
        std::string out_name_index = (out_name + ".index");
        std::cout << "Process split " << split  << " ... ";
        DBWriter write(out_name.c_str(), out_name_index.c_str(), 1);
        write.open();
        all_filenames.push_back(std::pair<std::string, std::string>(out_name,out_name_index));

        std::map<std::string, std::string * > swapMap;
        size_t startIndex = 0;
        size_t endIndex = 0;
        Util::decompose_domain(qdbr.getSize(), split, splitSize, &startIndex, &endIndex);
        for(size_t i = startIndex; i < endIndex; i++){
            char * outerKey = qdbr.getDbKey(i);
            char * data = qdbr.getData(i);
            if(*data == '\0'){ // check if file contains entry
                Debug(Debug::ERROR) << "ERROR: Sequence " << outerKey
                        << " does not containe any sequence!\n";
                continue;
            }

            while (*data != '\0')
            {
                words[0] = data;
                words[1] = data + Util::skipNoneWhitespace(data);

                ptrdiff_t keySize =  (words[1] - words[0]);
                strncpy(dbKey, data, keySize);
                dbKey[keySize] = '\0';
                //std::string entryKey(dbKey);

                std::string * entry = NULL;
//                std::cout << dbKey << std::endl;
                if(swapMap.find(dbKey) == swapMap.end()){
                    entry = new std::string();
                    entry->reserve(1620);
                    swapMap[dbKey] = entry;
                }else{
                    entry = swapMap[dbKey];
                }
                // write data to map
                entry->append(outerKey);
                //entry->push_back('\t');
                // next db key
                data = Util::skipLine(data);
                entry->append(words[1], data);
                //std::cout << swapMap[entryKey]->c_str() ;
                //std::flush(std::cout);
            }
        }
        // write results and delete memory
        size_t offset_sequence = 0;
        typedef std::map<std::string, std::string * >::iterator SwapIt;
        for(SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
            write.write((char *) iterator->second->c_str(),  iterator->second->size(), (char *) iterator->first.c_str() ,0);
            entries_num++;
            delete iterator->second;
        }
        swapMap.clear();
        write.close();
        std::cout << "Done." << std::endl;
    }

    DBWriter writer(outDB.c_str(), std::string( outDB +".index").c_str());
    writer.open();
    //TODO works only in a symetric setup. Means all all IDs in data are also in the index
    writer.mergeFiles(&qdbr, all_filenames, 1000000);
    for (size_t i = 0; i < all_filenames.size(); i++) {
        remove(all_filenames[i].first.c_str());
        remove(all_filenames[i].second.c_str());
    }
    writer.close();
    qdbr.close();

    return 0;
}
