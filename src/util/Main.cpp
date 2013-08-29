#include <iostream>
#include <time.h>
#include <unistd.h>
#include <string>
#include <sstream>
#include <map>
#include <set>

#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"


void printUsage(){

    std::string usage("\nCheck for inconsistent sets. Create a new consistent Database.\n");
    usage.append("Run time is O(n)*Log(n) space is O(n) \n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de)\n\n");
    usage.append("USAGE: check ffindexInputDBBase ffindexOutDBBase\n");
    std::cout << usage;
}

int getIdForKey(std::string key,
                std::map<std::string,int> &string_to_id,
                std::map<int,std::string> &id_to_string ){
    static int unique_id=0;
    if(string_to_id[key]==0){
        string_to_id[key]        = ++unique_id;
        id_to_string[unique_id]  = key;
    }
    // -1 to set it to 0
    return string_to_id[key]-1;
}

void readStructure(DBReader* qdbr,
                   int ** id_matche_lookup, int * element_sizes, int dbSize,
                   std::map<std::string,int> &string_to_id,
                   std::map<int,std::string> &id_to_string
                  ){
    const int MAX_ENTRY_SIZE=800000;
    int tmp[MAX_ENTRY_SIZE];
    // reserve first dbSize ids for index keys
    for (int i = 0; i < dbSize; i++){
        char* db_key  = qdbr->getDbKey(i);
        std::string db_key_str(db_key);
        getIdForKey(db_key_str,string_to_id,id_to_string);
    }
    for (int i = 0; i < dbSize; i++){
        std::cout << qdbr->getDbKey(i) << "\n";
        char* db_key  = qdbr->getDbKey(i);
        std::string db_key_str(db_key);
        int set_id = getIdForKey(db_key_str,string_to_id,id_to_string);
        
        char* seqData = qdbr->getData(i);
        std::stringstream seqDataStream(seqData);
        std::string    data_line;
        int count_elements=0;
        while(std::getline(seqDataStream, data_line))
        {
            std::istringstream line_stream(data_line);
            std::string token;
            // fist token is id
            while(std::getline(line_stream, token, '\t')){
                int element_id = getIdForKey(token,string_to_id,id_to_string);

                std::cout << string_to_id[token] << "\n";
                std::cout << id_to_string[element_id] << "\n";
                tmp[count_elements] = element_id;
                break;
            }
            count_elements++;
            if(count_elements >= MAX_ENTRY_SIZE){
                std::cerr << "One entry Contains more than"  <<MAX_ENTRY_SIZE << std::endl;
                exit(-1);
            }
            
        }
        id_matche_lookup[set_id] = new int[count_elements];
        element_sizes[set_id] = count_elements;
        memcpy ( &id_matche_lookup[set_id][0], &tmp[0], sizeof(int) *count_elements);
        std::sort(id_matche_lookup[set_id],id_matche_lookup[set_id]+count_elements);
    }
}

void addElementIfNotInArray(const int * array,const int size,
                            int element, std::set<int> * &inconsistenceSet){
    if(std::binary_search(array,array+size, element)==false){
        if(inconsistenceSet ==NULL)
            inconsistenceSet = new std::set<int>();
        inconsistenceSet->insert(element);
    }
}

void calculateInconsistenceSets(const int ** id_matche_lookup,
                                const int * element_sizes,
                                int dbSize,
                                std::set<int> ** inconsistenceSets){
    for (int id = 0; id < dbSize; id++){
        //check if element is inside them self
        addElementIfNotInArray(id_matche_lookup[id],element_sizes[id],id,inconsistenceSets[id]);
        // iterate over all elements and check if they are consistant 
        for(int j = 0; j < element_sizes[id];j++){
            const int element_id = id_matche_lookup[id][j];
            if(element_id < dbSize )
                addElementIfNotInArray(id_matche_lookup[element_id],
                                       element_sizes[element_id],
                                       id,
                                       inconsistenceSets[element_id]);
        }
    }
}


int main (int argc, const char * argv[])
{

    std::string inputDB = "../../test_data/example";
    std::string inputDBIndex = "";
    std::string outDB = "../../test_data/merge_out_test";
    std::string outDBIndex = "";
    std::string scoringMatrixFile = "";

    inputDBIndex = inputDB + ".index";
    outDBIndex = outDB + ".index";

    DBReader* qdbr = new DBReader(inputDB.c_str(), inputDBIndex.c_str());
    qdbr->open();
    int dbSize = qdbr->getSize();
    std::map<std::string,int> string_to_id;
    std::map<int,std::string> id_to_string;
    int ** id_matche_lookup = new int *[dbSize];
    int *  element_sizes    = new int [dbSize];

    readStructure(qdbr,id_matche_lookup,element_sizes,dbSize,string_to_id,id_to_string);
    
    std::set<int> ** inconsistenceSets = new std::set<int> *[dbSize];
    memset(inconsistenceSets,0,sizeof(std::set<int> **)*dbSize);
    calculateInconsistenceSets((const int **)id_matche_lookup,element_sizes,dbSize,inconsistenceSets);
    for(int i = 0; i<dbSize; i++){
        std::set<int> * insonsistent_set = inconsistenceSets[i];
        if(insonsistent_set != NULL){
            std::cout << "Inconsistence in " << i << std::endl;

        }
    }
    
//    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1);
//    dbw->open();

    // init the substitution matrices

    qdbr->close();
    
    delete [] element_sizes;
    for(int i = 0; i < dbSize; i++){
        delete id_matche_lookup[i];
    }
    delete [] id_matche_lookup;
//    dbw->close();
    std::cout << "done.\n";

    return 0;
}
