
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <istream>
#include <cstring>
#include <stdlib.h> 
#include <time.h>

#include "LinearMultiArray.h"
#include "SetCover.h"
#include "SimpleClustering.h"
#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"

void printUsage(){

    std::string usage("\nCalculates clustering of a sequence database based on Smith Waterman alignment scores with set cover algorithm.\n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE: kClust2_setcover ffindexAlnResultsDB ffindexOutDB [opts]\n"
             "-s              \t[file]\tffindex sequence database for the greedy clustering by sequence length.\n");
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexAlnDBBase, std::string* ffindexOutDBBase, std::string* ffindexSeqDBBase){
    if (argc < 3){
        printUsage();
        exit(EXIT_FAILURE);
    }
    ffindexAlnDBBase->assign(argv[1]);
    ffindexOutDBBase->assign(argv[2]);
    int i = 3;
    while (i < argc){
        if (strcmp(argv[i], "-s") == 0){
            if (++i < argc){
                ffindexSeqDBBase->assign(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else {
            printUsage();
            std::cerr << "Wrong argument: " << argv[i] << "\n";
            exit(EXIT_FAILURE);
        }
    }
}

struct set_data {
    // one set contains pointers to the cluster member ids
    unsigned int ** sets;
    unsigned short ** weights;
    unsigned int * set_sizes;
    unsigned int * element_size_lookup;
    unsigned int set_count;
    unsigned int uniqu_element_count;
    unsigned int all_element_count;
    unsigned int max_weight;
};

// check if every element is member in only one cluster
bool validate_result(std::list<set *> * ret,unsigned int uniqu_element_count){
    std::list<set *>::const_iterator iterator;
    unsigned int result_element_count=0;
    for (iterator = ret->begin(); iterator != ret->end(); ++iterator) {
        set::element * element =(*iterator)->elements;
        do{
            result_element_count++;
        }while((element=element->next)!=NULL);
    }
    std::cout << "uniqu_element_count: " << uniqu_element_count << ", result_element_count: " << result_element_count << "\n";
    return (uniqu_element_count==result_element_count)? true : false;
}


set_data read_in_set_data_set_cover(DBReader* alnDbr){
    
    set_data ret_struct;
    
    // n = overall sequence count
    int n = alnDbr->getSize();
    ret_struct.uniqu_element_count = n;
    ret_struct.set_count = n;
    
    unsigned int * element_buffer=(unsigned int *)malloc(n * sizeof(unsigned int));
    unsigned int * element_size=new unsigned int[n+2];
    memset(&element_size[0], 0, sizeof(unsigned int)*(n+2));
    ret_struct.element_size_lookup = element_size;
    unsigned int * set_size=new unsigned int[n];
    
    ret_struct.set_sizes = set_size;
    unsigned int ** sets   = new unsigned int*[n];
    ret_struct.sets = sets;
    unsigned short ** weights = new unsigned short*[n];
    ret_struct.weights = weights;
    ret_struct.max_weight = 0;
    ret_struct.all_element_count=0;

    char* buf = new char[1000000];

    int* idsCnt = new int[n];
    memset(idsCnt, 0, sizeof(int) * n);

    for(int i = 0; i<n; i++){
        char* data = alnDbr->getData(i);
        strcpy(buf, data);
        unsigned int element_counter=0;
        char* dbKey = strtok(buf, "\t");
        // ensure that each set contains the query sequence itself
        idsCnt[i]++;
        element_buffer[element_counter++] = i;
        element_size[i]++;
        ret_struct.all_element_count++;

        while (dbKey != 0)
        {
            int curr_element = alnDbr->getId(dbKey);
            if (curr_element != i){
                idsCnt[curr_element]++;
                element_buffer[element_counter++]=curr_element;
                element_size[curr_element]++;
                ret_struct.all_element_count++;
            }
            // skip the rest
            int score = atoi(strtok(NULL, "\t")); // score
            float qcov = atof(strtok(NULL, "\t")); // query sequence coverage
            float dbcov = atof(strtok(NULL, "\t")); // db sequence coverage
            float eval = atof(strtok(NULL, "\n")); // e-value
            dbKey = strtok(NULL, "\t"); // next db key
        }

        if(ret_struct.max_weight < element_counter){
            ret_struct.max_weight = element_counter;
        }

        unsigned int * elements = new unsigned int[element_counter];
        memcpy(elements,element_buffer,sizeof(int)*element_counter);
        unsigned short * weight = new unsigned short[element_counter];
        std::fill_n(weight, element_counter, 1);

        weights[i] = weight;
        sets[i] = elements;
        set_size[i]=element_counter;

    }

    return ret_struct;
}


set_data read_in_set_data_simple_clustering(DBReader* alnDbr, DBReader* seqDbr){
    
    set_data ret_struct;
    
    // n = overall sequence count
    int n = seqDbr->getSize();
    ret_struct.uniqu_element_count = n;
    ret_struct.set_count = n;
    
    unsigned int * element_buffer=(unsigned int *)malloc(n * sizeof(unsigned int));
    unsigned int * element_size=new unsigned int[n+2];
    memset(&element_size[0], 0, sizeof(unsigned int)*(n+2));
    ret_struct.element_size_lookup = element_size;
    unsigned int * set_size=new unsigned int[n];
    
    ret_struct.set_sizes = set_size;
    unsigned int ** sets   = new unsigned int*[n];
    ret_struct.sets = sets;
    ret_struct.all_element_count=0;

    char* buf = new char[1000000];

    int* idsCnt = new int[n];
    memset(idsCnt, 0, sizeof(int) * n);

    for(int i = 0; i < n; i++){
        char* qDbKey = seqDbr->getDbKey(i);
        char* data = alnDbr->getDataByDBKey(qDbKey);
        strcpy(buf, data);
        unsigned int element_counter=0;
        char* dbKey = strtok(buf, "\t");
        // ensure that each set contains the query sequence itself
        idsCnt[i]++;
        element_buffer[element_counter++]=i;
        element_size[i]++;
        ret_struct.all_element_count++;

        while (dbKey != 0)
        {
            int curr_element = alnDbr->getId(dbKey);
            if (curr_element != i){
                idsCnt[curr_element]++;
                element_buffer[element_counter++]=curr_element;
                element_size[curr_element]++;
                ret_struct.all_element_count++;
            }
            // skip the rest
            int score = atoi(strtok(NULL, "\t")); // score
            float qcov = atof(strtok(NULL, "\t")); // query sequence coverage
            float dbcov = atof(strtok(NULL, "\t")); // db sequence coverage
            float eval = atof(strtok(NULL, "\n")); // e-value
            dbKey = strtok(NULL, "\t"); // next db key
        }

        unsigned int * elements = new unsigned int[element_counter];
        memcpy(elements,element_buffer,sizeof(int)*element_counter);

        sets[i] = elements;
        set_size[i]=element_counter;

    }

    return ret_struct;
}


int main(int argc, const char * argv[])
{

    std::string seqDB = "";
    std::string seqDBIndex = "";
    std::string alnDB = "";
    std::string alnDBIndex = "";
    std::string outDB = "";
    std::string outDBIndex = "";

    parseArgs(argc, argv, &alnDB, &outDB, &seqDB);

    alnDBIndex = alnDB + ".index";
    outDBIndex = outDB + ".index";
    if (seqDB.compare("") != 0)
        seqDBIndex = seqDB + ".index";


    std::cout << "Opening databases...\n";

    DBReader* alnDbr = new DBReader(alnDB.c_str(), alnDBIndex.c_str());
    alnDbr->open(DBReader::NOSORT);

    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str());
    dbw->open();

    DBReader* seqDbr = NULL;
    if (seqDB.compare("") != 0){
        seqDbr = new DBReader(seqDB.c_str(), seqDBIndex.c_str());
        seqDbr->open(DBReader::SORT);
    }
   
    std::list<set *> ret;
    set_data set_data;
    if (seqDB.compare("") == 0){
        std::cout << "Reading the data...\n";
        set_data = read_in_set_data_set_cover(alnDbr);

        std::cout << "Init set cover...\n";
        SetCover setcover(set_data.set_count,
                set_data.uniqu_element_count,
                set_data.max_weight,
                set_data.all_element_count,
                set_data.element_size_lookup
                );

        for(int i = 0; i < set_data.set_count; i++){
            setcover.add_set(i+1, set_data.set_sizes[i]
                    ,(const unsigned int*)set_data.sets[i],
                    (const unsigned short*)set_data.weights[i],
                    set_data.set_sizes[i]);
            delete[] set_data.sets[i];
            delete[] set_data.weights[i];

        }
        delete[] set_data.sets;
        delete[] set_data.weights;


        std::cout.flush();

        std::cout << "Clustering...\n";
        ret = setcover.execute_set_cover();
        std::cout << "done.\n";

    }
    else {
        std::cout << "Reading the data...\n";
        set_data = read_in_set_data_simple_clustering(alnDbr, seqDbr);

        std::cout << "Init simple clustering...\n";
        SimpleClustering simpleClustering(set_data.set_count,
                set_data.uniqu_element_count,
                set_data.all_element_count,
                set_data.element_size_lookup);

        for(int i = 0; i < set_data.set_count; i++){
            simpleClustering.add_set((const unsigned int*)set_data.sets[i],
                    set_data.set_sizes[i]);
            delete[] set_data.sets[i];
        }
        delete[] set_data.sets;

        std::cout.flush();

        std::cout << "Clustering...\n";
        ret = simpleClustering.execute();
        std::cout << "done.\n";

    }
    std::cout << "Validating results...\n";
    if(validate_result(&ret,set_data.uniqu_element_count))
        std::cout << " VALID\n";
    else
        std::cout << " NOT VALID\n";

    std::cout.flush();

    std::cout << "Writing results...\n" << std::endl;

    size_t BUFFER_SIZE = 1000000;
    char* outBuffer = new char[BUFFER_SIZE];
    std::list<set *>::const_iterator iterator;
    for (iterator = ret.begin(); iterator != ret.end(); ++iterator) {
        std::stringstream res;
        set::element * element =(*iterator)->elements;
        // first entry is the representative sequence
        char* dbKey = alnDbr->getDbKey(element->element_id);
        do{
            res << alnDbr->getDbKey(element->element_id) << "\n";
        }while((element=element->next)!=NULL);

        std::string cluResultsOutString = res.str();
        const char* cluResultsOutData = cluResultsOutString.c_str();
        if (BUFFER_SIZE < strlen(cluResultsOutData)){
            std::cerr << "Tried to process the clustering list for the query " << dbKey << " , the length of the list = " << ret.size() << "\n";
            std::cerr << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " << cluResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters -> output buffer is already huge ;-)\n";
            continue;
        }
        memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length()*sizeof(char));
        dbw->write(outBuffer, cluResultsOutString.length(), dbKey);
    }
    dbw->close();
    delete[] outBuffer;
    std::cout << "done\n";

}

