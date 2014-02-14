#include "Clustering.h"

Clustering::Clustering(std::string seqDB, std::string seqDBIndex,
        std::string alnDB, std::string alnDBIndex,
        std::string outDB, std::string outDBIndex,
        float seqIdThr){

    seqDbr = new DBReader(seqDB.c_str(), seqDBIndex.c_str());
    seqDbr->open(DBReader::SORT);

    alnDbr = new DBReader(alnDB.c_str(), alnDBIndex.c_str());
    alnDbr->open(DBReader::NOSORT);
    std::cout << alnDB << " " << alnDBIndex << "\n";

    dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str());
    dbw->open();

    this->seqIdThr = seqIdThr;
}

void Clustering::run(int mode){

    struct timeval start, end;
    gettimeofday(&start, NULL);

    std::list<set *> ret;
    Clustering::set_data set_data;
    if (mode == SET_COVER){
        std::cout << "Clustering mode: SET COVER\n";
        std::cout << "Reading the data...\n";
        set_data = read_in_set_data_set_cover();

        std::cout << "Init set cover...\n";
        SetCover setcover(set_data.set_count,
                set_data.uniqu_element_count,
                set_data.max_weight,
                set_data.all_element_count,
                set_data.element_size_lookup
                );

        for(size_t i = 0; i < set_data.set_count; i++){
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
    else if (mode == GREEDY){

        std::cout << "Clustering mode: GREEDY\n";
        std::cout << "Reading the data...\n";
        set_data = read_in_set_data_simple_clustering();

        std::cout << "Init simple clustering...\n";
        SimpleClustering simpleClustering(set_data.set_count,
                set_data.uniqu_element_count,
                set_data.all_element_count,
                set_data.element_size_lookup);

        for(size_t i = 0; i < set_data.set_count; i++){
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
    else{
        std::cerr << "ERROR: Wrong clustering mode!\n";
        exit(EXIT_FAILURE);
    }
    std::cout << "Validating results...\n";
    if(validate_result(&ret,set_data.uniqu_element_count))
        std::cout << " VALID\n";
    else
        std::cout << " NOT VALID\n";

    std::cout.flush();

    int dbSize = alnDbr->getSize();
    int cluNum = ret.size();

    std::cout << "Writing results...\n";
    writeData(ret);
    seqDbr->close();
    alnDbr->close();
    dbw->close();
    std::cout << "...done.\n";

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for clustering: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";

    std::cout << "\nSize of the database was: " << dbSize << "\n";
    std::cout << "Number of clusters: " << cluNum << "\n";

}

void Clustering::writeData(std::list<set *> ret){

    size_t BUFFER_SIZE = 1000000;
    char* outBuffer = new char[BUFFER_SIZE];
    std::list<set *>::const_iterator iterator;
    for (iterator = ret.begin(); iterator != ret.end(); ++iterator) {
        std::stringstream res;
        set::element * element =(*iterator)->elements;
        // first entry is the representative sequence
        char* dbKey = alnDbr->getDbKey(element->element_id);
        do{ 
            char* nextDbKey = alnDbr->getDbKey(element->element_id);
            res << nextDbKey << "\n";
        }while((element=element->next)!=NULL);

        std::string cluResultsOutString = res.str();
        const char* cluResultsOutData = cluResultsOutString.c_str();
        if (BUFFER_SIZE < strlen(cluResultsOutData)){
            std::cerr << "Tried to process the clustering list for the query " << dbKey << " , length of the list = " << ret.size() << "\n";
            std::cerr << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " << cluResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters -> output buffer is already huge ;-)\n";
            continue;
        }
        memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length()*sizeof(char));
        dbw->write(outBuffer, cluResultsOutString.length(), dbKey);
    }
    delete[] outBuffer;
}


bool Clustering::validate_result(std::list<set *> * ret,unsigned int uniqu_element_count){
    std::list<set *>::const_iterator iterator;
    unsigned int result_element_count=0;
    for (iterator = ret->begin(); iterator != ret->end(); ++iterator) {
        set::element * element =(*iterator)->elements;
        do{ 
            result_element_count++;
        }while((element=element->next)!=NULL);
    }
    return (uniqu_element_count==result_element_count)? true : false;
}


Clustering::set_data Clustering::read_in_set_data_set_cover(){

    Clustering::set_data ret_struct;

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
            size_t curr_element = alnDbr->getId(dbKey);
            if (curr_element == UINT_MAX){
                std::cerr << "ERROR: Element " << dbKey << " not contained in the alignment index!\n";
                exit(1);
            }
            int score = atoi(strtok(NULL, "\t")); // score
            float qcov = atof(strtok(NULL, "\t")); // query sequence coverage
            float dbcov = atof(strtok(NULL, "\t")); // db sequence coverage
            float seqId = atof(strtok(NULL, "\t")); // sequence identity
            double eval = atof(strtok(NULL, "\n")); // e-value
            // add an edge if it meets the thresholds
            // the sequence itself has already beed added
            if (curr_element != i && seqId >= seqIdThr){
                idsCnt[curr_element]++;
                element_buffer[element_counter++]=curr_element;
                element_size[curr_element]++;
                ret_struct.all_element_count++;
            }
            // next db key
            dbKey = strtok(NULL, "\t");
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

Clustering::set_data Clustering::read_in_set_data_simple_clustering(){

    Clustering::set_data ret_struct;

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
            size_t curr_element = alnDbr->getId(dbKey);
            if (curr_element == UINT_MAX){
                std::cerr << "ERROR: Element " << dbKey << " not contained in the alignment index!\n";
                exit(1);
            }
            // skip the rest
            int score = atoi(strtok(NULL, "\t")); // score
            float qcov = atof(strtok(NULL, "\t")); // query sequence coverage
            float dbcov = atof(strtok(NULL, "\t")); // db sequence coverage
            float seqId = atof(strtok(NULL, "\t")); // sequence identity
            double eval = atof(strtok(NULL, "\n")); // e-value
            // add an edge if it meets the thresholds
            // the sequence itself has already beed added
            if (curr_element != i && seqId > seqIdThr){
                idsCnt[curr_element]++;
                element_buffer[element_counter++]=curr_element;
                element_size[curr_element]++;
                ret_struct.all_element_count++;
            }
           dbKey = strtok(NULL, "\t"); // next db key
        }

        unsigned int * elements = new unsigned int[element_counter];
        memcpy(elements,element_buffer,sizeof(int)*element_counter);

        sets[i] = elements;
        set_size[i]=element_counter;

    }
    return ret_struct;
}

