#include "Clustering.h"
#include "ClusteringAlgorithms.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"
#include <sys/time.h>

Clustering::Clustering(std::string seqDB, std::string seqDBIndex,
                       std::string alnDB, std::string alnDBIndex,
                       std::string outDB, std::string outDBIndex,
                       int validateClustering,
                       unsigned int maxIteration,
                       int similarityScoreType,  int threads){

    Debug(Debug::WARNING) << "Init...\n";
    Debug(Debug::INFO) << "Opening sequence database...\n";
    seqDbr = new DBReader<unsigned int>(seqDB.c_str(), seqDBIndex.c_str());
    seqDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);

    Debug(Debug::INFO) << "Opening alignment database...\n";
    alnDbr = new DBReader<unsigned int>(alnDB.c_str(), alnDBIndex.c_str());
    alnDbr->open(DBReader<unsigned int>::NOSORT);
    this->validate = validateClustering;
    Debug(Debug::INFO) << "done.\n";
    this->maxIteration=maxIteration;
    this->similarityScoreType=similarityScoreType;
    this->threads = threads;
    this->outDB = outDB;
    this->outDBIndex = outDBIndex;
}

Clustering::~Clustering() {
    delete seqDbr;
    delete alnDbr;
}


void Clustering::run(int mode) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    DBWriter * dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1);
    dbw->open();

    std::map<unsigned int, std::vector<unsigned int> > ret;
    if  (mode == Parameters::GREEDY){
        ClusteringAlgorithms* greedyincremental= new ClusteringAlgorithms(seqDbr,alnDbr,threads,similarityScoreType,maxIteration);
        ret =greedyincremental->execute(2);
        Debug(Debug::INFO) << "Writing results...\n";
        writeData(dbw, ret);
        Debug(Debug::INFO) << "...done.\n";
        gettimeofday(&end, NULL);
        int sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime for clustering: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
        delete greedyincremental;
    }else if (mode == Parameters::SET_COVER){
        Debug(Debug::INFO) << "Clustering mode: Set cover\n";
        ClusteringAlgorithms* setCover= new ClusteringAlgorithms(seqDbr,alnDbr,threads,similarityScoreType,maxIteration);
        ret = setCover->execute(1);
        writeData(dbw, ret);
        gettimeofday(&end, NULL);
        int sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime for clustering: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
        delete setCover;
    } else if (mode == Parameters::CONNECTED_COMPONENT){
        Debug(Debug::INFO) << "Clustering mode: connected component\n";
        Debug(Debug::INFO) << "Maxiteration " <<maxIteration<<"\n";
        ClusteringAlgorithms* connctedComponent= new ClusteringAlgorithms(seqDbr,alnDbr,threads,similarityScoreType,maxIteration);
        ret =connctedComponent->execute(3);
        writeData(dbw, ret);
        gettimeofday(&end, NULL);
        int sec = end.tv_sec - start.tv_sec;
        Debug(Debug::INFO) << "\nTime for clustering: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";
        delete connctedComponent;
    }else{
        Debug(Debug::ERROR)  << "ERROR: Wrong clustering mode!\n";
        EXIT(EXIT_FAILURE);
    }

//    if (validate == 1){
//        Debug(Debug::INFO) << "Validating results...\n";
//        if(validate_result(&ret,ret.size()))
//            Debug(Debug::INFO) << " VALID\n";
//        else
//            Debug(Debug::INFO) << " NOT VALID\n";
//    }

    unsigned int dbSize = alnDbr->getSize();
    unsigned int seqDbSize = seqDbr->getSize();
    unsigned int cluNum = ret.size();

    seqDbr->close();
    alnDbr->close();
    dbw->close();
    delete dbw;

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "\nTime for clustering: " << (sec / 60) << " m " << (sec % 60) << "s\n\n";

    Debug(Debug::INFO) << "\nSize of the sequence database: " << seqDbSize << "\n";
    Debug(Debug::INFO) << "Size of the alignment database: " << dbSize << "\n";
    Debug(Debug::INFO) << "Number of clusters: " << cluNum << "\n";
}

void Clustering::writeData(DBWriter *dbw, std::map<unsigned int, std::vector<unsigned int> > ret){
    std::map<unsigned int, std::vector<unsigned int> >::const_iterator iterator;
    for (iterator = ret.begin(); iterator != ret.end(); ++iterator) {
        std::stringstream res;
        std::vector<unsigned int> elements  = (*iterator).second;
        // first entry is the representative sequence
        for(size_t i = 0; i < elements.size(); i++){
            unsigned int nextDbKey = seqDbr->getDbKey(elements[i] );
            res << nextDbKey << "\n";
        }
        unsigned int dbKey = seqDbr->getDbKey((*iterator).first);
        std::string cluResultsOutString = res.str();
        dbw->writeData(cluResultsOutString.c_str(), cluResultsOutString.length(), SSTR(dbKey).c_str());
    }


}


bool Clustering::validate_result(std::list<set *> * ret,unsigned int uniqu_element_count){
    std::list<set *>::const_iterator iterator;
    unsigned int result_element_count=0;
    int* control = new int [uniqu_element_count+1];
    memset(control, 0, sizeof(int)*(uniqu_element_count+1));
    for (iterator = ret->begin(); iterator != ret->end(); ++iterator) {
        set::element * element =(*iterator)->elements;
        do{
            control[element->element_id]++;
            result_element_count++;
        }while((element=element->next)!=NULL);
    }
    int notin = 0;
    int toomuch = 0;
    for (unsigned int i = 0; i < uniqu_element_count; i++){
        if (control[i] == 0){
            Debug(Debug::INFO) << "id " << i << " (key " << seqDbr->getDbKey(i) << ") is missing in the clustering!\n";
            Debug(Debug::INFO) << "len = " <<  seqDbr->getSeqLens()[i] << "\n";
            Debug(Debug::INFO) << "alignment results len = " << strlen(alnDbr->getDataByDBKey(seqDbr->getDbKey(i))) << "\n";
            notin++;
        }
        else if (control[i] > 1){
            Debug(Debug::INFO) << "id " << i << " (key " << seqDbr->getDbKey(i) << ") is " << control[i] << " times in the clustering!\n";
            toomuch = toomuch + control[i];
        }
    }
    if (notin > 0)
        Debug(Debug::INFO) << "not in the clustering: " << notin << "\n";
    if (toomuch > 0)
        Debug(Debug::INFO) << "multiple times in the clustering: " << toomuch << "\n";
    delete[] control;
    if (uniqu_element_count==result_element_count)
        return true;
    else{
        Debug(Debug::ERROR) << "unique_element_count: " << uniqu_element_count
        << ", result_element_count: " << result_element_count << "\n";
        return false;
    }
}


