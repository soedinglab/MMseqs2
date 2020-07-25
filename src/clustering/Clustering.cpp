#include "Clustering.h"
#include "ClusteringAlgorithms.h"
#include "Debug.h"
#include "Util.h"
#include "itoa.h"
#include "Timer.h"

Clustering::Clustering(const std::string &seqDB, const std::string &seqDBIndex,
                       const std::string &alnDB, const std::string &alnDBIndex,
                       const std::string &outDB, const std::string &outDBIndex,
                       unsigned int maxIteration, int similarityScoreType, int threads, int compressed) : maxIteration(maxIteration),
                                                               similarityScoreType(similarityScoreType),
                                                               threads(threads),
                                                               compressed(compressed),
                                                               outDB(outDB),
                                                               outDBIndex(outDBIndex) {

    seqDbr = new DBReader<unsigned int>(seqDB.c_str(), seqDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX);
    seqDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);

    alnDbr = new DBReader<unsigned int>(alnDB.c_str(), alnDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnDbr->open(DBReader<unsigned int>::NOSORT);

}

Clustering::~Clustering() {
    delete seqDbr;
    delete alnDbr;
}


void Clustering::run(int mode) {
    Timer timer;
    DBWriter *dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1, compressed, Parameters::DBTYPE_CLUSTER_RES);
    dbw->open();

    std::pair<unsigned int, unsigned int> * ret;
    ClusteringAlgorithms *algorithm = new ClusteringAlgorithms(seqDbr, alnDbr,
                                                               threads, similarityScoreType,
                                                               maxIteration);

    if (mode == Parameters::GREEDY) {
        Debug(Debug::INFO) << "Clustering mode: Greedy\n";
        ret = algorithm->execute(4);
    } else if (mode == Parameters::GREEDY_MEM) {
        Debug(Debug::INFO) << "Clustering mode: Greedy Low Mem\n";
        ret = algorithm->execute(4);
    } else if (mode == Parameters::SET_COVER) {
        Debug(Debug::INFO) << "Clustering mode: Set Cover\n";
        ret = algorithm->execute(1);
    } else if (mode == Parameters::CONNECTED_COMPONENT) {
        Debug(Debug::INFO) << "Clustering mode: Connected Component\n";
        ret = algorithm->execute(3);
    } else {
        Debug(Debug::ERROR) << "Wrong clustering mode!\n";
        EXIT(EXIT_FAILURE);
    }

    Timer timerWrite;

    size_t dbSize = alnDbr->getSize();
    size_t seqDbSize = seqDbr->getSize();
    size_t cluNum = (dbSize > 0) ? 1 : 0;
    for(size_t i = 1; i < dbSize; i++){
        cluNum += (ret[i].first != ret[i-1].first);
    }
    Debug(Debug::INFO) << "Total time: " << timer.lap() << "\n";
    Debug(Debug::INFO) << "\nSize of the sequence database: " << seqDbSize << "\n";
    Debug(Debug::INFO) << "Size of the alignment database: " << dbSize << "\n";
    Debug(Debug::INFO) << "Number of clusters: " << cluNum << "\n\n";

    Debug(Debug::INFO) << "Writing results ";
    writeData(dbw, ret, dbSize);
    Debug(Debug::INFO) << timerWrite.lap() << "\n";
    delete [] ret;
    delete algorithm;

    dbw->close(false, false);
    seqDbr->close();
    alnDbr->close();
    delete dbw;

}

void Clustering::writeData(DBWriter *dbw, const std::pair<unsigned int, unsigned int> * ret, size_t dbSize) {
    std::string resultStr;
    resultStr.reserve(1024*1024*1024);
    char buffer[32];
    unsigned int prevRepresentativeKey = UINT_MAX;
    for(size_t i = 0; i < dbSize; i++){
        unsigned int currRepresentativeKey = ret[i].first;
        // write query key first
        if(prevRepresentativeKey != currRepresentativeKey) {
            if(prevRepresentativeKey != UINT_MAX){ // skip first
                dbw->writeData(resultStr.c_str(), resultStr.length(), prevRepresentativeKey);
            }
            resultStr.clear();
            char *outpos = Itoa::u32toa_sse2(currRepresentativeKey, buffer);
            resultStr.append(buffer, (outpos - buffer - 1));
            resultStr.push_back('\n');
        }
        unsigned int memberKey = ret[i].second;
        if(memberKey != currRepresentativeKey){
            char * outpos = Itoa::u32toa_sse2(memberKey, buffer);
            resultStr.append(buffer, (outpos - buffer - 1) );
            resultStr.push_back('\n');
        }

        prevRepresentativeKey = currRepresentativeKey;
    }
    if(prevRepresentativeKey != UINT_MAX){
        dbw->writeData(resultStr.c_str(), resultStr.length(), prevRepresentativeKey);
    }
}
