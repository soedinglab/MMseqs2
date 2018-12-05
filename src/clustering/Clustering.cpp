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
    Debug(Debug::INFO) << "Init...\n";
    Debug(Debug::INFO) << "Opening sequence database...\n";
    seqDbr = new DBReader<unsigned int>(seqDB.c_str(), seqDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX);
    seqDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);

    Debug(Debug::INFO) << "Opening alignment database...\n";
    alnDbr = new DBReader<unsigned int>(alnDB.c_str(), alnDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnDbr->open(DBReader<unsigned int>::NOSORT);

    Debug(Debug::INFO) << "done.\n";
}

Clustering::~Clustering() {
    delete seqDbr;
    delete alnDbr;
}


void Clustering::run(int mode) {
    Timer timer;

    DBWriter *dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1, compressed, Parameters::DBTYPE_CLUSTER_RES);
    dbw->open();

    std::unordered_map<unsigned int, std::vector<unsigned int>> ret;
    ClusteringAlgorithms *algorithm = new ClusteringAlgorithms(seqDbr, alnDbr,
                                                               threads, similarityScoreType,
                                                               maxIteration);

    if (mode == Parameters::GREEDY) {
        Debug(Debug::INFO) << "Clustering mode: Greedy\n";
        ret = algorithm->execute(2);
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
        Debug(Debug::ERROR) << "ERROR: Wrong clustering mode!\n";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Writing results...\n";
    writeData(dbw, ret);
    Debug(Debug::INFO) << "...done.\n";
    Debug(Debug::INFO) << "Time for clustering: " << timer.lap() << "\n";

    delete algorithm;

    size_t dbSize = alnDbr->getSize();
    size_t seqDbSize = seqDbr->getSize();
    size_t cluNum = ret.size();

    seqDbr->close();
    alnDbr->close();
    dbw->close();
    delete dbw;

    Debug(Debug::INFO) << "Total time: " << timer.lap() << "\n";
    Debug(Debug::INFO) << "\nSize of the sequence database: " << seqDbSize << "\n";
    Debug(Debug::INFO) << "Size of the alignment database: " << dbSize << "\n";
    Debug(Debug::INFO) << "Number of clusters: " << cluNum << "\n";
}

void Clustering::writeData(DBWriter *dbw, const std::unordered_map<unsigned int, std::vector<unsigned int>> &ret) {
    std::unordered_map<unsigned int, std::vector<unsigned int> >::const_iterator iterator;
    std::string resultStr;
    resultStr.reserve(1024*1024*1024);
    char buffer[32];
    for (iterator = ret.begin(); iterator != ret.end(); ++iterator) {
        std::vector<unsigned int> elements = (*iterator).second;
        // first entry is the representative sequence
        for (size_t i = 0; i < elements.size(); i++) {
            unsigned int nextDbKey = seqDbr->getDbKey(elements[i]);
            char * outpos = Itoa::u32toa_sse2(nextDbKey, buffer);
            resultStr.append(buffer, (outpos - buffer - 1) );
            resultStr.push_back('\n');
        }
        unsigned int dbKey = seqDbr->getDbKey((*iterator).first);
        dbw->writeData(resultStr.c_str(), resultStr.length(), dbKey);
        resultStr.clear();
    }
}
