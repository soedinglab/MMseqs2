#include "Clustering.h"
#include "ClusteringAlgorithms.h"
#include "AlignmentSymmetry.h"
#include "Debug.h"
#include "Util.h"
#include "itoa.h"
#include "Timer.h"
#include "SequenceWeights.h"
#include <fstream>

Clustering::Clustering(const std::string &seqDB, const std::string &seqDBIndex,
                       const std::string &alnDB, const std::string &alnDBIndex,
                       const std::string &outDB, const std::string &outDBIndex,
                       const std::string &sequenceWeightFile,
                       unsigned int maxIteration, int similarityScoreType, int threads, int compressed, bool needSET) : needSET(needSET),
                                                               maxIteration(maxIteration),
                                                               similarityScoreType(similarityScoreType),
                                                               threads(threads),
                                                               compressed(compressed),
                                                               outDB(outDB),
                                                               outDBIndex(outDBIndex) {

    seqDbr = new DBReader<unsigned int>(seqDB.c_str(), seqDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX);
    alnDbr = new DBReader<unsigned int>(alnDB.c_str(), alnDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnDbr->open(DBReader<unsigned int>::NOSORT);
    if (!sequenceWeightFile.empty()) {
        seqDbr->open(DBReader<unsigned int>::SORT_BY_ID);
        SequenceWeights *sequenceWeights = new SequenceWeights(sequenceWeightFile.c_str());
        float *localid2weight = new float[seqDbr->getSize()];
        for (size_t id = 0; id < seqDbr->getSize(); id++) {
            size_t key = seqDbr->getDbKey(id);
            localid2weight[id] = sequenceWeights->getWeightById(key);
        }
        seqDbr->sortIndex(localid2weight);
        delete[] localid2weight;
        delete sequenceWeights;

    } else {
        if (needSET == false) {
            seqDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);
        } else {
            DBReader<unsigned int> *originalseqDbr = new DBReader<unsigned int>(seqDB.c_str(), seqDBIndex.c_str(), threads, DBReader<unsigned int>::USE_INDEX);
            originalseqDbr->open(DBReader<unsigned int>::NOSORT);
            DBReader<unsigned int>::Index * seqIndex = originalseqDbr->getIndex();
            
            std::ifstream mappingStream(seqDB + ".lookup");
            std::string line;
            unsigned int setkey = 0;
            unsigned int maxsetkey = 0;
            unsigned int maxkey = 0;
            while (std::getline(mappingStream, line)) {
                std::vector<std::string> split = Util::split(line, "\t");
                unsigned int key = strtoul(split[0].c_str(), NULL, 10);
                setkey = strtoul(split[2].c_str(), NULL, 10);
                if (maxsetkey < setkey) {
                    maxsetkey = setkey;
                }
                maxkey = key;
            }
            unsigned int lastKey = maxkey;
            keyToSet = new unsigned int[lastKey+1];
            std::vector<bool> keysInSeq(lastKey+1, false);
            std::map<unsigned int, unsigned int> setToLength;

            mappingStream.close();
            mappingStream.open(seqDB + ".lookup");
            line = "";
            while (std::getline(mappingStream, line)) {
                std::vector<std::string> split = Util::split(line, "\t");
                unsigned int key = strtoul(split[0].c_str(), NULL, 10);
                setkey = strtoul(split[2].c_str(), NULL, 10);
                keyToSet[key] = setkey;
            }

            for (size_t id = 0; id < originalseqDbr->getSize(); id++) {
                setToLength[keyToSet[seqIndex[id].id]] += seqIndex[id].length;
                keysInSeq[seqIndex[id].id] = 1;
            }
            unsigned int sourceLen = maxsetkey + 1;
            seqnum = setToLength.size();
            sourceList = new(std::nothrow) unsigned int[lastKey];
            sourceOffsets = new(std::nothrow) size_t[sourceLen + 1]();
            sourceLookupTable = new(std::nothrow) unsigned int *[sourceLen];
            size_t * sourceOffsetsDecrease = new(std::nothrow) size_t[sourceLen + 1]();

            mappingStream.close();
            mappingStream.open(seqDB + ".lookup");

            line = "";
            while (std::getline(mappingStream, line)) {
                std::vector<std::string> split = Util::split(line, "\t");
                setkey = strtoul(split[2].c_str(), NULL, 10);
                sourceOffsets[setkey]++;
                sourceOffsetsDecrease[setkey]++;
            }
            AlignmentSymmetry::computeOffsetFromCounts(sourceOffsets, sourceLen);
            AlignmentSymmetry::setupPointers<unsigned int>(sourceList, sourceLookupTable, sourceOffsets, sourceLen, lastKey);
            
            mappingStream.close();
            mappingStream.open(seqDB + ".lookup");

            line = "";
            while (std::getline(mappingStream, line)) {
                std::vector<std::string> split = Util::split(line, "\t");
                unsigned int key = strtoul(split[0].c_str(), NULL, 10);
                setkey = strtoul(split[2].c_str(), NULL, 10);
                size_t order = sourceOffsets[setkey + 1] - sourceOffsetsDecrease[setkey];
                if(keysInSeq[key] == 1) {
                    sourceList[order] = key;
                } else {
                    sourceList[order] = UINT_MAX;
                }
                sourceOffsetsDecrease[setkey]--;
            }
            char* data = (char*)malloc(
                sizeof(size_t) +
                sizeof(size_t) +
                sizeof(unsigned int) +
                sizeof(int) +
                sizeof(unsigned int) +
                sizeof(DBReader<unsigned int>::Index) * seqnum
            );

            std::vector<DBReader<unsigned int>::Index*> indexStorage(seqnum);

            size_t n = 0;
            for (const auto& pairs : setToLength) {
                indexStorage[n] = new DBReader<unsigned int>::Index;
                indexStorage[n]->id = pairs.first;
                indexStorage[n]->length = pairs.second;
                indexStorage[n]->offset = 0;
                n++;
            }

            char* p = data;
            *((size_t*)p) = seqnum;
            p += sizeof(size_t);
            *((size_t*)p) = 0;
            p += sizeof(size_t);
            *((unsigned int*)p) = indexStorage[seqnum-1]->id;
            p += sizeof(unsigned int);
            *((int*)p) = originalseqDbr->getDbtype();
            p += sizeof(int);
            *((unsigned int*)p) = indexStorage[0]->length;
            p += sizeof(unsigned int);
            for (size_t i = 0; i < seqnum; ++i) {
                memcpy(
                    p + i * sizeof(DBReader<unsigned int>::Index),
                    indexStorage[i],
                    sizeof(DBReader<unsigned int>::Index)
                );
            }
            p += sizeof(DBReader<unsigned int>::Index) * seqnum;
            seqDbr = DBReader<unsigned int>::unserialize(data, threads);
            seqDbr->open(DBReader<unsigned int>::SORT_BY_LENGTH);
            for (auto* ptr : indexStorage) {
                delete ptr;
            }
        }
    }


}

Clustering::~Clustering() {
    delete seqDbr;
    delete alnDbr;
    if(needSET){
        delete keyToSet;
        delete sourceOffsets;
        delete sourceList;
        delete[] sourceLookupTable;
    }
}


void Clustering::run(int mode) {
    Timer timer;
    
    unsigned int dbType = Parameters::DBTYPE_CLUSTER_RES;
    unsigned int dbTypeSet = DBReader<unsigned int>::setExtendedDbtype(dbType, Parameters::DBTYPE_EXTENDED_SET);
    DBWriter *dbw;
    if(needSET) {
        dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1, compressed, dbTypeSet);
    } else {
        dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1, compressed, dbType);
    }
    dbw->open();

    std::pair<unsigned int, unsigned int> * ret;
    ClusteringAlgorithms *algorithm = new ClusteringAlgorithms(seqDbr, alnDbr,
                                                               threads, similarityScoreType,
                                                               maxIteration, keyToSet, sourceOffsets, sourceLookupTable, sourceList, seqnum, needSET);

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
    for(size_t i = 1; i < seqDbSize; i++){
        cluNum += (ret[i].first != ret[i-1].first);
    }
    Debug(Debug::INFO) << "Total time: " << timer.lap() << "\n";
    Debug(Debug::INFO) << "\nSize of the sequence database: " << seqDbSize << "\n";
    Debug(Debug::INFO) << "Size of the alignment database: " << dbSize << "\n";
    Debug(Debug::INFO) << "Number of clusters: " << cluNum << "\n\n";

    Debug(Debug::INFO) << "Writing results ";
    writeData(dbw, ret, seqDbSize);
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
