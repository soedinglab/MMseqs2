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

    seqDbr = new DBReader<DBKeyType>(seqDB.c_str(), seqDBIndex.c_str(), threads, DBReader<DBKeyType>::USE_INDEX);
    alnDbr = new DBReader<DBKeyType>(alnDB.c_str(), alnDBIndex.c_str(), threads, DBReader<DBKeyType>::USE_DATA|DBReader<DBKeyType>::USE_INDEX);
    alnDbr->open(DBReader<DBKeyType>::NOSORT);
    if (!sequenceWeightFile.empty()) {
        seqDbr->open(DBReader<DBKeyType>::SORT_BY_ID);
        SequenceWeights *sequenceWeights = new SequenceWeights(sequenceWeightFile.c_str());
        float *localid2weight = new float[seqDbr->getSize()];
        for (size_t id = 0; id < seqDbr->getSize(); id++) {
            DBKeyType key = seqDbr->getDbKey(id);
            localid2weight[id] = sequenceWeights->getWeightById(key);
        }
        seqDbr->sortIndex(localid2weight);
        delete[] localid2weight;
        delete sequenceWeights;

    } else {
        if (needSET == false) {
            seqDbr->open(DBReader<DBKeyType>::SORT_BY_LENGTH);
        } else {
            DBReader<DBKeyType> *originalseqDbr = new DBReader<DBKeyType>(seqDB.c_str(), seqDBIndex.c_str(), threads, DBReader<DBKeyType>::USE_INDEX);
            originalseqDbr->open(DBReader<DBKeyType>::NOSORT);
            DBReader<DBKeyType>::Index * seqIndex = originalseqDbr->getIndex();
            
            std::ifstream mappingStream(seqDB + ".lookup");
            std::string line;
            DBKeyType setkey = 0;
            DBKeyType maxsetkey = 0;
            DBKeyType maxkey = 0;
            while (std::getline(mappingStream, line)) {
                std::vector<std::string> split = Util::split(line, "\t");
                DBKeyType key = Util::fast_atoi<DBKeyType>(split[0].c_str());
                setkey = Util::fast_atoi<DBKeyType>(split[2].c_str());
                if (maxsetkey < setkey) {
                    maxsetkey = setkey;
                }
                maxkey = key;
            }
            DBKeyType lastKey = maxkey;
            keyToSet = new DBKeyType[lastKey+1];
            std::vector<bool> keysInSeq(lastKey+1, false);
            std::map<DBKeyType, size_t> setToLength;

            mappingStream.close();
            mappingStream.open(seqDB + ".lookup");
            line = "";
            while (std::getline(mappingStream, line)) {
                std::vector<std::string> split = Util::split(line, "\t");
                DBKeyType key = Util::fast_atoi<DBKeyType>(split[0].c_str());
                setkey = Util::fast_atoi<DBKeyType>(split[2].c_str());
                keyToSet[key] = setkey;
            }

            for (size_t id = 0; id < originalseqDbr->getSize(); id++) {
                setToLength[keyToSet[seqIndex[id].id]] += seqIndex[id].length;
                keysInSeq[seqIndex[id].id] = 1;
            }
            size_t sourceLen = maxsetkey + 1;
            seqnum = setToLength.size();
            sourceList = new(std::nothrow) DBKeyType[lastKey + 1];
            sourceOffsets = new(std::nothrow) size_t[sourceLen + 1]();
            sourceLookupTable = new(std::nothrow) DBKeyType *[sourceLen];
            size_t * sourceOffsetsDecrease = new(std::nothrow) size_t[sourceLen + 1]();

            mappingStream.close();
            mappingStream.open(seqDB + ".lookup");

            line = "";
            while (std::getline(mappingStream, line)) {
                std::vector<std::string> split = Util::split(line, "\t");
                setkey = Util::fast_atoi<DBKeyType>(split[2].c_str());
                sourceOffsets[setkey]++;
                sourceOffsetsDecrease[setkey]++;
            }
            AlignmentSymmetry::computeOffsetFromCounts(sourceOffsets, sourceLen);
            AlignmentSymmetry::setupPointers<DBKeyType>(sourceList, sourceLookupTable, sourceOffsets, sourceLen, lastKey + 1);
            
            mappingStream.close();
            mappingStream.open(seqDB + ".lookup");

            line = "";
            while (std::getline(mappingStream, line)) {
                std::vector<std::string> split = Util::split(line, "\t");
                DBKeyType key = Util::fast_atoi<DBKeyType>(split[0].c_str());
                setkey = Util::fast_atoi<DBKeyType>(split[2].c_str());
                size_t order = sourceOffsets[setkey + 1] - sourceOffsetsDecrease[setkey];
                if(keysInSeq[key] == 1) {
                    sourceList[order] = key;
                } else {
                    sourceList[order] = DB_KEY_INVALID;
                }
                sourceOffsetsDecrease[setkey]--;
            }
            char* data = (char*)malloc(
                sizeof(size_t) +
                sizeof(size_t) +
                sizeof(DBKeyType) +
                sizeof(int) +
                sizeof(unsigned int) +
                sizeof(DBReader<DBKeyType>::Index) * seqnum
            );

            std::vector<DBReader<DBKeyType>::Index*> indexStorage(seqnum);

            size_t n = 0;
            for (const auto& pairs : setToLength) {
                indexStorage[n] = new DBReader<DBKeyType>::Index;
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
            *((DBKeyType*)p) = indexStorage[seqnum-1]->id;
            p += sizeof(DBKeyType);
            *((int*)p) = originalseqDbr->getDbtype();
            p += sizeof(int);
            *((unsigned int*)p) = indexStorage[0]->length;
            p += sizeof(unsigned int);
            for (size_t i = 0; i < seqnum; ++i) {
                memcpy(
                    p + i * sizeof(DBReader<DBKeyType>::Index),
                    indexStorage[i],
                    sizeof(DBReader<DBKeyType>::Index)
                );
            }
            p += sizeof(DBReader<DBKeyType>::Index) * seqnum;
            seqDbr = DBReader<DBKeyType>::unserialize(data, threads);
            seqDbr->open(DBReader<DBKeyType>::SORT_BY_LENGTH);
            for (auto* ptr : indexStorage) {
                delete ptr;
            }
            delete[] sourceOffsetsDecrease;
            delete originalseqDbr;
        }
    }


}

Clustering::~Clustering() {
    delete seqDbr;
    delete alnDbr;
    if(needSET){
        delete[] keyToSet;
        delete[] sourceOffsets;
        delete[] sourceList;
        delete[] sourceLookupTable;
    }
}


void Clustering::run(int mode) {
    Timer timer;
    
    unsigned int dbType = Parameters::DBTYPE_CLUSTER_RES;
    unsigned int dbTypeSet = DBReader<DBKeyType>::setExtendedDbtype(dbType, Parameters::DBTYPE_EXTENDED_SET);
    DBWriter *dbw;
    if(needSET) {
        dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1, compressed, dbTypeSet);
    } else {
        dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1, compressed, dbType);
    }
    dbw->open();

    std::pair<DBKeyType, DBKeyType> * ret;
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

void Clustering::writeData(DBWriter *dbw, const std::pair<DBKeyType, DBKeyType> * ret, size_t dbSize) {
    std::string resultStr;
    resultStr.reserve(1024*1024*1024);
    char buffer[32];
    DBKeyType prevRepresentativeKey = DB_KEY_INVALID;
    for(size_t i = 0; i < dbSize; i++){
        DBKeyType currRepresentativeKey = ret[i].first;
        // write query key first
        if(prevRepresentativeKey != currRepresentativeKey) {
            if(prevRepresentativeKey != DB_KEY_INVALID){ // skip first
                dbw->writeData(resultStr.c_str(), resultStr.length(), prevRepresentativeKey);
            }
            resultStr.clear();
            char *outpos = Itoa::u64toa_sse2(static_cast<uint64_t>(currRepresentativeKey), buffer);
            resultStr.append(buffer, (outpos - buffer - 1));
            resultStr.push_back('\n');
        }
        DBKeyType memberKey = ret[i].second;
        if(memberKey != currRepresentativeKey){
            char * outpos = Itoa::u64toa_sse2(static_cast<uint64_t>(memberKey), buffer);
            resultStr.append(buffer, (outpos - buffer - 1) );
            resultStr.push_back('\n');
        }

        prevRepresentativeKey = currRepresentativeKey;
    }
    if(prevRepresentativeKey != DB_KEY_INVALID){
        dbw->writeData(resultStr.c_str(), resultStr.length(), prevRepresentativeKey);
    }
}
