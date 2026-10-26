#include "ClusteringAlgorithms.h"
#include "Util.h"
#include "Debug.h"
#include "AlignmentSymmetry.h"
#include "Timer.h"

#include <queue>
#include <algorithm>
#include <climits>
#include <unordered_map>
#include <FastSort.h>

#ifdef OPENMP
#include <omp.h>
#endif

ClusteringAlgorithms::ClusteringAlgorithms(DBReader<DBKeyType>* seqDbr, DBReader<DBKeyType>* alnDbr,
                                           int threads, int scoretype, int maxiterations,
                                           DBKeyType *keyToSet, size_t *sourceOffsets, DBKeyType **sourceLookupTable, DBKeyType *sourceList, size_t sourceLen, bool needSET){
    this->seqDbr=seqDbr;
    if(seqDbr->getSize() != alnDbr->getSize() && needSET == false){
        Debug(Debug::ERROR) << "Sequence db size != result db size\n";
        EXIT(EXIT_FAILURE);
    }
    this->alnDbr=alnDbr;
    this->dbSize=alnDbr->getSize();
    this->threads=threads;
    this->scoretype=scoretype;
    this->maxiterations=maxiterations;
    this->keyToSet=keyToSet;
    this->sourceOffsets=sourceOffsets;
    this->sourceLookupTable=sourceLookupTable;
    this->sourceList=sourceList;
    this->sourceLen=sourceLen;
    this->needSET=needSET;
    ///time
#ifdef MMSEQS_INT64_IDS
    this->clustersizes=new int64_t[dbSize];
#else
    this->clustersizes=new int32_t[dbSize];
#endif
    std::fill_n(clustersizes, dbSize, 0);
}

ClusteringAlgorithms::~ClusteringAlgorithms(){
    delete [] clustersizes;
}

std::pair<DBKeyType, DBKeyType> * ClusteringAlgorithms::execute(int mode) {
    // init data
    if(needSET){
        dbSize = sourceLen;
    }
    DBLocalId *assignedcluster = new(std::nothrow) DBLocalId[dbSize];
    Util::checkAllocation(assignedcluster, "Can not allocate assignedcluster memory in ClusteringAlgorithms::execute");
    std::fill_n(assignedcluster, dbSize, DB_LOCAL_ID_INVALID);

    //time
    if (mode==4 || mode==2) {
        greedyIncrementalLowMem(assignedcluster);
    } else {
        size_t elementCount = 0;
#pragma omp parallel reduction (+:elementCount)
        {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 10)
            for (size_t i = 0; i < alnDbr->getSize(); i++) {
                const char *data = alnDbr->getData(i, thread_idx);
                const size_t dataSize = alnDbr->getEntryLen(i);
                elementCount += (*data == '\0') ? 1 : Util::countLines(data, dataSize);
            }
        }
        DBLocalId * elements = new(std::nothrow) DBLocalId[elementCount];
        Util::checkAllocation(elements, "Can not allocate elements memory in ClusteringAlgorithms::execute");
        DBLocalId ** elementLookupTable = new(std::nothrow) DBLocalId*[dbSize];
        Util::checkAllocation(elementLookupTable, "Can not allocate elementLookupTable memory in ClusteringAlgorithms::execute");
        unsigned short **scoreLookupTable = new(std::nothrow) unsigned short *[dbSize];
        Util::checkAllocation(scoreLookupTable, "Can not allocate scoreLookupTable memory in ClusteringAlgorithms::execute");
        unsigned short *score = NULL;
        size_t *elementOffsets = new(std::nothrow) size_t[dbSize + 1];
        Util::checkAllocation(elementOffsets, "Can not allocate elementOffsets memory in ClusteringAlgorithms::execute");
        elementOffsets[dbSize] = 0;
        short *bestscore = new(std::nothrow) short[dbSize];
        Util::checkAllocation(bestscore, "Can not allocate bestscore memory in ClusteringAlgorithms::execute");
        std::fill_n(bestscore, dbSize, SHRT_MIN);

        readInClusterData(elementLookupTable, elements, scoreLookupTable, score, elementOffsets, elementCount);
        ClusteringAlgorithms::initClustersizes();
        if (mode == 1) {
            setCover(elementLookupTable, scoreLookupTable, assignedcluster, bestscore, elementOffsets);
        } else if (mode == 3) {
            Debug(Debug::INFO) << "connected component mode" << "\n";
            for (size_t cl_size = dbSize; cl_size > 0; cl_size--) {
                DBLocalId representative = sorted_clustersizes[cl_size - 1];
                if (assignedcluster[representative] == DB_LOCAL_ID_INVALID) {
                    assignedcluster[representative] = representative;
                    std::queue<DBLocalId> myqueue;
                    myqueue.push(representative);
                    std::queue<int> iterationcutoffs;
                    iterationcutoffs.push(0);
                    //delete clusters of members;
                    while (!myqueue.empty()) {
                        DBLocalId currentid = myqueue.front();
                        int iterationcutoff = iterationcutoffs.front();
                        assignedcluster[currentid] = representative;
                        myqueue.pop();
                        iterationcutoffs.pop();
                        size_t elementSize = (elementOffsets[currentid + 1] - elementOffsets[currentid]);
                        for (size_t elementId = 0; elementId < elementSize; elementId++) {
                            DBLocalId elementtodelete = elementLookupTable[currentid][elementId];
                            if (assignedcluster[elementtodelete] == DB_LOCAL_ID_INVALID && iterationcutoff < maxiterations) {
                                myqueue.push(elementtodelete);
                                iterationcutoffs.push((iterationcutoff + 1));
                            }
                            assignedcluster[elementtodelete] = representative;
                        }
                    }

                }
            }
        }
        //delete unnecessary datastructures
        delete [] sorted_clustersizes;
        delete [] clusterid_to_arrayposition;
        delete [] borders_of_set;


        delete [] elementLookupTable;
        delete [] elements;
        delete [] elementOffsets;
        delete [] scoreLookupTable;
        delete [] score;
        delete [] bestscore;
    }


    std::pair<DBKeyType, DBKeyType> * assignment = new std::pair<DBKeyType, DBKeyType> [dbSize];
#pragma omp parallel
    {

#pragma omp for schedule(static)
        for (size_t i = 0; i < dbSize; i++) {
            if (assignedcluster[i] == DB_LOCAL_ID_INVALID) {
                Debug(Debug::ERROR) << "there must be an error: " << seqDbr->getDbKey(i) << "\t" << i <<
                                    "\tis not assigned to a cluster\n";
                continue;
            }

            assignment[i].first = seqDbr->getDbKey(assignedcluster[i]);
            assignment[i].second = seqDbr->getDbKey(i);
        }
    }
    SORT_PARALLEL(assignment,assignment+dbSize);
    delete [] assignedcluster;
    return assignment;
}

void ClusteringAlgorithms::initClustersizes(){
    DBLocalId * setsize_abundance = new DBLocalId[maxClustersize+1];

    std::fill_n(setsize_abundance, maxClustersize+1, 0);
    //count how often a set size occurs
    for (size_t i = 0; i < dbSize; ++i) {
        setsize_abundance[clustersizes[i]]++;
    }
    //compute offsets
    borders_of_set = new size_t[maxClustersize+1];
    borders_of_set[0] = 0;
    for (size_t i = 1; i < maxClustersize+1; ++i) {
        borders_of_set[i] = borders_of_set[i-1] + setsize_abundance[i-1];
    }
    //fill array
    sorted_clustersizes = new(std::nothrow) DBLocalId[dbSize + 1];
    Util::checkAllocation(sorted_clustersizes, "Can not allocate sorted_clustersizes memory in ClusteringAlgorithms::initClustersizes");

    std::fill_n(sorted_clustersizes, dbSize+1, 0);
    clusterid_to_arrayposition = new(std::nothrow) DBLocalId[dbSize + 1];
    Util::checkAllocation(clusterid_to_arrayposition, "Can not allocate sorted_clustersizes memory in ClusteringAlgorithms::initClustersizes");

    std::fill_n(clusterid_to_arrayposition, dbSize + 1, 0);
    //reuse setsize_abundance as offset counter
    std::fill_n(setsize_abundance, maxClustersize + 1, 0);
    for (size_t i = 0; i < dbSize; ++i) {
        size_t position = borders_of_set[clustersizes[i]] + setsize_abundance[clustersizes[i]];
        sorted_clustersizes[position] = i;
        clusterid_to_arrayposition[i] = position;
        setsize_abundance[clustersizes[i]]++;
    }
    delete [] setsize_abundance;
}


void ClusteringAlgorithms::removeClustersize(DBLocalId clusterid){
    clustersizes[clusterid]=0;
    sorted_clustersizes[clusterid_to_arrayposition[clusterid]] = DB_LOCAL_ID_INVALID;
    clusterid_to_arrayposition[clusterid]=DB_LOCAL_ID_INVALID;
}

void ClusteringAlgorithms::decreaseClustersize(DBLocalId clusterid){
    const size_t oldposition=clusterid_to_arrayposition[clusterid];
    const size_t newposition=borders_of_set[clustersizes[clusterid]];
    const DBLocalId swapid=sorted_clustersizes[newposition];
    if(swapid != DB_LOCAL_ID_INVALID){
        clusterid_to_arrayposition[swapid]=oldposition;
    }
    sorted_clustersizes[oldposition]=swapid;

    sorted_clustersizes[newposition]=clusterid;
    clusterid_to_arrayposition[clusterid]=newposition;
    borders_of_set[clustersizes[clusterid]]++;
    clustersizes[clusterid]--;
}

void ClusteringAlgorithms::setCover(DBLocalId **elementLookupTable, unsigned short ** elementScoreLookupTable,
                                    DBLocalId *assignedcluster, short *bestscore, size_t *newElementOffsets) {
    for (int64_t cl_size = dbSize - 1; cl_size >= 0; cl_size--) {
        const DBLocalId representative = sorted_clustersizes[cl_size];
        if (representative == DB_LOCAL_ID_INVALID) {
            continue;
        }
//          Debug(Debug::INFO)<<alnDbr->getDbKey(representative)<<"\n";
        removeClustersize(representative);
        assignedcluster[representative] = representative;
        //delete clusters of members;
        size_t elementSize = (newElementOffsets[representative + 1] - newElementOffsets[representative]);
        for (size_t elementId = 0; elementId < elementSize; elementId++) {
            const DBLocalId elementtodelete = elementLookupTable[representative][elementId];
            // float seqId = elementScoreTable[representative][elementId];
            const short seqId = elementScoreLookupTable[representative][elementId];
            //  Debug(Debug::INFO)<<seqId<<"\t"<<bestscore[elementtodelete]<<"\n";
            // becareful of this criteria
            if (seqId > bestscore[elementtodelete]) {
                assignedcluster[elementtodelete] = representative;
                bestscore[elementtodelete] = seqId;
            }
            //Debug(Debug::INFO)<<bestscore[elementtodelete]<<"\n";
            if (elementtodelete == representative) {
                continue;
            }
            if (clustersizes[elementtodelete] < 1) {
                continue;
            }
            removeClustersize(elementtodelete);
        }

        for (size_t elementId = 0; elementId < elementSize; elementId++) {
            bool representativefound = false;
            const DBLocalId elementtodelete = elementLookupTable[representative][elementId];
            const size_t currElementSize = (newElementOffsets[elementtodelete + 1] -
                                                newElementOffsets[elementtodelete]);
            if (elementtodelete == representative) {
                clustersizes[elementtodelete] = -1;
                continue;
            }
            if (clustersizes[elementtodelete] < 0) {
                continue;
            }
            clustersizes[elementtodelete] = -1;
            //decrease clustersize of sets that contain the element
            for (size_t elementId2 = 0; elementId2 < currElementSize; elementId2++) {
                const DBLocalId elementtodecrease = elementLookupTable[elementtodelete][elementId2];
                if (representative == elementtodecrease) {
                    representativefound = true;
                }
                if (clustersizes[elementtodecrease] == 1) {
                    Debug(Debug::ERROR) << "there must be an error: " << seqDbr->getDbKey(elementtodelete) <<
                                        " deleted from " << seqDbr->getDbKey(elementtodecrease) <<
                                        " that now is empty, but not assigned to a cluster\n";
                } else if (clustersizes[elementtodecrease] > 0) {
                    decreaseClustersize(elementtodecrease);
                }
            }
            if (!representativefound) {
                Debug(Debug::ERROR) << "error with cluster:\t" << seqDbr->getDbKey(representative) <<
                                    "\tis not contained in set:\t" << seqDbr->getDbKey(elementtodelete) << ".\n";
            }
        }
    }
}

void ClusteringAlgorithms::greedyIncrementalLowMem( DBLocalId *assignedcluster) {

    const long BUFFER_SIZE = 100000; // Set this to a suitable value.
    const long numBuffers = (dbSize + BUFFER_SIZE - 1) / BUFFER_SIZE;

    // Pre-allocate buffer outside the loop to reuse it
    std::vector<std::pair<DBLocalId, std::vector<DBKeyType>>> buffer(BUFFER_SIZE);

    for (long bufferIndex = 0; bufferIndex < numBuffers; bufferIndex++) {
        long start = bufferIndex * BUFFER_SIZE;
        long end = std::min(start + BUFFER_SIZE, static_cast<long>(dbSize));

        // Clear the vectors within the buffer, but don't deallocate
        for (std::pair<DBLocalId, std::vector<DBKeyType>>& entry : buffer) {
            entry.second.clear();
        }

        // Parallel reading and parsing into buffer
#pragma omp parallel
        {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 4)
            for (long i = start; i < end; i++) {
                DBKeyType clusterKey = seqDbr->getDbKey(i);
                std::vector<DBKeyType>& keys = buffer[i - start].second;
                if(needSET) {
                    size_t start1 = sourceOffsets[clusterKey];
                    size_t end1 = sourceOffsets[clusterKey+1];
                    size_t len = end1 - start1;
                    for (size_t j = 0; j < len; ++j) {
                        DBKeyType value = sourceLookupTable[clusterKey][j];
                        if (value != DB_KEY_INVALID) {
                            const size_t alnId = alnDbr->getId(value);
                            char *data = alnDbr->getData(alnId, thread_idx);
                            while (*data != '\0') {
                                char dbKey[255 + 1];
                                Util::parseKey(data, dbKey);
                                const DBKeyType key = keyToSet[Util::fast_atoi<DBKeyType>(dbKey)];
                                keys.push_back(key);
                                data = Util::skipLine(data);
                            }
                        }
                    }
                } else {
                    const size_t alnId = alnDbr->getId(clusterKey);
                    char* data = alnDbr->getData(alnId, thread_idx);
                    while (*data != '\0') {
                        char dbKey[255 + 1];
                        Util::parseKey(data, dbKey);
                        const DBKeyType key = Util::fast_atoi<DBKeyType>(dbKey);
                        keys.push_back(key);
                        data = Util::skipLine(data);
                    }
                }

                buffer[i - start].first = i;
            }
        }

        // Sequential processing of the buffer
        for (long j = 0; j < (end - start); j++) {
            DBLocalId clusterId = buffer[j].first;
            const std::vector<DBKeyType>& keys = buffer[j].second;

            if (assignedcluster[clusterId] != DB_LOCAL_ID_INVALID) {
                continue;
            }

            if (keys.size() <= 1) {
                continue;
            }

            for (DBKeyType key : keys) {
                size_t currElement = seqDbr->getId(key);

                if (currElement != DB_ENTRY_NOT_FOUND && assignedcluster[currElement] == DB_LOCAL_ID_INVALID) {
                    assignedcluster[currElement] = clusterId;
                }
            }
        }
    }

    // correct edges that are not assigned properly
    for (size_t id = 0; id < dbSize; ++id) {
        // check if the assigned clusterid is a rep. sequence
        // if not, make it a rep. seq. again
        if(assignedcluster[id] == DB_LOCAL_ID_INVALID){
            assignedcluster[id] = id;
        }
    }
}

void ClusteringAlgorithms::readInClusterData(DBLocalId **elementLookupTable, DBLocalId *&elements,
                                             unsigned short **scoreLookupTable, unsigned short *&scores,
                                             size_t *elementOffsets, size_t totalElementCount) {
    Timer timer;
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < seqDbr->getSize(); i++) {
            const DBKeyType clusterId = seqDbr->getDbKey(i);
            if(needSET) {
                size_t start = sourceOffsets[clusterId];
                size_t end = sourceOffsets[clusterId+1];
                size_t len = end - start;
                size_t lineCounts = 0;
                for (size_t j = 0; j < len; ++j) {
                    DBKeyType value = sourceLookupTable[clusterId][j];
                    if (value != DB_KEY_INVALID) {
                        const size_t alnId = alnDbr->getId(value);
                        const char *data = alnDbr->getData(alnId, thread_idx);
                        const size_t dataSize = alnDbr->getEntryLen(alnId);
                        size_t lineCount = (*data == '\0') ? 1 : Util::countLines(data, dataSize);
                        lineCounts += lineCount;
                    }
                }
                elementOffsets[i] = lineCounts;
            } else {
                const size_t alnId = alnDbr->getId(clusterId);
                const char *data = alnDbr->getData(alnId, thread_idx);
                const size_t dataSize = alnDbr->getEntryLen(alnId);
                elementOffsets[i] = (*data == '\0') ? 1 : Util::countLines(data, dataSize);
            }
        }
    }

    // make offset table
    AlignmentSymmetry::computeOffsetFromCounts(elementOffsets, dbSize);
    // set element edge pointers by using the offset table
    AlignmentSymmetry::setupPointers<DBLocalId>(elements, elementLookupTable, elementOffsets, dbSize,
                                                   totalElementCount);
    // fill elements
    if(needSET) {
        AlignmentSymmetry::readInDataSet(alnDbr, seqDbr, elementLookupTable, NULL, 0, elementOffsets, sourceOffsets, sourceLookupTable, keyToSet, 1);
        AlignmentSymmetry::computeOffsetFromCounts(elementOffsets, dbSize);
    } else {
        AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable, NULL, 0, elementOffsets);
    }
    Debug(Debug::INFO) << "Sort entries\n";
    AlignmentSymmetry::sortElements(elementLookupTable, elementOffsets, dbSize);
    Debug(Debug::INFO) << "Find missing connections\n";

    size_t *newElementOffsets = new size_t[dbSize + 1];
    memcpy(newElementOffsets, elementOffsets, sizeof(size_t) * (dbSize + 1));
    
    // findMissingLinks detects new possible connections and updates the elementOffsets with new sizes
    const size_t symmetricElementCount = AlignmentSymmetry::findMissingLinks(elementLookupTable,
                                                                             newElementOffsets, dbSize,
                                                                             threads);
    // resize elements
    delete[] elements;
    elements = new(std::nothrow) DBLocalId[symmetricElementCount];
    Util::checkAllocation(elements, "Can not allocate elements memory in readInClusterData");
    std::fill_n(elements, symmetricElementCount, DB_LOCAL_ID_INVALID);
    // init score vector
    scores = new(std::nothrow) unsigned short[symmetricElementCount];
    Util::checkAllocation(scores, "Can not allocate scores memory in readInClusterData");
    std::fill_n(scores, symmetricElementCount, 0);
    if(needSET == false){
        Debug(Debug::INFO) << "Found " << symmetricElementCount - totalElementCount << " new connections.\n";
    }
    AlignmentSymmetry::setupPointers<DBLocalId>  (elements, elementLookupTable, newElementOffsets, dbSize, symmetricElementCount);
    AlignmentSymmetry::setupPointers<unsigned short>(scores, scoreLookupTable, newElementOffsets, dbSize, symmetricElementCount);
    //time
    Debug(Debug::INFO) << "Reconstruct initial order\n";
    alnDbr->remapData(); // need to free memory
    if(needSET) {
        AlignmentSymmetry::readInDataSet(alnDbr, seqDbr, elementLookupTable, scoreLookupTable, scoretype, elementOffsets, sourceOffsets, sourceLookupTable, keyToSet, 0);
    } else {
        AlignmentSymmetry::readInData(alnDbr, seqDbr, elementLookupTable, scoreLookupTable, scoretype, elementOffsets);
    }
    alnDbr->remapData(); // need to free memory
    Debug(Debug::INFO) << "Add missing connections\n";
    AlignmentSymmetry::addMissingLinks(elementLookupTable, elementOffsets, newElementOffsets, dbSize, scoreLookupTable);
    maxClustersize = 0;
    for (size_t i = 0; i < dbSize; i++) {
        size_t elementCount = newElementOffsets[i + 1] - newElementOffsets[i];
        maxClustersize = std::max(elementCount, maxClustersize);
        clustersizes[i] = elementCount;
    }

    memcpy(elementOffsets, newElementOffsets, sizeof(size_t) * (dbSize + 1));
    delete[] newElementOffsets;
    Debug(Debug::INFO) << "\nTime for read in: " << timer.lap() << "\n";
}
