#include <algorithm>
#include <string>
#include <vector>

#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/peer_access.cuh"

#include "dbdata.hpp"
#include "cudasw4.cuh"
#include "config.hpp"
namespace b64 {
#include "base64.h"
}
#include "marv.h"


size_t getMaxTempBytes(int maxSubjectLength) {
    int deviceId = 0;
    int numSMs = 0;
    cudaGetDevice(&deviceId); CUERR
    cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, deviceId); CUERR

    constexpr int maxGroupSize = 8; // from getMultiTileGroupRegConfigForPSSM_Gapless
    constexpr int threadBlockSize = 512;
    const int numGroupsPerBlock = threadBlockSize / maxGroupSize;
    const size_t tempStorageElementsPerGroup = SDIV(maxSubjectLength, 4);
    const size_t tempStorageElementsPerBlock = tempStorageElementsPerGroup * numGroupsPerBlock;
    return tempStorageElementsPerBlock * sizeof(float2) * numSMs;
}

cudasw4::ScanType mapCswScanType(Marv::AlignmentType type) {
    switch (type) {
        case Marv::AlignmentType::GAPLESS:
            return cudasw4::ScanType::Gapless;
        case Marv::AlignmentType::SMITH_WATERMAN:
            return cudasw4::ScanType::SW_Endpos;
        case Marv::AlignmentType::GAPLESS_SMITH_WATERMAN:
            return cudasw4::ScanType::GaplessPlusSW_Endpos;
        default:
            return cudasw4::ScanType::Gapless;
    }
}

Marv::Marv(size_t dbEntries, int alphabetSize, int maxSeqLength, size_t maxSeqs, Marv::AlignmentType alignmentType) : dbEntries(dbEntries), alphabetSize(alphabetSize), dbmanager(NULL), alignmentType(alignmentType) {
    std::vector<int> deviceIds = getDeviceIds();
    helpers::PeerAccess peerAccess(deviceIds, false);
    peerAccess.enableAllPeerAccesses();

    cudasw4::MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = 128ull * 1024ull * 1024ull;
    memoryConfig.maxBatchSequences = 10'000'000;
    // memoryConfig.maxTempBytes = 4ull * 1024ull * 1024ull * 1024ull;
    memoryConfig.maxTempBytes = getMaxTempBytes(maxSeqLength);
    memoryConfig.maxGpuMem = std::numeric_limits<size_t>::max();

    const int maxResults = std::min((int)maxSeqs, cudasw4::MaxNumberOfResults::value());
    const bool verbose = false;
    KernelConfigFilenames kernelConfigFilenames;
    //set the following to overwrite the hardcoded config
    // kernelConfigFilenames.gapless = "configfileA.txt";
    // kernelConfigFilenames.sw = "configfileB.txt";

    cudasw = static_cast<void*>(new cudasw4::CudaSW4(
        deviceIds,
        maxResults,
        cudasw4::BlosumType::BLOSUM62_20,
        memoryConfig,
        verbose,
        kernelConfigFilenames
    ));
    cudasw4::CudaSW4* sw = static_cast<cudasw4::CudaSW4*>(cudasw);
    sw->setScanType(mapCswScanType(alignmentType));
}
std::vector<std::shared_ptr<GpuDatabaseAllocationBase>> allocations_all;

Marv::~Marv() {
    allocations_all.clear();
    delete static_cast<cudasw4::CudaSW4*>(cudasw);
    // delete static_cast<cudasw4::MMseqsDB*>(db);
}

std::vector<int> Marv::getDeviceIds() {
    std::vector<int> deviceIds;
    int num = 0;
    cudaGetDeviceCount(&num); CUERR
    for (int i = 0; i < num; i++) {
        deviceIds.push_back(i);
    }
    return deviceIds;
}

void* Marv::loadDb(char* data, size_t* offset, int32_t* length, size_t dbByteSize) {
    return static_cast<void*>(new cudasw4::MMseqsDB(cudasw4::loadMMseqsDB(
        dbEntries, data, offset, length, dbByteSize
    )));
}

void* Marv::loadDb(char* data, size_t dbByteSize, void* otherdb) {
    cudasw4::MMseqsDB* db = static_cast<cudasw4::MMseqsDB*>(otherdb);
    const cudasw4::DBdataMetaData& meta = db->getData().getMetaData();
    return static_cast<void*>(new cudasw4::ExternalDB(cudasw4::loadExternalDB(
        dbEntries, dbByteSize, meta
    )));
}

std::vector<std::string> split(const std::string &str, const std::string &sep) {
    std::vector<std::string> arr;

    char *cstr = strdup(str.c_str());
    const char* csep = sep.c_str();
    char *rest;
    char *current = strtok_r(cstr, csep, &rest);
    while (current != NULL) {
        arr.emplace_back(current);
        current = strtok_r(NULL, csep, &rest);
    }
    free(cstr);

    return arr;
}

void Marv::setDbWithAllocation(void* dbhandle, const std::string& allocationinfo) {
    auto parts = split(allocationinfo, ":");

    cudaIpcMemHandle_t h1, h2, h3;
    char* charData;
    cudasw4::SequenceLengthT* lengths;
    size_t* offsets;
    size_t numChars, numSubjects;

    std::string decode;

    decode = b64::base64_decode(parts[0].data(), parts[0].length());
    memcpy((unsigned char *)(&h1), (unsigned char *)decode.data(), decode.length());
    cudaIpcOpenMemHandle((void **)&charData, h1, cudaIpcMemLazyEnablePeerAccess);
    CUERR
    decode = b64::base64_decode(parts[1].data(), parts[1].length());
    memcpy((unsigned char *)(&h2), (unsigned char *)decode.data(), decode.length());
    cudaIpcOpenMemHandle((void **)&lengths, h2, cudaIpcMemLazyEnablePeerAccess);
    CUERR
    decode = b64::base64_decode(parts[2].data(), parts[2].length());
    memcpy((unsigned char *)(&h3), (unsigned char *)decode.data(), decode.length());
    cudaIpcOpenMemHandle((void **)&offsets, h3, cudaIpcMemLazyEnablePeerAccess);
    CUERR
    numChars = strtoull(parts[3].c_str(), NULL, 10);
    numSubjects = strtoull(parts[4].c_str(), NULL, 10);
    std::vector<std::shared_ptr<GpuDatabaseAllocationBase>> allocations_remote;
    allocations_remote.emplace_back(std::make_shared<GpuDatabaseAllocationView>(GpuDatabaseAllocationView(charData, lengths, offsets, numChars, numSubjects)));

    cudasw4::MMseqsDB* db = static_cast<cudasw4::MMseqsDB*>(dbhandle);
    auto doNothingDeleter = [](cudasw4::MMseqsDB* ptr){ /* do nothing */ };
    std::shared_ptr<cudasw4::MMseqsDB> dbPtr(static_cast<cudasw4::MMseqsDB*>(db), doNothingDeleter);

    cudasw4::CudaSW4* sw = static_cast<cudasw4::CudaSW4*>(cudasw);
    // OpaqueAllocationManager* manager = static_cast<OpaqueAllocationManager*>(allocationhandle);
    sw->setDatabase(dbPtr, allocations_remote);
}

void Marv::setDb(void* dbhandle) {
    cudasw4::MMseqsDB* db = static_cast<cudasw4::MMseqsDB*>(dbhandle);
    auto doNothingDeleter = [](cudasw4::MMseqsDB* ptr){ /* do nothing */ };
    std::shared_ptr<cudasw4::MMseqsDB> dbPtr(static_cast<cudasw4::MMseqsDB*>(db), doNothingDeleter);
    cudasw4::CudaSW4* sw = static_cast<cudasw4::CudaSW4*>(cudasw);
    sw->setDatabase(dbPtr);
}

std::string Marv::getDbMemoryHandle() {
    cudasw4::CudaSW4* sw = static_cast<cudasw4::CudaSW4*>(cudasw);
    // sw->printDBInfo();
    // sw->printDBLengthPartitions();
    sw->prefetchDBToGpus();
    allocations_all = sw->getFullGpuDBAllocations();
    cudaIpcMemHandle_t h1, h2, h3;
    // char* charData;
    // cudasw4::SequenceLengthT* lengths;
    // size_t* offsets;
    size_t numChars, numSubjects;
    // for(const auto& alloc : allocations_all){
        const auto& alloc = allocations_all[0];
        cudaIpcGetMemHandle(&h1, alloc->getCharData());
        CUERR
        cudaIpcGetMemHandle(&h2, alloc->getLengthData());
        CUERR
        cudaIpcGetMemHandle(&h3, alloc->getOffsetData());
        CUERR
        // charData = alloc->getCharData();
        // lengths = alloc->getLengthData();
        // offsets = alloc->getOffsetData();
        numChars = alloc->getNumChars();
        numSubjects = alloc->getNumSubjects();
    // }

    std::vector<std::string> handles;
    std::string enc1 = b64::base64_encode(&h1, sizeof(cudaIpcMemHandle_t));
    std::string enc2 = b64::base64_encode(&h2, sizeof(cudaIpcMemHandle_t));
    std::string enc3 = b64::base64_encode(&h3, sizeof(cudaIpcMemHandle_t));

    std::string res;
    res.append(enc1);
    res.append(1, ':');
    res.append(enc2);
    res.append(1, ':');
    res.append(enc3);
    res.append(1, ':');
    res.append(std::to_string(numChars));
    res.append(1, ':');
    res.append(std::to_string(numSubjects));
    res.append(1, '\n');

    return res;
}

void Marv::printInfo() {
    cudasw4::CudaSW4* sw = static_cast<cudasw4::CudaSW4*>(cudasw);
    sw->printDBInfo();
    sw->printDBLengthPartitions();
}

void Marv::prefetch() {
    cudasw4::CudaSW4* sw = static_cast<cudasw4::CudaSW4*>(cudasw);
    sw->prefetchDBToGpus();
}

void Marv::startTimer() {
    cudasw4::CudaSW4* sw = static_cast<cudasw4::CudaSW4*>(cudasw);
    sw->totalTimerStart();
}

void Marv::stopTimer() {
    cudasw4::CudaSW4* sw = static_cast<cudasw4::CudaSW4*>(cudasw);
    auto totalBenchmarkStats = sw->totalTimerStop();
}

//sequence must be encoded
Marv::Stats Marv::scan(const char* sequence, size_t sequenceLength, int8_t* pssm, Result* results) {
    cudasw4::CudaSW4* sw = static_cast<cudasw4::CudaSW4*>(cudasw);
    cudasw4::EncodedQueryView queryView(sequence, sequenceLength);
    cudasw4::ScanResult scanResult = sw->scan(queryView, pssm);
    for (size_t i = 0; i < scanResult.scores.size(); i++) {
        results[i] = Result(
            scanResult.referenceIds[i],
            scanResult.scores[i],
            alignmentType != GAPLESS ? scanResult.endPositions[i].x : -1,
            alignmentType != GAPLESS ? scanResult.endPositions[i].y : -1
        );
    }
    Stats stats;
    stats.results = scanResult.scores.size();
    stats.numOverflows = scanResult.stats.numOverflows;
    stats.seconds = scanResult.stats.seconds;
    stats.gcups = scanResult.stats.gcups;
    return stats;
}
