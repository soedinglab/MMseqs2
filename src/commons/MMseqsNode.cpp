#include "MMseqsNode.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Parameters.h"
#include "Util.h"
#include "FastSort.h"
#include "DBReader.h"

#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <unistd.h>
#include <sys/stat.h>

#ifdef OPENMP
#include <omp.h>
#endif

unsigned int MMseqsNode::rank = 0;
unsigned int MMseqsNode::count = 1;
// the first line of every barrier file, so a file that a killed run left is told apart from one of this run
static std::string runId;

void MMseqsNode::init() {
    const char *procId = getenv("SLURM_PROCID");
    const char *taskCount = getenv("SLURM_STEP_NUM_TASKS");
    if (procId != NULL && taskCount != NULL) {
        rank = strtoul(procId, NULL, 10);
        count = strtoul(taskCount, NULL, 10);
    }
    const char *runVariables[] = { "SLURM_JOB_ID", "SLURM_STEP_ID", "SLURM_RESTART_COUNT" };
    runId.clear();
    for (size_t i = 0; i < sizeof(runVariables) / sizeof(runVariables[0]); i++) {
        const char *value = getenv(runVariables[i]);
        runId += std::string(value != NULL ? value : "") + "\t";
    }
    runId += "\n";
    if (count == 0 || rank >= count) {
        Debug(Debug::ERROR) << "Invalid SLURM_PROCID " << rank << " for " << count << " tasks\n";
        EXIT(EXIT_FAILURE);
    }
    if (isMaster() == false) {
        Parameters &par = Parameters::getInstance();
        par.verbosity = Debug::ERROR;
        Debug::setDebugLevel(Debug::ERROR);
    }
}

void MMseqsNode::publishFile(const std::string &file, const std::string &content) {
    const std::string tmpFile = file + ".tmp";
    FILE *handle = FileUtil::openAndDelete(tmpFile.c_str(), "w");
    if (fwrite(content.c_str(), sizeof(char), content.size(), handle) != content.size()) {
        Debug(Debug::ERROR) << "Cannot write " << tmpFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fflush(handle) != 0 || fsync(fileno(handle)) != 0) {
        Debug(Debug::ERROR) << "Cannot flush " << tmpFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(handle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmpFile << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (rename(tmpFile.c_str(), file.c_str()) != 0) {
        Debug(Debug::ERROR) << "Cannot rename " << tmpFile << " to " << file << "\n";
        EXIT(EXIT_FAILURE);
    }
}

std::string MMseqsNode::waitForFile(const std::string &file) {
    while (FileUtil::fileExists(file.c_str()) == false) {
        sleep(1);
    }
    FILE *handle = FileUtil::openFileOrDie(file.c_str(), "r", true);
    struct stat sb;
    if (fstat(fileno(handle), &sb) != 0) {
        Debug(Debug::ERROR) << "Cannot read " << file << "\n";
        EXIT(EXIT_FAILURE);
    }
    std::string content(sb.st_size, '\0');
    if (fread(&content[0], sizeof(char), content.size(), handle) != content.size()) {
        Debug(Debug::ERROR) << "Cannot read " << file << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(handle) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << file << "\n";
        EXIT(EXIT_FAILURE);
    }
    return content;
}

std::string MMseqsNode::gatherFiles(const std::string &prefix) {
    std::string content;
    for (unsigned int node = 0; node < count; node++) {
        const std::string file = prefix + SSTR(node);
        content += waitForFile(file);
        FileUtil::remove(file.c_str());
    }
    return content;
}

// a file of another run is read again until this run replaces it
static std::string waitForRunFile(const std::string &file) {
    std::string content = MMseqsNode::waitForFile(file);
    while (content.compare(0, runId.size(), runId) != 0) {
        sleep(1);
        content = MMseqsNode::waitForFile(file);
    }
    return content.substr(runId.size());
}

static std::string gatherRunFiles(const std::string &prefix) {
    std::string content;
    for (unsigned int node = 0; node < MMseqsNode::count; node++) {
        const std::string file = prefix + SSTR(node);
        content += waitForRunFile(file);
        FileUtil::remove(file.c_str());
    }
    return content;
}

std::string MMseqsNode::barrier(const std::string &file, const std::string &content) {
    publishFile(file + "." + SSTR(rank), runId + content);
    if (isMaster()) {
        publishFile(file, runId + gatherRunFiles(file + "."));
    }
    return waitForRunFile(file);
}

void MMseqsNode::addExchangeFlush(std::string &indexData, size_t &writtenBytes, const std::vector<size_t> &nodeOffsets) {
    indexData += SSTR(rank);
    for (size_t node = 0; node < nodeOffsets.size(); node++) {
        indexData += "\t" + SSTR(writtenBytes + nodeOffsets[node]);
    }
    indexData += "\n";
    writtenBytes += nodeOffsets.back();
}

std::string MMseqsNode::closeExchange(FILE *file, const std::string &fileName, const std::string &indexFile, const std::string &indexData) {
    if (fflush(file) != 0 || fsync(fileno(file)) != 0 || fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    return barrier(indexFile, indexData);
}

void MMseqsNode::openExchange(const std::string &index, const std::string &prefix, std::vector<ExchangeBlock> &blocks, std::vector<FILE *> &files) {
    files.assign(count, NULL);
    std::vector<const char *> columns(count + 2);
    for (char *line = (char *) index.c_str(); *line != '\0'; line = Util::skipLine(line)) {
        Util::getWordsOfLine(line, columns.data(), columns.size());
        const unsigned int node = Util::fast_atoi<unsigned int>(columns[0]);
        const size_t pos = Util::fast_atoi<size_t>(columns[1 + rank]);
        const size_t end = Util::fast_atoi<size_t>(columns[2 + rank]);
        if (pos == end) {
            continue;
        }
        if (files[node] == NULL) {
            files[node] = FileUtil::openFileOrDie((prefix + "." + SSTR(node)).c_str(), "r", true);
        }
        const ExchangeBlock block = { files[node], pos, end, std::vector<char>(), 0, 0, 0, 0 };
        blocks.push_back(block);
    }
}

size_t MMseqsNode::readExchangeRange(std::vector<ExchangeBlock> &blocks, size_t bufferMemory, size_t (*recordSize)(const char *, size_t),
                                     bool (*compare)(const ClusterPair &, const ClusterPair &), const ClusterPair &last,
                                     std::vector<char> &records, std::vector<size_t> *recordOffsets, const std::string &prefix) {
    static const size_t EXCHANGE_BUFFER_SIZE = 1024 * 1024;
    const size_t bufferSize = std::max(EXCHANGE_BUFFER_SIZE, bufferMemory / std::max<size_t>(1, blocks.size()));
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < blocks.size(); i++) {
        ExchangeBlock &block = blocks[i];
        if (block.pos == block.end || (block.bufferEnd - block.bufferPos >= bufferSize / 2 && block.recordEnd > block.bufferPos)) {
            continue;
        }
        memmove(block.buffer.data(), block.buffer.data() + block.bufferPos, block.bufferEnd - block.bufferPos);
        block.bufferEnd -= block.bufferPos;
        block.recordEnd -= block.bufferPos;
        block.lastRecord -= block.bufferPos;
        block.bufferPos = 0;
        // a record larger than the buffer doubles it
        do {
            if (block.bufferEnd == block.buffer.size()) {
                block.buffer.resize(std::min(std::max(bufferSize, 2 * block.buffer.size()), block.bufferEnd + block.end - block.pos));
            }
            const size_t readBytes = std::min(block.buffer.size() - block.bufferEnd, block.end - block.pos);
            for (size_t bytesRead = 0; bytesRead < readBytes;) {
                const ssize_t count = pread(fileno(block.file), block.buffer.data() + block.bufferEnd + bytesRead, readBytes - bytesRead, block.pos + bytesRead);
                if (count <= 0) {
                    Debug(Debug::ERROR) << "Cannot read " << prefix << "\n";
                    EXIT(EXIT_FAILURE);
                }
                bytesRead += count;
            }
            block.pos += readBytes;
            block.bufferEnd += readBytes;
            for (size_t size; (size = recordSize(block.buffer.data() + block.recordEnd, block.bufferEnd - block.recordEnd)) > 0; block.recordEnd += size) {
                block.lastRecord = block.recordEnd;
            }
        } while (block.recordEnd == 0 && block.pos < block.end);
        if (block.pos == block.end && block.recordEnd != block.bufferEnd) {
            Debug(Debug::ERROR) << "Cannot read " << prefix << ": a record ends after its block\n";
            EXIT(EXIT_FAILURE);
        }
    }
    // the records up to the last buffered record of every block that is not read to its end are all buffered
    ClusterPair rangeLast = last;
    for (size_t i = 0; i < blocks.size(); i++) {
        if (blocks[i].pos < blocks[i].end) {
            ClusterPair blockLast;
            memcpy(&blockLast, blocks[i].buffer.data() + blocks[i].lastRecord, sizeof(ClusterPair));
            if (compare(blockLast, rangeLast)) {
                rangeLast = blockLast;
            }
        }
    }
    std::vector<size_t> rangeEnds(blocks.size());
    std::vector<size_t> blockByteOffsets(blocks.size() + 1, 0);
    std::vector<size_t> blockRecordOffsets(blocks.size() + 1, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < blocks.size(); i++) {
        const ExchangeBlock &block = blocks[i];
        size_t pos = block.bufferPos;
        size_t recordCount = 0;
        for (ClusterPair pair; pos < block.recordEnd; pos += recordSize(block.buffer.data() + pos, block.recordEnd - pos), recordCount++) {
            memcpy(&pair, block.buffer.data() + pos, sizeof(ClusterPair));
            if (compare(rangeLast, pair)) {
                break;
            }
        }
        rangeEnds[i] = pos;
        blockByteOffsets[i + 1] = pos - block.bufferPos;
        blockRecordOffsets[i + 1] = recordCount;
    }
    for (size_t i = 0; i < blocks.size(); i++) {
        blockByteOffsets[i + 1] += blockByteOffsets[i];
        blockRecordOffsets[i + 1] += blockRecordOffsets[i];
    }
    if (records.size() < blockByteOffsets.back()) {
        records.resize(blockByteOffsets.back());
    }
    if (recordOffsets != NULL) {
        recordOffsets->resize(blockRecordOffsets.back() + 1);
        recordOffsets->back() = blockByteOffsets.back();
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < blocks.size(); i++) {
        ExchangeBlock &block = blocks[i];
        memcpy(records.data() + blockByteOffsets[i], block.buffer.data() + block.bufferPos, rangeEnds[i] - block.bufferPos);
        if (recordOffsets != NULL) {
            size_t pos = blockByteOffsets[i];
            for (size_t record = blockRecordOffsets[i]; record < blockRecordOffsets[i + 1]; record++) {
                (*recordOffsets)[record] = pos;
                pos += recordSize(records.data() + pos, blockByteOffsets[i + 1] - pos);
            }
        }
        block.bufferPos = rangeEnds[i];
    }
    return blockRecordOffsets.back();
}

void MMseqsNode::closeExchangeFiles(std::vector<FILE *> &files, const std::string &prefix) {
    for (size_t node = 0; node < files.size(); node++) {
        if (files[node] != NULL && fclose(files[node]) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << prefix << "." << node << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

void MMseqsNode::writePairs(ClusterPair *pairs, size_t pairCount, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets,
                            FILE *file, const std::string &fileName, std::string &indexData, size_t &writtenBytes) {
    SORT_PARALLEL(pairs, pairs + pairCount, compare);
    std::vector<size_t> nodeOffsets(count + 1, pairCount * sizeof(ClusterPair));
    for (size_t node = 0; node < count; node++) {
        const ClusterPair nodeFirstPair = { static_cast<DBKeyType>(nodeKeyOffsets[node]), static_cast<DBKeyType>(nodeKeyOffsets[node]) };
        nodeOffsets[node] = (std::lower_bound(pairs, pairs + pairCount, nodeFirstPair, compare) - pairs) * sizeof(ClusterPair);
    }
    if (fwrite(pairs, sizeof(ClusterPair), pairCount, file) != pairCount) {
        Debug(Debug::ERROR) << "Cannot write " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    addExchangeFlush(indexData, writtenBytes, nodeOffsets);
}

std::string MMseqsNode::sendClusterPairs(const std::string &clu, bool withRepresentative, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets,
                                         const std::string &exchange, const std::string &receivedIndex, const std::string &receivedPrefix) {
    std::vector<ExchangeBlock> receivedBlocks;
    std::vector<FILE *> receivedFiles;
    if (receivedPrefix.empty() == false) {
        openExchange(receivedIndex, receivedPrefix, receivedBlocks, receivedFiles);
    }
    const std::string exchangeFile = exchange + "." + SSTR(rank);
    FILE *handle = FileUtil::openAndDelete(exchangeFile.c_str(), "w");
    std::string indexData;
    size_t writtenBytes = 0;
    const size_t memoryLimit = Util::computeMemory(Parameters::getInstance().splitMemoryLimit);
    const size_t pairCapacity = memoryLimit / 2 / sizeof(ClusterPair);
    ClusterPair *pairs = new(std::nothrow) ClusterPair[pairCapacity];
    Util::checkAllocation(pairs, "Can not allocate pairs memory in sendClusterPairs");
    size_t pairCount = 0;

    DBReader<DBKeyType> cluDb((clu + "." + SSTR(rank)).c_str(), (clu + ".index." + SSTR(rank)).c_str(), Parameters::getInstance().threads,
                              DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_LINFINITY);
    cluDb.open(DBReader<DBKeyType>::NOSORT);
    // a loaded and a prefetched load take a quarter of the memory on average; an entry also holds its index cell, key and offset
    const size_t entryBytes = cluDb.getDataSize() / std::max<size_t>(1, cluDb.getSize()) + 2 * sizeof(uint64_t) + 2 * sizeof(size_t);
    const size_t entriesPerLoad = std::max<size_t>(1, memoryLimit / 8 / entryBytes);
    std::vector<DBKeyType> entryKeys;
    std::vector<DBKeyType> repKeys;
    std::vector<char> received;
    DBKeyType previousReceivedKey = DB_KEY_INVALID;
    for (size_t loadBegin = 0; loadBegin < cluDb.getSize(); loadBegin += entriesPerLoad) {
        const size_t loadEnd = std::min(loadBegin + entriesPerLoad, cluDb.getSize());
        cluDb.loadData(loadBegin, loadEnd);
        if (loadEnd < cluDb.getSize()) {
            cluDb.prefetchData(loadEnd, std::min(loadEnd + entriesPerLoad, cluDb.getSize()));
        }
        entryKeys.resize(loadEnd - loadBegin);
        repKeys.resize(loadEnd - loadBegin);
#pragma omp parallel for schedule(static)
        for (size_t i = loadBegin; i < loadEnd; i++) {
            entryKeys[i - loadBegin] = cluDb.getDbKey(i);
            repKeys[i - loadBegin] = entryKeys[i - loadBegin];
        }
        const ClusterPair loadLast = { entryKeys.back(), entryKeys.back() };
        size_t receivedCount;
        while ((receivedCount = readExchangeRange(receivedBlocks, memoryLimit / 64, ClusterPair::recordSize, ClusterPair::compareByMember, loadLast, received, NULL, receivedPrefix)) > 0) {
            ClusterPair *receivedPairs = reinterpret_cast<ClusterPair *>(received.data());
            SORT_PARALLEL(receivedPairs, receivedPairs + receivedCount, ClusterPair::compareByMember);
#pragma omp parallel for schedule(static)
            for (size_t i = 0; i < receivedCount; i++) {
                const DBKeyType memberKey = receivedPairs[i].memberKey;
                const size_t entry = std::lower_bound(entryKeys.begin(), entryKeys.end(), memberKey) - entryKeys.begin();
                // a received member is an entry of this node once, and the entries are in key order
                if (entry == entryKeys.size() || entryKeys[entry] != memberKey || memberKey == ((i > 0) ? receivedPairs[i - 1].memberKey : previousReceivedKey)) {
                    Debug(Debug::ERROR) << "Key " << memberKey << " of " << receivedPrefix << " is no entry of " << clu << "." << rank << " in key order\n";
                    EXIT(EXIT_FAILURE);
                }
                repKeys[entry] = receivedPairs[i].repKey;
            }
            previousReceivedKey = receivedPairs[receivedCount - 1].memberKey;
        }
        for (size_t firstEntry = loadBegin; firstEntry < loadEnd;) {
            // a line has at least a digit and a newline
            if (pairCount > 0 && pairCount + cluDb.getEntryLen(firstEntry) / 2 > pairCapacity) {
                writePairs(pairs, pairCount, compare, nodeKeyOffsets, handle, exchangeFile, indexData, writtenBytes);
                pairCount = 0;
            }
            // a cluster with more pairs than fit into memory is sent in several flushes
            if (cluDb.getEntryLen(firstEntry) / 2 > pairCapacity) {
                const DBKeyType entryKey = cluDb.getDbKey(firstEntry);
                char keyBuffer[255];
                for (char *data = cluDb.getData(firstEntry, 0); *data != '\0'; data = Util::skipLine(data)) {
                    Util::parseKey(data, keyBuffer);
                    const DBKeyType key = Util::fast_atoi<DBKeyType>(keyBuffer);
                    if (withRepresentative == false && key == entryKey) {
                        continue;
                    }
                    if (pairCount == pairCapacity) {
                        writePairs(pairs, pairCount, compare, nodeKeyOffsets, handle, exchangeFile, indexData, writtenBytes);
                        pairCount = 0;
                    }
                    pairs[pairCount].repKey = repKeys[firstEntry - loadBegin];
                    pairs[pairCount].memberKey = key;
                    pairCount++;
                }
                firstEntry++;
                continue;
            }
            size_t endEntry = firstEntry;
            for (size_t flushBytes = 0; endEntry < loadEnd && pairCount + (flushBytes + cluDb.getEntryLen(endEntry)) / 2 <= pairCapacity; endEntry++) {
                flushBytes += cluDb.getEntryLen(endEntry);
            }
#pragma omp parallel
            {
                int thread_idx = 0;
#ifdef OPENMP
                thread_idx = omp_get_thread_num();
#endif
                char keyBuffer[255];
#pragma omp for schedule(dynamic, 100)
                for (size_t i = firstEntry; i < endEntry; i++) {
                    const DBKeyType entryKey = cluDb.getDbKey(i);
                    char *data = cluDb.getData(i, thread_idx);
                    size_t pairOffset = __sync_fetch_and_add(&pairCount, Util::countLines(data, cluDb.getEntryLen(i)) - (withRepresentative ? 0 : 1));
                    while (*data != '\0') {
                        Util::parseKey(data, keyBuffer);
                        const DBKeyType key = Util::fast_atoi<DBKeyType>(keyBuffer);
                        if (withRepresentative || key != entryKey) {
                            pairs[pairOffset].repKey = repKeys[i - loadBegin];
                            pairs[pairOffset].memberKey = key;
                            pairOffset++;
                        }
                        data = Util::skipLine(data);
                    }
                }
            }
            firstEntry = endEntry;
        }
    }
    const ClusterPair lastPair = { DB_KEY_INVALID, DB_KEY_INVALID };
    if (readExchangeRange(receivedBlocks, memoryLimit / 64, ClusterPair::recordSize, ClusterPair::compareByMember, lastPair, received, NULL, receivedPrefix) > 0) {
        Debug(Debug::ERROR) << "Key " << reinterpret_cast<ClusterPair *>(received.data())->memberKey << " of " << receivedPrefix << " is no entry of " << clu << "." << rank << "\n";
        EXIT(EXIT_FAILURE);
    }
    writePairs(pairs, pairCount, compare, nodeKeyOffsets, handle, exchangeFile, indexData, writtenBytes);
    cluDb.close();
    delete[] pairs;
    closeExchangeFiles(receivedFiles, receivedPrefix);
    return closeExchange(handle, exchangeFile, exchange + ".index", indexData);
}

void MMseqsNode::writeClusterIndexEntry(FILE *file, const std::string &fileName, DBKeyType key, size_t offset) {
    const uint64_t cell[2] = { static_cast<uint64_t>(key), static_cast<uint64_t>(offset) };
    if (fwrite(cell, sizeof(uint64_t), 2, file) != 2) {
        Debug(Debug::ERROR) << "Cannot write " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void MMseqsNode::removeExchanges(const std::string &prefix, const std::string &lastExchange, const std::vector<std::string> &exchanges) {
    barrier(prefix + ".merged");
    FileUtil::remove((lastExchange + "." + SSTR(rank)).c_str());
    publishFile(prefix + ".done." + SSTR(rank), runId);
    if (isMaster()) {
        gatherRunFiles(prefix + ".done.");
        for (size_t i = 0; i < exchanges.size(); i++) {
            FileUtil::remove((exchanges[i] + ".index").c_str());
        }
        FileUtil::remove((prefix + ".merged").c_str());
    }
}
