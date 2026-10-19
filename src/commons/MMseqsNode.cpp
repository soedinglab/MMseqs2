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

void MMseqsNode::init() {
    const char *procId = getenv("SLURM_PROCID");
    const char *taskCount = getenv("SLURM_STEP_NUM_TASKS");
    if (procId != NULL && taskCount != NULL) {
        rank = strtoul(procId, NULL, 10);
        count = strtoul(taskCount, NULL, 10);
    }
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

std::string MMseqsNode::barrier(const std::string &file, const std::string &content) {
    publishFile(file + "." + SSTR(rank), content);
    if (isMaster()) {
        publishFile(file, gatherFiles(file + "."));
    }
    return waitForFile(file);
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
        const ExchangeBlock block = { files[node], pos, end, std::vector<char>(), 0, 0 };
        blocks.push_back(block);
    }
}

const char *MMseqsNode::readExchange(ExchangeBlock &block, size_t bytes, const std::string &prefix) {
    static const size_t EXCHANGE_BUFFER_SIZE = 1024 * 1024;
    if (block.bufferEnd - block.bufferPos < bytes) {
        if (block.pos == block.end) {
            return NULL;
        }
        const size_t kept = block.bufferEnd - block.bufferPos;
        const size_t readBytes = std::min(std::max(EXCHANGE_BUFFER_SIZE, bytes) - kept, block.end - block.pos);
        if (block.buffer.size() < kept + readBytes) {
            block.buffer.resize(kept + readBytes);
        }
        memmove(block.buffer.data(), block.buffer.data() + block.bufferPos, kept);
        if (fseek(block.file, block.pos, SEEK_SET) != 0 || fread(block.buffer.data() + kept, sizeof(char), readBytes, block.file) != readBytes) {
            Debug(Debug::ERROR) << "Cannot read " << prefix << "\n";
            EXIT(EXIT_FAILURE);
        }
        block.pos += readBytes;
        block.bufferPos = 0;
        block.bufferEnd = kept + readBytes;
        if (block.bufferEnd < bytes) {
            Debug(Debug::ERROR) << "Cannot read " << prefix << ": a record ends after its block\n";
            EXIT(EXIT_FAILURE);
        }
    }
    const char *data = block.buffer.data() + block.bufferPos;
    block.bufferPos += bytes;
    return data;
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

void MMseqsNode::openPairs(PairQueue &queue, std::vector<ExchangeBlock> &blocks, std::vector<FILE *> &files, const std::string &index, const std::string &prefix) {
    openExchange(index, prefix, blocks, files);
    for (size_t block = 0; block < blocks.size(); block++) {
        queueNextPair(queue, blocks, block, prefix);
    }
}

void MMseqsNode::queueNextPair(PairQueue &queue, std::vector<ExchangeBlock> &blocks, size_t block, const std::string &prefix) {
    const char *data = readExchange(blocks[block], sizeof(ClusterPair), prefix);
    if (data != NULL) {
        QueuedPair queuedPair;
        memcpy(&queuedPair.pair, data, sizeof(ClusterPair));
        queuedPair.block = block;
        queue.push(queuedPair);
    }
}

std::string MMseqsNode::sendClusterPairs(const std::string &clu, bool withRepresentative, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets,
                                         const std::string &exchange, const std::string &receivedIndex, const std::string &receivedPrefix) {
    const CompareQueuedPair compareQueuedByMember = { ClusterPair::compareByMember };
    PairQueue receivedQueue(compareQueuedByMember);
    std::vector<ExchangeBlock> receivedBlocks;
    std::vector<FILE *> receivedFiles;
    if (receivedPrefix.empty() == false) {
        openPairs(receivedQueue, receivedBlocks, receivedFiles, receivedIndex, receivedPrefix);
    }
    const std::string exchangeFile = exchange + "." + SSTR(rank);
    FILE *handle = FileUtil::openAndDelete(exchangeFile.c_str(), "w");
    std::string indexData;
    size_t writtenBytes = 0;
    const size_t memoryLimit = Util::computeMemory(0);
    const size_t pairCapacity = memoryLimit / 2 / sizeof(ClusterPair);
    ClusterPair *pairs = new(std::nothrow) ClusterPair[pairCapacity];
    Util::checkAllocation(pairs, "Can not allocate pairs memory in sendClusterPairs");
    size_t pairCount = 0;

    DBReader<DBKeyType> cluDb((clu + "." + SSTR(rank)).c_str(), (clu + ".index." + SSTR(rank)).c_str(), Parameters::getInstance().threads,
                              DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_LINFINITY);
    cluDb.open(DBReader<DBKeyType>::NOSORT);
    // a loaded and a prefetched load take a quarter of the memory on average
    const size_t entriesPerLoad = std::max<size_t>(1, memoryLimit / 8 / (cluDb.getDataSize() / std::max<size_t>(1, cluDb.getSize()) + 1));
    std::vector<DBKeyType> repKeys;
    for (size_t loadBegin = 0; loadBegin < cluDb.getSize(); loadBegin += entriesPerLoad) {
        const size_t loadEnd = std::min(loadBegin + entriesPerLoad, cluDb.getSize());
        cluDb.loadData(loadBegin, loadEnd);
        if (loadEnd < cluDb.getSize()) {
            cluDb.prefetchData(loadEnd, std::min(loadEnd + entriesPerLoad, cluDb.getSize()));
        }
        repKeys.resize(loadEnd - loadBegin);
        for (size_t i = loadBegin; i < loadEnd; i++) {
            const DBKeyType key = cluDb.getDbKey(i);
            while (receivedQueue.empty() == false && receivedQueue.top().pair.memberKey < key) {
                const size_t block = receivedQueue.top().block;
                receivedQueue.pop();
                queueNextPair(receivedQueue, receivedBlocks, block, receivedPrefix);
            }
            repKeys[i - loadBegin] = (receivedQueue.empty() == false && receivedQueue.top().pair.memberKey == key) ? receivedQueue.top().pair.repKey : key;
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
    publishFile(prefix + ".done." + SSTR(rank), "");
    if (isMaster()) {
        gatherFiles(prefix + ".done.");
        for (size_t i = 0; i < exchanges.size(); i++) {
            FileUtil::remove((exchanges[i] + ".index").c_str());
        }
        FileUtil::remove((prefix + ".merged").c_str());
    }
}
