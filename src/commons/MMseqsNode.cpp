#include "MMseqsNode.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Parameters.h"
#include "Util.h"
#include "FastSort.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "itoa.h"

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
    // mpirun (Open MPI, then PMI of MPICH/Intel MPI) before SLURM, whose variables mpirun tasks inherit from the allocation
    const char *rankVariables[][2] = { { "OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE" }, { "PMI_RANK", "PMI_SIZE" }, { "SLURM_PROCID", "SLURM_STEP_NUM_TASKS" } };
    for (size_t i = 0; i < sizeof(rankVariables) / sizeof(rankVariables[0]); i++) {
        const char *procId = getenv(rankVariables[i][0]);
        const char *taskCount = getenv(rankVariables[i][1]);
        if (procId != NULL && taskCount != NULL) {
            rank = strtoul(procId, NULL, 10);
            count = strtoul(taskCount, NULL, 10);
            break;
        }
    }
    const char *runVariables[] = { "SLURM_JOB_ID", "SLURM_STEP_ID", "SLURM_RESTART_COUNT", "OMPI_MCA_ess_base_jobid", "PMIX_NAMESPACE" };
    runId.clear();
    for (size_t i = 0; i < sizeof(runVariables) / sizeof(runVariables[0]); i++) {
        const char *value = getenv(runVariables[i]);
        runId += std::string(value != NULL ? value : "") + "\t";
    }
    runId += "\n";
    if (count == 0 || rank >= count) {
        Debug(Debug::ERROR) << "Invalid node rank " << rank << " for " << count << " nodes\n";
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

void MMseqsNode::sleepBeforePoll(unsigned int polls) {
    usleep(10000u << std::min(polls, 6u));
}

std::string MMseqsNode::waitForFile(const std::string &file) {
    for (unsigned int polls = 0; FileUtil::fileExists(file.c_str()) == false; polls++) {
        sleepBeforePoll(polls);
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
std::string MMseqsNode::waitForRunFile(const std::string &file) {
    std::string content = waitForFile(file);
    for (unsigned int polls = 0; content.compare(0, runId.size(), runId) != 0; polls++) {
        sleepBeforePoll(polls);
        content = waitForFile(file);
    }
    return content.substr(runId.size());
}

void MMseqsNode::publishRunFile(const std::string &file, const std::string &content) {
    publishFile(file, runId + content);
}

std::string MMseqsNode::gatherRunFiles(const std::string &prefix) {
    std::string content;
    for (unsigned int node = 0; node < MMseqsNode::count; node++) {
        const std::string file = prefix + SSTR(node);
        content += waitForRunFile(file);
        FileUtil::remove(file.c_str());
    }
    return content;
}

void MMseqsNode::publishLog(FILE *file, const std::string &log, size_t key, size_t bytes) {
    if (fflush(file) != 0 || fsync(fileno(file)) != 0) {
        Debug(Debug::ERROR) << "Cannot flush " << log << "\n";
        EXIT(EXIT_FAILURE);
    }
    publishRunFile(log + ".progress", SSTR(key) + "\t" + SSTR(bytes) + "\n");
}

bool MMseqsNode::readLogProgress(const std::string &log, size_t &key, size_t &bytes) {
    const std::string file = log + ".progress";
    if (FileUtil::fileExists(file.c_str()) == false) {
        return false;
    }
    const std::string content = waitForFile(file);
    if (content.compare(0, runId.size(), runId) != 0) {
        return false;
    }
    const char *columns[2];
    if (Util::getWordsOfLine((char *) content.c_str() + runId.size(), columns, 2) < 2) {
        Debug(Debug::ERROR) << "Cannot read " << file << "\n";
        EXIT(EXIT_FAILURE);
    }
    key = Util::fast_atoi<size_t>(columns[0]);
    bytes = Util::fast_atoi<size_t>(columns[1]);
    return true;
}

void MMseqsNode::waitForLog(const std::string &log, size_t &key, size_t &bytes) {
    for (unsigned int polls = 0; readLogProgress(log, key, bytes) == false; polls++) {
        sleepBeforePoll(polls);
    }
}

void MMseqsNode::readLog(const std::string &log, size_t from, size_t to, char *data) {
    FILE *file = FileUtil::openFileOrDie(log.c_str(), "r", true);
    FileUtil::preadAll(fileno(file), log.c_str(), data, to - from, from);
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << log << "\n";
        EXIT(EXIT_FAILURE);
    }
}

std::string MMseqsNode::barrier(const std::string &file, const std::string &content) {
    publishFile(file + "." + SSTR(rank), runId + content);
    if (isMaster()) {
        publishFile(file, runId + gatherRunFiles(file + "."));
    }
    return waitForRunFile(file);
}

size_t MMseqsNode::getSmallestMemoryLimit(const std::string &file, size_t memoryLimit, int threads) {
    const std::string nodeLimits = barrier(file, SSTR(memoryLimit) + ((threads > 0) ? "\t" + SSTR(threads) : "") + "\n");
    size_t smallestMemoryLimit = memoryLimit;
    for (char *line = (char *) nodeLimits.c_str(); *line != '\0'; line = Util::skipLine(line)) {
        const char *columns[2];
        Util::getWordsOfLine(line, columns, 2);
        if (threads > 0 && Util::fast_atoi<int>(columns[1]) != threads) {
            Debug(Debug::ERROR) << "Run all nodes with the same --threads\n";
            EXIT(EXIT_FAILURE);
        }
        smallestMemoryLimit = std::min(smallestMemoryLimit, Util::fast_atoi<size_t>(columns[0]));
    }
    return smallestMemoryLimit;
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
    // the buffers of the blocks and the copy of their records each take half of bufferMemory
    const size_t bufferSize = std::max(EXCHANGE_BUFFER_SIZE, bufferMemory / 2 / std::max<size_t>(1, blocks.size()));
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
            FileUtil::preadAll(fileno(block.file), prefix.c_str(), block.buffer.data() + block.bufferEnd, readBytes, block.pos);
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

MMseqsNode::ClusterPairWriter::ClusterPairWriter(const std::string &exchange, size_t pairCapacity, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets)
    : exchange(exchange), fileName(exchange + "." + SSTR(rank)), file(FileUtil::openAndDelete(fileName.c_str(), "w")), pairs(new(std::nothrow) ClusterPair[pairCapacity]),
      pairCapacity(pairCapacity), pairCount(0), compare(compare), nodeKeyOffsets(nodeKeyOffsets), writtenBytes(0) {
    Util::checkAllocation(pairs, "Can not allocate pairs memory in ClusterPairWriter");
}

MMseqsNode::ClusterPairWriter::~ClusterPairWriter() {
    delete[] pairs;
}

void MMseqsNode::ClusterPairWriter::flush() {
    writePairs(pairs, pairCount, compare, nodeKeyOffsets, file, fileName, indexData, writtenBytes);
    pairCount = 0;
}

std::string MMseqsNode::ClusterPairWriter::close() {
    if (pairCount > 0) {
        flush();
    }
    delete[] pairs;
    pairs = NULL;
    return closeExchange(file, fileName, exchange + ".index", indexData);
}

std::string MMseqsNode::sendClusterPairs(const std::string &clu, bool withRepresentative, bool (*compare)(const ClusterPair &, const ClusterPair &), const std::vector<size_t> &nodeKeyOffsets,
                                         const std::string &exchange, const std::string &receivedIndex, const std::string &receivedPrefix) {
    std::vector<ExchangeBlock> receivedBlocks;
    std::vector<FILE *> receivedFiles;
    if (receivedPrefix.empty() == false) {
        openExchange(receivedIndex, receivedPrefix, receivedBlocks, receivedFiles);
    }
    const size_t memoryLimit = Util::computeMemory(Parameters::getInstance().splitMemoryLimit);
    ClusterPairWriter pairWriter(exchange, memoryLimit / 2 / sizeof(ClusterPair), compare, nodeKeyOffsets);

    DBReader<DBKeyType> cluDb((clu + "." + SSTR(rank)).c_str(), (clu + ".index." + SSTR(rank)).c_str(), Parameters::getInstance().threads,
                              DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_LINFINITY);
    cluDb.open(DBReader<DBKeyType>::NOSORT);
    // a loaded and a prefetched split take a quarter of the memory, half for entry bytes and half for index cells; a larger entry is a split of its own
    const size_t splitBytes = memoryLimit / 16;
    const size_t splitEntries = std::max<size_t>(1, memoryLimit / 16 / (2 * sizeof(uint64_t) + 2 * sizeof(size_t)));
    const auto findSplitEnd = [&](size_t splitBegin, size_t splitOffset) {
        const size_t bytesEnd = (splitOffset + splitBytes >= cluDb.getDataSize()) ? cluDb.getSize() : std::max(splitBegin + 1, cluDb.findEntryByOffset(splitOffset + splitBytes + 1) - 1);
        return std::min(bytesEnd, splitBegin + splitEntries);
    };
    std::vector<DBKeyType> entryKeys;
    std::vector<DBKeyType> repKeys;
    std::vector<char> received;
    DBKeyType previousReceivedKey = DB_KEY_INVALID;
    size_t splitOffset = 0;
    for (size_t splitBegin = 0, splitEnd = findSplitEnd(0, 0), nextSplitEnd; splitBegin < cluDb.getSize(); splitBegin = splitEnd, splitEnd = nextSplitEnd) {
        cluDb.loadData(splitBegin, splitEnd);
        entryKeys.resize(splitEnd - splitBegin);
        repKeys.resize(splitEnd - splitBegin);
#pragma omp parallel for schedule(static) reduction(+:splitOffset)
        for (size_t i = splitBegin; i < splitEnd; i++) {
            entryKeys[i - splitBegin] = cluDb.getDbKey(i);
            repKeys[i - splitBegin] = entryKeys[i - splitBegin];
            splitOffset += cluDb.getEntryLen(i);
        }
        nextSplitEnd = findSplitEnd(splitEnd, splitOffset);
        if (splitEnd < cluDb.getSize()) {
            cluDb.prefetchData(splitEnd, nextSplitEnd);
        }
        const ClusterPair splitLast = { entryKeys.back(), entryKeys.back() };
        size_t receivedCount;
        while ((receivedCount = readExchangeRange(receivedBlocks, memoryLimit / 64, ClusterPair::recordSize, ClusterPair::compareByMember, splitLast, received, NULL, receivedPrefix)) > 0) {
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
        for (size_t firstEntry = splitBegin; firstEntry < splitEnd;) {
            // a line has at least a digit and a newline
            if (pairWriter.pairCount > 0 && pairWriter.pairCount + cluDb.getEntryLen(firstEntry) / 2 > pairWriter.pairCapacity) {
                pairWriter.flush();
            }
            // a cluster with more pairs than fit into memory is sent in several flushes
            if (cluDb.getEntryLen(firstEntry) / 2 > pairWriter.pairCapacity) {
                const DBKeyType entryKey = cluDb.getDbKey(firstEntry);
                char keyBuffer[255];
                for (char *data = cluDb.getData(firstEntry, 0); *data != '\0'; data = Util::skipLine(data)) {
                    Util::parseKey(data, keyBuffer);
                    const DBKeyType key = Util::fast_atoi<DBKeyType>(keyBuffer);
                    if (withRepresentative == false && key == entryKey) {
                        continue;
                    }
                    pairWriter.add(repKeys[firstEntry - splitBegin], key);
                }
                firstEntry++;
                continue;
            }
            size_t endEntry = firstEntry;
            for (size_t flushBytes = 0; endEntry < splitEnd && pairWriter.pairCount + (flushBytes + cluDb.getEntryLen(endEntry)) / 2 <= pairWriter.pairCapacity; endEntry++) {
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
                    size_t pairOffset = __sync_fetch_and_add(&pairWriter.pairCount, Util::countLines(data, cluDb.getEntryLen(i)) - (withRepresentative ? 0 : 1));
                    while (*data != '\0') {
                        Util::parseKey(data, keyBuffer);
                        const DBKeyType key = Util::fast_atoi<DBKeyType>(keyBuffer);
                        if (withRepresentative || key != entryKey) {
                            pairWriter.pairs[pairOffset].repKey = repKeys[i - splitBegin];
                            pairWriter.pairs[pairOffset].memberKey = key;
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
    cluDb.close();
    closeExchangeFiles(receivedFiles, receivedPrefix);
    return pairWriter.close();
}

void MMseqsNode::readKeyBits(int fd, const std::string &fileName, size_t fileOffset, size_t fileFirstKey, size_t from, size_t to, size_t firstKey, uint64_t *bits) {
    if (from >= to) {
        return;
    }
    const size_t firstWord = (from - fileFirstKey) / 64;
    std::vector<uint64_t> fileBits((to - 1 - fileFirstKey) / 64 + 1 - firstWord);
    FileUtil::preadAll(fd, fileName.c_str(), reinterpret_cast<char *>(fileBits.data()), fileBits.size() * sizeof(uint64_t), fileOffset + firstWord * sizeof(uint64_t));
    // a word of the file goes to bit (its first key - firstKey) of bits, without its keys outside [from, to)
    for (size_t word = 0; word < fileBits.size(); word++) {
        const size_t wordKey = fileFirstKey + (firstWord + word) * 64;
        uint64_t value = fileBits[word];
        value &= (wordKey < from) ? ~0ULL << (from - wordKey) : ~0ULL;
        value &= (to - wordKey < 64) ? (1ULL << (to - wordKey)) - 1 : ~0ULL;
        if (wordKey < firstKey) {
            value >>= firstKey - wordKey;
        }
        const size_t bit = (wordKey < firstKey) ? 0 : wordKey - firstKey;
        bits[bit / 64] |= value << (bit % 64);
        if (bit % 64 != 0) {
            bits[bit / 64 + 1] |= value >> (64 - bit % 64);
        }
    }
}

uint64_t *MMseqsNode::readRepresentativeBits(const std::string &db, const std::vector<size_t> &nodeKeyOffsets, size_t firstKey, size_t keyCount) {
    const size_t words = (keyCount + 63) / 64 + 1;
    uint64_t *bits = new(std::nothrow) uint64_t[words];
    Util::checkAllocation(bits, "Can not allocate representative bit memory");
    memset(bits, 0, words * sizeof(uint64_t));
    for (size_t node = 0; node + 1 < nodeKeyOffsets.size(); node++) {
        const size_t from = std::max(firstKey, nodeKeyOffsets[node]);
        const size_t to = std::min(firstKey + keyCount, nodeKeyOffsets[node + 1]);
        if (from >= to) {
            continue;
        }
        const std::string bitsFile = db + ".bits." + SSTR(node);
        FILE *file = FileUtil::openFileOrDie(bitsFile.c_str(), "r", true);
        uint64_t keyRange[2];
        FileUtil::preadAll(fileno(file), bitsFile.c_str(), reinterpret_cast<char *>(keyRange), sizeof(keyRange), 0);
        if (keyRange[0] > from || keyRange[0] + keyRange[1] < to) {
            Debug(Debug::ERROR) << "Cannot read " << bitsFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        readKeyBits(fileno(file), bitsFile, sizeof(keyRange), keyRange[0], from, to, firstKey, bits);
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << bitsFile << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    return bits;
}

void MMseqsNode::writeClusterIndexEntry(FILE *file, const std::string &fileName, DBKeyType key, size_t offset) {
    const uint64_t cell[2] = { static_cast<uint64_t>(key), static_cast<uint64_t>(offset) };
    if (fwrite(cell, sizeof(uint64_t), 2, file) != 2) {
        Debug(Debug::ERROR) << "Cannot write " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void MMseqsNode::writeClusterPairs(const std::string &index, const std::string &prefix, const std::string &clu, const uint64_t *singletonBits, size_t firstKey, size_t keyCount) {
    const Parameters &par = Parameters::getInstance();
    std::vector<ExchangeBlock> repBlocks;
    std::vector<FILE *> repFiles;
    openExchange(index, prefix, repBlocks, repFiles);
    const size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);
    const std::string resultDb = clu + "." + SSTR(rank);
    const std::string resultIndex = clu + ".index." + SSTR(rank);
    FILE *resultIndexFile = FileUtil::openAndDelete(resultIndex.c_str(), "w");
    DBWriter dbw(resultDb.c_str(), (resultDb + ".index").c_str(), 1, 0, Parameters::DBTYPE_OMIT_FILE);
    dbw.open();
    std::vector<char> records;
    std::vector<ClusterPair> singletons;
    std::vector<std::string> threadData(par.threads);
    // (key, offset) of the entries that start in the part of a thread
    std::vector<std::vector<uint64_t> > threadIndex(par.threads);
    std::vector<size_t> threadDataOffsets(par.threads + 1);
    DBKeyType previousRepKey = DB_KEY_INVALID;
    // the keys go in splits, so the pairs of their singletons fit into memory
    const size_t splitKeys = std::max<size_t>(1, memoryLimit / 16 / sizeof(ClusterPair));
    const size_t keyEnd = firstKey + keyCount;
    for (size_t splitBegin = firstKey; splitBegin < keyEnd; splitBegin += splitKeys) {
        const size_t splitEnd = std::min(keyEnd, splitBegin + splitKeys);
        const ClusterPair lastPair = { static_cast<DBKeyType>(splitEnd - 1), DB_KEY_INVALID };
        size_t nextSingletonKey = splitBegin;
        size_t pairCount;
        do {
            pairCount = readExchangeRange(repBlocks, memoryLimit / 8, ClusterPair::recordSize, ClusterPair::compareByRep, lastPair, records, NULL, prefix);
            singletons.clear();
            if (singletonBits != NULL) {
                // the singletons up to the largest rep of the pairs go with them in rep order, and the rest of the split after the last pairs
                size_t singletonEnd = (pairCount > 0) ? nextSingletonKey : splitEnd;
                for (size_t i = 0; i < pairCount; i++) {
                    singletonEnd = std::max<size_t>(singletonEnd, reinterpret_cast<ClusterPair *>(records.data())[i].repKey + 1);
                }
                const size_t bitBegin = nextSingletonKey - firstKey;
                const size_t bitEnd = singletonEnd - firstKey;
                for (size_t word = bitBegin / 64; word * 64 < bitEnd; word++) {
                    uint64_t bits = singletonBits[word] & ((word == bitBegin / 64) ? ~0ULL << (bitBegin % 64) : ~0ULL);
                    for (; bits != 0 && word * 64 + __builtin_ctzll(bits) < bitEnd; bits &= bits - 1) {
                        const DBKeyType key = firstKey + word * 64 + __builtin_ctzll(bits);
                        const ClusterPair singleton = { key, key };
                        singletons.push_back(singleton);
                    }
                }
                nextSingletonKey = std::max(nextSingletonKey, singletonEnd);
            }
            if (records.size() < (pairCount + singletons.size()) * sizeof(ClusterPair)) {
                records.resize((pairCount + singletons.size()) * sizeof(ClusterPair));
            }
            ClusterPair *pairs = reinterpret_cast<ClusterPair *>(records.data());
            std::copy(singletons.begin(), singletons.end(), pairs + pairCount);
            const size_t writeCount = pairCount + singletons.size();
            if (writeCount == 0) {
                break;
            }
            SORT_PARALLEL(pairs, pairs + writeCount, ClusterPair::compareByRep);
#pragma omp parallel for schedule(dynamic, 1)
            for (int thread = 0; thread < par.threads; thread++) {
                std::string &data = threadData[thread];
                std::vector<uint64_t> &threadEntries = threadIndex[thread];
                data.clear();
                threadEntries.clear();
                char buffer[32];
                const size_t begin = writeCount * thread / par.threads;
                DBKeyType repKey = (begin > 0) ? pairs[begin - 1].repKey : previousRepKey;
                for (size_t i = begin; i < writeCount * (thread + 1) / par.threads; i++) {
                    if (pairs[i].repKey != repKey) {
                        if (repKey != DB_KEY_INVALID) {
                            data.push_back('\0');
                        }
                        repKey = pairs[i].repKey;
                        threadEntries.push_back(static_cast<uint64_t>(repKey));
                        threadEntries.push_back(data.size());
                        char *tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(repKey), buffer);
                        *(tmpBuff - 1) = '\n';
                        data.append(buffer, tmpBuff - buffer);
                    }
                    if (pairs[i].memberKey != repKey) {
                        char *tmpBuff = Itoa::u64toa_sse2(static_cast<uint64_t>(pairs[i].memberKey), buffer);
                        *(tmpBuff - 1) = '\n';
                        data.append(buffer, tmpBuff - buffer);
                    }
                }
            }
            threadDataOffsets[0] = dbw.getOffset(0);
            for (int thread = 0; thread < par.threads; thread++) {
                threadDataOffsets[thread + 1] = threadDataOffsets[thread] + threadData[thread].size();
            }
#pragma omp parallel for schedule(dynamic, 1)
            for (int thread = 0; thread < par.threads; thread++) {
                for (size_t i = 1; i < threadIndex[thread].size(); i += 2) {
                    threadIndex[thread][i] += threadDataOffsets[thread];
                }
            }
            for (int thread = 0; thread < par.threads; thread++) {
                dbw.writeAdd(threadData[thread].data(), threadData[thread].size(), 0);
                if (fwrite(threadIndex[thread].data(), sizeof(uint64_t), threadIndex[thread].size(), resultIndexFile) != threadIndex[thread].size()) {
                    Debug(Debug::ERROR) << "Cannot write " << resultIndex << "\n";
                    EXIT(EXIT_FAILURE);
                }
            }
            previousRepKey = pairs[writeCount - 1].repKey;
        } while (pairCount > 0);
    }
    if (previousRepKey != DB_KEY_INVALID) {
        dbw.writeAdd("\0", 1, 0);
    }
    closeClusterWriter(dbw, resultIndexFile, resultIndex);
    closeExchangeFiles(repFiles, prefix);
}

void MMseqsNode::closeClusterWriter(DBWriter &writer, FILE *indexFile, const std::string &indexFileName) {
    writeClusterIndexEntry(indexFile, indexFileName, DB_KEY_INVALID, writer.getOffset(0));
    if (fclose(indexFile) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << indexFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    writer.close(false, false);
    FileUtil::remove(writer.getIndexFileName());
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
