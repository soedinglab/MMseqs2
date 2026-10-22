#include "DBWriter.h"
#include "Util.h"
#include "Parameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "ReducedMatrix.h"
#include "SubstitutionMatrix.h"
#include "DistanceCalculator.h"
#include "FastSort.h"
#include "MMseqsNode.h"
#include "itoa.h"
#include "FileUtil.h"

#include <limits>

#ifdef OPENMP
#include <omp.h>
#endif

struct __attribute__((__packed__)) HashEntry {
    size_t hash;
    DBLocalId id;

    static bool compareByHashAndId(const HashEntry &first, const HashEntry &second) {
        if (first.hash != second.hash) {
            return first.hash < second.hash;
        }
        return first.id < second.id;
    }

    static bool compareById(const HashEntry &first, const HashEntry &second) {
        return first.id < second.id;
    }
};

// keys of one length, starting at startKey and at offset bytes into the data file of the node, entries of entryLength bytes
struct LengthRange {
    size_t startKey;
    size_t offset;
    size_t entryLength;
};

static size_t hashWithLength(size_t hash, size_t length) {
    return hash ^ (length + static_cast<size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) + (hash >> 2));
}

// positions 64 * i .. 64 * i + 63 of a split length: which are members of the current hash range or block, and how many members come before them
struct MemberBits {
    uint64_t isMember;
    size_t membersBefore;
};

// index of the sequence at position pos of a split length among the members of the current hash range or block
static size_t memberIndex(const MemberBits *memberBits, size_t pos) {
    return memberBits[pos / 64].membersBefore + __builtin_popcountll(memberBits[pos / 64].isMember & ((1ull << (pos % 64)) - 1));
}

// a load fits when its sequence buffer is at most a third of the limit and its hashIdPairs fit next to two buffers
static bool fitsWithPrefetch(const LengthRange *begin, const LengthRange *end, bool isLinfinity, size_t seqBufferLimit, size_t hashIdPairsLimit) {
    const size_t seqBufferBytes = isLinfinity ? end->offset - begin->offset + 2 * Util::getPageSize() : 0;
    return seqBufferBytes <= seqBufferLimit && (end->startKey - begin->startKey) * sizeof(HashEntry) <= hashIdPairsLimit;
}

static size_t seqHash(Sequence &seq, size_t id, const char *data, size_t length) {
    seq.mapSequence(id, 0, data, length);
    return hashWithLength(Util::hash(seq.numSequence, seq.L), length);
}

// hashes each loaded sequence of ranges [0, rangeCount), ended by ranges[rangeCount].startKey, into hashIdPairs
static size_t fillHashIdPairs(DBReader<DBKeyType> &reader, const BaseMatrix *subMat, size_t maxSeqLen, const LengthRange *ranges, size_t rangeCount, HashEntry *hashIdPairs) {
    const size_t firstKey = ranges[0].startKey;
    const size_t seqCount = ranges[rangeCount].startKey - firstKey;
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence seq(maxSeqLen, reader.getDbtype(), subMat, 0, false, false);
        for (size_t range = 0; range < rangeCount; range++) {
            const size_t entryLength = ranges[range].entryLength;
            const char *rangeData = (entryLength > 0) ? reader.getDataFromLoadedBuffer(ranges[range].startKey) : NULL;
#pragma omp for schedule(dynamic, 100) nowait
            for (size_t key = ranges[range].startKey; key < ranges[range + 1].startKey; key++) {
                const char *data = (entryLength > 0) ? rangeData + (key - ranges[range].startKey) * entryLength : reader.getData(key, thread_idx);
                const size_t length = (entryLength > 0) ? entryLength - 2 : reader.getSeqLen(key);
                hashIdPairs[key - firstKey].hash = seqHash(seq, key, data, length);
                hashIdPairs[key - firstKey].id = key;
            }
        }
    }
    return seqCount;
}

// counts the sequences of each hash bucket of keys [firstKey, firstKey + seqCount), read readChunkSeqCount at a time,
// and cuts the buckets into hash ranges of at most maxRangeSeqCount sequences, like kmermatcher setupKmerSplits
static std::vector<std::pair<size_t, size_t>> setupHashRanges(DBReader<DBKeyType> &reader, const BaseMatrix *subMat, size_t maxSeqLen,
                                                              size_t firstKey, size_t seqCount, size_t readChunkSeqCount, size_t maxRangeSeqCount) {
    size_t *hashDist = new(std::nothrow) size_t[USHRT_MAX + 1];
    Util::checkAllocation(hashDist, "Can not allocate hashDist memory");
    memset(hashDist, 0, sizeof(size_t) * (USHRT_MAX + 1));
    for (size_t readChunkBegin = 0; readChunkBegin < seqCount; readChunkBegin += readChunkSeqCount) {
        const size_t readChunkEnd = std::min(readChunkBegin + readChunkSeqCount, seqCount);
        reader.loadData(firstKey + readChunkBegin, firstKey + readChunkEnd);
        if (readChunkEnd < seqCount) {
            reader.prefetchData(firstKey + readChunkEnd, firstKey + std::min(readChunkEnd + readChunkSeqCount, seqCount));
        }
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            Sequence seq(maxSeqLen, reader.getDbtype(), subMat, 0, false, false);
#pragma omp for schedule(dynamic, 100)
            for (size_t i = readChunkBegin; i < readChunkEnd; i++) {
                __sync_fetch_and_add(&hashDist[static_cast<unsigned short>(seqHash(seq, firstKey + i, reader.getData(firstKey + i, thread_idx), reader.getSeqLen(firstKey + i)))], 1);
            }
        }
    }

    std::vector<std::pair<size_t, size_t>> hashRanges;
    size_t rangeStart = 0;
    size_t rangeSeqCount = 0;
    for (size_t bucket = 0; bucket <= USHRT_MAX; bucket++) {
        if (hashDist[bucket] > maxRangeSeqCount) {
            Debug(Debug::ERROR) << "Not enough memory to hash " << hashDist[bucket] << " sequences of one hash bucket\n"
                                << "Raise RAM size or --split-memory-limit.\n";
            EXIT(EXIT_FAILURE);
        }
        if (rangeSeqCount + hashDist[bucket] > maxRangeSeqCount) {
            hashRanges.emplace_back(rangeStart, bucket);
            rangeStart = bucket;
            rangeSeqCount = 0;
        }
        rangeSeqCount += hashDist[bucket];
    }
    hashRanges.emplace_back(rangeStart, USHRT_MAX + 1);
    delete[] hashDist;
    return hashRanges;
}

// keeps the sequences of keys [firstKey, firstKey + seqCount) whose hash bucket is in hashRange
static size_t fillHashRange(DBReader<DBKeyType> &reader, const BaseMatrix *subMat, size_t maxSeqLen, size_t firstKey, size_t seqCount, size_t readChunkSeqCount,
                            const std::pair<size_t, size_t> &hashRange, HashEntry *hashIdPairs, MemberBits *memberBits) {
    size_t rangeSeqCount = 0;
    for (size_t readChunkBegin = 0; readChunkBegin < seqCount; readChunkBegin += readChunkSeqCount) {
        const size_t readChunkEnd = std::min(readChunkBegin + readChunkSeqCount, seqCount);
        reader.loadData(firstKey + readChunkBegin, firstKey + readChunkEnd);
        if (readChunkEnd < seqCount) {
            reader.prefetchData(firstKey + readChunkEnd, firstKey + std::min(readChunkEnd + readChunkSeqCount, seqCount));
        }
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            Sequence seq(maxSeqLen, reader.getDbtype(), subMat, 0, false, false);
            // a thread marks whole groups of 64 positions, so read chunks start at a multiple of 64
#pragma omp for schedule(dynamic, 1)
            for (size_t i = readChunkBegin / 64; i < (readChunkEnd + 63) / 64; i++) {
                memberBits[i].isMember = 0;
                for (size_t pos = i * 64; pos < std::min(i * 64 + 64, readChunkEnd); pos++) {
                    const unsigned short bucket = static_cast<unsigned short>(seqHash(seq, firstKey + pos, reader.getData(firstKey + pos, thread_idx), reader.getSeqLen(firstKey + pos)));
                    if (bucket >= hashRange.first && bucket < hashRange.second) {
                        memberBits[i].isMember |= 1ull << (pos % 64);
                    }
                }
            }
#pragma omp single
            for (size_t i = readChunkBegin / 64; i < (readChunkEnd + 63) / 64; i++) {
                memberBits[i].membersBefore = rangeSeqCount;
                rangeSeqCount += __builtin_popcountll(memberBits[i].isMember);
            }
#pragma omp for schedule(dynamic, 1)
            for (size_t i = readChunkBegin / 64; i < (readChunkEnd + 63) / 64; i++) {
                size_t slot = memberBits[i].membersBefore;
                for (uint64_t isMember = memberBits[i].isMember; isMember != 0; isMember &= isMember - 1) {
                    const size_t id = firstKey + i * 64 + __builtin_ctzll(isMember);
                    hashIdPairs[slot].hash = seqHash(seq, id, reader.getData(id, thread_idx), reader.getSeqLen(id));
                    hashIdPairs[slot].id = id;
                    slot++;
                }
            }
        }
    }
    return rangeSeqCount;
}

// hashes and sorts keys [firstKey, firstKey + seqCount) of one length, read readChunkSeqCount at a time; hash bins of one sequence become representatives,
// the other bins are moved to the front of hashIdPairs and their count is returned
static size_t hashAndSortLength(DBReader<DBKeyType> &reader, const BaseMatrix *subMat, size_t maxSeqLen, size_t firstKey, size_t seqCount, size_t readChunkSeqCount,
                                HashEntry *hashIdPairs, uint64_t *representativeBits, size_t localFirstKey) {
    const size_t length = reader.getSeqLen(firstKey);
    for (size_t readChunkBegin = 0; readChunkBegin < seqCount; readChunkBegin += readChunkSeqCount) {
        const size_t readChunkEnd = std::min(readChunkBegin + readChunkSeqCount, seqCount);
        reader.loadData(firstKey + readChunkBegin, firstKey + readChunkEnd);
        if (readChunkEnd < seqCount) {
            reader.prefetchData(firstKey + readChunkEnd, firstKey + std::min(readChunkEnd + readChunkSeqCount, seqCount));
        }
        const LengthRange readChunkRanges[2] = { { firstKey + readChunkBegin, 0, length + 2 }, { firstKey + readChunkEnd, 0, 0 } };
        fillHashIdPairs(reader, subMat, maxSeqLen, readChunkRanges, 1, hashIdPairs + readChunkBegin);
    }
    SORT_PARALLEL(hashIdPairs, hashIdPairs + seqCount, HashEntry::compareByHashAndId);
    size_t writePos = 0;
    for (size_t hashBinBegin = 0; hashBinBegin < seqCount;) {
        size_t hashBinEnd = hashBinBegin + 1;
        while (hashBinEnd < seqCount && hashIdPairs[hashBinEnd].hash == hashIdPairs[hashBinBegin].hash) {
            hashBinEnd++;
        }
        if (hashBinEnd - hashBinBegin == 1) {
            const size_t bit = hashIdPairs[hashBinBegin].id - localFirstKey;
            representativeBits[bit / 64] |= 1ull << (bit % 64);
        } else {
            for (size_t i = hashBinBegin; i < hashBinEnd; i++) {
                hashIdPairs[writePos++] = hashIdPairs[i];
            }
        }
        hashBinBegin = hashBinEnd;
    }
    return writePos;
}

// copies the sequences of blockHashIdPairs, keys of [firstKey, firstKey + seqCount), to their slots in blockSeqs, read readChunkSeqCount at a time
static void loadBlockSeqs(DBReader<DBKeyType> &reader, size_t firstKey, size_t seqCount, size_t readChunkSeqCount, const HashEntry *blockHashIdPairs, size_t blockSeqCount,
                          MemberBits *memberBits, char *blockSeqs) {
    const size_t length = reader.getSeqLen(firstKey);
    memset(memberBits, 0, sizeof(MemberBits) * (seqCount / 64 + 1));
#pragma omp parallel for
    for (size_t i = 0; i < blockSeqCount; i++) {
        const size_t pos = blockHashIdPairs[i].id - firstKey;
        __sync_fetch_and_or(&memberBits[pos / 64].isMember, 1ull << (pos % 64));
    }
    size_t memberCount = 0;
    for (size_t i = 0; i < seqCount / 64 + 1; i++) {
        memberBits[i].membersBefore = memberCount;
        memberCount += __builtin_popcountll(memberBits[i].isMember);
    }
    for (size_t readChunkBegin = 0; readChunkBegin < seqCount; readChunkBegin += readChunkSeqCount) {
        const size_t readChunkEnd = std::min(readChunkBegin + readChunkSeqCount, seqCount);
        reader.loadData(firstKey + readChunkBegin, firstKey + readChunkEnd);
        if (readChunkEnd < seqCount) {
            reader.prefetchData(firstKey + readChunkEnd, firstKey + std::min(readChunkEnd + readChunkSeqCount, seqCount));
        }
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(dynamic, 1)
            for (size_t i = readChunkBegin / 64; i < (readChunkEnd + 63) / 64; i++) {
                size_t member = memberBits[i].membersBefore;
                for (uint64_t isMember = memberBits[i].isMember; isMember != 0; isMember &= isMember - 1) {
                    const size_t id = firstKey + i * 64 + __builtin_ctzll(isMember);
                    memcpy(blockSeqs + member * length, reader.getData(id, thread_idx), length);
                    member++;
                }
            }
        }
    }
}

int clusthashfast(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE, 5));
    par.seqIdThr = (float)Parameters::CLUST_HASH_DEFAULT_MIN_SEQ_ID / 100.0f;
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<DBKeyType> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_LINFINITY);
    reader.open(DBReader<DBKeyType>::NOSORT);
    if (Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_AMINO_ACIDS) == false) {
        Debug(Debug::ERROR) << "clusthashfast only supports amino acid databases\n";
        EXIT(EXIT_FAILURE);
    }
    // a compressed or padded entry is decoded into one buffer per thread, so the representative is copied before its members are read
    const bool copyRepSeq = reader.isCompressed() || (DBReader<DBKeyType>::getExtendedDbtype(reader.getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU);
    const bool isLinfinity = DBReader<DBKeyType>::getExtendedDbtype(reader.getDbtype()) & Parameters::DBTYPE_EXTENDED_LINFINITY;
    if (isLinfinity == false && par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        reader.readMmapedDataInMemory();
    }

    std::vector<LengthRange> nodeSeqRange;
    if (isLinfinity) {
        MMseqsNode::init();
        if (reader.getDataFileCnt() != MMseqsNode::count) {
            Debug(Debug::ERROR) << "Database " << par.db1 << " has " << reader.getDataFileCnt() << " data files, but clusthashfast runs on " << MMseqsNode::count << " nodes\n";
            EXIT(EXIT_FAILURE);
        }
        const DBReader<DBKeyType>::Index *index = reader.getIndex();
        size_t endKey = reader.getSize();
        size_t endOffset = 0;
        for (size_t range = 0; range < reader.getIndexSize(); range++) {
            if (reader.getIndexFile(range) == MMseqsNode::rank) {
                nodeSeqRange.push_back({index[range].id, index[range].offset, index[range].length});
                endKey = (range + 1 < reader.getIndexSize()) ? index[range + 1].id : reader.getSize();
                endOffset = index[range].offset + (endKey - index[range].id) * index[range].length;
            }
        }
        nodeSeqRange.push_back({endKey, endOffset, 0});
    } else {
        nodeSeqRange.push_back({0, 0, 0});
        nodeSeqRange.push_back({reader.getSize(), reader.getDataSize(), 0});
    }
    const size_t localFirstKey = nodeSeqRange.front().startKey;
    const size_t localKeySize = nodeSeqRange.back().startKey - localFirstKey;

    if (isLinfinity && par.compressed) {
        Debug(Debug::WARNING) << "Compressed output cannot be combined with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
        Debug(Debug::WARNING) << "We recompute with --compressed 0\n";
        par.compressed = 0;
    }
    const std::string resultDb = isLinfinity ? par.db2 + "." + SSTR(MMseqsNode::rank) : par.db2;
    DBWriter resultWriter(resultDb.c_str(), (resultDb + ".index").c_str(), 1, par.compressed, isLinfinity ? Parameters::DBTYPE_OMIT_FILE : Parameters::DBTYPE_CLUSTER_RES);
    resultWriter.open();
    const std::string resultIndex = par.db2 + ".index." + SSTR(MMseqsNode::rank);
    FILE *resultIndexFile = isLinfinity ? FileUtil::openAndDelete(resultIndex.c_str(), "w") : NULL;
    const std::string representativeDb = par.db2 + "_redundancy";
    DBWriter representativeWriter(representativeDb.c_str(), (representativeDb + ".index").c_str(), 1, 0, Parameters::DBTYPE_OMIT_FILE);
    if (isLinfinity == false) {
        representativeWriter.open();
    }
    const size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);
    // a mode 3 node keeps its representatives as one bit per key
    const size_t representativeBitsSize = isLinfinity ? (localKeySize + 63) / 64 : 0;
    uint64_t *representativeBits = new(std::nothrow) uint64_t[representativeBitsSize]();
    Util::checkAllocation(representativeBits, "Can not allocate representative bit memory");
    if (representativeBitsSize * sizeof(uint64_t) >= memoryLimit) {
        Debug(Debug::ERROR) << "Not enough memory to run clusthashfast. Memory limit: " << memoryLimit
                            << " bytes, the representative bits need: " << representativeBitsSize * sizeof(uint64_t) << " bytes.\n"
                            << "Raise RAM size or --split-memory-limit.\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t splitMemoryLimit = static_cast<size_t>(static_cast<double>(memoryLimit - representativeBitsSize * sizeof(uint64_t)) * 0.95);

    // consecutive lengths that fit together are loaded at once;
    // a length that does not fit is loaded alone: a mode 3 length is read one read chunk at a time to hash it and once per block of candidates, a standard DB length once per hash range that fits, after 2 bits per key
    const size_t readChunkBytes = std::min(MAX_INFLIGHT_READ_SIZE, splitMemoryLimit / 4);
    const size_t seqBufferLimit = isLinfinity ? splitMemoryLimit / 3 : 0;
    const size_t hashIdPairsLimit = splitMemoryLimit - 2 * seqBufferLimit;
    std::vector<size_t> hashIdPairsCapacity(nodeSeqRange.size() - 1, 0);
    std::vector<size_t> loadOffsets;
    for (size_t i = 0; i + 1 < nodeSeqRange.size(); i++) {
        const size_t lengthSeqCount = nodeSeqRange[i + 1].startKey - nodeSeqRange[i].startKey;
        const size_t lengthBytes = (isLinfinity ? nodeSeqRange[i + 1].offset - nodeSeqRange[i].offset : 0) + lengthSeqCount * sizeof(HashEntry);
        if (lengthBytes > splitMemoryLimit) {
            const size_t readChunkSeqCount = isLinfinity ? std::min(readChunkBytes / nodeSeqRange[i].entryLength / 64 * 64, lengthSeqCount) : 0;
            const size_t readChunkBufferBytes = 2 * reader.getBufferSizeToLoad(nodeSeqRange[i].startKey, nodeSeqRange[i].startKey + readChunkSeqCount);
            const size_t keptBytes = isLinfinity ? (lengthSeqCount / 64 + 1) * sizeof(MemberBits) : lengthSeqCount / 4;
            const size_t hashBytes = sizeof(HashEntry) + (isLinfinity ? nodeSeqRange[i].entryLength - 2 : 0);
            // a split mode 3 length keeps a pair for every member to write its clusters in key order
            const size_t pairBytes = isLinfinity ? lengthSeqCount * sizeof(ClusterPair) : 0;
            if ((isLinfinity && readChunkSeqCount == 0) || splitMemoryLimit <= keptBytes + readChunkBufferBytes + pairBytes + hashBytes) {
                Debug(Debug::ERROR) << "Not enough memory to hash " << lengthSeqCount << " sequences\n"
                                    << "Raise RAM size or --split-memory-limit.\n";
                EXIT(EXIT_FAILURE);
            }
            // a split mode 3 length keeps its hashIdPairs and, after its candidates, the sequences of a block in hashIdPairs
            hashIdPairsCapacity[i] = (splitMemoryLimit - keptBytes - readChunkBufferBytes - pairBytes) / sizeof(HashEntry);
        }
        if (loadOffsets.empty() || fitsWithPrefetch(&nodeSeqRange[loadOffsets.back()], &nodeSeqRange[i + 1], isLinfinity, seqBufferLimit, hashIdPairsLimit) == false) {
            loadOffsets.push_back(i);
        }
    }
    loadOffsets.push_back(nodeSeqRange.size() - 1);
    HashEntry *hashIdPairs = NULL;
    size_t hashIdPairsSize = 0;

    std::string result;
    char buffer[32];
    ClusterPair *loadPairs = NULL;
    size_t loadPairCount = 0;
    bool needsSort = false;
    size_t clusters = 0;
    SubstitutionMatrix sMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, -0.2);
    ReducedMatrix subMat(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2num, sMat.num2aa, sMat.alphabetSize, par.alphabetSize.values.aminoacid(), 2.0);

    Debug(Debug::INFO) << "Hashing " << localKeySize << " sequences in " << loadOffsets.size() - 1 << " loads\n";
    Debug::Progress progress(loadOffsets.size() - 1);
    for (size_t load = 0; load + 1 < loadOffsets.size(); load++) {
        progress.updateProgress();
        const size_t loadFirstKey = nodeSeqRange[loadOffsets[load]].startKey;
        const size_t loadSeqCount = nodeSeqRange[loadOffsets[load + 1]].startKey - loadFirstKey;
        const bool isSplitLength = hashIdPairsCapacity[loadOffsets[load]] > 0;
        const size_t length = isLinfinity ? nodeSeqRange[loadOffsets[load]].entryLength - 2 : 0;
        const bool copiesSeqs = isSplitLength && isLinfinity;
        if (copiesSeqs) {
            loadPairs = new(std::nothrow) ClusterPair[loadSeqCount];
            Util::checkAllocation(loadPairs, "Can not allocate loadPairs memory in clusthashfast");
            loadPairCount = 0;
        }
        const size_t readChunkSeqCount = copiesSeqs ? std::min(readChunkBytes / (length + 2) / 64 * 64, loadSeqCount) : loadSeqCount;
        const size_t hashedSeqCapacity = isSplitLength ? hashIdPairsCapacity[loadOffsets[load]] : loadSeqCount;
        bool hashIdPairsFit = hashedSeqCapacity >= loadSeqCount;
        const bool prefetchesNextLoad = load + 2 < loadOffsets.size()
                                        && fitsWithPrefetch(&nodeSeqRange[loadOffsets[load]], &nodeSeqRange[loadOffsets[load + 1]], isLinfinity, seqBufferLimit, hashIdPairsLimit)
                                        && fitsWithPrefetch(&nodeSeqRange[loadOffsets[load + 1]], &nodeSeqRange[loadOffsets[load + 2]], isLinfinity, seqBufferLimit, hashIdPairsLimit);
        if (hashIdPairsSize < hashedSeqCapacity) {
            delete[] hashIdPairs;
            hashIdPairs = new(std::nothrow) HashEntry[hashedSeqCapacity];
            Util::checkAllocation(hashIdPairs, "Can not allocate hash entry memory");
            hashIdPairsSize = hashedSeqCapacity;
        }

        // [first, second) of hash buckets of a split length, or of hashIdPairs (blocks) of a split mode 3 length whose hashIdPairs fit
        std::vector<std::pair<size_t, size_t>> hashRanges(1, std::make_pair(0, USHRT_MAX + 1));
        MemberBits *memberBits = NULL;
        char *blockSeqs = NULL;
        if (isSplitLength) {
            memberBits = new(std::nothrow) MemberBits[loadSeqCount / 64 + 1];
            Util::checkAllocation(memberBits, "Can not allocate member bit memory");
        }
        if (copiesSeqs && hashIdPairsFit) {
            const size_t candidateCount = hashAndSortLength(reader, &subMat, par.maxSeqLen, loadFirstKey, loadSeqCount, readChunkSeqCount, hashIdPairs, representativeBits, localFirstKey);
            // the candidates are cut at hash bins into blocks whose copies fit
            const size_t copyCapacity = (hashedSeqCapacity - candidateCount) * sizeof(HashEntry) / length;
            hashRanges.clear();
            size_t blockBegin = 0;
            for (size_t hashBinBegin = 0; hashBinBegin < candidateCount;) {
                size_t hashBinEnd = hashBinBegin + 1;
                while (hashBinEnd < candidateCount && hashIdPairs[hashBinEnd].hash == hashIdPairs[hashBinBegin].hash) {
                    hashBinEnd++;
                }
                // the copies of one hash bin do not fit next to all hashIdPairs, so the length is hashed once per hash range
                if (hashBinEnd - hashBinBegin > copyCapacity) {
                    hashIdPairsFit = false;
                    break;
                }
                if (hashBinEnd - blockBegin > copyCapacity) {
                    hashRanges.emplace_back(blockBegin, hashBinBegin);
                    blockBegin = hashBinBegin;
                }
                hashBinBegin = hashBinEnd;
            }
            if (candidateCount > blockBegin) {
                hashRanges.emplace_back(blockBegin, candidateCount);
            }
            blockSeqs = reinterpret_cast<char *>(hashIdPairs + candidateCount);
        }
        if (isSplitLength && hashIdPairsFit == false) {
            // the sequences of a hash range of a mode 3 length are copied after its hashIdPairs
            const size_t maxRangeSeqCount = hashedSeqCapacity * sizeof(HashEntry) / (sizeof(HashEntry) + length);
            hashRanges = setupHashRanges(reader, &subMat, par.maxSeqLen, loadFirstKey, loadSeqCount, readChunkSeqCount, maxRangeSeqCount);
            blockSeqs = copiesSeqs ? reinterpret_cast<char *>(hashIdPairs + maxRangeSeqCount) : NULL;
        } else if (isSplitLength == false) {
            reader.loadData(loadFirstKey, loadFirstKey + loadSeqCount);
            if (prefetchesNextLoad) {
                reader.prefetchData(nodeSeqRange[loadOffsets[load + 1]].startKey, nodeSeqRange[loadOffsets[load + 2]].startKey);
            }
        }
        needsSort |= hashRanges.size() > 1;
        for (size_t range = 0; range < hashRanges.size(); range++) {
            HashEntry *rangeHashIdPairs = hashIdPairs;
            size_t hashedSeqCount = 0;
            if (copiesSeqs && hashIdPairsFit) {
                // hashIdPairs of a block are sorted by hash and id already
                rangeHashIdPairs = hashIdPairs + hashRanges[range].first;
                hashedSeqCount = hashRanges[range].second - hashRanges[range].first;
            } else {
                hashedSeqCount = isSplitLength
                               ? fillHashRange(reader, &subMat, par.maxSeqLen, loadFirstKey, loadSeqCount, readChunkSeqCount, hashRanges[range], hashIdPairs, memberBits)
                               : fillHashIdPairs(reader, &subMat, par.maxSeqLen, nodeSeqRange.data() + loadOffsets[load], loadOffsets[load + 1] - loadOffsets[load], hashIdPairs);
                SORT_PARALLEL(hashIdPairs, hashIdPairs + hashedSeqCount, HashEntry::compareByHashAndId);
            }
            if (copiesSeqs) {
                loadBlockSeqs(reader, loadFirstKey, loadSeqCount, readChunkSeqCount, rangeHashIdPairs, hashedSeqCount, memberBits, blockSeqs);
            }

            // bins never cross threads; threadOffsets are fixed before hash is overwritten with the representative id
            std::vector<size_t> threadOffsets(par.threads + 1, hashedSeqCount);
            for (int thread = 0; thread < par.threads; thread++) {
                const size_t pos = hashedSeqCount / par.threads * thread;
                threadOffsets[thread] = (pos == 0) ? 0 : std::upper_bound(rangeHashIdPairs + pos, rangeHashIdPairs + hashedSeqCount, HashEntry{rangeHashIdPairs[pos - 1].hash, std::numeric_limits<DBLocalId>::max()}, HashEntry::compareByHashAndId) - rangeHashIdPairs;
            }
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (int thread = 0; thread < par.threads; thread++) {
                const size_t threadEnd = threadOffsets[thread + 1];
                std::string repSeqCopy;
                for (size_t hashBinBegin = threadOffsets[thread]; hashBinBegin < threadEnd;) {
                    size_t hashBinEnd = hashBinBegin + 1;
                    while (hashBinEnd < threadEnd && rangeHashIdPairs[hashBinEnd].hash == rangeHashIdPairs[hashBinBegin].hash) {
                        hashBinEnd++;
                    }
                    // [hashBinBegin, i) is assigned, [i, hashBinEnd) is not: the smallest id left is the next representative, members are swapped in front
                    for (size_t i = hashBinBegin; i < hashBinEnd;) {
                        std::iter_swap(rangeHashIdPairs + i, std::min_element(rangeHashIdPairs + i, rangeHashIdPairs + hashBinEnd, HashEntry::compareById));
                        const DBLocalId repId = rangeHashIdPairs[i].id;
                        rangeHashIdPairs[i++].hash = repId;
                        if (i == hashBinEnd) {
                            break;
                        }
                        const unsigned int repLength = reader.getSeqLen(repId);
                        const char *repSeq = (blockSeqs == NULL) ? reader.getData(repId, thread) : blockSeqs + memberIndex(memberBits, repId - loadFirstKey) * length;
                        if (copyRepSeq) {
                            repSeqCopy.assign(repSeq, repLength);
                            repSeq = repSeqCopy.data();
                        }
                        for (size_t j = i; j < hashBinEnd; j++) {
                            const DBLocalId memberId = rangeHashIdPairs[j].id;
                            if (reader.getSeqLen(memberId) != repLength) {
                                continue;
                            }
                            const unsigned int distance = DistanceCalculator::computeInverseHammingDistance(repSeq, (blockSeqs == NULL) ? reader.getData(memberId, thread) : blockSeqs + memberIndex(memberBits, memberId - loadFirstKey) * length, repLength);
                            if (static_cast<float>(distance) / static_cast<float>(repLength) >= par.seqIdThr) {
                                std::swap(rangeHashIdPairs[i], rangeHashIdPairs[j]);
                                rangeHashIdPairs[i++].hash = repId;
                            }
                        }
                    }
                    hashBinBegin = hashBinEnd;
                }
            }

            // hash now holds the representative id, so sorting again gathers each cluster with its representative first
            SORT_PARALLEL(rangeHashIdPairs, rangeHashIdPairs + hashedSeqCount, HashEntry::compareByHashAndId);
            for (size_t clusterBegin = 0; clusterBegin < hashedSeqCount;) {
                size_t clusterEnd = clusterBegin + 1;
                while (clusterEnd < hashedSeqCount && rangeHashIdPairs[clusterEnd].hash == rangeHashIdPairs[clusterBegin].hash) {
                    clusterEnd++;
                }
                const DBLocalId repId = rangeHashIdPairs[clusterBegin].id;
                const DBKeyType repKey = reader.getDbKey(repId);
                if (isLinfinity) {
                    representativeBits[(repId - localFirstKey) / 64] |= 1ull << ((repId - localFirstKey) % 64);
                } else {
                    representativeWriter.writeIndexEntry(repKey, reader.getOffset(repId), reader.getEntryLen(repId), 0);
                    clusters++;
                }
                // the clusters of a split mode 3 length are written in key order after its last hash range
                if (copiesSeqs) {
                    for (size_t i = clusterBegin + 1; i < clusterEnd; i++) {
                        const ClusterPair pair = { repKey, reader.getDbKey(rangeHashIdPairs[i].id) };
                        loadPairs[loadPairCount++] = pair;
                    }
                } else {
                    for (size_t i = clusterBegin; i < clusterEnd; i++) {
                        char *end = Itoa::u64toa_sse2(static_cast<uint64_t>(reader.getDbKey(rangeHashIdPairs[i].id)), buffer);
                        *(end - 1) = '\n';
                        result.append(buffer, end - buffer);
                    }
                    if (isLinfinity) {
                        MMseqsNode::writeClusterIndexEntry(resultIndexFile, resultIndex, repKey, resultWriter.getOffset(0));
                    }
                    resultWriter.writeData(result.c_str(), result.length(), repKey, 0, true, isLinfinity == false);
                    result.clear();
                }
                clusterBegin = clusterEnd;
            }
        }
        if (copiesSeqs) {
            // every representative of the length starts an entry, also one that is alone in its hash
            SORT_PARALLEL(loadPairs, loadPairs + loadPairCount, ClusterPair::compareByRep);
            size_t pair = 0;
            for (size_t key = loadFirstKey; key < loadFirstKey + loadSeqCount; key++) {
                if ((representativeBits[(key - localFirstKey) / 64] & (1ull << ((key - localFirstKey) % 64))) == 0) {
                    continue;
                }
                MMseqsNode::writeClusterIndexEntry(resultIndexFile, resultIndex, key, resultWriter.getOffset(0));
                resultWriter.writeStart(0);
                char *end = Itoa::u64toa_sse2(static_cast<uint64_t>(key), buffer);
                *(end - 1) = '\n';
                resultWriter.writeAdd(buffer, end - buffer, 0);
                for (; pair < loadPairCount && loadPairs[pair].repKey == key; pair++) {
                    end = Itoa::u64toa_sse2(static_cast<uint64_t>(loadPairs[pair].memberKey), buffer);
                    *(end - 1) = '\n';
                    resultWriter.writeAdd(buffer, end - buffer, 0);
                }
                resultWriter.writeEnd(key, 0, true, false);
            }
            delete[] loadPairs;
            loadPairs = NULL;
        }
        delete[] memberBits;
        // only a load and the next one it prefetched keep both sequence buffers and hashIdPairs
        if (isLinfinity && prefetchesNextLoad == false) {
            reader.releaseDataPages();
            delete[] hashIdPairs;
            hashIdPairs = NULL;
            hashIdPairsSize = 0;
        }
    }

    if (isLinfinity) {
        MMseqsNode::closeClusterWriter(resultWriter, resultIndexFile, resultIndex);
        // the bits are on disk before the index part that tells the master this node is done
        const std::string bitsFile = representativeDb + ".bits." + SSTR(MMseqsNode::rank);
        FILE *bits = FileUtil::openAndDelete(bitsFile.c_str(), "w");
        const uint64_t header[2] = { localFirstKey, localKeySize };
        if (fwrite(header, sizeof(uint64_t), 2, bits) != 2 || fwrite(representativeBits, sizeof(uint64_t), representativeBitsSize, bits) != representativeBitsSize
            || fflush(bits) != 0 || fsync(fileno(bits)) != 0 || fclose(bits) != 0) {
            Debug(Debug::ERROR) << "Cannot write " << bitsFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        // the ranges of this node in the sequence index, with the number of its representatives as a sixth column
        std::string representativeIndex;
        for (size_t i = 0; i + 1 < nodeSeqRange.size(); i++) {
            const size_t rangeEndKey = nodeSeqRange[i + 1].startKey;
            const size_t rangeRepresentatives = Util::countSetBits(representativeBits, nodeSeqRange[i].startKey - localFirstKey, rangeEndKey - localFirstKey);
            representativeIndex += SSTR(nodeSeqRange[i].startKey) + "\t" + SSTR(nodeSeqRange[i].offset) + "\t" + SSTR(nodeSeqRange[i].entryLength) + "\t"
                                 + SSTR(rangeEndKey - nodeSeqRange[i].startKey) + "\t" + SSTR(MMseqsNode::rank) + "\t" + SSTR(rangeRepresentatives) + "\n";
        }
        MMseqsNode::publishRunFile(representativeDb + ".index." + SSTR(MMseqsNode::rank), representativeIndex);
        if (MMseqsNode::isMaster()) {
            MMseqsNode::publishFile(representativeDb + ".index", MMseqsNode::gatherRunFiles(representativeDb + ".index."));
        }
    } else {
        // clusters are written in key order unless a length was split into hash ranges
        resultWriter.close(false, needsSort);
        representativeWriter.close(false, needsSort);
    }
    if (MMseqsNode::isMaster()) {
        DBReader<DBKeyType>::softlinkDb(par.db1, representativeDb, DBFiles::DATA);
        DBWriter::writeDbtypeFile(representativeDb.c_str(), reader.getDbtype(), reader.isCompressed());
        if (isLinfinity) {
            DBReader<DBKeyType> representativeReader(representativeDb.c_str(), (representativeDb + ".index").c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_LINFINITY);
            representativeReader.open(DBReader<DBKeyType>::NOSORT);
            for (size_t range = 0; range < representativeReader.getIndexSize(); range++) {
                clusters += representativeReader.getRangeRepresentativeCount(range);
            }
            representativeReader.close();
            // a mode 3 subset only feeds the next pipeline steps, which read no headers, lookup or source
            DBWriter::writeDbtypeFile(par.db2.c_str(), DBReader<DBKeyType>::setExtendedDbtype(Parameters::DBTYPE_CLUSTER_RES, Parameters::DBTYPE_EXTENDED_LINFINITY), false);
        } else {
            DBReader<DBKeyType>::softlinkDb(par.db1, representativeDb, DBFiles::SEQUENCE_ANCILLARY);
        }
    }
    Debug(Debug::INFO) << "Size of the sequence database: " << reader.getSize() << "\n";
    Debug(Debug::INFO) << "Number of clusters: " << clusters << "\n";
    delete[] representativeBits;
    delete[] hashIdPairs;
    reader.close();
    return EXIT_SUCCESS;
}
