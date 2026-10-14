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
#include "MemoryMapped.h"
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

// keys of one length, starting at startKey and at offset bytes into the data
struct LengthStart {
    size_t startKey;
    size_t offset;
};

static size_t hashWithLength(size_t hash, size_t length) {
    return hash ^ (length + static_cast<size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) + (hash >> 2));
}

// positions 64 * i .. 64 * i + 63 of a length split into hash resultIndex: which belong to the current range, and how many copies come before them
struct HashRangeMembers {
    uint64_t inRange;
    size_t copiesBefore;
};

// slot of the copy of the sequence at position pos of a length, in the current hash range
static size_t copySlot(const HashRangeMembers *hashRangeMembers, size_t pos) {
    return hashRangeMembers[pos / 64].copiesBefore + __builtin_popcountll(hashRangeMembers[pos / 64].inRange & ((1ull << (pos % 64)) - 1));
}

static size_t seqHash(DBReader<DBKeyType> &reader, Sequence &seq, size_t id, unsigned int thread_idx) {
    const size_t length = reader.getSeqLen(id);
    seq.mapSequence(id, 0, reader.getData(id, thread_idx), length);
    return hashWithLength(Util::hash(seq.numSequence, seq.L), length);
}

// loads keys [firstKey, firstKey + seqCount) and hashes each sequence into hashIdPairs
static size_t fillHashIdPairs(DBReader<DBKeyType> &reader, const BaseMatrix *subMat, size_t maxSeqLen, size_t firstKey, size_t seqCount, HashEntry *hashIdPairs) {
    reader.loadData(firstKey, seqCount);
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence seq(maxSeqLen, reader.getDbtype(), subMat, 0, false, false);
#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < seqCount; i++) {
            hashIdPairs[i].hash = seqHash(reader, seq, firstKey + i, thread_idx);
            hashIdPairs[i].id = firstKey + i;
        }
    }
    return seqCount;
}

// counts the sequences of each hash bucket of keys [firstKey, firstKey + seqCount), read chunkSeqCount at a time,
// and cuts the buckets into hash resultIndex of at most maxRangeSeqCount sequences, like kmermatcher setupKmerSplits
static std::vector<std::pair<size_t, size_t>> setupHashRanges(DBReader<DBKeyType> &reader, const BaseMatrix *subMat, size_t maxSeqLen,
                                                              size_t firstKey, size_t seqCount, size_t chunkSeqCount, size_t maxRangeSeqCount) {
    size_t *hashDist = new(std::nothrow) size_t[USHRT_MAX + 1];
    Util::checkAllocation(hashDist, "Can not allocate hashDist memory");
    memset(hashDist, 0, sizeof(size_t) * (USHRT_MAX + 1));
    for (size_t chunkBegin = 0; chunkBegin < seqCount; chunkBegin += chunkSeqCount) {
        const size_t chunkEnd = std::min(chunkBegin + chunkSeqCount, seqCount);
        reader.loadData(firstKey + chunkBegin, chunkEnd - chunkBegin);
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            Sequence seq(maxSeqLen, reader.getDbtype(), subMat, 0, false, false);
#pragma omp for schedule(dynamic, 100)
            for (size_t i = chunkBegin; i < chunkEnd; i++) {
                __sync_fetch_and_add(&hashDist[static_cast<unsigned short>(seqHash(reader, seq, firstKey + i, thread_idx))], 1);
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
            hashRanges.emplace_back(rangeStart, bucket - 1);
            rangeStart = bucket;
            rangeSeqCount = 0;
        }
        rangeSeqCount += hashDist[bucket];
    }
    hashRanges.emplace_back(rangeStart, USHRT_MAX);
    delete[] hashDist;
    return hashRanges;
}

// keeps the sequences of keys [firstKey, firstKey + seqCount) whose hash bucket is in hashRange, copying them only if seqCopies is given
static size_t fillHashRange(DBReader<DBKeyType> &reader, const BaseMatrix *subMat, size_t maxSeqLen, size_t firstKey, size_t seqCount, size_t chunkSeqCount,
                            const std::pair<size_t, size_t> &hashRange, HashEntry *hashIdPairs, HashRangeMembers *hashRangeMembers, char *seqCopies) {
    const size_t length = reader.getSeqLen(firstKey);
    size_t rangeSeqCount = 0;
    for (size_t chunkBegin = 0; chunkBegin < seqCount; chunkBegin += chunkSeqCount) {
        const size_t chunkEnd = std::min(chunkBegin + chunkSeqCount, seqCount);
        reader.loadData(firstKey + chunkBegin, chunkEnd - chunkBegin);
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            Sequence seq(maxSeqLen, reader.getDbtype(), subMat, 0, false, false);
            // a thread marks whole groups of 64 positions, so chunks start at a multiple of 64
#pragma omp for schedule(dynamic, 1)
            for (size_t i = chunkBegin / 64; i < (chunkEnd + 63) / 64; i++) {
                hashRangeMembers[i].inRange = 0;
                for (size_t pos = i * 64; pos < std::min(i * 64 + 64, chunkEnd); pos++) {
                    const unsigned short bucket = static_cast<unsigned short>(seqHash(reader, seq, firstKey + pos, thread_idx));
                    if (bucket >= hashRange.first && bucket <= hashRange.second) {
                        hashRangeMembers[i].inRange |= 1ull << (pos % 64);
                    }
                }
            }
#pragma omp single
            for (size_t i = chunkBegin / 64; i < (chunkEnd + 63) / 64; i++) {
                hashRangeMembers[i].copiesBefore = rangeSeqCount;
                rangeSeqCount += __builtin_popcountll(hashRangeMembers[i].inRange);
            }
#pragma omp for schedule(dynamic, 1)
            for (size_t i = chunkBegin / 64; i < (chunkEnd + 63) / 64; i++) {
                size_t slot = hashRangeMembers[i].copiesBefore;
                for (uint64_t inRange = hashRangeMembers[i].inRange; inRange != 0; inRange &= inRange - 1) {
                    const size_t id = firstKey + i * 64 + __builtin_ctzll(inRange);
                    hashIdPairs[slot].hash = seqHash(reader, seq, id, thread_idx);
                    hashIdPairs[slot].id = id;
                    if (seqCopies != NULL) {
                        memcpy(seqCopies + slot * length, reader.getData(id, thread_idx), length);
                    }
                    slot++;
                }
            }
        }
    }
    return rangeSeqCount;
}

int clusthashfast(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE, 5));
    par.seqIdThr = (float)Parameters::CLUST_HASH_DEFAULT_MIN_SEQ_ID / 100.0f;
    par.parseParameters(argc, argv, command, true, 0, 0);
    MMseqsNode::init();

    DBReader<DBKeyType> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_LINFINITY);
    reader.open(DBReader<DBKeyType>::NOSORT);
    if (Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_AMINO_ACIDS) == false) {
        Debug(Debug::ERROR) << "clusthashfast only supports amino acid databases\n";
        EXIT(EXIT_FAILURE);
    }
    // a compressed or padded entry is decoded into one buffer per thread, so the representative is copied before its members are read
    const bool copyRepSeq = reader.isCompressed() || (DBReader<DBKeyType>::getExtendedDbtype(reader.getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU);
    const bool isDBLinfinity = DBReader<DBKeyType>::getExtendedDbtype(reader.getDbtype()) & Parameters::DBTYPE_EXTENDED_LINFINITY;
    if (isDBLinfinity == false && MMseqsNode::isMaster() == false) {
        reader.close();
        return EXIT_SUCCESS;
    }
    if (isDBLinfinity == false && par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        reader.readMmapedDataInMemory();
    }

    size_t localFirstKey = 0;
    size_t localKeySize = reader.getSize();
    std::vector<LengthStart> nodeSeqRange;
    if (isDBLinfinity) {
        reader.decomposeDomainByAminoAcid(MMseqsNode::rank, MMseqsNode::count, &localFirstKey, &localKeySize);
        const DBReader<DBKeyType>::Index *index = reader.getIndex();
        size_t offset = 0;
        for (size_t row = 0; row < reader.getIndexSize() && index[row].id < localFirstKey + localKeySize; row++) {
            if (index[row].id >= localFirstKey && (nodeSeqRange.empty() || index[row].length != index[row - 1].length)) {
                nodeSeqRange.push_back({index[row].id, offset});
            }
            const size_t rangeEnd = (row + 1 < reader.getIndexSize()) ? index[row + 1].id : reader.getSize();
            offset += (rangeEnd - index[row].id) * index[row].length;
        }
        nodeSeqRange.push_back({localFirstKey + localKeySize, offset});
    } else {
        nodeSeqRange.push_back({0, 0});
        nodeSeqRange.push_back({reader.getSize(), reader.getDataSize()});
    }

    if (isDBLinfinity && par.compressed) {
        Debug(Debug::WARNING) << "Compressed output cannot be combined with --createdb-mode " << Parameters::SEQUENCE_SPLIT_MODE_LINFINITY_LENGTH << " databases\n";
        Debug(Debug::WARNING) << "We recompute with --compressed 0\n";
        par.compressed = 0;
    }
    const size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);
    const std::string resultDb = isDBLinfinity ? par.db2 + "." + SSTR(MMseqsNode::rank) : par.db2;
    DBWriter writer(resultDb.c_str(), (resultDb + ".index").c_str(), 1, par.compressed, isDBLinfinity ? Parameters::DBTYPE_OMIT_FILE : Parameters::DBTYPE_CLUSTER_RES);
    writer.open();
    const std::string representativeDb = par.db2 + "_redundancy";
    DBWriter representativeWriter(representativeDb.c_str(), (representativeDb + ".index").c_str(), 1, 0, Parameters::DBTYPE_OMIT_FILE);
    if (isDBLinfinity == false) {
        representativeWriter.open();
    }
    // a mode 3 node keeps its representatives as one bit per key
    const size_t representativeBitsSize = isDBLinfinity ? (localKeySize + 63) / 64 : 0;
    uint64_t *representativeBits = new(std::nothrow) uint64_t[representativeBitsSize]();
    Util::checkAllocation(representativeBits, "Can not allocate representative bit memory");
    if (representativeBitsSize * sizeof(uint64_t) >= memoryLimit) {
        Debug(Debug::ERROR) << "Not enough memory to run clusthashfast. Memory limit: " << memoryLimit
                            << " bytes, the representative bits need: " << representativeBitsSize * sizeof(uint64_t) << " bytes.\n"
                            << "Raise RAM size or --split-memory-limit.\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t splitMemoryLimit = static_cast<size_t>(static_cast<double>(memoryLimit - representativeBitsSize * sizeof(uint64_t)) * 0.95);

    // consecutive lengths that fit together are loaded at once (a standard DB is mmapped, so only a mode 3 DB counts its sequences);
    // a length that does not fit is loaded alone and processed once per hash range that fits, after 2 bits per key
    std::vector<size_t> hashRangeSeqCount(nodeSeqRange.size() - 1, 0);
    std::vector<size_t> loadStarts;
    size_t loadBytes = 0;
    for (size_t i = 0; i + 1 < nodeSeqRange.size(); i++) {
        const size_t lengthSeqCount = nodeSeqRange[i + 1].startKey - nodeSeqRange[i].startKey;
        const size_t lengthBytes = (isDBLinfinity ? nodeSeqRange[i + 1].offset - nodeSeqRange[i].offset : 0) + lengthSeqCount * sizeof(HashEntry);
        if (lengthBytes > splitMemoryLimit) {
            // a mode 3 length is read in chunks of half the limit and copied, a standard DB is mmapped and its resultIndex get the whole limit
            const size_t length = isDBLinfinity ? reader.getSeqLen(nodeSeqRange[i].startKey) : 0;
            const size_t rangeBytes = isDBLinfinity ? splitMemoryLimit / 2 : splitMemoryLimit;
            const size_t keptBytes = lengthSeqCount / 4;
            if ((isDBLinfinity && rangeBytes / (length + 2) < 64) || rangeBytes <= keptBytes + length + sizeof(HashEntry)) {
                Debug(Debug::ERROR) << "Not enough memory to hash " << lengthSeqCount << " sequences\n"
                                    << "Raise RAM size or --split-memory-limit.\n";
                EXIT(EXIT_FAILURE);
            }
            hashRangeSeqCount[i] = (rangeBytes - keptBytes) / (length + sizeof(HashEntry));
        }
        if (loadStarts.empty() || loadBytes + lengthBytes > splitMemoryLimit) {
            loadStarts.push_back(i);
            loadBytes = 0;
        }
        loadBytes += lengthBytes;
    }
    loadStarts.push_back(nodeSeqRange.size() - 1);
    // the sequence buffer of reader.loadData and hashIdPairs are kept while the next load fits next to them
    HashEntry *hashIdPairs = NULL;
    size_t hashIdPairsSize = 0;
    size_t loadedDataBytes = 0;

    std::string result;
    char buffer[32];
    bool needsSort = false;
    size_t clusters = 0;
    std::string resultIndex;
    SubstitutionMatrix sMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, -0.2);
    ReducedMatrix subMat(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2num, sMat.num2aa, sMat.alphabetSize, par.alphabetSize.values.aminoacid(), 2.0);

    Debug(Debug::INFO) << "Hashing " << localKeySize << " sequences in " << loadStarts.size() - 1 << " loads\n";
    Debug::Progress progress(loadStarts.size() - 1);
    for (size_t load = 0; load + 1 < loadStarts.size(); load++) {
        progress.updateProgress();
        const size_t loadFirstKey = nodeSeqRange[loadStarts[load]].startKey;
        const size_t loadSeqCount = nodeSeqRange[loadStarts[load + 1]].startKey - loadFirstKey;
        const bool isLengthGroupExceedLimits = hashRangeSeqCount[loadStarts[load]] > 0;
        const size_t length = reader.getSeqLen(loadFirstKey);
        const bool copiesSeqs = isLengthGroupExceedLimits && isDBLinfinity;
        const size_t chunkSeqCount = copiesSeqs ? splitMemoryLimit / 2 / (length + 2) / 64 * 64 : loadSeqCount;
        const size_t loadDataBytes = isDBLinfinity == false ? 0 : (copiesSeqs ? chunkSeqCount * (length + 2) : nodeSeqRange[loadStarts[load + 1]].offset - nodeSeqRange[loadStarts[load]].offset);
        const size_t hashedSeqCapacity = isLengthGroupExceedLimits ? hashRangeSeqCount[loadStarts[load]] : loadSeqCount;
        // besides the sequence buffer a load needs hashIdPairs; a split mode 3 length also its copies in the other half of the limit
        const size_t hashBytes = isLengthGroupExceedLimits ? (copiesSeqs ? splitMemoryLimit / 2 : splitMemoryLimit) : hashedSeqCapacity * sizeof(HashEntry);
        // the sequence buffer and hashIdPairs are kept while this load fits next to them, otherwise they get the size of this load
        if (loadedDataBytes + hashBytes > splitMemoryLimit) {
            reader.releaseDataPages();
            loadedDataBytes = 0;
        }
        loadedDataBytes = std::max(loadedDataBytes, loadDataBytes);
        if (hashIdPairsSize < hashedSeqCapacity || loadedDataBytes + hashBytes + (hashIdPairsSize - hashedSeqCapacity) * sizeof(HashEntry) > splitMemoryLimit) {
            delete[] hashIdPairs;
            hashIdPairs = new(std::nothrow) HashEntry[hashedSeqCapacity];
            Util::checkAllocation(hashIdPairs, "Can not allocate hash entry memory");
            hashIdPairsSize = hashedSeqCapacity;
        }

        std::vector<std::pair<size_t, size_t>> hashRanges(1, std::make_pair(0, USHRT_MAX));
        HashRangeMembers *hashRangeMembers = NULL;
        char *seqCopies = NULL;
        if (isLengthGroupExceedLimits) {
            hashRanges = setupHashRanges(reader, &subMat, par.maxSeqLen, loadFirstKey, loadSeqCount, chunkSeqCount, hashedSeqCapacity);
            needsSort |= hashRanges.size() > 1;
            hashRangeMembers = new(std::nothrow) HashRangeMembers[loadSeqCount / 64 + 1];
            Util::checkAllocation(hashRangeMembers, "Can not allocate hash range memory");
            if (copiesSeqs) {
                seqCopies = new(std::nothrow) char[hashedSeqCapacity * length];
                Util::checkAllocation(seqCopies, "Can not allocate hash range memory");
            }
        }
        for (size_t range = 0; range < hashRanges.size(); range++) {
            const size_t hashedSeqCount = isLengthGroupExceedLimits
                                        ? fillHashRange(reader, &subMat, par.maxSeqLen, loadFirstKey, loadSeqCount, chunkSeqCount, hashRanges[range], hashIdPairs, hashRangeMembers, seqCopies)
                                        : fillHashIdPairs(reader, &subMat, par.maxSeqLen, loadFirstKey, loadSeqCount, hashIdPairs);
            SORT_PARALLEL(hashIdPairs, hashIdPairs + hashedSeqCount, HashEntry::compareByHashAndId);

            // bins never cross threads; resultIndex are fixed before hash is overwritten with the representative id
            std::vector<size_t> threadOffsets(par.threads + 1, hashedSeqCount);
            for (int thread = 0; thread < par.threads; thread++) {
                const size_t pos = hashedSeqCount / par.threads * thread;
                threadOffsets[thread] = (pos == 0) ? 0 : std::upper_bound(hashIdPairs + pos, hashIdPairs + hashedSeqCount, HashEntry{hashIdPairs[pos - 1].hash, std::numeric_limits<DBLocalId>::max()}, HashEntry::compareByHashAndId) - hashIdPairs;
            }
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (int thread = 0; thread < par.threads; thread++) {
                const size_t threadEnd = threadOffsets[thread + 1];
                std::string repSeqCopy;
                for (size_t hashBinBegin = threadOffsets[thread]; hashBinBegin < threadEnd;) {
                    size_t hashBinEnd = hashBinBegin + 1;
                    while (hashBinEnd < threadEnd && hashIdPairs[hashBinEnd].hash == hashIdPairs[hashBinBegin].hash) {
                        hashBinEnd++;
                    }
                    // [hashBinBegin, i) is assigned, [i, hashBinEnd) is not: the smallest id left is the next representative, members are swapped in front
                    for (size_t i = hashBinBegin; i < hashBinEnd;) {
                        std::iter_swap(hashIdPairs + i, std::min_element(hashIdPairs + i, hashIdPairs + hashBinEnd, HashEntry::compareById));
                        const DBLocalId repId = hashIdPairs[i].id;
                        hashIdPairs[i++].hash = repId;
                        if (i == hashBinEnd) {
                            break;
                        }
                        const unsigned int repLength = reader.getSeqLen(repId);
                        const char *repSeq = (seqCopies == NULL) ? reader.getData(repId, thread) : seqCopies + copySlot(hashRangeMembers, repId - loadFirstKey) * length;
                        if (copyRepSeq) {
                            repSeqCopy.assign(repSeq, repLength);
                            repSeq = repSeqCopy.data();
                        }
                        for (size_t j = i; j < hashBinEnd; j++) {
                            const DBLocalId memberId = hashIdPairs[j].id;
                            if (reader.getSeqLen(memberId) != repLength) {
                                continue;
                            }
                            const unsigned int distance = DistanceCalculator::computeInverseHammingDistance(repSeq, (seqCopies == NULL) ? reader.getData(memberId, thread) : seqCopies + copySlot(hashRangeMembers, memberId - loadFirstKey) * length, repLength);
                            if (static_cast<float>(distance) / static_cast<float>(repLength) >= par.seqIdThr) {
                                std::swap(hashIdPairs[i], hashIdPairs[j]);
                                hashIdPairs[i++].hash = repId;
                            }
                        }
                    }
                    hashBinBegin = hashBinEnd;
                }
            }

            // hash now holds the representative id, so sorting again gathers each cluster with its representative first
            SORT_PARALLEL(hashIdPairs, hashIdPairs + hashedSeqCount, HashEntry::compareByHashAndId);
            const size_t rangeOffset = writer.getOffset(0);
            for (size_t clusterBegin = 0; clusterBegin < hashedSeqCount;) {
                size_t clusterEnd = clusterBegin + 1;
                while (clusterEnd < hashedSeqCount && hashIdPairs[clusterEnd].hash == hashIdPairs[clusterBegin].hash) {
                    clusterEnd++;
                }
                const DBLocalId repId = hashIdPairs[clusterBegin].id;
                const DBKeyType repKey = reader.getDbKey(repId);
                if (isDBLinfinity) {
                    representativeBits[(repId - localFirstKey) / 64] |= 1ull << ((repId - localFirstKey) % 64);
                } else {
                    representativeWriter.writeIndexEntry(repKey, reader.getOffset(repId), reader.getEntryLen(repId), 0);
                    clusters++;
                }
                // a mode 3 database leaves singletons out: in the representative subset but not in a cluster
                if (isDBLinfinity == false || clusterEnd - clusterBegin > 1) {
                    for (size_t i = clusterBegin; i < clusterEnd; i++) {
                        char *end = Itoa::u64toa_sse2(static_cast<uint64_t>(reader.getDbKey(hashIdPairs[i].id)), buffer);
                        *(end - 1) = '\n';
                        result.append(buffer, end - buffer);
                    }
                    writer.writeData(result.c_str(), result.length(), repKey, 0, true, isDBLinfinity == false);
                    result.clear();
                }
                clusterBegin = clusterEnd;
            }
            // a mode 3 cluster DB is indexed by one line per processing unit: offset, length and node of its clusters
            if (isDBLinfinity && writer.getOffset(0) > rangeOffset) {
                resultIndex += SSTR(rangeOffset) + "\t" + SSTR(writer.getOffset(0) - rangeOffset) + "\t" + SSTR(MMseqsNode::rank) + "\n";
            }
        }
        delete[] hashRangeMembers;
        delete[] seqCopies;
    }

    // clusters are written in key order unless a length was split into hash resultIndex
    writer.close(false, needsSort);
    if (isDBLinfinity) {
        FileUtil::remove((resultDb + ".index").c_str());
        // the bits are on disk before the index part that tells the master this node is done
        const std::string bitsFile = representativeDb + ".bits." + SSTR(MMseqsNode::rank);
        FILE *bits = FileUtil::openAndDelete(bitsFile.c_str(), "w");
        const uint64_t header[2] = { localFirstKey, localKeySize };
        if (fwrite(header, sizeof(uint64_t), 2, bits) != 2 || fwrite(representativeBits, sizeof(uint64_t), representativeBitsSize, bits) != representativeBitsSize
            || fflush(bits) != 0 || fsync(fileno(bits)) != 0 || fclose(bits) != 0) {
            Debug(Debug::ERROR) << "Cannot write " << bitsFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        // the rows of this node in the sequence index, with the number of its representatives as a sixth column
        std::string representativeIndex;
        MemoryMapped indexData(par.db1Index, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        const char *columns[6];
        for (char *row = (char *) indexData.getData(); row < (char *) indexData.getData() + indexData.size(); ) {
            char *rowEnd = Util::skipLine(row);
            Util::getWordsOfLine(row, columns, 4);
            const size_t rowFirstKey = Util::fast_atoi<size_t>(columns[0]);
            const size_t rowEndKey = rowFirstKey + Util::fast_atoi<size_t>(columns[3]);
            if (rowFirstKey >= localFirstKey + localKeySize) {
                break;
            }
            if (rowFirstKey >= localFirstKey) {
                size_t rowRepresentatives = 0;
                for (size_t bit = rowFirstKey - localFirstKey; bit < rowEndKey - localFirstKey;) {
                    const size_t take = std::min<size_t>(64 - bit % 64, rowEndKey - localFirstKey - bit);
                    const uint64_t word = representativeBits[bit / 64] >> (bit % 64);
                    rowRepresentatives += __builtin_popcountll(take < 64 ? word & ((1ull << take) - 1) : word);
                    bit += take;
                }
                representativeIndex.append(row, rowEnd - row - 1);
                representativeIndex += "\t" + SSTR(rowRepresentatives) + "\n";
            }
            row = rowEnd;
        }
        indexData.close();
        MMseqsNode::publishFile(par.db2Index + "." + SSTR(MMseqsNode::rank), resultIndex);
        MMseqsNode::publishFile(representativeDb + ".index." + SSTR(MMseqsNode::rank), representativeIndex);
        if (MMseqsNode::isMaster()) {
            MMseqsNode::publishFile(par.db2Index, MMseqsNode::gatherFiles(par.db2Index + "."));
            representativeIndex = MMseqsNode::gatherFiles(representativeDb + ".index.");
            for (char *row = (char *) representativeIndex.c_str(); *row != '\0'; row = Util::skipLine(row)) {
                Util::getWordsOfLine(row, columns, 6);
                clusters += Util::fast_atoi<size_t>(columns[5]);
            }
            MMseqsNode::publishFile(representativeDb + ".index", representativeIndex);
        }
    } else {
        representativeWriter.close(false, needsSort);
    }
    if (MMseqsNode::isMaster()) {
        DBReader<DBKeyType>::softlinkDb(par.db1, representativeDb, DBFiles::DATA);
        DBWriter::writeDbtypeFile(representativeDb.c_str(), reader.getDbtype(), reader.isCompressed());
        if (isDBLinfinity) {
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
