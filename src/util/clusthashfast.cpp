#include "DBWriter.h"
#include "FileUtil.h"
#include "Util.h"
#include "Parameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "ReducedMatrix.h"
#include "SubstitutionMatrix.h"
#include "DistanceCalculator.h"
#include "Orf.h"
#include "FastSort.h"
#include "itoa.h"
#include "Timer.h"

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <mutex>

#ifdef OPENMP
#include <omp.h>
#endif


static const size_t CLUSTHASH_CHUNK_RECORDS = 65536;
static const size_t CLUSTHASH_RESULT_RESERVE = 4096;
static const size_t CLUSTHASH_MAX_RETAINED_RESULT = 4 * 1024 * 1024;
static const size_t CLUSTHASH_MAX_RETAINED_CLAIMS = 1024 * 1024;
static const unsigned int CLUSTHASH_MAX_PARTITIONS = 4096;
static const size_t CLUSTHASH_PENDING_BYTES = 8 * 1024 * 1024;
static const size_t CLUSTHASH_PROGRESS_STEP = 1ull << 20;

static size_t progressSteps(size_t dbSize) {
    const size_t steps = (dbSize + CLUSTHASH_PROGRESS_STEP - 1) / CLUSTHASH_PROGRESS_STEP;
    return (steps == 0) ? 1 : steps;
}

struct ClusterCounts {
    size_t clusters;
    size_t merged;
    size_t runs;

    ClusterCounts() : clusters(0), merged(0), runs(0) {}

    void add(const ClusterCounts &other) {
        clusters += other.clusters;
        merged += other.merged;
        runs += other.runs;
    }
};

struct __attribute__((packed)) HashEntry {
    size_t hash;
    DBLocalId id;

    static bool compareByHashAndId(const HashEntry &first, const HashEntry &second) {
        if (first.hash != second.hash) {
            return first.hash < second.hash;
        }
        return first.id < second.id;
    }
};

static void appendKey(std::string &buffer, DBKeyType key, char *lineBuffer) {
    char *outPos = Itoa::u64toa_sse2(static_cast<uint64_t>(key), lineBuffer);
    buffer.append(lineBuffer, (outPos - lineBuffer - 1));
    buffer.push_back('\n');
}

struct ClusterWorker {
    std::string result;
    std::string querySeq;
    std::vector<unsigned char> claimed;
    std::vector<unsigned int> memberLength;
    char lineBuffer[32];
    ClusterCounts counts;

    ClusterWorker() {
        result.reserve(CLUSTHASH_RESULT_RESERVE);
    }

    void beginCluster(DBKeyType representativeKey) {
        if (result.capacity() > CLUSTHASH_MAX_RETAINED_RESULT) {
            std::string().swap(result);
            result.reserve(CLUSTHASH_RESULT_RESERVE);
        } else {
            result.clear();
        }
        appendKey(result, representativeKey, lineBuffer);
    }

    void addMember(DBKeyType memberKey) {
        appendKey(result, memberKey, lineBuffer);
        counts.merged++;
    }

    void writeCluster(DBWriter &writer, DBWriter &representativeWriter, DBReader<DBKeyType> &reader,
                      DBLocalId representativeId, DBKeyType representativeKey, unsigned int thread_idx) {
        writer.writeData(result.c_str(), result.length(), representativeKey, thread_idx);
        representativeWriter.writeIndexEntry(representativeKey, reader.getOffset(representativeId),
                                             reader.getEntryLen(representativeId), thread_idx);
        counts.clusters++;
    }

    void trimClaimed() {
        if (claimed.capacity() > CLUSTHASH_MAX_RETAINED_CLAIMS) {
            std::vector<unsigned char>().swap(claimed);
        }
        if (memberLength.capacity() > CLUSTHASH_MAX_RETAINED_CLAIMS) {
            std::vector<unsigned int>().swap(memberLength);
        }
    }
};

static size_t hashNucleotideSequence(const char *data, size_t length) {
    const size_t A = 31;
    size_t h1 = 0;
    size_t h2 = 0;
    for (size_t i = 0; i < length; ++i) {
        h1 = ((h1 * A) + data[i]);
        h2 = ((h2 * A) + Orf::complement(data[length - i - 1]));
    }
    return std::min(h1, h2);
}

static size_t hashWithLength(size_t hash, size_t length) {
    return hash ^ (length + static_cast<size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) + (hash >> 2));
}

static bool hashScanIsSequential(DBReader<DBKeyType> &reader) {
    return reader.isSortedByOffset();
}

static size_t idScanBlock(DBReader<DBKeyType> &reader, size_t dbSize, int threads) {
    const size_t DEFAULT_BLOCK = 1000;
    const size_t MIN_BLOCKS_PER_THREAD = 16;
    const size_t TARGET_SPAN = 64 * 1024 * 1024;
    if (threads < 2 || dbSize == 0 || hashScanIsSequential(reader) == false) {
        return DEFAULT_BLOCK;
    }
    size_t fileSize = 0;
    for (size_t fileIdx = 0; fileIdx < reader.getDataFileCnt(); fileIdx++) {
        fileSize += reader.getDataSizeForFile(fileIdx);
    }
    size_t lastEnd = std::min(reader.getIndex(dbSize - 1)->offset + reader.getEntryLen(dbSize - 1), fileSize);
    if (fileSize == 0 || lastEnd == 0) {
        return DEFAULT_BLOCK;
    }
    size_t stride = std::max<size_t>(lastEnd / dbSize, 1);
    size_t wanted = (TARGET_SPAN + stride - 1) / stride;
    size_t balanced = dbSize / (static_cast<size_t>(threads) * MIN_BLOCKS_PER_THREAD);
    return std::max(std::min(wanted, balanced), DEFAULT_BLOCK);
}

static size_t hashOf(DBReader<DBKeyType> &reader, size_t id, Sequence *seq, unsigned int thread_idx) {
    const size_t length = reader.getSeqLen(id);
    const char *data = reader.getData(id, thread_idx);
    if (seq == NULL) {
        return hashWithLength(hashNucleotideSequence(data, length), length);
    }
    seq->mapSequence(id, 0, data, length);
    return hashWithLength(Util::hash(seq->numSequence, seq->L), length);
}

static size_t hashRunEnd(const HashEntry *entries, size_t entryCount, size_t runBegin) {
    const size_t hash = entries[runBegin].hash;
    size_t runEnd = runBegin + 1;
    while (runEnd < entryCount && entries[runEnd].hash == hash) {
        runEnd++;
    }
    return runEnd;
}

static size_t firstOwnedRun(const HashEntry *entries, size_t chunkBegin, size_t chunkEnd) {
    size_t runBegin = chunkBegin;
    while (runBegin < chunkEnd && runBegin > 0 && entries[runBegin - 1].hash == entries[runBegin].hash) {
        runBegin++;
    }
    return runBegin;
}

static void clusterHashRun(DBReader<DBKeyType> &reader, DBWriter &writer, DBWriter &representativeWriter,
                           const HashEntry *run, size_t runSize, float seqIdThr, ClusterWorker &worker,
                           unsigned int thread_idx) {
    worker.counts.runs++;
    if (runSize == 1) {
        const DBKeyType representativeKey = reader.getDbKey(run[0].id);
        worker.beginCluster(representativeKey);
        worker.writeCluster(writer, representativeWriter, reader, run[0].id, representativeKey, thread_idx);
        return;
    }

    worker.memberLength.resize(runSize);
    for (size_t k = 0; k < runSize; ++k) {
        const size_t memberLen = reader.getSeqLen(run[k].id);
        if (memberLen > UINT_MAX) {
            Debug(Debug::ERROR) << "Sequence " << reader.getDbKey(run[k].id) << " is " << memberLen
                                << " residue long; clusthashfast stores lengths in 32 bit.\n";
            EXIT(EXIT_FAILURE);
        }
        worker.memberLength[k] = static_cast<unsigned int>(memberLen);
    }

    worker.claimed.assign(runSize, 0);
    for (size_t i = 0; i < runSize; i++) {
        if (worker.claimed[i]) {
            continue;
        }
        worker.claimed[i] = 1;
        const DBLocalId queryId = run[i].id;
        const unsigned int queryLength = worker.memberLength[i];
        const DBKeyType representativeKey = reader.getDbKey(queryId);
        worker.beginCluster(representativeKey);

        const char *querySeq = NULL;
        for (size_t j = i + 1; j < runSize; j++) {
            if (worker.claimed[j] || worker.memberLength[j] != queryLength) {
                continue;
            }
            if (querySeq == NULL) {
                worker.querySeq.assign(reader.getData(queryId, thread_idx), queryLength);
                querySeq = worker.querySeq.data();
            }
            const char *targetSeq = reader.getData(run[j].id, thread_idx);
            const unsigned int distance =
                DistanceCalculator::computeInverseHammingDistance(querySeq, targetSeq, queryLength);
            const float seqId = static_cast<float>(distance) / static_cast<float>(queryLength);
            if (seqId >= seqIdThr) {
                worker.addMember(reader.getDbKey(run[j].id));
                worker.claimed[j] = 1;
            }
        }
        worker.writeCluster(writer, representativeWriter, reader, queryId, representativeKey, thread_idx);
    }
    worker.trimClaimed();
}

static ClusterCounts clusterHashRuns(DBReader<DBKeyType> &reader, DBWriter &writer,
                                     DBWriter &representativeWriter, const HashEntry *entries,
                                     size_t entryCount, float seqIdThr, bool showProgress, int threads) {
    const size_t chunkCount = (entryCount + CLUSTHASH_CHUNK_RECORDS - 1) / CLUSTHASH_CHUNK_RECORDS;
    Debug::Progress progress(chunkCount);
    size_t totalClusters = 0;
    size_t totalMerged = 0;
    size_t totalRuns = 0;
#pragma omp parallel num_threads(threads) reduction(+:totalClusters) reduction(+:totalMerged) reduction(+:totalRuns)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        ClusterWorker worker;

#pragma omp for schedule(dynamic, 1)
        for (size_t chunk = 0; chunk < chunkCount; ++chunk) {
            const size_t chunkBegin = chunk * CLUSTHASH_CHUNK_RECORDS;
            const size_t chunkEnd = std::min(entryCount, chunkBegin + CLUSTHASH_CHUNK_RECORDS);
            if (showProgress) {
                progress.updateProgress(chunk);
            }
            size_t pos = firstOwnedRun(entries, chunkBegin, chunkEnd);
            while (pos < chunkEnd) {
                const size_t runEnd = hashRunEnd(entries, entryCount, pos);
                clusterHashRun(reader, writer, representativeWriter, entries + pos, runEnd - pos,
                               seqIdThr, worker, thread_idx);
                pos = runEnd;
            }
        }
        totalClusters += worker.counts.clusters;
        totalMerged += worker.counts.merged;
        totalRuns += worker.counts.runs;
    }

    ClusterCounts counts;
    counts.clusters = totalClusters;
    counts.merged = totalMerged;
    counts.runs = totalRuns;
    return counts;
}



struct HashPartitions {
    std::vector<FILE *> files;
    std::vector<std::string> names;
    std::vector<size_t> counts;
    unsigned int bits;

    HashPartitions() : bits(0) {}

    ~HashPartitions() {
        close();
        for (size_t part = 0; part < names.size(); part++) {
            if (FileUtil::fileExists(names[part].c_str())) {
                FileUtil::remove(names[part].c_str());
            }
        }
    }

    bool open(const std::string &prefix, unsigned int partitions) {
        bits = 0;
        while ((1u << bits) < partitions) {
            bits++;
        }
        files.assign(partitions, NULL);
        counts.assign(partitions, 0);
        for (unsigned int part = 0; part < partitions; part++) {
            names.push_back(prefix + "." + SSTR(part));
            files[part] = fopen(names[part].c_str(), "wb");
            if (files[part] == NULL) {
                return false;
            }
        }
        return true;
    }

    void close() {
        for (size_t part = 0; part < files.size(); part++) {
            if (files[part] != NULL) {
                fclose(files[part]);
                files[part] = NULL;
            }
        }
    }

    // a short sequence leaves the top bits of Util::hash zero, so mix before shifting
    unsigned int partitionOf(size_t hash) const {
        if (bits == 0) {
            return 0u;
        }
        size_t mixed = hash + 0x9e3779b97f4a7c15ULL;
        mixed = (mixed ^ (mixed >> 30)) * 0xbf58476d1ce4e5b9ULL;
        mixed = (mixed ^ (mixed >> 27)) * 0x94d049bb133111ebULL;
        mixed ^= mixed >> 31;
        return static_cast<unsigned int>(mixed >> (sizeof(size_t) * 8 - bits));
    }
};

static bool flushPartition(HashPartitions &parts, std::vector<std::mutex> &locks, unsigned int partition,
                           std::vector<HashEntry> &pending) {
    std::lock_guard<std::mutex> lock(locks[partition]);
    const bool ok = fwrite(pending.data(), sizeof(HashEntry), pending.size(), parts.files[partition])
                    == pending.size();
    parts.counts[partition] += pending.size();
    pending.clear();
    return ok;
}

static bool hashIntoPartitions(DBReader<DBKeyType> &reader, HashPartitions &parts, bool isNuclInput,
                               BaseMatrix *subMat, size_t maxSeqLen, bool showProgress, int threads) {
    const size_t dbSize = reader.getSize();
    const size_t scanChunk = idScanBlock(reader, dbSize, threads);
    const unsigned int partitionCount = static_cast<unsigned int>(parts.files.size());
    const size_t flush = std::max<size_t>(64, CLUSTHASH_PENDING_BYTES
                                              / (partitionCount * sizeof(HashEntry)));
    std::vector<std::mutex> locks(partitionCount);
    Debug::Progress progress(progressSteps(dbSize));
    bool ok = true;
#pragma omp parallel num_threads(threads) reduction(&&:ok)
    {
        bool threadOk = true;
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<std::vector<HashEntry> > pending(partitionCount);
        Sequence *seq = isNuclInput ? NULL : new Sequence(maxSeqLen, reader.getDbtype(), subMat, 0, false, false);
#pragma omp for schedule(static, scanChunk)
        for (size_t id = 0; id < dbSize; ++id) {
            if (showProgress && (id & (CLUSTHASH_PROGRESS_STEP - 1)) == 0) {
                progress.updateProgress();
            }
            HashEntry entry;
            entry.hash = hashOf(reader, id, seq, thread_idx);
            entry.id = static_cast<DBLocalId>(id);
            const unsigned int partition = parts.partitionOf(entry.hash);
            pending[partition].push_back(entry);
            if (pending[partition].size() >= flush) {
                threadOk = flushPartition(parts, locks, partition, pending[partition]) && threadOk;
            }
        }
        for (unsigned int partition = 0; partition < partitionCount; partition++) {
            if (pending[partition].empty() == false) {
                threadOk = flushPartition(parts, locks, partition, pending[partition]) && threadOk;
            }
        }
        delete seq;
        ok = ok && threadOk;
    }
    return ok;
}

static unsigned int hashPartitionCount(size_t dbSize, size_t budget) {
    const size_t needed = dbSize * sizeof(HashEntry);
    const size_t usable = std::max<size_t>(budget, 1);
    unsigned int partitions = 1;
    while (partitions < CLUSTHASH_MAX_PARTITIONS && needed / partitions > usable) {
        partitions *= 2;
    }
    return partitions;
}

static void hashAllSequences(DBReader<DBKeyType> &reader, HashEntry *entries, size_t dbSize,
                             bool isNuclInput, BaseMatrix *subMat, size_t maxSeqLen, bool showProgress,
                             int threads) {
    const size_t scanChunk = idScanBlock(reader, dbSize, threads);
    Debug::Progress progress(progressSteps(dbSize));
#pragma omp parallel num_threads(threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence *seq = isNuclInput ? NULL : new Sequence(maxSeqLen, reader.getDbtype(), subMat, 0, false, false);
#pragma omp for schedule(static, scanChunk)
        for (size_t id = 0; id < dbSize; ++id) {
            if (showProgress && (id & (CLUSTHASH_PROGRESS_STEP - 1)) == 0) {
                progress.updateProgress();
            }
            entries[id].hash = hashOf(reader, id, seq, thread_idx);
            entries[id].id = static_cast<DBLocalId>(id);
        }
        delete seq;
    }
}


static ClusterCounts clusterSequences(const Parameters &par, DBReader<DBKeyType> &reader, DBWriter &writer,
                                      DBWriter &representativeWriter, BaseMatrix *subMat, bool isNuclInput,
                                      size_t memoryLimit) {
    ClusterCounts counts;
    const size_t dbSize = reader.getSize();
    if (dbSize == 0) {
        return counts;
    }

    const bool showProgress = (Debug::debugLevel >= Debug::INFO);
    const std::string tmpPrefix = std::string(writer.getDataFileName());
    const size_t entryBudget = std::max<size_t>(memoryLimit / 2, 1);
    const unsigned int partitionCount = hashPartitionCount(dbSize, entryBudget);
    Debug(Debug::INFO) << "Memory limit " << memoryLimit << " byte, entries "
                       << dbSize * sizeof(HashEntry) / partitionCount << " byte per partition\n";

    if (reader.isSortedByOffset()) {
        reader.setSequentialAdvice();
    }

    if (partitionCount == 1) {
        Debug(Debug::INFO) << "Hashing " << dbSize << " sequences...\n";
        HashEntry *entries = new(std::nothrow) HashEntry[dbSize];
        Util::checkAllocation(entries, "Cannot allocate hash entry memory in clusthashfast");
        hashAllSequences(reader, entries, dbSize, isNuclInput, subMat, par.maxSeqLen, showProgress,
                         par.threads);
        Debug(Debug::INFO) << "Sort sequence hashes...\n";
        Timer sortTimer;
        SORT_PARALLEL(entries, entries + dbSize, HashEntry::compareByHashAndId);
        Debug(Debug::INFO) << "Time for sorting hashes: " << sortTimer.lap() << "\n";
        Debug(Debug::INFO) << "Cluster equal length sequences...\n";
        counts = clusterHashRuns(reader, writer, representativeWriter, entries, dbSize, par.seqIdThr,
                                 showProgress, par.threads);
        delete[] entries;
        Debug(Debug::INFO) << "Found " << counts.runs << " unique hash/length groups\n";
        return counts;
    }

    Debug(Debug::INFO) << "Hashing " << dbSize << " sequences into " << partitionCount << " partitions...\n";
    HashPartitions parts;
    if (parts.open(tmpPrefix + ".hashpart", partitionCount) == false) {
        Debug(Debug::ERROR) << "Cannot open hash partitions under " << tmpPrefix << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (hashIntoPartitions(reader, parts, isNuclInput, subMat, par.maxSeqLen, showProgress,
                           par.threads) == false) {
        Debug(Debug::ERROR) << "Cannot write hash partitions under " << tmpPrefix << "\n";
        EXIT(EXIT_FAILURE);
    }
    parts.close();

    std::vector<HashEntry> part;
    for (unsigned int partition = 0; partition < partitionCount; partition++) {
        const size_t partSize = parts.counts[partition];
        if (partSize == 0) {
            continue;
        }
        part.resize(partSize);
        Timer partReadTimer;
        FILE *in = fopen(parts.names[partition].c_str(), "rb");
        if (in == NULL || fread(part.data(), sizeof(HashEntry), partSize, in) != partSize) {
            Debug(Debug::ERROR) << "Cannot read hash partition " << parts.names[partition] << "\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(in);
        Debug(Debug::INFO) << "Time for reading partition " << (partition + 1) << ": "
                           << partReadTimer.lap() << "\n";
        FileUtil::remove(parts.names[partition].c_str());

        Timer partSortTimer;
        SORT_PARALLEL(part.data(), part.data() + partSize, HashEntry::compareByHashAndId);
        Debug(Debug::INFO) << "Time for sorting partition " << (partition + 1) << ": "
                           << partSortTimer.lap() << "\n";
        const ClusterCounts partCounts = clusterHashRuns(reader, writer, representativeWriter, part.data(),
                                                        partSize, par.seqIdThr, showProgress, par.threads);
        counts.add(partCounts);
        Debug(Debug::INFO) << "Partition " << (partition + 1) << "/" << partitionCount << ": "
                           << partCounts.runs << " groups, " << partCounts.clusters << " clusters\n";
    }
    Debug(Debug::INFO) << "Found " << counts.runs << " unique hash/length groups\n";
    return counts;
}

int clusthashfast(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.alphabetSize = MultiParam<NuclAA<int> >(NuclAA<int>(Parameters::CLUST_HASH_DEFAULT_ALPH_SIZE, 5));
    par.seqIdThr = static_cast<float>(Parameters::CLUST_HASH_DEFAULT_MIN_SEQ_ID) / 100.0f;
    par.parseParameters(argc, argv, command, true, 0, 0);

    // NOSORT still orders by key, which is what makes the representative the smallest key
    DBReader<DBKeyType> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                               DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX);
    reader.open(DBReader<DBKeyType>::NOSORT);
    const bool dataFitsInMemory = reader.getDataSize() < Util::getTotalSystemMemory() / 2;
    if (par.preloadMode == Parameters::PRELOAD_MODE_MMAP_TOUCH
        || (par.preloadMode != Parameters::PRELOAD_MODE_MMAP && dataFitsInMemory)) {
        reader.readMmapedDataInMemory();
    }

    const bool isNuclInput = Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    BaseMatrix *subMat = NULL;
    if (isNuclInput == false) {
        SubstitutionMatrix sMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, -0.2);
        subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2num, sMat.num2aa,
                                   sMat.alphabetSize, par.alphabetSize.values.aminoacid(), 2.0);
    }

    const size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed,
                    Parameters::DBTYPE_CLUSTER_RES);
    writer.open();

    const std::string representativeDb = par.db2 + "_redundancy";
    DBWriter representativeWriter(representativeDb.c_str(), (representativeDb + ".index").c_str(),
                                  par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    representativeWriter.open();

    const ClusterCounts counts = clusterSequences(par, reader, writer, representativeWriter, subMat,
                                                  isNuclInput, memoryLimit);

    if (subMat != NULL) {
        delete subMat;
    }
    writer.close();

    const bool mergeRepresentativeIndex = Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_HMM_PROFILE)
        || Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_AMINO_ACIDS)
        || Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    representativeWriter.close(mergeRepresentativeIndex, true);
    DBReader<DBKeyType>::softlinkDb(par.db1, representativeDb, DBFiles::DATA);
    DBWriter::writeDbtypeFile(representativeDb.c_str(), reader.getDbtype(), reader.isCompressed());
    DBReader<DBKeyType>::softlinkDb(par.db1, representativeDb, DBFiles::SEQUENCE_ANCILLARY);
    reader.close();

    Debug(Debug::INFO) << counts.clusters << " clusters, " << counts.merged << " sequences merged by hamming distance\n";
    return EXIT_SUCCESS;
}
