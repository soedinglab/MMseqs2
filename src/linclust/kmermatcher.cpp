// include xxhash early to avoid incompatibilites with SIMDe
#define XXH_INLINE_ALL
#include "xxhash.h"


#include "kmermatcher.h"
#include "Debug.h"
#include "Indexer.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "NucleotideMatrix.h"
#include "QueryMatcher.h"
#include "KmerGenerator.h"
#include "MarkovKmerScore.h"
#include "FileUtil.h"
#include "FastSort.h"
#include "Timer.h"
#include "SequenceWeights.h"
#include "Masker.h"

#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/statvfs.h>
#include <fcntl.h>
#include <unistd.h>

#include <limits>
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <zstd.h>

#ifdef OPENMP
#include <omp.h>
#endif
#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) -1)
#endif

// SEQUENTIAL only widens the fault window; this starts the next chunk's read so the fault does not wait
static void prefetchScanChunk(DBReader<DBKeyType> &reader, size_t startId, size_t count) {
    if (count == 0 || startId >= reader.getSize()) {
        return;
    }
    if (count > reader.getSize() - startId) {
        count = reader.getSize() - startId;
    }
    const size_t lastId = startId + count - 1;
    char *begin = reader.getDataUncompressed(startId);
    char *end = reader.getDataUncompressed(lastId) + reader.getEntryLen(lastId);
    if (end <= begin) {
        return;
    }
    // madvise needs a page aligned address, and widening outward costs nothing for a readahead hint
    const size_t page = Util::getPageSize();
    const uintptr_t first = reinterpret_cast<uintptr_t>(begin) & ~(uintptr_t)(page - 1);
    const uintptr_t last = (reinterpret_cast<uintptr_t>(end) + page - 1) & ~(uintptr_t)(page - 1);
    if (last > first) {
        const int rc = posix_madvise(reinterpret_cast<void *>(first),
                                     static_cast<size_t>(last - first), POSIX_MADV_WILLNEED);
        if (rc != 0) {
            static bool warned = false;
            if (warned == false) {
                warned = true;
                Debug(Debug::WARNING) << "posix_madvise(WILLNEED) failed: " << rc << "\n";
            }
        }
    }
}

static size_t kmerScanChunkSize(DBReader<DBKeyType> &reader, size_t startId, size_t cnt, int threads) {
    const size_t DEFAULT_CHUNK = 100;
    const size_t MIN_CHUNKS_PER_THREAD = 16;
    if (threads < 2 || cnt == 0 || reader.getDataFileCnt() != 1 || reader.isSortedByOffset() == false) {
        return DEFAULT_CHUNK;
    }
    size_t fileSize = reader.getDataSizeForFile(0);
    if (fileSize == 0) {
        return DEFAULT_CHUNK;
    }
    size_t firstOffset = reader.getIndex(startId)->offset;
    size_t lastId = startId + cnt - 1;
    // a soft-linked createdb-mode 1 db has no room for the entry terminator, so clamp to the file
    size_t lastEnd = std::min(reader.getIndex(lastId)->offset + reader.getEntryLen(lastId), fileSize);
    if (lastEnd <= firstOffset) {
        return DEFAULT_CHUNK;
    }
    // one page table maps pageSize/sizeof(void*) pages, so derive its span instead of hardcoding 2 MB
    size_t pageSize = Util::getPageSize();
    size_t pmdSpan = pageSize * (pageSize / sizeof(void *));
    size_t stride = std::max<size_t>((lastEnd - firstOffset) / cnt, 1);
    size_t wanted = (pmdSpan + stride - 1) / stride;
    size_t balanced = cnt / (static_cast<size_t>(threads) * MIN_CHUNKS_PER_THREAD);
    return std::max(std::min(wanted, balanced), DEFAULT_CHUNK);
}

uint64_t hashUInt64(uint64_t in, uint64_t seed) {
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    in = __builtin_bswap64(in);
#endif
    return XXH64(&in, sizeof(uint64_t), seed);
}


template <typename T, bool includeAdjacency, bool IncludeSeqLen>
KmerPosition<T, includeAdjacency, IncludeSeqLen> *initKmerPositionMemory(size_t size) {
    KmerPosition<T, includeAdjacency, IncludeSeqLen> * hashSeqPair = new(std::nothrow) KmerPosition<T, includeAdjacency, IncludeSeqLen>[size + 1];
    Util::checkAllocation(hashSeqPair, "Can not allocate memory");
    size_t pageSize = Util::getPageSize()/sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>);
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (size_t page = 0; page < size+1; page += pageSize) {
            size_t readUntil = std::min(size+1, page + pageSize) - page;
            memset(hashSeqPair+page, 0xFF, sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>)* readUntil);
        }
    }
    return hashSeqPair;
}

// batches the atomic reservation into the shared k-mer array; a contention/memory trade-off only
static const size_t KMER_STAGING_BUFFER_SIZE = 65536;

static void removeKmerTmpFileIfExists(const std::string &fileName);
static FILE *openKmerTmpFileForOverwriteOrDie(const std::string &fileName, const char *mode);
static void printKmerTmpFsState(const std::string &fileName);
static void kmerTmpIoErrorDie(const char *op, const std::string &fileName);

// --compress-kmer-tmp-files selects the codec for every k-mer spill file: 0 raw, 1 zstd
static const int KMER_TMP_CODEC_NONE = 0;
static const int KMER_TMP_CODEC_ZSTD = 1;

// the suffix is the only on-disk marker of the codec, so readers recover it from the file name
static const char *kmerTmpCodecSuffix(int codec) {
    if (codec == KMER_TMP_CODEC_ZSTD) {
        return ".zst";
    }
    return "";
}

// buckets are transient spill, so they compress at a lower level than the split-result tmp files
static const int KMER_BUCKET_ZSTD_LEVEL = 2;
// non-template so they can be declared here and defined after the stream writers
static void *bucketWriterOpen(const std::string &fileName, int codec, int zstdLevel);
static void bucketWriterAppend(void *writer, int codec, const void *data, size_t byteSize);
static void bucketWriterClose(void *writer, int codec);
static size_t bucketReadFile(const std::string &fileName, int codec, void *dst, size_t maxBytes,
                             size_t ioBufferBytes, bool &overflow);

// a handler without SA_RESTART turns a signal into EINTR, and the read itself is idempotent
static ssize_t kmerPreadRetry(int fd, void *buf, size_t len, size_t offset) {
    ssize_t got;
    do {
        got = pread(fd, buf, len, static_cast<off_t>(offset));
    } while (got < 0 && errno == EINTR);
    return got;
}

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
static void flushKmerBuffer(KmerPosition<T, includeAdjacency, IncludeSeqLen> *kmerArray,
                            size_t kmerArraySize,
                            KmerPosition<T, includeAdjacency, IncludeSeqLen> *threadKmerBuffer,
                            size_t bufferPos,
                            size_t *offset) {
    size_t writeOffset = __sync_fetch_and_add(offset, bufferPos);
    if (writeOffset > kmerArraySize || bufferPos > kmerArraySize - writeOffset) {
        Debug(Debug::ERROR) << "Kmer array overflow. currKmerArrayOffset=" << writeOffset
                            << ", kmerBufferPos=" << bufferPos
                            << ", kmerArraySize=" << kmerArraySize << ".\n";
        EXIT(EXIT_FAILURE);
    }
    if (kmerArray != NULL) {
        memcpy(kmerArray + writeOffset, threadKmerBuffer,
               sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) * bufferPos);
    }
}

// one pass routes every k-mer to a per-thread file for its split's range, so no split has to rescan
static std::string kmerBucketCountFileName(const std::string &base, size_t bucket) {
    return base + "_" + SSTR(bucket) + ".cnt";
}

static std::string kmerBucketFileName(const std::string &base, size_t bucket, int tid, int codec) {
    return base + "_" + SSTR(bucket) + "_" + SSTR(tid) + kmerTmpCodecSuffix(codec);
}

// hash (unsigned short) -> bucket index, from the contiguous split hash ranges.
static unsigned int *buildHashToBucketLookup(const std::vector<std::pair<size_t, size_t>> &ranges) {
    unsigned int *lut = new(std::nothrow) unsigned int[USHRT_MAX + 1];
    Util::checkAllocation(lut, "Can not allocate hash-to-bucket lookup");
    for (size_t b = 0; b < ranges.size(); b++) {
        size_t hi = std::min(ranges[b].second, static_cast<size_t>(USHRT_MAX));
        for (size_t h = ranges[b].first; h <= hi; h++) {
            lut[h] = static_cast<unsigned int>(b);
        }
    }
    return lut;
}

static const size_t KMERMATCHER_PROGRESS_STEP = 1ull << 20;

static const size_t KMER_SPOOL_BUFFER_BUDGET = 2ull * 1024 * 1024 * 1024;
static const size_t KMER_SPOOL_BUFFER_MIN = 64ull * 1024;
static const size_t KMER_SPOOL_BUFFER_MAX = 8ull * 1024 * 1024;

static size_t kmerSpoolBufferBytes(size_t slots) {
    if (slots == 0) {
        return KMER_SPOOL_BUFFER_MAX;
    }
    size_t perSlot = KMER_SPOOL_BUFFER_BUDGET / slots;
    return std::min(std::max(perSlot, KMER_SPOOL_BUFFER_MIN), KMER_SPOOL_BUFFER_MAX);
}

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
struct KmerPartitionSink {
    typedef KmerPosition<T, includeAdjacency, IncludeSeqLen> KP;
    size_t bucketBuffer;
    int numThreads;
    size_t numBuckets;
    int codec;
    const unsigned int *hashToBucket;   // [USHRT_MAX+1]
    std::string base;
    std::vector<void *> writers;        // indexed [tid * numBuckets + bucket]: FILE* or stream writer
    std::vector<KP *> buffers;
    std::vector<size_t> bufPos;
    std::vector<size_t> recordsWritten;

    KmerPartitionSink(const std::string &base, int numThreads, size_t numBuckets,
                      const unsigned int *hashToBucket, int codec)
        : numThreads(numThreads), numBuckets(numBuckets), codec(codec), hashToBucket(hashToBucket), base(base) {
        size_t slots = static_cast<size_t>(numThreads) * numBuckets;
        bucketBuffer = std::max<size_t>(kmerSpoolBufferBytes(slots) / sizeof(KP), 1);
        writers.assign(slots, NULL);
        buffers.assign(slots, NULL);
        bufPos.assign(slots, 0);
        recordsWritten.assign(slots, 0);
        for (int tid = 0; tid < numThreads; tid++) {
            for (size_t b = 0; b < numBuckets; b++) {
                size_t k = static_cast<size_t>(tid) * numBuckets + b;
                writers[k] = bucketWriterOpen(kmerBucketFileName(base, b, tid, codec), codec,
                                              KMER_BUCKET_ZSTD_LEVEL);
                buffers[k] = new(std::nothrow) KP[bucketBuffer];
                Util::checkAllocation(buffers[k], "Can not allocate k-mer bucket buffer");
            }
        }
    }

    inline void flushSlot(size_t k) {
        if (bufPos[k] == 0) {
            return;
        }
        bucketWriterAppend(writers[k], codec, buffers[k], sizeof(KP) * bufPos[k]);
        recordsWritten[k] += bufPos[k];
        bufPos[k] = 0;
    }

    // concurrent, but each thread touches only its own [tid] slots, so no locking is needed
    inline void emit(int tid, const KP &rec, unsigned int hash) {
        size_t k = static_cast<size_t>(tid) * numBuckets + hashToBucket[hash];
        buffers[k][bufPos[k]++] = rec;
        if (bufPos[k] >= bucketBuffer) {
            flushSlot(k);
        }
    }

    // slots share no state, so closing them in parallel keeps this off the threads x buckets path
    void finish() {
        size_t slots = static_cast<size_t>(numThreads) * numBuckets;
#pragma omp parallel for schedule(static)
        for (size_t k = 0; k < slots; k++) {
            flushSlot(k);
            bucketWriterClose(writers[k], codec);
            delete[] buffers[k];
        }
        // per-thread record counts let loadKmerBucket place each file without decompressing first
        for (size_t b = 0; b < numBuckets; b++) {
            std::string countFile = kmerBucketCountFileName(base, b);
            FILE *cf = fopen(countFile.c_str(), "w");
            if (cf == NULL) {
                kmerTmpIoErrorDie("open", countFile);
            }
            for (int tid = 0; tid < numThreads; tid++) {
                if (fprintf(cf, "%zu\n", recordsWritten[static_cast<size_t>(tid) * numBuckets + b]) < 0) {
                    kmerTmpIoErrorDie("write", countFile);
                }
            }
            if (fclose(cf) != 0) {
                kmerTmpIoErrorDie("close", countFile);
            }
        }
    }
};

// reads and deletes one bucket's per-thread files into arr[0..count); the sentinel tail stays
template <typename T, bool includeAdjacency, bool IncludeSeqLen>
size_t loadKmerBucket(const std::string &base, size_t bucket, int numThreads,
                      KmerPosition<T, includeAdjacency, IncludeSeqLen> *arr, size_t cap, int codec) {
    typedef KmerPosition<T, includeAdjacency, IncludeSeqLen> KP;
    // the partition pass recorded per-thread counts, so every file has a fixed slot and reads run in parallel
    std::vector<size_t> counts;
    std::string countFile = kmerBucketCountFileName(base, bucket);
    FILE *cf = fopen(countFile.c_str(), "r");
    if (cf != NULL) {
        size_t v;
        while (fscanf(cf, "%zu", &v) == 1) {
            counts.push_back(v);
        }
        fclose(cf);
    }
    // one bucket file per thread is in flight at a time, so they share the budget between them
    const size_t ioBufferBytes = kmerSpoolBufferBytes(static_cast<size_t>(numThreads));
    if (counts.size() == static_cast<size_t>(numThreads)) {
        std::vector<size_t> offset(static_cast<size_t>(numThreads) + 1, 0);
        for (int tid = 0; tid < numThreads; tid++) {
            offset[tid + 1] = offset[tid] + counts[tid];
        }
        if (offset[numThreads] > cap) {
            Debug(Debug::ERROR) << "k-mer bucket " << bucket << " exceeds the allocated split array\n";
            EXIT(EXIT_FAILURE);
        }
        bool failed = false;
#pragma omp parallel for schedule(dynamic, 1)
        for (int tid = 0; tid < numThreads; tid++) {
            std::string fileName = kmerBucketFileName(base, bucket, tid, codec);
            if (counts[tid] == 0) {
                removeKmerTmpFileIfExists(fileName);
                continue;
            }
            // skipping here would leave sentinel records in the array, silently dropping the k-mers
            if (FileUtil::fileExists(fileName.c_str()) == false) {
                Debug(Debug::ERROR) << "k-mer bucket file " << fileName << " is missing but its count file records "
                                    << counts[tid] << " entries\n";
                failed = true;
                continue;
            }
            bool overflow = false;
            size_t want = counts[tid] * sizeof(KP);
            size_t bytes = bucketReadFile(fileName, codec, arr + offset[tid], want, ioBufferBytes, overflow);
            if (overflow || bytes != want) {
                failed = true;
            }
            removeKmerTmpFileIfExists(fileName);
        }
        if (failed) {
            Debug(Debug::ERROR) << "k-mer bucket " << bucket << " does not match its recorded record count\n";
            EXIT(EXIT_FAILURE);
        }
        removeKmerTmpFileIfExists(countFile);
        return offset[numThreads];
    }
    size_t count = 0;
    for (int tid = 0; tid < numThreads; tid++) {
        std::string fileName = kmerBucketFileName(base, bucket, tid, codec);
        if (FileUtil::fileExists(fileName.c_str()) == false) {
            continue;
        }
        bool overflow = false;
        size_t bytes = bucketReadFile(fileName, codec, arr + count, (cap - count) * sizeof(KP), ioBufferBytes, overflow);
        // whole KmerPosition records only, so a partial one means the file or the cap is inconsistent
        if (overflow || (bytes % sizeof(KP)) != 0) {
            Debug(Debug::ERROR) << "k-mer bucket " << bucket << " exceeds the allocated split array\n";
            EXIT(EXIT_FAILURE);
        }
        count += bytes / sizeof(KP);
        removeKmerTmpFileIfExists(fileName);
    }
    return count;
}

template <int TYPE, typename T, bool includeAdjacency, bool IncludeSeqLen>
std::pair<size_t, size_t> fillKmerPositionArray(KmerPosition<T, includeAdjacency, IncludeSeqLen> * kmerArray, size_t kmerArraySize, DBReader<DBKeyType> &seqDbr,
                                                Parameters & par, BaseMatrix * subMat, bool hashWholeSequence,
                                                size_t hashStartRange, size_t hashEndRange, size_t * hashDistribution,
                                                KmerPartitionSink<T, includeAdjacency, IncludeSeqLen> *partitionSink){
    size_t offset = 0;
    int querySeqType  =  seqDbr.getDbtype();
    size_t longestKmer = par.kmerSize;
    const unsigned char xIndex = subMat->aa2num[static_cast<int>('X')];


    ScoreMatrix two;
    ScoreMatrix three;
    if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
        two = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        three = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
    }

    // Debug::Progress paints the id-1..id delta, so it has to count updates, not sequences
    const size_t progressSteps = (seqDbr.getSize() + KMERMATCHER_PROGRESS_STEP - 1) / KMERMATCHER_PROGRESS_STEP;
    Debug::Progress progress(progressSteps == 0 ? 1 : progressSteps);
#pragma omp parallel num_threads(par.threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        size_t threadLongestKmer = par.kmerSize;
        unsigned short * scoreDist= new(std::nothrow) unsigned short[65536];
        Util::checkAllocation(scoreDist, "Can not allocate scoreDist memory in fillKmerPositionArray");
        unsigned int * hierarchicalScoreDist= new(std::nothrow) unsigned int[128];
        Util::checkAllocation(hierarchicalScoreDist, "Can not allocate hierarchicalScoreDist memory in fillKmerPositionArray");
        // Zero once; kept clean by clearing only touched bins after each sequence (below).
        memset(scoreDist, 0, sizeof(unsigned short) * 65536);
        memset(hierarchicalScoreDist, 0, sizeof(unsigned int) * 128);

        Masker *masker = NULL;
        if (par.maskMode == 1) {
            masker = new Masker(*subMat);
        }
        const int adjustedKmerSize = (par.adjustKmerLength) ? std::min( par.kmerSize+5, 23) :   par.kmerSize;
        Sequence seq(par.maxSeqLen, querySeqType, subMat, adjustedKmerSize, par.spacedKmer, false, true, par.spacedKmerPattern);
        KmerGenerator* generator;
        if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
            generator = new KmerGenerator( par.kmerSize, subMat->alphabetSize, 150);
            generator->setDivideStrategy(&three, &two);
        }
        Indexer idxer(subMat->alphabetSize - 1,  par.kmerSize);
        const unsigned int BUFFER_SIZE = static_cast<unsigned int>(KMER_STAGING_BUFFER_SIZE);
        size_t bufferPos = 0;
        KmerPosition<T, includeAdjacency, IncludeSeqLen> * threadKmerBuffer = NULL;
        if (hashDistribution == NULL) {
            threadKmerBuffer = new(std::nothrow) KmerPosition<T, includeAdjacency, IncludeSeqLen>[BUFFER_SIZE];
            Util::checkAllocation(threadKmerBuffer, "Can not allocate threadKmerBuffer memory in fillKmerPositionArray");
        }
        SequencePosition * kmers = (SequencePosition *) malloc((par.pickNbest * (par.maxSeqLen + 1) + 1) * sizeof(SequencePosition) * 2);
        Util::checkAllocation(kmers, "Can not allocate kmers memory in fillKmerPositionArray");
        size_t kmersArraySize = par.maxSeqLen;
        const size_t flushSize = 100000000;
        size_t iterations = static_cast<size_t>(ceil(static_cast<double>(seqDbr.getSize()) / static_cast<double>(flushSize)));
        for (size_t i = 0; i < iterations; i++) {
            size_t start = (i * flushSize);
            size_t bucketSize = std::min(seqDbr.getSize() - (i * flushSize), flushSize);

            size_t scanChunk = kmerScanChunkSize(seqDbr, start, bucketSize, par.threads);
            // the pointer arithmetic above needs one mapping, and the span needs offset-ordered ids
            const bool canPrefetch = seqDbr.getDataFileCnt() == 1 && seqDbr.isSortedByOffset()
                                     && seqDbr.isDirectIo() == false;
// every thread in the team computes the same chunk from shared read-only state, as the schedule needs
#pragma omp for schedule(dynamic, scanChunk)
            for (size_t id = start; id < (start + bucketSize); id++) {
                // one atomic per 2^20 sequences, so the contention the sparse guard avoids stays avoided
                if ((id & (KMERMATCHER_PROGRESS_STEP - 1)) == 0) {
                    progress.updateProgress();
                }
                // the team sweeps forward, so ask for the next chunk while this one is still hashing
                if (canPrefetch && ((id - start) % scanChunk) == 0) {
                    prefetchScanChunk(seqDbr, id + scanChunk, scanChunk);
                }

                seq.mapSequence(id, seqDbr.getDbKey(id), seqDbr.getData(id, thread_idx), seqDbr.getSeqLen(id));

                size_t seqHash =  SIZE_T_MAX;
                //TODO, how to handle this in reverse?
                if(hashWholeSequence){
                    seqHash = Util::hash(seq.numSequence, seq.L);
                    seqHash = hashUInt64(seqHash, par.hashShift);
                }
                if(masker != NULL){
                    masker->maskSequence(seq, par.maskMode,  par.maskProb, par.maskLowerCaseMode, par.maskNrepeats);
                }
                size_t seqKmerCount = 0;
                DBKeyType seqId = (par.kmerMatcherMode == Parameters::KMERMATCHER_MODE_LOCAL)
                    ? static_cast<DBKeyType>(id) : seq.getDbKey();
                while (seq.hasNextKmer()) {
                    unsigned char *kmer = (unsigned char*) seq.nextKmer();
                    if(seq.kmerContainsX()){
                        continue;
                    }
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                        NucleotideMatrix * nuclMatrix = (NucleotideMatrix*)subMat;
                        size_t kmerLen =  par.kmerSize;
                        size_t kmerIdx = Indexer::computeKmerIdx(kmer, kmerLen);
                        size_t revkmerIdx = Util::revComplement(kmerIdx, kmerLen);
                        // skip forward and rev. identical k-mers.
                        // We can not know how to align these afterwards
                        if(revkmerIdx == kmerIdx){
                            continue;
                        }
                        bool pickReverseKmer = (revkmerIdx<kmerIdx);
                        kmerIdx = (pickReverseKmer) ? revkmerIdx : kmerIdx;
                        const unsigned short hash = hashUInt64(kmerIdx, par.hashShift);

                        if(par.adjustKmerLength) {
                            unsigned char revKmer[32];
                            unsigned char * kmerToHash = kmer;
                            if(pickReverseKmer){
                                for(int pos = static_cast<int>(adjustedKmerSize)-1; pos > -1; pos--){
                                    revKmer[(adjustedKmerSize - 1) - pos]=nuclMatrix->reverseResidue(kmer[pos]);
                                }
                                kmerToHash = revKmer;
                            }
                            kmerLen = MarkovKmerScore::adjustedLength(kmerToHash, adjustedKmerSize,
                                                                      (par.kmerSize - MarkovScores::MARKOV_ORDER) * MarkovScores::MEDIAN_SCORE);
                            threadLongestKmer = std::max(kmerLen, threadLongestKmer);
                            kmerIdx = Indexer::computeKmerIdx(kmerToHash, kmerLen);
                        }

                        // set signed bit for normal kmers to make the  SIZE_T_MAX logic easier
                        // reversed kmers do not have a signed bit
                        size_t kmerRev = (pickReverseKmer) ? BIT_CLEAR(kmerIdx, 63) : BIT_SET(kmerIdx, 63);
                        (kmers + seqKmerCount)->kmer = kmerRev;
                        int pos = seq.getCurrentPosition();
                        (kmers + seqKmerCount)->pos = (pickReverseKmer) ? (seq.L) - pos - kmerLen : pos;
                        (kmers + seqKmerCount)->score = hash;
                        scoreDist[hash]++;
                        hierarchicalScoreDist[hash >> 9]++;
                        seqKmerCount++;
                    } else if(TYPE == Parameters::DBTYPE_HMM_PROFILE) {
                        std::pair<size_t*, size_t>  scoreMat = generator->generateKmerList(kmer, true);
                        for(size_t kmerPos = 0; kmerPos < scoreMat.second && kmerPos < static_cast<size_t >(par.pickNbest); kmerPos++){
                            size_t kmerIdx = scoreMat.first[kmerPos];
                            (kmers + seqKmerCount)->kmer = kmerIdx;
                            (kmers + seqKmerCount)->pos = seq.getCurrentPosition();
                            const unsigned short hash = hashUInt64(kmerIdx, par.hashShift);
                            (kmers + seqKmerCount)->score = hash;
                            scoreDist[hash]++;
                            hierarchicalScoreDist[hash >> 9]++;
                            seqKmerCount++;
                        }
                    } else {
                        size_t kmerIdx = idxer.int2index(kmer, 0, par.kmerSize);
                        (kmers + seqKmerCount)->kmer = kmerIdx;
                        (kmers + seqKmerCount)->pos = seq.getCurrentPosition();
                        const unsigned short hash = hashUInt64(kmerIdx, par.hashShift);
                        (kmers + seqKmerCount)->score = hash;
                        scoreDist[hash]++;
                        hierarchicalScoreDist[hash >> 9]++;
                        seqKmerCount++;
                    }
                    if(seqKmerCount >= kmersArraySize){
                        kmersArraySize = seq.getMaxLen();
                        SequencePosition *newKmers = (SequencePosition *) realloc(kmers, (par.pickNbest * (kmersArraySize + 1) + 1) * sizeof(SequencePosition) * 2);
                        Util::checkAllocation(newKmers, "Can not reallocate kmers memory in fillKmerPositionArray");
                        kmers = newKmers;
                    }

                }
                float kmersPerSequenceScale = (TYPE == Parameters::DBTYPE_NUCLEOTIDES) ? par.kmersPerSequenceScale.values.nucleotide()
                                                                                       : par.kmersPerSequenceScale.values.aminoacid();
                size_t kmerConsidered = std::min(static_cast<size_t >(par.kmersPerSequence  - 1 + (kmersPerSequenceScale * seq.L)), seqKmerCount);

                unsigned int threshold = 0;
                size_t kmerInBins = 0;
                if (seqKmerCount > 0) {
                    size_t hierarchicaThreshold = 0;
                    for(hierarchicaThreshold = 0; hierarchicaThreshold < 128 && kmerInBins < kmerConsidered; hierarchicaThreshold++){
                        kmerInBins += hierarchicalScoreDist[hierarchicaThreshold];
                    }
                    hierarchicaThreshold -= (hierarchicaThreshold > 0) ? 1: 0;
                    kmerInBins -= hierarchicalScoreDist[hierarchicaThreshold];
                    for(threshold = hierarchicaThreshold*512; threshold <= USHRT_MAX && kmerInBins < kmerConsidered; threshold++){
                        kmerInBins += scoreDist[threshold];
                    }
                }
                int tooMuchElemInLastBin = (kmerInBins - kmerConsidered);

                // add k-mer to represent the identity
                if (static_cast<unsigned short>(seqHash) >= hashStartRange && static_cast<unsigned short>(seqHash) <= hashEndRange) {
                    if(hashDistribution != NULL){
                        __sync_fetch_and_add(&hashDistribution[static_cast<unsigned short>(seqHash)], 1);
                    }
                    else{
                        threadKmerBuffer[bufferPos].kmer = seqHash;
                        threadKmerBuffer[bufferPos].id = seqId;
                        threadKmerBuffer[bufferPos].pos = 0;
                        threadKmerBuffer[bufferPos].sl.setSeqLen(static_cast<T>(seq.L));
                        if (includeAdjacency) {
                            for (size_t i = 0; i < 6; i++) {
                                threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                            }
                        }
                        if (partitionSink != NULL) {
                            partitionSink->emit(thread_idx, threadKmerBuffer[bufferPos], static_cast<unsigned short>(seqHash));
                        } else {
                            bufferPos++;
                            if (bufferPos >= BUFFER_SIZE) {
                                flushKmerBuffer(kmerArray, kmerArraySize, threadKmerBuffer, bufferPos, &offset);
                                bufferPos = 0;
                            }
                        }
                    }
                }

                if(par.ignoreMultiKmer){
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        SORT_SERIAL(kmers, kmers + seqKmerCount, SequencePosition::compareByScoreReverse);
                    }else{
                        SORT_SERIAL(kmers, kmers + seqKmerCount, SequencePosition::compareByScore);
                    }
                }
                size_t selectedKmer = 0;
                for (size_t kmerIdx = 0; kmerIdx < seqKmerCount && selectedKmer < kmerConsidered; kmerIdx++) {

                    /* skip repeated kmer */
                    if (par.ignoreMultiKmer) {
                        size_t kmer = (kmers + kmerIdx)->kmer;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(kmer, 63);
                        }
                        if (kmerIdx + 1 < seqKmerCount) {
                            size_t nextKmer = (kmers + kmerIdx + 1)->kmer;
                            if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                                nextKmer = BIT_SET(nextKmer, 63);
                            }
                            if (kmer == nextKmer) {
                                while (kmer == nextKmer && kmerIdx < seqKmerCount) {
                                    kmerIdx++;
                                    if(kmerIdx >= seqKmerCount)
                                        break;
                                    nextKmer = (kmers + kmerIdx)->kmer;
                                    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                                        nextKmer = BIT_SET(nextKmer, 63);
                                    }
                                }
                            }
                        }
                        if(kmerIdx >= seqKmerCount)
                            break;
                    }

                    if ((kmers + kmerIdx)->score < threshold ){
                        if((kmers + kmerIdx)->score == (threshold - 1) && tooMuchElemInLastBin){
                            tooMuchElemInLastBin--;
                            threshold -= (tooMuchElemInLastBin == 0) ? 1 : 0;
                        }

                        selectedKmer++;
                        if ((kmers + kmerIdx)->score >= hashStartRange && (kmers + kmerIdx)->score <= hashEndRange)
                        {
                            if(hashDistribution != NULL){
                                __sync_fetch_and_add(&hashDistribution[(kmers + kmerIdx)->score], 1);
                                continue;
                            }
                            threadKmerBuffer[bufferPos].kmer = (kmers + kmerIdx)->kmer;
                            threadKmerBuffer[bufferPos].id = seqId;
                            threadKmerBuffer[bufferPos].pos = (kmers + kmerIdx)->pos;
                            threadKmerBuffer[bufferPos].sl.setSeqLen(static_cast<T>(seq.L));
                            if (includeAdjacency) {
                                unsigned int startPos = (kmers + kmerIdx)->pos;
                                unsigned int endPos = (kmers + kmerIdx)->pos + seq.getEffectiveKmerSize() - 1;
                                for (size_t i = 0; i < 6; i++) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(i, xIndex);
                                }

                                if (startPos >= 3) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(0, seq.numSequence[startPos - 3]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(1, seq.numSequence[startPos - 2]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(2, seq.numSequence[startPos - 1]);
                                } else if (startPos == 2) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(1, seq.numSequence[startPos - 2]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(2, seq.numSequence[startPos - 1]);
                                } else if (startPos == 1) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(2, seq.numSequence[startPos - 1]);
                                }

                                if (endPos + 3 <= static_cast<unsigned int>(seq.L) - 1) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(3, seq.numSequence[endPos + 1]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(4, seq.numSequence[endPos + 2]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(5, seq.numSequence[endPos + 3]);
                                } else if (endPos + 2 == static_cast<unsigned int>(seq.L) - 1) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(3, seq.numSequence[endPos + 1]);
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(4, seq.numSequence[endPos + 2]);
                                } else if (endPos + 1 == static_cast<unsigned int>(seq.L) - 1) {
                                    threadKmerBuffer[bufferPos].setAdjacentSeq(3, seq.numSequence[endPos + 1]);
                                }
                            }
                            if (partitionSink != NULL) {
                                partitionSink->emit(thread_idx, threadKmerBuffer[bufferPos], (kmers + kmerIdx)->score);
                            } else {
                                bufferPos++;
                                if (bufferPos >= BUFFER_SIZE) {
                                    flushKmerBuffer(kmerArray, kmerArraySize, threadKmerBuffer, bufferPos, &offset);
                                    bufferPos = 0;
                                }
                            }
                        }
                    }
                }
                // clear only the bins this sequence touched, not a 128 KB memset per sequence
                for (size_t k = 0; k < seqKmerCount; k++) {
                    scoreDist[(kmers + k)->score] = 0;
                    hierarchicalScoreDist[(kmers + k)->score >> 9] = 0;
                }
            }
#pragma omp barrier
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            if (thread_idx == 0) {
                // keeps the mapping: remapData reopens by name, which dies on an unlinked file and drops the advice
                seqDbr.dropCacheAll();
            }
#pragma omp barrier
        }
        if (masker != NULL) {
            delete masker;
        }

        if(threadKmerBuffer != NULL && bufferPos > 0){
            flushKmerBuffer(kmerArray, kmerArraySize, threadKmerBuffer, bufferPos, &offset);
        }
        free(kmers);
        delete[] threadKmerBuffer;
        delete[] hierarchicalScoreDist;
        delete[] scoreDist;
        if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
            delete generator;
        }
#pragma omp critical
        {
            longestKmer = std::max(threadLongestKmer, longestKmer);
        }
    }

    if (TYPE == Parameters::DBTYPE_HMM_PROFILE) {
        ExtendedSubstitutionMatrix::freeScoreMatrix(three);
        ExtendedSubstitutionMatrix::freeScoreMatrix(two);
    }

    return std::make_pair(offset, longestKmer);
}


template <int TYPE, typename T, bool includeAdjacency, bool IncludeSeqLen>
void swapCenterSequence(KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair, size_t splitKmerCount, SequenceWeights &seqWeights) {

    size_t prevHash = hashSeqPair[0].kmer;
    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
        prevHash = BIT_SET(prevHash, 63);
    }

    size_t repSeqPos = 0;
    size_t prevHashStart = 0;
    float repSeqWeight = seqWeights.getWeightById(hashSeqPair[repSeqPos].id);
    for (size_t elementIdx = 0; elementIdx < splitKmerCount; elementIdx++) {

        size_t currKmer = hashSeqPair[elementIdx].kmer;
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            currKmer = BIT_SET(currKmer, 63);
        }
        if (prevHash != currKmer) {

            if (repSeqPos != prevHashStart)
                std::swap(hashSeqPair[repSeqPos],hashSeqPair[prevHashStart]);

            prevHashStart = elementIdx;
            prevHash = hashSeqPair[elementIdx].kmer;
            if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                prevHash = BIT_SET(prevHash, 63);
            }
            repSeqPos = elementIdx;
            repSeqWeight = seqWeights.getWeightById(hashSeqPair[repSeqPos].id);
        }
        else {
            float currWeight = seqWeights.getWeightById(hashSeqPair[elementIdx].id);
            if (currWeight > repSeqWeight) {
                repSeqWeight = currWeight;
                repSeqPos = elementIdx;
            }
        }

        if (hashSeqPair[elementIdx].kmer == SIZE_T_MAX) {
            break;
        }

    }
}

template <int TYPE, typename T, bool includeAdjacency, bool IncludeSeqLen>
size_t assignGroup(KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair, KmerPosition<T, false, IncludeSeqLen> *writeSeqPair,
                    bool includeOnlyExtendable, int covMode, float covThr,
                    SequenceWeights *sequenceWeights, float weightThr, int threads,
                    std::vector<size_t> &threadOffsets, BaseMatrix *subMat,
                    AssignGroupMask assignGroupMask, ComputationPhase phase, uint8_t *countTable) {

    // Current assign group mode based on assignGroupMask
    const bool useAdjacentSeq = (includeAdjacency && hasFeature(assignGroupMask, AssignGroupFeature::AdjacentSeq));
    const bool useCountTable = hasFeature(assignGroupMask, AssignGroupFeature::CountTable);
    const bool isSetupCountTable = (phase == ComputationPhase::SetupCountTable);

    if (isSetupCountTable) {
        Debug(Debug::INFO) << "Assign group Mode: SetupCountTable: ";
    } else if (useAdjacentSeq) {
        Debug(Debug::INFO) << "Assign group Mode: Adjacent sequence: ";
    } else if (useCountTable) {
        Debug(Debug::INFO) << "Assign group Mode: CountTable: ";
    } else {
        Debug(Debug::INFO) << "Assign group Mode: Longest center: ";
    }

    std::vector<size_t> localWritePos;
    localWritePos.resize(threads);
    for (int thread = 0; thread < threads; thread++) {
        localWritePos[thread] = threadOffsets[thread];
    }

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int thread = 0; thread < threads; thread++) {
        size_t startIdx = threadOffsets[thread];
        size_t endIdx = threadOffsets[thread + 1];
        if (startIdx >= endIdx) {
            continue;
        }

        size_t prevHash = hashSeqPair[startIdx].kmer;
        size_t repSeqKey = hashSeqPair[startIdx].id;
        size_t repSeqId = repSeqKey;
        bool repIsReverse = false;

        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
            repIsReverse = (BIT_CHECK(hashSeqPair[startIdx].kmer, 63) == false);
            repSeqId = (repIsReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
            prevHash = BIT_SET(prevHash, 63);
        }

        size_t prevHashStart = startIdx;
        size_t prevSetSize = 0;
        size_t skipByWeightCount = 0;
        T queryLen = hashSeqPair[startIdx].sl.getSeqLen(hashSeqPair[startIdx].id);
        T repSeq_i_pos = hashSeqPair[startIdx].pos;

        short *subMatPos[6] = {NULL, NULL, NULL, NULL, NULL, NULL};

        // prepare subMatPos for adj mode
        if (useAdjacentSeq && hashSeqPair[prevHashStart].getAdjacentSeq(0) != UCHAR_MAX) {
            for (size_t i = 0; i < 6; i++) {
                subMatPos[i] = subMat->subMatrix[hashSeqPair[prevHashStart].getAdjacentSeq(i)];
            }
        }

        for (size_t elementIdx = startIdx; elementIdx <= endIdx; elementIdx++) {
            size_t currKmer = hashSeqPair[elementIdx].kmer;
            if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                currKmer = BIT_SET(currKmer, 63);
            }

            if (prevHash != currKmer) {

                // Phase 1: find and swap in the best representative for this mode
                if (useAdjacentSeq && subMatPos[0] != NULL) {
                    // find member with lowest adj score → swap to prevHashStart
                    size_t bestPos = prevHashStart;
                    int minAdjScore = INT_MAX;
                    T bestLen = 0;
                    for (size_t i = prevHashStart; i < elementIdx; i++) {
                        if (i > prevHashStart && sequenceWeights != nullptr &&
                            sequenceWeights->getWeightById(hashSeqPair[i].id) > weightThr) {
                            continue;
                        }
                        size_t kmer = hashSeqPair[i].kmer;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(kmer, 63);
                        }
                        if (kmer == SIZE_T_MAX) continue;

                        if (hashSeqPair[i].id == repSeqKey) {
                            hashSeqPair[i].setAdjacentSeq(0, UCHAR_MAX);
                        }
                        if (hashSeqPair[i].getAdjacentSeq(0) != UCHAR_MAX) {
                            int currAdjScore = 0;
                            for (size_t j = 0; j < 6; j++) {
                                currAdjScore += subMatPos[j][hashSeqPair[i].getAdjacentSeq(j)];
                            }
                            // an equal adjacency score goes to the longer member, which greedy makes the centre
                            const T currLen = hashSeqPair[i].sl.getSeqLen(hashSeqPair[i].id);
                            if (currAdjScore < minAdjScore
                                || (currAdjScore == minAdjScore && currLen > bestLen)) {
                                minAdjScore = currAdjScore;
                                bestLen = currLen;
                                bestPos = i;
                            }
                        }
                    }
                    if (bestPos != prevHashStart &&
                        hashSeqPair[bestPos].kmer != SIZE_T_MAX &&
                        hashSeqPair[bestPos].getAdjacentSeq(0) != UCHAR_MAX) {
                        std::swap(hashSeqPair[bestPos], hashSeqPair[prevHashStart]);
                    }
                } else if (useCountTable && countTable != NULL) {
                    // find member with highest count → swap to prevHashStart
                    size_t bestPos = prevHashStart;
                    int maxCount = -1;
                    for (size_t i = prevHashStart + 1; i < elementIdx; i++) {
                        if (sequenceWeights != nullptr &&
                            sequenceWeights->getWeightById(hashSeqPair[i].id) > weightThr) {
                            continue;
                        }
                        size_t kmer = hashSeqPair[i].kmer;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(kmer, 63);
                        }
                        if (kmer == SIZE_T_MAX) continue;

                        const size_t memid = hashSeqPair[i].id;
                        if (memid != repSeqKey) {
                            int cnt = countTable[memid];
                            if (cnt >= maxCount) {
                                if (includeAdjacency == false || hashSeqPair[i].getAdjacentSeq(0) != UCHAR_MAX) {
                                    maxCount = cnt;
                                    bestPos = i;
                                }
                            }
                        }
                    }
                    if (bestPos != prevHashStart &&
                        hashSeqPair[bestPos].kmer != SIZE_T_MAX &&
                        (includeAdjacency == false || hashSeqPair[bestPos].getAdjacentSeq(0) != UCHAR_MAX)) {
                        std::swap(hashSeqPair[bestPos], hashSeqPair[prevHashStart]);
                    }
                }

                // After swap, update rep info
                repSeqKey = hashSeqPair[prevHashStart].id;
                repSeqId = repSeqKey;
                repIsReverse = false;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    repIsReverse = (BIT_CHECK(hashSeqPair[prevHashStart].kmer, 63) == false);
                    repSeqId = (repIsReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
                }
                queryLen = hashSeqPair[prevHashStart].sl.getSeqLen(hashSeqPair[prevHashStart].id);
                repSeq_i_pos = hashSeqPair[prevHashStart].pos;

                // Phase 2: assign group members to representative
                bool skipProcessing = false;
                if (useAdjacentSeq) {
                    skipProcessing = (hashSeqPair[prevHashStart].getAdjacentSeq(0) == UCHAR_MAX);
                }

                if (skipProcessing == false) {
                    for (size_t i = prevHashStart; i < elementIdx; i++) {
                        if (i > prevHashStart && sequenceWeights != nullptr &&
                            sequenceWeights->getWeightById(hashSeqPair[i].id) > weightThr) {
                            continue;
                        }

                        size_t kmer = hashSeqPair[i].kmer;
                        if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                            kmer = BIT_SET(hashSeqPair[i].kmer, 63);
                        }

                        size_t rId = (kmer != SIZE_T_MAX) ? ((prevSetSize - skipByWeightCount == 1) ? SIZE_T_MAX : repSeqId) : SIZE_T_MAX;

                        if (rId != SIZE_T_MAX) {
                            int diagonal = repSeq_i_pos - hashSeqPair[i].pos;
                            if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                                bool targetIsReverse = (BIT_CHECK(hashSeqPair[i].kmer, 63) == false);
                                bool queryNeedsToBeRev = false;

                                T queryPos = 0;
                                T targetPos = 0;
                                T targetLen = hashSeqPair[i].sl.getSeqLen(hashSeqPair[i].id);

                                if (repIsReverse == true && targetIsReverse == false) {
                                    queryPos = repSeq_i_pos;
                                    targetPos = hashSeqPair[i].pos;
                                    queryNeedsToBeRev = true;
                                } else if (repIsReverse == true && targetIsReverse == true) {
                                    queryPos = (queryLen - 1) - repSeq_i_pos;
                                    targetPos = (targetLen - 1) - hashSeqPair[i].pos;
                                    queryNeedsToBeRev = false;
                                } else if (repIsReverse == false && targetIsReverse == true) {
                                    queryPos = (queryLen - 1) - repSeq_i_pos;
                                    targetPos = (targetLen - 1) - hashSeqPair[i].pos;
                                    queryNeedsToBeRev = true;
                                } else {
                                    queryPos = repSeq_i_pos;
                                    targetPos = hashSeqPair[i].pos;
                                    queryNeedsToBeRev = false;
                                }

                                diagonal = queryPos - targetPos;
                                rId = (queryNeedsToBeRev) ? BIT_CLEAR(rId, 63) : BIT_SET(rId, 63);
                            }

                            T targetLen = hashSeqPair[i].sl.getSeqLen(hashSeqPair[i].id);
                            bool canBeExtended = (diagonal < 0) || (diagonal > (queryLen - targetLen));
                            bool canBeCovered = Util::canBeCovered(covThr, covMode, static_cast<float>(queryLen), static_cast<float>(targetLen));

                            if ((includeOnlyExtendable == false && canBeCovered) ||
                                (canBeExtended && includeOnlyExtendable == true)) {
                                if (isSetupCountTable) {
                                    if (countTable != NULL) {
                                        uint8_t *slot = &countTable[hashSeqPair[i].id];
                                        // saturating, so the byte counter keeps the group ranking instead of wrapping
                                        uint8_t seen = *slot;
                                        while (seen < UCHAR_MAX
                                               && __sync_bool_compare_and_swap(slot, seen, seen + 1) == false) {
                                            seen = *slot;
                                        }
                                    }
                                } else {
                                    if (writeSeqPair != NULL) {
                                        if (queryLen < hashSeqPair[i].sl.getSeqLen(hashSeqPair[i].id) && covMode == Parameters::COV_MODE_TARGET) {
                                            writeSeqPair[localWritePos[thread]].kmer = hashSeqPair[i].id;
                                            writeSeqPair[localWritePos[thread]].pos = static_cast<short>(-diagonal);
                                            writeSeqPair[localWritePos[thread]].sl.setSeqLen(targetLen);
                                            writeSeqPair[localWritePos[thread]].id = rId;
                                        } else {
                                            writeSeqPair[localWritePos[thread]].kmer = rId;
                                            writeSeqPair[localWritePos[thread]].pos = static_cast<short>(diagonal);
                                            writeSeqPair[localWritePos[thread]].sl.setSeqLen(targetLen);
                                            writeSeqPair[localWritePos[thread]].id = hashSeqPair[i].id;
                                        }
                                    } else {
                                        if (queryLen < hashSeqPair[i].sl.getSeqLen(hashSeqPair[i].id) && covMode == Parameters::COV_MODE_TARGET) {
                                            hashSeqPair[localWritePos[thread]].kmer = hashSeqPair[i].id;
                                            hashSeqPair[localWritePos[thread]].pos = static_cast<short>(-diagonal);
                                            hashSeqPair[localWritePos[thread]].sl.setSeqLen(targetLen);
                                            hashSeqPair[localWritePos[thread]].id = rId;
                                        } else {
                                            hashSeqPair[localWritePos[thread]].kmer = rId;
                                            hashSeqPair[localWritePos[thread]].pos = static_cast<short>(diagonal);
                                            hashSeqPair[localWritePos[thread]].sl.setSeqLen(targetLen);
                                            hashSeqPair[localWritePos[thread]].id = hashSeqPair[i].id;
                                        }
                                    }
                                    localWritePos[thread]++;
                                }
                            }
                        }
                    }
                }

                if (elementIdx == endIdx || hashSeqPair[elementIdx].kmer == SIZE_T_MAX) {
                    break;
                }

                prevSetSize = 0;
                skipByWeightCount = 0;
                prevHash = currKmer;
                prevHashStart = elementIdx;

                // Reset rep info for next hash group (will be updated after Phase 1 swap)
                repSeqKey = hashSeqPair[elementIdx].id;
                repSeqId = repSeqKey;
                repIsReverse = false;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    prevHash = BIT_SET(prevHash, 63);
                    repIsReverse = (BIT_CHECK(hashSeqPair[elementIdx].kmer, 63) == 0);
                    repSeqId = (repIsReverse) ? BIT_CLEAR(repSeqId, 63) : BIT_SET(repSeqId, 63);
                }

                queryLen = hashSeqPair[elementIdx].sl.getSeqLen(hashSeqPair[elementIdx].id);
                repSeq_i_pos = hashSeqPair[elementIdx].pos;

                if (useAdjacentSeq && hashSeqPair[prevHashStart].getAdjacentSeq(0) != UCHAR_MAX) {
                    for (size_t i = 0; i < 6; i++) {
                        subMatPos[i] = subMat->subMatrix[hashSeqPair[prevHashStart].getAdjacentSeq(i)];
                    }
                } else {
                    for (size_t i = 0; i < 6; i++) {
                        subMatPos[i] = NULL;
                    }
                }
            }

            if (hashSeqPair[elementIdx].kmer == SIZE_T_MAX) {
                break;
            }

            prevSetSize++;
            if (prevSetSize > 1 && sequenceWeights != nullptr &&
                sequenceWeights->getWeightById(hashSeqPair[elementIdx].id) > weightThr) {
                skipByWeightCount++;
            }
        }
    }

    if (isSetupCountTable) {
        return 0;
    }

    size_t writePos = localWritePos[0];
    if (writeSeqPair != nullptr) {
        for (int thread = 1; thread < threads; thread++) {
            size_t startIdx = threadOffsets[thread];
            size_t endIdx = localWritePos[thread];

            for (size_t cpid = startIdx; cpid < endIdx; cpid++) {
                writeSeqPair[writePos++] = writeSeqPair[cpid];
            }
        }
        writeSeqPair[writePos].kmer = SIZE_T_MAX;
    } else {
        for (int thread = 1; thread < threads; thread++) {
            size_t startIdx = threadOffsets[thread];
            size_t endIdx = localWritePos[thread];
            for (size_t cpid = startIdx; cpid < endIdx; cpid++) {
                hashSeqPair[writePos++] = hashSeqPair[cpid];
            }
        }
        hashSeqPair[writePos].kmer = SIZE_T_MAX;
    }
    return writePos;
}


// A lane merges its files from every iteration as one sorted run, so it must own the same query ids in
// all of them: the cuts stay query ids and only their spacing changes. Weighting by sequence length
// undoes the skew a length-sorted db otherwise puts on lane 0, and the cuts stay a pure function of the
// db and the thread count, so every iteration recomputes the same ones.
static void computeWeightedQueryCuts(DBReader<DBKeyType> &seqDbr, bool localIdSpace, size_t idSpaceSize,
                                     int threads, std::vector<size_t> &cuts) {
    const size_t BUCKETS = 1u << 16;
    // the cuts are a pure function of the db and the thread count, so the later iterations reuse them
    static std::vector<size_t> cached;
    static size_t cachedSize = 0, cachedSpace = 0;
    static int cachedThreads = 0;
    if (cached.empty() == false && seqDbr.getSize() == cachedSize
        && idSpaceSize == cachedSpace && threads == cachedThreads) {
        cuts = cached;
        return;
    }
    cuts.assign(static_cast<size_t>(threads) + 1, 0);
    cuts[threads] = idSpaceSize;
    if (threads < 2 || idSpaceSize == 0) {
        cached = cuts; cachedSize = seqDbr.getSize(); cachedSpace = idSpaceSize; cachedThreads = threads;
        return;
    }
    std::vector<size_t> bucket(BUCKETS, 0);
    const size_t dbSize = seqDbr.getSize();
#pragma omp parallel for schedule(static) num_threads(threads)
    for (size_t i = 0; i < dbSize; i++) {
        const size_t pos = localIdSpace ? i : static_cast<size_t>(seqDbr.getDbKey(i));
        if (pos >= idSpaceSize) {
            continue;
        }
        const size_t b = pos * BUCKETS / idSpaceSize;
        const size_t w = seqDbr.getSeqLen(i);
#pragma omp atomic
        bucket[b] += w;
    }
    size_t total = 0;
    for (size_t b = 0; b < BUCKETS; b++) {
        const size_t w = bucket[b];
        bucket[b] = total;
        total += w;
    }
    if (total == 0) {
        for (int t = 1; t < threads; t++) {
            cuts[t] = idSpaceSize / static_cast<size_t>(threads) * static_cast<size_t>(t);
        }
        cached = cuts; cachedSize = seqDbr.getSize(); cachedSpace = idSpaceSize; cachedThreads = threads;
        return;
    }
    size_t b = 0;
    for (int t = 1; t < threads; t++) {
        const size_t target = total / static_cast<size_t>(threads) * static_cast<size_t>(t);
        while (b + 1 < BUCKETS && bucket[b + 1] <= target) {
            b++;
        }
        const size_t cut = b * idSpaceSize / BUCKETS;
        cuts[t] = std::max(cut, cuts[t - 1]);
    }
    cached = cuts; cachedSize = seqDbr.getSize(); cachedSpace = idSpaceSize; cachedThreads = threads;
}

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
static void runIteration(
    AssignGroupMask mask, int &iteration, size_t &writePos,
    size_t hashEndRange, const std::string &splitFile,
    KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair,
    KmerPosition<T, false, IncludeSeqLen> *writeSeqPair,
    DBReader<DBKeyType> &seqDbr, Parameters &par,
    BaseMatrix *subMat, uint8_t *countTable,
    SequenceWeights *sequenceWeights,
    std::vector<size_t> &threadOffsets, Timer &timer) {

    timer.reset();
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        writePos = assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
            hashSeqPair, writeSeqPair,
            par.includeOnlyExtendable, par.covMode, par.covThr,
            sequenceWeights, par.weightThr,
            par.threads, threadOffsets, subMat,
            mask, ComputationPhase::Main, countTable);
    } else {
        writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
            hashSeqPair, writeSeqPair,
            par.includeOnlyExtendable, par.covMode, par.covThr,
            sequenceWeights, par.weightThr,
            par.threads, threadOffsets, subMat,
            mask, ComputationPhase::Main, countTable);
    }
    Debug(Debug::INFO) << "Time for assign: " << timer.lap() << "\n";

    Debug(Debug::INFO) << "Sort by rep. sequence ";
    timer.reset();
    if (par.needWriteBuffer) {
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            SORT_PARALLEL(
                writeSeqPair, writeSeqPair + writePos,
                KmerPosition<T, false, IncludeSeqLen>::compareRepSequenceAndIdAndDiagReverse);
        } else {
            SORT_PARALLEL(
                writeSeqPair, writeSeqPair + writePos,
                KmerPosition<T, false, IncludeSeqLen>::compareRepSequenceAndIdAndDiag);
        }
    } else {
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            SORT_PARALLEL(
                hashSeqPair, hashSeqPair + writePos,
                KmerPosition<T, includeAdjacency, IncludeSeqLen>::compareRepSequenceAndIdAndDiagReverse);
        } else {
            SORT_PARALLEL(
                hashSeqPair, hashSeqPair + writePos,
                KmerPosition<T, includeAdjacency, IncludeSeqLen>::compareRepSequenceAndIdAndDiag);
        }
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    if (hashEndRange != SIZE_T_MAX || par.needWriteBuffer) {
        std::vector<size_t> threadQueryOffsets(par.threads + 1);
        // .kmer is a db key in mode 1 and a local id in mode 2, and sparse keys would crowd lane 0
        const size_t idSpaceSize = (par.kmerMatcherMode == Parameters::KMERMATCHER_MODE_LOCAL)
            ? seqDbr.getSize() : (static_cast<size_t>(seqDbr.getLastKey()) + 1);
        std::vector<size_t> queryCuts;
        computeWeightedQueryCuts(seqDbr, par.kmerMatcherMode == Parameters::KMERMATCHER_MODE_LOCAL,
                                 idSpaceSize, par.threads, queryCuts);
        threadQueryOffsets[0] = 0;

        if (par.needWriteBuffer) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (int thread = 1; thread < par.threads; thread++) {
                size_t startqid = queryCuts[thread];
                KmerPosition<T, false, IncludeSeqLen> *it = std::lower_bound(
                    writeSeqPair, writeSeqPair + writePos, startqid,
                    [](const KmerPosition<T, false, IncludeSeqLen> &elem, size_t k) {
                        return BIT_SET(elem.kmer, 63) < BIT_SET(k, 63);
                    });
                threadQueryOffsets[thread] = it - writeSeqPair;
            }
        } else {
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
            for (int thread = 1; thread < par.threads; thread++) {
                size_t startqid = queryCuts[thread];
                KmerPosition<T, includeAdjacency, IncludeSeqLen> *it = std::lower_bound(
                    hashSeqPair, hashSeqPair + writePos, startqid,
                    [](const KmerPosition<T, includeAdjacency, IncludeSeqLen> &elem, size_t k) {
                        return BIT_SET(elem.kmer, 63) < BIT_SET(k, 63);
                    });
                threadQueryOffsets[thread] = it - hashSeqPair;
            }
        }
        threadQueryOffsets[par.threads] = writePos;
        {
            // the merge runs one lane per writer, so a skewed split here is a serial tail there
            size_t minCnt = SIZE_T_MAX, maxCnt = 0;
            for (int thread = 0; thread < par.threads; thread++) {
                const size_t cnt = threadQueryOffsets[thread + 1] - threadQueryOffsets[thread];
                minCnt = std::min(minCnt, cnt);
                maxCnt = std::max(maxCnt, cnt);
            }
            Debug(Debug::INFO) << "Writer split: " << writePos << " entries over " << par.threads
                               << " lanes, min " << minCnt << ", max " << maxCnt << "\n";
        }

        timer.reset();
        if (par.needWriteBuffer) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, false, IncludeSeqLen>(
                    splitFile, writeSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            } else {
                writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T, false, IncludeSeqLen>(
                    splitFile, writeSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            }
        } else {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, T, includeAdjacency, IncludeSeqLen>(
                    splitFile, hashSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            } else {
                writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, T, includeAdjacency, IncludeSeqLen>(
                    splitFile, hashSeqPair, writePos + 1,
                    par.threads, &threadQueryOffsets, iteration);
            }
        }
        Debug(Debug::INFO) << "Time for write: " << timer.lap() << "\n";
    }

    iteration++;
}


template <typename T, bool includeAdjacency, bool IncludeSeqLen>
KmerPosition<T, includeAdjacency, IncludeSeqLen> *doComputation(
    size_t totalKmers, size_t hashStartRange, size_t hashEndRange,
    const std::string &splitFile, AssignGroupMask assignGroupMask,
    ComputationPhase phase, DBReader<DBKeyType> &seqDbr,
    Parameters &par, BaseMatrix *subMat, uint8_t *countTable,
    const std::string *kmerWriteBase = NULL, size_t kmerWriteBucket = 0) {

    KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair =
        initKmerPositionMemory<T, includeAdjacency, IncludeSeqLen>(totalKmers);

    // The write buffer is only used by the assignGroup/runIteration step below, not by the
    // fill step. Allocate it after fill frees its per-thread scratch to reduce peak memory
    // by keeping the scratch and full-size write buffer from coexisting.
    KmerPosition<T, false, IncludeSeqLen> *writeSeqPair = NULL;

    size_t elementsToSort;
    if (kmerWriteBase != NULL) {
        // the partition pass already wrote this split's k-mers, so read the bucket instead of rescanning
        Debug(Debug::INFO) << "Load k-mer bucket ";
        Timer bucketTimer;
        elementsToSort = loadKmerBucket<T, includeAdjacency, IncludeSeqLen>(
            *kmerWriteBase, kmerWriteBucket, par.threads, hashSeqPair, totalKmers, par.compressKmerTmpFiles);
        Debug(Debug::INFO) << bucketTimer.lap() << "\n";
    } else if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        std::pair<size_t, size_t> ret =
            fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
                hashSeqPair, totalKmers, seqDbr, par,
                subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
        par.kmerSize = ret.second;
        Debug(Debug::INFO) << "\nAdjusted k-mer length " << par.kmerSize << "\n";
    } else {
        std::pair<size_t, size_t> ret =
            fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
                hashSeqPair, totalKmers, seqDbr, par,
                subMat, true, hashStartRange, hashEndRange, NULL);
        elementsToSort = ret.first;
    }

    if (hashEndRange == SIZE_T_MAX) {
        seqDbr.unmapData();
    }

    Debug(Debug::INFO) << "Sort kmer ";
    Timer timer;
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        SORT_PARALLEL(
            hashSeqPair, hashSeqPair + elementsToSort,
            KmerPosition<T, includeAdjacency, IncludeSeqLen>::compareRepSequenceAndIdAndPosReverse);
    } else {
        SORT_PARALLEL(
            hashSeqPair, hashSeqPair + elementsToSort,
            KmerPosition<T, includeAdjacency, IncludeSeqLen>::compareRepSequenceAndIdAndPos);
    }
    Debug(Debug::INFO) << timer.lap() << "\n";

    SequenceWeights *sequenceWeights = NULL;
    if (par.PARAM_WEIGHT_FILE.wasSet) {
        sequenceWeights = new SequenceWeights(par.weightFile.c_str());
        if (sequenceWeights != NULL) {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                swapCenterSequence<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
                    hashSeqPair, totalKmers, *sequenceWeights);
            } else {
                swapCenterSequence<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
                    hashSeqPair, totalKmers, *sequenceWeights);
            }
        }
    }

    std::vector<size_t> threadOffsets;
    size_t splitSize = elementsToSort / par.threads;

    threadOffsets.push_back(0);
    for (int thread = 1; thread < par.threads; thread++) {
        if (!par.useParallelism) {
            threadOffsets.push_back(elementsToSort);
            continue;
        }

        size_t prevHash = hashSeqPair[thread * splitSize].kmer;
        if (prevHash == SIZE_T_MAX) {
            for (int i = thread; i < par.threads; i++) {
                threadOffsets.push_back(elementsToSort);
            }
            break;
        }

        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            prevHash = BIT_SET(prevHash, 63);
        }

        bool wasSet = false;
        for (size_t pos = thread * splitSize; pos < elementsToSort; pos++) {
            size_t currKmer = hashSeqPair[pos].kmer;
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                currKmer = BIT_SET(currKmer, 63);
            }
            if (prevHash != currKmer) {
                wasSet = true;
                threadOffsets.push_back(pos);
                break;
            }
        }

        if (wasSet == false) {
            for (int i = thread; i < par.threads; i++) {
                threadOffsets.push_back(elementsToSort);
            }
            break;
        }
    }
    threadOffsets.push_back(elementsToSort);

    if (phase == ComputationPhase::SetupCountTable) {
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
                hashSeqPair, NULL,
                par.includeOnlyExtendable, par.covMode, par.covThr,
                sequenceWeights, par.weightThr,
                par.threads, threadOffsets, subMat,
                AssignGroupFeature::Default,
                ComputationPhase::SetupCountTable, countTable);
        } else {
            assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
                hashSeqPair, NULL,
                par.includeOnlyExtendable, par.covMode, par.covThr,
                sequenceWeights, par.weightThr,
                par.threads, threadOffsets, subMat,
                AssignGroupFeature::Default,
                ComputationPhase::SetupCountTable, countTable);
        }

        delete sequenceWeights;
        delete[] hashSeqPair;
        delete[] writeSeqPair;
        return NULL;
    }

    // Now that the fill step (and its per-thread scratch) is done, allocate the write
    // buffer used by the assignGroup iterations.
    if (par.needWriteBuffer) {
        writeSeqPair = initKmerPositionMemory<T, false, IncludeSeqLen>(totalKmers);
    }

    int iteration = 0;
    size_t writePos = 0;

    runIteration<T, includeAdjacency, IncludeSeqLen>(
        AssignGroupFeature::Default, iteration, writePos,
        hashEndRange, splitFile, hashSeqPair, writeSeqPair,
        seqDbr, par, subMat, countTable,
        sequenceWeights, threadOffsets, timer);

    if (hasFeature(assignGroupMask, AssignGroupFeature::AdjacentSeq)) {
        for (int iter = 0; iter < par.adjIteration; iter++) {
            runIteration<T, includeAdjacency, IncludeSeqLen>(
                AssignGroupFeature::AdjacentSeq, iteration, writePos,
                hashEndRange, splitFile, hashSeqPair, writeSeqPair,
                seqDbr, par, subMat, countTable,
                sequenceWeights, threadOffsets, timer);
        }
    }

    if (hasFeature(assignGroupMask, AssignGroupFeature::CountTable)) {
        for (int iter = 0; iter < par.countTableIteration; iter++) {
            runIteration<T, includeAdjacency, IncludeSeqLen>(
                AssignGroupFeature::CountTable, iteration, writePos,
                hashEndRange, splitFile, hashSeqPair, writeSeqPair,
                seqDbr, par, subMat, countTable,
                sequenceWeights, threadOffsets, timer);
        }
    }

    delete sequenceWeights;
    if (hashEndRange != SIZE_T_MAX || par.needWriteBuffer) {
        delete[] hashSeqPair;
        hashSeqPair = NULL;
    }
    delete[] writeSeqPair;
    writeSeqPair = NULL;

    return hashSeqPair;
}

void setLinearFilterDefault(Parameters *p) {
    p->covThr = 0.8;
    p->maskMode = 0;
    p->spacedKmer = 0;
    p->kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
    p->alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE, 5));
    p->kmersPerSequence = Parameters::CLUST_LINEAR_KMER_PER_SEQ;
}


size_t computeKmerCount(DBReader<DBKeyType> &reader, size_t KMER_SIZE, size_t chooseTopKmer, float chooseTopKmerScale) {
    size_t totalKmers = 0;
    for(size_t id = 0; id < reader.getSize(); id++ ){
        int seqLen = static_cast<int>(reader.getSeqLen(id));
        int kmerAdjustedSeqLen = std::max(1, seqLen  - static_cast<int>(KMER_SIZE ) + 2) ;
        totalKmers += std::min(kmerAdjustedSeqLen, static_cast<int>( chooseTopKmer + (chooseTopKmerScale * seqLen)));
    }
    return totalKmers;
}

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
size_t computeMemoryNeededLinearfilter(size_t totalKmer) {
    return sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) * totalKmer;
}


// the histogram mode never touches the record layout, so one scan serves every instantiation
template <typename T, bool includeAdjacency, bool IncludeSeqLen>
size_t *fillHashDist(Parameters &par, BaseMatrix *subMat, DBReader<DBKeyType> &seqDbr) {
    size_t *hashDist = new(std::nothrow) size_t[USHRT_MAX + 1];
    Util::checkAllocation(hashDist, "Can not allocate hashDist memory");
    memset(hashDist, 0, sizeof(size_t) * (USHRT_MAX + 1));
    if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
            NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
    } else {
        fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
            NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, hashDist);
    }
    seqDbr.dropCacheAll();
    return hashDist;
}

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
std::vector<std::pair<size_t, size_t>> setupCountTable(
    Parameters &par,
    BaseMatrix *subMat,
    DBReader<DBKeyType> &seqDbr,
    size_t splits,
    size_t totalKmersPerSplit,
    size_t totalKmersCountTable,
    size_t **sharedHashDist
) {
    std::vector<std::pair<size_t, size_t>> hashRanges;
    Debug(Debug::INFO) << "Initiating count table\n";
    size_t *hashDist = fillHashDist<T, includeAdjacency, IncludeSeqLen>(par, subMat, seqDbr);

    if (splits > 1) {
        Debug(Debug::INFO) << "Not enough memory to process at once need to split for initiating count table\n";
        size_t maxBucketSize = 0;
        for(size_t i = 0; i < (USHRT_MAX+1); i++) {
            if(maxBucketSize < hashDist[i]){
                maxBucketSize = hashDist[i];
            }
        }
        if(maxBucketSize > totalKmersPerSplit){
            Debug(Debug::INFO) << "Not enough memory to run the kmermatcher. Minimum is at least " << maxBucketSize* sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) << " bytes\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // find hashRange [0, totalKmersCountTable]
    size_t currBucketSize = 0;
    size_t currBucketStart = 0;
    size_t totalBucketSize = 0;
    for(size_t i = 0; i < (USHRT_MAX+1); i++){
        // if bucketsize exceeds subsampled countable range then break
        if(totalBucketSize+hashDist[i] >= totalKmersCountTable){
            hashRanges.emplace_back(currBucketStart, i - 1);
            break;
        }
        // if bucketsize exceeds memory limits -> multiple hashranges and reset startpos
        if(currBucketSize+hashDist[i] >= totalKmersPerSplit){
            hashRanges.emplace_back(currBucketStart, i - 1);
            currBucketSize = 0;
            currBucketStart = i;
        }
        currBucketSize+=hashDist[i];
        totalBucketSize+=hashDist[i];
    }
    if (sharedHashDist != NULL) {
        *sharedHashDist = hashDist;
    } else {
        delete [] hashDist;
    }
    return hashRanges;
}


template <typename T, bool includeAdjacency, bool IncludeSeqLen>
int kmermatcherInner(Parameters& par, DBReader<DBKeyType>& seqDbr) {
    int querySeqType = seqDbr.getDbtype();
    const bool localIds = par.kmerMatcherMode == Parameters::KMERMATCHER_MODE_LOCAL;
    const size_t idSpaceSize = localIds ? seqDbr.getSize() : seqDbr.getLastKey() + 1;
    if (localIds && par.weightFile.empty() == false) {
        Debug(Debug::ERROR) << "--kmermatcher-mode 2 cannot be combined with --weights yet; "
                            << "weights are keyed by DB key.\n";
        EXIT(EXIT_FAILURE);
    }
    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    }else {
        if (par.alphabetSize.values.aminoacid() == 21) {
            subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
        } else {
            SubstitutionMatrix sMat(par.scoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2f);
            subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2num, sMat.num2aa, sMat.alphabetSize, par.alphabetSize.values.aminoacid(), 2.0);
        }
    }

    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);

    Debug(Debug::INFO) << "\n";
    float kmersPerSequenceScale = (Parameters::isEqualDbtype(querySeqType, Parameters::DBTYPE_NUCLEOTIDES)) ?
                                        par.kmersPerSequenceScale.values.nucleotide() : par.kmersPerSequenceScale.values.aminoacid();
    size_t totalKmers = computeKmerCount(seqDbr, par.kmerSize, par.kmersPerSequence, kmersPerSequenceScale);
    // the write buffer carries no adjacency, so sizing it separately keeps splits from being over-counted
    size_t totalSizeNeeded = computeMemoryNeededLinearfilter<T, includeAdjacency, IncludeSeqLen>(totalKmers)
                             + (par.needWriteBuffer ? computeMemoryNeededLinearfilter<T, false, IncludeSeqLen>(totalKmers) : 0);

    // --split-memory-limit sizes the main k-mer/write buffers (hashSeqPair [+ writeSeqPair]).
    // Reserve the tables that stay resident next to them for the whole split loop first,
    // then keep 5% headroom for transient allocations:
    //   seqLenTable : per-sequence length lookup (skipped when lengths are stored inline)
    //   countTable  : per-sequence k-mer counts, used by the center-swapping iterations
    size_t seqLenTableMemory = (!IncludeSeqLen) ? (idSpaceSize + 1) * sizeof(T) : 0;
    size_t countTableMemory = 0;
    if (par.includeCountTable) {
        if (idSpaceSize > std::numeric_limits<size_t>::max() / sizeof(uint8_t)) {
            Debug(Debug::ERROR) << "Count table is too large to allocate.\n";
            EXIT(EXIT_FAILURE);
        }
        countTableMemory = idSpaceSize * sizeof(uint8_t);
    }
    size_t fixedMemory = seqLenTableMemory + countTableMemory;
    if (fixedMemory >= memoryLimit) {
        Debug(Debug::ERROR) << "Not enough memory to run the kmermatcher. Memory limit: " << memoryLimit
                            << " bytes, resident tables (seqkey_to_len + count table) need: " << fixedMemory << " bytes.\n"
                            << "Raise --split-memory-limit.\n";
        EXIT(EXIT_FAILURE);
    }
    size_t splitMemoryLimit = static_cast<size_t>(static_cast<double>(memoryLimit - fixedMemory) * 0.95);
    if (splitMemoryLimit == 0) {
        Debug(Debug::ERROR) << "Not enough memory to run the kmermatcher after reserving fixed memory\n";
        EXIT(EXIT_FAILURE);
    }

    size_t splits = std::max(static_cast<size_t>(1), static_cast<size_t>(std::ceil(static_cast<float>(totalSizeNeeded) / splitMemoryLimit)));
    size_t totalKmersPerSplit = std::max(
                                    static_cast<size_t>(1024 + 1),
                                    static_cast<size_t>(
                                        std::min(totalSizeNeeded, splitMemoryLimit) /
                                        (sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) +
                                        (par.needWriteBuffer ? sizeof(KmerPosition<T, false, IncludeSeqLen>) : 0))
                                    ) + 1
                                );

    // one partition has no bucket to reuse, and splitting to make one buys a second pass over the db
    if (par.kmerWriteToDisk && splits == 1) {
        par.kmerWriteToDisk = false;
        Debug(Debug::INFO) << "--kmer-write-to-disk has nothing to spill at one k-mer partition, ignoring it\n";
    }

    if (!IncludeSeqLen) {
        T * seqkey_to_len = new(std::nothrow) T[idSpaceSize+1];
        Util::checkAllocation(seqkey_to_len, "Can not allocate seqkey_to_len memory");
        memset(seqkey_to_len, 0, sizeof(T)*(idSpaceSize+1));
#pragma omp parallel
        {
#pragma omp for schedule(dynamic, 1000)
            for (size_t id = 0; id < seqDbr.getSize(); id++) {
                const DBKeyType resultId = localIds ? static_cast<DBKeyType>(id) : seqDbr.getDbKey(id);
                seqkey_to_len[resultId] = static_cast<T>(seqDbr.getSeqLen(id));
            }
        }
        SeqLenData<T, false>::seqkey_to_len = seqkey_to_len;
    }

    std::vector<uint8_t> countTable;
    // --adjust-kmer-len rewrites par.kmerSize between the two setups, so only share the scan without it
    size_t *sharedHashDist = NULL;
    if (par.includeCountTable) {
        // countTable is already reserved in fixedMemory above; its fill step allocates its own
        // k-mer buffers, which fit in the same splitMemoryLimit (seqkey_to_len + countTable stay
        // resident during the fill too).
        countTable.assign(idSpaceSize, 0);
        size_t countTableTotalKmers = static_cast<size_t>(totalKmers * par.countTableScale);
        // hashSeqPair + writeSeqPair for the count-table fill
        size_t countTableTotalSizeNeeded = computeMemoryNeededLinearfilter<T, false, false>(countTableTotalKmers) * 2;

        size_t countTableSplits = std::max(static_cast<size_t>(1), static_cast<size_t>(
            std::ceil(static_cast<double>(countTableTotalSizeNeeded) / splitMemoryLimit)
        ));

        size_t countTableKmersPerSplit = std::max(
            static_cast<size_t>(1024 + 1),
            static_cast<size_t>(
                std::min(countTableTotalSizeNeeded, splitMemoryLimit) /
                (sizeof(KmerPosition<T, false, false>) * 2)
            ) + 1
        );

        std::vector<std::pair<size_t, size_t>> countTableHashRanges;
        countTableHashRanges = setupCountTable<T, false, false>(
            par, subMat, seqDbr,
            countTableSplits, countTableKmersPerSplit,
            static_cast<size_t>(totalKmers * par.countTableScale),
            &sharedHashDist
        );
        for (size_t split = 0; split < countTableHashRanges.size(); split++) {
            Debug(Debug::INFO) << "Fill count table for " << (split + 1) << " split\n";

            doComputation<T, false, false>(
                countTableKmersPerSplit,
                countTableHashRanges[split].first,
                countTableHashRanges[split].second,
                "COUNT_TABLE",
                AssignGroupFeature::Default,
                ComputationPhase::SetupCountTable,
                seqDbr, par, subMat, countTable.data());
        }
    }
    // set assigngroup matching mode
    AssignGroupMask assignGroupMask = AssignGroupFeature::Default;
    if (par.includeAdjacency) {
        assignGroupMask |= AssignGroupFeature::AdjacentSeq;
    }
    if (par.includeCountTable) {
        assignGroupMask |= AssignGroupFeature::CountTable;
    }

    std::vector<std::pair<size_t, size_t>> hashRanges = setupKmerSplits<T, includeAdjacency, IncludeSeqLen>(
        par, subMat, seqDbr, totalKmersPerSplit, splits, par.adjustKmerLength ? NULL : sharedHashDist);
    delete [] sharedHashDist;
    if(splits > 1){
        Debug(Debug::INFO) << "Process file into " << hashRanges.size() << " parts\n";
    }
    std::vector<std::string> splitFiles;
    KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair = NULL;

    std::string kmerWriteBase;
    bool useKmerWrite = false;

    size_t mpiRank = 0;
#ifdef HAVE_MPI
    splits = hashRanges.size();
    size_t fromSplit = 0;
    size_t splitCount = 1;
    mpiRank = MMseqsMPI::rank;
    unsigned int * splitCntPerProc = new unsigned int[MMseqsMPI::numProc];
    memset(splitCntPerProc, 0, sizeof(unsigned int) * MMseqsMPI::numProc);
    for(size_t i = 0; i < splits; i++){
        splitCntPerProc[i % MMseqsMPI::numProc] += 1;
    }
    for(int i = 0; i < MMseqsMPI::rank; i++){
        fromSplit += splitCntPerProc[i];
    }
    splitCount = splitCntPerProc[MMseqsMPI::rank];
    delete[] splitCntPerProc;

    for(size_t split = fromSplit; split < fromSplit+splitCount; split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        hashSeqPair = doComputation<T, includeAdjacency, IncludeSeqLen>(
                        totalKmers, hashRanges[split].first, hashRanges[split].second, splitFileName,
                        assignGroupMask, ComputationPhase::Main,
                        seqDbr, par, subMat, countTable.empty() ? NULL : countTable.data());
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(mpiRank == 0){
        for(size_t split = 0; split < splits; split++) {
            std::string splitFileName = par.db2 + "_split_" +SSTR(split);
            splitFiles.push_back(splitFileName);
        }
    }
#else
    // --kmer-write-to-disk: extract once into per-split buckets so the split loop never rescans
    if ((splits > 1) && par.kmerWriteToDisk) {
        useKmerWrite = true;
        kmerWriteBase = par.db2 + "_bucket";
        size_t openFiles = static_cast<size_t>(par.threads) * hashRanges.size();
        long openMax = sysconf(_SC_OPEN_MAX);
        size_t fdBudget = (openMax > 128) ? static_cast<size_t>(openMax - 64) : 64;
        if (openFiles >= fdBudget) {
            // FileUtil::fixRlimitNoFile already raised the soft limit, so only the hard limit is actionable
            struct rlimit lim;
            unsigned long long hard = (getrlimit(RLIMIT_NOFILE, &lim) == 0) ? (unsigned long long) lim.rlim_max : 0;
            Debug(Debug::ERROR) << "The k-mer buckets need " << openFiles << " open bucket files (threads x splits) "
                                << "but the open-file limit is " << openMax << " (hard limit " << hard
                                << "). Lower --threads or raise the hard limit (ulimit -Hn).\n";
            EXIT(EXIT_FAILURE);
        }
        unsigned int *hashToBucket = buildHashToBucketLookup(hashRanges);
        Debug(Debug::INFO) << "Partition k-mers into " << hashRanges.size() << " buckets\n";
        // all buckets coexist until their split consumes them, so this is the peak tmp usage
        size_t bucketBytes = totalKmers * sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>);
        struct statvfs bucketVfs;
        size_t bucketFsFree = (statvfs(FileUtil::dirName(kmerWriteBase).c_str(), &bucketVfs) == 0)
                                  ? static_cast<size_t>(bucketVfs.f_bavail) * bucketVfs.f_frsize : SIZE_T_MAX;
        Debug(Debug::INFO) << "K-mer buckets hold up to " << (bucketBytes >> 30)
                           << " GB before compression\n";
        if (bucketFsFree != SIZE_T_MAX && bucketBytes > bucketFsFree) {
            Debug(Debug::WARNING) << "The tmp filesystem has " << (bucketFsFree >> 30)
                                  << " GB free for the k-mer buckets; compression may still fit them\n";
        }
        KmerPartitionSink<T, includeAdjacency, IncludeSeqLen> sink(kmerWriteBase, par.threads, hashRanges.size(), hashToBucket, par.compressKmerTmpFiles);
        if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
            std::pair<size_t, size_t> ret = fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
                NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, NULL, &sink);
            par.kmerSize = ret.second;
        } else {
            fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
                NULL, SIZE_T_MAX, seqDbr, par, subMat, true, 0, SIZE_T_MAX, NULL, &sink);
        }
        sink.finish();
        // last use of the sequence bodies, so hand their cache to the k-mer arrays before the split loop
        seqDbr.dropCacheAll();
        seqDbr.unmapData();
        delete[] hashToBucket;
    }
    for(size_t split = 0; split < hashRanges.size(); split++) {
        std::string splitFileName = par.db2 + "_split_" +SSTR(split);
        Debug(Debug::INFO) << "Generate k-mers list for " << (split+1) <<" split\n";

        hashSeqPair = doComputation<T, includeAdjacency, IncludeSeqLen>(
                    totalKmersPerSplit, hashRanges[split].first, hashRanges[split].second, splitFileName,
                    assignGroupMask, ComputationPhase::Main,
                    seqDbr, par, subMat, countTable.empty() ? NULL : countTable.data(),
                    useKmerWrite ? &kmerWriteBase : NULL, split);

        splitFiles.push_back(splitFileName);
    }
#endif
    // the merge reads neither table, so give the pages back before repSequence and the merge buffers
    countTable.clear();
    countTable.shrink_to_fit();
    if (!IncludeSeqLen && SeqLenData<T, false>::seqkey_to_len != NULL) {
        delete[] SeqLenData<T, false>::seqkey_to_len;
        SeqLenData<T, false>::seqkey_to_len = NULL;
    }
    if(mpiRank == 0){
        std::vector<char> repSequence(idSpaceSize);
        std::fill(repSequence.begin(), repSequence.end(), false);
        const int resultDbtype = localIds
            ? Parameters::DBTYPE_PREFILTER_LOCAL_RES
            : ((Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES))
                ? Parameters::DBTYPE_PREFILTER_REV_RES : Parameters::DBTYPE_PREFILTER_RES);
        DBWriter dbw(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, resultDbtype);
        dbw.open();

        Timer timer;
        if (splits > 1 || par.needWriteBuffer) {
            int maxIter = 1;
            if (par.includeAdjacency) {
                maxIter += par.adjIteration;
            }
            if (par.includeCountTable) {
                maxIter += par.countTableIteration;
            }

            seqDbr.unmapData();
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev>(dbw, splitFiles, repSequence,
                                                                                         par.threads, maxIter, &seqDbr, localIds);
            } else {
                mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry>(dbw, splitFiles, repSequence,
                                                                                      par.threads, maxIter, &seqDbr, localIds);
            }

        } else {
            if (Parameters::isEqualDbtype(seqDbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
                writeKmerMatcherResult<Parameters::DBTYPE_NUCLEOTIDES, T, includeAdjacency, IncludeSeqLen>(
                    dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1, &seqDbr, localIds);
            } else {
                writeKmerMatcherResult<Parameters::DBTYPE_AMINO_ACIDS, T, includeAdjacency, IncludeSeqLen>(
                    dbw, hashSeqPair, totalKmersPerSplit, repSequence, 1, &seqDbr, localIds);
            }
        }
        Debug(Debug::INFO) << "Time for fill: " << timer.lap() << "\n";

        // one entry per ungrouped sequence, and no consumer needs write order, so static is fine
#pragma omp parallel
        {
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
#pragma omp for schedule(static)
            for (size_t id = 0; id < seqDbr.getSize(); id++) {
                char buffer[100];
                const DBKeyType dbKey = seqDbr.getDbKey(id);
                const DBKeyType resultId = localIds ? static_cast<DBKeyType>(id) : dbKey;
                if (repSequence[resultId] == false) {
                    hit_t h;
                    h.prefScore = 0;
                    h.diagonal = 0;
                    h.seqId = resultId;
                    int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                    dbw.writeData(buffer, len, dbKey, thread_idx);
                }
            }
        }
        dbw.close(false, true);
    }
    delete subMat;
    if(hashSeqPair){
        delete [] hashSeqPair;
    }
    return EXIT_SUCCESS;
}

template <typename T, bool includeAdjacency, bool IncludeSeqLen>
std::vector<std::pair<size_t, size_t>> setupKmerSplits(Parameters &par, BaseMatrix * subMat, DBReader<DBKeyType> &seqDbr, size_t totalKmers, size_t splits, size_t *reuseHashDist){
    std::vector<std::pair<size_t, size_t>> hashRanges;
    if (splits > 1) {
        Debug(Debug::INFO) << "Not enough memory to process at once need to split\n";
        const bool owned = (reuseHashDist == NULL);
        size_t *hashDist = owned ? fillHashDist<T, includeAdjacency, IncludeSeqLen>(par, subMat, seqDbr)
                                 : reuseHashDist;
        size_t maxBucketSize = 0;
        for(size_t i = 0; i < (USHRT_MAX+1); i++) {
            if(maxBucketSize < hashDist[i]){
                maxBucketSize = hashDist[i];
            }
        }
        if(maxBucketSize > totalKmers){
            Debug(Debug::INFO) << "Not enough memory to run the kmermatcher. Minimum is at least " << maxBucketSize* sizeof(KmerPosition<T, includeAdjacency, IncludeSeqLen>) << " bytes\n";
            EXIT(EXIT_FAILURE);
        }
        size_t currBucketSize = 0;
        size_t currBucketStart = 0;
        for(size_t i = 0; i < (USHRT_MAX+1); i++){
            if(currBucketSize+hashDist[i] >= totalKmers){
                hashRanges.emplace_back(currBucketStart, i - 1);
                currBucketSize = 0;
                currBucketStart = i;
            }
            currBucketSize+=hashDist[i];
        }
        hashRanges.emplace_back(currBucketStart, (USHRT_MAX+1));
        if (owned) { delete [] hashDist; }
    }else{
        hashRanges.emplace_back(0, SIZE_T_MAX);
    }
    return hashRanges;
}

int kmermatcher(int argc, const char **argv, const Command &command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    setLinearFilterDefault(&par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    DBReader<DBKeyType> seqDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    // releasing the sequence cache after the partition pass needs a descriptor that outlives the mapping
    seqDbr.setIoCacheAdvice(true);
    // the whole-db scans below are sequential, which only NOSORT plus this hint can tell the reader
    seqDbr.open(DBReader<DBKeyType>::NOSORT);
    // NOSORT is key order and createdb writes offsets monotone in it, so ask for readahead too
    seqDbr.setSequentialAdvice();
    int querySeqType = seqDbr.getDbtype();

    setKmerLengthAndAlphabet(par, seqDbr.getAminoAcidDBSize(), querySeqType);
    std::vector<MMseqsParameter *> *params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    Debug(Debug::INFO) << "Database size: " << seqDbr.getSize() << " type: " << seqDbr.getDbTypeName() << "\n";

    
    if (par.includeCountTable || par.includeAdjacency) {
        par.needWriteBuffer = true;
        par.useParallelism = true;
    } else {
        // if user turns on parallelism, it lead to using writebuffer with extra memory usage
        // linclust 1 setting: useParallelism = false
        par.needWriteBuffer = par.useParallelism;
    }

    if (par.linclustVersion == 1) {
        par.needWriteBuffer = false;
        // in-place lanes share their boundary element, so linclust1 stays single-lane
        par.useParallelism = false;
        par.includeCountTable = false;
        par.includeAdjacency = false;
        par.adjIteration = 0;
        par.countTableIteration = 0;
    } else if (par.linclustVersion == 2){
        if (par.includeAdjacency && par.adjIteration == 0) {
            Debug(Debug::ERROR) << "Adjacent Iteration must be greater than 0 when include Adjacent is true\n";
            EXIT(EXIT_FAILURE);
        }
        if (par.includeCountTable && par.countTableIteration == 0) {
            Debug(Debug::ERROR) << "CountTable Iteration must be greater than 0 when include CountTable is true\n";
            EXIT(EXIT_FAILURE);
        }
    }
    
    if (seqDbr.getMaxSeqLen() < SHRT_MAX) {
        if (par.includeAdjacency) { 
            kmermatcherInner<short, true, false>(par, seqDbr);
        }
        else {
            kmermatcherInner<short, false, false>(par, seqDbr);
        }
    }
    else {
        if (par.includeAdjacency) {
            kmermatcherInner<int, true, false>(par, seqDbr);
        }
        else {
            kmermatcherInner<int, false, false>(par, seqDbr);
        }
    }

    seqDbr.close();

    return EXIT_SUCCESS;
}

template <int TYPE, typename T, bool includeAdjacency, bool IncludeSeqLen>
void writeKmerMatcherResult(DBWriter & dbw,
                            KmerPosition<T, includeAdjacency, IncludeSeqLen> *hashSeqPair, size_t totalKmers,
                            std::vector<char> &repSequence, size_t threads,
                            DBReader<DBKeyType> *sequenceDbr, bool localIds) {
    if (localIds && sequenceDbr == NULL) {
        Debug(Debug::ERROR) << "Local-ID kmermatcher output needs its source sequence DB.\n";
        EXIT(EXIT_FAILURE);
    }
    std::vector<size_t> threadOffsets;
    size_t splitSize = totalKmers/threads;
    threadOffsets.push_back(0);
    for(size_t thread = 1; thread < threads; thread++){
        size_t kmer = hashSeqPair[thread*splitSize].kmer;
        size_t repSeqId = static_cast<size_t>(kmer);
        repSeqId=BIT_SET(repSeqId, 63);
        bool wasSet = false;
        for(size_t pos = thread*splitSize; pos < totalKmers; pos++){
            size_t currSeqId = hashSeqPair[pos].kmer;
            currSeqId=BIT_SET(currSeqId, 63);
            if(repSeqId != currSeqId){
                wasSet = true;
                threadOffsets.push_back(pos);
                break;
            }
        }
        if(wasSet == false){
            threadOffsets.push_back(totalKmers);
        }
    }
    threadOffsets.push_back(totalKmers);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for(size_t thread = 0; thread < threads; thread++){
        std::string prefResultsOutString;
        prefResultsOutString.reserve(100000000);
        char buffer[100];
        size_t lastTargetId = SIZE_T_MAX;
        // counts every hit line a lane writes, so a 32 bit counter can wrap back to zero at 10^10
        size_t writeSets = 0;
        size_t kmerPos=0;
        size_t repSeqId = SIZE_T_MAX;
        for(kmerPos = threadOffsets[thread]; kmerPos < threadOffsets[thread+1] && hashSeqPair[kmerPos].kmer != SIZE_T_MAX; kmerPos++){
            size_t currKmer = hashSeqPair[kmerPos].kmer;
            int reverMask = 0;
            if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                reverMask  = BIT_CHECK(currKmer, 63)==false;
                currKmer = BIT_CLEAR(currKmer, 63);
            }
            if(repSeqId != currKmer) {
                if (writeSets > 0) {
                    repSequence[repSeqId] = true;
                    const DBKeyType outputKey = localIds ? sequenceDbr->getDbKey(repSeqId) : repSeqId;
                    dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), outputKey, thread);
                }else{
                    if(repSeqId != SIZE_T_MAX) {
                        repSequence[repSeqId] = false;
                    }
                }
                lastTargetId = SIZE_T_MAX;
                prefResultsOutString.clear();
                repSeqId = currKmer;
                hit_t h;
                h.seqId = repSeqId;
                h.prefScore = 0;
                h.diagonal = 0;
                int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                prefResultsOutString.append(buffer, len);
            }
            DBKeyType targetId = hashSeqPair[kmerPos].id;
            T diagonal = hashSeqPair[kmerPos].pos;
            size_t kmerOffset = 0;
            T prevDiagonal = diagonal;
            size_t maxDiagonal = 0;
            size_t diagonalCnt = 0;
            size_t topScore =0;
            int bestReverMask = reverMask;
            while(lastTargetId != targetId
                  && kmerPos+kmerOffset < threadOffsets[thread+1]
                  && hashSeqPair[kmerPos+kmerOffset].kmer == repSeqId
                  && hashSeqPair[kmerPos+kmerOffset].id == targetId){
                if(prevDiagonal == hashSeqPair[kmerPos+kmerOffset].pos){
                    diagonalCnt++;
                }else{
                    diagonalCnt = 1;
                }
                if(diagonalCnt >= maxDiagonal){
                    diagonal = hashSeqPair[kmerPos+kmerOffset].pos;
                    maxDiagonal = diagonalCnt;
                    if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
                        bestReverMask = BIT_CHECK(hashSeqPair[kmerPos+kmerOffset].kmer, 63) == false;
                    }
                }
                prevDiagonal = hashSeqPair[kmerPos+kmerOffset].pos;
                kmerOffset++;
                topScore++;
            }
            if(targetId == repSeqId || lastTargetId == targetId){
                lastTargetId = targetId;
                continue;
            }
            hit_t h;
            h.seqId = targetId;
            h.prefScore = (bestReverMask) ? -topScore : topScore;
            h.diagonal = diagonal;
            int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
            prefResultsOutString.append(buffer, len);
            lastTargetId = targetId;
            writeSets++;
        }
        if (writeSets > 0) {
            repSequence[repSeqId] = true;
            const DBKeyType outputKey = localIds ? sequenceDbr->getDbKey(repSeqId) : repSeqId;
            dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), outputKey, thread);
        }else{
            if(repSeqId != SIZE_T_MAX) {
                repSequence[repSeqId] = false;
            }
        }
    }
}

static const size_t KMER_TMP_ZSTD_INPUT_BUFFER_SIZE = 65536;
static const size_t KMER_TMP_ZSTD_OUTPUT_BUFFER_SIZE = 65536;
static const int KMER_TMP_ZSTD_COMPRESSION_LEVEL = 2;
static const size_t KMER_MERGE_RESULT_BUFFER_RESERVE = 1024 * 1024;
static const size_t KMER_MERGE_RESULT_BUFFER_MAX_KEEP = 16 * 1024 * 1024;
static const size_t KMER_MERGE_MIN_MEMORY_HEADROOM = 1024ull * 1024ull * 1024ull;

static std::string kmerTmpFileName(const std::string &tmpFile, int iteration, int threadIdx, int codec) {
    return tmpFile + "_iter_" + std::to_string(iteration) + "_thread_" + std::to_string(threadIdx)
           + kmerTmpCodecSuffix(codec);
}

// a rerun may inherit files written under a different codec, so try the requested one and then the rest
static std::string existingKmerTmpFileName(const std::string &baseName, int preferredCodec) {
    static const int codecs[] = { KMER_TMP_CODEC_NONE, KMER_TMP_CODEC_ZSTD };
    const std::string preferred = baseName + kmerTmpCodecSuffix(preferredCodec);
    if (FileUtil::fileExists(preferred.c_str())) {
        return preferred;
    }
    for (size_t i = 0; i < sizeof(codecs) / sizeof(codecs[0]); i++) {
        const std::string name = baseName + kmerTmpCodecSuffix(codecs[i]);
        if (FileUtil::fileExists(name.c_str())) {
            return name;
        }
    }
    return "";
}

static void removeKmerTmpFileIfExists(const std::string &fileName) {
    if (FileUtil::fileExists(fileName.c_str())) {
        FileUtil::remove(fileName.c_str());
    }
}

// the fs counters tell a full disk, a full quota and a full inode table apart from the log alone
static void printKmerTmpFsState(const std::string &fileName) {
    struct statvfs vfs;
    const std::string dir = FileUtil::dirName(fileName);
    if (statvfs(dir.c_str(), &vfs) == 0) {
        Debug(Debug::ERROR) << "Filesystem of " << dir << " has "
                            << static_cast<uint64_t>(vfs.f_bavail) * vfs.f_frsize << " bytes and "
                            << static_cast<uint64_t>(vfs.f_favail) << " inodes available to this user\n";
    }
}

static void kmerTmpIoErrorDie(const char *op, const std::string &fileName) {
    const int errsv = errno;
    Debug(Debug::ERROR) << "Can not " << op << " " << fileName << ". errno=" << errsv
                        << " (" << strerror(errsv) << ")\n";
    printKmerTmpFsState(fileName);
    EXIT(EXIT_FAILURE);
}

static FILE *openKmerTmpFileForOverwriteOrDie(const std::string &fileName, const char *mode) {
    FILE *file = fopen(fileName.c_str(), mode);
    if (file == NULL) {
        kmerTmpIoErrorDie("open", fileName);
    }
    return file;
}

static void writeAllOrDie(FILE *file, const void *data, size_t dataSize, const std::string &fileName) {
    if (dataSize == 0) {
        return;
    }
    size_t written = fwrite(data, sizeof(char), dataSize, file);
    if (written != dataSize) {
        const int errsv = errno;
        Debug(Debug::ERROR) << "Can not write to file " << fileName << ": wrote " << written
                            << " of " << dataSize << " bytes. errno=" << errsv
                            << " (" << strerror(errsv) << ")\n";
        printKmerTmpFsState(fileName);
        EXIT(EXIT_FAILURE);
    }
}

class ZstdKmerTmpFileWriter {
public:
    ZstdKmerTmpFileWriter(const std::string &fileName, int compressionLevel)
        : fileName(fileName), file(NULL), cstream(NULL), outBuffer(NULL), closed(true) {
        file = openKmerTmpFileForOverwriteOrDie(fileName, "wb");
        cstream = ZSTD_createCStream();
        if (cstream == NULL) {
            Debug(Debug::ERROR) << "ZSTD_createCStream() failed for " << fileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        outBuffer = static_cast<char *>(malloc(KMER_TMP_ZSTD_OUTPUT_BUFFER_SIZE));
        if (outBuffer == NULL) {
            Debug(Debug::ERROR) << "Cannot allocate zstd output buffer for " << fileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (compressionLevel < 1) {
            compressionLevel = 1;
        }
        size_t ret = ZSTD_initCStream(cstream, compressionLevel);
        if (ZSTD_isError(ret)) {
            Debug(Debug::ERROR) << "ZSTD_initCStream() error for " << fileName << ". Error "
                                << ZSTD_getErrorName(ret) << "\n";
            EXIT(EXIT_FAILURE);
        }
        closed = false;
    }

    ~ZstdKmerTmpFileWriter() {
        if (closed == false) {
            close();
        }
        if (outBuffer != NULL) {
            free(outBuffer);
        }
        if (cstream != NULL) {
            ZSTD_freeCStream(cstream);
        }
    }

    void write(const void *data, size_t dataSize) {
        ZSTD_inBuffer input = { data, dataSize, 0 };
        while (input.pos < input.size) {
            ZSTD_outBuffer output = { outBuffer, KMER_TMP_ZSTD_OUTPUT_BUFFER_SIZE, 0 };
            size_t ret = ZSTD_compressStream(cstream, &output, &input);
            if (ZSTD_isError(ret)) {
                Debug(Debug::ERROR) << "ZSTD_compressStream() error for " << fileName << ". Error "
                                    << ZSTD_getErrorName(ret) << "\n";
                EXIT(EXIT_FAILURE);
            }
            writeAllOrDie(file, outBuffer, output.pos, fileName);
        }
    }

    void close() {
        if (closed) {
            return;
        }
        size_t remainingToFlush = 0;
        do {
            ZSTD_outBuffer output = { outBuffer, KMER_TMP_ZSTD_OUTPUT_BUFFER_SIZE, 0 };
            remainingToFlush = ZSTD_endStream(cstream, &output);
            if (ZSTD_isError(remainingToFlush)) {
                Debug(Debug::ERROR) << "ZSTD_endStream() error for " << fileName << ". Error "
                                    << ZSTD_getErrorName(remainingToFlush) << "\n";
                EXIT(EXIT_FAILURE);
            }
            writeAllOrDie(file, outBuffer, output.pos, fileName);
        } while (remainingToFlush != 0);

        if (fclose(file) != 0) {
            kmerTmpIoErrorDie("close", fileName);
        }
        file = NULL;
        closed = true;
    }

private:
    std::string fileName;
    FILE *file;
    ZSTD_CStream *cstream;
    char *outBuffer;
    bool closed;
};

// same shape as ZstdKmerTmpFileWriter so the two codecs stay interchangeable at the call sites

// wraps the raw FILE* spill writer so its error paths can name the file, like the compressed writers
struct RawKmerTmpFileWriter {
    std::string fileName;
    FILE *file;
    RawKmerTmpFileWriter(const std::string &fileName)
        : fileName(fileName), file(openKmerTmpFileForOverwriteOrDie(fileName, "wb")) {}
    void write(const void *data, size_t dataSize) {
        writeAllOrDie(file, data, dataSize, fileName);
    }
    // an unchecked fclose here would truncate the buffered tail silently on a full filesystem
    void close() {
        if (file != NULL && fclose(file) != 0) {
            kmerTmpIoErrorDie("close", fileName);
        }
        file = NULL;
    }
};

// bucket file io wrappers, declared near the top and defined here after the stream writers
static void *bucketWriterOpen(const std::string &fileName, int codec, int zstdLevel) {
    if (codec == KMER_TMP_CODEC_ZSTD) {
        return new ZstdKmerTmpFileWriter(fileName, zstdLevel);
    }
    return new RawKmerTmpFileWriter(fileName);
}

static void bucketWriterAppend(void *writer, int codec, const void *data, size_t byteSize) {
    if (byteSize == 0) {
        return;
    }
    if (codec == KMER_TMP_CODEC_ZSTD) {
        static_cast<ZstdKmerTmpFileWriter *>(writer)->write(data, byteSize);
        } else {
        static_cast<RawKmerTmpFileWriter *>(writer)->write(data, byteSize);
    }
}

static void bucketWriterClose(void *writer, int codec) {
    if (writer == NULL) {
        return;
    }
    if (codec == KMER_TMP_CODEC_ZSTD) {
        ZstdKmerTmpFileWriter *z = static_cast<ZstdKmerTmpFileWriter *>(writer);
        z->close();
        delete z;
        } else {
        RawKmerTmpFileWriter *r = static_cast<RawKmerTmpFileWriter *>(writer);
        r->close();
        delete r;
    }
}

// reads one bucket file into dst; overflow means it held more than the split array, which should not happen
static size_t bucketReadFile(const std::string &fileName, int codec, void *dst, size_t maxBytes,
                             size_t ioBufferBytes, bool &overflow) {
    overflow = false;
    FILE *f = fopen(fileName.c_str(), "rb");
    // callers verified the file exists, so failing here is a real IO error, not an empty bucket
    if (f == NULL) {
        kmerTmpIoErrorDie("open", fileName);
    }
    size_t written = 0;
    if (codec == KMER_TMP_CODEC_NONE) {
        written = fread(dst, 1, maxBytes, f);
        char probe;
        if (fread(&probe, 1, 1, f) == 1) {
            overflow = true;
        }
        fclose(f);
        return written;
    }

    ZSTD_DStream *dstream = NULL;
    size_t inBufferSize = std::max(ioBufferBytes, KMER_TMP_ZSTD_INPUT_BUFFER_SIZE);
    if (codec == KMER_TMP_CODEC_ZSTD) {
        dstream = ZSTD_createDStream();
        if (dstream == NULL) {
            Debug(Debug::ERROR) << "ZSTD_createDStream() failed for " << fileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        ZSTD_initDStream(dstream);
        inBufferSize = std::max(ioBufferBytes, ZSTD_DStreamInSize());
    }
    // decompression already lands straight in dst, so the read size is all that is left to set
    char *inBuffer = static_cast<char *>(malloc(inBufferSize));
    Util::checkAllocation(inBuffer, "Can not allocate k-mer bucket input buffer");
    const int fd = fileno(f);
    size_t readOffset = 0;
    size_t inSize = 0;
    size_t inPos = 0;
    while (overflow == false) {
        if (inPos == inSize) {
            const ssize_t got = kmerPreadRetry(fd, inBuffer, inBufferSize, readOffset);
            if (got < 0) {
                kmerTmpIoErrorDie("read", fileName);
            }
            readOffset += static_cast<size_t>(got);
            inSize = static_cast<size_t>(got);
            inPos = 0;
            if (inSize == 0) {
                break;  // EOF
            }
        }
        if (codec == KMER_TMP_CODEC_ZSTD) {
            ZSTD_inBuffer input = { inBuffer, inSize, inPos };
            ZSTD_outBuffer output = { dst, maxBytes, written };
            const size_t ret = ZSTD_decompressStream(dstream, &output, &input);
            if (ZSTD_isError(ret)) {
                Debug(Debug::ERROR) << "ZSTD_decompressStream() error for " << fileName << ": "
                                    << ZSTD_getErrorName(ret) << "\n";
                EXIT(EXIT_FAILURE);
            }
            written = output.pos;
            inPos = input.pos;
        }
        if (written == maxBytes && inPos < inSize) {
            overflow = true;  // output full but input remains
        }
    }
    free(inBuffer);
    if (dstream != NULL) {
        ZSTD_freeDStream(dstream);
    }
    fclose(f);
    return written;
}

template <typename T>
class KmerTmpFileReader {
public:
    KmerTmpFileReader(const std::string &fileName, size_t ioBufferBytes)
        : fileName(fileName), file(NULL), codec(codecFromFileName(fileName)),
          entries(NULL), entrySize(0), offsetPos(0), released(0), dataSize(0),
          dstream(NULL), inBuffer(NULL), outBuffer(NULL), inBufferSize(0), outBufferSize(0),
          input(), readOffset(0), outPos(0), outLen(0),
          eof(false), frameComplete(false), closed(true), ioBufferBytes(ioBufferBytes) {
        open();
    }

    ~KmerTmpFileReader() {
        close();
    }

    bool next(T &entry) {
        if (codec != KMER_TMP_CODEC_NONE) {
            return nextCompressed(entry);
        }
        if (offsetPos >= entrySize) {
            return false;
        }
        entry = entries[offsetPos];
        offsetPos++;
        releaseConsumed();
        return true;
    }

    void close() {
        if (closed) {
            return;
        }
        if (dataSize > 0 && entries != NULL && munmap((void *) entries, dataSize) < 0) {
            Debug(Debug::ERROR) << "Failed to munmap memory dataSize=" << dataSize << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (dstream != NULL) {
            ZSTD_freeDStream(dstream);
            dstream = NULL;
        }
        if (inBuffer != NULL) {
            free(inBuffer);
            inBuffer = NULL;
        }
        if (outBuffer != NULL) {
            free(outBuffer);
            outBuffer = NULL;
        }
        if (file != NULL && fclose(file) != 0) {
            kmerTmpIoErrorDie("close", fileName);
        }
        file = NULL;
        closed = true;
    }

private:
    // walked strictly forward, so release behind offsetPos in 64 MB steps, madvise before fadvise
    void releaseConsumed() {
        static const size_t RELEASE_STEP = 64ull * 1024 * 1024;
        const size_t consumed = offsetPos * sizeof(T);
        if (entries == NULL || consumed < released + RELEASE_STEP) {
            return;
        }
        const size_t page = Util::getPageSize();
        const size_t end = consumed & ~(page - 1);
        if (end <= released) {
            return;
        }
        char *base = reinterpret_cast<char *>(entries);
        ::madvise(base + released, end - released, MADV_DONTNEED);
#if defined(HAVE_POSIX_FADVISE)
        posix_fadvise(fileno(file), static_cast<off_t>(released),
                      static_cast<off_t>(end - released), POSIX_FADV_DONTNEED);
#endif
        released = end;
    }

    // the codec is not stored in the file, so the suffix the writer chose is what identifies it
    static int codecFromFileName(const std::string &fileName) {
        if (Util::endsWith(".zst", fileName)) {
            return KMER_TMP_CODEC_ZSTD;
        }
        return KMER_TMP_CODEC_NONE;
    }

    void open() {
        file = FileUtil::openFileOrDie(fileName.c_str(), "rb", true);
        if (codec != KMER_TMP_CODEC_NONE) {
            // both codecs want at least their own stream sizes; beyond that the budget sets the read size
            inBufferSize = std::max(ioBufferBytes, KMER_TMP_ZSTD_INPUT_BUFFER_SIZE);
            outBufferSize = std::max(ioBufferBytes, KMER_TMP_ZSTD_OUTPUT_BUFFER_SIZE);
            if (codec == KMER_TMP_CODEC_ZSTD) {
                dstream = ZSTD_createDStream();
                if (dstream == NULL) {
                    Debug(Debug::ERROR) << "ZSTD_createDStream() failed for " << fileName << "\n";
                    EXIT(EXIT_FAILURE);
                }
                size_t ret = ZSTD_initDStream(dstream);
                if (ZSTD_isError(ret)) {
                    Debug(Debug::ERROR) << "ZSTD_initDStream() error for " << fileName << ". Error "
                                        << ZSTD_getErrorName(ret) << "\n";
                    EXIT(EXIT_FAILURE);
                }
                inBufferSize = std::max(ioBufferBytes, ZSTD_DStreamInSize());
                outBufferSize = std::max(ioBufferBytes, ZSTD_DStreamOutSize());
            }
            inBuffer = static_cast<char *>(malloc(inBufferSize));
            outBuffer = static_cast<char *>(malloc(outBufferSize));
            if (inBuffer == NULL || outBuffer == NULL) {
                Debug(Debug::ERROR) << "Cannot allocate decompression buffer for " << fileName << "\n";
                EXIT(EXIT_FAILURE);
            }
            input.src = inBuffer;
            input.size = 0;
            input.pos = 0;
            readOffset = 0;
            outPos = 0;
            outLen = 0;
            eof = false;
            frameComplete = false;
        } else {
            struct stat sb;
            if (fstat(fileno(file), &sb) == 0 && sb.st_size > 0) {
                entries = static_cast<T *>(FileUtil::mmapFile(file, &dataSize));
                if (dataSize % sizeof(T) != 0) {
                    Debug(Debug::ERROR) << "Malformed kmer temporary file " << fileName
                                        << ": size is not a multiple of entry size\n";
                    EXIT(EXIT_FAILURE);
                }
                Util::madviseLogged(entries, dataSize, POSIX_MADV_SEQUENTIAL, fileName.c_str());
            } else {
                entries = NULL;
                dataSize = 0;
            }
            entrySize = dataSize / sizeof(T);
            offsetPos = 0;
            released = 0;
        }
        closed = false;
    }

    // entries are served from the decompression output buffer; only one straddling two blocks is moved
    bool nextCompressed(T &entry) {
        while (outLen - outPos < sizeof(T)) {
            const size_t carry = outLen - outPos;
            if (carry > 0) {
                memmove(outBuffer, outBuffer + outPos, carry);
            }
            outPos = 0;
            outLen = carry;
            const size_t produced = decompressChunk(outBuffer + outLen, outBufferSize - outLen);
            if (produced == 0) {
                if (outLen != 0) {
                    Debug(Debug::ERROR) << "Malformed compressed kmer temporary file " << fileName
                                        << ": trailing partial entry\n";
                    EXIT(EXIT_FAILURE);
                }
                return false;
            }
            outLen += produced;
        }

        memcpy(&entry, outBuffer + outPos, sizeof(T));
        outPos += sizeof(T);
        return true;
    }

    // returns 0 only at end of stream; a decode call that yields nothing just pulls more input
    size_t decompressChunk(char *dst, size_t cap) {
        size_t produced = 0;
        while (produced == 0) {
            if (input.pos == input.size && eof == false) {
                const ssize_t got = kmerPreadRetry(fileno(file), inBuffer, inBufferSize, readOffset);
                if (got < 0) {
                    kmerTmpIoErrorDie("read", fileName);
                }
                if (got == 0) {
                    eof = true;
                }
                readOffset += static_cast<size_t>(got);
                input.src = inBuffer;
                input.size = static_cast<size_t>(got);
                input.pos = 0;
            }
            if (input.pos == input.size && eof) {
                if (frameComplete == false) {
                    Debug(Debug::ERROR) << "Malformed compressed kmer temporary file " << fileName
                                        << ": truncated frame\n";
                    EXIT(EXIT_FAILURE);
                }
                return 0;
            }

            if (codec == KMER_TMP_CODEC_ZSTD) {
                ZSTD_outBuffer output = { dst, cap, produced };
                const size_t ret = ZSTD_decompressStream(dstream, &output, &input);
                if (ZSTD_isError(ret)) {
                    Debug(Debug::ERROR) << "ZSTD_decompressStream() error for " << fileName << ". Error "
                                        << ZSTD_getErrorName(ret) << "\n";
                    EXIT(EXIT_FAILURE);
                }
                produced = output.pos;
                frameComplete = (ret == 0);
            }
        }
        return produced;
    }

    std::string fileName;
    FILE *file;
    int codec;
    T *entries;
    size_t entrySize;
    size_t offsetPos;
    size_t released;
    size_t dataSize;
    ZSTD_DStream *dstream;
    char *inBuffer;
    char *outBuffer;
    size_t inBufferSize;
    size_t outBufferSize;
    ZSTD_inBuffer input;
    size_t readOffset;
    size_t outPos;
    size_t outLen;
    bool eof;
    bool frameComplete;
    bool closed;
    size_t ioBufferBytes;
};

template <int TYPE, typename T>
static bool queueNextEntryStreaming(KmerPositionQueue &queue, int file,
                                    KmerTmpFileReader<T> &reader, DBKeyType &repSeqId) {
    T entry;
    while (reader.next(entry)) {
        if (entry.seqId == DB_KEY_INVALID) {
            repSeqId = DB_KEY_INVALID;
            continue;
        }
        if (repSeqId == DB_KEY_INVALID) {
            repSeqId = entry.seqId;
            continue;
        }
        if(TYPE == Parameters::DBTYPE_NUCLEOTIDES){
            queue.push(FileKmerPosition(repSeqId, entry.seqId, entry.diagonal,
                                        entry.score, entry.getRev(), file));
        }else{
            queue.push(FileKmerPosition(repSeqId, entry.seqId, entry.diagonal,
                                        entry.score, file));
        }
        return true;
    }
    return false;
}

template <int TYPE, typename T>
void mergeKmerFilesAndOutput(DBWriter &dbw,
                             std::vector<std::string> tmpFiles,
                             std::vector<char> &repSequence,
                             int numThreads, int maxIter,
                             DBReader<DBKeyType> *sequenceDbr, bool localIds) {
    if (localIds && sequenceDbr == NULL) {
        Debug(Debug::ERROR) << "Local-ID kmermatcher output needs its source sequence DB.\n";
        EXIT(EXIT_FAILURE);
    }
    Debug(Debug::INFO) << "Merge splits ... ";

    std::vector<std::vector<std::string>> threadedFiles;
    threadedFiles.resize(numThreads);

    const int tmpFileCodec = Parameters::getInstance().compressKmerTmpFiles;
    for (int threadIdx = 0; threadIdx < numThreads; threadIdx++) {
        for (int iter = 0; iter < maxIter; iter++) {
            for (size_t i = 0; i < tmpFiles.size(); ++i) {
                std::string splitFileName = existingKmerTmpFileName(
                    kmerTmpFileName(tmpFiles[i], iter, threadIdx, KMER_TMP_CODEC_NONE), tmpFileCodec);
                if (splitFileName.empty() == false) {
                    threadedFiles[threadIdx].push_back(splitFileName);
                }
            }
        }
    }

    size_t maxFilesPerThread = 0;
    for (int threadIdx = 0; threadIdx < numThreads; threadIdx++) {
        maxFilesPerThread = std::max(maxFilesPerThread, threadedFiles[threadIdx].size());
    }
    // every lane holds one in and one out buffer for each file it merges
    const size_t readerBufferBytes = kmerSpoolBufferBytes(
        static_cast<size_t>(numThreads) * std::max<size_t>(maxFilesPerThread, 1) * 2);

    int mergeThreads = numThreads;
    if (maxFilesPerThread > 0) {
        // both codecs hold one descriptor per file per lane, so the fd throttle applies to both
        long openMax = sysconf(_SC_OPEN_MAX);
        size_t fdReserve = 2 * static_cast<size_t>(numThreads) + 64;
        size_t fdBudget = (openMax > 0 && static_cast<size_t>(openMax) > fdReserve)
                              ? (static_cast<size_t>(openMax) - fdReserve) : 1;
        if (maxFilesPerThread >= fdBudget) {
            Debug(Debug::ERROR) << "One merge lane needs " << maxFilesPerThread
                                << " open k-mer temporary files but the open-file limit is "
                                << openMax << ". Reduce kmermatcher splits or raise the limit (ulimit -n).\n";
            EXIT(EXIT_FAILURE);
        }
        mergeThreads = std::min(mergeThreads,
                                std::max(1, static_cast<int>(fdBudget / maxFilesPerThread)));

        if (tmpFileCodec != KMER_TMP_CODEC_NONE) {
            const size_t perReaderBytes = (2 * readerBufferBytes) + 4096;
            const size_t memoryLimit = Util::computeMemory(Parameters::getInstance().splitMemoryLimit);
            const size_t mergeHeadroom = std::max(KMER_MERGE_MIN_MEMORY_HEADROOM, memoryLimit / 10);
            const size_t memBudget = (memoryLimit > mergeHeadroom) ? (memoryLimit - mergeHeadroom) : memoryLimit / 2;
            const size_t bytesPerLane = maxFilesPerThread * perReaderBytes;
            if (bytesPerLane > 0) {
                mergeThreads = std::min(mergeThreads,
                                        std::max(1, static_cast<int>(memBudget / bytesPerLane)));
            }
        }
    }

    if (mergeThreads < numThreads) {
        Debug(Debug::INFO) << "K-mer tmp merge uses " << mergeThreads << " concurrent lanes for "
                           << numThreads << " writer lanes to stay within file descriptor/memory bounds.\n";
    }

    Timer mergeTimer;
#pragma omp parallel for num_threads(mergeThreads)
    for (int threadIdx = 0; threadIdx < numThreads; threadIdx++) {
        const int fileCnt = threadedFiles[threadIdx].size();
        if (fileCnt == 0) {
            continue;
        }

        KmerTmpFileReader<T> **readers = new KmerTmpFileReader<T> *[fileCnt];
        DBKeyType *repSeqIds = new DBKeyType[fileCnt];

        for (size_t file = 0; file < threadedFiles[threadIdx].size(); file++) {
            readers[file] = new KmerTmpFileReader<T>(threadedFiles[threadIdx][file], readerBufferBytes);
            repSeqIds[file]  = DB_KEY_INVALID;
        }

        KmerPositionQueue queue;
        for (int file = 0; file < fileCnt; file++) {
            queueNextEntryStreaming<TYPE, T>(queue, file, *readers[file], repSeqIds[file]);
        }

        std::string prefResultsOutString;
        prefResultsOutString.reserve(KMER_MERGE_RESULT_BUFFER_RESERVE);
        char buffer[1024];
        bool hasRepSeq = (repSequence.size() > 0);
        DBKeyType currRepSeq = DB_KEY_INVALID;
        DBKeyType prevHitId = DB_KEY_INVALID;
        short prevDiagonal = 0;
        int diagonalScore = 0;
        int bestDiagonalCnt = 0;
        int bestRevertMask = 0;
        short bestDiagonal = 0;
        int topScore = 0;
        bool hasHit = false;

        const auto flushHit = [&]() {
            if(hasHit == false){
                return;
            }
            if(diagonalScore >= bestDiagonalCnt){
                bestDiagonalCnt = diagonalScore;
                bestDiagonal = prevDiagonal;
            }
            hit_t h;
            h.seqId = prevHitId;
            h.prefScore = (bestRevertMask) ? -topScore : topScore;
            h.diagonal = bestDiagonal;
            int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
            prefResultsOutString.append(buffer, len);

            prevHitId = DB_KEY_INVALID;
            prevDiagonal = 0;
            diagonalScore = 0;
            bestDiagonalCnt = 0;
            bestRevertMask = 0;
            bestDiagonal = 0;
            topScore = 0;
            hasHit = false;
        };

        const auto flushRepSeq = [&]() {
            if(currRepSeq == DB_KEY_INVALID){
                return;
            }
            flushHit();
            const DBKeyType outputKey = localIds ? sequenceDbr->getDbKey(currRepSeq) : currRepSeq;
            dbw.writeData(prefResultsOutString.c_str(), prefResultsOutString.length(), outputKey, threadIdx);
            if (hasRepSeq) {
                repSequence[currRepSeq] = true;
            }
            prefResultsOutString.clear();
            if (prefResultsOutString.capacity() > KMER_MERGE_RESULT_BUFFER_MAX_KEEP) {
                std::string tmp;
                tmp.reserve(KMER_MERGE_RESULT_BUFFER_RESERVE);
                prefResultsOutString.swap(tmp);
            }
        };

        const auto startRepSeq = [&](DBKeyType repSeq) {
            flushRepSeq();
            currRepSeq = repSeq;
            if (hasRepSeq) {
                hit_t h;
                h.seqId = currRepSeq;
                h.prefScore = 0;
                h.diagonal = 0;
                int len = QueryMatcher::prefilterHitToBuffer(buffer, h);
                prefResultsOutString.append(buffer, len);
            }
        };

        while (queue.empty() == false) {
            FileKmerPosition res = queue.top();
            queue.pop();

            if (currRepSeq != res.repSeq) {
                startRepSeq(res.repSeq);
            }

            bool hitIsRepSeq = (currRepSeq == res.id);
            if (hitIsRepSeq) {
                queueNextEntryStreaming<TYPE, T>(queue, res.file, *readers[res.file], repSeqIds[res.file]);
                continue;
            }

            if (hasHit == false || prevHitId != res.id) {
                flushHit();
                hasHit = true;
                prevHitId = res.id;
                prevDiagonal = res.pos;
                diagonalScore = 0;
            } else if (prevDiagonal != res.pos) {
                if(diagonalScore >= bestDiagonalCnt){
                    bestDiagonalCnt = diagonalScore;
                    bestDiagonal = prevDiagonal;
                }
                prevDiagonal = res.pos;
                diagonalScore = 0;
            }

            diagonalScore += res.score;
            if(diagonalScore >= bestDiagonalCnt){
                bestDiagonalCnt = diagonalScore;
                bestDiagonal = res.pos;
                bestRevertMask = res.reverse;
            }
            topScore += res.score;

            queueNextEntryStreaming<TYPE, T>(queue, res.file, *readers[res.file], repSeqIds[res.file]);
        }

        if(currRepSeq != DB_KEY_INVALID){
            flushRepSeq();
        }

        for (size_t file = 0; file < threadedFiles[threadIdx].size(); file++) {
            readers[file]->close();
            delete readers[file];
            // a lane owns its files, so unlinking here frees disk now, not after the slowest lane
            FileUtil::remove(threadedFiles[threadIdx][file].c_str());
        }

        delete[] repSeqIds;
        delete[] readers;
    }

    Debug(Debug::INFO) << "Time for merging k-mer splits: " << mergeTimer.lap() << "\n";

    for (int tid = 0; tid < numThreads; ++tid) {
        for (std::string file : threadedFiles[tid]) {
            // already unlinked by its lane unless that lane died before reaching the loop above
            if (FileUtil::fileExists(file.c_str())) {
                FileUtil::remove(file.c_str());
            }
        }
    }
}

template <int TYPE, typename T, typename seqLenType, bool includeAdjacency, bool IncludeSeqLen>
void writeKmersToDisk(std::string tmpFile, KmerPosition<seqLenType, includeAdjacency, IncludeSeqLen> *hashSeqPair, size_t totalKmers,
                      int numThreads, std::vector<size_t> *threadQueryOffsets, int iteration) {
    const size_t BUFFER_SIZE = 2048;
    const Parameters &par = Parameters::getInstance();
    // the result writer stays on the direct path; bucket files keep their own codec policy
    const bool compressTmpFiles = par.compressKmerTmpFiles == KMER_TMP_CODEC_ZSTD;
#ifndef OPENMP
    (void) numThreads;
#endif

#pragma omp parallel num_threads(numThreads)
    {
        int tid = 0;
#ifdef OPENMP
        tid = omp_get_thread_num();
#endif

        size_t startIdx, endIdx;

        if (threadQueryOffsets == nullptr) {
            startIdx = 0;
            endIdx = totalKmers;
        } else {
            startIdx = (*threadQueryOffsets)[tid];
            endIdx = (*threadQueryOffsets)[tid + 1];
        }
        
        const int tmpFileCodec = compressTmpFiles ? KMER_TMP_CODEC_ZSTD : KMER_TMP_CODEC_NONE;
        std::string tmpFileThread = kmerTmpFileName(tmpFile, iteration, tid, tmpFileCodec);
        // a rerun under a different codec would otherwise leave the old file for the merge to pick up
        removeKmerTmpFileIfExists(kmerTmpFileName(tmpFile, iteration, tid, KMER_TMP_CODEC_NONE));
        removeKmerTmpFileIfExists(kmerTmpFileName(tmpFile, iteration, tid, KMER_TMP_CODEC_ZSTD));
        if (startIdx < endIdx && hashSeqPair[startIdx].kmer != SIZE_T_MAX) {
            FILE *filePtr = NULL;
            ZstdKmerTmpFileWriter *zstdWriter = NULL;
            if (compressTmpFiles) {
                zstdWriter = new ZstdKmerTmpFileWriter(tmpFileThread, KMER_TMP_ZSTD_COMPRESSION_LEVEL);
            } else {
                filePtr = openKmerTmpFileForOverwriteOrDie(tmpFileThread, "wb");
            }

            const auto writeEntries = [&](const T *entries, size_t count) {
                if (count == 0) {
                    return;
                }
                if (compressTmpFiles) {
                    zstdWriter->write(entries, sizeof(T) * count);
                } else {
                    writeAllOrDie(filePtr, entries, sizeof(T) * count, tmpFileThread);
                }
            };

            size_t repSeqId = SIZE_T_MAX;
            size_t lastTargetId = SIZE_T_MAX;
            seqLenType lastDiagonal = 0;
            int diagonalScore = 0;
            // counts every hit line a lane writes, so a 32 bit counter can wrap back to zero at 10^10
            size_t writeSets = 0;
            size_t bufferPos = 0;
            size_t elementCnt = 0;

            T writeBuffer[BUFFER_SIZE];
            T nullEntry;
            nullEntry.seqId = DB_KEY_INVALID;
            nullEntry.diagonal = 0;

            for (size_t kmerPos = startIdx; kmerPos < endIdx && hashSeqPair[kmerPos].kmer != SIZE_T_MAX; kmerPos++) {
                size_t currKmer = hashSeqPair[kmerPos].kmer;
                if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                    currKmer = BIT_CLEAR(currKmer, 63);
                }

                if (repSeqId != currKmer) {
                    if (writeSets > 0 && elementCnt > 0) {
                        if (bufferPos > 0) {
                            writeEntries(writeBuffer, bufferPos);
                        }
                        writeEntries(&nullEntry, 1);
                    }
                    lastTargetId = SIZE_T_MAX;
                    bufferPos = 0;
                    elementCnt = 0;
                    repSeqId = currKmer;
                    writeBuffer[bufferPos].seqId = repSeqId;
                    writeBuffer[bufferPos].score = 0;
                    writeBuffer[bufferPos].diagonal = 0;
                    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        bool isReverse = (BIT_CHECK(hashSeqPair[kmerPos].kmer, 63) == false);
                        writeBuffer[bufferPos].setReverse(isReverse);
                    }
                    bufferPos++;
                }

                DBKeyType targetId = hashSeqPair[kmerPos].id;
                seqLenType diagonal = hashSeqPair[kmerPos].pos;
                int forward = 0;
                int reverse = 0;

                do {
                    diagonalScore += (diagonalScore == 0 || (lastTargetId == targetId && lastDiagonal == diagonal));
                    lastTargetId = hashSeqPair[kmerPos].id;
                    lastDiagonal = hashSeqPair[kmerPos].pos;
                    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        bool isReverse = (BIT_CHECK(hashSeqPair[kmerPos].kmer, 63) == false);
                        forward += (isReverse == false);
                        reverse += (isReverse == true);
                    }
                    kmerPos++;
                } while (kmerPos < endIdx &&
                        hashSeqPair[kmerPos].kmer != SIZE_T_MAX &&
                        repSeqId == hashSeqPair[kmerPos].kmer &&
                        targetId == hashSeqPair[kmerPos].id &&
                        hashSeqPair[kmerPos].pos == diagonal);
                kmerPos--;

                elementCnt++;
                // score is one byte, so a longer run is emitted as several records that the merge adds up
                int remainingScore = diagonalScore;
                diagonalScore = 0;
                do {
                    int chunk = (remainingScore > UCHAR_MAX) ? UCHAR_MAX : remainingScore;
                    writeBuffer[bufferPos].seqId = targetId;
                    writeBuffer[bufferPos].score = static_cast<unsigned char>(chunk);
                    writeBuffer[bufferPos].diagonal = diagonal;
                    if (TYPE == Parameters::DBTYPE_NUCLEOTIDES) {
                        bool isReverse = (reverse > forward) ? true : false;
                        writeBuffer[bufferPos].setReverse(isReverse);
                    }
                    bufferPos++;

                    if (bufferPos >= BUFFER_SIZE) {
                        writeEntries(writeBuffer, bufferPos);
                        bufferPos = 0;
                    }
                    remainingScore -= chunk;
                } while (remainingScore > 0);
                lastTargetId = targetId;
                writeSets++;
            }

            if (writeSets > 0 && elementCnt > 0) {
                if (bufferPos > 0) {
                    writeEntries(writeBuffer, bufferPos);
                }
                writeEntries(&nullEntry, 1);
            }
            if (compressTmpFiles) {
                zstdWriter->close();
                delete zstdWriter;
            } else if (fclose(filePtr) != 0) {
                kmerTmpIoErrorDie("close", tmpFileThread);
            }
        }
    }
    
}

void setKmerLengthAndAlphabet(Parameters &parameters, size_t aaDbSize, int seqTyp) {
    if(Parameters::isEqualDbtype(seqTyp, Parameters::DBTYPE_NUCLEOTIDES)){
        if(parameters.kmerSize == 0) {
            parameters.kmerSize = std::max(17, static_cast<int>(log(static_cast<float>(aaDbSize))/log(4)));
            parameters.spacedKmerPattern = "";
            parameters.alphabetSize = 5;
        }
        if(parameters.kmersPerSequence == 0){
            parameters.kmersPerSequence = 60;
        }
    }else{
        if(parameters.kmerSize == 0){
            if((parameters.seqIdThr+0.001)>=0.99){
                parameters.kmerSize = 14;
                parameters.alphabetSize = 21;
            }else if((parameters.seqIdThr+0.001)>=0.9){
                parameters.kmerSize = 14;
                parameters.alphabetSize = 13;
            }else{
                parameters.kmerSize = std::max(10, static_cast<int>(log(static_cast<float>(aaDbSize))/log(8.7)));
                parameters.alphabetSize = 13;
            }
            parameters.spacedKmerPattern = "";
        }
        if(parameters.kmersPerSequence == 0){
            parameters.kmersPerSequence = 20;
        }
    }
}

// only what other translation units need is explicit; kmersearch/kmerindexdb use the IncludeSeqLen=true variants
template std::pair<size_t, size_t>  fillKmerPositionArray<0, short, false, true>(KmerPosition<short, false, true> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *, KmerPartitionSink<short, false, true> *);
template std::pair<size_t, size_t>  fillKmerPositionArray<1, short, false, true>(KmerPosition<short, false, true> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *, KmerPartitionSink<short, false, true> *);
template std::pair<size_t, size_t>  fillKmerPositionArray<2, short, false, true>(KmerPosition<short, false, true> *, size_t, DBReader<DBKeyType> &, Parameters &, BaseMatrix *, bool, size_t, size_t, size_t *, KmerPartitionSink<short, false, true> *);

template KmerPosition<short, false, true> *initKmerPositionMemory(size_t size);

template size_t computeMemoryNeededLinearfilter<short, false, true>(size_t totalKmer);

template std::vector<std::pair<size_t, size_t>>  setupKmerSplits<short, false, true>(Parameters &, BaseMatrix *, DBReader<DBKeyType> &, size_t, size_t, size_t *);

template void writeKmersToDisk<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev, short, false, true>(std::string, KmerPosition<short, false, true> *, size_t, int, std::vector<size_t> *, int);
template void writeKmersToDisk<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry, short, false, true>(std::string, KmerPosition<short, false, true> *, size_t, int, std::vector<size_t> *, int);

template void mergeKmerFilesAndOutput<Parameters::DBTYPE_NUCLEOTIDES, KmerEntryRev>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int, DBReader<DBKeyType> *, bool);
template void mergeKmerFilesAndOutput<Parameters::DBTYPE_AMINO_ACIDS, KmerEntry>(DBWriter &, std::vector<std::string>, std::vector<char> &, int, int, DBReader<DBKeyType> *, bool);

#undef SIZE_T_MAX
