/*
 * createdbparallel
 *
 * Builds a sequence database from FASTA input with many independent worker
 * processes that coordinate only through a shared directory -- no MPI, no rank
 * argument, and no direct node-to-node communication. Every worker runs the same
 * command line; identity comes from a counter file.
 *
 * Two things make this different from stock createdb, and both exist to survive
 * a trillion sequences:
 *
 *   - Keys are *length-ranked and dense*: key 0 is the longest sequence and key
 *     == global order rank. That is what lets DenseIndex address an entry by key
 *     without a resident index, and what makes the distributed greedy in the
 *     reduce stage exact (order within a component is just sorted keys).
 *   - Nothing is merged. A first pass histograms sequence lengths per input byte
 *     range; those histograms alone determine every byte offset in the finished
 *     database (see LengthRankedPlan.h), so the second pass writes each sequence
 *     straight to its final position. Stock createdb's GPU mode achieves the same
 *     ordering, but only by sorting and then merging every shard through one
 *     process holding an entry per sequence in RAM.
 *
 * Memory per worker is bounded by --chunk-size times the thread count, not by
 * the size of the input, which is the property the whole design rests on.
 */
#include "Command.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "DenseIndex.h"
#include "FastSort.h"
#include "FileUtil.h"
#include "KSeqWrapper.h"
#include "KmerPartition.h"
#include "LengthRankedPlan.h"
#include "ParallelCoordination.h"
#include "Parameters.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {

// One work item: the records whose '>' falls in [begin, end) of one input file.
struct Chunk {
    size_t fileIdx;
    size_t begin;
    size_t end;
};

// Per-chunk coordination files live in a subdirectory of 1000, not all in one.
//
// At 1e12 sequences the input is ~350 TB, which is 1.37e6 chunks at the 256 MB
// default and therefore 2.74e6 entries -- a histogram and a plan each -- in a
// single directory. That is a directory no shared filesystem enjoys creating,
// listing or removing. Grouping them keeps any one directory at a thousand files
// while the path stays a pure function of the chunk index, which is what every
// worker relies on.
const size_t COORD_FILES_PER_DIR = 1000;

std::string chunkGroupDir(const std::string &coordDir, size_t chunkIdx) {
    return coordDir + "/c" + SSTR(chunkIdx / COORD_FILES_PER_DIR);
}

// Racing workers may both create it; only a failure that also leaves no directory
// behind is real.
void ensureChunkGroupDir(const std::string &coordDir, size_t chunkIdx) {
    const std::string path = chunkGroupDir(coordDir, chunkIdx);
    if (FileUtil::directoryExists(path.c_str()) == true) {
        return;
    }
    if (mkdir(path.c_str(), 0777) != 0 && FileUtil::directoryExists(path.c_str()) == false) {
        Debug(Debug::ERROR) << "Cannot create " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

std::string chunkHistPath(const std::string &coordDir, size_t chunkIdx) {
    return chunkGroupDir(coordDir, chunkIdx) + "/chunk." + SSTR(chunkIdx) + ".hist";
}

std::string chunkPlanPath(const std::string &coordDir, size_t chunkIdx) {
    return chunkGroupDir(coordDir, chunkIdx) + "/chunk." + SSTR(chunkIdx) + ".plan";
}

void writeAt(int fd, const void *data, size_t length, size_t offset, const char *what) {
    const char *cursor = static_cast<const char *>(data);
    size_t done = 0;
    while (done < length) {
        const ssize_t written = pwrite(fd, cursor + done, length - done, static_cast<off_t>(offset + done));
        if (written <= 0) {
            if (errno == EINTR) {
                continue;
            }
            Debug(Debug::ERROR) << "Cannot write " << what << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        done += static_cast<size_t>(written);
    }
}

int openForWrite(const std::string &path) {
    const int fd = open(path.c_str(), O_WRONLY);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << path << " for writing: " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    return fd;
}

// Offset of the first record that starts at or after `from`.
//
// A record starts at a '>' that is either the first byte of the file or directly
// preceded by a newline. Resolving boundaries this way is what makes chunks
// independent: a record is owned by exactly the chunk its '>' falls in, so no
// record is split between chunks or claimed by two of them.
size_t findRecordStart(int fd, size_t from, size_t fileSize) {
    if (from == 0) {
        return 0;
    }
    if (from >= fileSize) {
        return fileSize;
    }

    const size_t windowSize = 65536;
    std::vector<char> window(windowSize);
    // Start one byte early so a "\n>" straddling the requested offset is seen.
    size_t pos = from - 1;
    while (pos < fileSize) {
        const size_t want = std::min(windowSize, fileSize - pos);
        ssize_t got = pread(fd, window.data(), want, static_cast<off_t>(pos));
        if (got <= 0) {
            if (errno == EINTR) {
                continue;
            }
            Debug(Debug::ERROR) << "Cannot read input file: " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        for (ssize_t i = 0; i + 1 < got; i++) {
            if (window[i] == '\n' && window[i + 1] == '>') {
                return pos + static_cast<size_t>(i) + 1;
            }
        }
        if (static_cast<size_t>(got) < want) {
            continue;
        }
        // Overlap by one byte so a pair split across windows is not missed.
        pos += static_cast<size_t>(got) - 1;
    }
    return fileSize;
}

// How many chunks each input file is cut into, and where each file's chunks start
// in the global numbering. No I/O at all: a chunk is the piece of the fixed
// `chunkSize` grid it sits on, so the count follows from the file size.
//
// This replaces a planner that resolved *every* chunk boundary up front, in every
// worker, with a 64 KB pread each. At 1e12 (350 TB at the 256 MB default) that is
// 1.37e6 preads -- ~90 GB per worker and ~730 TB across 8192 of them -- paid
// before any worker starts work, every time the stage is entered.
//
// The numbering is still a pure function of (input, chunkSize) and still runs in
// file order, which is what the key assignment's tie-break depends on: ties go to
// the lower chunk index and then to position within the chunk, i.e. input order.
// The grid can leave a chunk empty where one record spans a whole boundary
// region; an empty chunk contributes nothing to the ordering, so the keys are
// unchanged.
struct ChunkLayout {
    std::vector<size_t> countPerFile;
    std::vector<size_t> firstChunkOfFile;
    size_t total;

    ChunkLayout() : total(0) {}
};

ChunkLayout planChunkLayout(const std::vector<std::string> &filenames, size_t chunkSize) {
    ChunkLayout layout;
    layout.countPerFile.assign(filenames.size(), 0);
    layout.firstChunkOfFile.assign(filenames.size(), 0);
    for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
        const size_t fileSize = FileUtil::getFileSize(filenames[fileIdx]);
        layout.firstChunkOfFile[fileIdx] = layout.total;
        layout.countPerFile[fileIdx] = (fileSize + chunkSize - 1) / chunkSize;
        layout.total += layout.countPerFile[fileIdx];
    }
    return layout;
}

size_t fileOfChunk(const ChunkLayout &layout, size_t chunkIdx) {
    // Few input files in practice, and this runs once per work item.
    for (size_t fileIdx = layout.countPerFile.size(); fileIdx > 0; fileIdx--) {
        if (chunkIdx >= layout.firstChunkOfFile[fileIdx - 1]
                && layout.countPerFile[fileIdx - 1] > 0) {
            return fileIdx - 1;
        }
    }
    return 0;
}

// Resolves one chunk's byte range: two preads, in the worker that owns it.
Chunk resolveChunk(const std::vector<std::string> &filenames, const ChunkLayout &layout,
                   size_t chunkSize, size_t chunkIdx) {
    const size_t fileIdx = fileOfChunk(layout, chunkIdx);
    const size_t local = chunkIdx - layout.firstChunkOfFile[fileIdx];
    const size_t fileSize = FileUtil::getFileSize(filenames[fileIdx]);
    const int fd = open(filenames[fileIdx].c_str(), O_RDONLY);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << filenames[fileIdx] << ": " << strerror(errno)
                            << "\n";
        EXIT(EXIT_FAILURE);
    }
    Chunk chunk;
    chunk.fileIdx = fileIdx;
    chunk.begin = findRecordStart(fd, local * chunkSize, fileSize);
    chunk.end = findRecordStart(fd, std::min((local + 1) * chunkSize, fileSize), fileSize);
    close(fd);
    if (chunk.end < chunk.begin) {
        chunk.end = chunk.begin;
    }
    return chunk;
}

// Keeps the chunk count within reach of the coordination machinery.
//
// The default 256 MB is right for the sizes this is usually run at, but at 1e12
// (350 TB) it derives 1.37e6 chunks -- two work queues of that size, 2.74e6
// coordination files, and a planner sweep that holds every histogram and plan
// resident (~550 GB) while it runs single-threaded on one node. Scaling the chunk
// size so the count lands near TARGET_CHUNKS makes most of that go away for free
// and costs only a larger per-item buffer, which is bounded by --chunk-size times
// the thread count either way.
//
// Only ever *raises* it: a user who asked for a specific --chunk-size gets it, and
// small inputs keep small chunks so a second worker still has something to claim.
const size_t TARGET_CHUNKS = 100000;

size_t deriveChunkSize(const std::vector<std::string> &filenames, size_t requested,
                       bool userSupplied) {
    if (userSupplied) {
        return requested;
    }
    uint64_t totalBytes = 0;
    for (size_t i = 0; i < filenames.size(); i++) {
        totalBytes += FileUtil::getFileSize(filenames[i]);
    }
    const uint64_t wanted = totalBytes / TARGET_CHUNKS;
    if (wanted <= requested) {
        return requested;
    }
    // Rounded up to a whole number of the requested unit, so the derived value
    // stays a recognisable multiple of the documented default.
    return static_cast<size_t>((wanted + requested - 1) / requested) * requested;
}

// Reads a chunk's bytes so they can be handed to the FASTA parser in one piece.
void readChunk(const std::string &filename, const Chunk &chunk, std::vector<char> &buffer) {
    const int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << filename << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t length = chunk.end - chunk.begin;
    buffer.resize(length);
    size_t done = 0;
    while (done < length) {
        const ssize_t got = pread(fd, buffer.data() + done, length - done,
                                  static_cast<off_t>(chunk.begin + done));
        if (got <= 0) {
            if (errno == EINTR) {
                continue;
            }
            Debug(Debug::ERROR) << "Cannot read " << filename << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        done += static_cast<size_t>(got);
    }
    close(fd);
}

// Builds the header exactly as createdb does, so a database built either way
// holds byte-identical header entries.
void buildHeader(const KSeqWrapper::KSeqEntry &entry, std::string &header) {
    header.clear();
    header.append(entry.name.s, entry.name.l);
    if (entry.comment.l > 0) {
        header.append(" ", 1);
        header.append(entry.comment.s, entry.comment.l);
    }
    header.push_back('\n');
}

// Pass 1: report what the chunk contains without writing any sequence data.
ChunkHistogram scanChunk(const std::string &filename, const Chunk &chunk, size_t chunkIdx) {
    std::vector<char> buffer;
    readChunk(filename, chunk, buffer);

    // Counted into parallel arrays kept ascending by length. The number of
    // distinct lengths is small (a few thousand for proteins) regardless of how
    // big the chunk is, which is what keeps a histogram tens of KB rather than
    // proportional to the sequence count.
    std::vector<uint64_t> lengthOf;
    std::vector<uint64_t> countOf;
    std::vector<uint64_t> headerBytesOf;
    std::vector<uint64_t> accessionBytesOf;

    KSeqBuffer kseq(buffer.data(), buffer.size());
    std::string header;
    uint64_t seqCount = 0;
    uint64_t nuclVotes = 0;
    uint64_t sampleCount = 0;
    const size_t testForNucSequence = 100;

    while (kseq.ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq.entry;
        if (entry.name.l == 0) {
            Debug(Debug::ERROR) << "Invalid FASTA entry in " << filename << " at byte "
                                << chunk.begin << "\n";
            EXIT(EXIT_FAILURE);
        }
        buildHeader(entry, header);

        const uint64_t length = entry.sequence.l;
        // Insertion sort into the ascending length list. Chunks hold few distinct
        // lengths, so a linear probe from the end is cheaper than a hash map.
        size_t slot = lengthOf.size();
        while (slot > 0 && lengthOf[slot - 1] > length) {
            slot--;
        }
        if (slot == 0 || lengthOf[slot - 1] != length) {
            lengthOf.insert(lengthOf.begin() + slot, length);
            countOf.insert(countOf.begin() + slot, 0);
            headerBytesOf.insert(headerBytesOf.begin() + slot, 0);
            accessionBytesOf.insert(accessionBytesOf.begin() + slot, 0);
        } else {
            slot--;
        }
        countOf[slot]++;
        // DBWriter appends a NUL after the header text.
        headerBytesOf[slot] += header.length() + 1;
        // The same extraction the emit pass will do, so the planner's byte count
        // and the bytes actually written cannot disagree.
        accessionBytesOf[slot] += Util::parseFastaHeader(header.c_str()).length();

        if (sampleCount < testForNucSequence) {
            size_t nucleotideLike = 0;
            for (size_t i = 0; i < entry.sequence.l; i++) {
                switch (toupper(entry.sequence.s[i])) {
                    case 'T': case 'A': case 'G': case 'C': case 'U': case 'N':
                        nucleotideLike++;
                        break;
                }
            }
            if (entry.sequence.l > 0 &&
                static_cast<float>(nucleotideLike) / static_cast<float>(entry.sequence.l) > 0.9f) {
                nuclVotes++;
            }
            sampleCount++;
        }
        seqCount++;
    }

    ChunkHistogram histogram;
    histogram.chunkIdx = chunkIdx;
    histogram.fileIdx = chunk.fileIdx;
    histogram.seqCount = seqCount;
    histogram.nuclVotes = nuclVotes;
    histogram.sampleCount = sampleCount;
    histogram.buckets.resize(lengthOf.size());
    for (size_t i = 0; i < lengthOf.size(); i++) {
        histogram.buckets[i].length = lengthOf[i];
        histogram.buckets[i].count = countOf[i];
        histogram.buckets[i].headerBytes = headerBytesOf[i];
        histogram.buckets[i].accessionBytes = accessionBytesOf[i];
    }
    return histogram;
}

// Everything one length bucket of one chunk writes. Buffered rather than written
// per sequence so each bucket becomes a single large, contiguous write instead of
// a scatter of ~250 byte writes over a parallel filesystem.
struct BucketOutput {
    uint64_t keyStart;
    uint64_t lookupOffset;
    std::string lookupText;
    uint64_t dataOffset;
    uint64_t hdrOffset;
    uint64_t written;
    std::vector<char> seqData;
    std::vector<char> hdrData;
    std::vector<DenseIndex::Entry> seqIndex;
    std::vector<DenseIndex::Entry> hdrIndex;
};

// Pass 2: rescan the chunk and write every sequence at its planned position.
void emitChunk(const std::string &filename, const Chunk &chunk, const ChunkPlan &plan,
               int seqFd, int hdrFd, int seqIdxFd, int hdrIdxFd, int lookupFd) {
    std::vector<char> buffer;
    readChunk(filename, chunk, buffer);

    std::vector<BucketOutput> buckets(plan.entries.size());
    for (size_t i = 0; i < plan.entries.size(); i++) {
        buckets[i].keyStart = plan.entries[i].keyStart;
        buckets[i].dataOffset = plan.entries[i].dataOffset;
        buckets[i].hdrOffset = plan.entries[i].hdrOffset;
        buckets[i].lookupOffset = plan.entries[i].lookupOffset;
        buckets[i].written = 0;
    }
    const std::string fileIdxField = "\t" + SSTR(plan.fileIdx) + "\n";

    KSeqBuffer kseq(buffer.data(), buffer.size());
    std::string header;
    while (kseq.ReadEntry()) {
        const KSeqWrapper::KSeqEntry &entry = kseq.entry;
        buildHeader(entry, header);
        const uint64_t length = entry.sequence.l;

        // The plan lists one entry per length, ascending, so this is a binary
        // search over a few thousand entries at most.
        size_t lo = 0;
        size_t hi = plan.entries.size();
        while (lo < hi) {
            const size_t mid = lo + (hi - lo) / 2;
            if (plan.entries[mid].length < length) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if (lo >= plan.entries.size() || plan.entries[lo].length != length) {
            Debug(Debug::ERROR) << "Chunk " << plan.chunkIdx << " holds a sequence of length "
                                << length << " that pass 1 did not record. The input changed "
                                << "between the two passes.\n";
            EXIT(EXIT_FAILURE);
        }
        BucketOutput &bucket = buckets[lo];

        DenseIndex::Entry seqEntry;
        seqEntry.offset = bucket.dataOffset + bucket.seqData.size();
        seqEntry.length = static_cast<uint32_t>(length + 2);
        bucket.seqIndex.push_back(seqEntry);
        bucket.seqData.insert(bucket.seqData.end(), entry.sequence.s, entry.sequence.s + length);
        bucket.seqData.push_back('\n');
        bucket.seqData.push_back('\0');

        DenseIndex::Entry hdrEntry;
        hdrEntry.offset = bucket.hdrOffset + bucket.hdrData.size();
        hdrEntry.length = static_cast<uint32_t>(header.length() + 1);
        bucket.hdrIndex.push_back(hdrEntry);
        bucket.hdrData.insert(bucket.hdrData.end(), header.begin(), header.end());
        bucket.hdrData.push_back('\0');

        // The bucket's keys run from keyStart in the order it consumes them, so
        // the key of this line is known without any shared counter.
        bucket.lookupText.append(SSTR(bucket.keyStart + bucket.written));
        bucket.lookupText.push_back('\t');
        bucket.lookupText.append(Util::parseFastaHeader(header.c_str()));
        bucket.lookupText.append(fileIdxField);

        bucket.written++;
    }

    for (size_t i = 0; i < buckets.size(); i++) {
        BucketOutput &bucket = buckets[i];
        if (bucket.written != plan.entries[i].count) {
            Debug(Debug::ERROR) << "Chunk " << plan.chunkIdx << " wrote " << bucket.written
                                << " sequences of length " << plan.entries[i].length
                                << " but pass 1 counted " << plan.entries[i].count << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (bucket.lookupText.size() != plan.entries[i].lookupBytes) {
            Debug(Debug::ERROR) << "Chunk " << plan.chunkIdx << " built "
                                << bucket.lookupText.size() << " lookup bytes for length "
                                << plan.entries[i].length << " but the plan reserved "
                                << plan.entries[i].lookupBytes << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (bucket.written == 0) {
            continue;
        }
        writeAt(seqFd, bucket.seqData.data(), bucket.seqData.size(), bucket.dataOffset, "sequence data");
        writeAt(hdrFd, bucket.hdrData.data(), bucket.hdrData.size(), bucket.hdrOffset, "header data");
        writeAt(seqIdxFd, bucket.seqIndex.data(), bucket.seqIndex.size() * sizeof(DenseIndex::Entry),
                DenseIndex::entryOffset(bucket.keyStart), "sequence index");
        writeAt(hdrIdxFd, bucket.hdrIndex.data(), bucket.hdrIndex.size() * sizeof(DenseIndex::Entry),
                DenseIndex::entryOffset(bucket.keyStart), "header index");
        if (lookupFd >= 0) {
            writeAt(lookupFd, bucket.lookupText.data(), bucket.lookupText.size(),
                    bucket.lookupOffset, "lookup");
        }
    }
}

// Creates a file of exactly `size` bytes for the workers to pwrite into.
void allocateFile(const std::string &path, size_t size) {
    FILE *file = FileUtil::openAndDelete(path.c_str(), "wb");
    if (ftruncate(fileno(file), static_cast<off_t>(size)) != 0) {
        Debug(Debug::ERROR) << "Cannot size " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
}

// Runs a queue to completion with every thread of this worker claiming its own
// items. Leases mean a worker that dies mid-item only delays that item; another
// worker picks it up once the lease expires.
//
// Every thread drains the queue, so each needs an id of its own.
//
// claim()'s lease-expiry test and renew()'s ownership test are both keyed on
// (item, worker). With one id shared by all threads of a process, a thread whose
// lease lapsed could have its item re-claimed by a *sibling* thread, and the first
// thread's heartbeat would then keep renewing a record the sibling now owns --
// both running the same item with the queue providing no exclusion at all. Only
// the idempotence of the bodies here made that harmless, which is not a property
// to rely on silently.
//
// workerId * threads + threadNum is unique across the run because worker ids come
// from a shared counter and every process uses the same thread count.
template <typename Body>
void runQueue(WorkQueue &queue, int threads, int64_t workerId, Body body) {
    bool stalled = false;
#pragma omp parallel num_threads(threads)
    {
        int threadNum = 0;
#ifdef OPENMP
        threadNum = omp_get_thread_num();
#endif
        const int64_t threadWorkerId = workerId * static_cast<int64_t>(threads) + threadNum;
        // drain() rather than a plain claim loop: it goes back to claiming after
        // waiting, so items abandoned by a crashed worker are picked up once their
        // leases expire instead of being waited on by everyone forever.
        if (queue.drain(threadWorkerId, body) == false) {
#pragma omp critical
            stalled = true;
        }
    }
    if (stalled) {
        Debug(Debug::ERROR) << "Work queue stalled: work remains but no item is claimable\n";
        EXIT(EXIT_FAILURE);
    }
}

}  // namespace

int createdbparallel(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_VARIADIC, 0);
    par.printParameters(command.cmd, argc, argv, *command.params);

    std::vector<std::string> filenames(par.filenames);
    const std::string dataFile = filenames.back();
    filenames.pop_back();

    // Same ordering rule as createdb, so file numbering matches between them.
    SORT_SERIAL(filenames.begin(), filenames.end(), [](const std::string &a, const std::string &b) {
        return FileUtil::baseName(a) < FileUtil::baseName(b);
    });
    for (size_t i = 0; i < filenames.size(); i++) {
        if (FileUtil::directoryExists(filenames[i].c_str()) == true) {
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory\n";
            EXIT(EXIT_FAILURE);
        }
        // Named here rather than left to the planner. Compressed input is read as
        // plain text, finds no '>' and plans zero sequences, so the run used to die
        // several stages later with "The input files have no entry" -- true, and
        // useless. createdb decompresses; this cannot, because both passes address
        // the file by byte offset and a gzip stream has no seekable record
        // boundaries.
        unsigned char magic[2] = {0, 0};
        FILE *probe = fopen(filenames[i].c_str(), "rb");
        if (probe != NULL) {
            const size_t got = fread(magic, 1, sizeof(magic), probe);
            fclose(probe);
            if (got == sizeof(magic)
                    && ((magic[0] == 0x1f && magic[1] == 0x8b)
                        || (magic[0] == 'B' && magic[1] == 'Z'))) {
                Debug(Debug::ERROR) << filenames[i] << " is compressed, which this command cannot "
                                    << "read: both passes address the input by byte offset, and a "
                                    << "compressed stream has no seekable record boundaries. "
                                    << "Decompress it first.\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }

    const std::string hdrDataFile = dataFile + "_h";
    const std::string lookupFile = dataFile + ".lookup";
    // Derived, not passed, so every worker runs a byte-identical command line.
    const std::string coordDir = dataFile + ".coord";
    if (FileUtil::directoryExists(coordDir.c_str()) == false) {
        FileUtil::makeDir(coordDir.c_str());
    }

    // Derived once and recorded, so every worker of the run -- including one that
    // joins after a restart -- cuts the input on the same grid. A different chunk
    // size is a different chunk numbering, and the numbering is what breaks length
    // ties when keys are assigned.
    const size_t chunkSize =
        deriveChunkSize(filenames, par.chunkSize, par.PARAM_CHUNK_SIZE.wasSet);
    const ChunkLayout layout = planChunkLayout(filenames, chunkSize);
    if (layout.total == 0) {
        Debug(Debug::ERROR) << "The input files have no entry\n";
        EXIT(EXIT_FAILURE);
    }
    if (chunkSize != par.chunkSize) {
        Debug(Debug::INFO) << "Using a chunk size of " << chunkSize << " B so the input is "
                           << layout.total << " chunks rather than "
                           << "one work item per " << par.chunkSize << " B\n";
    }

    SharedCounter workerCounter(coordDir + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();
    Debug(Debug::INFO) << "Worker " << workerId << " joined, " << layout.total << " chunks\n";

    // Pass 1: histogram every chunk.
    {
        WorkQueue scanQueue(coordDir + "/scan.queue", static_cast<int64_t>(layout.total));
        runQueue(scanQueue, par.threads, workerId, [&](size_t chunkIdx) {
            const std::string path = chunkHistPath(coordDir, chunkIdx);
            if (FileUtil::fileExists(path.c_str()) == true) {
                return;
            }
            // Boundaries resolved here, by the chunk's own owner, rather than for
            // every chunk in every worker before the stage starts.
            const Chunk chunk = resolveChunk(filenames, layout, chunkSize, chunkIdx);
            ChunkHistogram histogram = scanChunk(filenames[chunk.fileIdx], chunk, chunkIdx);
            histogram.fileIdx = chunk.fileIdx;
            ensureChunkGroupDir(coordDir, chunkIdx);
            histogram.write(path);
        });

    }
    Debug(Debug::INFO) << "Scan pass done\n";

    // Plan: one worker turns the histograms into placements and lays out the
    // output files. The rest wait on the sentinel.
    const std::string planDone = coordDir + "/plan.done";
    {
        FileLock planLock(coordDir + "/plan.lock");
        planLock.lock();
        if (FileUtil::fileExists(planDone.c_str()) == false) {
            std::vector<ChunkHistogram> histograms;
            histograms.reserve(layout.total);
            for (size_t i = 0; i < layout.total; i++) {
                histograms.push_back(ChunkHistogram::read(chunkHistPath(coordDir, i)));
            }

            std::vector<ChunkPlan> plans;
            const LengthRankedTotals totals = buildLengthRankedPlan(histograms, plans);
            if (totals.seqCount == 0) {
                Debug(Debug::ERROR) << "The input files have no entry\n";
                EXIT(EXIT_FAILURE);
            }
            // Refused rather than wrapped. Keys are dense, so the last one is
            // seqCount - 1, and two ceilings apply: DBKeyType, which is a uint32_t
            // unless the build sets MMSEQS_INT64_IDS -- the CMake default is 0 --
            // and the 48-bit id field of KmerRecord. Past either, distinct
            // sequences alias to one key and the reduce groups unrelated sequences
            // together, silently. The 32-bit ceiling bites three orders of
            // magnitude below the scale this command exists for, so it is worth
            // saying which limit was hit.
            const uint64_t keyCeiling =
                std::min<uint64_t>(KmerRecord::MAX_ID, static_cast<uint64_t>(DB_KEY_INVALID) - 1);
            if (totals.seqCount - 1 > keyCeiling) {
                Debug(Debug::ERROR)
                    << "The input holds " << totals.seqCount << " sequences, past the "
                    << (keyCeiling + 1) << " this build can key.\n"
                    << (static_cast<uint64_t>(DB_KEY_INVALID) - 1 < KmerRecord::MAX_ID
                            ? "Rebuild with -DMMSEQS_INT64_IDS=1; the default build uses 32-bit "
                              "database keys.\n"
                            : "The k-mer record carries a 48-bit id, which is the hard limit.\n");
                EXIT(EXIT_FAILURE);
            }
            for (size_t i = 0; i < plans.size(); i++) {
                ensureChunkGroupDir(coordDir, plans[i].chunkIdx);
                plans[i].write(chunkPlanPath(coordDir, plans[i].chunkIdx));
            }

            allocateFile(dataFile, totals.dataBytes);
            allocateFile(hdrDataFile, totals.headerBytes);
            if (par.writeLookup) {
                allocateFile(lookupFile, totals.lookupBytes);
            }
            DenseIndex::createEmpty(dataFile, totals.seqCount, 0, totals.dataBytes,
                                    static_cast<uint32_t>(totals.maxSeqLen + 2));
            DenseIndex::createEmpty(hdrDataFile, totals.seqCount, 0, totals.headerBytes, 0);

            const std::string sourceFile = dataFile + ".source";
            FILE *source = FileUtil::openAndDelete(sourceFile.c_str(), "w");
            for (size_t i = 0; i < filenames.size(); i++) {
                fprintf(source, "%zu\t%s\n", i, FileUtil::baseName(filenames[i]).c_str());
            }
            if (fclose(source) != 0) {
                Debug(Debug::ERROR) << "Cannot close " << sourceFile << "\n";
                EXIT(EXIT_FAILURE);
            }

            // The database type is decided from votes across the whole input, not
            // from whichever chunk happened to be scanned first.
            int dbType = Parameters::DBTYPE_AMINO_ACIDS;
            if (par.dbType == 2) {
                dbType = Parameters::DBTYPE_NUCLEOTIDES;
            } else if (par.dbType == 0 && totals.sampleCount > 0 &&
                       totals.nuclVotes == totals.sampleCount) {
                dbType = Parameters::DBTYPE_NUCLEOTIDES;
            }
            FileUtil::writeFile(coordDir + "/dbtype",
                                reinterpret_cast<const unsigned char *>(&dbType), sizeof(int));

            Debug(Debug::INFO) << "Planned " << totals.seqCount << " sequences, "
                               << totals.dataBytes << " data bytes, longest " << totals.maxSeqLen << "\n";
            FILE *sentinel = FileUtil::openAndDelete(planDone.c_str(), "w");
            fclose(sentinel);
        }
        planLock.unlock();
    }

    // Pass 2: write the sequences.
    {
        const int seqFd = openForWrite(dataFile);
        const int hdrFd = openForWrite(hdrDataFile);
        const int seqIdxFd = openForWrite(DenseIndex::fileName(dataFile));
        const int hdrIdxFd = openForWrite(DenseIndex::fileName(hdrDataFile));
        // Off by request only: at 1e12 sequences the lookup is ~30 TB, and the
        // clustering path never reads it -- keys translate through the header
        // database, which is addressed by the same dense keys.
        const int lookupFd = par.writeLookup ? openForWrite(lookupFile) : -1;

        WorkQueue emitQueue(coordDir + "/emit.queue", static_cast<int64_t>(layout.total));
        runQueue(emitQueue, par.threads, workerId, [&](size_t chunkIdx) {
            const ChunkPlan plan = ChunkPlan::read(chunkPlanPath(coordDir, chunkIdx));
            const Chunk chunk = resolveChunk(filenames, layout, chunkSize, chunkIdx);
            emitChunk(filenames[chunk.fileIdx], chunk, plan,
                      seqFd, hdrFd, seqIdxFd, hdrIdxFd, lookupFd);
        });
        // fsync before the sentinel, so a worker that finalises after a crash
        // cannot read a partially flushed database.
        //
        // Both the fsync and the close are checked. A delayed ENOSPC, an exhausted
        // quota or a shared-filesystem EIO surfaces at exactly these two calls and
        // nowhere earlier, and ignoring them let this worker go on to write
        // finalize.done over a database with holes in it -- which every later run
        // then treats as finished.
        static const char *names[] = {"sequence data", "header data", "sequence index",
                                      "header index", "lookup"};
        const int fds[] = {seqFd, hdrFd, seqIdxFd, hdrIdxFd, lookupFd};
        for (size_t i = 0; i < sizeof(fds) / sizeof(fds[0]); i++) {
            if (fds[i] < 0) {
                continue;
            }
            if (fsync(fds[i]) != 0) {
                const int err = errno;
                Debug(Debug::ERROR) << "Cannot flush the " << names[i] << " of " << dataFile << ": "
                                    << strerror(err) << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (close(fds[i]) != 0) {
                const int err = errno;
                Debug(Debug::ERROR) << "Cannot close the " << names[i] << " of " << dataFile << ": "
                                    << strerror(err) << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }
    Debug(Debug::INFO) << "Emit pass done\n";

    // Finalise: one worker writes the type files and the text indices that the
    // stock tools still expect. The distributed stages read the dense index.
    const std::string finalizeDone = coordDir + "/finalize.done";
    {
        FileLock finalizeLock(coordDir + "/finalize.lock");
        finalizeLock.lock();
        if (FileUtil::fileExists(finalizeDone.c_str()) == false) {
            int dbType = Parameters::DBTYPE_AMINO_ACIDS;
            FILE *typeFile = FileUtil::openFileOrDie((coordDir + "/dbtype").c_str(), "rb", true);
            if (fread(&dbType, sizeof(int), 1, typeFile) != 1) {
                Debug(Debug::ERROR) << "Cannot read the planned database type\n";
                EXIT(EXIT_FAILURE);
            }
            fclose(typeFile);

            // Never compressed. emitChunk pwrites raw bytes at offsets a plan
            // derived from uncompressed lengths fixed, so nothing here compresses
            // anything; passing par.compressed through only *labelled* those raw
            // files as compressed, and readers then tried to decode them. See the
            // note beside createdbparallel's parameter list.
            DBWriter::writeDbtypeFile(dataFile.c_str(), dbType, false);
            DBWriter::writeDbtypeFile(hdrDataFile.c_str(), Parameters::DBTYPE_GENERIC_DB, false);
            // Opt-in. No stage of this pipeline reads a text index -- they all
            // address entries through the dense .index.bin -- and it exists only so
            // stock MMseqs2 tools can open the database. Measured, the snprintf and
            // fwrite loop costs 177 ns per line writing to /dev/null, so at 1e12
            // sequences the two databases are 2e12 lines: ~98 h single-threaded on
            // one node, holding the finalize lock, plus ~36 TB of a 1 PB scratch
            // budget for output nothing in the run consumes.
            if (par.writeTextIndex) {
                DenseIndex::writeTextIndex(dataFile);
                DenseIndex::writeTextIndex(hdrDataFile);
            } else {
                Debug(Debug::INFO) << "Skipping the stock-compatible text indices "
                                   << "(--write-text-index 0); the dense .index.bin is what the "
                                   << "distributed stages read\n";
            }

            Debug(Debug::INFO) << "Database type: " << Parameters::getDbTypeName(dbType) << "\n";
            FILE *sentinel = FileUtil::openAndDelete(finalizeDone.c_str(), "w");
            fclose(sentinel);
        }
        finalizeLock.unlock();
    }

    return EXIT_SUCCESS;
}
