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

std::string chunkHistPath(const std::string &coordDir, size_t chunkIdx) {
    return coordDir + "/chunk." + SSTR(chunkIdx) + ".hist";
}

std::string chunkPlanPath(const std::string &coordDir, size_t chunkIdx) {
    return coordDir + "/chunk." + SSTR(chunkIdx) + ".plan";
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

// Splits every input file into chunk-size pieces aligned to record boundaries.
// Computed identically and independently by every worker, so the chunk numbering
// -- and therefore the key assignment that breaks length ties by chunk index --
// never depends on who is running.
std::vector<Chunk> planChunks(const std::vector<std::string> &filenames, size_t chunkSize) {
    std::vector<Chunk> chunks;
    for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
        const size_t fileSize = FileUtil::getFileSize(filenames[fileIdx]);
        if (fileSize == 0) {
            continue;
        }
        const int fd = open(filenames[fileIdx].c_str(), O_RDONLY);
        if (fd < 0) {
            Debug(Debug::ERROR) << "Cannot open " << filenames[fileIdx] << ": " << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        size_t begin = 0;
        while (begin < fileSize) {
            const size_t nominalEnd = std::min(begin + chunkSize, fileSize);
            const size_t end = findRecordStart(fd, nominalEnd, fileSize);
            if (end > begin) {
                Chunk chunk = {fileIdx, begin, end};
                chunks.push_back(chunk);
            }
            if (end <= begin) {
                break;
            }
            begin = end;
        }
        close(fd);
    }
    return chunks;
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
// The claim is recorded against the *process*, not the thread, because that is
// the unit that dies: if this worker is killed, every item its threads held must
// become re-claimable together once the lease expires. Passing a thread index
// here instead would also make identities collide across worker processes, since
// every process numbers its threads from zero.
template <typename Body>
void runQueue(WorkQueue &queue, int threads, int64_t workerId, Body body) {
    bool stalled = false;
#pragma omp parallel num_threads(threads)
    {
        // drain() rather than a plain claim loop: it goes back to claiming after
        // waiting, so items abandoned by a crashed worker are picked up once their
        // leases expire instead of being waited on by everyone forever.
        if (queue.drain(workerId, body) == false) {
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
    }

    const std::string hdrDataFile = dataFile + "_h";
    const std::string lookupFile = dataFile + ".lookup";
    // Derived, not passed, so every worker runs a byte-identical command line.
    const std::string coordDir = dataFile + ".coord";
    if (FileUtil::directoryExists(coordDir.c_str()) == false) {
        FileUtil::makeDir(coordDir.c_str());
    }

    const std::vector<Chunk> chunks = planChunks(filenames, par.chunkSize);
    if (chunks.empty()) {
        Debug(Debug::ERROR) << "The input files have no entry\n";
        EXIT(EXIT_FAILURE);
    }

    SharedCounter workerCounter(coordDir + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();
    Debug(Debug::INFO) << "Worker " << workerId << " joined, " << chunks.size() << " chunks\n";

    // Pass 1: histogram every chunk.
    {
        WorkQueue scanQueue(coordDir + "/scan.queue", static_cast<int64_t>(chunks.size()));
        runQueue(scanQueue, par.threads, workerId, [&](size_t chunkIdx) {
            const std::string path = chunkHistPath(coordDir, chunkIdx);
            if (FileUtil::fileExists(path.c_str()) == true) {
                return;
            }
            ChunkHistogram histogram = scanChunk(filenames[chunks[chunkIdx].fileIdx],
                                                 chunks[chunkIdx], chunkIdx);
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
            histograms.reserve(chunks.size());
            for (size_t i = 0; i < chunks.size(); i++) {
                histograms.push_back(ChunkHistogram::read(chunkHistPath(coordDir, i)));
            }

            std::vector<ChunkPlan> plans;
            const LengthRankedTotals totals = buildLengthRankedPlan(histograms, plans);
            if (totals.seqCount == 0) {
                Debug(Debug::ERROR) << "The input files have no entry\n";
                EXIT(EXIT_FAILURE);
            }
            for (size_t i = 0; i < plans.size(); i++) {
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

        WorkQueue emitQueue(coordDir + "/emit.queue", static_cast<int64_t>(chunks.size()));
        runQueue(emitQueue, par.threads, workerId, [&](size_t chunkIdx) {
            const ChunkPlan plan = ChunkPlan::read(chunkPlanPath(coordDir, chunkIdx));
            emitChunk(filenames[chunks[chunkIdx].fileIdx], chunks[chunkIdx], plan,
                      seqFd, hdrFd, seqIdxFd, hdrIdxFd, lookupFd);
        });
        // fsync before the sentinel, so a worker that finalises after a crash
        // cannot read a partially flushed database.
        fsync(seqFd);
        fsync(hdrFd);
        fsync(seqIdxFd);
        fsync(hdrIdxFd);
        if (lookupFd >= 0) {
            fsync(lookupFd);
        }
        close(seqFd);
        close(hdrFd);
        close(seqIdxFd);
        close(hdrIdxFd);
        if (lookupFd >= 0) {
            close(lookupFd);
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

            DBWriter::writeDbtypeFile(dataFile.c_str(), dbType, par.compressed);
            DBWriter::writeDbtypeFile(hdrDataFile.c_str(), Parameters::DBTYPE_GENERIC_DB, par.compressed);
            DenseIndex::writeTextIndex(dataFile);
            DenseIndex::writeTextIndex(hdrDataFile);

            Debug(Debug::INFO) << "Database type: " << Parameters::getDbTypeName(dbType) << "\n";
            FILE *sentinel = FileUtil::openAndDelete(finalizeDone.c_str(), "w");
            fclose(sentinel);
        }
        finalizeLock.unlock();
    }

    return EXIT_SUCCESS;
}
