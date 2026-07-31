/*
 * kmermatcherparallel -- the map half of the distributed linclust.
 *
 * Scans the sequence database once and writes every extracted k-mer into the
 * bucket file of the partition its hash selects. Many worker processes run the
 * same command line and coordinate only through files in a shared directory; a
 * worker may join late, die, or be restarted without losing or duplicating work.
 *
 * Why this exists instead of running stock kmermatcher with more splits: a split
 * re-reads the whole database and re-extracts every k-mer, then discards the
 * ones outside its hash range (kmermatcher.cpp:322). At 1e12 sequences on a 2 TB
 * node that is ~68 passes over the most expensive stage in the pipeline. Here the
 * extraction is paid for exactly once and the partitioning is what carries k-mers
 * to the reduce, so the cost is one scan plus one shuffle.
 *
 * The partitioning is lossless, not an approximation: the 16-bit score is a pure
 * function of the k-mer, so equal k-mers always land in the same partition and a
 * partition can be grouped in isolation (see KmerPartition.h).
 *
 * Memory per worker is bounded by the bucket buffers and the per-thread k-mer
 * scratch, not by the database size: the database is read through DenseIndex in
 * key ranges, so no worker ever holds an index entry per sequence.
 */
#include "Command.h"
#include "Debug.h"
#include "DBReader.h"
#include "DenseIndex.h"
#include "FileUtil.h"
#include "CandidateEdge.h"
#include "KmerPartition.h"
#include "NucleotideMatrix.h"
#include "ParallelCoordination.h"
#include "Parameters.h"
#include "ReducedMatrix.h"
#include "SubstitutionMatrix.h"
#include "Util.h"
#include "kmermatcher.h"

#include <climits>
#include <cstdint>
#include <cmath>
#include <string>

namespace {

// Every file a sequence database actually occupies, not just its data file.
// createdbparallel writes the data, the dense index, the text index, the headers
// with their own two indices, .lookup and .source; sizing against the data file
// alone undercounts by ~1.64x on real data.
static uint64_t databaseFootprint(const std::string &db) {
    static const char *suffixes[] = {"", ".index", ".index.bin", ".dbtype", ".lookup", ".source",
                                     "_h", "_h.index", "_h.index.bin", "_h.dbtype"};
    uint64_t total = 0;
    for (size_t i = 0; i < sizeof(suffixes) / sizeof(suffixes[0]); i++) {
        const std::string path = db + suffixes[i];
        if (FileUtil::fileExists(path.c_str())) {
            total += FileUtil::getFileSize(path);
        }
    }
    return total;
}

// Work items are contiguous key ranges, so an item is also one contiguous read
// of the data file. Their size is derived rather than fixed: a fixed size that
// suits 1e12 sequences leaves a 1e6-sequence database with a single item, which
// no second worker can ever claim.
//
// Aim for many more items than there will plausibly be workers, so late joiners
// find work and a dead worker costs one item; bound the size from below so the
// fixed per-item cost is amortised, and from above so the largest runs keep item
// counts and lease traffic sane.
//
// The floor is the load-bearing one. An item opens a DBReader over its key range,
// which re-maps the data file, and that costs ~50 ms independent of how big the
// item is. Measured on 1M sequences: 1000-sequence items took 61 s against 7.8 s
// for the same work in one item. At 50k the fixed cost is a few percent, and at
// the scales this stage is for, items sit at the ceiling anyway.
const uint64_t TARGET_ITEM_COUNT = 4096;
const uint64_t MIN_SEQUENCES_PER_ITEM = 50000;
const uint64_t MAX_SEQUENCES_PER_ITEM = 1000000;

uint64_t deriveSequencesPerItem(uint64_t entryCount) {
    const uint64_t target = entryCount / TARGET_ITEM_COUNT;
    if (target < MIN_SEQUENCES_PER_ITEM) {
        return MIN_SEQUENCES_PER_ITEM;
    }
    if (target > MAX_SEQUENCES_PER_ITEM) {
        return MAX_SEQUENCES_PER_ITEM;
    }
    return target;
}

// What the first worker records so every later worker can check it is joining
// the run it thinks it is. Repartitioning halfway through -- which is what a
// changed --split-memory-limit or --scratch-budget would silently do -- would
// scatter equal k-mers across different partitions and quietly lose edges.
struct ShuffleManifest {
    uint64_t entryCount;
    uint64_t partitionCount;
    uint64_t waveCount;
    uint64_t kmerSize;

    void write(const std::string &path) const {
        FILE *file = FileUtil::openAndDelete(path.c_str(), "w");
        fprintf(file, "entryCount\t%zu\n", (size_t)entryCount);
        fprintf(file, "partitionCount\t%zu\n", (size_t)partitionCount);
        fprintf(file, "waveCount\t%zu\n", (size_t)waveCount);
        fprintf(file, "kmerSize\t%zu\n", (size_t)kmerSize);
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    static ShuffleManifest read(const std::string &path) {
        ShuffleManifest manifest;
        FILE *file = FileUtil::openFileOrDie(path.c_str(), "r", true);
        char name[64];
        size_t value;
        manifest.entryCount = 0;
        manifest.partitionCount = 0;
        manifest.waveCount = 0;
        manifest.kmerSize = 0;
        while (fscanf(file, "%63s\t%zu\n", name, &value) == 2) {
            const std::string key = name;
            if (key == "entryCount") {
                manifest.entryCount = value;
            } else if (key == "partitionCount") {
                manifest.partitionCount = value;
            } else if (key == "waveCount") {
                manifest.waveCount = value;
            } else if (key == "kmerSize") {
                manifest.kmerSize = value;
            }
        }
        fclose(file);
        return manifest;
    }

    void requireMatches(const ShuffleManifest &other, const std::string &path) const {
        if (entryCount != other.entryCount || partitionCount != other.partitionCount ||
            waveCount != other.waveCount || kmerSize != other.kmerSize) {
            Debug(Debug::ERROR)
                << "This worker derived a different k-mer shuffle than the one already in "
                << "progress (" << path << "): " << partitionCount << " partitions, "
                << waveCount << " waves, k=" << kmerSize << " over " << entryCount
                << " sequences, against " << other.partitionCount << " partitions, "
                << other.waveCount << " waves, k=" << other.kmerSize << " over "
                << other.entryCount << " sequences.\n"
                << "Every worker must run the same command line on the same database.\n";
            EXIT(EXIT_FAILURE);
        }
    }
};

BaseMatrix *createSubstitutionMatrix(Parameters &par, int dbType) {
    if (Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES)) {
        return new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    }
    if (par.alphabetSize.values.aminoacid() == 21) {
        return new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    }
    SubstitutionMatrix sMat(par.scoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2f);
    return new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, sMat.aa2num, sMat.num2aa,
                             sMat.alphabetSize, par.alphabetSize.values.aminoacid(), 2.0);
}

// K-mers selected per sequence, by the same rule fillKmerPositionArray applies
// (kmersPerSequence - 1 + scale * L), plus the identity k-mer every sequence
// contributes. Computed from the average length rather than measured, because
// measuring would mean a pass over the database this stage does not otherwise
// need; it only sizes the shuffle, and the placement stays exact either way.
unsigned int estimateKmersPerSequence(Parameters &par, int dbType, uint64_t residues,
                                      uint64_t entryCount) {
    const float scale = Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES)
                            ? par.kmersPerSequenceScale.values.nucleotide()
                            : par.kmersPerSequenceScale.values.aminoacid();
    const double averageLength = static_cast<double>(residues) / static_cast<double>(entryCount);
    const double selected = static_cast<double>(par.kmersPerSequence) - 1.0 + scale * averageLength;
    return static_cast<unsigned int>(std::max(1.0, std::ceil(selected))) + 1;
}

// Extracts one key range into the bucket files.
//
// T is the position/length type of the staged KmerPosition, picked from the
// longest sequence exactly as stock kmermatcher picks it. It never reaches the
// bucket file, which always stores both fields as uint16.
template <typename T>
void scanKeyRange(const std::string &seqDb, DBKeyType keyFrom, DBKeyType keyTo, int dbType,
                  Parameters &par, BaseMatrix *subMat, KmerBucketWriter &writer,
                  const KmerPartitioner &partitioner) {
    DenseIndex::Info rangeInfo;
    DBReader<DBKeyType>::Index *index = DenseIndex::loadRange(seqDb, keyFrom, keyTo, &rangeInfo);
    if (rangeInfo.entryCount == 0) {
        delete[] index;
        return;
    }

    DBReader<DBKeyType> reader(index, rangeInfo.entryCount, rangeInfo.dataSize,
                               static_cast<DBKeyType>(keyFrom + rangeInfo.entryCount - 1), dbType,
                               rangeInfo.maxSeqLen, par.threads);
    reader.setDataFile(seqDb.c_str());
    reader.open(DBReader<DBKeyType>::NOSORT);

    // The k-mer array is never touched in bucketing mode -- every selected k-mer
    // goes straight to a bucket -- so there is nothing to allocate for it.
    if (Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES)) {
        fillKmerPositionArray<Parameters::DBTYPE_NUCLEOTIDES, T, true, true>(
            NULL, SIZE_MAX, reader, par, subMat, true, 0, SIZE_MAX, NULL, &writer, &partitioner);
    } else {
        fillKmerPositionArray<Parameters::DBTYPE_AMINO_ACIDS, T, true, true>(
            NULL, SIZE_MAX, reader, par, subMat, true, 0, SIZE_MAX, NULL, &writer, &partitioner);
    }

    reader.close();
    // DBReader does not own an externally supplied index.
    delete[] index;
}

}  // namespace

int kmermatcherparallel(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setLinearFilterDefault(&par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    const std::string seqDb = par.db1;
    const std::string kmerDir = par.db2;

    if (DenseIndex::exists(seqDb) == false) {
        Debug(Debug::ERROR) << "No dense index next to " << seqDb << ". The distributed stages "
                            << "address sequences by key without a resident index, which needs "
                            << "the dense keys createdbparallel produces.\n";
        EXIT(EXIT_FAILURE);
    }
    const DenseIndex::Info info = DenseIndex::readInfo(seqDb);
    if (info.entryCount == 0) {
        Debug(Debug::ERROR) << "Database " << seqDb << " is empty\n";
        EXIT(EXIT_FAILURE);
    }
    // The packed bucket record carries pos and seqLen in 16 bits, so a longer
    // sequence would silently wrap rather than merely lose precision.
    if (info.maxSeqLen > USHRT_MAX) {
        Debug(Debug::ERROR) << "Longest sequence is " << info.maxSeqLen
                            << " residues; the packed k-mer record holds positions up to "
                            << USHRT_MAX << ".\n";
        EXIT(EXIT_FAILURE);
    }
    if (par.adjustKmerLength) {
        // The adjusted length is a property of the whole database, and a worker
        // scanning one key range cannot arrive at the value the others will.
        Debug(Debug::ERROR) << "--adjust-kmer-len is not supported by the distributed map\n";
        EXIT(EXIT_FAILURE);
    }

    const int dbType = FileUtil::parseDbType(seqDb.c_str());
    // Every entry stores its residues plus a newline and DBWriter's terminating
    // null, so the residue count follows from the data size without a scan.
    const uint64_t residues = info.dataSize - 2 * info.entryCount;
    setKmerLengthAndAlphabet(par, residues, dbType);
    par.printParameters(command.cmd, argc, argv, *command.params);
    Debug(Debug::INFO) << "Database size: " << info.entryCount << " sequences, " << residues
                       << " residues\n";

    const unsigned int kmersPerSequence =
        estimateKmersPerSequence(par, dbType, residues, info.entryCount);

    // What the shuffle may NOT use, which is everything that is on the scratch
    // filesystem while it runs and outlives it.
    //
    // This used to be `info.dataSize` -- the sequence *data file* alone. Measured
    // on MGnify, the real database footprint is 1.64x that once the dense index,
    // the text index, the headers, .lookup and .source are counted, and none of
    // the pipeline's own output was counted at all. The result was that a 1e11 run
    // handed the true 100 TB ceiling derived a single wave and peaked at ~186 TB.
    //
    // Two corrections. First, prefer a *measured* occupied figure passed by the
    // caller (--scratch-used): the workflow knows what is already on disk, which
    // is the only way pass 2 can account for what pass 1 left behind, and is why
    // one budget knob can now size both passes. Second, reserve the candidate
    // edges this stage's own reduce will write. Their volume is bounded, not
    // guessed: a candidate edge needs a k-mer to have found it, so there can be no
    // more edge records than k-mer records, and each is sizeof(CandidateEdge)
    // against sizeof(KmerRecord).
    const uint64_t totalKmerBytes =
        static_cast<uint64_t>(info.entryCount) * kmersPerSequence * sizeof(KmerRecord);
    const uint64_t projectedEdgeBytes =
        totalKmerBytes / sizeof(KmerRecord) * sizeof(CandidateEdge);
    uint64_t occupied = par.scratchUsed > 0 ? static_cast<uint64_t>(par.scratchUsed)
                                            : databaseFootprint(seqDb);
    const uint64_t persistentBytes = occupied + projectedEdgeBytes;
    if (par.scratchBudget > 0) {
        Debug(Debug::INFO) << "Scratch budget " << par.scratchBudget << " B; already occupied "
                           << occupied << " B; reserving " << projectedEdgeBytes
                           << " B for candidate edges\n";
    }
    const KmerShuffleSizing sizing =
        deriveKmerShuffleSizing(info.entryCount, kmersPerSequence, par.scratchBudget,
                                persistentBytes, Util::computeMemory(par.splitMemoryLimit));
    Debug(Debug::INFO) << "K-mer shuffle: " << sizing.partitionCount << " partitions, "
                       << sizing.waveCount << " wave(s), " << sizing.totalKmerBytes
                       << " k-mer bytes, " << sizing.bytesPerPartition << " bytes per partition\n";
    // Which slice of the partition space this invocation writes. Waves exist
    // because the whole shuffle need not fit at once: each re-extracts every
    // k-mer but keeps only its own slice, so peak scratch is divided by the wave
    // count at the cost of that many extraction passes.
    unsigned int waveFrom = 0;
    unsigned int waveTo = sizing.partitionCount;
    if (sizing.waveCount > 1) {
        if (par.kmerWave < 0) {
            Debug(Debug::ERROR) << "The scratch budget of " << par.scratchBudget << " bytes needs "
                                << sizing.waveCount << " extraction waves. Run waves 0.."
                                << (sizing.waveCount - 1) << " with --kmer-wave, reducing and "
                                << "deleting each wave's buckets before starting the next, or "
                                << "raise --scratch-budget to fit a single wave.\n";
            EXIT(EXIT_FAILURE);
        }
        if (static_cast<unsigned int>(par.kmerWave) >= sizing.waveCount) {
            Debug(Debug::ERROR) << "--kmer-wave " << par.kmerWave << " is out of range; this run "
                                << "has " << sizing.waveCount << " waves\n";
            EXIT(EXIT_FAILURE);
        }
        // Exact: the sizing guarantees the wave count divides P.
        const unsigned int perWave = sizing.partitionCount / sizing.waveCount;
        waveFrom = static_cast<unsigned int>(par.kmerWave) * perWave;
        waveTo = waveFrom + perWave;
        Debug(Debug::INFO) << "Wave " << par.kmerWave << " of " << sizing.waveCount
                           << ": partitions [" << waveFrom << ", " << waveTo << ")\n";
    } else if (par.kmerWave > 0) {
        Debug(Debug::ERROR) << "--kmer-wave " << par.kmerWave << " given but this run needs only "
                            << "one wave\n";
        EXIT(EXIT_FAILURE);
    }

    if (FileUtil::directoryExists(kmerDir.c_str()) == false) {
        FileUtil::makeDir(kmerDir.c_str());
    }
    // Derived rather than passed, so every worker runs a byte-identical command line.
    const std::string coordDir = kmerDir + "/coord";
    if (FileUtil::directoryExists(coordDir.c_str()) == false) {
        FileUtil::makeDir(coordDir.c_str());
    }

    ShuffleManifest manifest;
    manifest.entryCount = info.entryCount;
    manifest.partitionCount = sizing.partitionCount;
    manifest.waveCount = sizing.waveCount;
    manifest.kmerSize = static_cast<uint64_t>(par.kmerSize);
    const std::string manifestPath = coordDir + "/shuffle.info";
    {
        FileLock manifestLock(coordDir + "/shuffle.lock");
        manifestLock.lock();
        if (FileUtil::fileExists(manifestPath.c_str()) == true) {
            manifest.requireMatches(ShuffleManifest::read(manifestPath), manifestPath);
        } else {
            KmerBucketWriter::createLayout(kmerDir, sizing.partitionCount);
            manifest.write(manifestPath);
        }
        manifestLock.unlock();
    }

    SharedCounter workerCounter(coordDir + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();
    const uint64_t sequencesPerItem = deriveSequencesPerItem(info.entryCount);
    const int64_t itemCount =
        static_cast<int64_t>((info.entryCount + sequencesPerItem - 1) / sequencesPerItem);
    Debug(Debug::INFO) << "Worker " << workerId << " joined, " << itemCount << " key ranges of "
                       << sequencesPerItem << " sequences\n";

    BaseMatrix *subMat = createSubstitutionMatrix(par, dbType);
    KmerPartitioner partitioner(sizing.partitionCount);
    // One writer for the whole process, shared by every thread and kept open
    // across items so bucket files are appended to rather than reopened.
    KmerBucketWriter writer(kmerDir, sizing.partitionCount, "w" + SSTR(workerId),
                            1024 * 1024 * 1024, waveFrom, waveTo);

    {
        WorkQueue queue(coordDir + "/scan." + SSTR(par.kmerWave < 0 ? 0 : par.kmerWave) +
                            ".queue", itemCount);
        // Claimed one item at a time by the process rather than by each thread:
        // the extraction inside an item is already threaded, and nesting a second
        // parallel region inside a claiming one would oversubscribe the node.
        const bool finished = queue.drain(workerId, [&](size_t item) {
            const DBKeyType keyFrom = static_cast<DBKeyType>(item * sequencesPerItem);
            const DBKeyType keyTo = static_cast<DBKeyType>(
                std::min<uint64_t>((item + 1) * sequencesPerItem, info.entryCount));
            if (info.maxSeqLen < SHRT_MAX) {
                scanKeyRange<short>(seqDb, keyFrom, keyTo, dbType, par, subMat, writer, partitioner);
            } else {
                scanKeyRange<int>(seqDb, keyFrom, keyTo, dbType, par, subMat, writer, partitioner);
            }
            // Before drain() records the item as done, so a worker that dies can
            // never leave an item marked complete whose k-mers were still buffered.
            // One flush per item keeps the writes large: item size and P grow
            // together, so a partition takes ~0.5 MB per flush at the 100B target.
            writer.flushAll();
        });
        if (finished == false) {
            Debug(Debug::ERROR) << "Map stage stalled: work remains but no item is claimable\n";
            EXIT(EXIT_FAILURE);
        }
    }

    writer.close();
    Debug(Debug::INFO) << "Worker " << workerId << " wrote " << writer.getRecordCount()
                       << " k-mers\n";
    delete subMat;

    return EXIT_SUCCESS;
}
