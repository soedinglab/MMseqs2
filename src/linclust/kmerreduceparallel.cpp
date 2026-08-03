/*
 * kmerreduceparallel -- the reduce half of the distributed linclust.
 *
 * Reads back one k-mer partition written by kmermatcherparallel, groups it, and
 * writes the resulting (representative, member) candidate edges. Many worker
 * processes run the same command line and claim partitions from a shared work
 * queue, exactly as the map stage claims key ranges.
 *
 * Why grouping a partition in isolation gives the same answer as grouping the
 * whole database: the greedy only ever compares *equal* k-mers, and the partition
 * of a k-mer is a pure function of the k-mer (see KmerPartition.h), so every
 * occurrence of a k-mer is in the partition being read and nowhere else. The
 * grouping call is stock's `assignGroup`, unchanged -- this command supplies it
 * with one partition where stock supplies it with one split.
 *
 * Two things here exist because the stock reduce cannot run at 1e12 sequences:
 *   - lengths come from the k-mer record, not from `seqkey_to_len[dbKeySize]`
 *     (kmermatcher.cpp:1236), which is sized by key space;
 *   - output is packed binary edges rather than a prefilter DB, whose per-key
 *     index is likewise per-key state no node can hold.
 * Neither the sequence database nor any per-key table is opened at all: a
 * partition is self-contained, which is the whole point of carrying seqLen and
 * the adjacent residues in the record.
 */
#include "Command.h"
#include "CandidateEdge.h"
#include "Debug.h"
#include "DenseIndex.h"
#include "FastSort.h"
#include "FileUtil.h"
#include "KmerPartition.h"
#include "NucleotideMatrix.h"
#include "ParallelCoordination.h"
#include "Parameters.h"
#include "ReducedMatrix.h"
#include "SubstitutionMatrix.h"
#include "Util.h"
#include "kmermatcher.h"

#ifdef OPENMP
#include <omp.h>
#endif

#include <cerrno>
#include <climits>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

#include <unistd.h>

namespace {

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

// Streams a partition's shards straight into the KmerPosition array.
//
// Reading into a vector<KmerRecord> first and converting afterwards would hold
// both representations at once, which is 24 extra bytes per k-mer at the moment
// of peak memory -- ~50% more for a stage whose partition size is chosen to just
// fit a node.
template <typename T>
size_t readPartitionAsPositions(const std::string &kmerDir, unsigned int partition,
                                KmerPosition<T, true, true> *out, size_t capacity) {
    const std::vector<std::string> shards = KmerBucketReader::shardFiles(kmerDir, partition);
    const size_t blockRecords = 1024 * 1024;
    std::vector<KmerRecord> block(blockRecords);
    size_t filled = 0;
    for (size_t i = 0; i < shards.size(); i++) {
        // Not openFileOrDie: a worker whose lease lapsed can still be here while
        // the workers that finished the wave unlink its shards. Its results are
        // discarded by block header regardless, so the shard vanishing under it is
        // a race it should survive rather than a reason to fail the whole stage.
        FILE *file = fopen(shards[i].c_str(), "rb");
        if (file == NULL) {
            if (errno == ENOENT) {
                continue;
            }
            Debug(Debug::ERROR) << "Cannot open k-mer bucket " << shards[i] << ": "
                                << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        while (true) {
            const size_t got = fread(block.data(), sizeof(KmerRecord), blockRecords, file);
            if (got == 0) {
                // A short read is EOF only if nothing went wrong. Treating an I/O
                // error as EOF silently groups the partition with fewer k-mers,
                // and the missing edges are indistinguishable from "this k-mer had
                // no partner" -- a wrong answer with no diagnostic.
                if (ferror(file)) {
                    Debug(Debug::ERROR) << "Cannot read k-mer bucket " << shards[i] << ": "
                                        << strerror(errno) << "\n";
                    EXIT(EXIT_FAILURE);
                }
                break;
            }
            if (filled + got > capacity) {
                Debug(Debug::ERROR) << "Partition " << partition << " holds more k-mers than the "
                                    << "shard sizes reported. A shard was written while this was "
                                    << "reading it.\n";
                EXIT(EXIT_FAILURE);
            }
            for (size_t r = 0; r < got; r++) {
                kmerRecordToPosition(block[r], out[filled + r]);
            }
            filled += got;
        }
        fclose(file);
    }
    return filled;
}

// Drops records that are byte-identical to their predecessor.
//
// The map is not idempotent across a crash: a worker killed mid-item may have
// flushed part of it, and the worker that redoes the item writes the same records
// into a different shard. Removing them here is exact rather than a heuristic --
// (kmer, id, pos) names one k-mer occurrence in one sequence, so identical
// records can only be the same occurrence written twice. Duplicates are adjacent
// after the sort because every field they are ordered by is equal.
template <typename T>
size_t dropDuplicates(KmerPosition<T, true, true> *positions, size_t count) {
    if (count == 0) {
        return 0;
    }
    size_t out = 1;
    for (size_t i = 1; i < count; i++) {
        if (memcmp(&positions[i], &positions[out - 1], sizeof(positions[0])) != 0) {
            positions[out] = positions[i];
            out++;
        }
    }
    return out;
}

// Splits the sorted array into per-thread ranges that never cut a k-mer group,
// so each thread's greedy sees whole groups. Same rule as doComputation
// (kmermatcher.cpp:1006-1047).
template <typename T>
void buildThreadOffsets(KmerPosition<T, true, true> *positions, size_t count, int threads,
                        bool isNucleotide, std::vector<size_t> &threadOffsets) {
    threadOffsets.clear();
    threadOffsets.push_back(0);
    const size_t splitSize = count / threads;
    for (int thread = 1; thread < threads; thread++) {
        size_t prevKmer = positions[thread * splitSize].kmer;
        if (prevKmer == SIZE_MAX) {
            for (int i = thread; i < threads; i++) {
                threadOffsets.push_back(count);
            }
            break;
        }
        if (isNucleotide) {
            prevKmer = BIT_SET(prevKmer, 63);
        }
        bool wasSet = false;
        for (size_t pos = thread * splitSize; pos < count; pos++) {
            size_t currKmer = positions[pos].kmer;
            if (isNucleotide) {
                currKmer = BIT_SET(currKmer, 63);
            }
            if (prevKmer != currKmer) {
                wasSet = true;
                threadOffsets.push_back(pos);
                break;
            }
        }
        if (wasSet == false) {
            for (int i = thread; i < threads; i++) {
                threadOffsets.push_back(count);
            }
            break;
        }
    }
    threadOffsets.push_back(count);
}


// Collapses one grouping round's output into (representative, member, diagonal)
// edges, counting how many k-mers support each diagonal.
//
// It deliberately does **not** pick a best diagonal here. Stock chooses the
// diagonal the most k-mers agree on *globally* (kmermatcher.cpp:1913-1925), and a
// partition only ever sees a fraction of a pair's k-mers. Collapsing to one
// diagonal per partition would discard the other diagonals' counts before the
// align stage can add them up -- measured, that reproduces stock's diagonal for
// only 98.64% of pairs. Emitting one edge per distinct diagonal keeps the counts
// intact so the global merge is exact, and costs ~12% more edges because a
// colinear pair puts all its k-mers on one diagonal anyway.
template <typename T>
void collectRoundEdges(KmerPosition<T, false, true> *grouped, size_t writePos, bool isNucleotide,
                       std::vector<CandidateEdge> &edges) {
    CandidateEdge edge;
    size_t i = 0;
    while (i < writePos && grouped[i].kmer != SIZE_MAX) {
        const size_t repRaw = grouped[i].kmer;
        const DBKeyType member = grouped[i].id;
        size_t j = i;
        while (j < writePos && grouped[j].kmer == repRaw && grouped[j].id == member) {
            const T diagonal = grouped[j].pos;
            size_t count = 0;
            while (j < writePos && grouped[j].kmer == repRaw && grouped[j].id == member &&
                   grouped[j].pos == diagonal) {
                count++;
                j++;
            }
            size_t rep = repRaw;
            // Stock keeps the strand in bit 63 of the representative key; a 48-bit
            // key has no room for it, so it moves into its own field.
            edge.reverseStrand = 0;
            if (isNucleotide) {
                edge.reverseStrand = BIT_CHECK(rep, 63) == false;
                rep = BIT_CLEAR(rep, 63);
            }
            // A sequence being its own representative carries no information: keys
            // are dense, so any key never appearing as a member is its own
            // representative by construction.
            if (static_cast<DBKeyType>(rep) == member) {
                continue;
            }
            edge.setRep(static_cast<uint64_t>(rep));
            edge.setMember(static_cast<uint64_t>(member));
            edge.diagonal = static_cast<int16_t>(diagonal);
            edge.score = static_cast<uint16_t>(std::min<size_t>(count, 65535));
            edges.push_back(edge);
        }
        i = j;
    }
}

// Orders edges so copies of the same (pair, diagonal) are adjacent.
bool compareEdge(const CandidateEdge &a, const CandidateEdge &b) {
    const uint64_t ra = a.getRep(), rb = b.getRep();
    if (ra != rb) return ra < rb;
    const uint64_t ma = a.getMember(), mb = b.getMember();
    if (ma != mb) return ma < mb;
    if (a.diagonal != b.diagonal) return a.diagonal < b.diagonal;
    // Strand belongs in the key: two opposite-strand edges on the same diagonal
    // are different alignments, and collapsing them would sum their support.
    // Only reachable for nucleotide input, which alignparallel rejects today.
    return a.reverseStrand < b.reverseStrand;
}



// Groups one partition and writes its edges.
template <typename T>
uint64_t reducePartition(const std::string &kmerDir,
                         unsigned int partition, int dbType, Parameters &par, BaseMatrix *subMat,
                         EdgeBucketWriter &writer, uint64_t bucketSpan) {
    const bool isNucleotide = Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);

    const uint64_t recordCount = KmerBucketReader::countRecords(kmerDir, partition);
    if (recordCount == 0) {
        return 0;
    }

    // One slot past the end for the SIZE_T_MAX sentinel assignGroup stops on,
    // matching initKmerPositionMemory (kmermatcher.cpp:43).
    KmerPosition<T, true, true> *positions =
        new (std::nothrow) KmerPosition<T, true, true>[recordCount + 1];
    Util::checkAllocation(positions, "Cannot allocate the k-mer partition");

    size_t count = readPartitionAsPositions<T>(kmerDir, partition, positions, recordCount);
    if (isNucleotide) {
        SORT_PARALLEL(positions, positions + count,
                      KmerPosition<T, true, true>::compareRepSequenceAndIdAndPosReverse);
    } else {
        SORT_PARALLEL(positions, positions + count,
                      KmerPosition<T, true, true>::compareRepSequenceAndIdAndPos);
    }
    count = dropDuplicates<T>(positions, count);
    memset(&positions[count], 0xFF, sizeof(positions[0]));

    std::vector<size_t> threadOffsets;
    buildThreadOffsets<T>(positions, count, par.threads, isNucleotide, threadOffsets);

    // assignGroup writes the grouped pairs here rather than over its input, so a
    // round can be re-run over the k-mers rather than over the previous result.
    KmerPosition<T, false, true> *grouped =
        new (std::nothrow) KmerPosition<T, false, true>[count + 1];
    Util::checkAllocation(grouped, "Cannot allocate the grouping output");

    // linclust v2's rounds, ported.
    //
    // Round 0 picks the longest sequence of each k-mer group as its centre; every
    // later round re-picks it by adjacency agreement. The rounds are *stateful*:
    // assignGroup swaps the chosen centre to the front of its group inside the
    // k-mer array (kmermatcher.cpp:607,638), so each round starts from the previous
    // round's arrangement and finds pairs the earlier ones did not. Measured on
    // stock, the rounds are strictly additive -- 3 rounds is a superset of 1,
    // adding 90,628 pairs and losing none.
    //
    // That mutation is the whole reason this ports for free: all the round-to-round
    // state lives in the k-mer array, which is partition-local. Nothing global is
    // needed, so the rounds cost the same here as they do in stock.
    const int rounds = par.includeAdjacency ? 1 + par.adjIteration : 1;
    std::vector<CandidateEdge> edges;
    for (int round = 0; round < rounds; round++) {
        const AssignGroupMask mask =
            (round == 0) ? AssignGroupFeature::Default : AssignGroupFeature::AdjacentSeq;
        size_t writePos = 0;
        if (isNucleotide) {
            writePos = assignGroup<Parameters::DBTYPE_NUCLEOTIDES, T, true, true>(
                positions, grouped, par.includeOnlyExtendable, par.covMode, par.covThr, NULL,
                par.weightThr, par.threads, threadOffsets, subMat, mask, ComputationPhase::Main,
                NULL);
            SORT_PARALLEL(grouped, grouped + writePos,
                          KmerPosition<T, false, true>::compareRepSequenceAndIdAndDiagReverse);
        } else {
            writePos = assignGroup<Parameters::DBTYPE_AMINO_ACIDS, T, true, true>(
                positions, grouped, par.includeOnlyExtendable, par.covMode, par.covThr, NULL,
                par.weightThr, par.threads, threadOffsets, subMat, mask, ComputationPhase::Main,
                NULL);
            SORT_PARALLEL(grouped, grouped + writePos,
                          KmerPosition<T, false, true>::compareRepSequenceAndIdAndDiag);
        }
        collectRoundEdges<T>(grouped, writePos, isNucleotide, edges);
    }
    delete[] positions;
    delete[] grouped;

    // Rounds overlap heavily by construction, so collapse the union here rather
    // than writing the same pair once per round. Keeping the highest score means
    // the round that agreed on the most k-mers wins.
    // The v2 rounds overlap heavily, so the same (pair, diagonal) is produced by
    // several of them. Sum the supporting k-mer counts rather than keeping one,
    // so the count the align stage accumulates globally stays meaningful.
    SORT_PARALLEL(edges.begin(), edges.end(), compareEdge);
    size_t unique = 0;
    for (size_t i = 0; i < edges.size(); i++) {
        if (unique > 0 && edges[i].getRep() == edges[unique - 1].getRep() &&
            edges[i].getMember() == edges[unique - 1].getMember() &&
            edges[i].diagonal == edges[unique - 1].diagonal &&
            edges[i].reverseStrand == edges[unique - 1].reverseStrand) {
            edges[unique - 1].score = static_cast<uint16_t>(
                std::min<int>(edges[unique - 1].score + edges[i].score, 65535));
            continue;
        }
        edges[unique] = edges[i];
        unique++;
    }
    edges.resize(unique);

    // Bucketed by representative key, which is what the align stage needs and
    // what brings every copy of a pair together (see CandidateEdge.h).
    for (size_t i = 0; i < edges.size(); i++) {
        writer.append(static_cast<unsigned int>(edges[i].getRep() / bucketSpan), edges[i]);
    }

    return edges.size();
}

}  // namespace

int kmerreduceparallel(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setLinearFilterDefault(&par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_CLUSTLINEAR);

    const std::string seqDb = par.db1;
    const std::string kmerDir = par.db2;
    const std::string edgeDir = par.db3;

    // Only the header is read, for the sequence type and the longest sequence.
    // The partitions themselves are self-contained.
    const int dbType = FileUtil::parseDbType(seqDb.c_str());
    const DenseIndex::Info info = DenseIndex::readInfo(seqDb);

    const std::string coordDir = kmerDir + "/coord";
    const std::string manifestPath = coordDir + "/shuffle.info";
    if (FileUtil::fileExists(manifestPath.c_str()) == false) {
        Debug(Debug::ERROR) << "No k-mer shuffle at " << kmerDir << " (" << manifestPath
                            << " is missing). Run kmermatcherparallel first.\n";
        EXIT(EXIT_FAILURE);
    }
    unsigned int partitionCount = 0;
    unsigned int kmerSize = 0;
    unsigned int waveCount = 1;
    {
        FILE *file = FileUtil::openFileOrDie(manifestPath.c_str(), "r", true);
        char name[64];
        size_t value;
        while (fscanf(file, "%63s\t%zu\n", name, &value) == 2) {
            const std::string key = name;
            if (key == "partitionCount") {
                partitionCount = static_cast<unsigned int>(value);
            } else if (key == "kmerSize") {
                kmerSize = static_cast<unsigned int>(value);
            } else if (key == "waveCount") {
                waveCount = static_cast<unsigned int>(value);
            }
        }
        fclose(file);
    }
    if (partitionCount == 0) {
        Debug(Debug::ERROR) << "Shuffle manifest " << manifestPath << " has no partition count\n";
        EXIT(EXIT_FAILURE);
    }
    // A wave's map wrote only its own slice of partition space, so its reduce
    // claims exactly that slice. The slicing must be identical to the map's, and
    // each wave needs its own queue: a shared one would record the whole
    // partition space done after wave 0 and skip every later wave outright.
    unsigned int waveFrom = 0;
    unsigned int waveTo = partitionCount;
    if (waveCount > 1) {
        if (par.kmerWave < 0) {
            Debug(Debug::ERROR) << "This shuffle was written in " << waveCount
                                << " waves. Reduce each one with --kmer-wave 0.."
                                << (waveCount - 1) << ", matching the wave its map wrote.\n";
            EXIT(EXIT_FAILURE);
        }
        if (static_cast<unsigned int>(par.kmerWave) >= waveCount) {
            Debug(Debug::ERROR) << "--kmer-wave " << par.kmerWave << " is out of range; this "
                                << "shuffle has " << waveCount << " waves\n";
            EXIT(EXIT_FAILURE);
        }
        // Exact: the sizing guarantees the wave count divides P.
        const unsigned int perWave = partitionCount / waveCount;
        waveFrom = static_cast<unsigned int>(par.kmerWave) * perWave;
        waveTo = waveFrom + perWave;
        Debug(Debug::INFO) << "Wave " << par.kmerWave << " of " << waveCount << ": partitions ["
                           << waveFrom << ", " << waveTo << ")\n";
    } else if (par.kmerWave > 0) {
        Debug(Debug::ERROR) << "--kmer-wave " << par.kmerWave << " given but this shuffle was "
                            << "written in one wave\n";
        EXIT(EXIT_FAILURE);
    }
    // The map decided k, so take it from the manifest rather than re-deriving it
    // here: the two must agree, and the map's value is the one on disk.
    par.kmerSize = static_cast<int>(kmerSize);
    par.printParameters(command.cmd, argc, argv, *command.params);
    Debug(Debug::INFO) << "Reducing " << partitionCount << " partitions of " << kmerDir << "\n";

    if (FileUtil::directoryExists(edgeDir.c_str()) == false) {
        FileUtil::makeDir(edgeDir.c_str());
    }
    const std::string reduceCoordDir = edgeDir + "/coord";
    if (FileUtil::directoryExists(reduceCoordDir.c_str()) == false) {
        FileUtil::makeDir(reduceCoordDir.c_str());
    }

    // Edge buckets are ranges of representative key. Sized so one bucket's
    // sequences are a comfortable slice for the align stage, which loads them
    // whole; the align workers' memory, not the reduce's, sets this.
    const uint64_t targetBucketBytes = Util::computeMemory(par.splitMemoryLimit) / 4;
    unsigned int bucketCount = 1;
    while (bucketCount < 65536 && info.dataSize / bucketCount > targetBucketBytes) {
        bucketCount *= 2;
    }
    uint64_t bucketSpan = (info.entryCount + bucketCount - 1) / bucketCount;

    const std::string edgeManifest = reduceCoordDir + "/edge.info";
    {
        FileLock lock(reduceCoordDir + "/edge.lock");
        lock.lock();
        if (FileUtil::fileExists(edgeManifest.c_str()) == false) {
            EdgeBucketWriter::createLayout(edgeDir, bucketCount);
            FILE *f = FileUtil::openAndDelete(edgeManifest.c_str(), "w");
            fprintf(f, "bucketCount\t%zu\n", (size_t)bucketCount);
            fprintf(f, "bucketSpan\t%zu\n", (size_t)bucketSpan);
            fprintf(f, "entryCount\t%zu\n", (size_t)info.entryCount);
            fclose(f);
        }
        lock.unlock();
    }
    // Re-read what the manifest actually says and refuse to disagree with it.
    //
    // Every worker derives bucketCount from Util::computeMemory(--split-memory-limit),
    // which with the workflow's default of 0 is *this node's* RAM. A heterogeneous
    // Slurm array, a restart on a different node, or a later wave would otherwise
    // route edges into a different bucketing than the run began with -- writing
    // into r<n> directories createLayout never made. The map already does exactly
    // this for shuffle.info.
    {
        unsigned int fileBucketCount = 0;
        uint64_t fileBucketSpan = 0;
        FILE *f = FileUtil::openFileOrDie(edgeManifest.c_str(), "r", true);
        char name[64];
        size_t value;
        while (fscanf(f, "%63s\t%zu\n", name, &value) == 2) {
            const std::string key = name;
            if (key == "bucketCount") fileBucketCount = static_cast<unsigned int>(value);
            else if (key == "bucketSpan") fileBucketSpan = value;
        }
        fclose(f);
        if (fileBucketCount != bucketCount || fileBucketSpan != bucketSpan) {
            Debug(Debug::ERROR) << "This worker derived " << bucketCount << " edge buckets of "
                                << bucketSpan << " keys, but " << edgeManifest << " says "
                                << fileBucketCount << " of " << fileBucketSpan
                                << ". The run was started with a different --split-memory-limit or "
                                << "on a node with different memory; edges would be routed into a "
                                << "different bucketing than the rest of the run. Re-run every "
                                << "worker of this stage with the same --split-memory-limit.\n";
            EXIT(EXIT_FAILURE);
        }
        bucketCount = fileBucketCount;
        bucketSpan = fileBucketSpan;
    }
    Debug(Debug::INFO) << "Writing edges into " << bucketCount << " representative-key buckets of "
                       << bucketSpan << " keys\n";

    SharedCounter workerCounter(reduceCoordDir + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();
    Debug(Debug::INFO) << "Worker " << workerId << " joined\n";

    BaseMatrix *subMat = createSubstitutionMatrix(par, dbType);

    EdgeBucketWriter *edgeWriter =
        new EdgeBucketWriter(edgeDir, bucketCount, "w" + SSTR(workerId));

    uint64_t edgeCount = 0;
    {
        WorkQueue queue(reduceCoordDir + "/reduce." + SSTR(par.kmerWave < 0 ? 0 : par.kmerWave) +
                            ".queue",
                        static_cast<int64_t>(waveTo - waveFrom));
        // One partition at a time per process: the sort and the greedy inside a
        // partition are already threaded, and a partition is sized to fill a node.
        const bool finished = queue.drain(workerId, [&](size_t item) {
            const size_t partition = waveFrom + item;
            // Stamps every block this partition writes with (partition, worker).
            // If this worker dies before the queue records the item done, another
            // redoes it and the align stage drops these blocks in favour of the
            // redo's -- without which both copies would be summed.
            edgeWriter->beginPartition(static_cast<unsigned int>(partition), workerId);
            if (info.maxSeqLen < SHRT_MAX) {
                edgeCount += reducePartition<short>(kmerDir, static_cast<unsigned int>(partition),
                                                    dbType, par, subMat, *edgeWriter, bucketSpan);
            } else {
                edgeCount += reducePartition<int>(kmerDir, static_cast<unsigned int>(partition), dbType,
                                                  par, subMat, *edgeWriter, bucketSpan);
            }
            // Before drain() records the bucket done, so a worker that dies never
            // leaves an item complete whose edges were still buffered.
            edgeWriter->flushAll();
        });
        if (finished == false) {
            Debug(Debug::ERROR) << "Reduce stage stalled: work remains but no partition is claimable\n";
            EXIT(EXIT_FAILURE);
        }
    }

    {
        // The k-mer buckets are dead once their partitions have been reduced --
        // nothing downstream reads them again -- so they go here, at every wave
        // count. Gating this on waveCount > 1, as it used to be, meant a
        // single-wave run kept the whole shuffle on disk for the rest of the job,
        // including all of pass 2: measured at 100M, kmer1 (49.8 GB) and kmer2
        // (52.9 GB) were both resident at peak, and the shuffle is the largest
        // intermediate the pipeline writes.
        //
        // Safe here because drain() returns only once every partition of the wave
        // is recorded done, so no worker holding a live lease is still reading
        // one. Every worker that observes the queue finished reaches this and
        // unlinks, hence the ENOENT tolerance below.
        //
        // A worker whose lease lapsed while still inside reducePartition can be
        // reading a shard as another unlinks it. That costs nothing: its item was
        // redone and recorded by someone else, so its edges are discarded by block
        // header anyway. What it must not do is turn into a stage failure, so
        // readPartitionAsPositions treats a shard that disappears under it as
        // empty rather than as an error.
        for (unsigned int p = waveFrom; p < waveTo; p++) {
            const std::vector<std::string> shards = KmerBucketReader::shardFiles(kmerDir, p);
            for (size_t i = 0; i < shards.size(); i++) {
                if (unlink(shards[i].c_str()) != 0 && errno != ENOENT) {
                    Debug(Debug::ERROR) << "Cannot remove reduced bucket " << shards[i] << ": "
                                        << strerror(errno) << "\n";
                    EXIT(EXIT_FAILURE);
                }
            }
        }
        Debug(Debug::INFO) << "Removed the consumed k-mer buckets"
                           << (waveCount > 1 ? " of wave " + SSTR(par.kmerWave) : "") << "\n";
    }

    Debug(Debug::INFO) << "Worker " << workerId << " wrote " << edgeCount
                       << " candidate edges\n";
    edgeWriter->close();
    delete edgeWriter;
    delete subMat;

    return EXIT_SUCCESS;
}
