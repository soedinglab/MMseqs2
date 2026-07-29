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

#include <climits>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

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
        FILE *file = FileUtil::openFileOrDie(shards[i].c_str(), "rb", true);
        while (true) {
            const size_t got = fread(block.data(), sizeof(KmerRecord), blockRecords, file);
            if (got == 0) {
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


// Collapses one grouping round's output into (representative, member) edges.
//
// The same pair is produced once for every k-mer the two sequences share, so
// emitting them raw duplicates heavily -- measured 4.1 M records for 1.1 M
// distinct pairs on 1M sequences. Stock collapses them the same way in
// writeKmerMatcherResult (kmermatcher.cpp:1789): keep the diagonal that the most
// k-mers agree on, and use that count as the prefilter score. The array is sorted
// by (rep, member, diagonal), so both runs are contiguous.
template <typename T>
void collectRoundEdges(KmerPosition<T, false, true> *grouped, size_t writePos, bool isNucleotide,
                       std::vector<CandidateEdge> &edges) {
    CandidateEdge edge;
    size_t i = 0;
    while (i < writePos && grouped[i].kmer != SIZE_MAX) {
        const size_t repRaw = grouped[i].kmer;
        const DBKeyType member = grouped[i].id;
        T bestDiagonal = grouped[i].pos;
        size_t bestCount = 0;
        T runDiagonal = grouped[i].pos;
        size_t runCount = 0;
        size_t j = i;
        while (j < writePos && grouped[j].kmer == repRaw && grouped[j].id == member) {
            if (grouped[j].pos == runDiagonal) {
                runCount++;
            } else {
                runDiagonal = grouped[j].pos;
                runCount = 1;
            }
            if (runCount > bestCount) {
                bestCount = runCount;
                bestDiagonal = runDiagonal;
            }
            j++;
        }
        i = j;

        size_t rep = repRaw;
        // Stock keeps the strand in bit 63 of the representative key; a 48-bit key
        // has no room for it, so it moves into its own field.
        edge.reverseStrand = 0;
        if (isNucleotide) {
            edge.reverseStrand = BIT_CHECK(rep, 63) == false;
            rep = BIT_CLEAR(rep, 63);
        }
        // A sequence being its own representative carries no information here:
        // keys are dense, so any key that never appears as a member is its own
        // representative by construction, and the final greedy derives singletons
        // from that rather than from a marker edge.
        if (static_cast<DBKeyType>(rep) == member) {
            continue;
        }
        edge.setRep(static_cast<uint64_t>(rep));
        edge.setMember(static_cast<uint64_t>(member));
        edge.diagonal = static_cast<int16_t>(bestDiagonal);
        edge.score = static_cast<uint8_t>(std::min<size_t>(bestCount, 255));
        edges.push_back(edge);
    }
}

// Orders edges so duplicates of a pair are adjacent, best score first.
bool compareEdge(const CandidateEdge &a, const CandidateEdge &b) {
    const uint64_t ra = a.getRep(), rb = b.getRep();
    if (ra != rb) return ra < rb;
    const uint64_t ma = a.getMember(), mb = b.getMember();
    if (ma != mb) return ma < mb;
    return a.score > b.score;
}

// Groups one partition and writes its edges.
template <typename T>
uint64_t reducePartition(const std::string &kmerDir, const std::string &edgeDir,
                         unsigned int partition, int dbType, Parameters &par, BaseMatrix *subMat) {
    const bool isNucleotide = Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);
    EdgeWriter writer(EdgeWriter::partitionPath(edgeDir, partition));

    const uint64_t recordCount = KmerBucketReader::countRecords(kmerDir, partition);
    if (recordCount == 0) {
        writer.close();
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
    SORT_PARALLEL(edges.begin(), edges.end(), compareEdge);
    for (size_t i = 0; i < edges.size(); i++) {
        if (i > 0 && edges[i].getRep() == edges[i - 1].getRep() &&
            edges[i].getMember() == edges[i - 1].getMember()) {
            continue;
        }
        writer.append(edges[i]);
    }

    writer.close();
    return writer.getEdgeCount();
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
            }
        }
        fclose(file);
    }
    if (partitionCount == 0) {
        Debug(Debug::ERROR) << "Shuffle manifest " << manifestPath << " has no partition count\n";
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

    SharedCounter workerCounter(reduceCoordDir + "/worker.counter");
    const int64_t workerId = workerCounter.fetchAdd();
    Debug(Debug::INFO) << "Worker " << workerId << " joined\n";

    BaseMatrix *subMat = createSubstitutionMatrix(par, dbType);
    uint64_t edgeCount = 0;
    {
        WorkQueue queue(reduceCoordDir + "/reduce.queue", static_cast<int64_t>(partitionCount));
        // One partition at a time per process: the sort and the greedy inside a
        // partition are already threaded, and a partition is sized to fill a node.
        const bool finished = queue.drain(workerId, [&](size_t partition) {
            if (info.maxSeqLen < SHRT_MAX) {
                edgeCount += reducePartition<short>(kmerDir, edgeDir,
                                                    static_cast<unsigned int>(partition), dbType,
                                                    par, subMat);
            } else {
                edgeCount += reducePartition<int>(kmerDir, edgeDir,
                                                  static_cast<unsigned int>(partition), dbType, par,
                                                  subMat);
            }
        });
        if (finished == false) {
            Debug(Debug::ERROR) << "Reduce stage stalled: work remains but no partition is claimable\n";
            EXIT(EXIT_FAILURE);
        }
    }

    Debug(Debug::INFO) << "Worker " << workerId << " wrote " << edgeCount << " candidate edges\n";
    delete subMat;

    return EXIT_SUCCESS;
}
