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
#include "LengthRankTable.h"
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
                                KmerPosition<T, true, true> *out, size_t capacity,
                                const LengthRankTable *lengths, unsigned int sliceCount,
                                unsigned int slice) {
    const std::vector<std::string> shards = KmerBucketReader::shardFiles(kmerDir, partition);
    std::vector<KmerRecord> block;
    size_t filled = 0;
    for (size_t i = 0; i < shards.size(); i++) {
        // A worker whose lease lapsed can still be here while the workers that
        // finished the wave unlink its shards. KmerShardReader treats a shard that
        // will not open as empty, which is the race this stage should survive
        // rather than a reason to fail.
        KmerShardReader reader(shards[i]);
        while (true) {
            block.clear();
            if (reader.next(block, lengths) == false) {
                break;
            }
            // Bound before filtering: a slice keeps at most the whole block, and
            // checking after would already have written past the end.
            if (sliceCount <= 1 && filled + block.size() > capacity) {
                Debug(Debug::ERROR) << "Partition " << partition << " holds more k-mers than the "
                                    << "shard sizes reported. A shard was written while this was "
                                    << "reading it.\n";
                EXIT(EXIT_FAILURE);
            }
            for (size_t r = 0; r < block.size(); r++) {
                if (sliceCount > 1 && kmerSliceOf(block[r].kmer, sliceCount) != slice) {
                    continue;
                }
                kmerRecordToPosition(block[r], out[filled]);
                filled++;
            }
        }
        if (reader.endedTorn()) {
            // Stopped at the last whole block rather than failing. A torn tail is
            // what an interrupted worker leaves; its item is redone into a
            // *different* shard, so nothing is lost. Making it fatal meant every
            // later reduce of that partition died on it forever, with no recovery
            // but deleting the file by hand.
            Debug(Debug::WARNING) << "K-mer bucket " << shards[i]
                                  << " ends on a partial block, as an interrupted worker leaves "
                                     "it; using the whole blocks it holds.\n";
        }
    }
    return filled;
}

// Never slice further than this. Past it the passes over the partition cost more
// than the grouping, and a budget this far below one partition is a configuration
// error worth surfacing rather than working around.
const unsigned int MAX_REDUCE_SLICES = 1024;

// Records per slice, so each slice's array is allocated for what it actually
// holds rather than for an assumed even split.
//
// One extra decode pass over the partition. Against the sliceCount passes the
// slicing already costs it is nothing, and the alternative -- allocating
// recordCount / sliceCount with a safety factor -- either wastes memory in the
// case where memory is already short or overflows and has to start again.
void countPartitionSlices(const std::string &kmerDir, unsigned int partition,
                          const LengthRankTable *lengths, unsigned int sliceCount,
                          std::vector<uint64_t> &counts) {
    counts.assign(sliceCount, 0);
    const std::vector<std::string> shards = KmerBucketReader::shardFiles(kmerDir, partition);
    std::vector<KmerRecord> block;
    for (size_t i = 0; i < shards.size(); i++) {
        KmerShardReader reader(shards[i]);
        while (true) {
            block.clear();
            if (reader.next(block, lengths) == false) {
                break;
            }
            for (size_t r = 0; r < block.size(); r++) {
                counts[kmerSliceOf(block[r].kmer, sliceCount)]++;
            }
        }
    }
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



// Groups one slice of one partition and writes its edges. sliceCount == 1 is the
// whole partition, which is the normal case.
template <typename T>
uint64_t reduceSlice(const std::string &kmerDir, unsigned int partition, int dbType,
                     Parameters &par, BaseMatrix *subMat, EdgeBucketWriter &writer,
                     uint64_t bucketSpan, const LengthRankTable *lengths, unsigned int sliceCount,
                     unsigned int slice, uint64_t sliceRecords) {
    const bool isNucleotide = Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_NUCLEOTIDES);

    const uint64_t recordCount = sliceRecords;
    if (recordCount == 0) {
        return 0;
    }

    // One slot past the end for the SIZE_T_MAX sentinel assignGroup stops on,
    // matching initKmerPositionMemory (kmermatcher.cpp:43).
    KmerPosition<T, true, true> *positions =
        new (std::nothrow) KmerPosition<T, true, true>[recordCount + 1];
    Util::checkAllocation(positions, "Cannot allocate the k-mer partition");

    size_t count = readPartitionAsPositions<T>(kmerDir, partition, positions, recordCount, lengths,
                                               sliceCount, slice);
    if (isNucleotide) {
        SORT_PARALLEL(positions, positions + count,
                      KmerPosition<T, true, true>::compareRepSequenceAndIdAndPosReverse);
    } else {
        SORT_PARALLEL(positions, positions + count,
                      KmerPosition<T, true, true>::compareRepSequenceAndIdAndPos);
    }
    count = dropDuplicates<T>(positions, count);
    memset(&positions[count], 0xFF, sizeof(positions[0]));

    // The largest k-mer group in this partition.
    //
    // This is the quantity that separates the two kinds of skew, and they have
    // different fixes. If partitions are uneven but every group is small, the
    // partition can be cut into k-mer-hash sub-slices and each grouped
    // independently -- exact, because every occurrence of a k-mer shares a hash
    // and so stays together. If one *group* is itself too large, no sub-slicing
    // helps: the group is the atomic unit, since assignGroup picks a centre and
    // emits pairs within one group's index range and buildThreadOffsets refuses
    // to split a group across threads.
    //
    // Measured here rather than assumed, because the two need very different
    // amounts of work and nobody has reported which one real data produces.
    {
        size_t largestGroup = 0;
        size_t groupCount = 0;
        size_t groupStart = 0;
        for (size_t i = 1; i <= count; i++) {
            const bool boundary = (i == count) || (positions[i].kmer != positions[groupStart].kmer);
            if (boundary) {
                largestGroup = std::max(largestGroup, i - groupStart);
                groupCount++;
                groupStart = i;
            }
        }
        Debug(Debug::INFO) << "Partition " << partition
                           << (sliceCount > 1 ? (" slice " + SSTR(slice) + "/" + SSTR(sliceCount))
                                              : std::string())
                           << ": " << count << " k-mers, "
                           << groupCount << " groups, largest group " << largestGroup << " ("
                           << (count > 0 ? (100.0 * static_cast<double>(largestGroup) /
                                            static_cast<double>(count))
                                         : 0.0)
                           << "% of the partition)\n";
    }

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

// Groups one partition, slicing it if it will not fit the worker's budget.
template <typename T>
uint64_t reducePartition(const std::string &kmerDir, unsigned int partition, int dbType,
                         Parameters &par, BaseMatrix *subMat, EdgeBucketWriter &writer,
                         uint64_t bucketSpan, const LengthRankTable *lengths,
                         uint64_t expectedRecords, uint64_t memoryLimit, uint64_t *recordsOut) {
    const uint64_t recordCount = KmerBucketReader::countRecords(kmerDir, partition);
    // Handed back so the caller can keep a running mean without re-scanning every
    // block header of the partition a second time.
    if (recordsOut != NULL) {
        *recordsOut = recordCount;
    }
    if (recordCount == 0) {
        return 0;
    }

    // The two arrays the grouping holds: the k-mers it reads and the pairs
    // assignGroup writes. The edges accumulate alongside them, which is why the
    // arrays are sized against a fraction of the budget rather than all of it --
    // the same split deriveKmerShuffleSizing's factor of 3 encodes.
    const size_t bytesPerRecord =
        sizeof(KmerPosition<T, true, true>) + sizeof(KmerPosition<T, false, true>);
    const uint64_t residentBytes = recordCount * bytesPerRecord;

    // P is derived from the *average* partition, so a partition well above it is
    // not a sizing mistake -- it is skew, and the two need different responses.
    // Saying so matters because raising --split-memory-limit cannot help: the
    // partition of a k-mer is a pure function of the k-mer, so a k-mer carrying an
    // outsized share of the database stays in one partition however large P is.
    if (expectedRecords > 0 && recordCount > 2 * expectedRecords) {
        Debug(Debug::WARNING)
            << "Partition " << partition << " holds " << recordCount << " k-mers, "
            << (static_cast<double>(recordCount) / static_cast<double>(expectedRecords))
            << "x the average of " << expectedRecords << " over the partitions this worker has "
            << "already grouped. This is k-mer skew, not a partition count that is too low.\n";
    }

    // Sliced only when the two arrays alone exceed the whole budget.
    //
    // Not a fraction of it, deliberately. deriveKmerShuffleSizing already aims a
    // partition at roughly limit/3 on-disk bytes, which is 0.64 * limit once the
    // 24-byte records become 46 bytes of KmerPosition -- so a threshold of, say,
    // 0.6 * limit would slice every partition of every ordinary run and double
    // its passes for nothing. At 1.0 the trigger is a partition ~56% above what
    // the sizing intended, which is the case slicing exists for.
    unsigned int sliceCount = 1;
    if (par.reduceSlices > 0) {
        sliceCount = static_cast<unsigned int>(par.reduceSlices);
    } else if (memoryLimit > 0) {
        while (sliceCount < MAX_REDUCE_SLICES && residentBytes / sliceCount > memoryLimit) {
            sliceCount *= 2;
        }
    }
    if (sliceCount == 1) {
        return reduceSlice<T>(kmerDir, partition, dbType, par, subMat, writer, bucketSpan, lengths,
                              1, 0, recordCount);
    }

    Debug(Debug::INFO) << "Partition " << partition << " needs about "
                       << (residentBytes / (1024 * 1024)) << " MB to group whole, past this "
                       << "worker's budget; grouping it in " << sliceCount
                       << " k-mer slices instead. Every occurrence of a k-mer shares a slice, so "
                       << "the edges are the same ones grouping it whole would produce.\n";

    std::vector<uint64_t> sliceRecords;
    countPartitionSlices(kmerDir, partition, lengths, sliceCount, sliceRecords);

    uint64_t largest = 0;
    for (size_t i = 0; i < sliceRecords.size(); i++) {
        largest = std::max(largest, sliceRecords[i]);
    }
    if (memoryLimit > 0 && largest * bytesPerRecord > memoryLimit) {
        // Only reachable when one k-mer *group* is itself too big: slices split
        // groups apart, never within one, so no slice count can shrink it. Said
        // plainly, because the obvious responses -- more slices, more partitions,
        // a bigger limit -- are all useless against it.
        Debug(Debug::WARNING)
            << "Slice " << sliceCount << "-way still leaves " << largest << " k-mers ("
            << (largest * bytesPerRecord / (1024 * 1024)) << " MB) in one slice. A single k-mer "
            << "group cannot be divided without changing the answer, so this is the floor for "
            << "this partition.\n";
    }

    uint64_t edges = 0;
    for (unsigned int slice = 0; slice < sliceCount; slice++) {
        edges += reduceSlice<T>(kmerDir, partition, dbType, par, subMat, writer, bucketSpan,
                                lengths, sliceCount, slice, sliceRecords[slice]);
        // Flushed between slices, so no block ever holds two of them.
        //
        // A slice's edges are sorted by representative, but the next slice starts
        // over at low keys, and the writer buffers across calls -- so without this
        // a bucket's buffer would hold a descending step in the middle and the
        // packed encoder's delta would wrap. The encoder checks for exactly that
        // and refused the run, which is how this was found rather than shipped.
        writer.flushAll();
    }
    return edges;
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

    // The k-mer record does not carry seqLen; the key does, because keys are
    // length-ranked. Opened once per process -- it is a few tens of kilobytes and
    // every decoded record consults it.
    LengthRankTable lengths;
    lengths.open(seqDb);
    if (lengths.getEntryCount() != info.entryCount) {
        Debug(Debug::ERROR) << "The length-rank table describes " << lengths.getEntryCount()
                            << " sequences but the database has " << info.entryCount
                            << ". They are from different builds.\n";
        EXIT(EXIT_FAILURE);
    }

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
    unsigned int alphabetSizeAa = 0;
    unsigned int alphabetSizeNucl = 0;
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
            } else if (key == "alphabetSizeAa") {
                alphabetSizeAa = static_cast<unsigned int>(value);
            } else if (key == "alphabetSizeNucl") {
                alphabetSizeNucl = static_cast<unsigned int>(value);
            }
        }
        fclose(file);
    }
    if (partitionCount == 0) {
        Debug(Debug::ERROR) << "Shuffle manifest " << manifestPath << " has no partition count\n";
        EXIT(EXIT_FAILURE);
    }
    if (alphabetSizeAa == 0) {
        Debug(Debug::ERROR) << "Shuffle manifest " << manifestPath << " has no alphabet size. It "
                            << "was written by an older build, whose reduce scored the k-mers with "
                            << "a different alphabet than the map encoded them in. Re-run the map "
                            << "for this shuffle with this build.\n";
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
    // The map decided k *and* the alphabet, so both come from the manifest rather
    // than being re-derived or re-passed here. setKmerLengthAndAlphabet picks them
    // together from --min-seq-id (kmermatcher.cpp:2098-2110), and this command
    // cannot run it -- it has no sequence database to size it from. Taking only k
    // and leaving the alphabet at setLinearFilterDefault's 13 meant that at
    // --min-seq-id >= 0.99 the map encoded adjacent residues with the 21-letter
    // aa2num while the reduce scored them through a ReducedMatrix(13) whose rows
    // 13..20 are zero. The alphabet is a property of the k-mers on disk, so the
    // manifest is where it belongs.
    par.kmerSize = static_cast<int>(kmerSize);
    par.alphabetSize.values.aminoacid(static_cast<int>(alphabetSizeAa));
    if (alphabetSizeNucl > 0) {
        par.alphabetSize.values.nucleotide(static_cast<int>(alphabetSizeNucl));
    }
    par.printParameters(command.cmd, argc, argv, *command.params);
    Debug(Debug::INFO) << "Reducing " << partitionCount << " partitions of " << kmerDir << "\n";

    if (FileUtil::directoryExists(edgeDir.c_str()) == false) {
        FileUtil::makeDir(edgeDir.c_str());
    }
    const std::string reduceCoordDir = edgeDir + "/coord";
    if (FileUtil::directoryExists(reduceCoordDir.c_str()) == false) {
        FileUtil::makeDir(reduceCoordDir.c_str());
    }

    // Edge buckets are ranges of representative key, and they are two things at
    // once:
    //
    //   1. a slice that has to fit an align worker, which loads its sequences
    //      whole -- the align workers' memory, not the reduce's, sets this;
    //   2. the unit of work the align stage's queue hands out, which makes the
    //      count a hard ceiling on how many workers that stage can use at all.
    //
    // Sizing it from (1) alone -- which is what this did -- means a *generous*
    // --split-memory-limit yields *fewer* buckets and less parallelism. That is
    // backwards, and silent: the starved workers find nothing claimable, sleep in
    // the heartbeat and exit 0. Measured at 1e9 with 800G per node, the edges are
    // 125 GB against a 200 GB target, so bucketCount stayed 1 and alignparallel
    // ran 3h16m on a single node while three others idled -- half the wall clock
    // of the whole run.
    //
    // So the memory figure becomes a ceiling and the target is taken well below
    // it. The cost is bounded and already handled: EdgeBucketWriter caps its open
    // descriptors (deriveEdgeDescriptorBudget) and floors the per-bucket write
    // buffer at 64 records, so a few thousand buckets is a few MB of buffering
    // rather than descriptor exhaustion.
    // Derived in CandidateEdge so it can be tested without a database; see the
    // comment on deriveAlignBucketCount for why the memory limit is a ceiling
    // rather than the target.
    unsigned int bucketCount = deriveAlignBucketCount(
        info.dataSize, info.entryCount, Util::computeMemory(par.splitMemoryLimit),
        static_cast<unsigned int>(std::max(par.workerCount, 0)));
    uint64_t bucketSpan = (info.entryCount + bucketCount - 1) / bucketCount;

    // The manifest is authoritative once it exists; the derivation above only
    // creates it. This is the same rule the map applies to shuffle.info, and it is
    // load-bearing for the same reason: bucketCount comes from
    // Util::computeMemory(--split-memory-limit), which with the workflow default
    // of 0 is *this node's* physical RAM (and ignores cgroup limits besides).
    //
    // Refusing to run when the two disagree -- which is what this used to do --
    // meant one differently-sized node in a Slurm array, a restart on another
    // partition, or a container memory limit killed the whole stage. Adopting the
    // recorded layout instead is both safer and correct: the bucketing has to be
    // fixed for the life of the edge directory, and the run that created it is the
    // one that decided.
    const std::string edgeManifest = reduceCoordDir + "/edge.info";
    {
        FileLock lock(reduceCoordDir + "/edge.lock");
        lock.lock();
        if (FileUtil::fileExists(edgeManifest.c_str()) == false) {
            EdgeBucketWriter::createLayout(edgeDir, bucketCount);
            // Private temporary then rename, so a worker killed mid-write cannot
            // leave an empty manifest that every later worker reads as
            // authoritative and no restart can repair.
            const std::string tmp = edgeManifest + ".tmp." + SSTR(getpid());
            FILE *f = FileUtil::openAndDelete(tmp.c_str(), "w");
            fprintf(f, "bucketCount\t%zu\n", (size_t)bucketCount);
            fprintf(f, "bucketSpan\t%zu\n", (size_t)bucketSpan);
            fprintf(f, "entryCount\t%zu\n", (size_t)info.entryCount);
            // The shape of the reduce itself, so the align stage can tell "every
            // partition has been reduced" from "the queues I happened to find were
            // complete". Without it, starting the align before the last wave
            // existed produced a short authority vector that looked finished.
            fprintf(f, "partitionCount\t%zu\n", (size_t)partitionCount);
            fprintf(f, "waveCount\t%zu\n", (size_t)waveCount);
            if (fclose(f) != 0) {
                Debug(Debug::ERROR) << "Cannot close " << tmp << ": " << strerror(errno) << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (rename(tmp.c_str(), edgeManifest.c_str()) != 0) {
                Debug(Debug::ERROR) << "Cannot publish " << edgeManifest << ": " << strerror(errno)
                                    << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        lock.unlock();
    }
    {
        unsigned int fileBucketCount = 0;
        uint64_t fileBucketSpan = 0;
        uint64_t fileEntryCount = 0;
        FILE *f = FileUtil::openFileOrDie(edgeManifest.c_str(), "r", true);
        char name[64];
        size_t value;
        while (fscanf(f, "%63s\t%zu\n", name, &value) == 2) {
            const std::string key = name;
            if (key == "bucketCount") fileBucketCount = static_cast<unsigned int>(value);
            else if (key == "bucketSpan") fileBucketSpan = value;
            else if (key == "entryCount") fileEntryCount = value;
        }
        fclose(f);
        if (fileBucketCount == 0 || fileBucketSpan == 0) {
            Debug(Debug::ERROR) << "Edge manifest " << edgeManifest << " has no bucket layout\n";
            EXIT(EXIT_FAILURE);
        }
        // entryCount still *is* checked, because unlike the bucketing it describes
        // the input rather than the node: a different database means these workers
        // are not running the same job at all, and the key ranges would not line up.
        if (fileEntryCount != info.entryCount) {
            Debug(Debug::ERROR) << "The edge buckets at " << edgeDir << " were laid out for "
                                << fileEntryCount << " sequences, but this worker was given "
                                << info.entryCount << ". Every worker must run on the same "
                                << "database.\n";
            EXIT(EXIT_FAILURE);
        }
        if (fileBucketCount != bucketCount || fileBucketSpan != bucketSpan) {
            Debug(Debug::WARNING) << "This worker would have derived " << bucketCount
                                  << " edge buckets of " << bucketSpan << " keys from its own "
                                  << "memory, but " << edgeManifest << " says " << fileBucketCount
                                  << " of " << fileBucketSpan
                                  << "; using the recorded layout. This is expected on a "
                                  << "heterogeneous allocation or under a cgroup memory limit.\n";
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
        new EdgeBucketWriter(edgeDir, bucketCount, "w" + SSTR(workerId),
                             256 * 1024 * 1024, par.rawRecords != 0);

    uint64_t edgeCount = 0;
    {
        WorkQueue queue(reduceCoordDir + "/reduce." + SSTR(par.kmerWave < 0 ? 0 : par.kmerWave) +
                            ".queue",
                        static_cast<int64_t>(waveTo - waveFrom));
    // Say so when this worker cannot possibly get work.
    //
    // The stage sizes its work units from --split-memory-limit and --workers, not
    // from how many workers actually turn up, so an allocation larger than the
    // unit count leaves the surplus with nothing to claim -- they sleep in the
    // heartbeat and exit 0, and the stage looks like it succeeded on every node
    // while running on one. That is exactly how this went unnoticed until a 1e9
    // run was profiled. It cannot be made an error, because a worker joining a
    // nearly finished stage legitimately finds nothing left; but it can be said
    // out loud when the worker never had a chance.
    if (workerId >= queue.getItemCount()) {
        Debug(Debug::WARNING)
            << "Worker " << workerId << " has no work: this stage was sized for "
            << queue.getItemCount() << " work unit(s) and cannot use more workers than that. "
            << "Pass --workers <n> so the unit count is sized for the allocation.\n";
    }

        // One partition at a time per process: the sort and the greedy inside a
        // partition are already threaded, and a partition is sized to fill a node.
        // A running mean over the partitions this worker has already grouped.
        //
        // Compared against its peers rather than against a figure derived from the
        // manifest: the reduce does not know --kmer-per-seq (the workflow does not
        // pass it here), and "unusual next to the other partitions of this run" is
        // the comparison that actually identifies skew anyway.
        uint64_t seenPartitions = 0;
        uint64_t seenRecords = 0;
        const bool finished = queue.drain(workerId, [&](size_t item) {
            const size_t partition = waveFrom + item;
            const uint64_t expectedRecords =
                (seenPartitions > 0) ? (seenRecords / seenPartitions) : 0;
            // Util::computeMemory, as every other reader of this parameter does:
            // --split-memory-limit is a byte-parsed string whose 0 means "derive
            // from the machine", and the workflow leaves it at 0 by default. Taken
            // raw, the limit is 0, the slicing below never triggers however skewed
            // a partition is, and the stage has no way to fit one that does not
            // fit -- every worker in turn allocates, dies, and hands the item to
            // the next.
            const uint64_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);
            uint64_t partitionRecords = 0;
            // Stamps every block this partition writes with (partition, worker).
            // If this worker dies before the queue records the item done, another
            // redoes it and the align stage drops these blocks in favour of the
            // redo's -- without which both copies would be summed.
            edgeWriter->beginPartition(static_cast<unsigned int>(partition), workerId);
            if (info.maxSeqLen < SHRT_MAX) {
                edgeCount += reducePartition<short>(kmerDir, static_cast<unsigned int>(partition),
                                                    dbType, par, subMat, *edgeWriter, bucketSpan,
                                                    &lengths, expectedRecords, memoryLimit,
                                                    &partitionRecords);
            } else {
                edgeCount += reducePartition<int>(kmerDir, static_cast<unsigned int>(partition), dbType,
                                                  par, subMat, *edgeWriter, bucketSpan, &lengths,
                                                  expectedRecords, memoryLimit, &partitionRecords);
            }
            seenRecords += partitionRecords;
            seenPartitions++;
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
        // is recorded done, so no worker holding a live lease is still reading one
        // -- and, just as importantly, so that a worker cannot delete an item's
        // input before the item is recorded complete. Deleting inside the work
        // item would be cheaper still, but a worker that unlinked and then died
        // before completing would leave the redo to read an empty partition, emit
        // no edges, and record *that* as success.
        //
        // Done by exactly one worker, under a lock and behind a sentinel. Every
        // worker used to walk every partition of the wave and unlink every shard:
        // at 1e12 that is 8192 workers x 1024 partitions x 8192 shards, roughly
        // 7e10 unlink syscalls of which all but 1/8192 return ENOENT, plus 8.4e6
        // opendir calls per worker. The others block only for as long as the one
        // sweep takes, which is work they were all doing anyway.
        //
        // The sentinel rather than a counter, so an interrupted sweep is redone by
        // the next run of the stage instead of being remembered as finished.
        //
        // A worker whose lease lapsed while still inside reducePartition can be
        // reading a shard as the sweeper unlinks it. That costs nothing: its item
        // was redone and recorded by someone else, so its edges are discarded by
        // block header anyway. What it must not do is turn into a stage failure, so
        // readPartitionAsPositions treats a shard that disappears under it as
        // empty rather than as an error.
        const std::string waveTag = SSTR(par.kmerWave < 0 ? 0 : par.kmerWave);
        const std::string cleanupDone = reduceCoordDir + "/cleanup." + waveTag + ".done";
        FileLock cleanupLock(reduceCoordDir + "/cleanup." + waveTag + ".lock");
        cleanupLock.lock();
        if (FileUtil::fileExists(cleanupDone.c_str()) == false) {
            for (unsigned int p = waveFrom; p < waveTo; p++) {
                const std::vector<std::string> shards = KmerBucketReader::shardFiles(kmerDir, p);
                for (size_t i = 0; i < shards.size(); i++) {
                    if (unlink(shards[i].c_str()) != 0 && errno != ENOENT) {
                        Debug(Debug::ERROR) << "Cannot remove reduced bucket " << shards[i] << ": "
                                            << strerror(errno) << "\n";
                        cleanupLock.unlock();
                        EXIT(EXIT_FAILURE);
                    }
                }
            }
            FILE *sentinel = FileUtil::openAndDelete(cleanupDone.c_str(), "w");
            if (fclose(sentinel) != 0) {
                Debug(Debug::ERROR) << "Cannot close " << cleanupDone << "\n";
                cleanupLock.unlock();
                EXIT(EXIT_FAILURE);
            }
            Debug(Debug::INFO) << "Removed the consumed k-mer buckets"
                               << (waveCount > 1 ? " of wave " + SSTR(par.kmerWave) : "") << "\n";
        }
        cleanupLock.unlock();
    }

    Debug(Debug::INFO) << "Worker " << workerId << " wrote " << edgeCount
                       << " candidate edges\n";
    edgeWriter->close();
    delete edgeWriter;
    delete subMat;

    return EXIT_SUCCESS;
}
