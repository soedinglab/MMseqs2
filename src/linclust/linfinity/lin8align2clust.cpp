#include "Lin8Db.h"
#include "Lin8DbReader.h"
#include "Parameters.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include "Timer.h"
#include "NodePlacement.h"
#include "Alignment.h"
#include "Matcher.h"

#include <sstream>
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "EvalueComputation.h"
#include "BlockAligner.h"
#include "FastSort.h"

#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

static const size_t STREAM_ROWS = 1u << 16;
// the lanes carry slice reads, so a lane arena only has to hold the longest sequence
static const size_t ARENA_BYTES = 1u << 20;

static const size_t MEMBERS_PER_ALIGN_BATCH = 1024;
// a batch is read in slices sorted by file position, the next one in flight while this one aligns
static const size_t READ_SLICES = 8;
static const size_t SLICE_MIN_ROWS = 1u << 16;

struct MemberBatch {
    size_t group;
    size_t from;
    size_t count;
};

static void readPipelineShape(const std::string &path, unsigned int &nodes, size_t &repRankBlocks, uint64_t &ranks) {
    FILE *in = fopen(path.c_str(), "r");
    if (in == NULL) {
        Debug(Debug::ERROR) << "Cannot open " << path << ". Run lin8pref first\n";
        EXIT(EXIT_FAILURE);
    }
    nodes = 0;
    repRankBlocks = 0;
    size_t seen = 0;
    const bool read = fscanf(in, "nodes\t%u\nrepRankBlocks\t%zu\nranks\t%zu", &nodes, &repRankBlocks, &seen) == 3;
    fclose(in);
    ranks = seen;
    if (read == false || nodes == 0 || repRankBlocks == 0 || repRankBlocks > PairRecord::MAX_REP_RANK_BLOCKS) {
        Debug(Debug::ERROR) << path << " names " << nodes << " nodes and " << repRankBlocks
                            << " repRankBlocks, which is not a set of prefilter pairs\n";
        EXIT(EXIT_FAILURE);
    }
}

class RepRankBlockReader {
public:
    RepRankBlockReader(const std::string &prefix, unsigned int prefNodes, size_t repRankBlock, size_t bufferRows,
                uint64_t skipRows)
        : buffer(bufferRows), at(0), filled(0), stopped(false) {
        const std::string path =
            prefix + "." + SSTR(repRankBlock % prefNodes) + "." + SSTR(repRankBlock);
        file = fopen(path.c_str(), "r");
        if (file == NULL) {
            Debug(Debug::ERROR) << "Cannot open " << path << ", which node " << (repRankBlock % prefNodes)
                                << " should have written for repRankBlock " << repRankBlock << "\n";
            EXIT(EXIT_FAILURE);
        }
        const off_t skip = static_cast<off_t>(skipRows * PairRecord::DISK_BYTES);
        if (skip > 0 && fseeko(file, skip, SEEK_SET) != 0) {
            Debug(Debug::ERROR) << "Cannot seek " << skip << " byte into " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
        refill();
    }

    ~RepRankBlockReader() {
        if (file != NULL) {
            fclose(file);
        }
    }

    bool next(PairRecord &row) {
        if (at >= filled) {
            refill();
            if (filled == 0) {
                return false;
            }
        }
        row = buffer[at++];
        return true;
    }

    bool fillBatch(std::vector<PairRecord> &into, size_t want, uint64_t until) {
        into.clear();
        into.swap(carry);
        PairRecord row;
        while (stopped == false && into.size() < want) {
            if (next(row) == false || row.rep() >= until) {
                stopped = true;
                break;
            }
            into.push_back(row);
        }
        if (stopped == false && into.empty() == false) {
            const uint64_t last = into.back().rep();
            size_t keep = into.size();
            while (keep > 0 && into[keep - 1].rep() == last) {
                keep--;
            }
            if (keep == 0) {
                while (true) {
                    if (next(row) == false || row.rep() >= until) {
                        stopped = true;
                        break;
                    }
                    if (row.rep() != last) {
                        carry.push_back(row);
                        break;
                    }
                    into.push_back(row);
                }
            } else {
                carry.assign(into.begin() + keep, into.end());
                into.resize(keep);
            }
        }
        return into.empty() == false;
    }

private:
    RepRankBlockReader(const RepRankBlockReader &);
    RepRankBlockReader &operator=(const RepRankBlockReader &);

    void refill() {
        filled = readRecords(buffer.data(), buffer.size(), file);
        at = 0;
    }

    std::vector<PairRecord> buffer;
    FILE *file;
    size_t at;
    size_t filled;
    bool stopped;
    std::vector<PairRecord> carry;
};

struct Candidates {
    std::vector<uint64_t> members;
    std::vector<int> diagonals;
    uint64_t queryRank;

    Candidates() : queryRank(0) {}
    void clear() {
        members.clear();
        diagonals.clear();
    }
};

// one thread's rank range of a slice, laid out and read by that thread
struct SlicePart {
    std::vector<uint64_t> ranks;
    std::vector<const char *> at;
    std::vector<IoRing::Read> reads;
    size_t bytes;
};

// the rows of one slice and every sequence they need, read once in file order
struct ReadSlice {
    size_t firstItem;
    size_t lastItem;
    std::vector<Candidates> items;
    std::vector<std::vector<uint64_t> > gathered;
    std::vector<uint64_t> splitters;
    std::vector<SlicePart> parts;
    std::vector<size_t> partOffset;
    std::vector<char> arena;

    const char *sequenceOf(uint64_t rank) const {
        const SlicePart &part = parts[std::upper_bound(splitters.begin(), splitters.end(), rank) - splitters.begin()];
        return part.at[std::lower_bound(part.ranks.begin(), part.ranks.end(), rank) - part.ranks.begin()];
    }
};

struct AlignWorker {
    AlignWorker(size_t maxLen, SubstitutionMatrix &subMat,
                SubstitutionMatrix::FastMatrix &fastMatrix, EvalueComputation &evaluer,
                const Parameters &par)
        : query(maxLen, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, par.compBiasCorrection),
          target(maxLen, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, par.compBiasCorrection),
          aligner(Parameters::DBTYPE_AMINO_ACIDS, maxLen, &subMat, &fastMatrix, &evaluer,
                  par.compBiasCorrection, par.compBiasCorrectionScale,
                  -par.gapOpen.values.aminoacid(), -par.gapExtend.values.aminoacid(),
                  BlockAligner::CLUSTER_MAX_BAND) {}
    Sequence query;
    Sequence target;
    BlockAligner aligner;
};

struct GateCounts {
    GateCounts() : seen(0), rejected(0), rescued(0), kept(0), stale(0),
                   selfRow(0), targetSkip(0), covSkip(0), querySkip(0) {}
    uint64_t seen;
    uint64_t rejected;
    uint64_t rescued;
    uint64_t kept;
    uint64_t stale;
    uint64_t selfRow;
    uint64_t targetSkip;
    uint64_t covSkip;
    uint64_t querySkip;
};

static float lookupScorePerColumn(const std::string &table, double seqId, double cov,
                                  double recall) {
    std::stringstream in(table);
    std::string line;
    const int hundredths = static_cast<int>((seqId + 0.0001) * 100);
    const float wantSeqId = static_cast<float>(hundredths - hundredths % 5) / 100;
    const float wantCov = static_cast<float>(static_cast<int>((cov + 0.0001) * 10)) / 10;
    while (std::getline(in, line)) {
        const std::vector<std::string> values = Util::split(line, " ");
        if (values.size() < 4) {
            continue;
        }
        if (MathUtil::AreSame((float) strtod(values[0].c_str(), NULL), wantCov)
            && MathUtil::AreSame((float) strtod(values[1].c_str(), NULL), wantSeqId)
            && strtod(values[3].c_str(), NULL) >= recall) {
            return (float) strtod(values[2].c_str(), NULL);
        }
    }
    Debug(Debug::WARNING) << "No calibrated score per column for coverage " << wantCov
                          << " and identity " << wantSeqId << ", gapped rescue off\n";
    return 0;
}

// the same fields Matcher::resultToBuffer writes, with a rank where a 32 bit key would not fit
static void appendAlnTextLine(std::string &into, uint64_t member, int score, float seqId, double eval,
                              int qStart, int qEnd, unsigned int qLen, int dbStart, int dbEnd,
                              unsigned int dbLen, const std::string &backtrace) {
    char buffer[192];
    const int wrote = snprintf(buffer, sizeof(buffer),
                               "%llu\t%d\t%1.3f\t%.3E\t%d\t%d\t%u\t%d\t%d\t%u\t",
                               (unsigned long long) member, score, seqId, eval, qStart, qEnd, qLen,
                               dbStart, dbEnd, dbLen);
    into.append(buffer, wrote);
    into.append(backtrace);
    into.push_back('\n');
}

// the representative is the center of its own alignment, so its line is full length and identical
static void appendSelfLine(std::string &into, uint64_t rep, uint32_t seqLen) {
    appendAlnTextLine(into, rep, 0, 1.0f, 0.0, 0, (int) seqLen - 1, seqLen, 0, (int) seqLen - 1,
                      seqLen, SSTR(seqLen) + "M");
}

static void appendRepMark(std::string &into, uint64_t rep, uint32_t seqLen) {
    into.push_back(ALN_TEXT_REP_MARK);
    into.push_back('\t');
    into.append(SSTR(rep));
    into.push_back('\n');
    appendSelfLine(into, rep, seqLen);
}

static void writeText(FILE *out, const std::string &path, const std::string &text) {
    if (text.empty() == false && fwrite(text.c_str(), 1, text.size(), out) != text.size()) {
        Debug(Debug::ERROR) << "Cannot write " << text.size() << " byte to " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
}

static bool rescueWithGaps(uint64_t member, uint32_t queryLen, uint32_t targetLen,
                           const char *querySeq, const char *targetSeq,
                           const BlockAligner::UngappedAln_res &hit, const ClusterAssignmentBitmap &assignedCluster,
                           Sequence &target, BlockAligner &aligner, const Parameters &par,
                           float scorePerColThreshold, int xDrop, GateCounts &gate,
                           std::string *line) {
    if (hit.diagonalLen <= 0
        || (float) hit.score / (float) hit.diagonalLen < scorePerColThreshold) {
        return false;
    }
    if (hit.qStart < 0 || hit.tStart < 0 || hit.alnLen < 3 || assignedCluster.isAssigned(member)) {
        return false;
    }
    int queryFrom = -1;
    int targetFrom = -1;
    for (int at = 0; at <= hit.alnLen - 3; at++) {
        const int q = hit.qStart + at;
        const int t = hit.tStart + at;
        if (querySeq[q] == targetSeq[t] && querySeq[q + 1] == targetSeq[t + 1]
            && querySeq[q + 2] == targetSeq[t + 2]) {
            queryFrom = q + 1;
            targetFrom = t + 1;
            break;
        }
    }
    if (queryFrom < 0) {
        return false;
    }
    gate.rescued++;
    std::string backtrace;
    const s_align gapped = aligner.bandedalign(&target, (size_t) queryFrom, (size_t) targetFrom,
                                               backtrace, xDrop, par.covThr, par.covMode);
    if (gapped.evalue < 0 || backtrace.empty()) {
        return false;
    }
    const unsigned int alnLen = static_cast<unsigned int>(backtrace.size());
    const float seqId = Util::computeSeqId(par.seqIdMode, gapped.identicalAACnt, queryLen, targetLen,
                                           alnLen);
    Matcher::result_t result(0, gapped.score1, gapped.qCov, gapped.tCov, seqId, gapped.evalue, alnLen,
                             gapped.qStartPos1, gapped.qEndPos1, queryLen, gapped.dbStartPos1,
                             gapped.dbEndPos1, targetLen, backtrace);
    if (Alignment::checkCriteria(result, false, par.evalThr, par.seqIdThr, par.alnLenThr, par.covMode,
                                 par.covThr) == false) {
        return false;
    }
    gate.kept++;
    if (line != NULL) {
        line->clear();
        appendAlnTextLine(*line, member, gapped.score1, seqId, gapped.evalue, gapped.qStartPos1,
                          gapped.qEndPos1, queryLen, gapped.dbStartPos1, gapped.dbEndPos1, targetLen,
                          Matcher::compressAlignment(backtrace));
    }
    return true;
}

static void collectCandidates(const Lin8DbReader &reader, uint64_t rep, const PairRecord *rows, size_t count,
                              const ClusterAssignmentBitmap &assignedCluster, const Parameters &par,
                              Candidates &candidates, GateCounts &gate) {
    candidates.clear();
    candidates.queryRank = rep;
    if (assignedCluster.isAssigned(rep)) {
        return;
    }
    const uint32_t queryLen = reader.getSeqLen(rep);
    for (size_t i = 0; i < count; i++) {
        const uint64_t member = rows[i].member();
        if (member == rep) {
            gate.selfRow++;
            continue;
        }
        if (assignedCluster.isAssigned(member)) {
            gate.targetSkip++;
            continue;
        }
        if (Util::canBeCovered(par.covThr, par.covMode, queryLen, reader.getSeqLen(member)) == false) {
            gate.covSkip++;
            continue;
        }
        candidates.members.push_back(member);
        candidates.diagonals.push_back(rows[i].diagonal());
    }
}

// each thread sorts the ranks of every threads-th item, and thread 0 picks the range splitters from its share
static void gatherRanks(ReadSlice &slice, unsigned int thread, unsigned int threads) {
    std::vector<uint64_t> &ranks = slice.gathered[thread];
    ranks.clear();
    for (size_t k = thread; k < slice.items.size(); k += threads) {
        const Candidates &item = slice.items[k];
        if (item.members.empty()) {
            continue;
        }
        ranks.push_back(item.queryRank);
        ranks.insert(ranks.end(), item.members.begin(), item.members.end());
    }
    SORT_SERIAL(ranks.begin(), ranks.end());
    ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
    if (thread == 0) {
        slice.splitters.resize(threads - 1);
        for (unsigned int t = 0; t + 1 < threads; t++) {
            slice.splitters[t] = ranks.empty() ? 0 : ranks[ranks.size() * (t + 1) / threads];
        }
    }
}

// merges this thread's rank range out of every gathered list and measures its reads
static void measurePart(const Lin8DbReader &reader, ReadSlice &slice, unsigned int thread) {
    SlicePart &part = slice.parts[thread];
    part.ranks.clear();
    for (size_t g = 0; g < slice.gathered.size(); g++) {
        const std::vector<uint64_t> &ranks = slice.gathered[g];
        std::vector<uint64_t>::const_iterator from = thread == 0 ? ranks.begin()
            : std::lower_bound(ranks.begin(), ranks.end(), slice.splitters[thread - 1]);
        std::vector<uint64_t>::const_iterator until = thread == slice.parts.size() - 1 ? ranks.end()
            : std::lower_bound(from, ranks.end(), slice.splitters[thread]);
        part.ranks.insert(part.ranks.end(), from, until);
    }
    SORT_SERIAL(part.ranks.begin(), part.ranks.end());
    part.ranks.erase(std::unique(part.ranks.begin(), part.ranks.end()), part.ranks.end());
    part.bytes = reader.layoutReads(part.ranks.data(), part.ranks.size(), NULL, part.reads, part.at);
}

static void placeParts(ReadSlice &slice) {
    slice.partOffset.resize(slice.parts.size() + 1);
    slice.partOffset[0] = 0;
    for (size_t t = 0; t < slice.parts.size(); t++) {
        slice.partOffset[t + 1] = slice.partOffset[t] + slice.parts[t].bytes;
    }
    if (slice.arena.size() < slice.partOffset.back()) {
        slice.arena.resize(slice.partOffset.back());
    }
}

static void alignMemberBatch(const Lin8DbReader &reader, uint64_t rep,
                             const ClusterAssignmentBitmap &assignedCluster, Sequence &query, Sequence &target,
                             BlockAligner &aligner, const Parameters &par, const Candidates &candidates,
                             const ReadSlice &slice, std::vector<uint64_t> &out, GateCounts &gate, float scorePerColThreshold,
                             int xDrop, std::vector<std::string> *lines) {
    const uint32_t queryLen = reader.getSeqLen(rep);
    std::string line;
    const char *querySeq = slice.sequenceOf(rep);
    query.mapSequence(0, 0, (char *) querySeq, queryLen);
    aligner.initQuery(&query);
    for (size_t k = 0; k < candidates.members.size(); k++) {
        const uint64_t member = candidates.members[k];
        if (assignedCluster.isAssigned(member)) {
            gate.stale++;
            continue;
        }
        const uint32_t targetLen = reader.getSeqLen(member);
        const char *targetSeq = slice.sequenceOf(member);
        target.mapSequence(0, 0, (char *) targetSeq, targetLen);
        const BlockAligner::UngappedAln_res hit = aligner.ungappedAlign(
            &target, static_cast<unsigned short>(candidates.diagonals[k]));
        gate.seen++;
        if (hit.eval > par.evalThr || hit.alnLen < par.alnLenThr
            || Util::hasCoverage(par.covThr, par.covMode, hit.qcov, hit.tcov) == false) {
            gate.rejected++;
            if (rescueWithGaps(member, queryLen, targetLen, querySeq, targetSeq, hit, assignedCluster,
                               target, aligner, par, scorePerColThreshold, xDrop, gate,
                               lines == NULL ? NULL : &line)) {
                out.push_back(member);
                if (lines != NULL) {
                    lines->push_back(line);
                }
            }
            continue;
        }
        int identical = 0;
        for (int q = hit.qStart; q <= hit.qEnd; q++) {
            const char a = querySeq[q] & static_cast<unsigned char>(~0x20);
            const char b =
                targetSeq[hit.tStart + (q - hit.qStart)] & static_cast<unsigned char>(~0x20);
            identical += (a == b);
        }
        const float seqId =
            Util::computeSeqId(par.seqIdMode, identical, queryLen, targetLen, hit.alnLen);
        if (seqId < par.seqIdThr - FLT_EPSILON) {
            gate.rejected++;
            if (rescueWithGaps(member, queryLen, targetLen, querySeq, targetSeq, hit, assignedCluster,
                               target, aligner, par, scorePerColThreshold, xDrop, gate,
                               lines == NULL ? NULL : &line)) {
                out.push_back(member);
                if (lines != NULL) {
                    lines->push_back(line);
                }
            }
            continue;
        }
        out.push_back(member);
        if (lines != NULL) {
            line.clear();
            appendAlnTextLine(line, member, hit.bitScore, seqId, hit.eval, hit.qStart, hit.qEnd,
                              queryLen, hit.tStart, hit.tEnd, targetLen,
                              SSTR(hit.qEnd - hit.qStart + 1) + "M");
            lines->push_back(line);
        }
    }
}

int lin8align2clust(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    FileUtil::fixRlimitNoFile();
    const NodePlacement node = NodePlacement::resolve(par);
    unsigned int writerNodes = 0;
    size_t repRankBlocks = 0;
    uint64_t ranks = 0;
    readPipelineShape(par.db2, writerNodes, repRankBlocks, ranks);
    requireEveryNodeDone(par.db2, writerNodes);

    Lin8DbReader reader(par.db1);
    reader.open();
    if (reader.getSize() != ranks) {
        Debug(Debug::ERROR) << "The database holds " << reader.getSize()
                            << " sequences and the pairs were made for " << ranks << "\n";
        EXIT(EXIT_FAILURE);
    }
    const unsigned int threads = par.threads;
    // a batch fills every lane four times over, so the fork and join are lost in the aligning
    const size_t batchRows = (size_t) threads * Lin8DbReader::LANES * MEMBERS_PER_ALIGN_BATCH * 4;
    Debug(Debug::INFO) << "Batches of " << batchRows << " rows\n";
    reader.openBatch(threads, ARENA_BYTES, Util::computeMemory(par.splitMemoryLimit),
                     Lin8DbReader::READ_AGAIN, Lin8DbReader::ACCESS_RANDOM);

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    SubstitutionMatrix::FastMatrix fastMatrix = SubstitutionMatrix::createAsciiSubMat(subMat);
    EvalueComputation evaluer(reader.getDataSize(), &subMat);
    const size_t maxLen = std::max<size_t>(reader.getIndex().getMaxSeqLen(), 1);

    size_t firstRepRankBlock = par.lin8RepRankBlock < 0 ? 0 : (size_t) par.lin8RepRankBlock;
    const size_t lastRepRankBlock =
        par.lin8RepRankBlock < 0
            ? repRankBlocks
            : (par.lin8RepRankBlockCount <= 0
                   ? repRankBlocks
                   : std::min(firstRepRankBlock + (size_t) par.lin8RepRankBlockCount, repRankBlocks));

    const bool decideHere = node.count == 1;
    if (decideHere && par.lin8RepRankBlock < 0) {
        while (firstRepRankBlock < lastRepRankBlock
               && FileUtil::fileExists((par.db4 + ".0." + SSTR(firstRepRankBlock)).c_str())) {
            firstRepRankBlock++;
        }
        if (firstRepRankBlock > 0) {
            Debug(Debug::INFO) << "Resuming at rank block " << firstRepRankBlock << "\n";
        }
    }
    Debug(Debug::INFO) << node.says("takes")
                       << node.share(lastRepRankBlock - firstRepRankBlock, repRankBlocks)
                       << " rank blocks\n";

    ClusterAssignmentBitmap assignedCluster;
    assignedCluster.open(par.db4 + ".align_assigned_" + SSTR(node.index), ranks);
    if (decideHere) {
        assignedCluster.catchUpTo(par.db4, firstRepRankBlock);
    } else {
        const size_t lookahead = (size_t) par.lin8RepRankBlockLookahead;
        const size_t floor = firstRepRankBlock > lookahead ? firstRepRankBlock - lookahead : 0;
        assignedCluster.catchUpToAvailable(par.db4, firstRepRankBlock, floor);
        assignedCluster.save(0);
    }
    const SubBucketCounts prefCounts(par.db2, writerNodes, PairRecord::REP_RANK_SUB_BLOCKS, repRankBlocks);

    Timer timer;
    uint64_t aligned = 0;
    uint64_t querySkipped = 0;
    uint64_t passed = 0;
    std::vector<PairRecord> batch;
    std::vector<size_t> starts;
    std::vector<std::vector<uint64_t> > survivors;
    std::vector<MemberBatch> work;
    std::vector<std::vector<uint64_t> > batchSurvivors;
    std::vector<size_t> workFirst;
    const bool wantText = par.includeAlignFiles;
    std::vector<std::vector<std::string> > survivorLines;
    std::vector<std::vector<std::string> > batchSurvivorLines;
    std::string text;
    double spentReading = 0, spentAligning = 0, spentDeciding = 0, spentWriting = 0;
    uint64_t wastedPasses = 0;
    double aligningThreadSeconds = 0;
    Debug::Progress progress(repRankBlocks);
    ReadSlice slices[2];
    for (int i = 0; i < 2; i++) {
        slices[i].gathered.resize(threads);
        slices[i].parts.resize(threads);
    }
    std::vector<PairRecord> outBuffer;
    std::vector<AlignWorker *> workers(threads, NULL);
    std::vector<GateCounts> gate(threads);
    const float scorePerColThreshold =
        lookupScorePerColumn(par.covMode == Parameters::COV_MODE_BIDIRECTIONAL
                                 ? getCovSeqidQscPercMinDiag()
                                 : getCovSeqidQscPercMinDiagTargetCov(),
                             par.seqIdThr, par.covThr, par.rescueRecall);
    const int xDrop = (int) BlockAligner::MIN_BAND * par.gapExtend.values.aminoacid()
                      + par.gapOpen.values.aminoacid();
    Debug(Debug::INFO) << "Gapped rescue score per column: " << scorePerColThreshold
                       << " at recall " << par.rescueRecall << "\n";

    uint64_t clusters = 0;
    uint64_t assigned = 0;

    std::vector<std::pair<std::string, std::string> > pendingOut;
    size_t pendingFirst = firstRepRankBlock;
    uint64_t pendingBytes = 0;
    for (size_t repRankBlock = firstRepRankBlock; repRankBlock < lastRepRankBlock; repRankBlock++) {
        if (decideHere == false) {
            if (FileUtil::fileExists(nodeDonePath(par.db3 + "." + SSTR(repRankBlock), node.index).c_str())) {
                progress.updateProgress(repRankBlock);
                continue;
            }
        }
        const std::vector<uint64_t> subRows = prefCounts.of(repRankBlock, repRankBlock % writerNodes);
        uint64_t blockRows = 0;
        for (size_t sub = 0; sub < subRows.size(); sub++) {
            blockRows += subRows[sub];
        }
        // cut the block where the rows are, not where the ranks are, and weigh the first node by --first-node-share
        size_t subFrom = 0, subUntil = 0;
        uint64_t skipRows = 0, seenRows = 0;
        for (size_t sub = 0, at = 0; at <= node.count; ) {
            const double before = at == 0 ? 0.0 : par.lin8FirstNodeShare + (at - 1);
            const uint64_t want = (uint64_t) (blockRows * before / (par.lin8FirstNodeShare + (node.count - 1)));
            if (seenRows >= want) {
                if (at == node.index) { subFrom = sub; skipRows = seenRows; }
                if (at == node.index + 1) { subUntil = sub; }
                at++;
                continue;
            }
            seenRows += subRows[sub++];
        }
        const size_t fineTotal = repRankBlocks * PairRecord::REP_RANK_SUB_BLOCKS;
        const uint64_t myFrom =
            ((repRankBlock * PairRecord::REP_RANK_SUB_BLOCKS + subFrom) * ranks + fineTotal - 1) / fineTotal;
        const uint64_t myUntil =
            ((repRankBlock * PairRecord::REP_RANK_SUB_BLOCKS + subUntil) * ranks + fineTotal - 1) / fineTotal;

        const std::string outPath = decideHere ? par.db4 + ".0." + SSTR(repRankBlock)
                                              : par.db3 + "." + SSTR(node.index) + "." + SSTR(repRankBlock);
        const std::string outTmp = outPath + ".tmp";
        FILE *out = FileUtil::openAndDelete(outTmp.c_str(), "w");
        const std::string textPath =
            alnTextPath(decideHere ? par.db4 : par.db3, node.index, repRankBlock);
        const std::string textTmp = textPath + ".tmp";
        FILE *textOut = wantText ? FileUtil::openAndDelete(textTmp.c_str(), "w") : NULL;
        if (textOut != NULL) {
            setvbuf(textOut, NULL, _IOFBF, 1u << 20);
        }
        outBuffer.clear();

        RepRankBlockReader stream(par.db2, writerNodes, repRankBlock, STREAM_ROWS, skipRows);

        while (true) {
            double mark = omp_get_wtime();
            if (decideHere == false) {
                const size_t lookahead = (size_t) par.lin8RepRankBlockLookahead;
                assignedCluster.catchUpToAvailable(par.db4, repRankBlock,
                                                   repRankBlock > lookahead ? repRankBlock - lookahead : 0);
            }
            const bool more = stream.fillBatch(batch, batchRows, myUntil);
            spentReading += omp_get_wtime() - mark;
            if (more == false) {
                break;
            }
            starts.clear();
            starts.push_back(0);
            for (size_t i = 1; i < batch.size(); i++) {
                if (batch[i].rep() != batch[i - 1].rep()) {
                    starts.push_back(i);
                }
            }
            starts.push_back(batch.size());
            const size_t groups = starts.size() - 1;
            if (survivors.size() < groups) {
                survivors.resize(groups);
            }
            if (wantText && survivorLines.size() < groups) {
                survivorLines.resize(groups);
            }

            work.clear();
            workFirst.assign(groups + 1, 0);
            for (size_t g = 0; g < groups; g++) {
                workFirst[g] = work.size();
                survivors[g].clear();
                if (wantText) {
                    survivorLines[g].clear();
                }
                const uint64_t rep = batch[starts[g]].rep();
                if (rep < myFrom || rep >= myUntil) {
                    continue;
                }
                if (assignedCluster.isAssigned(rep)) {
                    querySkipped += starts[g + 1] - starts[g];
                    continue;
                }
                const size_t rows = starts[g + 1] - starts[g];
                for (size_t from = 0; from < rows; from += MEMBERS_PER_ALIGN_BATCH) {
                    MemberBatch item;
                    item.group = g;
                    item.from = from;
                    item.count = std::min(MEMBERS_PER_ALIGN_BATCH, rows - from);
                    work.push_back(item);
                }
            }
            workFirst[groups] = work.size();
            if (batchSurvivors.size() < work.size()) {
                batchSurvivors.resize(work.size());
            }
            if (wantText && batchSurvivorLines.size() < work.size()) {
                batchSurvivorLines.resize(work.size());
            }

            mark = omp_get_wtime();
            const size_t slicesWanted = std::max<size_t>(1, std::min(READ_SLICES, batch.size() / SLICE_MIN_ROWS));
            const size_t sliceItems = (work.size() + slicesWanted - 1) / slicesWanted;
            const size_t sliceCount = sliceItems == 0 ? 0 : (work.size() + sliceItems - 1) / sliceItems;
#pragma omp parallel num_threads(threads)
            {
                unsigned int thread = 0;
#ifdef OPENMP
                thread = static_cast<unsigned int>(omp_get_thread_num());
#endif
                if (workers[thread] == NULL) {
                    workers[thread] = new AlignWorker(maxLen, subMat, fastMatrix, evaluer, par);
                }
                AlignWorker &worker = *workers[thread];
                double aligning = 0;
                for (size_t s = 0; s <= sliceCount; s++) {
                    const double began = omp_get_wtime();
                    double waited = 0;
                    if (s < sliceCount) {
                        ReadSlice &slice = slices[s % 2];
#pragma omp single
                        {
                            slice.firstItem = s * sliceItems;
                            slice.lastItem = std::min(work.size(), slice.firstItem + sliceItems);
                            slice.items.resize(slice.lastItem - slice.firstItem);
                        }
#pragma omp for schedule(dynamic, 16)
                        for (size_t w = slice.firstItem; w < slice.lastItem; w++) {
                            const MemberBatch &item = work[w];
                            const size_t at = starts[item.group];
                            collectCandidates(reader, batch[at].rep(), &batch[at + item.from], item.count,
                                              assignedCluster, par, slice.items[w - slice.firstItem], gate[thread]);
                        }
                        gatherRanks(slice, thread, threads);
#pragma omp barrier
                        measurePart(reader, slice, thread);
#pragma omp barrier
#pragma omp single
                        placeParts(slice);
                        SlicePart &part = slice.parts[thread];
                        reader.layoutReads(part.ranks.data(), part.ranks.size(),
                                           slice.arena.data() + slice.partOffset[thread], part.reads, part.at);
                        reader.submitReads(thread, s % 2, part.reads.data(), part.reads.size());
                    }
                    if (s > 0) {
                        ReadSlice &slice = slices[(s - 1) % 2];
                        const double stalled = omp_get_wtime();
                        reader.awaitBatch(thread, (s - 1) % 2);
#pragma omp barrier
                        waited = omp_get_wtime() - stalled;
#pragma omp for schedule(dynamic, 1)
                        for (size_t w = slice.firstItem; w < slice.lastItem; w++) {
                            const Candidates &item = slice.items[w - slice.firstItem];
                            batchSurvivors[w].clear();
                            if (wantText) {
                                batchSurvivorLines[w].clear();
                            }
                            if (item.members.empty()) {
                                continue;
                            }
                            alignMemberBatch(reader, item.queryRank, assignedCluster, worker.query, worker.target,
                                             worker.aligner, par, item, slice, batchSurvivors[w], gate[thread],
                                             scorePerColThreshold, xDrop, wantText ? &batchSurvivorLines[w] : NULL);
                        }
                    }
                    aligning += omp_get_wtime() - began - waited;
                }
#pragma omp atomic
                aligningThreadSeconds += aligning;
            }
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
            for (size_t g = 0; g < groups; g++) {
                for (size_t w = workFirst[g]; w < workFirst[g + 1]; w++) {
                    survivors[g].insert(survivors[g].end(), batchSurvivors[w].begin(),
                                        batchSurvivors[w].end());
                    if (wantText) {
                        survivorLines[g].insert(survivorLines[g].end(),
                                                batchSurvivorLines[w].begin(),
                                                batchSurvivorLines[w].end());
                    }
                }
            }
            spentAligning += omp_get_wtime() - mark;
            mark = omp_get_wtime();
            for (size_t g = 0; g < groups; g++) {
                const uint64_t rep = batch[starts[g]].rep();
                if (rep < myFrom || rep >= myUntil) {
                    continue;
                }
                aligned += starts[g + 1] - starts[g];
                passed += survivors[g].size();
                if (decideHere) {
                    if (assignedCluster.isAssigned(rep)) {
                        wastedPasses += survivors[g].size();
                    } else {
                        for (size_t k = 0; k < survivors[g].size(); k++) {
                            wastedPasses += survivors[g][k] != rep && assignedCluster.isAssigned(survivors[g][k]);
                        }
                    }
                    const size_t before = outBuffer.size();
                    const bool made = assignCluster(rep, survivors[g].data(), survivors[g].size(),
                                                    assignedCluster, outBuffer, assigned);
                    clusters += made ? 1 : 0;
                    if (wantText && made) {
                        text.clear();
                        appendRepMark(text, rep, reader.getSeqLen(rep));
                        size_t j = before + 1;
                        for (size_t k = 0; k < survivors[g].size() && j < outBuffer.size(); k++) {
                            if (survivors[g][k] == outBuffer[j].member()) {
                                text.append(survivorLines[g][k]);
                                j++;
                            }
                        }
                        writeText(textOut, textTmp, text);
                    }
                } else {
                    PairRecord line;
                    line.set(rep, rep, 0);
                    outBuffer.push_back(line);
                    if (wantText) {
                        text.clear();
                        appendSelfLine(text, rep, reader.getSeqLen(rep));
                    }
                    for (size_t k = 0; k < survivors[g].size(); k++) {
                        line.set(rep, survivors[g][k], 0);
                        outBuffer.push_back(line);
                        if (wantText) {
                            text.append(survivorLines[g][k]);
                        }
                    }
                    if (wantText) {
                        writeText(textOut, textTmp, text);
                    }
                }
                if (outBuffer.size() >= STREAM_ROWS) {
                    const double put = omp_get_wtime();
                    if (writeRecords(outBuffer.data(), outBuffer.size(), out) != outBuffer.size()) {
                        Debug(Debug::ERROR) << "Cannot write " << outTmp << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    outBuffer.clear();
                    spentWriting += omp_get_wtime() - put;
                }
            }
            spentDeciding += omp_get_wtime() - mark;
        }

        if (outBuffer.empty() == false
            && writeRecords(outBuffer.data(), outBuffer.size(), out)
                   != outBuffer.size()) {
            Debug(Debug::ERROR) << "Cannot write " << outTmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        const uint64_t outBytes = static_cast<uint64_t>(ftello(out));
        if (fclose(out) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << outTmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        // published first, because the pair file is what a rerun takes as the block being done
        if (textOut != NULL) {
            if (fclose(textOut) != 0) {
                Debug(Debug::ERROR) << "Cannot close " << textTmp << "\n";
                EXIT(EXIT_FAILURE);
            }
            FileUtil::publishAtomically(textTmp, textPath);
        }
        if (decideHere) {
            pendingOut.push_back(std::make_pair(outTmp, outPath));
            pendingBytes += outBytes;
            if (pendingBytes >= PUBLISH_BATCH_BYTES || pendingOut.size() >= PUBLISH_BATCH_FILES
                || repRankBlock + 1 == lastRepRankBlock) {
                publishAllAtomically(pendingOut, threads);
                if (par.removeTmpFiles) {
                    removeConsumedBuckets(par.db2, writerNodes, pendingFirst, repRankBlock + 1, 1);
                }
                pendingFirst = repRankBlock + 1;
                pendingBytes = 0;
            }
        } else {
            FileUtil::publishAtomically(outTmp, outPath);
            markNodeDone(par.db3 + "." + SSTR(repRankBlock), node.index);
        }
        progress.updateProgress(repRankBlock);
    }

    if (decideHere) {
        assignedCluster.save(lastRepRankBlock);
    }

    const std::string shapeTmp = (decideHere ? par.db4 : par.db3)
                                 + "." + SSTR(node.index) + ".shape.tmp";
    FILE *shape = FileUtil::openAndDelete(shapeTmp.c_str(), "w");
    if (decideHere) {
        fprintf(shape, "repRankBlocks\t%zu\nranks\t%zu\n", repRankBlocks, (size_t) ranks);
    } else {
        fprintf(shape, "nodes\t%u\nrepRankBlocks\t%zu\nranks\t%zu\n", node.count, repRankBlocks, (size_t) ranks);
    }
    if (fclose(shape) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << shapeTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(shapeTmp, decideHere ? par.db4 : par.db3);

    GateCounts all;
    for (size_t i = 0; i < gate.size(); i++) {
        all.seen += gate[i].seen;
        all.rejected += gate[i].rejected;
        all.rescued += gate[i].rescued;
        all.kept += gate[i].kept;
        all.stale += gate[i].stale;
    }
    Debug(Debug::INFO) << "Ungapped filter: " << all.seen << " pairs, "
                       << (all.seen - all.rejected + all.kept) << " accepted; " << all.rescued
                       << " of the failures retried with gaps and " << all.kept << " passed\n";
    all.selfRow = 0; all.targetSkip = 0; all.covSkip = 0; all.querySkip = 0;
    for (int i = 0; i < par.threads; i++) {
        all.selfRow += gate[i].selfRow; all.targetSkip += gate[i].targetSkip;
        all.covSkip += gate[i].covSkip; all.querySkip += gate[i].querySkip;
    }
    all.querySkip = querySkipped;
    Debug(Debug::INFO) << "Before aligning: " << all.querySkip << " rows whose representative was already clustered, "
                       << all.targetSkip << " members already clustered, " << all.selfRow
                       << " self rows, " << all.covSkip << " that no alignment could cover\n";
    if (all.stale > 0) {
        Debug(Debug::INFO) << "Raced with another node on " << all.stale << " candidate pairs it clustered first\n";
    }
    Debug(Debug::INFO) << "Time for reading: " << (uint64_t) spentReading << "s aligning: "
                       << (uint64_t) spentAligning << "s deciding: " << (uint64_t) spentDeciding
                       << "s writing: " << (uint64_t) spentWriting << "s\n";
    Debug(Debug::INFO) << "Aligning used "
                       << (spentAligning > 0 ? aligningThreadSeconds / spentAligning : 0)
                       << " of " << threads << " threads on average, "
                       << (uint64_t) (spentAligning * threads - aligningThreadSeconds)
                       << " thread seconds idle\n";
    Debug(Debug::INFO) << "Read " << aligned << " candidate rows, " << passed << " passed, in "
                       << timer.lap() << "\n";
    if (decideHere) {
        Debug(Debug::INFO) << "Made " << clusters << " clusters holding " << (clusters + assigned)
                           << " sequences\n";
        Debug(Debug::INFO) << wastedPasses << " passed alignments were for members another group in the same batch had already taken\n";
    }
    for (size_t i = 0; i < workers.size(); i++) {
        delete workers[i];
    }
    reader.close();
    return EXIT_SUCCESS;
}
