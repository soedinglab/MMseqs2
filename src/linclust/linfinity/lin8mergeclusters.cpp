
#include "Lin8Db.h"

#include <cstring>
#include "Parameters.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include "Timer.h"
#include "FastSort.h"
#include <algorithm>
#include <cstdio>
#include <vector>
#include "DBWriter.h"
#include "itoa.h"

static const uint64_t LINCLUSTHASH_MAGIC = 0x4C494E4348504153ull;
struct PairFileHeader {
    uint64_t magic;
    uint64_t version;
    uint64_t keyWidth;
    uint64_t pairs;
};

struct JoinRow {
    uint64_t on;
    uint64_t carried;

    static const size_t DISK_BYTES = 16;
    void pack(unsigned char *to) const { memcpy(to, this, DISK_BYTES); }
    void unpack(const unsigned char *from) { memcpy(this, from, DISK_BYTES); }

    static bool byJoinKey(const JoinRow &first, const JoinRow &second) {
        if (first.on != second.on) {
            return first.on < second.on;
        }
        return first.carried < second.carried;
    }
};

static void readJoinRowsFromFile(const std::string &path, std::vector<JoinRow> &into) {
    into.clear();
    FILE *in = fopen(path.c_str(), "r");
    if (in == NULL) {
        return;
    }
    std::vector<JoinRow> buffer(1u << 16);
    size_t read = 0;
    while ((read = fread(buffer.data(), sizeof(JoinRow), buffer.size(), in)) > 0) {
        into.insert(into.end(), buffer.begin(), buffer.begin() + read);
    }
    if (ferror(in) != 0) {
        Debug(Debug::ERROR) << "Cannot read " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(in);
}

int lin8mergehashredundancy(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    size_t repRankBlocks = 0;
    size_t ranks = 0;
    FILE *shape = fopen(par.db1.c_str(), "r");
    if (shape == NULL || fscanf(shape, "repRankBlocks\t%zu\nranks\t%zu", &repRankBlocks, &ranks) != 2
        || repRankBlocks == 0) {
        Debug(Debug::ERROR) << "Cannot read " << par.db1 << ". Run lin8align2clustmulti first\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(shape);

    const size_t budget = static_cast<size_t>(Util::computeMemory(par.splitMemoryLimit) * 0.95);
    const std::string byMember = par.db3 + ".member";
    const std::string byKept = par.db3 + ".kept";
    const std::string extra = par.db3 + ".extra";
    Timer timer;

    {
        BucketWriter<JoinRow> members(byMember, repRankBlocks, par.threads, budget);
        BucketWriter<JoinRow> kept(byKept, repRankBlocks, par.threads, budget);
        std::vector<uint64_t> nothing(repRankBlocks, 0);
        members.openAt(nothing, nothing);
        kept.openAt(nothing, nothing);

        Debug(Debug::INFO) << "Indexing clusters and their near-duplicates in " << repRankBlocks << " rank groups\n";
        Debug::Progress routeProgress(repRankBlocks);
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (size_t repRankBlock = 0; repRankBlock < repRankBlocks; repRankBlock++) {
            unsigned int thread = 0;
#ifdef OPENMP
            thread = static_cast<unsigned int>(omp_get_thread_num());
#endif
            std::vector<PairRecord> pairs(1u << 16);
            const std::string path = par.db1 + ".0." + SSTR(repRankBlock);
            FILE *in = fopen(path.c_str(), "r");
            if (in == NULL) {
                Debug(Debug::ERROR) << "Cannot open " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            size_t read = 0;
            while ((read = readRecords(pairs.data(), pairs.size(), in)) > 0) {
                for (size_t k = 0; k < read; k++) {
                    JoinRow row;
                    row.on = pairs[k].member();
                    row.carried = pairs[k].rep();
                    members.add(thread, row, PairRecord::repRankBlockOf(row.on, ranks, repRankBlocks));
                }
            }
            if (ferror(in) != 0) {
                Debug(Debug::ERROR) << "Cannot read " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            fclose(in);
            routeProgress.updateProgress();
        }
        Debug(Debug::INFO) << "  indexed cluster members: " << timer.lap() << "\n";

        FILE *probe = fopen(par.db2.c_str(), "r");
        if (probe == NULL) {
            Debug(Debug::ERROR) << "Cannot open " << par.db2 << ". Run lin8clusthash first\n";
            EXIT(EXIT_FAILURE);
        }
        PairFileHeader header;
        if (fread(&header, sizeof(PairFileHeader), 1, probe) != 1
            || header.magic != LINCLUSTHASH_MAGIC) {
            Debug(Debug::ERROR) << par.db2 << " is not a set of redundancy pairs\n";
            EXIT(EXIT_FAILURE);
        }
        fclose(probe);
        const size_t redundantCount = header.pairs;
        const size_t REDUNDANCY_CHUNK = 1u << 22;
        const size_t redundancyChunks =
            redundantCount == 0 ? 0 : (redundantCount + REDUNDANCY_CHUNK - 1) / REDUNDANCY_CHUNK;
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads)
        for (size_t c = 0; c < redundancyChunks; c++) {
            unsigned int thread = 0;
#ifdef OPENMP
            thread = static_cast<unsigned int>(omp_get_thread_num());
#endif
            const size_t from = c * REDUNDANCY_CHUNK;
            const size_t take = std::min<size_t>(REDUNDANCY_CHUNK, redundantCount - from);
            FILE *in = fopen(par.db2.c_str(), "r");
            if (in == NULL
                || fseeko(in, (off_t) sizeof(PairFileHeader) + (off_t) (from * sizeof(JoinRow)), SEEK_SET) != 0) {
                Debug(Debug::ERROR) << "Cannot read " << par.db2 << "\n";
                EXIT(EXIT_FAILURE);
            }
            std::vector<JoinRow> redundancy(take);
            if (fread(redundancy.data(), sizeof(JoinRow), take, in) != take) {
                Debug(Debug::ERROR) << "Cannot read " << par.db2 << "\n";
                EXIT(EXIT_FAILURE);
            }
            fclose(in);
            for (size_t k = 0; k < take; k++) {
                JoinRow row;
                row.on = redundancy[k].carried;
                row.carried = redundancy[k].on;
                kept.add(thread, row, PairRecord::repRankBlockOf(row.on, ranks, repRankBlocks));
            }
        }
        members.flushAll(par.threads);
        kept.flushAll(par.threads);
        members.close();
        kept.close();
        Debug(Debug::INFO) << "  indexed near-duplicates: " << timer.lap() << "\n";
    }

    uint64_t added = 0;
    {
        BucketWriter<PairRecord> writer(extra, repRankBlocks, par.threads, budget);
        std::vector<uint64_t> nothing(repRankBlocks, 0);
        writer.openAt(nothing, nothing);
        Debug(Debug::INFO) << "Matching near-duplicates to their clusters\n";
        Debug::Progress backProgress(repRankBlocks);
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads) reduction(+ : added)
        for (size_t repRankBlock = 0; repRankBlock < repRankBlocks; repRankBlock++) {
            unsigned int thread = 0;
#ifdef OPENMP
            thread = static_cast<unsigned int>(omp_get_thread_num());
#endif
            std::vector<JoinRow> left;
            std::vector<JoinRow> right;
            readJoinRowsFromFile(byKept + "." + SSTR(repRankBlock), right);
            if (right.empty()) {
                backProgress.updateProgress();
                continue;
            }
            readJoinRowsFromFile(byMember + "." + SSTR(repRankBlock), left);
            SORT_SERIAL(left.begin(), left.end(), JoinRow::byJoinKey);
            SORT_SERIAL(right.begin(), right.end(), JoinRow::byJoinKey);
            size_t at = 0;
            for (size_t i = 0; i < right.size(); i++) {
                while (at < left.size() && left[at].on < right[i].on) {
                    at++;
                }
                for (size_t k = at; k < left.size() && left[k].on == right[i].on; k++) {
                    PairRecord row;
                    row.set(left[k].carried, right[i].carried, 0);
                    writer.add(thread, row, PairRecord::repRankBlockOf(left[k].carried, ranks, repRankBlocks));
                    added++;
                }
            }
            backProgress.updateProgress();
        }
        writer.flushAll(par.threads);
        writer.close();
        Debug(Debug::INFO) << "  matched near-duplicates to clusters: " << timer.lap() << "\n";
    }

    uint64_t rows = 0;
    Debug(Debug::INFO) << "Adding near-duplicates back into their clusters\n";
    Debug::Progress mergeProgress(repRankBlocks);
#pragma omp parallel for schedule(dynamic, 1) num_threads(par.threads) reduction(+ : rows)
    for (size_t repRankBlock = 0; repRankBlock < repRankBlocks; repRankBlock++) {
        // the additions are small next to clu_accepted, so sort them by rank
        std::vector<PairRecord> add;
        std::vector<PairRecord> buffer(1u << 16);
        FILE *extraIn = fopen((extra + "." + SSTR(repRankBlock)).c_str(), "r");
        if (extraIn != NULL) {
            size_t read = 0;
            while ((read = readRecords(buffer.data(), buffer.size(), extraIn)) > 0) {
                add.insert(add.end(), buffer.begin(), buffer.begin() + read);
            }
            fclose(extraIn);
        }
        SORT_SERIAL(add.begin(), add.end(), PairRecord::byRepAndMember);

        const std::string outPath = par.db3 + ".0." + SSTR(repRankBlock);
        const std::string outTmp = outPath + ".tmp";
        FILE *out = FileUtil::openAndDelete(outTmp.c_str(), "w");
        std::vector<PairRecord> outBuf;
        outBuf.reserve(1u << 16);
        uint64_t written = 0;
        size_t at = 0;
        bool haveRep = false;
        uint64_t lastRep = 0;
        // clu_accepted is already grouped by representative rank, so splice the sorted additions in
        FILE *in = fopen((par.db1 + ".0." + SSTR(repRankBlock)).c_str(), "r");
        if (in != NULL) {
            size_t read = 0;
            while ((read = readRecords(buffer.data(), buffer.size(), in)) > 0) {
                for (size_t k = 0; k < read; k++) {
                    const uint64_t rep = buffer[k].rep();
                    if (haveRep && rep < lastRep) {
                        Debug(Debug::ERROR) << "clu_accepted is not in representative order at block "
                                            << repRankBlock << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    while (at < add.size() && add[at].rep() < rep) {
                        outBuf.push_back(add[at++]);
                    }
                    outBuf.push_back(buffer[k]);
                    lastRep = rep;
                    haveRep = true;
                    if (outBuf.size() >= (1u << 16)) {
                        if (writeRecords(outBuf.data(), outBuf.size(), out) != outBuf.size()) {
                            Debug(Debug::ERROR) << "Cannot write " << outTmp << "\n";
                            EXIT(EXIT_FAILURE);
                        }
                        written += outBuf.size();
                        outBuf.clear();
                    }
                }
            }
            fclose(in);
        }
        while (at < add.size()) {
            outBuf.push_back(add[at++]);
        }
        if (outBuf.empty() == false) {
            if (writeRecords(outBuf.data(), outBuf.size(), out) != outBuf.size()) {
                Debug(Debug::ERROR) << "Cannot write " << outTmp << "\n";
                EXIT(EXIT_FAILURE);
            }
            written += outBuf.size();
        }
        if (fclose(out) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << outTmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        FileUtil::publishAtomically(outTmp, outPath);
        rows += written;
        mergeProgress.updateProgress();
    }

    const std::string shapeTmp = par.db3 + ".shape.tmp";
    FILE *shapeOut = FileUtil::openAndDelete(shapeTmp.c_str(), "w");
    fprintf(shapeOut, "repRankBlocks\t%zu\nranks\t%zu\n", repRankBlocks, ranks);
    if (fclose(shapeOut) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << shapeTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(shapeTmp, par.db3);

    Debug(Debug::INFO) << "Added " << added << " near-duplicate sequences back, " << rows
                       << " cluster rows, in " << timer.lap() << "\n";
    return EXIT_SUCCESS;
}

static void appendDecimalKey(std::string &into, uint64_t key) {
    char buffer[32];
    char *end = Itoa::u64toa_sse2(key, buffer);
    into.append(buffer, end - buffer - 1);
    into.push_back('\n');
}

int lin8createclusterdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    size_t repRankBlocks = 0;
    size_t ranks = 0;
    FILE *shape = fopen(par.db1.c_str(), "r");
    if (shape == NULL || fscanf(shape, "repRankBlocks\t%zu\nranks\t%zu", &repRankBlocks, &ranks) != 2
        || repRankBlocks == 0) {
        Debug(Debug::ERROR) << "Cannot read " << par.db1 << ". Run lin8align2clustmulti first\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(shape);
    if (ranks > (size_t) DB_KEY_INVALID) {
        Debug(Debug::ERROR) << ranks << " sequences do not fit a cluster key. "
                            << "Rebuild with -DMMSEQS_INT64_IDS=1\n";
        EXIT(EXIT_FAILURE);
    }

    const unsigned int threads = std::max<unsigned int>(1, par.threads);
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), threads, par.compressed,
                    Parameters::DBTYPE_CLUSTER_RES);
    writer.open();

    Timer timer;
    uint64_t clusters = 0;
    uint64_t members = 0;
    Debug::Progress progress(repRankBlocks);
#pragma omp parallel for schedule(static, 1) num_threads(threads) reduction(+ : clusters, members)
    for (unsigned int part = 0; part < threads; part++) {
        const size_t from = repRankBlocks * part / threads;
        const size_t until = repRankBlocks * (part + 1) / threads;
        uint64_t lastRep = 0;
        bool haveRep = false;
        std::string entry;
        std::vector<PairRecord> buffer(1u << 16);
        for (size_t repRankBlock = from; repRankBlock < until; repRankBlock++) {
            const std::string path = par.db1 + ".0." + SSTR(repRankBlock);
            FILE *in = fopen(path.c_str(), "r");
            if (in == NULL) {
                Debug(Debug::ERROR) << "Cannot open " << path << ", which repRankBlock " << repRankBlock
                                    << " should have decided\n";
                EXIT(EXIT_FAILURE);
            }
            size_t read = 0;
            while ((read = readRecords(buffer.data(), buffer.size(), in)) > 0) {
                for (size_t k = 0; k < read; k++) {
                    const uint64_t rep = buffer[k].rep();
                    const uint64_t member = buffer[k].member();
                    if (haveRep == false || rep != lastRep) {
                        if (haveRep) {
                            writer.writeData(entry.c_str(), entry.length(), lastRep, part);
                        }
                        if (haveRep && rep < lastRep) {
                            Debug(Debug::ERROR) << "The pairs go back from representative " << lastRep
                                                << " to " << rep << ", so they are not in order\n";
                            EXIT(EXIT_FAILURE);
                        }
                        entry.clear();
                        appendDecimalKey(entry, rep);
                        lastRep = rep;
                        haveRep = true;
                        clusters++;
                    }
                    if (member != rep) {
                        appendDecimalKey(entry, member);
                        members++;
                    }
                }
            }
            if (ferror(in) != 0) {
                Debug(Debug::ERROR) << "Cannot read " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            fclose(in);
            progress.updateProgress();
        }
        if (haveRep) {
            writer.writeData(entry.c_str(), entry.length(), lastRep, part);
        }
    }
    writer.close(false, false);

    Debug(Debug::INFO) << "Wrote " << clusters << " clusters holding " << (clusters + members)
                       << " sequences in " << timer.lap() << "\n";
    return EXIT_SUCCESS;
}
