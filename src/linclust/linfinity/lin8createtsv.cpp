#include "Lin8DbReader.h"
#include "Parameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "FileUtil.h"
#include "Util.h"
#include "Timer.h"

#include <cstdio>
#include <string>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

static void writeAt(int fd, const std::string &path, std::string &line, uint64_t &at) {
    if (line.empty()) {
        return;
    }
    const ssize_t wrote = pwrite(fd, line.c_str(), line.size(), (off_t) at);
    if (wrote < 0 || (size_t) wrote != line.size()) {
        Debug(Debug::ERROR) << "Cannot write " << line.size() << " byte to " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    at += line.size();
    line.clear();
}

static const unsigned int NAME_CHUNK_BITS = 30;
static const uint64_t NAME_CHUNK = uint64_t(1) << NAME_CHUNK_BITS;

struct PackedNames {
    std::vector<std::vector<char> > chunk;
    std::vector<uint64_t> at;
    uint64_t end;
    PackedNames(size_t ranks) : at(ranks + 1, 0), end(0) {}

    void append(const char *name, size_t length) {
        uint64_t wrote = 0;
        while (wrote < length) {
            const size_t which = (size_t) (end >> NAME_CHUNK_BITS);
            if (which >= chunk.size()) {
                chunk.push_back(std::vector<char>(NAME_CHUNK));
            }
            const size_t off = (size_t) (end & (NAME_CHUNK - 1));
            const size_t take = std::min<uint64_t>(length - wrote, NAME_CHUNK - off);
            memcpy(&chunk[which][off], name + wrote, take);
            wrote += take;
            end += take;
        }
    }

    void appendTo(std::string &line, uint64_t rank) const {
        uint64_t from = at[rank];
        const uint64_t until = at[rank + 1];
        while (from < until) {
            const size_t which = (size_t) (from >> NAME_CHUNK_BITS);
            const size_t off = (size_t) (from & (NAME_CHUNK - 1));
            const size_t take = std::min<uint64_t>(until - from, NAME_CHUNK - off);
            line.append(&chunk[which][off], take);
            from += take;
        }
    }

    size_t lengthOf(uint64_t rank) const { return (size_t) (at[rank + 1] - at[rank]); }
    size_t bytes() const { return chunk.size() * NAME_CHUNK + at.size() * sizeof(uint64_t); }
};

int lin8createtsv(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    RunDbReader reader(par.db1, true);
    reader.open();

    const size_t budget = Util::computeMemory(par.splitMemoryLimit);
    const size_t need = reader.getSize() * (sizeof(uint64_t) + 24);
    if (need > budget) {
        // no room to name every sequence, so skip the tsv instead of ending the run
        Debug(Debug::WARNING) << "Naming " << reader.getSize() << " sequences would need about "
                              << (need >> 30) << " GB and the limit is " << (budget >> 30)
                              << " GB, so no tsv was written. Use the cluster database (by rank) or "
                              << "lin8-createrepseqfasta instead\n";
        reader.close();
        return EXIT_SUCCESS;
    }

    Timer timer;
    const unsigned int threads = std::max<unsigned int>(1, par.threads);
    PackedNames nameOfRank(reader.getSize());
    const size_t NAME_BATCH = 1u << 16;
    std::vector<std::string> parsed(NAME_BATCH);
    RunDbReader::HeaderStream headers(reader);
    Debug(Debug::INFO) << "Naming " << reader.getSize() << " sequences\n";
    Debug::Progress nameProgress(reader.getSize() / NAME_BATCH + 1);
    const char *begin = NULL;
    size_t length = 0;
    uint64_t rank = 0;
    while (true) {
        size_t got = 0;
        while (got < NAME_BATCH && headers.next(begin, length)) {
            if (rank + got >= reader.getSize()) {
                Debug(Debug::ERROR) << "The headers hold more entries than the " << reader.getSize()
                                    << " the sequence locator names\n";
                EXIT(EXIT_FAILURE);
            }
            parsed[got].assign(begin, length);
            got++;
        }
        if (got == 0) {
            break;
        }
        nameProgress.updateProgress();
#pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t i = 0; i < got; i++) {
            parsed[i] = Util::parseFastaHeader(parsed[i].c_str());
        }
        for (size_t i = 0; i < got; i++) {
            nameOfRank.at[rank + i] = nameOfRank.end;
            nameOfRank.append(parsed[i].data(), parsed[i].size());
        }
        rank += got;
    }
    if (rank != reader.getSize()) {
        Debug(Debug::ERROR) << "The headers hold " << rank << " entries and the sequence locator names "
                            << reader.getSize() << "\n";
        EXIT(EXIT_FAILURE);
    }
    nameOfRank.at[reader.getSize()] = nameOfRank.end;
    Debug(Debug::INFO) << "Read " << rank << " names in " << timer.lap() << ", " << (nameOfRank.bytes() >> 30) << " GB\n";

    DBReader<DBKeyType> clusters(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                 DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    clusters.open(DBReader<DBKeyType>::LINEAR_ACCCESS);

    std::vector<size_t> edge(threads + 1, 0);
    for (unsigned int t = 0; t <= threads; t++) {
        edge[t] = clusters.getSize() * t / threads;
    }
    std::vector<uint64_t> bytesOf(threads, 0);
    std::vector<uint64_t> rowsOf(threads, 0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (unsigned int t = 0; t < threads; t++) {
        uint64_t sum = 0;
        uint64_t mine = 0;
        for (size_t i = edge[t]; i < edge[t + 1]; i++) {
            const uint64_t rep = clusters.getDbKey(i);
            if (rep >= reader.getSize()) {
                Debug(Debug::ERROR) << "The clustering names rank " << rep << ", past the database\n";
                EXIT(EXIT_FAILURE);
            }
            char *data = clusters.getData(i, t);
            while (data != NULL && *data != '\0') {
                const uint64_t member = strtoull(data, NULL, 10);
                if (member >= reader.getSize()) {
                    Debug(Debug::ERROR) << "The clustering names rank " << member
                                        << ", past the database\n";
                    EXIT(EXIT_FAILURE);
                }
                sum += nameOfRank.lengthOf(rep) + nameOfRank.lengthOf(member) + 2;
                mine++;
                data = Util::skipLine(data);
            }
        }
        bytesOf[t] = sum;
        rowsOf[t] = mine;
    }
    std::vector<uint64_t> startOf(threads + 1, 0);
    for (unsigned int t = 0; t < threads; t++) {
        startOf[t + 1] = startOf[t] + bytesOf[t];
    }

    const std::string tmp = par.db3 + ".tmp";
    const int out = open(tmp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
    Debug(Debug::INFO) << "Writing " << (startOf[threads] >> 20) << " MB of rows in "
                       << threads << " parts\n";
    Debug::Progress writeProgress(threads);
    if (out < 0 || ftruncate(out, (off_t) startOf[threads]) != 0) {
        Debug(Debug::ERROR) << "Cannot make " << tmp << " of " << startOf[threads] << " byte\n";
        EXIT(EXIT_FAILURE);
    }
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (unsigned int t = 0; t < threads; t++) {
        uint64_t at = startOf[t];
        std::string line;
        line.reserve(1u << 20);
        for (size_t i = edge[t]; i < edge[t + 1]; i++) {
            const uint64_t rep = clusters.getDbKey(i);
            char *data = clusters.getData(i, t);
            while (data != NULL && *data != '\0') {
                const uint64_t member = strtoull(data, NULL, 10);
                nameOfRank.appendTo(line, rep);
                line.push_back('\t');
                nameOfRank.appendTo(line, member);
                line.push_back('\n');
                data = Util::skipLine(data);
                if (line.size() >= (1u << 20)) {
                    writeAt(out, tmp, line, at);
                }
            }
        }
        writeAt(out, tmp, line, at);
        if (at != startOf[t + 1]) {
            Debug(Debug::ERROR) << "Thread " << t << " wrote to " << at << " and was sized to "
                                << startOf[t + 1] << "\n";
            EXIT(EXIT_FAILURE);
        }
        writeProgress.updateProgress();
    }
    if (::close(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(tmp, par.db3);
    uint64_t rows = 0;
    for (unsigned int t = 0; t < threads; t++) {
        rows += rowsOf[t];
    }

    clusters.close();
    reader.close();
    Debug(Debug::INFO) << "Wrote " << rows << " rows in " << timer.lap() << "\n";
    return EXIT_SUCCESS;
}
