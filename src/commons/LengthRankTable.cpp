#include "LengthRankTable.h"

#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>

#include <unistd.h>

namespace {

struct Header {
    uint64_t magic;
    uint32_t version;
    uint32_t reserved;
    uint64_t maxLength;
    uint64_t entryCount;
    uint64_t runCount;
};

}  // namespace

// Out-of-class definitions: C++14 has no inline variables, so an in-class
// initialiser alone is not a definition and any use that takes a reference --
// std::min, a container insert -- fails to link.
const uint64_t LengthRankTable::MAX_SEQUENCE_LENGTH;
const uint64_t LengthRankTable::MAGIC;
const uint32_t LengthRankTable::VERSION;

std::string LengthRankTable::fileName(const std::string &dbPath) { return dbPath + ".lenrank"; }

bool LengthRankTable::exists(const std::string &dbPath) {
    return FileUtil::fileExists(fileName(dbPath).c_str());
}

void LengthRankTable::write(const std::string &dbPath, const std::vector<Run> &runs,
                            uint64_t entryCount) {
    std::vector<Entry> table;
    table.reserve(runs.size());
    uint64_t maxLength = 0;
    uint64_t seen = 0;
    for (size_t i = 0; i < runs.size(); i++) {
        // Strictly descending lengths tiling key space with no hole is what makes
        // the lookup a plain binary search. Both are guaranteed by the planner;
        // checking here turns a planner regression into an immediate failure
        // rather than a subtly wrong length two stages later.
        if (i > 0 && runs[i].length >= runs[i - 1].length) {
            Debug(Debug::ERROR) << "Length runs are not strictly descending at " << i << ": "
                                << runs[i - 1].length << " then " << runs[i].length << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (i > 0 && runs[i].firstKey != runs[i - 1].firstKey + runs[i - 1].count) {
            Debug(Debug::ERROR) << "Length runs leave a hole in key space at " << i << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (runs[i].count == 0) {
            Debug(Debug::ERROR) << "Length run " << i << " is empty\n";
            EXIT(EXIT_FAILURE);
        }
        maxLength = std::max(maxLength, runs[i].length);
        seen += runs[i].count;

        Entry entry;
        entry.firstKey = runs[i].firstKey;
        // Widened to 32 bits so the table itself stays usable for databases the
        // k-mer stages would refuse; the 65535 cap is enforced where the 16-bit
        // fields actually live, not here.
        if (runs[i].length > 0xFFFFFFFFULL) {
            Debug(Debug::ERROR) << "Sequence length " << runs[i].length
                                << " exceeds what the length-rank table can record\n";
            EXIT(EXIT_FAILURE);
        }
        entry.length = static_cast<uint32_t>(runs[i].length);
        entry.reserved = 0;
        table.push_back(entry);
    }
    if (seen != entryCount) {
        Debug(Debug::ERROR) << "Length runs cover " << seen << " sequences but the database has "
                            << entryCount << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (runs.empty() == false && runs[0].firstKey != 0) {
        Debug(Debug::ERROR) << "Length runs do not start at key 0\n";
        EXIT(EXIT_FAILURE);
    }

    Header header;
    header.magic = MAGIC;
    header.version = VERSION;
    header.reserved = 0;
    header.maxLength = maxLength;
    header.entryCount = entryCount;
    header.runCount = table.size();

    const std::string finalPath = fileName(dbPath);
    const std::string tmpPath = finalPath + ".tmp." + SSTR(getpid());
    FILE *file = fopen(tmpPath.c_str(), "wb");
    if (file == NULL) {
        Debug(Debug::ERROR) << "Cannot open " << tmpPath << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    bool ok = fwrite(&header, sizeof(Header), 1, file) == 1;
    if (ok && table.empty() == false) {
        ok = fwrite(table.data(), sizeof(Entry), table.size(), file) == table.size();
    }
    if (ok == false || fflush(file) != 0 || fsync(fileno(file)) != 0) {
        Debug(Debug::ERROR) << "Cannot write " << tmpPath << ": " << strerror(errno) << "\n";
        fclose(file);
        EXIT(EXIT_FAILURE);
    }
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << tmpPath << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::move(tmpPath.c_str(), finalPath.c_str());
}

void LengthRankTable::open(const std::string &dbPath) {
    const std::string path = fileName(dbPath);
    FILE *file = fopen(path.c_str(), "rb");
    if (file == NULL) {
        Debug(Debug::ERROR) << "Cannot open length-rank table " << path << ": " << strerror(errno)
                            << "\nIt is written by createdbparallel; a database built before this "
                               "existed must be rebuilt.\n";
        EXIT(EXIT_FAILURE);
    }
    Header header;
    if (fread(&header, sizeof(Header), 1, file) != 1) {
        Debug(Debug::ERROR) << "Length-rank table " << path << " is truncated\n";
        fclose(file);
        EXIT(EXIT_FAILURE);
    }
    if (header.magic != MAGIC || header.version != VERSION) {
        Debug(Debug::ERROR) << "Length-rank table " << path << " has a bad magic or version\n";
        fclose(file);
        EXIT(EXIT_FAILURE);
    }
    maxLength = header.maxLength;
    entryCount = header.entryCount;
    entries.resize(static_cast<size_t>(header.runCount));
    if (entries.empty() == false &&
        fread(entries.data(), sizeof(Entry), entries.size(), file) != entries.size()) {
        Debug(Debug::ERROR) << "Length-rank table " << path << " is short: expected "
                            << entries.size() << " runs\n";
        fclose(file);
        EXIT(EXIT_FAILURE);
    }
    // Nothing may follow. A longer file means this is not the format we think it
    // is, and reading it as though it were would give plausible wrong lengths.
    char extra = 0;
    const bool trailing = fread(&extra, 1, 1, file) == 1;
    fclose(file);
    if (trailing) {
        Debug(Debug::ERROR) << "Length-rank table " << path << " has trailing bytes\n";
        EXIT(EXIT_FAILURE);
    }
    // The invariants the binary search depends on. A table that violates them
    // returns plausible wrong lengths rather than failing, so they are checked on
    // open rather than assumed -- it costs one pass over at most a few thousand
    // entries, once per process.
    for (size_t i = 0; i < entries.size(); i++) {
        if (i > 0 && (entries[i].firstKey <= entries[i - 1].firstKey ||
                      entries[i].length >= entries[i - 1].length)) {
            Debug(Debug::ERROR) << "Length-rank table " << path << " is not ordered at run " << i
                                << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (entries[i].firstKey >= entryCount) {
            Debug(Debug::ERROR) << "Length-rank table " << path << " run " << i
                                << " starts past the end of the database\n";
            EXIT(EXIT_FAILURE);
        }
    }
    if (entries.empty() == false && entries[0].firstKey != 0) {
        Debug(Debug::ERROR) << "Length-rank table " << path << " does not start at key 0\n";
        EXIT(EXIT_FAILURE);
    }
}

bool LengthRankTable::tryLengthOf(uint64_t key, unsigned int &length) const {
    if (entries.empty() || key >= entryCount) {
        return false;
    }
    // Runs are ordered by ascending firstKey, so the run owning a key is the last
    // one that starts at or below it.
    size_t lo = 0;
    size_t hi = entries.size();
    while (lo + 1 < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        if (entries[mid].firstKey <= key) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    length = static_cast<unsigned int>(entries[lo].length);
    return true;
}

unsigned int LengthRankTable::lengthOf(uint64_t key) const {
    unsigned int length = 0;
    if (tryLengthOf(key, length) == false) {
        Debug(Debug::ERROR) << "Key " << key << " has no length in a table of " << entryCount
                            << " sequences\n";
        EXIT(EXIT_FAILURE);
    }
    return length;
}
