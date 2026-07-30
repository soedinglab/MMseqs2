#include "PartitionSequences.h"

#include "Debug.h"
#include "DenseIndex.h"
#include "Util.h"

#include <algorithm>
#include <cerrno>
#include <cstring>

#include <fcntl.h>
#include <unistd.h>

namespace {
// Index entries within this distance are fetched in one read rather than two.
// The dense index is 12 B per entry, so this coalesces anything within ~5000
// keys -- which, because keys are length-ranked, is a very short stretch of the
// length distribution and therefore extremely common between neighbouring hits.
const size_t INDEX_COALESCE_BYTES = 64 * 1024;
// Same idea on the data file, where entries are ~250 B, so this bridges gaps of
// a few thousand sequences.
const size_t DATA_COALESCE_BYTES = 1024 * 1024;
}  // namespace

PartitionSequences::PartitionSequences(const std::string &dbName)
    : dbName(dbName), dataFd(-1), indexFd(-1), bytesRead(0) {
    const DenseIndex::Info info = DenseIndex::readInfo(dbName);
    entryCount = info.entryCount;
    firstKey = info.firstKey;

    dataFd = open(dbName.c_str(), O_RDONLY);
    if (dataFd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << dbName << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    const std::string indexName = DenseIndex::fileName(dbName);
    indexFd = open(indexName.c_str(), O_RDONLY);
    if (indexFd < 0) {
        Debug(Debug::ERROR) << "Cannot open " << indexName << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
}

PartitionSequences::~PartitionSequences() {
    if (dataFd >= 0) {
        close(dataFd);
    }
    if (indexFd >= 0) {
        close(indexFd);
    }
}

void PartitionSequences::readAt(int fd, void *dst, size_t length, size_t offset, const char *what) {
    char *cursor = static_cast<char *>(dst);
    size_t done = 0;
    while (done < length) {
        const ssize_t got = pread(fd, cursor + done, length - done, static_cast<off_t>(offset + done));
        if (got < 0) {
            if (errno == EINTR) {
                continue;
            }
            Debug(Debug::ERROR) << "Cannot read " << what << " from " << dbName << ": "
                                << strerror(errno) << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (got == 0) {
            Debug(Debug::ERROR) << "Unexpected end of file reading " << what << " from " << dbName
                                << " at offset " << (offset + done) << "\n";
            EXIT(EXIT_FAILURE);
        }
        done += static_cast<size_t>(got);
    }
    bytesRead += length;
}

void PartitionSequences::load(const std::vector<uint64_t> &sortedKeys) {
    keys = sortedKeys;
    offsets.assign(keys.size(), 0);
    lengths.assign(keys.size(), 0);
    arena.clear();
    bytesRead = 0;
    if (keys.empty()) {
        return;
    }

    // Pass 1: the index entries, read in coalesced ascending runs.
    std::vector<DenseIndex::Entry> entries(keys.size());
    size_t i = 0;
    std::vector<char> block;
    while (i < keys.size()) {
        const uint64_t firstRow = keys[i] - firstKey;
        size_t j = i;
        while (j + 1 < keys.size() &&
               (keys[j + 1] - keys[i]) * sizeof(DenseIndex::Entry) < INDEX_COALESCE_BYTES) {
            j++;
        }
        const uint64_t lastRow = keys[j] - firstKey;
        const size_t span = static_cast<size_t>(lastRow - firstRow + 1) * sizeof(DenseIndex::Entry);
        block.resize(span);
        readAt(indexFd, block.data(), span, DenseIndex::entryOffset(firstRow), "index entries");
        for (size_t k = i; k <= j; k++) {
            const size_t at = static_cast<size_t>(keys[k] - firstKey - firstRow) * sizeof(DenseIndex::Entry);
            memcpy(&entries[k], block.data() + at, sizeof(DenseIndex::Entry));
        }
        i = j + 1;
    }

    // Pass 2: the residues themselves, again in coalesced ascending runs. Dense
    // keys make ascending key mean ascending offset, so this is a forward scan.
    uint64_t total = 0;
    for (size_t k = 0; k < keys.size(); k++) {
        lengths[k] = entries[k].length > 0 ? entries[k].length - 1 : 0;  // drop the trailing newline
        offsets[k] = total;
        total += lengths[k];
    }
    arena.resize(total);

    i = 0;
    while (i < keys.size()) {
        size_t j = i;
        while (j + 1 < keys.size() &&
               entries[j + 1].offset >= entries[j].offset &&
               (entries[j + 1].offset - entries[i].offset) < DATA_COALESCE_BYTES) {
            j++;
        }
        const uint64_t from = entries[i].offset;
        const uint64_t to = entries[j].offset + entries[j].length;
        block.resize(static_cast<size_t>(to - from));
        readAt(dataFd, block.data(), block.size(), from, "sequence data");
        for (size_t k = i; k <= j; k++) {
            memcpy(arena.data() + offsets[k], block.data() + (entries[k].offset - from), lengths[k]);
        }
        i = j + 1;
    }
}

const char *PartitionSequences::get(uint64_t key, unsigned int *length) const {
    const std::vector<uint64_t>::const_iterator it =
        std::lower_bound(keys.begin(), keys.end(), key);
    if (it == keys.end() || *it != key) {
        return NULL;
    }
    const size_t idx = static_cast<size_t>(it - keys.begin());
    *length = lengths[idx];
    return arena.data() + offsets[idx];
}
