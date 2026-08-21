#include "DenseIndex.h"

#include "Debug.h"
#include "FileUtil.h"
#include "MemoryMapped.h"
#include "Util.h"

#include <cerrno>
#include <cstring>

#include <unistd.h>

std::string DenseIndex::fileName(const std::string &dbName) {
    return dbName + ".index.bin";
}

// Validated, not merely present.
//
// build() and createEmpty() used to write straight to the final path, so a worker
// killed part-way left a file that exists() accepted and every reader then failed
// on -- a zero header, a truncated entry array, or both. They now publish by
// atomic rename, but a database built by an older build (or by a killed run of
// one) can still be sitting there, so this checks what it is looking at: the
// magic, the version, and that the file is exactly as long as its own entry count
// says it should be.
bool DenseIndex::exists(const std::string &dbName) {
    const std::string path = fileName(dbName);
    if (FileUtil::fileExists(path.c_str()) == false) {
        return false;
    }
    FILE *file = fopen(path.c_str(), "rb");
    if (file == NULL) {
        return false;
    }
    Header header;
    const bool readHeader = fread(&header, sizeof(header), 1, file) == 1;
    fclose(file);
    if (readHeader == false || header.magic != MAGIC || header.version != VERSION) {
        return false;
    }
    return FileUtil::getFileSize(path) == entryOffset(header.entryCount);
}

void DenseIndex::build(const std::string &dbName) {
    const std::string textIndexName = dbName + ".index";
    MemoryMapped indexData(textIndexName, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (indexData.isValid() == false) {
        Debug(Debug::ERROR) << "Cannot open index file " << textIndexName << "\n";
        EXIT(EXIT_FAILURE);
    }

    const std::string outputName = fileName(dbName);
    // Written to a private sibling and renamed on success. Writing the final path
    // directly meant an interrupted build left a file that exists() accepted --
    // header zeroed, because the real header is only written at the end -- and no
    // restart could repair, since the path was there.
    const std::string tmpName = outputName + ".tmp." + SSTR(getpid());
    FILE *out = fopen(tmpName.c_str(), "wb");
    if (out == NULL) {
        Debug(Debug::ERROR) << "Cannot write dense index " << tmpName << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }

    // Reserve the header; it is rewritten at the end once the totals are known.
    Header header;
    memset(&header, 0, sizeof(header));
    if (fwrite(&header, sizeof(header), 1, out) != 1) {
        Debug(Debug::ERROR) << "Cannot write dense index header " << outputName << "\n";
        EXIT(EXIT_FAILURE);
    }

    char *cursor = (char *) indexData.getData();
    const size_t indexDataSize = indexData.size();
    size_t consumed = 0;

    const size_t bufferCapacity = 65536;
    Entry *buffer = new Entry[bufferCapacity];
    size_t buffered = 0;

    uint64_t entryCount = 0;
    uint64_t firstKey = 0;
    uint64_t dataSize = 0;
    uint32_t maxSeqLen = 0;

    const char *cols[3];
    while (consumed < indexDataSize) {
        Util::getWordsOfLine(cursor, cols, 3);
        const DBKeyType key = Util::fast_atoi<DBKeyType>(cols[0]);
        const uint64_t offset = Util::fast_atoi<size_t>(cols[1]);
        const uint32_t length = Util::fast_atoi<size_t>(cols[2]);

        if (entryCount == 0) {
            firstKey = key;
        } else if (key != firstKey + entryCount) {
            // Row number and key must agree, otherwise a range load would return
            // the wrong sequences without any way to notice.
            Debug(Debug::ERROR) << "Database " << dbName << " does not have dense keys: expected key "
                                << (firstKey + entryCount) << " at row " << entryCount << ", found "
                                << key << ". A dense index requires keys numbered consecutively.\n";
            EXIT(EXIT_FAILURE);
        }

        buffer[buffered].offset = offset;
        buffer[buffered].length = length;
        buffered++;
        if (buffered == bufferCapacity) {
            if (fwrite(buffer, sizeof(Entry), buffered, out) != buffered) {
                Debug(Debug::ERROR) << "Cannot write dense index " << outputName << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffered = 0;
        }

        dataSize += length;
        maxSeqLen = std::max(maxSeqLen, length);
        entryCount++;

        cursor = Util::skipLine(cursor);
        consumed = cursor - (char *) indexData.getData();
    }

    if (buffered > 0 && fwrite(buffer, sizeof(Entry), buffered, out) != buffered) {
        Debug(Debug::ERROR) << "Cannot write dense index " << outputName << "\n";
        EXIT(EXIT_FAILURE);
    }
    delete[] buffer;
    indexData.close();

    header.magic = MAGIC;
    header.version = VERSION;
    header.entryCount = entryCount;
    header.firstKey = firstKey;
    header.dataSize = dataSize;
    header.maxSeqLen = maxSeqLen;
    if (fseek(out, 0, SEEK_SET) != 0 || fwrite(&header, sizeof(header), 1, out) != 1) {
        Debug(Debug::ERROR) << "Cannot finalise dense index " << outputName << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close dense index " << tmpName << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (rename(tmpName.c_str(), outputName.c_str()) != 0) {
        Debug(Debug::ERROR) << "Cannot publish dense index " << outputName << ": "
                            << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Dense index for " << dbName << ": " << entryCount << " entries, first key "
                       << firstKey << "\n";
}

DenseIndex::Info DenseIndex::readInfo(const std::string &dbName) {
    const std::string path = fileName(dbName);
    FILE *file = fopen(path.c_str(), "rb");
    if (file == NULL) {
        Debug(Debug::ERROR) << "Cannot open dense index " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }

    Header header;
    if (fread(&header, sizeof(header), 1, file) != 1) {
        Debug(Debug::ERROR) << "Cannot read dense index header " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(file);

    if (header.magic != MAGIC) {
        Debug(Debug::ERROR) << "File " << path << " is not a dense index\n";
        EXIT(EXIT_FAILURE);
    }
    if (header.version != VERSION) {
        Debug(Debug::ERROR) << "Dense index " << path << " has version " << header.version
                            << ", this build reads version " << VERSION << "\n";
        EXIT(EXIT_FAILURE);
    }

    Info info;
    info.entryCount = header.entryCount;
    info.firstKey = header.firstKey;
    info.dataSize = header.dataSize;
    info.maxSeqLen = header.maxSeqLen;
    return info;
}

DBReader<DBKeyType>::Index *DenseIndex::loadRange(const std::string &dbName, DBKeyType keyFrom,
                                                   DBKeyType keyTo, Info *rangeInfo) {
    const Info whole = readInfo(dbName);

    if (keyFrom < whole.firstKey || keyTo < keyFrom
            || keyTo > whole.firstKey + whole.entryCount) {
        Debug(Debug::ERROR) << "Key range [" << keyFrom << ", " << keyTo << ") is outside "
                            << dbName << ", which holds keys [" << whole.firstKey << ", "
                            << (whole.firstKey + whole.entryCount) << ")\n";
        EXIT(EXIT_FAILURE);
    }

    const size_t count = static_cast<size_t>(keyTo - keyFrom);
    const size_t rowFrom = static_cast<size_t>(keyFrom - whole.firstKey);

    const std::string path = fileName(dbName);
    FILE *file = fopen(path.c_str(), "rb");
    if (file == NULL) {
        Debug(Debug::ERROR) << "Cannot open dense index " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }

    DBReader<DBKeyType>::Index *index = new DBReader<DBKeyType>::Index[count == 0 ? 1 : count];

    uint64_t dataSize = 0;
    uint32_t maxSeqLen = 0;

    // Read in chunks rather than one record at a time so a large range costs a
    // handful of sequential reads instead of one syscall per key.
    const size_t bufferCapacity = 65536;
    Entry *buffer = new Entry[bufferCapacity];
    size_t done = 0;
    while (done < count) {
        const size_t batch = std::min(bufferCapacity, count - done);
        // fseeko/off_t, not fseek/long. A 32-bit long caps the seek at 2 GB, which
        // a dense index passes at 180M entries -- three orders of magnitude below
        // what this exists for -- and the overflow is silent.
        const off_t offset = static_cast<off_t>(sizeof(Header) + (rowFrom + done) * sizeof(Entry));
        if (fseeko(file, offset, SEEK_SET) != 0
                || fread(buffer, sizeof(Entry), batch, file) != batch) {
            Debug(Debug::ERROR) << "Cannot read dense index " << path << " at row "
                                << (rowFrom + done) << "\n";
            EXIT(EXIT_FAILURE);
        }
        for (size_t i = 0; i < batch; i++) {
            index[done + i].id = static_cast<DBKeyType>(keyFrom + done + i);
            index[done + i].offset = buffer[i].offset;
            index[done + i].length = buffer[i].length;
            dataSize += buffer[i].length;
            maxSeqLen = std::max(maxSeqLen, buffer[i].length);
        }
        done += batch;
    }
    delete[] buffer;
    fclose(file);

    if (rangeInfo != NULL) {
        rangeInfo->entryCount = count;
        rangeInfo->firstKey = keyFrom;
        rangeInfo->dataSize = dataSize;
        rangeInfo->maxSeqLen = maxSeqLen;
    }
    return index;
}

size_t DenseIndex::entryOffset(uint64_t row) {
    return sizeof(Header) + row * sizeof(Entry);
}

void DenseIndex::createEmpty(const std::string &dbName, uint64_t entryCount, uint64_t firstKey,
                             uint64_t dataSize, uint32_t maxSeqLen) {
    const std::string path = fileName(dbName);
    FILE *file = fopen(path.c_str(), "wb");
    if (file == NULL) {
        Debug(Debug::ERROR) << "Cannot create dense index " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }

    Header header;
    memset(&header, 0, sizeof(header));
    header.magic = MAGIC;
    header.version = VERSION;
    header.entryCount = entryCount;
    header.firstKey = firstKey;
    header.dataSize = dataSize;
    header.maxSeqLen = maxSeqLen;
    if (fwrite(&header, sizeof(header), 1, file) != 1) {
        Debug(Debug::ERROR) << "Cannot write dense index header " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fflush(file) != 0) {
        Debug(Debug::ERROR) << "Cannot flush dense index " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    // Size the file up front so the workers' pwrites land inside it rather than
    // each extending a sparse file from a different direction.
    if (ftruncate(fileno(file), static_cast<off_t>(entryOffset(entryCount))) != 0) {
        Debug(Debug::ERROR) << "Cannot size dense index " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close dense index " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
}

void DenseIndex::writeTextIndex(const std::string &dbName) {
    const Info info = readInfo(dbName);
    const std::string path = fileName(dbName);
    FILE *in = fopen(path.c_str(), "rb");
    if (in == NULL) {
        Debug(Debug::ERROR) << "Cannot open dense index " << path << ": " << strerror(errno) << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fseeko(in, static_cast<off_t>(entryOffset(0)), SEEK_SET) != 0) {
        Debug(Debug::ERROR) << "Cannot seek dense index " << path << "\n";
        EXIT(EXIT_FAILURE);
    }

    // Same publish-by-rename rule as build(). This one runs for hours at scale and
    // is the last thing a build does, so an interrupted run leaving a half-written
    // `.index` where a stock tool would read it is the likeliest way to get a
    // silently short database.
    const std::string textIndexName = dbName + ".index";
    const std::string textTmpName = textIndexName + ".tmp." + SSTR(getpid());
    FILE *out = FileUtil::openAndDelete(textTmpName.c_str(), "w");
    setvbuf(out, NULL, _IOFBF, 1024 * 1024 * 4);

    const size_t bufferCapacity = 65536;
    Entry *buffer = new Entry[bufferCapacity];
    char line[128];
    uint64_t done = 0;
    while (done < info.entryCount) {
        const size_t batch = std::min(bufferCapacity, static_cast<size_t>(info.entryCount - done));
        if (fread(buffer, sizeof(Entry), batch, in) != batch) {
            Debug(Debug::ERROR) << "Cannot read dense index " << path << " at row " << done << "\n";
            EXIT(EXIT_FAILURE);
        }
        for (size_t i = 0; i < batch; i++) {
            const int written = snprintf(line, sizeof(line), "%llu\t%llu\t%u\n",
                                         (unsigned long long) (info.firstKey + done + i),
                                         (unsigned long long) buffer[i].offset,
                                         buffer[i].length);
            if (fwrite(line, 1, static_cast<size_t>(written), out) != static_cast<size_t>(written)) {
                Debug(Debug::ERROR) << "Cannot write index " << textTmpName << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        done += batch;
    }
    delete[] buffer;
    fclose(in);
    if (fclose(out) != 0) {
        Debug(Debug::ERROR) << "Cannot close index " << textTmpName << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (rename(textTmpName.c_str(), textIndexName.c_str()) != 0) {
        Debug(Debug::ERROR) << "Cannot publish index " << textIndexName << ": " << strerror(errno)
                            << "\n";
        EXIT(EXIT_FAILURE);
    }
}
