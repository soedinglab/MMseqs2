#ifndef MMSEQS_DENSEINDEX_H
#define MMSEQS_DENSEINDEX_H

#include "DBReader.h"
#include "IndexTypes.h"

#include <cstdint>
#include <string>

// A fixed-width companion index that makes a *key range* of a database loadable
// without touching the rest of it.
//
// MMseqs2's native .index is text with variable-length lines, so locating the
// entries for keys [from, to) means parsing every preceding line, and DBReader
// materialises the whole thing as a 24 byte Index record per entry. At 1e12
// sequences that is ~34 TB of text parsed into ~24 TB of RAM -- per worker. No
// node can hold it, and every distributed stage would pay it just to read its own
// slice.
//
// This companion file stores one fixed-width record per key, so a key range maps
// straight to a byte range and a worker reads only what it needs. Loading a range
// costs 12 bytes per key in that range and nothing for the rest of the database.
//
// The scheme requires keys to be *dense*: key == firstKey + row. That is what the
// length-ranked createdb produces, and it is also why the length ranking matters
// beyond ordering -- it makes the key its own array offset. Databases with sparse
// keys are unaffected and keep using the text index.
class DenseIndex {
public:
    // Packed so the file is 12 bytes per entry rather than 16. At 1e12 sequences
    // that difference is 4 TB of scratch, and the reads are copied into aligned
    // structs anyway so unaligned access never happens.
    struct __attribute__((__packed__)) Entry {
        uint64_t offset;
        uint32_t length;
    };

    struct Info {
        uint64_t entryCount;
        uint64_t firstKey;
        uint64_t dataSize;
        uint32_t maxSeqLen;
    };

    static std::string fileName(const std::string &dbName);
    static bool exists(const std::string &dbName);

    // Streams the database's text index and writes the companion file. Runs in
    // constant memory -- it never materialises the index -- so it is usable at any
    // scale, including as a one-off converter for a database built by stock
    // createdb.
    //
    // Exits if the keys are not dense. Falling back silently would leave callers
    // reading a companion file whose row numbers do not mean what they think.
    static void build(const std::string &dbName);

    // Creates the companion file for a database that many workers are about to
    // write concurrently. The header is final from the start (the planner knows
    // the totals before any sequence is written) and the rows are left zeroed;
    // each worker then fills only its own rows with pwrite at entryOffset(row).
    // Rows are fixed width, so those writes never overlap and need no locking.
    static void createEmpty(const std::string &dbName, uint64_t entryCount, uint64_t firstKey,
                            uint64_t dataSize, uint32_t maxSeqLen);

    // Byte offset of a row within the companion file.
    static size_t entryOffset(uint64_t row);

    // Writes the text .index from the companion file, streaming in constant
    // memory. The inverse of build(), for the stock tools that still read the
    // text index; the distributed stages read the companion file directly.
    static void writeTextIndex(const std::string &dbName);

    static Info readInfo(const std::string &dbName);

    // Loads entries for keys [keyFrom, keyTo) into a newly allocated array that
    // the caller owns. The array is laid out exactly as DBReader expects, so the
    // caller can hand it to DBReader's external-index constructor together with
    // setDataFile()/open() and get a reader scoped to that key range:
    //
    //   DBReader<DBKeyType>::Index *idx = DenseIndex::loadRange(db, from, to, &info);
    //   DBReader<DBKeyType> reader(idx, info.entryCount, info.dataSize,
    //                              from + info.entryCount - 1, dbType,
    //                              info.maxSeqLen, threads);
    //   reader.setDataFile(db.c_str());
    //   reader.open(DBReader<DBKeyType>::NOSORT);
    //
    // rangeInfo receives the entry count, data size and max length *of the range*,
    // which are the values that constructor needs.
    static DBReader<DBKeyType>::Index *loadRange(const std::string &dbName,
                                                 DBKeyType keyFrom, DBKeyType keyTo,
                                                 Info *rangeInfo);

private:
    struct __attribute__((__packed__)) Header {
        uint64_t magic;
        uint64_t version;
        uint64_t entryCount;
        uint64_t firstKey;
        uint64_t dataSize;
        uint32_t maxSeqLen;
        uint32_t reserved;
    };

    static const uint64_t MAGIC = 0x4d4d44454e534931ULL;  // "MMDENSI1"
    static const uint64_t VERSION = 1;
};

#endif
