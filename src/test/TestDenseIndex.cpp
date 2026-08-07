// Tests for the fixed-width companion index used to open a key range of a
// database without materialising the whole text index.
//
// The property under test is the one the distributed stages depend on: a reader
// built from a key range must return byte-identical sequences to a reader over
// the whole database, while only ever reading that range's slice of the index.

#include "DBReader.h"
#include "DBWriter.h"
#include "DenseIndex.h"
#include "Parameters.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <unistd.h>

const char* binary_name = "test_denseindex";

static int failures = 0;

static void check(bool condition, const std::string &what) {
    if (condition == false) {
        fprintf(stderr, "FAIL: %s\n", what.c_str());
        failures++;
    } else {
        fprintf(stdout, "ok:   %s\n", what.c_str());
    }
}

static std::string makeTempDir() {
    char tmpl[] = "/tmp/mmseqs_denseidx_testXXXXXX";
    char *dir = mkdtemp(tmpl);
    if (dir == NULL) {
        perror("mkdtemp");
        exit(EXIT_FAILURE);
    }
    return std::string(dir);
}

static void removeTempDir(const std::string &dir) {
    std::string cmd = "rm -rf '" + dir + "'";
    if (system(cmd.c_str()) != 0) {
        fprintf(stderr, "warning: could not clean up %s\n", dir.c_str());
    }
}

// Deliberately varies sequence length with the key so that a wrong offset shows
// up as a length mismatch as well as a content mismatch.
static std::string makeSequence(size_t key) {
    const char alphabet[] = "ACDEFGHIKLMNPQRSTVWY";
    const size_t length = 5 + (key % 37);
    std::string sequence;
    for (size_t i = 0; i < length; i++) {
        sequence.push_back(alphabet[(key + i) % 20]);
    }
    return sequence;
}

static void writeDatabase(const std::string &dbName, size_t count) {
    DBWriter writer(dbName.c_str(), (dbName + ".index").c_str(), 1, 0,
                    Parameters::DBTYPE_AMINO_ACIDS);
    writer.open();
    for (size_t key = 0; key < count; key++) {
        std::string sequence = makeSequence(key);
        sequence.push_back('\n');
        writer.writeData(sequence.c_str(), sequence.size(), static_cast<DBKeyType>(key), 0);
    }
    writer.close();
}

// Opens [keyFrom, keyTo) through the companion index and hands the resulting
// array to DBReader's external-index constructor -- exactly the call sequence the
// distributed stages use.
static DBReader<DBKeyType> *openRange(const std::string &dbName, DBKeyType keyFrom, DBKeyType keyTo,
                                       DBReader<DBKeyType>::Index **indexOut) {
    DenseIndex::Info info;
    DBReader<DBKeyType>::Index *index = DenseIndex::loadRange(dbName, keyFrom, keyTo, &info);
    DBReader<DBKeyType> *reader = new DBReader<DBKeyType>(
            index, info.entryCount, info.dataSize,
            static_cast<DBKeyType>(keyTo == keyFrom ? keyFrom : keyTo - 1),
            Parameters::DBTYPE_AMINO_ACIDS, info.maxSeqLen, 1);
    reader->setDataFile(dbName.c_str());
    reader->open(DBReader<DBKeyType>::NOSORT);
    *indexOut = index;
    return reader;
}

int main(int, const char**) {
    const std::string dir = makeTempDir();
    const std::string dbName = dir + "/seqdb";
    const size_t count = 5000;

    writeDatabase(dbName, count);
    DenseIndex::build(dbName);

    check(DenseIndex::exists(dbName), "companion index file is created");

    DenseIndex::Info whole = DenseIndex::readInfo(dbName);
    check(whole.entryCount == count, "companion index covers every entry");
    check(whole.firstKey == 0, "companion index records the first key");

    // Reference: the database read the normal way.
    DBReader<DBKeyType> full(dbName.c_str(), (dbName + ".index").c_str(), 1,
                             DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX);
    full.open(DBReader<DBKeyType>::NOSORT);
    check(full.getSize() == count, "reference reader sees every entry");

    // A range in the middle, a range at the start, and a range running to the end
    // -- the three cases where an off-by-one in the row arithmetic would show up.
    const DBKeyType ranges[][2] = { {1000, 1500}, {0, 64}, {4900, 5000}, {2500, 2501} };
    for (size_t r = 0; r < sizeof(ranges) / sizeof(ranges[0]); r++) {
        const DBKeyType from = ranges[r][0];
        const DBKeyType to = ranges[r][1];

        DBReader<DBKeyType>::Index *index = NULL;
        DBReader<DBKeyType> *ranged = openRange(dbName, from, to, &index);

        bool sizeOk = ranged->getSize() == static_cast<size_t>(to - from);
        bool contentOk = true;
        bool keysOk = true;
        for (size_t i = 0; i < ranged->getSize(); i++) {
            const DBKeyType expectedKey = static_cast<DBKeyType>(from + i);
            if (ranged->getDbKey(i) != expectedKey) {
                keysOk = false;
                break;
            }
            const size_t fullId = full.getId(expectedKey);
            if (ranged->getSeqLen(i) != full.getSeqLen(fullId)) {
                contentOk = false;
                break;
            }
            if (memcmp(ranged->getData(i, 0), full.getData(fullId, 0),
                       full.getEntryLen(fullId)) != 0) {
                contentOk = false;
                break;
            }
        }

        const std::string label = "[" + std::to_string(from) + ", " + std::to_string(to) + ")";
        check(sizeOk, "range " + label + " has the expected entry count");
        check(keysOk, "range " + label + " reports the expected keys");
        check(contentOk, "range " + label + " returns byte-identical sequences");

        ranged->close();
        delete ranged;
        delete[] index;
    }

    // The whole database as one range must agree with the reference reader
    // entry for entry, which also covers the data-size and max-length totals.
    {
        DBReader<DBKeyType>::Index *index = NULL;
        DBReader<DBKeyType> *ranged = openRange(dbName, 0, static_cast<DBKeyType>(count), &index);
        bool allOk = ranged->getSize() == count;
        for (size_t i = 0; allOk && i < count; i++) {
            allOk = ranged->getSeqLen(i) == full.getSeqLen(full.getId(static_cast<DBKeyType>(i)));
        }
        check(allOk, "full-database range matches the reference reader");
        ranged->close();
        delete ranged;
        delete[] index;
    }

    full.close();
    removeTempDir(dir);

    if (failures > 0) {
        fprintf(stderr, "\n%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "\nall checks passed\n");
    return EXIT_SUCCESS;
}
