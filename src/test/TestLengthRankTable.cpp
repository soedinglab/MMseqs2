// Tests for the key -> length table that lets the k-mer record drop its seqLen
// field.
//
// The property under test is simple to state and expensive to get wrong: for
// every key in the database, the table must return exactly the length the
// length-ranked key assignment gave that key. A table that is merely *close* does
// not fail visibly -- it mis-ranks the centre of a k-mer group, which produces a
// different clustering with no error message. So the round-trip here is checked
// exhaustively against a brute-force expectation rather than sampled.

#include "LengthRankTable.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <functional>

#include <sys/wait.h>
#include <unistd.h>

const char* binary_name = "test_lengthranktable";

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
    char tmpl[] = "/tmp/mmseqs_lenrank_testXXXXXX";
    char *dir = mkdtemp(tmpl);
    if (dir == NULL) {
        perror("mkdtemp");
        exit(EXIT_FAILURE);
    }
    return std::string(dir);
}

static uint64_t nextRandom(uint64_t &state) {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

// Turns a multiset of sequence lengths into the runs the planner would produce:
// keys handed out longest first, so a length's run starts at the number of
// strictly longer sequences.
static std::vector<LengthRankTable::Run> runsFor(std::vector<uint64_t> lengths) {
    std::sort(lengths.begin(), lengths.end(), std::greater<uint64_t>());
    std::vector<LengthRankTable::Run> runs;
    size_t i = 0;
    while (i < lengths.size()) {
        size_t j = i;
        while (j < lengths.size() && lengths[j] == lengths[i]) {
            j++;
        }
        LengthRankTable::Run run;
        run.length = lengths[i];
        run.firstKey = i;
        run.count = j - i;
        runs.push_back(run);
        i = j;
    }
    return runs;
}

// The brute-force expectation: sort descending, and key k has the length at
// position k. This is the definition of a length-ranked key.
static void checkRoundTrip(const std::string &dir, const std::string &name,
                           std::vector<uint64_t> lengths, const std::string &what) {
    const std::string db = dir + "/" + name;
    std::vector<uint64_t> sorted = lengths;
    std::sort(sorted.begin(), sorted.end(), std::greater<uint64_t>());

    LengthRankTable::write(db, runsFor(lengths), sorted.size());
    LengthRankTable table;
    table.open(db);

    bool ok = table.getEntryCount() == sorted.size();
    ok = ok && table.getMaxLength() == (sorted.empty() ? 0 : sorted[0]);
    for (size_t k = 0; k < sorted.size(); k++) {
        unsigned int got = 0;
        if (table.tryLengthOf(k, got) == false || got != sorted[k]) {
            ok = false;
            break;
        }
    }
    // One past the end must be refused, not clamped: a key outside the database
    // is a caller bug, and returning the last length would hide it.
    unsigned int ignored = 0;
    ok = ok && table.tryLengthOf(sorted.size(), ignored) == false;
    check(ok, what);
}

static void testRoundTrips(const std::string &dir) {
    // One sequence.
    checkRoundTrip(dir, "single", std::vector<uint64_t>(1, 42), "one sequence round-trips");

    // Every sequence the same length: one run covering the whole key space.
    checkRoundTrip(dir, "uniform", std::vector<uint64_t>(1000, 300),
                   "1000 sequences of one length round-trip");

    // Every sequence a distinct length: one run per key, the worst case for the
    // run count.
    std::vector<uint64_t> distinct;
    for (uint64_t l = 1; l <= 2000; l++) {
        distinct.push_back(l);
    }
    checkRoundTrip(dir, "distinct", distinct, "2000 distinct lengths round-trip");

    // A realistic protein-ish distribution: heavy repetition, wide gaps, a long
    // tail. Gaps matter because a length with no sequences must never be returned.
    std::vector<uint64_t> realistic;
    uint64_t state = 0x243F6A8885A308D3ULL;
    for (int i = 0; i < 20000; i++) {
        const uint64_t bucket = nextRandom(state) % 100;
        if (bucket < 70) {
            realistic.push_back(100 + nextRandom(state) % 300);
        } else if (bucket < 95) {
            realistic.push_back(400 + nextRandom(state) % 2000);
        } else {
            realistic.push_back(3000 + nextRandom(state) % 60000);
        }
    }
    checkRoundTrip(dir, "realistic", realistic,
                   "20000 sequences over a skewed length distribution round-trip");

    // Length 1 present, and a maximum at the 16-bit ceiling the k-mer record
    // imposes.
    std::vector<uint64_t> extremes;
    extremes.push_back(1);
    extremes.push_back(1);
    extremes.push_back(LengthRankTable::MAX_SEQUENCE_LENGTH);
    extremes.push_back(2);
    checkRoundTrip(dir, "extremes", extremes, "length 1 and the 65535 ceiling round-trip");
}

static void testCorruptionRejected(const std::string &dir) {
    const std::string db = dir + "/corrupt";
    std::vector<uint64_t> lengths;
    for (uint64_t l = 1; l <= 50; l++) {
        lengths.push_back(l);
        lengths.push_back(l);
    }
    LengthRankTable::write(db, runsFor(lengths), lengths.size());
    const std::string path = LengthRankTable::fileName(db);

    // Baseline: the untouched file opens.
    {
        LengthRankTable table;
        table.open(db);
        check(table.getRunCount() == 50, "a table of 50 distinct lengths has 50 runs");
    }

    std::vector<char> good;
    {
        FILE *f = fopen(path.c_str(), "rb");
        char buffer[4096];
        size_t got = 0;
        while ((got = fread(buffer, 1, sizeof(buffer), f)) > 0) {
            good.insert(good.end(), buffer, buffer + got);
        }
        fclose(f);
    }

    // Each of these must be rejected. They are checked by running a child
    // process, because the reader's contract is to exit rather than return an
    // error -- a silently wrong length table is the one outcome that must be
    // impossible.
    struct Case {
        const char *what;
        int kind;  // 0 truncate, 1 magic, 2 trailing
    };
    const Case cases[] = {
        {"a truncated length-rank table is rejected", 0},
        {"a bad magic is rejected", 1},
        {"trailing bytes are rejected", 2},
    };
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
        std::vector<char> broken = good;
        if (cases[c].kind == 0) {
            broken.resize(broken.size() - 8);
        } else if (cases[c].kind == 1) {
            broken[0] = static_cast<char>(broken[0] ^ 0xFF);
        } else {
            broken.push_back('x');
        }
        FILE *f = fopen(path.c_str(), "wb");
        fwrite(broken.data(), 1, broken.size(), f);
        fclose(f);

        fflush(stdout);
        fflush(stderr);
        const pid_t pid = fork();
        if (pid == 0) {
            // Silence the expected error message so the log stays readable.
            freopen("/dev/null", "w", stderr);
            LengthRankTable table;
            table.open(db);
            _exit(0);  // reached only if the corruption was accepted
        }
        int status = 0;
        waitpid(pid, &status, 0);
        const bool rejected = !(WIFEXITED(status) && WEXITSTATUS(status) == 0);
        check(rejected, cases[c].what);
    }
}

int main(int, const char**) {
    const std::string dir = makeTempDir();

    testRoundTrips(dir);
    testCorruptionRejected(dir);

    // Best-effort cleanup; the test owns a fresh mkdtemp directory.
    std::string cmd = "rm -rf '" + dir + "'";
    if (system(cmd.c_str()) != 0) {
        fprintf(stderr, "warning: could not remove %s\n", dir.c_str());
    }

    if (failures > 0) {
        fprintf(stderr, "\n%d check(s) failed\n", failures);
        return EXIT_FAILURE;
    }
    fprintf(stdout, "\nall checks passed\n");
    return EXIT_SUCCESS;
}
