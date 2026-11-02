// Exactness test of the gapless score-only search through the public Marv API, for every backend: deterministic
// inputs, every score compared with a CPU reference, and a digest of all results that every backend must reproduce.
//   test_gapless                  run all checks (exit code 0 = passed)
//   test_gapless --print-digest   print the result digest only
//   test_gapless --gpu-compatible only the padded, length-sorted, 21-letter cases the CUDA/HIP kernels accept
//   test_gapless --case N         run one case
#include "marv.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <stdexcept>
#include <vector>

// recorded from the Metal and SIMD backends; -DEXPECTED_DIGEST=0 disables the comparison
#ifndef EXPECTED_DIGEST
#define EXPECTED_DIGEST 0x491a6f2bcba8cbedULL
#endif

struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next() { uint64_t z = (s += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
    int uniform(int lo, int hi) { return lo + (int)(next() % (uint64_t)(hi - lo + 1)); }
};

struct Digest {
    uint64_t h = 0xcbf29ce484222325ULL;
    void add(uint64_t v) { for (int i = 0; i < 8; i++) { h ^= (v >> (8 * i)) & 0xff; h *= 0x100000001b3ULL; } }
};

struct Db {
    std::vector<char> data;
    std::vector<char> codes;
    std::vector<size_t> offsets;
    std::vector<int32_t> lengths;
    const char* ref() const { return codes.empty() ? data.data() : codes.data(); }
};

static unsigned char asciiTable[256];
static void initAsciiTable(int alphabet) { for (int b = 0; b < 256; b++) asciiTable[b] = (unsigned char)(alphabet - 1); for (int c = 0; c < alphabet; c++) asciiTable['A' + c] = (unsigned char) c; }

static Db makeDb(Rng& rng, size_t n, int minLen, int maxLen, int alphabet, bool aligned, bool sorted, bool ascii,
                 const std::vector<std::vector<char>>& extra = {}) {
    std::vector<std::vector<char>> seqs;
    for (size_t i = 0; i < n; i++) {
        std::vector<char> q(rng.uniform(minLen, maxLen));
        for (auto& c : q) c = (char) rng.uniform(0, alphabet - 1);
        seqs.push_back(q);
    }
    seqs.insert(seqs.end(), extra.begin(), extra.end());
    if (sorted) std::stable_sort(seqs.begin(), seqs.end(), [](const std::vector<char>& a, const std::vector<char>& b) { return a.size() < b.size(); });
    Db db;
    size_t pos = 0;
    for (auto& q : seqs) {
        db.offsets.push_back(pos); db.lengths.push_back((int32_t) q.size());
        db.data.insert(db.data.end(), q.begin(), q.end());
        pos += q.size();
        if (aligned) { while (pos % 4) { db.data.push_back((char)(alphabet - 1)); pos++; } }
    }
    db.offsets.push_back(pos);
    if (db.data.empty()) db.data.push_back(0);
    if (ascii) { db.codes = db.data; for (auto& c : db.data) c = (char)('A' + (unsigned char) c); }
    return db;
}

static std::vector<int8_t> makePssm(Rng& rng, int alphabet, int L, int lo, int hi, const std::vector<char>* query, int matchScore) {
    std::vector<int8_t> p((size_t) alphabet * L);
    for (int s = 0; s < alphabet; s++) for (int c = 0; c < L; c++) p[(size_t) s * L + c] = (int8_t) rng.uniform(lo, hi);
    if (query) for (int c = 0; c < L; c++) p[(size_t)(unsigned char)(*query)[c] * L + c] = (int8_t) matchScore;
    return p;
}

static int referenceScore(const int8_t* pssm, int alphabet, int L, const char* subject, int subjectLen) {
    std::vector<int> prev(L + 1, 0), cur(L + 1, 0);
    int best = 0;
    for (int j = 0; j < subjectLen; j++) {
        int letter = std::min((int)(unsigned char) subject[j], alphabet - 1);
        cur[0] = 0;
        for (int i = 1; i <= L; i++) {
            int v = prev[i - 1] + pssm[(size_t) letter * L + (i - 1)];
            cur[i] = v > 0 ? v : 0;
            if (cur[i] > best) best = cur[i];
        }
        std::swap(prev, cur);
    }
    return best;
}

static int failures = 0;
static bool translationSupported = true;
static bool gpuCompatible = false;
static int selectedCase = -1;
static Digest digest;
static void check(bool ok, const std::string& what) { if (!ok) { failures++; fprintf(stderr, "FAIL: %s\n", what.c_str()); } }

struct Case {
    size_t numSubjects; int minLen, maxLen, alphabet; bool aligned, sorted;
    std::vector<int> queryLens; bool plantQueries; int topN; bool subset;
    int matchScore = 11;
    size_t budget = 0; // database memory budget forcing the streaming path, 0 = unlimited
};

static std::vector<Marv::Result> runCase(uint64_t seed, const Case& cs, bool int8, bool ascii, const std::string& tag) {
    Rng rng(seed);
    std::vector<std::vector<char>> queries; std::vector<std::vector<int8_t>> pssms;
    for (int L : cs.queryLens) {
        std::vector<char> q(L); for (auto& c : q) c = (char) rng.uniform(0, cs.alphabet - 1);
        pssms.push_back(makePssm(rng, cs.alphabet, L, -6, 11, cs.plantQueries ? &q : nullptr, cs.matchScore));
        queries.push_back(q);
    }
    Db db = makeDb(rng, cs.numSubjects, cs.minLen, cs.maxLen, cs.alphabet, cs.aligned, cs.sorted, ascii, cs.plantQueries ? queries : std::vector<std::vector<char>>{});
    const size_t numSubjects = db.lengths.size();
    Marv::MemoryConfig memory;
    if (cs.budget) { memory.maxGpuMem = cs.budget; memory.maxBatchBytes = 64 * 1024; memory.maxBatchSequences = 256; }
    Marv marv(numSubjects, cs.alphabet, memory);
    void* h = marv.loadDb(db.data.data(), db.offsets.data(), db.lengths.data(), db.data.size());
    marv.setDb(h);
    if (ascii) { initAsciiTable(cs.alphabet); marv.setLetterTranslation(asciiTable); }

    std::vector<const char*> qp; std::vector<int> ql; std::vector<int8_t*> pp;
    for (size_t i = 0; i < queries.size(); i++) { qp.push_back(queries[i].data()); ql.push_back((int) queries[i].size()); pp.push_back(pssms[i].data()); }
    std::vector<int32_t> subsetIds;
    if (cs.subset) { for (size_t s = 0; s < numSubjects; s++) if (rng.uniform(0, 2) == 0) subsetIds.push_back((int32_t) s); }
    const size_t universe = cs.subset ? subsetIds.size() : numSubjects;
    const int N = (int) std::min<size_t>(cs.topN > 0 ? (size_t) cs.topN : universe, universe);

    std::vector<Marv::Result> res((size_t) queries.size() * std::max(N, 1), Marv::Result(0, 0, 0, 0));
    auto stats = marv.search_scoreOnly(res.data(), qp, ql, pp, N, Marv::GAPLESS, 0, 0, int8, true, cs.subset ? &subsetIds : nullptr);
    check(stats.size() == queries.size(), tag + "stats size");
    std::vector<Marv::Result> again(res.size(), Marv::Result(0, 0, 0, 0)); // second search reuses the retained workspace
    marv.search_scoreOnly(again.data(), qp, ql, pp, N, Marv::GAPLESS, 0, 0, int8, true, cs.subset ? &subsetIds : nullptr);
    bool same = true;
    for (size_t i = 0; same && i < res.size(); i++) same = again[i].id == res[i].id && again[i].score == res[i].score;
    check(same, tag + "repeated search with the retained workspace gives different results");

    for (size_t qi = 0; qi < queries.size(); qi++) {
        const int L = (int) queries[qi].size();
        std::vector<std::pair<int, uint32_t>> ref;
        auto score = [&](size_t s) { return referenceScore(pssms[qi].data(), cs.alphabet, L, db.ref() + db.offsets[s], db.lengths[s]); };
        if (cs.subset) for (int32_t s : subsetIds) ref.emplace_back(score(s), (uint32_t) s);
        else for (size_t s = 0; s < numSubjects; s++) ref.emplace_back(score(s), (uint32_t) s);
        std::sort(ref.begin(), ref.end(), [](auto& a, auto& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
        const std::string t = tag + " L=" + std::to_string(L) + ": ";
        check(stats[qi].results == (size_t) N, t + "result count " + std::to_string(stats[qi].results) + " != " + std::to_string(N));
        for (int i = 0; i < N; i++) {
            const Marv::Result& r = res[(size_t) qi * N + i];
            if (r.id != ref[i].second || r.score != ref[i].first) {
                check(false, t + "rank " + std::to_string(i) + " got (id " + std::to_string(r.id) + ", score " + std::to_string(r.score) + ") expected (id " + std::to_string(ref[i].second) + ", score " + std::to_string(ref[i].first) + ")");
                break;
            }
            check(r.qEndPos == -1 && r.dbEndPos == -1, t + "end positions must be -1");
        }
    }
    return res;
}

static void runAll() {
    const std::vector<int> tiledLens{1, 3, 7, 16, 79, 100, 191, 192, 193, 300, 511, 512, 513, 800, 1000};
    const std::vector<int> untiledLens{1537, 2500, 4096, 4097, 6000, 8192};
    std::vector<Case> cases;
    for (int alphabet : {21, 25}) {
        for (bool aligned : {true, false}) {
            cases.push_back({3000, 1, 400, alphabet, aligned, true, tiledLens, false, 0, false});
            cases.push_back({3000, 1, 400, alphabet, aligned, true, tiledLens, false, 50, false});
            cases.push_back({2000, 1, 400, alphabet, aligned, true, tiledLens, false, 40, true});
        }
        cases.push_back({400, 1, 3000, alphabet, true, true, untiledLens, false, 0, false});
        cases.push_back({300, 5, 300, alphabet, true, true, {201, 401, 1001}, true, 0, false});
        cases.push_back({300, 5, 300, alphabet, true, false, {201, 401, 1001}, true, 0, false});
        cases.push_back({300, 5, 300, alphabet, true, true, {3001}, true, 0, false});
    }
    cases.push_back({500, 0, 3, 21, true, true, {1, 20, 79}, false, 0, false});
    if (selectedCase >= int(cases.size())) { check(false, "case index is out of range"); return; }

    for (size_t ci = 0; ci < cases.size(); ci++) {
        const Case& cs = cases[ci];
        // CUDA/HIP kernels consume the padded, length-sorted 21-letter GPU database format.
        // Keep case indices/seeds unchanged so both backends can check identical inputs.
        if (selectedCase >= 0 && ci != size_t(selectedCase)) continue;
        if (gpuCompatible && (cs.alphabet != 21 || !cs.aligned || !cs.sorted)) continue;
        const uint64_t seed = 12345 + 7919 * ci;
        for (bool int8 : {false, true}) {
            if (int8 && cs.queryLens.size() == 1 && cs.queryLens[0] == 3001) continue;
            const std::string tag = "case " + std::to_string(ci) + (int8 ? " int8" : " int16");
            std::vector<Marv::Result> numeric = runCase(seed, cs, int8, false, tag);
            for (const Marv::Result& r : numeric) { digest.add(r.id); digest.add((uint64_t)(int64_t) r.score); }
            if (translationSupported) {
                try {
                    std::vector<Marv::Result> translated = runCase(seed, cs, int8, true, tag + " ascii");
                    bool same = translated.size() == numeric.size();
                    for (size_t i = 0; same && i < numeric.size(); i++) same = translated[i].id == numeric[i].id && translated[i].score == numeric[i].score;
                    check(same, tag + ": ASCII-encoded database gives different results than the numeric one");
                } catch (const std::runtime_error& e) {
                    translationSupported = false;
                    fprintf(stderr, "note: letter translation not available on this backend (%s); ASCII cases skipped\n", e.what());
                }
            }
        }
    }
    if (selectedCase < 0) {
        // Regression checks outside the historical digest: 8-bit retries may also overflow,
        // signed packed addition can wrap, and a long FP16 alignment can reach infinity.
        const Case highScores{0, 1, 1, 21, true, true, {260, 3001, 4096}, true, 0, false, 127};
        for (bool int8 : {false, true})
            runCase(98765, highScores, int8, false, int8 ? "wide retry int8" : "wide retry int16");
        // streamed database (small budget, partial cache) with a subset that is gathered from host memory
        const Case streamed{2048, 100, 300, 21, true, true, {48, 79, 161}, true, 50, true, 11, 384 * 1024};
        for (bool int8 : {false, true}) runCase(4242, streamed, int8, false, int8 ? "streamed int8" : "streamed int16");
        // queries longer than one 8192-column pass take the multi-tile scan; planted full-length hits cross every tile
        const Case multiTile{100, 1, 1000, 21, true, true, {8193, 12000, 20000}, true, 0, false};
        for (bool int8 : {false, true}) runCase(5150, multiTile, int8, false, int8 ? "multi-tile int8" : "multi-tile int16");
        const Case multiTileSubset{300, 1, 1500, 25, false, false, {9000}, false, 40, true};
        if (!gpuCompatible) runCase(5151, multiTileSubset, false, false, "multi-tile subset");
        if (translationSupported) runCase(5151, multiTileSubset, false, true, "multi-tile subset ascii");
        // a 16385-column hit far above 16 bits forces the 32-bit multi-tile retry
        const Case multiTileWide{0, 1, 1, 21, true, true, {16385}, true, 0, false, 127};
        for (bool int8 : {false, true}) runCase(5152, multiTileWide, int8, false, int8 ? "multi-tile wide int8" : "multi-tile wide int16");
    }
}

int main(int argc, char** argv) {
    bool printOnly = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--print-digest") printOnly = true;
        else if (arg == "--gpu-compatible") { gpuCompatible = true; translationSupported = false; }
        else if (arg == "--case" && i + 1 < argc) selectedCase = std::stoi(argv[++i]);
        else { fprintf(stderr, "usage: test_gapless [--print-digest] [--gpu-compatible] [--case N]\n"); return 2; }
    }
    if (Marv::getDeviceIds().empty()) { fprintf(stderr, "no device\n"); return 2; }
    runAll();
    if (printOnly) { printf("digest 0x%016llx\n", (unsigned long long) digest.h); return failures ? 1 : 0; }
    const uint64_t expected = selectedCase >= 0 || EXPECTED_DIGEST == 0 ? 0 :
        gpuCompatible ? 0xcba3a00c27ffe2dfULL : EXPECTED_DIGEST;
    if (expected != 0 && digest.h != expected) {
        char buf[160]; snprintf(buf, sizeof buf, "result digest 0x%016llx differs from the recorded 0x%016llx: results differ from the other backends", (unsigned long long) digest.h, (unsigned long long) expected); check(false, buf);
    }
    if (failures == 0) { printf("libmarv gapless test: all checks passed, digest 0x%016llx%s\n", (unsigned long long) digest.h, translationSupported ? "" : " (letter translation not tested)"); return 0; }
    fprintf(stderr, "libmarv gapless test: %d failures\n", failures);
    return 1;
}
