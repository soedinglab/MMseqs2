#include "Parameters.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "MemoryMapped.h"
#include "FastSort.h"

#include <dirent.h>

#ifdef OPENMP
#include <omp.h>
#endif

static uint32_t be32(const unsigned char *p) {
    uint32_t v = ((uint32_t)p[0] << 24) | ((uint32_t)p[1] << 16) | ((uint32_t)p[2] << 8) | (uint32_t)p[3];
    return v;
}

static inline uint64_t le64(const unsigned char *p) {
    return ((uint64_t)p[0]) |
           ((uint64_t)p[1] << 8) |
           ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) |
           ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) |
           ((uint64_t)p[7] << 56);
}

static bool readBe32String(const unsigned char *idx, size_t idxSize, size_t &pos) {
    if (pos + 4 > idxSize) {
        return false;
    }

    uint32_t len = be32(idx + pos);
    pos += 4;
    if (pos + (size_t)len > idxSize) {
        return false;
    }
    pos += (size_t)len;
    return true;
}



static inline char iupacFromMask(unsigned char m) {
    static const char LUT[16] = {'N', 'A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
    if (m < 16) {
        return LUT[m];
    }
    return 'N';
}

// 2-bit (NCBI2na) -> DNA base
static inline char baseFrom2na(unsigned char code) {
    static const char LUT[4] = {'A', 'C', 'G', 'T'};

    return LUT[code & 3];
}

// compute number of residues between [startByte, endByte) in .nsq 2na block
static size_t countNucResidues(const unsigned char *nsq, uint32_t startByte, uint32_t endByte) {
    if (endByte <= startByte) {
        return 0;
    }

    size_t bytes = (size_t)(endByte - startByte);
    unsigned char last = nsq[endByte - 1];
    size_t lastCount = (size_t)(last & 0x03);
    return (bytes - 1U) * 4U + lastCount;
}

// decode NCBI2na from [startByte, endByte)
static std::string decodeNucSlice(const unsigned char *nsq, uint32_t startByte, uint32_t endByte) {
    size_t nres = countNucResidues(nsq, startByte, endByte);

    std::string out;
    out.resize(nres);
    for (size_t r = 0; r < nres; ++r) {
        size_t byteIdx = (size_t)startByte + (r >> 2);
        unsigned char b = nsq[byteIdx];
        unsigned int shift = 6U - 2U * (unsigned int)(r & 3U);
        unsigned char code = (unsigned char)((b >> shift) & 0x03U);
        out[r] = baseFrom2na(code);
    }

    return out;
}

// Apply ambiguity table to an ASCII sequence
// ambStart points at the table,
// ambEnd is the start of the next sequence (.nin: S[i+1]).
static void applyAmbiguityPatches(
    std::string &seqStr,
    const unsigned char *nsq,
    uint32_t ambStart,
    uint32_t ambEnd
) {
    if (ambStart >= ambEnd) {
        return;
    }

    if ((size_t)ambEnd - (size_t)ambStart < 4U) {
        return;
    }

    uint32_t countWords = be32(nsq + ambStart);
    bool is64 = ((countWords & 0x80000000U) != 0);
    countWords &= 0x7FFFFFFFU;

    size_t tableBytes = (size_t)countWords * 4U;
    size_t avail = (size_t)ambEnd - (size_t)ambStart;
    if (4U + tableBytes > avail) {
        return;
    }

    size_t numEntries = is64 ? (tableBytes / 8U) : (tableBytes / 4U);
    size_t p = (size_t)ambStart + 4U;
    for (size_t e = 0; e < numEntries; ++e) {
        if (!is64) {
            uint32_t w = be32(nsq + p);
            p += 4U;
            unsigned int sym = (unsigned int)(w >> 28);
            unsigned int rep = ((unsigned int)((w >> 24) & 0x0FU)) + 1U;
            unsigned int off = (unsigned int)(w & 0x00FFFFFFU);
            char ch = iupacFromMask((unsigned char)sym);
            size_t a = (size_t)off;
            size_t b = a + (size_t)rep;

            if (a >= seqStr.size()) {
                continue;
            }

            if (b > seqStr.size()) {
                b = seqStr.size();
            }

            for (size_t i = a; i < b; ++i) {
                seqStr[i] = ch;
            }
        } else {
            uint32_t hi = be32(nsq + p);
            uint32_t lo = be32(nsq + p + 4U);
            p += 8U;
            uint64_t word = (((uint64_t)hi) << 32) | (uint64_t)lo;
            unsigned int sym = (unsigned int)((word >> 60) & 0x0FULL);
            unsigned int rep = (unsigned int)(((word >> 48) & 0x0FFFULL) + 1ULL);
            uint64_t off = (word & 0x0000FFFFFFFFFFFFULL);
            char ch = iupacFromMask((unsigned char)sym);

            size_t a = (size_t)off;
            size_t b = a + (size_t)rep;
            if (a >= seqStr.size()) {
                continue;
            }

            if (b > seqStr.size()) {
                b = seqStr.size();
            }

            for (size_t i = a; i < b; ++i) {
                seqStr[i] = ch;
            }
        }
    }
}

static const char STDAA[28] = {
    '-', 'A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','U','*','O','J'
};
static std::string decodePsqSlice(const unsigned char *psq, uint32_t s, uint32_t e) {
    if (e > s && psq[e - 1] == 0) {
        --e;
    }

    std::string out;
    out.resize(e - s);
    for (uint32_t i = s; i < e; ++i) {
        unsigned char b = psq[i];
        out[i - s] = (b < 28) ? STDAA[b] : 'X';
    }

    return out;
}

static void pushNumberedVols(
    const std::string &base,
    const char *ext,
    std::vector<std::string> &out) {
    std::string dir = FileUtil::dirName(base);

    size_t pos = base.find_last_of("\\/");
    std::string leaf = (pos == std::string::npos) ? base : base.substr(pos + 1);

    std::string scanDir = dir;
    if (scanDir.empty()) {
        scanDir = ".";
    }

    std::vector< std::pair<int, std::string> > numbered;
    bool hasPlain = false;

    DIR* dp = opendir(scanDir.c_str());
    if (dp != NULL) {
        struct dirent* de = NULL;
        while ((de = readdir(dp)) != NULL) {
#if defined(_DIRENT_HAVE_D_TYPE)
            if (de->d_type == DT_DIR) {
                continue;
            }
#endif
            const char *name = de->d_name;
            size_t len = std::strlen(name);
            size_t extLen = std::strlen(ext);

            if (len < extLen) {
                continue;
            }

            if (std::strcmp(name + (len - extLen), ext) != 0) {
                continue;
            }

            std::string fname(name, name + len);
            if (fname == leaf + ext) {
                hasPlain = true;
            } else {
                if (fname.size() > leaf.size() + 1 + extLen) {
                    bool prefixOk = (fname.compare(0, leaf.size(), leaf) == 0);
                    bool dotOk = (fname[leaf.size()] == '.');

                    if (prefixOk && dotOk) {
                        size_t start = leaf.size() + 1;
                        size_t end = fname.size() - extLen;
                        std::string suf = fname.substr(start, end - start);

                        bool allDigits = !suf.empty();
                        for (size_t i = 0; i < suf.size(); ++i) {
                            if (suf[i] < '0' || suf[i] > '9') {
                                allDigits = false;
                                break;
                            }
                        }

                        if (allDigits) {
                            int v = std::atoi(suf.c_str());

                            std::string baseNoExt;
                            if (dir.empty() || dir == ".") {
                                baseNoExt = leaf + "." + suf;
                            } else {
                                baseNoExt = dir + "/" + leaf + "." + suf;
                            }

                            numbered.emplace_back(v, baseNoExt);
                        }
                    }
                }
            }
        }

        closedir(dp);
    }

    if (!numbered.empty()) {
        // pair orders by int then string
        SORT_SERIAL(numbered.begin(), numbered.end());

        for (size_t i = 0; i < numbered.size(); ++i) {
            const std::string &p = numbered[i].second;
            if (FileUtil::fileExists((p + ext).c_str())) {
                out.push_back(p);
            }
        }
        return;
    }

    std::string plainPath;
    if (dir.empty() || dir == ".") {
        plainPath = leaf;
    } else {
        plainPath = dir + "/" + leaf;
    }

    if (hasPlain && FileUtil::fileExists((plainPath + ext).c_str())) {
        out.push_back(plainPath);
    }
}

static unsigned char asciiToLower(unsigned char c) {
    if (c >= 'A' && c <= 'Z') {
        return (unsigned char)(c + ('a' - 'A'));
    }
    return c;
}

static void skipWsBuf(const char *buf, size_t n, size_t &k) {
    while (k < n && (unsigned char)buf[k] <= 32) {
        ++k;
    }
}

static bool parseAliasFileVolumes(
    const std::string &aliasPath,
    const char *dataExt,
    std::vector<std::string> &volumes) {
    MemoryMapped pal(aliasPath, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (!pal.isValid()) {
        Debug(Debug::ERROR) << "Cannot mmap alias file " << aliasPath << "\n";
        EXIT(EXIT_FAILURE);
    }

    const char *buf = reinterpret_cast<const char *>(pal.getData());
    size_t n = pal.mappedSize();

    const char needle[] = "dblist";
    size_t tlen = 6;

    size_t pos = (size_t)-1;
    if (n >= tlen) {
        for (size_t i = 0; i <= n - tlen; ++i) {
            bool match = true;

            for (size_t k = 0; k < tlen; ++k) {
                if (asciiToLower((unsigned char)buf[i + k]) != (unsigned char)needle[k]) {
                    match = false;
                    break;
                }
            }

            if (match) {
                pos = i;
                break;
            }
        }
    }

    if (pos == (size_t)-1) {
        pal.close();
        return false;
    }

    std::string baseDir = FileUtil::dirName(aliasPath);

    size_t i = pos + tlen;
    skipWsBuf(buf, n, i);

    std::vector<std::string> toks;
    while (i < n) {
        skipWsBuf(buf, n, i);
        if (i >= n) {
            break;
        }

        if (buf[i] == '"') {
            ++i;
            size_t j = i;
            while (j < n && buf[j] != '"') {
                ++j;
            }
            toks.emplace_back(std::string(buf + i, (j > i) ? (j - i) : 0));
            i = (j < n) ? j + 1 : j;
        } else {
            size_t j = i;
            while (j < n && (unsigned char)buf[j] > 32) {
                ++j;
            }
            toks.emplace_back(std::string(buf + i, (j > i) ? (j - i) : 0));
            i = j;
        }
    }
    pal.close();

    for (size_t t = 0; t < toks.size(); ++t) {
        const std::string &tok = toks[t];
        if (tok.size() >= 6 && strncasecmp(tok.c_str(), "DBLIST", 6) == 0) {
            continue;
        }

        std::string cand = tok;
        if (!(cand.size() && (cand[0] == '/' || (cand.size() > 1 && cand[1] == ':')))) {
            cand = baseDir.empty() ? cand : (baseDir + "/" + cand);
        }

        if (FileUtil::fileExists((cand + dataExt).c_str())) {
            volumes.push_back(cand);
        } else {
            pushNumberedVols(cand, dataExt, volumes);
        }
    }

    return !volumes.empty();
}

static void findVolumes(
    const std::string &db,
    int &kind,
    std::vector<std::string> &volumes) {
    // direct data files
    if (FileUtil::fileExists((db + ".psq").c_str())) {
        kind = Parameters::DBTYPE_AMINO_ACIDS;
        volumes.push_back(db);
        return;
    }
    if (FileUtil::fileExists((db + ".nsq").c_str())) {
        kind = Parameters::DBTYPE_NUCLEOTIDES;
        volumes.push_back(db);
        return;
    }

    // numbered siblings
    pushNumberedVols(db, ".psq", volumes);
    if (!volumes.empty()) {
        kind = Parameters::DBTYPE_AMINO_ACIDS;
        return;
    }
    pushNumberedVols(db, ".nsq", volumes);
    if (!volumes.empty()) {
        kind = Parameters::DBTYPE_NUCLEOTIDES;
        return;
    }

    // alias files
    std::string palPath = db + ".pal";
    std::string nalPath = db + ".nal";
    bool tried = false;

    if (FileUtil::fileExists(palPath.c_str())) {
        tried = true;
        if (parseAliasFileVolumes(palPath, ".psq", volumes)) {
            kind = Parameters::DBTYPE_AMINO_ACIDS;
            return;
        }
    }
    if (FileUtil::fileExists(nalPath.c_str())) {
        tried = true;
        if (parseAliasFileVolumes(nalPath, ".nsq", volumes)) {
            kind = Parameters::DBTYPE_NUCLEOTIDES;
            return;
        }
    }

    if (tried == false) {
        Debug(Debug::ERROR) << "No .psq/.nsq or .pal/.nal found for '" << db << "'\n";
    } else {
        Debug(Debug::ERROR) << "No volumes listed in alias for '" << db << "'\n";
    }
    EXIT(EXIT_FAILURE);
}

static bool parseIndexHeader(
    const unsigned char *idx, size_t idxSize,
    int expectedDbType,
    size_t &arraysPos,
    uint32_t &nseq
) {
    // require minimal header
    if (idxSize < 16) {
        return false;
    }

    // read format version
    size_t pos = 0;
    uint32_t ver = be32(idx + pos);
    pos += 4;
    // read seq type: nucleotide (0) protein (1)
    if (pos + 4 > idxSize) {
        return false;
    }

    uint32_t seqType = be32(idx + pos);
    pos += 4;
    int expectedKind = (expectedDbType == Parameters::DBTYPE_AMINO_ACIDS) ? 1 : 0;
    if (expectedKind != 0 && (int)seqType != expectedKind) {
        return false;
    }

    if (ver == 5U) {
        // read volume index
        if (pos + 4 > idxSize) {
            return false;
        }
        (void)be32(idx + pos);
        pos += 4;
        // read title
        if (!readBe32String(idx, idxSize, pos)) {
            return false;
        }
        // read LMDB name
        if (!readBe32String(idx, idxSize, pos)) {
            return false;
        }
        // read padded date
        if (!readBe32String(idx, idxSize, pos)) {
            return false;
        }
    } else if (ver == 4U) {
        // read title
        if (!readBe32String(idx, idxSize, pos)) {
            return false;
        }
        // read padded date
        if (!readBe32String(idx, idxSize, pos)) {
            return false;
        }
    } else {
        // unknown version
        return false;
    }
    // read #OIDs (nseq)
    if (pos + 4 > idxSize) {
        return false;
    }
    nseq = be32(idx + pos);
    pos += 4;
    // read letters (u64 LE)
    if (pos + 8 > idxSize) {
        return false;
    }
    (void)le64(idx + pos);
    pos += 8;
    // read maxlength (u32 BE)
    if (pos + 4 > idxSize) {
        return false;
    }
    (void)be32(idx + pos);
    pos += 4;

    // arrays begin here
    arraysPos = pos;

    return true;
}

static bool getVolumeNseq(const std::string& base, int kind, uint32_t& nseqOut) {
    const char *idxExt = (kind == Parameters::DBTYPE_AMINO_ACIDS) ? ".pin" : ".nin";
    MemoryMapped idx(base + idxExt, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (!idx.isValid()) {
        return false;
    }
    size_t arraysPos = 0;
    if (!parseIndexHeader(idx.getData(), idx.mappedSize(), kind, arraysPos, nseqOut)) {
        return false;
    }
    idx.close();
    return true;
}

// load (nseq+1) u32 be offsets; validates monotonic non-decreasing and <= fileSize
static bool loadOffsetArray(
    const unsigned char *idx, size_t idxSize, size_t &pos,
    size_t entries,
    size_t fileSize,
    bool requireFirstValue,
    uint32_t requiredFirst,
    std::vector<uint32_t> &out
) {
    // check bounds
    if (pos + entries * 4ULL > idxSize) {
        return false;
    }
    // read values
    out.resize(entries);
    for (size_t i = 0; i < entries; ++i) {
        out[i] = be32(idx + pos);
        pos += 4;
    }
    // validate first element if requested
    if (requireFirstValue) {
        if (out[0] != requiredFirst) {
            return false;
        }
    }
    // validate monotonic non-decreasing
    for (size_t i = 1; i < entries; ++i) {
        if (out[i] < out[i - 1]) {
            return false;
        }
    }
    // validate last element bound
    if ((size_t)out.back() > fileSize) {
        return false;
    }

    return true;
}

static bool readIdxOffsets(
    const unsigned char *idx, size_t idxSize,
    int expectedDbType,
    size_t seqFileSize, size_t hdrFileSize,
    std::vector<uint32_t> &hdr, std::vector<uint32_t> &seq,
    std::vector<uint32_t> *amb
) {
    // parse header for v4 or v5
    size_t arraysPos = 0;
    uint32_t nseq = 0;

    if (!parseIndexHeader(idx, idxSize, expectedDbType, arraysPos, nseq)) {
        return false;
    }

    size_t entries = (size_t)nseq + 1;
    size_t pos = arraysPos;
    // read header offsets H[0..nseq] with H[0] == 0
    if (!loadOffsetArray(idx, idxSize, pos, entries, hdrFileSize, true, 0, hdr)) {
        return false;
    }

    // read sequence offsets S[0..nseq] with S[0] == 1 (see writer: m_Seq.push_back(1))
    if (!loadOffsetArray(idx, idxSize, pos, entries, seqFileSize, true, 1, seq)) {
        // some legacy producers used 0; allow as a fallback but keep strong bounds
        // rewind to after H
        pos = arraysPos + entries * 4;
        seq.clear();
        if (!loadOffsetArray(idx, idxSize, pos, entries, seqFileSize, false, 0, seq)) {
            return false;
        }
    }

    if (amb != NULL && expectedDbType == Parameters::DBTYPE_NUCLEOTIDES) {
        // check if an ambiguity block is present at all
        if (pos + entries * 4ULL <= idxSize) {
            // read A[0..nseq] (writer emits nseq entries + extra S.back())
            std::vector<uint32_t> A;
            size_t apos = pos;
            if (!loadOffsetArray(idx, idxSize, apos, entries, seqFileSize, false, 0, A)) {
                amb->clear();
            } else {
                // validate per-entry bounds
                bool ok = true;
                for (size_t i = 0; i + 1 < entries; ++i) {
                    if (A[i] < seq[i] || A[i] > seq[i + 1]) {
                        ok = false;
                        break;
                    }
                }
                // validate last element equals S.back()
                if (ok) {
                    if (A.back() != seq.back()) {
                        ok = false;
                    }
                }
                if (ok) {
                    *amb = std::move(A);
                    pos = apos;
                } else {
                    amb->clear();
                }
            }
        } else {
            amb->clear();
        }
    }
    return true;
}

// static const unsigned char CLS_UNIV = 0x00;
static const unsigned char CLS_CTX  = 0x80;

struct Tlv {
    unsigned char tag = 0;
    bool constructed = false;
    unsigned char cls = 0;
    long len = 0;
    size_t vpos = 0;
};

static inline unsigned char tagNum(unsigned char tag) {
    return tag & 0x1F;
}

static bool isEoc(const unsigned char *b, size_t n, size_t p) {
    if (p + 1 >= n) {
        return false;
    }
    return (b[p] == 0x00 && b[p + 1] == 0x00);
}

static bool readLen(const unsigned char *b, size_t n, size_t &i, long &L) {
    if (i >= n) {
        return false;
    }

    unsigned char first = b[i++];
    if (first == 0x80) {
        L = -1;
        return true;
    }

    if ((first & 0x80) == 0) {
        L = first;
        return true;
    }

    int cnt = first & 0x7F;
    if (cnt == 0 || i + (size_t)cnt > n) {
        return false;
    }

    uint32_t v = 0;
    for (int k = 0; k < cnt; ++k) {
        v = (v << 8) | b[i++];
    }
    L = (long)v;

    return true;
}

static bool readTlv(const unsigned char *b, size_t n, size_t &i, Tlv &out) {
    if (i >= n) {
        return false;
    }

    out.tag = b[i++];
    out.constructed = (out.tag & 0x20) != 0;
    out.cls = out.tag & 0xC0;
    if (!readLen(b, n, i, out.len)) {
        return false;
    }
    out.vpos = i;

    return true;
}

static bool skipValue(const unsigned char *b, size_t n, size_t start, Tlv &tlv, size_t &nextPos) {
    size_t i = start;

    if (!readTlv(b, n, i, tlv)) {
        return false;
    }

    if (tlv.len >= 0) {
        nextPos = tlv.vpos + (size_t)tlv.len;

        return (nextPos <= n);
    }

    size_t p = tlv.vpos;

    while (p < n && !isEoc(b, n, p)) {
        Tlv ch;
        size_t cn;

        if (!skipValue(b, n, p, ch, cn)) {
            return false;
        }

        p = cn;
    }

    if (!isEoc(b, n, p)) {
        return false;
    }

    nextPos = p + 2;

    return true;
}

static bool getVisibleUtf8String(const unsigned char *b, size_t n, size_t nodePos, std::string &s) {
    Tlv t;
    size_t nxt;

    if (!skipValue(b, n, nodePos, t, nxt)) {
        return false;
    }

    if ((t.tag == 0x1A || t.tag == 0x0C) && t.len >= 0) {
        s.assign(reinterpret_cast<const char*>(b + t.vpos), (size_t)t.len);

        size_t a = 0;
        size_t z = s.size();

        while (a < z && (unsigned char)s[a] <= 32) {
            ++a;
        }

        while (z > a && (unsigned char)s[z - 1] <= 32) {
            --z;
        }

        s = s.substr(a, z - a);

        return true;
    }

    if (t.constructed) {
        size_t p = t.vpos;

        while (true) {
            if (t.len >= 0 && p >= t.vpos + (size_t)t.len) {
                break;
            }

            if (t.len < 0 && isEoc(b, n, p)) {
                break;
            }

            std::string inner;

            if (getVisibleUtf8String(b, n, p, inner)) {
                s = inner;
                return true;
            }

            Tlv ch;
            size_t cn;

            if (!skipValue(b, n, p, ch, cn)) {
                break;
            }

            p = cn;
        }
    }

    if ((t.cls == CLS_CTX) && !t.constructed && t.len >= 0) {
        s.assign(reinterpret_cast<const char*>(b + t.vpos), (size_t)t.len);
        return true;
    }

    return false;
}

static bool readPrim(const unsigned char *b, const struct Tlv &tt, long &val) {
    if (tt.len < 0) {
        return false;
    }

    long v = 0;
    for (long k = 0; k < tt.len; ++k) {
        v = (v << 8) | b[tt.vpos + (size_t)k];
    }
    val = v;

    return true;
}

static bool getInteger(const unsigned char *b, size_t n, size_t nodePos, long &val) {
    Tlv t;
    size_t nxt;

    if (!skipValue(b, n, nodePos, t, nxt)) {
        return false;
    }

    if (t.tag == 0x02) {
        return readPrim(b, t, val);
    }

    if (t.constructed) {
        size_t p = t.vpos;

        while (true) {
            if (t.len >= 0 && p >= t.vpos + (size_t)t.len) {
                break;
            }

            if (t.len < 0 && isEoc(b, n, p)) {
                break;
            }

            long v = 0;
            if (getInteger(b, n, p, v)) {
                val = v;
                return true;
            }

            Tlv ch;
            size_t cn;
            if (!skipValue(b, n, p, ch, cn)) {
                break;
            }

            p = cn;
        }
    }

    if ((t.cls == CLS_CTX) && !t.constructed) {
        return readPrim(b, t, val);
    }

    return false;
}

struct SeqId {
    int which = -1;
    int type  = 0;
    std::string accession;
    std::string name;
    std::string release;
    std::string version;
    std::string db;
    std::string tag;
};

static int typeFromChoice(int which) {
    switch (which) {
        case 7:  { return 1; }
        case 9:  { return 2; }
        case 4:  { return 3; }
        case 5:  { return 4; }
        case 12: { return 5; }
        case 6:  { return 6; }
        case 13: { return 7; }
        case 10: { return 8; }
        case 11: { return 9; }
        default: { return 10; }
    }
}

static bool parseTextseqId(const unsigned char *b, size_t n, size_t nodePos, SeqId &out) {
    Tlv t;
    size_t nxt;

    if (!skipValue(b, n, nodePos, t, nxt)) {
        return false;
    }

    size_t p = t.vpos;
    while (true) {
        if (t.len >= 0 && p >= t.vpos + (size_t)t.len) {
            break;
        }

        if (t.len < 0 && isEoc(b, n, p)) {
            break;
        }

        Tlv ch;
        size_t cn;
        if (!skipValue(b, n, p, ch, cn)) {
            break;
        }

        if (ch.cls == CLS_CTX) {
            unsigned char ntag = tagNum(ch.tag);

            if (ntag == 0) {
                std::string s;
                if (getVisibleUtf8String(b, n, p, s)) {
                    out.name = s;
                }
            } else if (ntag == 1) {
                std::string s;
                if (getVisibleUtf8String(b, n, p, s)) {
                    out.accession = s;
                }
            } else if (ntag == 2) {
                std::string s;
                if (getVisibleUtf8String(b, n, p, s)) {
                    out.release = s;
                }
            } else if (ntag == 3) {
                long v = 0;
                if (getInteger(b, n, p, v)) {
                    out.version = SSTR(v);
                }
            }
        }

        p = cn;
    }

    return true;
}

static bool parseDbtag(const unsigned char *b, size_t n, size_t nodePos, SeqId &out) {
    Tlv t;
    size_t nxt;

    if (!skipValue(b, n, nodePos, t, nxt)) {
        return false;
    }

    if (t.tag != 0x30) {
        return false;
    }

    bool gotDb = false;
    size_t p = t.vpos;
    while (true) {
        if (t.len >= 0 && p >= t.vpos + (size_t)t.len) {
            break;
        }

        if (t.len < 0 && isEoc(b, n, p)) {
            break;
        }

        Tlv ch;
        size_t cn;
        if (!skipValue(b, n, p, ch, cn)) {
            break;
        }

        if ((ch.tag == 0x1A || ch.tag == 0x0C) && ch.len >= 0 && !gotDb) {
            out.db.assign(reinterpret_cast<const char*>(b + ch.vpos), (size_t)ch.len);
            gotDb = true;
        } else if (ch.cls == CLS_CTX) {
            unsigned char nt = tagNum(ch.tag);
            if (nt == 0) {
                long v = 0;
                if (getInteger(b, n, p, v)) {
                    out.tag = SSTR(v);
                }
            } else if (nt == 1) {
                std::string s;
                if (getVisibleUtf8String(b, n, p, s)) {
                    out.tag = s;
                }
            }
        }

        p = cn;
    }

    return (gotDb || !out.tag.empty());
}

// parse Seq-id::pdb { mol, chain, rel }
static bool parsePdbSeqId(const unsigned char *b, size_t n, size_t nodePos, SeqId &out) {
    Tlv t;
    size_t nxt;
    if (!skipValue(b, n, nodePos, t, nxt)) {
        return false;
    }

    std::string mol;
    std::string chainId;
    std::string chainLegacy;

    // walk children of PDB-seq-id
    size_t p = t.vpos;
    while (true) {
        if ((t.len >= 0 && p >= t.vpos + (size_t)t.len) || (t.len < 0 && isEoc(b, n, p))) {
            break;
        }

        Tlv ch;
        size_t cn;
        if (!skipValue(b, n, p, ch, cn)) {
            break;
        }

        if (ch.cls == CLS_CTX) {
            unsigned char nt = tagNum(ch.tag);
            if (nt == 0 && mol.empty()) {
                // [0] mol-id (VisibleString)
                std::string s;
                if (getVisibleUtf8String(b, n, p, s)) {
                    mol = s;
                }
            } else if (nt == 1 && chainLegacy.empty()) {
                // [1] legacy chain (CHOICE { char VisibleString | num INTEGER })
                std::string s;
                if (getVisibleUtf8String(b, n, p, s) && !s.empty()) {
                    chainLegacy = s;
                } else {
                    long v = 0;
                    if (getInteger(b, n, p, v)) {
                        chainLegacy.assign(1, static_cast<char>(v));
                    }
                }
            } else if (nt == 3 && chainId.empty()) {
                // [2] ignored rel (Date-std)
                // [3] any length chain-id (VisibleString)
                std::string s;
                if (getVisibleUtf8String(b, n, p, s)) {
                    chainId = s;
                }
            }
        }
        p = cn;
    }

    if (mol.empty()) {
        return false;
    }

    for (size_t i = 0; i < mol.size(); ++i) {
        mol[i] = std::toupper((unsigned char)mol[i]);
    }

    std::string chain = !chainId.empty() ? chainId : chainLegacy;
    out.accession = chain.empty() ? mol : (mol + "_" + chain);
    out.name.clear();
    out.release.clear();
    out.version.clear();
    out.db.clear();
    out.tag.clear();

    return true;
}


static bool parseSeqId(const unsigned char *b, size_t n, unsigned char ctag, size_t vpos, long vlen, SeqId &out) {
    (void)vlen;

    out.which = (int)tagNum(ctag);
    out.type  = typeFromChoice(out.which);

    if (out.which == 14) {
        return parsePdbSeqId(b, n, vpos, out);
    }

    bool isTextseq =
        (out.which == 4 || out.which == 5 || out.which == 6 || out.which == 7 || out.which == 9 ||
         out.which == 12 || out.which == 13 || out.which == 15 || out.which == 16 ||
         out.which == 17 || out.which == 18 || out.which == 19);

    if (isTextseq) {
        size_t i = vpos;
        Tlv inner;
        size_t nx;

        if (!skipValue(b, n, i, inner, nx)) {
            return false;
        }

        if (inner.tag == 0x30 || inner.constructed) {
            parseTextseqId(b, n, vpos, out);
            return (!out.accession.empty() || !out.name.empty());
        }

        std::string acc;
        if (inner.tag == 0x1A || inner.tag == 0x0C) {
            acc.assign(reinterpret_cast<const char*>(b + inner.vpos), (size_t)inner.len);
        } else if (inner.cls == CLS_CTX && !inner.constructed && inner.len >= 0) {
            acc.assign(reinterpret_cast<const char*>(b + inner.vpos), (size_t)inner.len);
        }

        size_t a = 0;
        size_t z = acc.size();
        while (a < z && (unsigned char)acc[a] <= 32) {
            ++a;
        }

        while (z > a && (unsigned char)acc[z - 1] <= 32) {
            --z;
        }

        if (z > a) {
            out.accession = acc.substr(a, z - a);
            return true;
        }

        return false;
    }

    if (out.which == 10) {
        return parseDbtag(b, n, vpos, out);
    }

    if (out.which == 11) {
        long gi = 0;
        if (getInteger(b, n, vpos, gi)) {
            out.tag = SSTR(gi);
            return true;
        }
        return false;
    }

    return false;
}

static void parseSeqidList(const unsigned char *b, size_t n, size_t nodePos, std::vector<SeqId> &out) {
    Tlv t;
    size_t nxt;

    if (!skipValue(b, n, nodePos, t, nxt)) {
        return;
    }

    std::vector<size_t> stack;
    size_t p = t.vpos;
    while (true) {
        if (t.len >= 0 && p >= t.vpos + (size_t)t.len) {
            break;
        }

        if (t.len < 0 && isEoc(b, n, p)) {
            break;
        }

        stack.push_back(p);

        Tlv ch;
        size_t cn;
        if (!skipValue(b, n, p, ch, cn)) {
            break;
        }

        p = cn;
    }

    while (!stack.empty()) {
        size_t node = stack.back();

        stack.pop_back();

        Tlv x;
        size_t xn;
        if (!skipValue(b, n, node, x, xn)) {
            continue;
        }

        if (x.cls == CLS_CTX) {
            SeqId sid;
            if (parseSeqId(b, n, x.tag, x.vpos, x.len, sid)) {
                out.push_back(sid);
            }
        } else if (x.constructed) {
            size_t q = x.vpos;
            while (true) {
                if (x.len >= 0 && q >= x.vpos + (size_t)x.len) {
                    break;
                }

                if (x.len < 0 && isEoc(b, n, q)) {
                    break;
                }

                stack.push_back(q);

                Tlv ch;
                size_t cn;
                if (!skipValue(b, n, q, ch, cn)) {
                    break;
                }

                q = cn;
            }
        }
    }
}

// rank order: accession+version, pir|..., prf||..., accession, db|tag, gi|..., name, empty
static std::pair<int, std::string> formatId(const SeqId &sid) {
    if (sid.which == 6 /* PIR */ && !sid.accession.empty()) {
        std::string s = "pir|";
        s += sid.accession;
        s += "|";
        if (!sid.name.empty()) {
            s += sid.name;
        }
        return std::make_pair(1, s);
    }

    if (sid.which == 13 /* PRF */) {
        if (!sid.name.empty()) {
            return {4, std::string("prf||") + sid.name};
        }
        if (!sid.accession.empty()) {
            return {1, std::string("prf||") + sid.accession};
        }
    }

    if (!sid.accession.empty()) {
        if (!sid.version.empty()) {
            size_t dot = sid.accession.rfind('.');
            if (dot == std::string::npos || sid.accession.substr(dot + 1) != sid.version) {
                return std::make_pair(0, sid.accession + "." + sid.version);
            } else {
                return std::make_pair(0, sid.accession);
            }
        }
        return std::make_pair(1, sid.accession);
    }

    if (sid.type == 8 && !sid.db.empty() && !sid.tag.empty()) {
        return std::make_pair(2, sid.db + "|" + sid.tag);
    }

    if (sid.type == 9 && !sid.tag.empty()) {
        return std::make_pair(3, std::string("gi|") + sid.tag);
    }

    if (!sid.name.empty()) {
        return std::make_pair(4, sid.name);
    }

    return std::make_pair(5, std::string());
}

struct DefInfo {
    std::string title;
    long taxid = -1;
    std::vector<SeqId> seqids;
    long pig = -1;
};

static std::string parseBlastDefline(const unsigned char *blob, size_t blobSize, DefInfo& first) {
    std::string header;
    bool hasFirst = false;

    size_t i = 0;
    while (i < blobSize) {
        Tlv t;
        size_t nxt;
        if (!skipValue(blob, blobSize, i, t, nxt)) {
            break;
        }

        if (t.tag == 0x30) {
            size_t p = t.vpos;
            while (true) {
                if (t.len >= 0 && p >= t.vpos + (size_t)t.len) {
                    break;
                }

                if (t.len < 0 && isEoc(blob, blobSize, p)) {
                    break;
                }

                Tlv def;
                size_t dn;
                if (!skipValue(blob, blobSize, p, def, dn)) {
                    break;
                }

                if (def.tag == 0x30) {
                    DefInfo out;
                    size_t f = def.vpos;
                    while (true) {
                        if (def.len >= 0 && f >= def.vpos + (size_t)def.len) {
                            break;
                        }

                        if (def.len < 0 && isEoc(blob, blobSize, f)) {
                            break;
                        }

                        Tlv fld;
                        size_t fn;
                        if (!skipValue(blob, blobSize, f, fld, fn)) {
                            break;
                        }

                        if (fld.cls == CLS_CTX) {
                            unsigned char ntag = tagNum(fld.tag);
                            if (ntag == 0 && out.title.empty()) {
                                std::string s;
                                if (getVisibleUtf8String(blob, blobSize, f, s)) {
                                    out.title = s;
                                }
                            } else if (ntag == 1) {
                                parseSeqidList(blob, blobSize, f, out.seqids);
                            } else if (ntag == 2 && out.taxid < 0) {
                                long v = 0;
                                if (getInteger(blob, blobSize, f, v)) {
                                    out.taxid = v;
                                }
                            } else if (ntag == 4 && out.pig < 0) {
                                long v = 0;
                                if (getInteger(blob, blobSize, f, v)) {
                                    out.pig = v;
                                }
                            }
                        }
                        f = fn;
                    }

                    int bestRank = INT_MAX;
                    std::string id;
                    for (size_t si = 0; si < out.seqids.size(); ++si) {
                        std::pair<int, std::string> cand = formatId(out.seqids[si]);
                        // Debug(Debug::INFO) << "which =" << out.seqids[si].which
                        //                   << " type=" << out.seqids[si].type
                        //                   << " acc='" << out.seqids[si].accession << "'"
                        //                   << " name='" << out.seqids[si].name << "'"
                        //                   << " name='" << out.seqids[si].release << "'"
                        //                   << " ver='" << out.seqids[si].version << "'"
                        //                   << " db='" << out.seqids[si].db << "'"
                        //                   << " tag='" << out.seqids[si].tag << "'"
                        //                   << " => rank=" << cand.first
                        //                   << " id='" << cand.second << "'\n";
                        if (!cand.second.empty() && cand.first < bestRank) {
                            bestRank = cand.first;
                            id = cand.second;
                            if (bestRank == 0) {
                                break;
                            }
                        }
                    }
                    // EXIT(EXIT_SUCCESS);

                    std::string part;
                    if (!id.empty() && !out.title.empty()) {
                        part = id + " " + out.title;
                    } else if (!id.empty()) {
                        part = id;
                    } else {
                        part = out.title;
                    }

                    if (!hasFirst) {
                        first = out;
                        hasFirst = true;
                    }

                    if (!part.empty()) {
                        if (header.empty()) {
                            header = part;
                        } else {
                            header += " >";
                            header += part;
                        }
                    }
                }
                p = dn;
            }
            break;
        }
        i = nxt;
    }
    return header;
}

static void dumpVolumeToDb(const std::string &base,
                           int kind,
                           DBWriter &seqWriter,
                           DBWriter &hdrWriter,
                           DBWriter &lookupWriter,
                           DBWriter &mappingWriter,
                           size_t baseOid,
                           int threadIdx) {
    const char *seqExt = (kind == Parameters::DBTYPE_AMINO_ACIDS) ? ".psq" : ".nsq";
    const char *hdrExt = (kind == Parameters::DBTYPE_AMINO_ACIDS) ? ".phr" : ".nhr";
    const char *idxExt = (kind == Parameters::DBTYPE_AMINO_ACIDS) ? ".pin" : ".nin";

    std::string seqPath = base + seqExt;
    std::string hdrPath = base + hdrExt;
    std::string idxPath = base + idxExt;

    if (!FileUtil::fileExists(seqPath.c_str()) || !FileUtil::fileExists(hdrPath.c_str()) || !FileUtil::fileExists(idxPath.c_str())) {
        Debug(Debug::ERROR) << "Missing " << seqExt << "/" << hdrExt << "/" << idxExt << " for '" << base << "'\n";
        EXIT(EXIT_FAILURE);
    }

    MemoryMapped seq(seqPath, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    MemoryMapped hdr(hdrPath, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    MemoryMapped idx(idxPath, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    if (!seq.isValid() || !hdr.isValid() || !idx.isValid()) {
        Debug(Debug::ERROR) << "mmap failed for volume '" << base << "'\n";
        EXIT(EXIT_FAILURE);
    }

    std::vector<uint32_t> hdrOff;
    std::vector<uint32_t> seqOff;
    std::vector<uint32_t> ambOff;
    std::vector<uint32_t> *ambPtr = (kind == Parameters::DBTYPE_NUCLEOTIDES) ? &ambOff : NULL;
    if (!readIdxOffsets(idx.getData(), idx.mappedSize(), kind, seq.mappedSize(), hdr.mappedSize(), hdrOff, seqOff, ambPtr)) {
        Debug(Debug::ERROR) << "Cannot locate offset tables in '" << idxPath << "'\n";
        EXIT(EXIT_FAILURE);
    }

    const unsigned char *seqData = seq.getData();
    const unsigned char *hdrData = hdr.getData();

    size_t nseq = (hdrOff.size() ? hdrOff.size() - 1 : 0);
    for (size_t i = 0; i < nseq; ++i) {
        uint32_t h0 = hdrOff[i];
        uint32_t h1 = hdrOff[i + 1];
        uint32_t s0 = seqOff[i];
        uint32_t s1 = seqOff[i + 1];

        const unsigned char *blob = (h1 >= h0 && (size_t)h1 <= hdr.mappedSize()) ? (hdrData + h0) : NULL;
        size_t blobSize = (h1 >= h0 && (size_t)h1 <= hdr.mappedSize()) ? (size_t)(h1 - h0) : 0;

        unsigned int key = baseOid++;

        std::string header;
        DefInfo info;
        if (blob && blobSize) {
            header = parseBlastDefline(blob, blobSize, info);
        }
        if (header.empty()) {
            header = "OID:" + SSTR(baseOid);
        }
        if (header.back() != '\n') {
            header.append(1, '\n');
        }

        std::string seqStr;
        if (kind == Parameters::DBTYPE_AMINO_ACIDS) {
            seqStr = decodePsqSlice(seqData, s0, s1);
        } else {
            uint32_t a0 = s1;
            if (!ambOff.empty() && ambOff.size() == seqOff.size()) {
                a0 = ambOff[i];
                if (a0 < s0 || a0 > s1) {
                    a0 = s1;
                }
            }

            seqStr = decodeNucSlice(seqData, s0, a0);

            if (!ambOff.empty() && ambOff[i] < seqOff[i + 1]) {
                if (ambOff[i] < seqOff[i + 1]) {
                    applyAmbiguityPatches(seqStr, seqData, ambOff[i], seqOff[i + 1]);
                }
            }
        }
        hdrWriter.writeData(header.c_str(), header.size(), key, threadIdx);

        seqWriter.writeStart(threadIdx);
        seqWriter.writeAdd(seqStr.c_str(), seqStr.size(), threadIdx);
        seqWriter.writeAdd("\n", 1, threadIdx);
        seqWriter.writeEnd(key, threadIdx);

        if (info.taxid >= 0) {
            char buf[4096];
            int len = std::snprintf(buf, sizeof(buf), "%u\t%ld\n", key, info.taxid);
            if (len < 0 || (size_t)len >= sizeof(buf)) {
                Debug(Debug::ERROR) << "mapping line overflow\n";
                EXIT(EXIT_FAILURE);
            }
            mappingWriter.writeData(buf, (size_t)len, key, threadIdx, false, false);
        }

        std::string accession = Util::parseFastaHeader(header.c_str());
        char buf[4096];
        int len = std::snprintf(buf, sizeof(buf), "%u\t%s\t%s\n",
                                key,
                                accession.c_str(),
                                (info.pig >= 0 ? SSTR(info.pig).c_str() : "0"));
        if (len < 0 || (size_t)len >= sizeof(buf)) {
            Debug(Debug::ERROR) << "lookup line overflow\n";
            EXIT(EXIT_FAILURE);
        }
        lookupWriter.writeData(buf, (size_t)len, key, threadIdx, false, false);
    }

    seq.close();
    hdr.close();
    idx.close();
}

int convertblastdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    int kind = Parameters::DBTYPE_AMINO_ACIDS;
    std::vector<std::string> volumes;
    findVolumes(par.db1, kind, volumes);

    Debug(Debug::INFO) << "Found " << volumes.size() << " volume(s) (" << (kind == Parameters::DBTYPE_AMINO_ACIDS ? "protein" : "nucleotide") << ")\n";

    DBWriter seqWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, kind);
    seqWriter.open();

    DBWriter hdrWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    hdrWriter.open();

    std::string lookupPath = par.db2 + ".lookup";
    std::string lookupIndex = lookupPath + ".index";
    DBWriter lookupWriter(lookupPath.c_str(), lookupIndex.c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    lookupWriter.open();

    std::string mappingPath = par.db2 + "_mapping";
    std::string mappingIndex = mappingPath + ".index";
    DBWriter mappingWriter(mappingPath.c_str(), mappingIndex.c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    mappingWriter.open();

    std::vector<size_t> baseOid(volumes.size(), 0);
    size_t total = 0;
    for (size_t i = 0; i < volumes.size(); ++i) {
        uint32_t nseq = 0;
        if (!getVolumeNseq(volumes[i], kind, nseq)) {
            Debug(Debug::ERROR) << "Failed reading nseq from index for '" << volumes[i] << "'\n";
            EXIT(EXIT_FAILURE);
        }
        baseOid[i] = total;
        total += nseq;
    }

    Debug::Progress progress(volumes.size());
#pragma omp parallel
{
    unsigned int thread_idx = 0;
#ifdef OPENMP
    thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(static)
    for (size_t vi = 0; vi < volumes.size(); ++vi) {
        dumpVolumeToDb(volumes[vi], kind, seqWriter, hdrWriter, lookupWriter, mappingWriter, baseOid[vi], thread_idx);
        progress.updateProgress();
    }
}
    mappingWriter.close(true);
    FileUtil::remove(mappingIndex.c_str());

    lookupWriter.close(true);
    FileUtil::remove(lookupIndex.c_str());

    hdrWriter.close(true);
    seqWriter.close(true);

    Debug(Debug::INFO) << "Wrote " << total << " sequences\n";

    return EXIT_SUCCESS;
}
