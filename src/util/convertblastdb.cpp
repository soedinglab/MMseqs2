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

static inline char iupacFromMask(unsigned char m) {
    static const char LUT[16] = {'-', 'A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
    if (m < 16) {
        return LUT[m];
    }
    return 'N';
}

static const char STDAA[28] = {
    '-', 'A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','U','*','O','J'
};

static std::string decodeNucSlice(const unsigned char *nsq, uint32_t s, uint32_t e) {
    if (e > s && nsq[e - 1] == 0) {
        --e;
    }

    std::string out;
    out.resize(e - s);
    for (uint32_t i = s; i < e; ++i) {
        unsigned char b = nsq[i];
        out[i - s] = iupacFromMask((unsigned char)(b & 0x0F));
    }

    return out;
}

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
    if (dp != nullptr) {
        struct dirent* de = nullptr;
        while ((de = readdir(dp)) != nullptr) {
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

static bool tryPinAt(const unsigned char *pin, size_t pinSize, size_t pos,
                     size_t psqSize, size_t phrSize,
                     std::vector<uint32_t> &hdr, std::vector<uint32_t> &seq) {
    if (pos + 4 > pinSize) {
        return false;
    }

    uint32_t nseq = be32(pin + pos);

    pos += 4;

    if (nseq == 0) {
        return false;
    }

    if (pos + 8 + 4 > pinSize) {
        return false;
    }

    pos += 8 + 4;

    size_t remain = pinSize - pos;

    if (remain < 8) {
        return false;
    }

    size_t maxEntries = remain / 8;

    if ((size_t)nseq + 1 > maxEntries) {
        return false;
    }

    size_t entries = (size_t)nseq + 1;

    std::vector<uint32_t> H(entries);
    std::vector<uint32_t> S(entries);

    for (size_t i = 0; i < entries; ++i) {
        if (pos + 4 > pinSize) {
            return false;
        }
        H[i] = be32(pin + pos);
        pos += 4;
    }

    for (size_t i = 0; i < entries; ++i) {
        if (pos + 4 > pinSize) {
            return false;
        }
        S[i] = be32(pin + pos);
        pos += 4;
    }

    for (size_t i = 0; i + 1 < entries; ++i) {
        if (H[i] > H[i + 1]) {
            return false;
        }
    }

    for (size_t i = 0; i + 1 < entries; ++i) {
        if (S[i] > S[i + 1]) {
            return false;
        }
    }

    if ((size_t)H.back() > phrSize || (size_t)S.back() > psqSize) {
        return false;
    }

    hdr.swap(H);
    seq.swap(S);

    return true;
}


static bool readIdxOffsets(
    const unsigned char *idx, size_t idxSize,
    int expectedDbType,
    size_t seqFileSize, size_t hdrFileSize,
    std::vector<uint32_t> &hdr, std::vector<uint32_t> &seq) {
    if (idxSize < 12) {
        Debug(Debug::ERROR) << "BlastDB idx too small: idxSize " << idxSize << "\n";
        return false;
    }

    size_t o = 0;

    uint32_t ver = be32(idx + o);
    (void)ver;
    o += 4;

    uint32_t dbt = be32(idx + o);
    o += 4;

    int expectedKind = (expectedDbType == Parameters::DBTYPE_AMINO_ACIDS) ? 1 : 0;
    if (expectedKind != 0 && (int)dbt != expectedKind) {
        Debug(Debug::ERROR) << "BlastDB dbtype mismatch (saw " << dbt << ")\n";
        return false;
    }

    uint32_t tlen = be32(idx + o);
    o += 4;

    if (o + tlen + 4 > idxSize) {
        Debug(Debug::ERROR) << "Bad BlastDB header: titleLen " << tlen << " exceeds idx\n";
        return false;
    }

    o += tlen;

    uint32_t slen = be32(idx + o);
    o += 4;

    if (o + slen > idxSize) {
        Debug(Debug::ERROR) << "Bad BlastDB timestamp: tsLen " << slen << " exceeds idx\n";
        return false;
    }

    size_t tsEnd = o + slen;

    for (int pad = 0; pad < 8; ++pad) {
        if (tryPinAt(idx, idxSize, tsEnd + (size_t)pad, seqFileSize, hdrFileSize, hdr, seq)) {
            return true;
        }
    }

    size_t maxProbe = std::min(idxSize, tsEnd + (size_t)4096);
    for (size_t pos = tsEnd; pos + 12 < maxProbe; ++pos) {
        if (tryPinAt(idx, idxSize, pos, seqFileSize, hdrFileSize, hdr, seq)) {
            return true;
        }
    }

    return false;
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

static bool getInteger(const unsigned char *b, size_t n, size_t nodePos, long &val) {
    Tlv t;
    size_t nxt;

    if (!skipValue(b, n, nodePos, t, nxt)) {
        return false;
    }

    auto readPrim = [&](const Tlv &tt) -> bool {
        if (tt.len < 0) {
            return false;
        }

        long v = 0;
        for (long k = 0; k < tt.len; ++k) {
            v = (v << 8) | b[tt.vpos + (size_t)k];
        }
        val = v;

        return true;
    };

    if (t.tag == 0x02) {
        return readPrim(t);
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
        return readPrim(t);
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

        if (ch.tag == 0x1A && ch.len >= 0 && !gotDb) {
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

static bool parseSeqId(const unsigned char *b, size_t n, unsigned char ctag, size_t vpos, long vlen, SeqId &out) {
    (void)vlen;

    out.which = (int)tagNum(ctag);
    out.type  = typeFromChoice(out.which);

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

struct DefInfo {
    std::string title;
    long taxid = -1;
    std::vector<SeqId> seqids;
    long pig = -1;
};

static void parseBlastDeflineSet(const unsigned char *blob, size_t blobSize, DefInfo &out) {
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
                }
                p = dn;
            }
            break;
        }
        i = nxt;
    }
}

static std::string choosePrimaryAccession(const std::vector<SeqId> &sids) {
    std::vector<SeqId> by[11];

    for (const auto &s : sids) {
        int t = (s.type < 0 || s.type > 10) ? 0 : s.type;
        by[t].push_back(s);
    }

    auto pickAccVer = [](const std::vector<SeqId> &v) -> std::string {
        for (const auto &s : v) {
            if (!s.accession.empty()) {
                if (!s.version.empty()) {
                    size_t dot = s.accession.rfind('.');

                    if (dot == std::string::npos || s.accession.substr(dot + 1) != s.version) {
                        return s.accession + "." + s.version;
                    } else {
                        return s.accession;
                    }
                }

                return s.accession;
            }
        }

        for (const auto &s : v) {
            if (!s.name.empty()) {
                return s.name;
            }
        }

        return std::string();
    };

    const int order[] = {1, 2, 3, 4, 5, 6, 7};
    for (int o : order) {
        std::string s = pickAccVer(by[o]);
        if (!s.empty()) {
            return s;
        }
    }

    for (const auto &s : by[8]) {
        if (!s.db.empty() && !s.tag.empty()) {
            return s.db + "|" + s.tag;
        }
    }

    for (const auto &s : by[9]) {
        if (!s.tag.empty()) {
            return std::string("gi|") + s.tag;
        }
    }

    return std::string();
}

static std::string makeDbKey(const std::vector<SeqId> &sids, const std::string &primary) {
    for (const auto &s : sids) {
        if (s.type == 8 && !s.db.empty() && !s.tag.empty()) {
            return s.db + "|" + s.tag;
        }
    }

    for (const auto &s : sids) {
        if (s.type == 9 && !s.tag.empty()) {
            return std::string("gi|") + s.tag;
        }
    }

    return primary;
}

static void dumpVolumeToDb(const std::string &base,
                           int kind,
                           DBWriter &seqWriter,
                           DBWriter &hdrWriter,
                           DBWriter &lookupWriter,
                           FILE *mapFp,
                           size_t &globalOid,
                           int threadIdx = 0) {
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
    if (!readIdxOffsets(idx.getData(), idx.mappedSize(), kind, seq.mappedSize(), hdr.mappedSize(), hdrOff, seqOff)) {
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

        const unsigned char *blob = (h1 >= h0 && (size_t)h1 <= hdr.mappedSize()) ? (hdrData + h0) : nullptr;
        size_t blobSize = (h1 >= h0 && (size_t)h1 <= hdr.mappedSize()) ? (size_t)(h1 - h0) : 0;

        DefInfo info;
        if (blob && blobSize) {
            parseBlastDeflineSet(blob, blobSize, info);
        }

        std::string primary = choosePrimaryAccession(info.seqids);
        std::string header;
        if (!primary.empty() && !info.title.empty()) {
            header = primary + " " + info.title;
        } else if (!primary.empty()) {
            header = primary;
        } else if (!info.title.empty()) {
            header = info.title;
        } else {
            header = "OID:" + SSTR(globalOid);
        }
        if (header.back() != '\n') {
            header.append(1, '\n');
        }

        std::string seqStr;
        if (kind == Parameters::DBTYPE_AMINO_ACIDS) {
            seqStr = decodePsqSlice(seqData, s0, s1);
        } else {
            seqStr = decodeNucSlice(seqData, s0, s1);
        }
        unsigned int key = (unsigned int)globalOid;
        hdrWriter.writeData(header.c_str(), header.size(), key, threadIdx);

        seqWriter.writeStart(threadIdx);
        seqWriter.writeAdd(seqStr.c_str(), seqStr.size(), threadIdx);
        seqWriter.writeAdd("\n", 1, threadIdx);
        seqWriter.writeEnd(key, threadIdx);

        if (mapFp && info.taxid >= 0) {
            if (std::fprintf(mapFp, "%u\t%ld\n", key, info.taxid) < 0) {
                Debug(Debug::ERROR) << "Failed writing mapping line\n";
                EXIT(EXIT_FAILURE);
            }
        }

        std::string accession = makeDbKey(info.seqids, primary);
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

        __sync_fetch_and_add(&globalOid, 1);
    }

    seq.close();
    hdr.close();
    idx.close();
}

int convertblastdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.threads = 1; // TODO
    par.parseParameters(argc, argv, command, true, 0, 0);
    par.threads = 1;

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

    std::string mapPath = par.db2 + "_mapping";
    FILE *mapFp = std::fopen(mapPath.c_str(), "w");
    if (!mapFp) {
        Debug(Debug::ERROR) << "Cannot open " << mapPath << " for writing\n";
        return EXIT_FAILURE;
    }

    size_t globalOid = 0;
    Debug::Progress progress(volumes.size());
// #pragma omp parallel for schedule(static)
    for (size_t vi = 0; vi < volumes.size(); ++vi) {
        progress.updateProgress();
        dumpVolumeToDb(volumes[vi], kind, seqWriter, hdrWriter, lookupWriter, mapFp, globalOid, 0);
    }

    if (std::fclose(mapFp) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << mapPath << "\n";
        return EXIT_FAILURE;
    }

    lookupWriter.close(true);
    FileUtil::remove(lookupIndex.c_str());

    hdrWriter.close(true);
    seqWriter.close(true);

    Debug(Debug::INFO) << "Wrote " << globalOid << " sequences\n";

    return EXIT_SUCCESS;
}
