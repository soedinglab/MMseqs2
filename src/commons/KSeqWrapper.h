#ifndef MMSEQS_KSEQWRAPPER_H
#define MMSEQS_KSEQWRAPPER_H

#include <string>
#include <kseq/kseq.h>

class KSeqWrapper {
public:
    struct KSeqEntry {
        kstring_t name;
        kstring_t sequence;
        kstring_t comment;
        kstring_t qual;
    } entry;

    virtual bool ReadEntry() = 0;
    virtual ~KSeqWrapper() {};

protected:
    void* seq;
};

class KSeqFile : public KSeqWrapper {
public:
    KSeqFile(const char* file);
    bool ReadEntry();
    ~KSeqFile();
private:
    FILE* file;
};

#ifdef HAVE_ZLIB
#include <zlib.h>

class KSeqGzip : public KSeqWrapper {
public:
    KSeqGzip(const char* file);
    bool ReadEntry();
    ~KSeqGzip();
private:
    gzFile file;
};
#endif

#ifdef HAVE_BZLIB
#include <bzlib.h>

class KSeqBzip : public KSeqWrapper {
public:
    KSeqBzip(const char* file);
    bool ReadEntry();
    ~KSeqBzip();
private:
    BZFILE *file;
};
#endif


KSeqWrapper* KSeqFactory(const char* file);

#endif //MMSEQS_KSEQWRAPPER_H
