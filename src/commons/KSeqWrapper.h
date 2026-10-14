#ifndef MMSEQS_KSEQWRAPPER_H
#define MMSEQS_KSEQWRAPPER_H

#include "kseq.h"
#include "KSeqBufferReader.h"

#include <string>

class KSeqWrapper {
public:
    struct KSeqEntry {
        kstring_t name;
        kstring_t sequence;
        kstring_t comment;
        kstring_t qual;
        size_t headerOffset;
        size_t sequenceOffset;
        int newlineCount;
    } entry;

    enum kseq_type {
        KSEQ_FILE,
        KSEQ_STREAM,
        KSEQ_GZIP,
        KSEQ_BZIP,
        KSEQ_BUFFER,
        KSEQ_ZSTD
    };
    kseq_type type;

    virtual bool ReadEntry() = 0;
    virtual ~KSeqWrapper() {};

protected:
    void* seq;
};

class KSeqFile : public KSeqWrapper {
public:
    KSeqFile(const char* file, size_t start = 0, size_t end = SIZE_MAX);
    bool ReadEntry();
    ~KSeqFile();
private:
    FILE* file;
    size_t end;
};


class KSeqStream : public KSeqWrapper {
public:
    KSeqStream();
    bool ReadEntry();
    ~KSeqStream();
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

class KSeqBuffer : public KSeqWrapper {
public:
    KSeqBuffer(const char* buffer, size_t length);
    bool ReadEntry();
    ~KSeqBuffer();
private:
    kseq_buffer_t d;
};

#include <zstd.h>

class KSeqZstd : public KSeqWrapper {
public:
    KSeqZstd(const char* file, size_t start = 0, size_t end = SIZE_MAX);
    bool ReadEntry();
    ~KSeqZstd();

    struct ZstdFile {
        FILE* file;
        ZSTD_DStream* dstream;
        ZSTD_inBuffer input;
        size_t toRead;
        size_t position;
        size_t end;
        size_t outputPosition;
        size_t outputEnd;
    };
private:
    ZstdFile file;
};

KSeqWrapper* KSeqFactory(const char* file, size_t start = 0, size_t end = SIZE_MAX);

#endif //MMSEQS_KSEQWRAPPER_H
