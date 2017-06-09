#ifndef MMSEQS_KSEQWRAPPER_H
#define MMSEQS_KSEQWRAPPER_H

#include <string>

class KSeqWrapper {
public:
    struct KSeqEntry {
        std::string name;
        std::string comment;
        std::string sequence;
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

KSeqWrapper* KSeqFactory(const char* file);

#endif //MMSEQS_KSEQWRAPPER_H
