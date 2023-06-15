#include "KSeqWrapper.h"
#include "kseq.h"
#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"
#include <unistd.h>

namespace KSEQFILE {
    KSEQ_INIT(int, read)
}

KSeqFile::KSeqFile(const char* fileName) {
    file = FileUtil::openFileOrDie(fileName, "r", true);
    seq = (void*) KSEQFILE::kseq_init(fileno(file));
    type = KSEQ_FILE;
}

bool KSeqFile::ReadEntry() {
    KSEQFILE::kseq_t* s = (KSEQFILE::kseq_t*) seq;
    if (KSEQFILE::kseq_read(s) < 0)
        return false;
    entry.headerOffset = s->headerOffset;
    entry.sequenceOffset = s->sequenceOffset;
    entry.newlineCount = s->newlineCount;
    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;

    return true;
}

KSeqFile::~KSeqFile() {
    kseq_destroy((KSEQFILE::kseq_t*)seq);
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close KSeq input file\n";
        EXIT(EXIT_FAILURE);
    }
}


namespace KSEQSTREAM {
    KSEQ_INIT(int, read)
}

KSeqStream::KSeqStream() {
    seq = (void*) KSEQSTREAM::kseq_init(STDIN_FILENO);
    type = KSEQ_STREAM;
}

bool KSeqStream::ReadEntry() {
    KSEQSTREAM::kseq_t* s = (KSEQSTREAM::kseq_t*) seq;
    if (KSEQSTREAM::kseq_read(s) < 0)
        return false;

    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;

    return true;
}

KSeqStream::~KSeqStream() {
    kseq_destroy((KSEQSTREAM::kseq_t*)seq);
}

#ifdef HAVE_ZLIB
namespace KSEQGZIP {
    KSEQ_INIT(gzFile, gzread)
}

KSeqGzip::KSeqGzip(const char* fileName) {
    if(FileUtil::fileExists(fileName) == false) {
        errno = ENOENT;
        perror(fileName);
        EXIT(EXIT_FAILURE);
    }

    file = gzopen(fileName, "r");
    if(file == NULL) {
        perror(fileName); EXIT(EXIT_FAILURE);
    }

    seq = (void*) KSEQGZIP::kseq_init(file);
    type = KSEQ_GZIP;
}

bool KSeqGzip::ReadEntry() {
    KSEQGZIP::kseq_t* s = (KSEQGZIP::kseq_t*) seq;
    if (KSEQGZIP::kseq_read(s) < 0)
        return false;

    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;
    entry.headerOffset = 0;
    entry.sequenceOffset = 0;
    entry.newlineCount = s->newlineCount;

    return true;
}

KSeqGzip::~KSeqGzip() {
    kseq_destroy((KSEQGZIP::kseq_t*)seq);
    gzclose(file);
}
#endif


#ifdef HAVE_BZLIB
namespace KSEQBZIP {
    KSEQ_INIT(BZFILE *, BZ2_bzread)
}

KSeqBzip::KSeqBzip(const char* fileName) {
    if(FileUtil::fileExists(fileName) == false) {
        errno = ENOENT;
        perror(fileName);
        EXIT(EXIT_FAILURE);
    }
    FILE *fp = FileUtil::openFileOrDie(fileName, "r+b", true);
    int bzError;
    file = BZ2_bzReadOpen(&bzError, fp, 0, 0, NULL, 0);
    if(bzError != 0){
        perror(fileName); EXIT(EXIT_FAILURE);
    }
    seq = (void*) KSEQBZIP::kseq_init(file);
    type = KSEQ_BZIP;
}

bool KSeqBzip::ReadEntry() {
    KSEQBZIP::kseq_t* s = (KSEQBZIP::kseq_t*) seq;
    if (KSEQBZIP::kseq_read(s) < 0)
        return false;

    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;
    entry.headerOffset = 0;
    entry.sequenceOffset = 0;
    entry.newlineCount = s->newlineCount;

    return true;
}

KSeqBzip::~KSeqBzip() {
    kseq_destroy((KSEQBZIP::kseq_t*)seq);
    int bzError;
    BZ2_bzReadClose(&bzError, file);
}
#endif

KSeqWrapper* KSeqFactory(const char* file) {
    KSeqWrapper* kseq = NULL;
    if( strcmp(file, "stdin") == 0 ){
        kseq = new KSeqStream();
        return kseq;
    }

    if(Util::endsWith(".gz", file) == false && Util::endsWith(".bz2", file) == false ) {
        kseq = new KSeqFile(file);
        return kseq;
    }
#ifdef HAVE_ZLIB
    else if(Util::endsWith(".gz", file) == true) {
        kseq = new KSeqGzip(file);
        return kseq;
    }
#else
    else if(Util::endsWith(".gz", file) == true) {
        Debug(Debug::ERROR) << "MMseqs was not compiled with zlib support. Can not read compressed input!\n";
        EXIT(EXIT_FAILURE);
    }
#endif

#ifdef HAVE_BZLIB
    else if(Util::endsWith(".bz2", file) == true) {
        kseq = new KSeqBzip(file);
        return kseq;
    }
#else
    else if(Util::endsWith(".bz2", file) == true) {
        Debug(Debug::ERROR) << "MMseqs was not compiled with bz2lib support. Can not read compressed input!\n";
        EXIT(EXIT_FAILURE);
    }
#endif

    return kseq;
}

namespace KSEQBUFFER {
    KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)
}

KSeqBuffer::KSeqBuffer(const char* buffer, size_t length) {
    d.buffer = (char*)buffer;
    d.length = length;
    d.position = 0;
    seq = (void*) KSEQBUFFER::kseq_init(&d);
    type = KSEQ_BUFFER;
}

bool KSeqBuffer::ReadEntry() {
    KSEQBUFFER::kseq_t* s = (KSEQBUFFER::kseq_t*) seq;
    if (KSEQBUFFER::kseq_read(s) < 0)
        return false;
    entry.headerOffset = s->headerOffset;
    entry.sequenceOffset = s->sequenceOffset;
    entry.newlineCount = s->newlineCount;
    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;

    return true;
}

KSeqBuffer::~KSeqBuffer() {
    kseq_destroy((KSEQBUFFER::kseq_t*)seq);
}
