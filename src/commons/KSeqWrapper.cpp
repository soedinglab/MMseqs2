#include "KSeqWrapper.h"
#include "kseq.h"
#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"
#include <unistd.h>

namespace KSEQFILE {
    KSEQ_INIT(int, read)
}

KSeqFile::KSeqFile(const char* fileName, size_t start, size_t end) : end(end) {
    file = FileUtil::openFileOrDie(fileName, "r", true);
    seq = (void*) KSEQFILE::kseq_init(fileno(file));
    type = KSEQ_FILE;
    if (start > 0) {
        if (lseek(fileno(file), start, SEEK_SET) != (off_t) start) {
            Debug(Debug::ERROR) << "Cannot seek in file " << fileName << "\n";
            EXIT(EXIT_FAILURE);
        }
        KSEQFILE::kseq_t* s = (KSEQFILE::kseq_t*) seq;
        s->f->cur_buf_pos = start;
        int c;
        while ((c = KSEQFILE::ks_getc(s->f)) >= 0 && c != '\n');
    }
}

bool KSeqFile::ReadEntry() {
    KSEQFILE::kseq_t* s = (KSEQFILE::kseq_t*) seq;
    if (KSEQFILE::kseq_read(s) < 0 || s->headerOffset - 1 > end)
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

namespace KSEQZSTD {
    static int zstdread(KSeqZstd::ZstdFile *file, char *buffer, int length) {
        ZSTD_outBuffer output = {buffer, (size_t) length, 0};
        while (output.pos == 0) {
            if (file->input.pos == file->input.size) {
                size_t readSize = ZSTD_DStreamInSize();
                if (file->position < file->end) {
                    readSize = std::min(readSize, file->end - file->position);
                } else if (file->outputEnd == SIZE_MAX) {
                    readSize = 0;
                }
                file->input.size = fread((void*) file->input.src, sizeof(char), readSize, file->file);
                file->input.pos = 0;
                file->position += file->input.size;
                if (file->input.size == 0 && file->toRead == 0 && ferror(file->file) == 0) {
                    break;
                }
            }
            file->toRead = ZSTD_decompressStream(file->dstream, &output, &file->input);
            if (ZSTD_isError(file->toRead)) {
                Debug(Debug::ERROR) << "ZSTD_decompressStream " << ZSTD_getErrorName(file->toRead) << "\n";
                EXIT(EXIT_FAILURE);
            }
            if (file->input.size == 0 && output.pos == 0) {
                Debug(Debug::ERROR) << "Cannot read truncated zstd input\n";
                EXIT(EXIT_FAILURE);
            }
            if (file->outputEnd == SIZE_MAX && file->position >= file->end && file->input.pos == file->input.size && file->toRead == 0) {
                file->outputEnd = file->outputPosition + output.pos;
            }
        }
        file->outputPosition += output.pos;
        return output.pos;
    }

    KSEQ_INIT(KSeqZstd::ZstdFile*, zstdread)
}

KSeqZstd::KSeqZstd(const char* fileName, size_t start, size_t end) {
    file.file = FileUtil::openFileOrDie(fileName, "r", true);
    if (start > 0 && fseek(file.file, start, SEEK_SET) != 0) {
        Debug(Debug::ERROR) << "Cannot seek in file " << fileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    file.dstream = ZSTD_createDStream();
    if (file.dstream == NULL) {
        Debug(Debug::ERROR) << "ZSTD_createDStream() error\n";
        EXIT(EXIT_FAILURE);
    }
    file.input.src = malloc(ZSTD_DStreamInSize());
    file.input.size = 0;
    file.input.pos = 0;
    file.toRead = 0;
    file.position = start;
    file.end = end;
    file.outputPosition = 0;
    file.outputEnd = SIZE_MAX;
    seq = (void*) KSEQZSTD::kseq_init(&file);
    type = KSEQ_ZSTD;
    if (start > 0) {
        KSEQZSTD::kseq_t* s = (KSEQZSTD::kseq_t*) seq;
        int c;
        while ((c = KSEQZSTD::ks_getc(s->f)) >= 0 && c != '\n');
    }
}

bool KSeqZstd::ReadEntry() {
    KSEQZSTD::kseq_t* s = (KSEQZSTD::kseq_t*) seq;
    if (KSEQZSTD::kseq_read(s) < 0 || s->headerOffset - 1 > file.outputEnd)
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

KSeqZstd::~KSeqZstd() {
    kseq_destroy((KSEQZSTD::kseq_t*)seq);
    ZSTD_freeDStream(file.dstream);
    free((void*) file.input.src);
    if (fclose(file.file) != 0) {
        Debug(Debug::ERROR) << "Cannot close KSeq input file\n";
        EXIT(EXIT_FAILURE);
    }
}

KSeqWrapper* KSeqFactory(const char* file, size_t start, size_t end) {
    KSeqWrapper* kseq = NULL;
    if( strcmp(file, "stdin") == 0 ){
        kseq = new KSeqStream();
        return kseq;
    }

    if(Util::endsWith(".gz", file) == false && Util::endsWith(".bz2", file) == false && Util::endsWith(".zst", file) == false) {
        kseq = new KSeqFile(file, start, end);
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
    else if(Util::endsWith(".zst", file) == true) {
        kseq = new KSeqZstd(file, start, end);
        return kseq;
    }

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
