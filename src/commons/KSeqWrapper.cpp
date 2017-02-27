#include "KSeqWrapper.h"
#include "kseq.h"
#include "FileUtil.h"
#include "Util.h"

namespace KSEQFILE {
    KSEQ_INIT(int, read)
}

KSeqFile::KSeqFile(const char* fileName) {
    file = FileUtil::openFileOrDie(fileName, "r", true);
    seq = (void*) KSEQFILE::kseq_init(fileno(file));
}

bool KSeqFile::ReadEntry() {
    KSEQFILE::kseq_t* s = (KSEQFILE::kseq_t*) seq;
    int result = KSEQFILE::kseq_read(s);
    if (result < 0)
        return false;

    entry.name = std::string(s->name.s, s->name.l);
    entry.comment = std::string(s->comment.s, s->comment.l);
    entry.sequence = std::string(s->seq.s, s->seq.l);

    return true;
}

KSeqFile::~KSeqFile() {
    kseq_destroy((KSEQFILE::kseq_t*)seq);
    fclose(file);
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
}

bool KSeqGzip::ReadEntry() {
    KSEQGZIP::kseq_t* s = (KSEQGZIP::kseq_t*) seq;
    int result = KSEQGZIP::kseq_read(s);
    if (result < 0)
        return false;

    entry.name = std::string(s->name.s, s->name.l);
    entry.comment = std::string(s->comment.s, s->comment.l);
    entry.sequence = std::string(s->seq.s, s->seq.l);

    return true;
}

KSeqGzip::~KSeqGzip() {
    kseq_destroy((KSEQGZIP::kseq_t*)seq);
    gzclose(file);
}
#endif

KSeqWrapper* KSeqFactory(const char* file) {
    KSeqWrapper* kseq;
    if(Util::endsWith(".gz", file) == false) {
        kseq = new KSeqFile(file);
    }
#ifdef HAVE_ZLIB
    else {
        kseq = new KSeqGzip(file);
    }
#else
    else {
        Debug(Debug::ERROR) << "MMseqs was not compiled with zlib support. Can not read compressed input!\n";
        EXIT(EXIT_FAILURE);
    }
#endif

    return kseq;
}

