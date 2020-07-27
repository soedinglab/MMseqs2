#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "PatternCompiler.h"

#include "microtar.h"

#ifdef HAVE_ZLIB
#include <zlib.h>
static int file_gzread(mtar_t *tar, void *data, size_t size) {
    size_t res = gzread((gzFile)tar->stream, data, size);
    return (res == size) ? MTAR_ESUCCESS : MTAR_EREADFAIL;
}

static int file_gzseek(mtar_t *tar, long offset, int whence) {
    int res = gzseek((gzFile)tar->stream, offset, whence);
    return (res != -1) ? MTAR_ESUCCESS : MTAR_ESEEKFAIL;
}

static int file_gzclose(mtar_t *tar) {
    gzclose((gzFile)tar->stream);
    return MTAR_ESUCCESS;
}

int mtar_gzopen(mtar_t *tar, const char *filename) {
    // Init tar struct and functions
    memset(tar, 0, sizeof(*tar));
    tar->read = file_gzread;
    tar->seek = file_gzseek;
    tar->close = file_gzclose;

    // Open file
    tar->stream = gzopen(filename, "rb");
    if (!tar->stream) {
        return MTAR_EOPENFAIL;
    }

    // Return ok
    return MTAR_ESUCCESS;
}
#endif

#ifdef HAVE_BZLIB
#include <bzlib.h>
#endif

int tar2db(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::vector<std::string> filenames(par.filenames);
    for (size_t i = 0; i < filenames.size(); i++) {
        if (FileUtil::directoryExists(filenames[i].c_str()) == true) {
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory.\n";
            EXIT(EXIT_FAILURE);
        }
    }

    PatternCompiler include(par.tarInclude.c_str());
    PatternCompiler exclude(par.tarExclude.c_str());

    std::string dataFile = filenames.back();
    filenames.pop_back();
    std::string indexFile = dataFile + ".index";

    std::string sourceFile = dataFile + ".source";
    FILE *source = FileUtil::openAndDelete(sourceFile.c_str(), "w");

    std::string lookupFile = dataFile + ".lookup";
    FILE *lookup = FileUtil::openAndDelete(lookupFile.c_str(), "w");

    DBWriter writer(dataFile.c_str(), indexFile.c_str(), 1, par.compressed, par.outputDbType);
    writer.open();
    Debug::Progress progress;
    char buffer[4096];

#ifdef HAVE_ZLIB
    const unsigned int CHUNK = 128 * 1024;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];
    z_stream strm;
    memset(&strm, 0, sizeof(z_stream));
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.next_in = in;
    strm.avail_in = 0;
    int status = inflateInit2(&strm, 15 | 32);
    if (status < 0) {
        Debug(Debug::ERROR) << "Cannot initialize zlib stream\n";
        EXIT(EXIT_FAILURE);
    }
#endif

    size_t key = 0;
    for (size_t i = 0; i < filenames.size(); i++) {
        size_t len = snprintf(buffer, sizeof(buffer), "%zu\t%s\n", i, FileUtil::baseName(filenames[i]).c_str());
        int written = fwrite(buffer, sizeof(char), len, source);
        if (written != (int) len) {
            Debug(Debug::ERROR) << "Cannot write to source file " << sourceFile << "\n";
            EXIT(EXIT_FAILURE);
        }

        mtar_t tar;
        if (Util::endsWith(".tar.gz", filenames[i]) || Util::endsWith(".tgz", filenames[i])) {
#ifdef HAVE_ZLIB
            if (mtar_gzopen(&tar, filenames[i].c_str()) != MTAR_ESUCCESS) {
                Debug(Debug::ERROR) << "Cannot open file " << filenames[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
#else
            Debug(Debug::ERROR) << "MMseqs2 was not compiled with zlib support. Cannot read compressed input.\n";
            EXIT(EXIT_FAILURE);
#endif
        } else {
            if (mtar_open(&tar, filenames[i].c_str()) != MTAR_ESUCCESS) {
                Debug(Debug::ERROR) << "Cannot open file " << filenames[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
        }

        size_t bufferSize = 10 * 1024;
        char* dataBuffer = (char*) malloc(bufferSize);

        size_t inflateSize = 10 * 1024;
        char* inflateBuffer = (char*) malloc(inflateSize);

        mtar_header_t header;
        while ((mtar_read_header(&tar, &header)) != MTAR_ENULLRECORD ) {
            if (header.type != MTAR_TREG) {
                continue;
            }
            progress.updateProgress();
            if (include.isMatch(header.name) == false || exclude.isMatch(header.name) == true) {
                key++;
                continue;
            }
            if (header.size > bufferSize) {
                bufferSize = header.size * 1.5;
                dataBuffer = (char*)realloc(dataBuffer, bufferSize);
            }
            if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                Debug(Debug::ERROR) << "Cannot read entry " << header.name << "\n";
                EXIT(EXIT_FAILURE);
            }

            if (Util::endsWith(".gz", header.name)) {
#ifdef HAVE_ZLIB
                inflateReset(&strm);
                writer.writeStart(0);
                strm.avail_in = header.size;
                strm.next_in = (unsigned char*)dataBuffer;
                do {
                    unsigned have;
                    strm.avail_out = CHUNK;
                    strm.next_out = out;
                    int err = inflate(&strm, Z_NO_FLUSH);
                    switch (err) {
                        case Z_OK:
                        case Z_STREAM_END:
                        case Z_BUF_ERROR:
                            break;
                        default:
                            inflateEnd(&strm);
                            Debug(Debug::ERROR) << "Gzip error " << err << " entry " << header.name << "\n";
                            EXIT(EXIT_FAILURE);
                    }
                    have = CHUNK - strm.avail_out;
                    writer.writeAdd((const char*)out, have, 0);
                } while (strm.avail_out == 0);
                writer.writeEnd(key, 0);
#else
                Debug(Debug::ERROR) << "MMseqs2 was not compiled with zlib support. Cannot read compressed input.\n";
                EXIT(EXIT_FAILURE);
#endif
            } else if (Util::endsWith(".bz2", header.name)) {
#ifdef HAVE_BZLIB
                unsigned int entrySize = inflateSize;
                int err;
                while ((err = BZ2_bzBuffToBuffDecompress(inflateBuffer, &entrySize, dataBuffer, header.size, 0, 0) == BZ_OUTBUFF_FULL)) {
                    entrySize = inflateSize = inflateSize * 1.5;
                    inflateBuffer = (char*)realloc(inflateBuffer, inflateSize);
                }
                if (err != BZ_OK) {
                    Debug(Debug::ERROR) << "Could not decompress " << header.name  << "\n";
                    EXIT(EXIT_FAILURE);
                }
                writer.writeData(inflateBuffer, entrySize, key, 0);
#else
                Debug(Debug::ERROR) << "MMseqs2 was not compiled with bzlib support. Cannot read compressed input.\n";
                EXIT(EXIT_FAILURE);
#endif
            } else {
                writer.writeData(dataBuffer, header.size, key, 0);
            }
            size_t len = snprintf(buffer, sizeof(buffer), "%zu\t%s\t%zu\n", key, FileUtil::baseName(header.name).c_str(), i);
            int written = fwrite(buffer, sizeof(char), len, lookup);
            if (written != (int) len) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            key++;
        }

        free(inflateBuffer);
        free(dataBuffer);

        mtar_close(&tar);
    }
    fclose(lookup);
    fclose(source);
    writer.close();

#ifdef HAVE_ZLIB
    inflateEnd(&strm);
#endif

    return EXIT_SUCCESS;
}
