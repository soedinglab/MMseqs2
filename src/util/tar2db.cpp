#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "PatternCompiler.h"

#include "microtar.h"

#ifdef OPENMP
#include <omp.h>
#endif

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
    tar->isFinished = 0;
    // Open file
    tar->stream = gzopen(filename, "rb");
    if (!tar->stream) {
        return MTAR_EOPENFAIL;
    }

#if defined(ZLIB_VERNUM) && ZLIB_VERNUM >= 0x1240
    if (gzbuffer((gzFile)tar->stream, 1 * 1024 * 1024) != 0) {
        Debug(Debug::WARNING) << "Could not set gzbuffer size, performance might be bad\n";
    }
#endif

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

    DBWriter writer(dataFile.c_str(), indexFile.c_str(), par.threads, par.compressed, par.outputDbType);
    writer.open();

    std::string lookupFile = dataFile + ".lookup";
    DBWriter lookupWriter(lookupFile.c_str(), (lookupFile + ".index").c_str(), par.threads, 0, Parameters::DBTYPE_OMIT_FILE);
    lookupWriter.open();

    Debug::Progress progress;

    size_t globalKey = 0;
    for (size_t i = 0; i < filenames.size(); i++) {
        char buffer[4096];
        size_t len = snprintf(buffer, sizeof(buffer), "%zu\t%s\n", i, FileUtil::baseName(filenames[i]).c_str());
        int written = fwrite(buffer, sizeof(char), len, source);
        if (written != (int) len) {
            Debug(Debug::ERROR) << "Cannot write to source file " << sourceFile << "\n";
            EXIT(EXIT_FAILURE);
        }

#ifdef OPENMP
        int localThreads = par.threads;
#endif
        mtar_t tar;
        if (Util::endsWith(".tar.gz", filenames[i]) || Util::endsWith(".tgz", filenames[i])) {
#ifdef HAVE_ZLIB
            if (mtar_gzopen(&tar, filenames[i].c_str()) != MTAR_ESUCCESS) {
                Debug(Debug::ERROR) << "Cannot open file " << filenames[i] << "\n";
                EXIT(EXIT_FAILURE);
            }
#ifdef OPENMP
            localThreads = 1;
#endif
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

#pragma omp parallel shared(tar) num_threads(localThreads)
        {
            char buffer[4096];
            size_t bufferSize = 1024 * 1024;
            char *dataBuffer = (char *) malloc(bufferSize);
            size_t inflateSize = 1024 * 1024;
            char *inflateBuffer = (char *) malloc(inflateSize);
            mtar_header_t header;
            size_t currentKey = 0;
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
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            std::string name;
            bool proceed = true;
            while (proceed) {
                bool writeEntry = true;
#pragma omp critical
                {
                    if (tar.isFinished == 0 && (mtar_read_header(&tar, &header)) != MTAR_ENULLRECORD) {
                        // GNU tar has special blocks for long filenames
                        if (header.type == MTAR_TGNU_LONGNAME || header.type == MTAR_TGNU_LONGLINK) {
                            if (header.size > bufferSize) {
                                bufferSize = header.size * 1.5;
                                dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                            }
                            if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                Debug(Debug::ERROR) << "Cannot read entry " << header.name << "\n";
                                EXIT(EXIT_FAILURE);
                            }
                            name.assign(dataBuffer, header.size);
                            // skip to next record
                            if (mtar_read_header(&tar, &header) == MTAR_ENULLRECORD) {
                                Debug(Debug::ERROR) << "Tar truncated after entry " << name << "\n";
                                EXIT(EXIT_FAILURE);
                            }
                        } else {
                            name = header.name;
                        }
                        if (header.type == MTAR_TREG || header.type == MTAR_TCONT || header.type == MTAR_TOLDREG) {
                            if (include.isMatch(name.c_str()) == false || exclude.isMatch(name.c_str()) == true) {
                                __sync_fetch_and_add(&(globalKey), 1);
                                proceed = true;
                                writeEntry = false;
                            } else {
                                if (header.size > bufferSize) {
                                    bufferSize = header.size * 1.5;
                                    dataBuffer = (char *) realloc(dataBuffer, bufferSize);
                                }
                                if (mtar_read_data(&tar, dataBuffer, header.size) != MTAR_ESUCCESS) {
                                    Debug(Debug::ERROR) << "Cannot read entry " << name << "\n";
                                    EXIT(EXIT_FAILURE);
                                }
                                proceed = true;
                                writeEntry = true;
                                currentKey = __sync_fetch_and_add(&(globalKey), 1);
                            }
                        } else {
                            proceed = true;
                            writeEntry = false;
                        }
                    } else {
                        tar.isFinished = 1;
                        proceed = false;
                        writeEntry = false;
                    }
                }
                if (proceed && writeEntry) {
                    progress.updateProgress();
                    if (Util::endsWith(".gz", name)) {
#ifdef HAVE_ZLIB
                        inflateReset(&strm);
                        writer.writeStart(thread_idx);
                        strm.avail_in = header.size;
                        strm.next_in = (unsigned char *) dataBuffer;
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
                                    Debug(Debug::ERROR) << "Gzip error " << err << " entry " << name << "\n";
                                    EXIT(EXIT_FAILURE);
                            }
                            have = CHUNK - strm.avail_out;
                            writer.writeAdd((const char *) out, have, thread_idx);
                        } while (strm.avail_out == 0);
                        writer.writeEnd(currentKey, thread_idx);
#else
                        Debug(Debug::ERROR) << "MMseqs2 was not compiled with zlib support. Cannot read compressed input.\n";
                EXIT(EXIT_FAILURE);
#endif
                    } else if (Util::endsWith(".bz2", name)) {
#ifdef HAVE_BZLIB
                        unsigned int entrySize = inflateSize;
                        int err;
                        while ((err = BZ2_bzBuffToBuffDecompress(inflateBuffer, &entrySize, dataBuffer, header.size, 0,
                                                                 0) == BZ_OUTBUFF_FULL)) {
                            entrySize = inflateSize = inflateSize * 1.5;
                            inflateBuffer = (char *) realloc(inflateBuffer, inflateSize);
                        }
                        if (err != BZ_OK) {
                            Debug(Debug::ERROR) << "Could not decompress " << name << "\n";
                            EXIT(EXIT_FAILURE);
                        }
                        writer.writeData(inflateBuffer, entrySize, currentKey, thread_idx);
#else
                        Debug(Debug::ERROR) << "MMseqs2 was not compiled with bzlib support. Cannot read compressed input.\n";
                EXIT(EXIT_FAILURE);
#endif
                    } else {
                        writer.writeData(dataBuffer, header.size, currentKey, thread_idx);
                    }
                    size_t len = snprintf(buffer, sizeof(buffer), "%zu\t%s\t%zu\n", currentKey, FileUtil::baseName(name).c_str(), i);
                    lookupWriter.writeData(buffer, len, thread_idx, false, false);
                }
            }

#ifdef HAVE_ZLIB
            inflateEnd(&strm);
#endif
            free(inflateBuffer);
            free(dataBuffer);
        } // end omp

        mtar_close(&tar);
    } // filename for
    writer.close();
    lookupWriter.close(true);
    FileUtil::remove(lookupWriter.getIndexFileName());
    if (fclose(source) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << sourceFile << "\n";
        EXIT(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
