#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"

#include "microtar.h"

#ifdef HAVE_ZLIB
#include <zlib.h>
static int file_gzwrite(mtar_t *tar, const void *data, size_t size) {
    size_t res = gzwrite((gzFile)tar->stream, data, size);
    return (res == size) ? MTAR_ESUCCESS : MTAR_EREADFAIL;
    return MTAR_ESUCCESS;
}

static int file_gzclose(mtar_t *tar) {
    gzclose((gzFile)tar->stream);
    return MTAR_ESUCCESS;
}

int mtar_gzopenw(mtar_t *tar, const char *filename) {
    // Init tar struct and functions
    memset(tar, 0, sizeof(*tar));
    tar->read = NULL;
    tar->seek = NULL;
    tar->write = file_gzwrite;
    tar->close = file_gzclose;
    tar->isFinished = 0;
    // Open file
    tar->stream = gzopen(filename, "wb");
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

int db2tar(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_LOOKUP);
    reader.open(DBReader<unsigned int>::NOSORT);

    int err;
    mtar_t tar;
    if (Util::endsWith(par.db2, ".gz") || Util::endsWith(par.db2, ".tgz")) {
#ifdef HAVE_ZLIB
        err = mtar_gzopenw(&tar, par.db2.c_str());
#else
        Debug(Debug::ERROR) << "MMseqs2 was not compiled with zlib support. Cannot create gzipped tar file.\n";
        return EXIT_FAILURE;
#endif
    } else {
        err = mtar_open(&tar, par.db2.c_str(), "w");
    }
    if (err) {
        Debug(Debug::ERROR) << "Could not open tar file " << par.db2 << " for writing\n";
        return EXIT_FAILURE;
    }

    for (size_t i = 0; i < reader.getSize(); ++i) {
        // unsigned int key = reader.getDbKey(i);
        char* data = reader.getData(i, 0);
        size_t length = std::max(reader.getEntryLen(i), (size_t)1) - 1;

        std::string name = reader.getLookupEntryName(i);

        err = mtar_write_file_header(&tar, name.c_str(), length);
        if (err) {
            Debug(Debug::ERROR) << "Could not write tar header for entry " << i << "\n";
            return EXIT_FAILURE;
        }
        err = mtar_write_data(&tar, data, length);
        if (err) {
            Debug(Debug::ERROR) << "Could not write tar data for entry " << i << "\n";
            return EXIT_FAILURE;
        }
    }
    err = mtar_write_finalize(&tar);
    if (err) {
        Debug(Debug::ERROR) << "Could not finalize tar file " << par.db2 << "\n";
        return EXIT_FAILURE;
    }
    err = mtar_close(&tar);
    if (err) {
        Debug(Debug::ERROR) << "Could not close tar file " << par.db2 << "\n";
        return EXIT_FAILURE;
    }

    reader.close();

    return EXIT_SUCCESS;
}
