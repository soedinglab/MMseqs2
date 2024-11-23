#include "Debug.h"

#include <cstdio>
#include <string>
#include <cstring>

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

class GzReader {
public:
    enum Mode {
        FILE_MODE,
#ifdef HAVE_ZLIB
        GZ_MODE
#endif
    };

    GzReader(const std::string &filename) {
        if (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz") {
#ifdef HAVE_ZLIB
            mode = GZ_MODE;
            gzHandle = gzopen(filename.c_str(), "r");
            openFailed = !gzHandle;
            return;
#else
            Debug(Debug::ERROR) << "MMseqs2 was not compiled with zlib support. Cannot read compressed input\n";
            EXIT(EXIT_FAILURE);
#endif
        }
        mode = FILE_MODE;
        file = fopen(filename.c_str(), "r");
        openFailed = !file;
    }

    ~GzReader() {
        if (mode == FILE_MODE && file) fclose(file);
#ifdef HAVE_ZLIB
        else if (mode == GZ_MODE && gzHandle) gzclose(gzHandle);
#endif
    }

    bool fail() const {
        return openFailed;
    }

    bool getline(std::string &line) {
        line.clear();
        if (openFailed) return false;

        char buffer[4096];
        bool complete = false;
        while (!complete) {
            if (mode == FILE_MODE) {
                if (fgets(buffer, sizeof(buffer), file) != NULL) {
                    if (char *newline = strchr(buffer, '\n')) {
                        line.append(buffer, newline - buffer);
                        complete = true;
                    } else {
                        line.append(buffer);
                    }
                } else {
                    return !line.empty();
                }
            }
    #ifdef HAVE_ZLIB
            else if (mode == GZ_MODE) {
                if (gzgets(gzHandle, buffer, sizeof(buffer)) != NULL) {
                    if (char *newline = strchr(buffer, '\n')) {
                        line.append(buffer, newline - buffer);
                        complete = true;
                    } else {
                        line.append(buffer);
                    }
                } else {
                    return !line.empty();
                }
            }
    #endif
        }

        return true;
    }

private:
    Mode mode;
    bool openFailed = false;
    FILE *file = NULL;
#ifdef HAVE_ZLIB
    gzFile gzHandle = NULL;
#endif
};
