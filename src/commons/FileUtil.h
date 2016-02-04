#ifndef MMSEQS_FILEUTIL_H
#define MMSEQS_FILEUTIL_H

#include <cstdio>

class FileUtil {

public:
    static void errorIfFileExist(const char * file);

    static bool fileExists(const char* fileName);

    static FILE* openFileOrDie(const char * fileName, const char * mode, bool shouldExist);
};


#endif //MMSEQS_FILEUTIL_H
