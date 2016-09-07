#ifndef MMSEQS_FILEUTIL_H
#define MMSEQS_FILEUTIL_H

#include <cstdio>
#include <list>
#include <string>

class FileUtil {

public:
    static void errorIfFileExist(const char * file);

    static bool fileExists(const char* fileName);

    static bool directoryExists(const char* directoryName);

    static FILE* openFileOrDie(const char * fileName, const char * mode, bool shouldExist);

    static size_t countLines(const char* name);

    static void deleteTempFiles(std::list<std::string> tmpFiles);

    static void writeFile(std::string pathToFile, unsigned char *sh, size_t len);
};


#endif //MMSEQS_FILEUTIL_H
