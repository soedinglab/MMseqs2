#ifndef MMSEQS_FILEUTIL_H
#define MMSEQS_FILEUTIL_H

#include <cstdio>
#include <list>
#include <string>
#include <vector>
#include <cstddef>
#include <utility>
#include "Parameters.h"

class FileUtil {

public:
    static bool fileExists(const char *fileName);

    static bool directoryExists(const char *directoryName);

    static FILE* openFileOrDie(const char *fileName, const char *mode, bool shouldExist);

    static size_t countLines(const char *name);

    static bool makeDir(const char *dirName, const int mode = 0777);

    static void deleteTempFiles(const std::list<std::string> &tmpFiles);

    static std::string getRealPathFromSymLink(const std::string path);

    static std::string getHashFromSymLink(const std::string path);

    static void* mmapFile(FILE * file, size_t *dataSize);

    static void munmapData(void * ptr, size_t dataSize);

    static void writeFile(const std::string &pathToFile, const unsigned char *sh, size_t len);

    static std::string dirName(const std::string &file);

    static std::string baseName(const std::string &file);

    static size_t getFreeSpace(const char *dir);

    static std::string getCurrentWorkingDirectory();

    static void symlinkAlias(const std::string &file, const std::string &alias);
    static void symlinkAbs(const std::string &target, const std::string &link);

    static size_t getFileSize(const std::string &fileName);

    static bool symlinkExists(const std::string &path);

    static void copyFile(const char *src, const char *dst);

    static FILE *openAndDelete(const char *fileName, const char *mode);

    static std::vector<std::string> findDatafiles(const char * datafiles);

    static void remove(const char * file);

    static void move(const char * src, const char * dst);

    static int parseDbType(const char *name);

    static std::string createTemporaryDirectory(const std::string& basePath, const std::string& subDirectory);

    static void fixRlimitNoFile();
};


#endif //MMSEQS_FILEUTIL_H
