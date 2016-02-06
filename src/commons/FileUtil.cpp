#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"
#include <sys/stat.h>

void FileUtil::errorIfFileExist(const char * file){
    struct stat st;
    if(stat(file, &st) == 0) { errno = EEXIST; perror(file); EXIT(EXIT_FAILURE); }
}

bool FileUtil::fileExists(const char* fileName) {
    struct stat st;
    return stat(fileName, &st) == 0;
}

FILE* FileUtil::openFileOrDie(const char * fileName, const char * mode, bool shouldExist) {
    bool exists = FileUtil::fileExists(fileName);
    if(exists && !shouldExist) {
        errno = EEXIST;
        perror(fileName);
        EXIT(EXIT_FAILURE);
    }
    if(!exists && shouldExist) {
        errno = ENOENT;
        perror(fileName);
        EXIT(EXIT_FAILURE);
    }

    FILE* file;
    file = fopen(fileName, mode);
    if(file == NULL) { perror(fileName); EXIT(EXIT_FAILURE); }
    return file;
}

size_t FileUtil::countLines(const char* name) {
    std::ifstream index(name);
    if (index.fail()) {
        Debug(Debug::ERROR) << "File " << name << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    size_t cnt = 0;
    std::vector<char> buffer(1024 * 1024);
    index.read(buffer.data(), buffer.size());
    while (size_t r = index.gcount()) {
        for (size_t i = 0; i < r; i++) {
            const char *p = buffer.data();
            if (p[i] == '\n') {
                cnt++;
            }
        }
        index.read(buffer.data(), buffer.size());
    }
    index.close();

    return cnt;
}
