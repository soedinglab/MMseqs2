#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"
#include <sys/stat.h>
#include <stdio.h>
#include <fstream>

void FileUtil::errorIfFileExist(const char * file){
    struct stat st;
    if(stat(file, &st) == 0) { errno = EEXIST; perror(file); EXIT(EXIT_FAILURE); }
}

bool FileUtil::fileExists(const char* fileName) {
    struct stat st;
    return stat(fileName, &st) == 0;
}

bool FileUtil::directoryExists(const char* directoryName) {
    struct stat st;
    return stat(directoryName, &st) == 0 && S_ISDIR(st.st_mode);
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
    FILE *file1 = fopen(name, "r");
#if HAVE_POSIX_FADVISE
    if (posix_fadvise (fileno(file1), 0, 0, POSIX_FADV_SEQUENTIAL) != 0){
       Debug(Debug::ERROR) << "posix_fadvise returned an error\n";
    }
#endif
    size_t cnt = 0;
    int c1;
    while((c1=getc_unlocked(file1)) != EOF) {
        cnt += (c1 == '\n') ? 1 : 0;
    }
    fclose(file1);
    return cnt;
}

void FileUtil::deleteTempFiles(std::list<std::string> tmpFiles) {
    for (std::list<std::string>::const_iterator it = tmpFiles.begin(); it != tmpFiles.end(); it++) {
        Debug(Debug::INFO) << "Deleting " << *it << "\n";
        if (remove((*it).c_str()) != 0) {
            Debug(Debug::WARNING) << "Error deleting file " << *it << "\n";
        }
    }
}
