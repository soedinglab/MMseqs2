#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"
#include "MemoryMapped.h"
#include <sys/stat.h>
#include <stdio.h>
#include <fstream>
#include <fcntl.h>
#include <unistd.h>

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
    MemoryMapped indexData(name, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    size_t cnt = 0;
    char* indexDataChar = (char *) indexData.getData();
    for(size_t pos = 0; pos < indexData.size(); pos++) {
        cnt += (indexDataChar[pos] == '\n') ? 1 : 0;
    }
    indexData.close();
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

void FileUtil::writeFile(std::string pathToFile, const unsigned char *data, size_t len) {
    int fd = open(pathToFile.c_str(), O_RDWR|O_CREAT, 0700);
    if (fd == -1) {
        Debug(Debug::ERROR) << "Could not write file " << pathToFile << "!\n";
        EXIT(EXIT_FAILURE);
    }

    ssize_t res = write(fd, data, len);
    if (res == -1) {
        Debug(Debug::ERROR) << "Error writing file " << pathToFile << "!\n";
        EXIT(EXIT_FAILURE);
    }

    if (fsync(fd) != 0) {
        Debug(Debug::ERROR) << "Error syncing file " << pathToFile << "!\n";
        EXIT(EXIT_FAILURE);
    }

    if (close(fd) != 0) {
        Debug(Debug::ERROR) << "Error closing file " << pathToFile << "!\n";
        EXIT(EXIT_FAILURE);
    }
}
