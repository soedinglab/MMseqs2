#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"
#include "MemoryMapped.h"
#include <sys/stat.h>
#include <stdio.h>
#include <fstream>
#include <fcntl.h>
#include <unistd.h>
#include <sys/statvfs.h>
#include <sys/types.h>
#include <dirent.h>

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

void FileUtil::deleteFile(std::string tmpFiles) {
    if (remove(tmpFiles.c_str()) != 0) {
        Debug(Debug::WARNING) << "Error deleting file " << tmpFiles << "\n";
    }
}

void FileUtil::deleteTempFiles(std::list<std::string> tmpFiles) {
    for (std::list<std::string>::const_iterator it = tmpFiles.begin(); it != tmpFiles.end(); it++) {
        Debug(Debug::INFO) << "Deleting " << *it << "\n";
        deleteFile(*it);
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

std::string FileUtil::dirName(const std::string &fileName) {
        size_t pos = fileName.find_last_of("\\/");
        return (std::string::npos == pos)
               ? ""
               : fileName.substr(0, pos);
}

size_t FileUtil::getFreeSpace(const char *path) {
        struct statvfs stat;
        if (statvfs(path, &stat) != 0) {
            // error happens, just quits here
            return -1;
        }

        // the available size is f_bsize * f_bavail
        return stat.f_bfree * stat.f_frsize;
}

void FileUtil::symlinkAlias(const std::string &file, const std::string &alias) {
    char *p = realpath(file.c_str(), NULL);
    std::string path = dirName(p);
    free(p);

    DIR *dir = opendir(path.c_str());
    if (dir == NULL) {
        Debug(Debug::ERROR) << "Error opening directory " << path << "!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string pathToAlias = (path + "/" + alias);
    if (FileUtil::fileExists(pathToAlias.c_str())) {
        remove(pathToAlias.c_str());
    }

    if (symlinkat(file.c_str(), dirfd(dir), alias.c_str()) != 0) {
        Debug(Debug::ERROR) << "Could not create symlink of " << file << "!\n";
        EXIT(EXIT_FAILURE);
    }

    if (closedir(dir) != 0) {
        Debug(Debug::ERROR) << "Error closing directory " << path << "!\n";
        EXIT(EXIT_FAILURE);
    }
}
