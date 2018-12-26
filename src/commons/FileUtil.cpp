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
#include <sys/mman.h>

bool FileUtil::fileExists(const char* fileName) {
    struct stat st;
    return stat(fileName, &st) == 0;
}

bool FileUtil::directoryExists(const char* directoryName) {
    struct stat st;
    return stat(directoryName, &st) == 0 && S_ISDIR(st.st_mode);
}

bool FileUtil::makeDir(const char* directoryName, const int mode ) {
    return mkdir(directoryName, mode) == 0;
}

void* FileUtil::mmapFile(FILE * file, size_t *dataSize){
    struct stat sb;
    if (fstat(fileno(file), &sb) < 0) {
        int errsv = errno;
        Debug(Debug::ERROR) << "Failed to fstat." << ". Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }
    *dataSize = sb.st_size;
    int mode = PROT_READ;
    int fd =  fileno(file);

    void * ret = mmap(NULL, *dataSize, mode, MAP_PRIVATE, fd, 0);
    if(ret == MAP_FAILED){
        int errsv = errno;
        Debug(Debug::ERROR) << "Failed to mmap memory dataSize=" << *dataSize <<". Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }
    return ret;
}

FILE* FileUtil::openFileOrDie(const char * fileName, const char * mode, bool shouldExist) {
    bool exists = FileUtil::fileExists(fileName);
    if (exists && !shouldExist) {
        errno = EEXIST;
        perror(fileName);
        EXIT(EXIT_FAILURE);
    }
    if (!exists && shouldExist) {
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
    for (size_t pos = 0; pos < indexData.size(); pos++) {
        cnt += (indexDataChar[pos] == '\n') ? 1 : 0;
    }
    indexData.close();
    return cnt;
}

void FileUtil::deleteFile(const std::string &file) {
    if (remove(file.c_str()) != 0) {
        Debug(Debug::WARNING) << "Error deleting file " << file << "\n";
    }
}

void FileUtil::deleteTempFiles(const std::list<std::string> &tmpFiles) {
    for (std::list<std::string>::const_iterator it = tmpFiles.begin(); it != tmpFiles.end(); it++) {
        Debug(Debug::INFO) << "Deleting " << *it << "\n";
        deleteFile(*it);
    }
}

void FileUtil::writeFile(const std::string &pathToFile, const unsigned char *data, size_t len) {
    int fd = open(pathToFile.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IXUSR);
    if (fd == -1) {
        Debug(Debug::ERROR) << "Could not write file " << pathToFile << "!\n";
        EXIT(EXIT_FAILURE);
    }

    ssize_t res = write(fd, data, len);
    if (res == -1) {
        Debug(Debug::ERROR) << "Error writing file " << pathToFile << "!\n";
        EXIT(EXIT_FAILURE);
    }

    if (close(fd) != 0) {
        Debug(Debug::ERROR) << "Error closing file " << pathToFile << "!\n";
        EXIT(EXIT_FAILURE);
    }
}

std::string FileUtil::dirName(const std::string &file) {
        size_t pos = file.find_last_of("\\/");
        return (std::string::npos == pos)
               ? "."
               : file.substr(0, pos);
}

std::string FileUtil::baseName(const std::string &file) {
    size_t pos = file.find_last_of("\\/");
    return (std::string::npos == pos)
           ? file
           : file.substr(pos+1, file.length());
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

std::string FileUtil::getHashFromSymLink(const std::string path){
    char *p = realpath(path.c_str(), NULL);
    if (p == NULL) {
        Debug(Debug::ERROR) << "Could not get path of " << path << "!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string base = baseName(p);
    free(p);

    return base;
}

void FileUtil::symlinkAlias(const std::string &file, const std::string &alias) {
    char *p = realpath(file.c_str(), NULL);
    if (p == NULL) {
        Debug(Debug::ERROR) << "Could not get path of " << file << "!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string path = dirName(p);
    std::string base = baseName(p);
    free(p);

    DIR *dir = opendir(path.c_str());
    if (dir == NULL) {
        Debug(Debug::ERROR) << "Error opening directory " << path << "!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string pathToAlias = (path + "/" + alias);
    if (symlinkExists(pathToAlias) == true && remove(pathToAlias.c_str()) != 0){
        Debug(Debug::ERROR) << "Could not remove old symlink " << pathToAlias << "!\n";
        EXIT(EXIT_FAILURE);
    }

    if (symlinkat(base.c_str(), dirfd(dir), alias.c_str()) != 0) {
        Debug(Debug::ERROR) << "Could not create symlink of " << file << "!\n";
        EXIT(EXIT_FAILURE);
    }

    if (closedir(dir) != 0) {
        Debug(Debug::ERROR) << "Error closing directory " << path << "!\n";
        EXIT(EXIT_FAILURE);
    }
}

void FileUtil::symlinkAbs(const std::string &target, const std::string &link) {
    char *t = realpath(target.c_str(), NULL);
    if (t == NULL) {
        Debug(Debug::ERROR) << "Could not get realpath of " << target << "!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string realLink;
    char *l = realpath(link.c_str(), NULL);
    if (l == NULL) {
        std::string path = dirName(link);
        std::string base = baseName(link);
        l = realpath(path.c_str(), NULL);
        if (l == NULL) {
            Debug(Debug::ERROR) << "Could not get realpath of " << link << "!\n";
            EXIT(EXIT_FAILURE);
        } else {
            realLink = (std::string(l) + "/" + base);
        }
    } else {
        realLink = l;
        if (symlinkExists(realLink) == true && remove(realLink.c_str()) != 0){
            Debug(Debug::ERROR) << "Could not remove old symlink " << link << "!\n";
            EXIT(EXIT_FAILURE);
        }
    }

    if (symlink(t, realLink.c_str()) != 0) {
        Debug(Debug::ERROR) << "Could not create symlink of " << target << "!\n";
        EXIT(EXIT_FAILURE);
    }

    free(t);
    free(l);
}

size_t FileUtil::getFileSize(const std::string &fileName) {
    struct stat stat_buf;
    int rc = stat(fileName.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}


bool FileUtil::symlinkExists(const std::string &path)  {
    struct stat buf;
    int result = lstat(path.c_str(), &buf);
    return (result == 0);
}

void FileUtil::copyFile(const char *src, const char *dst) {
    //https://stackoverflow.com/questions/10195343/copy-a-file-in-a-sane-safe-and-efficient-way
    char buf[BUFSIZ];
    size_t size;

    int source = open(src, O_RDONLY, 0);
    if (source == -1) {
        Debug(Debug::ERROR) << "Could not open file " << src << "!\n";
        EXIT(EXIT_FAILURE);
    }
    int dest = open(dst, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (dest == -1) {
        Debug(Debug::ERROR) << "Could not open file " << dst << "!\n";
        EXIT(EXIT_FAILURE);
    }
    while ((size = read(source, buf, BUFSIZ)) > 0) {
        write(dest, buf, size);
    }
    close(source);
    close(dest);
}

FILE * FileUtil::openAndDelete(const char *fileName, const char *mode) {
    if(FileUtil::fileExists(fileName) == true){
        FileUtil::deleteFile(fileName);
    }
    FILE * file = fopen(fileName, mode);
    if (file == NULL) {
        Debug(Debug::ERROR) << "Could not open " << fileName << " for writing!\n";
        EXIT(EXIT_FAILURE);
    }
    return file;
}
