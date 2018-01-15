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
    if(mkdir(directoryName, mode) == 0){
        return true;
    }else{
        return false;
    }
}

void* FileUtil::mmapFile(FILE * file, size_t *dataSize){
    struct stat sb;
    if (fstat(fileno(file), &sb) < 0)
    {
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

size_t FileUtil::getFileSize(std::string fileName) {
    struct stat stat_buf;
    int rc = stat(fileName.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}


bool FileUtil::symlinkExists(const std::string  path)
{
    struct stat buf;
    int result;

    result = lstat(path.c_str(), &buf);

    return (result == 0);
}

bool FileUtil::symlinkCreateOrRepleace(const std::string linkname, const std::string linkdest) {
    if(symlinkExists(linkname)==true){
        if(remove(linkname.c_str()) != 0){
            return false;
        }
    }
    char *abs_in_header_filename = realpath(linkdest.c_str(), NULL);
    symlink(abs_in_header_filename, linkname.c_str());
    free(abs_in_header_filename);
    return true;
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
