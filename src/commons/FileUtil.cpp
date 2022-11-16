#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"

#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/simde-common.h>

#include <cstddef>
#include <cstring>
#include <climits>
#include <algorithm>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

#include <sys/stat.h>
#include <sys/statvfs.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/resource.h>

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

void FileUtil::munmapData(void * ptr, size_t dataSize){
    if(munmap(ptr, dataSize) < 0){
        Debug(Debug::ERROR) << "Failed to munmap memory\n";
        EXIT(EXIT_FAILURE);
    }
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
    FILE *fp = FileUtil::openFileOrDie(name, "r", true);
    size_t cnt = 0;
    while (!feof(fp)) {
        char ch = fgetc(fp);
        cnt += (ch == '\n') ? 1 : 0;
    }
    if (fclose(fp) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << name << "\n";
        EXIT(EXIT_FAILURE);
    }
    return cnt;
}


void FileUtil::deleteTempFiles(const std::list<std::string> &tmpFiles) {
    for (std::list<std::string>::const_iterator it = tmpFiles.begin(); it != tmpFiles.end(); it++) {
        Debug(Debug::INFO) << "Deleting " << *it << "\n";
        std::string file = *it;
        FileUtil::remove(file.c_str());
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
        return SIZE_MAX;
    }

    // the available size is f_bsize * f_bavail
    return stat.f_bfree * stat.f_frsize;
}


std::string FileUtil::getRealPathFromSymLink(const std::string path){
    char *p = realpath(path.c_str(), NULL);
    if (p == NULL) {
        Debug(Debug::ERROR) << "Could not get path of " << path << "!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string name(p);
    free(p);
    return name;
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

    std::string pathToAlias = (path + "/" + alias);
    if (symlinkExists(pathToAlias) == true){
        FileUtil::remove(pathToAlias.c_str());
    }
    // symlinkat is not available in Conda macOS
    // Conda uses the macOS 10.9 SDK, and symlinkat was introduced in 10.10
    // We emulate symlinkat by manipulating the CWD instead
    std::string oldWd = FileUtil::getCurrentWorkingDirectory();
    if (chdir(path.c_str()) != 0) {
        Debug(Debug::ERROR) << "Could not change working directory to " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (symlink(base.c_str(), alias.c_str()) != 0) {
        Debug(Debug::ERROR) << "Could not create symlink of " << file << "!\n";
        EXIT(EXIT_FAILURE);
    }
    if (chdir(oldWd.c_str()) != 0) {
        Debug(Debug::ERROR) << "Could not change working directory to " << oldWd << "\n";
        EXIT(EXIT_FAILURE);
    }
}

std::string FileUtil::getCurrentWorkingDirectory() {
    // CWD can be larger than PATH_MAX and allocating enough memory is somewhat tricky
    char* wd = NULL;
#ifdef PATH_MAX
    size_t bufferSize = PATH_MAX;
#else
    size_t bufferSize = 1024;
#endif
    do {
        if (wd != NULL) {
            free(wd);
            bufferSize *= 2;
        }
        errno = 0;
        wd = getcwd(NULL, bufferSize);
        if (wd == NULL && errno != ERANGE && errno != 0) {
            Debug(Debug::ERROR) << "Could not get current working directory\n";
            EXIT(EXIT_FAILURE);
        }
    } while (wd == NULL && errno == ERANGE);
    std::string cwd(wd);
    free(wd);
    return cwd;
}

void FileUtil::symlinkAbs(const std::string &target, const std::string &link) {
    if(FileUtil::fileExists(link.c_str())){
        FileUtil::remove(link.c_str());
    }
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
        if (symlinkExists(realLink) == true){
            FileUtil::remove(realLink.c_str());
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
        size_t res = write(dest, buf, size);
        if (res != size) {
            Debug(Debug::ERROR) << "Error writing file " << dst << "!\n";
            EXIT(EXIT_FAILURE);
        }
    }
    close(source);
    close(dest);
}

void FileUtil::copyFile(const std::string& src, const std::string& dst) {
    copyFile(src.c_str(), dst.c_str());
}

FILE * FileUtil::openAndDelete(const char *fileName, const char *mode) {
    if(FileUtil::fileExists(fileName) == true){
        if(FileUtil::directoryExists(fileName)){
            Debug(Debug::ERROR) << "Can not open " << fileName << " for writing. It is a directory.\n";
            EXIT(EXIT_FAILURE);
        }else {
            FileUtil::remove(fileName);
        }
    }
    FILE * file = fopen(fileName, mode);
    if (file == NULL) {
        Debug(Debug::ERROR) << "Could not open " << fileName << " for writing!\n";
        EXIT(EXIT_FAILURE);
    }
    return file;
}

std::vector<std::string> FileUtil::findDatafiles(const char * datafiles){
    std::string baseName = std::string(datafiles);
    std::string checkName = baseName + ".0";
    std::vector<std::string> filenames;
    size_t cnt = 0;
    while(FileUtil::fileExists(checkName.c_str()) == true){
        filenames.push_back(checkName);
        cnt++;
        checkName = baseName + "." + SSTR(cnt);
    }
    if(cnt == 0){
        if(FileUtil::fileExists(baseName.c_str())){
            filenames.push_back(baseName);
        }
    }
    return filenames;
}

void FileUtil::remove(const char * file ) {
    if (std::remove(file) != 0){
        Debug(Debug::ERROR) << "Could not delete " << file << "!\n";
        EXIT(EXIT_FAILURE);
    }
}

void FileUtil::move(const char * src, const char * dst) {
    struct stat srcFileInfo;
    FILE * srcFile = FileUtil::openFileOrDie(src, "rw", true);
    if (fstat(fileno(srcFile), &srcFileInfo) < 0) {
        int errsv = errno;
        Debug(Debug::ERROR) << "Failed to fstat File=" << src << ". Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }
    struct stat srcDirInfo;
    std::string dirName = FileUtil::dirName(dst);
    FILE * dstDir = FileUtil::openFileOrDie(dirName.c_str(), "r", true);
    if (fstat(fileno(dstDir), &srcDirInfo) < 0) {
        int errsv = errno;
        Debug(Debug::ERROR) << "Failed to fstat File=" << dirName << ". Error " << errsv << ".\n";
        EXIT(EXIT_FAILURE);
    }
    bool sameFileSystem = (srcDirInfo.st_dev == srcFileInfo.st_dev);
    if (fclose(srcFile) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << src << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(dstDir) != 0) {
        Debug(Debug::ERROR) << "Cannot close directory " << dirName << "\n";
        EXIT(EXIT_FAILURE);
    }
    if(sameFileSystem){
        if (std::rename(src, dst) != 0){
            Debug(Debug::ERROR) << "Could not copy file " << src << " to " << dst << "!\n";
            EXIT(EXIT_FAILURE);
        }
    }else{
        FileUtil::copyFile(src, dst);
        FileUtil::remove(src);
    }
}

int FileUtil::parseDbType(const char *name) {
    std::string dbTypeFile = std::string(name) + ".dbtype";
    if (FileUtil::fileExists(dbTypeFile.c_str()) == false) {
        return Parameters::DBTYPE_GENERIC_DB;
    }

    size_t fileSize = FileUtil::getFileSize(dbTypeFile);
    if (fileSize != sizeof(int)) {
        Debug(Debug::ERROR) << "File size of " << dbTypeFile << " seems to be wrong!\n";
        Debug(Debug::ERROR) << "It should have 4 bytes but it has " <<  fileSize << " bytes.";
        EXIT(EXIT_FAILURE);
    }
    FILE *file = fopen(dbTypeFile.c_str(), "r");
    if (file == NULL) {
        Debug(Debug::ERROR) << "Could not open data file " << dbTypeFile << "!\n";
        EXIT(EXIT_FAILURE);
    }
    int dbtype;
    size_t result = fread(&dbtype, 1, fileSize, file);
    if (result != fileSize) {
        Debug(Debug::ERROR) << "Could not read " << dbTypeFile << "!\n";
        EXIT(EXIT_FAILURE);
    }
    if (fclose(file) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << dbTypeFile << "\n";
        EXIT(EXIT_FAILURE);
    }
#if SIMDE_ENDIAN_ORDER == SIMDE_ENDIAN_BIG
    dbtype = __builtin_bswap32(dbtype);
#endif
    return dbtype;
}

std::string FileUtil::createTemporaryDirectory(const std::string& basePath, const std::string& subDirectory) {
    std::string tmpDir(basePath);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        Debug(Debug::INFO) << "Path " << tmpDir << " does not exist or is not a directory.\n";
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Cannot create temporary folder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        } else {
            Debug(Debug::INFO) << "Created directory " << tmpDir << "\n";
        }
    }
    tmpDir += "/" + subDirectory;
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Cannot create temporary subfolder " << tmpDir << ".\n";
            EXIT(EXIT_FAILURE);
        }
    }
    FileUtil::symlinkAlias(tmpDir, "latest");
    return tmpDir;
}

void FileUtil::fixRlimitNoFile() {
    static bool increasedRlimitNoFile(false);
    if (increasedRlimitNoFile == false) {
        increasedRlimitNoFile = true;
        struct rlimit limit;
        if (getrlimit(RLIMIT_NOFILE, &limit) != 0) {
            Debug(Debug::WARNING) << "Could not increase maximum number of open files (getrlimit " << errno << "). Use ulimit manually\n";
            return;
        }
        limit.rlim_cur = std::min(std::max((rlim_t)8192, limit.rlim_cur), limit.rlim_max);
        limit.rlim_max = std::min(RLIM_INFINITY, limit.rlim_max);
        if (setrlimit(RLIMIT_NOFILE, &limit) != 0) {
            Debug(Debug::WARNING) << "Could not increase maximum number of open files (setrlimit " << errno << "). Use ulimit manually\n";
        }
    }
}

std::string FileUtil::sanitizeFilename(std::string name){
    static const std::vector<std::pair<char, char>> symbolTable =
            {{'\\', '@'},
             {'/', '@'},
             {':', '@'},
             {'*', '@'},
             {'?', '@'},
             {'<', '@'},
             {'>', '@'},
             {'|', '!'}};

    std::vector<std::pair<char, char>>::const_iterator it;
    for (it = symbolTable.begin(); it != symbolTable.end(); ++it) {
        std::replace(name.begin(), name.end(), it->first, it->second);
    }
    return name;
}
