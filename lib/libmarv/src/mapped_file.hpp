#ifndef MAPPED_FILE_HPP
#define MAPPED_FILE_HPP

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <string>
#include <stdexcept>
#include <future>
#include <chrono>
#include <iostream>
#include <cassert>

namespace cudasw4{

class MappedFileException : public std::exception{
    std::string message;
public:
    MappedFileException() : MappedFileException("MMapException"){}
    MappedFileException(const std::string& msg) : message(msg){}

    const char* what() const noexcept override{
        return message.c_str();
    }
};

struct MappedFile{
    struct Options{
        bool readaccess = true;
        bool writeaccess = false;
        bool prefault = false;
    };
    
    int fd;
    size_t filesize;
    void* rawMmapPtr;
    std::string filename;
    Options options;

    MappedFile(const std::string& filename_, MappedFile::Options options_)
        :   filesize(getFileSizeInBytes(filename_)), 
            filename(filename_),
            options(options_){

        int openflags = 0;
        if(options.readaccess && options.writeaccess){
            openflags = O_RDWR;
        }else if(options.readaccess && !options.writeaccess){
            openflags = O_RDONLY;
        }else if(!options.readaccess && options.writeaccess){
            openflags = O_WRONLY;
        }else{
            throw MappedFileException("Invalid options for MappedFile");
        }

        fd = open(filename.c_str(), openflags);
        if(fd == -1){
            perror("open");
            throw MappedFileException("Could not open file " + filename);
        }

        int mmapprot = 0;
        if(options.readaccess){
            mmapprot |= PROT_READ;
        }
        if(options.writeaccess){
            mmapprot |= PROT_WRITE; 
        }

        int mmapflags = MAP_PRIVATE;
        if(options.writeaccess){
            mmapflags = MAP_SHARED;
        }
        if(options.prefault){
            mmapflags |= MAP_POPULATE; //load the file into memory immediately
        }

        rawMmapPtr = mmap(nullptr, filesize, mmapprot, mmapflags, fd, 0);
        if(rawMmapPtr == MAP_FAILED){
            close(fd);
            throw MappedFileException("Could not map file " + filename);
        }
    }

    ~MappedFile(){
        munmap(rawMmapPtr, filesize);
        close(fd);
    }

    char* data() noexcept{
        return (char*)rawMmapPtr;
    }
    const char* data() const noexcept{
        return (const char*)rawMmapPtr;
    }
    size_t size() const noexcept{
        return filesize;
    }

    template<class T>
    size_t numElements() const noexcept{
        return size() / sizeof(T);
    }
private:
    size_t getFileSizeInBytes(const std::string& filename){
        struct stat stat_buf;
        int rc = stat(filename.c_str(), &stat_buf);
        if(rc == 0){
            return stat_buf.st_size;
        }else{
            throw MappedFileException("Could not determine file size of file " + filename);
        }
    }
};


} //namespace cudasw4

#endif