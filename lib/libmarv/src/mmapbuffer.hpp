#ifndef MMAP_BUFFER_HPP
#define MMAP_BUFFER_HPP

#include "hpc_helpers/all_helpers.cuh"

#include <sys/mman.h>
#include <unistd.h>
#include <cstdint>
#include <cstdio>
#include <cassert>
#include <string>
#include <stdexcept>
#include <utility>
#include <iostream>

#include <type_traits>

namespace cudasw4{

class FileBackedMMapBuffer{    
private:
    class OpenCFile{
    private:
        FILE* file = nullptr;

    public:
        OpenCFile() = default;

        OpenCFile(const char* filename, const char* mode){
            file = fopen(filename, mode);
            if(file == nullptr){
                perror("OpenCFile fopen");
                throw std::runtime_error("Cannot open file " + std::string(filename));
            }
        }

        OpenCFile(const OpenCFile&) = delete;
        OpenCFile(OpenCFile&& rhs){
            file = std::exchange(rhs.file, nullptr);
        }
        
        OpenCFile& operator=(OpenCFile rhs){
            std::swap(*this, rhs);
            return *this;
        }
        
        ~OpenCFile(){
            if(file != nullptr){
                fclose(file);
            }
        }

        friend void swap(OpenCFile& l, OpenCFile& r) noexcept{
            using std::swap;

            swap(l.file, r.file);
        }

        FILE* getFile() const noexcept{
            return file;
        }

        int getFd() const noexcept{
            int ret = fileno(file);
            if(ret == -1){
                perror("OpenCFile fileno");
            }
            return ret;
        }
    };

    void* rawtotaldata = nullptr;
    std::size_t size = 0;
    std::size_t capacity = 0;

    std::size_t memoryCapacity = 0;
    std::size_t fileCapacity = 0;

    std::size_t memoryLimit = 0;

    OpenCFile filehandle{};
    std::string filename{};

    //undo all mappings
    int unmapMemoryAndFile(){
        int ret = unmapFile();
        if(ret != 0){
            return ret;
        }

        if(memoryCapacity > 0){
            ret = munmap(rawtotaldata, memoryCapacity);
            if(ret == 0){
                rawtotaldata = nullptr;
            }
        };
        return ret;
    }

    //undo mappings of file
    int unmapFile(){
        if(fileCapacity > 0){
            int ret = msync(((char*)rawtotaldata) + memoryCapacity, fileCapacity, MS_SYNC);
            if(ret != 0){
                return ret;
            }

            ret = munmap(((char*)rawtotaldata) + memoryCapacity, fileCapacity);
            if(ret != 0){
                return ret;
            }
            if(memoryCapacity == 0){
                rawtotaldata = nullptr;
            }
        }
       
        return 0;
    }

    //create virtual address range of size newcapacity which is backed by anonymous mapping
    void* remap(std::size_t newcapacity){
        void* newptr = nullptr;

        if(rawtotaldata != nullptr){           
            newptr = mremap(rawtotaldata, memoryCapacity, newcapacity, MREMAP_MAYMOVE);
        }else{
            //neither memory mapping nor file mapping. make a fresh range
            newptr = mmap(0, newcapacity, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
        }

        return newptr;
    }

    //create virtual address range of size newcapacity which is backed by anonymous mapping for at most memoryLimit bytes
    //After memoryLimit bytes, the adress range is backed by a file mapping
    void remapWithFile(std::size_t newcapacity){
        //remove file mappings, if any
        int ret = unmapFile();
        if(ret != 0){
            perror("FileBackedMMapBuffer::remapWithFile unmapFile");
            throw std::runtime_error("remapWithFile failed");
        }

        #ifndef NDEBUG
        const std::size_t pagesize = getpagesize();
        assert(newcapacity % pagesize == 0);
        #endif

        //get new virtual adress range of size newcapacity with anonymous mapping
        void* newptr = remap(newcapacity);
        if(newptr == MAP_FAILED){
            perror("FileBackedMMapBuffer::reserve remap");
            throw std::runtime_error("Reserve failed");
        }

        rawtotaldata = newptr;

        memoryCapacity = std::min(newcapacity, memoryLimit);
        fileCapacity = newcapacity - memoryCapacity;

        #ifndef NDEBUG
        assert(memoryCapacity % pagesize == 0);
        assert(fileCapacity % pagesize == 0);
        #endif

        if(fileCapacity > 0){
            //update file size
            int ret = ftruncate(filehandle.getFd(), fileCapacity);
            if(ret != 0){
                perror("FileBackedMMapBuffer::reserve ftruncate");
                throw std::runtime_error("Reserve failed");
            }

            //update file mappings to virtual adress range
            //this overwrites the previous anonymous mapping.
            newptr = mmap(((char*)rawtotaldata) + memoryCapacity, fileCapacity, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_FIXED | MAP_POPULATE | MAP_NONBLOCK, filehandle.getFd(), 0);
            if(newptr == MAP_FAILED){
                perror("FileBackedMMapBuffer::reserve file mmap");
                throw std::runtime_error("Reserve failed");
            }
        }
    }

public:    

    FileBackedMMapBuffer() = delete;

    FileBackedMMapBuffer(std::size_t size_, std::size_t memoryLimit_, std::string filename_)
        : memoryLimit(memoryLimit_), filehandle(filename_.c_str(), "a+"), filename(filename_){

        const std::size_t pagesize = getpagesize();
        //round down memoryLimit to pagesize
        memoryLimit = (memoryLimit / pagesize) * pagesize;

        std::size_t numPagesForSize = SDIV(size_, pagesize);

        reserve(numPagesForSize * pagesize);

        size = size_;
    }

    FileBackedMMapBuffer(const FileBackedMMapBuffer&) = delete;
    FileBackedMMapBuffer(FileBackedMMapBuffer&& rhs){
        rawtotaldata = std::exchange(rhs.rawtotaldata, nullptr);
        size = std::exchange(rhs.size, 0);
        capacity = std::exchange(rhs.capacity, 0);
        memoryCapacity = std::exchange(rhs.memoryCapacity, 0);
        fileCapacity = std::exchange(rhs.fileCapacity, 0);
        memoryLimit = std::exchange(rhs.memoryLimit, 0);
        filehandle = std::move(rhs.filehandle);
        filename = std::exchange(rhs.filename, "");
    }

    ~FileBackedMMapBuffer(){
        int ret = unmapMemoryAndFile();
        if(ret != 0){
            perror("~FileBackedMMapBuffer unmapMemoryAndFile");
        }else{
            rawtotaldata = nullptr;
            size = 0;
            capacity = 0;
            memoryCapacity = 0;
            memoryLimit = 0;
        }

        remove(filename.c_str());
    }

    friend void swap(FileBackedMMapBuffer& l, FileBackedMMapBuffer& r) noexcept{
        using std::swap;

        std::swap(l.rawtotaldata, r.rawtotaldata);
        std::swap(l.size, r.size);
        std::swap(l.capacity, r.capacity);
        std::swap(l.memoryCapacity, r.memoryCapacity);
        std::swap(l.fileCapacity, r.fileCapacity);
        std::swap(l.memoryLimit, r.memoryLimit);
        std::swap(l.filehandle, r.filehandle);
        std::swap(l.filename, r.filename);
    }

    void destroy(){
        int ret = unmapMemoryAndFile();
        if(ret != 0){
            perror("FileBackedMMapBuffer::destroy unmapMemoryAndFile");
        }else{
            rawtotaldata = nullptr;
            size = 0;
            capacity = 0;
            memoryCapacity = 0;
            memoryLimit = 0;
        }
    }

    void reserve(std::size_t newcapacity){
        if(newcapacity > capacity){
            const std::size_t pagesize = getpagesize();
            newcapacity = SDIV(newcapacity, pagesize) * pagesize;
            remapWithFile(newcapacity);
            capacity = newcapacity;
        }        
    }

    void resize(std::size_t newsize){        
        if(newsize <= capacity){
            size = newsize;
        }else{
            reserve(newsize);
            size = newsize;
        }        
    }

    void shrink_to_fit(){
        if(size == 0){
            int ret = unmapMemoryAndFile();
            if(ret != 0){
                perror("FileBackedMMapBuffer::shrink_to_fit unmapMemoryAndFile");
                throw std::runtime_error("shrink_to_fit failed");
            }else{
                rawtotaldata = nullptr;
                size = 0;
                capacity = 0;
                memoryCapacity = 0;
                memoryLimit = 0;
            }
        }else{
            if(size < capacity){
                const std::size_t pagesize = getpagesize();
                std::size_t newcapacity = SDIV(size, pagesize) * pagesize;
                remapWithFile(newcapacity);
                capacity = newcapacity;
            }
        }
    }

    void clear(){
        size = 0;
    }

    void* get() noexcept{
        return rawtotaldata;
    }

    const void* get() const noexcept{
        return rawtotaldata;
    }

    std::size_t getSize() const noexcept{
        return size;
    }
    
    std::size_t getCapacity() const noexcept{
        return capacity;
    }

    std::size_t getCapacityInMemory() const noexcept{
        return memoryCapacity;
    }

    std::size_t getCapacityInFile() const noexcept{
        return fileCapacity;
    }

    void printStatus(std::ostream& os) const{
        os << "size: " << getSize() << ", capacity: " << getCapacity();
        os << ", memoryCapacity: " << getCapacityInMemory() << ", fileCapacity: " << getCapacityInFile();
        os << ", memoryLimit: " << memoryLimit;
    }
};


template<class T, int growth = 150>
class FileBackedUVector{
    static_assert(std::is_trivial<T>::value, "FileBackedUVector: T must be trivial");
    static_assert(growth >= 100, "growth must be >= 100");

private:
    FileBackedMMapBuffer buffer;
    std::size_t size_ = 0;
    std::size_t capacity_ = 0;

    void grow(std::size_t mincapacity){
        constexpr double growthFactor = double(growth) / 100.0;

        const size_t growthCapacity = std::max(capacity() + 1, size_t(capacity() * growthFactor));
        if(mincapacity <= growthCapacity){
            buffer.resize(sizeof(T) * growthCapacity);
            capacity_ = growthCapacity;
        }else{
            buffer.resize(sizeof(T) * mincapacity);
            capacity_ = mincapacity;
        }
        
    }
public:
    using value_type = T;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using iterator = value_type*;
    using const_iterator = const value_type*;
    using reference = value_type&;
    using const_reference = const value_type&;


    FileBackedUVector(std::size_t elements, std::size_t maxBytesInMemory, const std::string& backingfile)
        : buffer(sizeof(T) * elements, maxBytesInMemory, backingfile){

        size_ = elements;
        capacity_ = elements;
    }

    void push_back(T obj){
        if(size() >= capacity()){
            reserve(size() + 1);
        }

        data()[size_++] = std::move(obj);
    }

    template<class InputIt>
    iterator insert(const_iterator pos, InputIt first, InputIt last){
        const std::size_t insertsize = std::distance(first, last);
        const std::size_t newsize = size() + insertsize;
        const std::size_t where = std::distance(cbegin(), pos);
        assert(where <= size());

        if(newsize > capacity()){
            reserve(newsize);
        }

        std::copy_backward(cbegin() + where, cend(), begin() + newsize);
        std::copy(first, last, begin() + where);

        size_ = newsize;

        return begin() + where;
    }

    iterator insert(const_iterator pos, size_t count, const T& value){
        const std::size_t insertsize = count;
        const std::size_t newsize = size() + insertsize;
        const std::size_t where = std::distance(cbegin(), pos);
        assert(where <= size());

        if(newsize > capacity()){
            reserve(newsize);
        }

        std::copy_backward(cbegin() + where, cend(), begin() + newsize);
        std::fill(begin() + where, begin() + where + insertsize, value);

        size_ = newsize;

        return begin() + where;
    }

    std::size_t size() const noexcept{
        return size_;
    }

    std::size_t capacity() const noexcept{
        return capacity_;
    }


    pointer data() noexcept{
        return reinterpret_cast<pointer>(buffer.get());
    }

    const_pointer data() const noexcept{
        return reinterpret_cast<const_pointer>(buffer.get());
    }

    iterator begin() noexcept{
        return data();
    }

    const_iterator cbegin() const noexcept{
        return data();
    }

    iterator end() noexcept{
        return data() + size();
    }

    const_iterator cend() const noexcept{
        return data() + size();
    }    

    bool empty() const noexcept{
        return size() == 0;
    }

    reference operator[](size_t i){
        return data()[i];
    }

    const_reference operator[](size_t i) const{
        return data()[i];
    }

    reference front(){
        return data()[0];
    }

    const_reference front() const{
        return data()[0];
    }

    reference back(){
        return data()[size() - 1];
    }

    const_reference back() const{
        return data()[size() - 1];
    }

    std::size_t getCapacityInMemoryInBytes() const noexcept{
        return buffer.getCapacityInMemory();
    }

    std::size_t getCapacityInFileInBytes() const noexcept{
        return buffer.getCapacityInFile();
    }

    void clear(){
        buffer.clear();
    }

    void resize(std::size_t newsize){
        if(newsize > capacity()){
            grow(newsize);
        }
        size_ = newsize;
    }

    void reserve(std::size_t newcapacity){
        if(newcapacity > capacity()){
            grow(newcapacity);
        }
    }
};

} //namespace cudasw4



#endif