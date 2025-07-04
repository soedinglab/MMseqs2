#ifdef HAVE_CUDA
#include "GpuUtil.h"

#include "Debug.h"
#include "FileUtil.h"
#include "PrefilteringIndexReader.h"

#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <sys/mman.h>
#include <unistd.h>
#include <utility>
#include <thread>

extern const char* version;

std::string GPUSharedMemory::getShmHash(const std::string& db) {
    std::string dbpath = FileUtil::getRealPathFromSymLink(PrefilteringIndexReader::dbPathWithoutIndex(db));
    char* visibleDevices = getenv("CUDA_VISIBLE_DEVICES");
    if (visibleDevices) {
        dbpath.append(visibleDevices);
    }
    dbpath.append(version);
    size_t hash = Util::hash(dbpath.c_str(), dbpath.length());
    return SSTR(hash);
}

// Allocate and initialize shared memory
GPUSharedMemory* GPUSharedMemory::alloc(const std::string& name, unsigned int maxSeqLen, unsigned int maxResListLen) {
    size_t shm_size = calculateSize(maxSeqLen, maxResListLen);
    int fd = shm_open(name.c_str(), O_CREAT | O_RDWR, 0666);
    if (fd == -1) {
        Debug(Debug::ERROR) << "Failed to open shared memory\n";
        EXIT(EXIT_FAILURE);
    }
    if (ftruncate(fd, shm_size) == -1) {
        close(fd);
        Debug(Debug::ERROR) << "Failed to size shared memory\n";
        EXIT(EXIT_FAILURE);
    }
    void* ptr = mmap(0, shm_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);  // Close the file descriptor early as it's no longer needed after mmap
    if (ptr == MAP_FAILED) {
        Debug(Debug::ERROR) << "Failed to map shared memory\n";
        EXIT(EXIT_FAILURE);
    }

    GPUSharedMemory* layout = new (ptr) GPUSharedMemory;
    layout->maxSeqLen = maxSeqLen;
    layout->maxResListLen = maxResListLen;
    layout->queryOffset = sizeof(GPUSharedMemory);
    layout->resultsOffset = layout->queryOffset + sizeof(char) * maxSeqLen;
    layout->profileOffset = layout->resultsOffset + sizeof(Marv::Result) * maxResListLen;
    return layout;
}

// Deallocate shared memory
void GPUSharedMemory::dealloc(GPUSharedMemory* layout, const std::string& name) {
    if (layout) {
        size_t shm_size = GPUSharedMemory::calculateSize(layout->maxSeqLen, layout->maxResListLen);
        if (munmap(layout, shm_size) == -1) {
            Debug(Debug::ERROR) << "Error unmapping shared memory\n";
        }
        if (shm_unlink(name.c_str()) == -1) {
            Debug(Debug::ERROR) << "Error unlinking shared memory\n";
        }
    }
}

void GPUSharedMemory::unmap(GPUSharedMemory* layout) {
    if (layout) {
        size_t shm_size = GPUSharedMemory::calculateSize(layout->maxSeqLen, layout->maxResListLen);
        if (munmap(layout, shm_size) == -1) {
            Debug(Debug::ERROR) << "Error unmapping shared memory\n";
        }
    }
}

// Function to open and map existing shared memory and automatically determine sizes
GPUSharedMemory* GPUSharedMemory::openSharedMemory(const std::string& name) {
    int fd = shm_open(name.c_str(), O_RDWR, 0666);
    if (fd == -1) {
        Debug(Debug::ERROR) << "Failed to open shared memory\n";
        EXIT(EXIT_FAILURE);
    }

    // Map enough memory to access the first part of the structure
    void* ptr = mmap(0, sizeof(GPUSharedMemory), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (ptr == MAP_FAILED) {
        close(fd);
        Debug(Debug::ERROR) << "Failed to map shared memory\n";
        EXIT(EXIT_FAILURE);
    }

    // Now read maxSeqLen and maxResListLen from the mapped memory
    unsigned int maxSeqLen = *(reinterpret_cast<unsigned int*>(ptr));
    unsigned int maxResListLen = *(reinterpret_cast<unsigned int*>(ptr) + 1);

    // Correctly calculate the total size of the shared memory using read values
    size_t shm_size = GPUSharedMemory::calculateSize(maxSeqLen, maxResListLen);

    // Re-map with the full size now that we know it
    munmap(ptr, sizeof(GPUSharedMemory));  // Unmap the initial small mapping
    ptr = mmap(0, shm_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);  // Close the file descriptor as it's no longer needed after mmap
    if (ptr == MAP_FAILED) {
        Debug(Debug::ERROR) << "Failed to remap shared memory\n";
        EXIT(EXIT_FAILURE);
    }
    return reinterpret_cast<GPUSharedMemory*>(ptr);
}

#endif
