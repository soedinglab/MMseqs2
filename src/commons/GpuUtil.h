#ifndef GPUUTIL_H
#define GPUUTIL_H

#include "Debug.h"
#include "FileUtil.h"
#include "PrefilteringIndexReader.h"
#include "marv.h"
#include <atomic>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <string>
#include <sys/mman.h>
#include <unistd.h>
#include <utility>
#include <thread>

struct GPUSharedMemory {
    enum State {
        IDLE,
        RESERVED,
        READY,
        DONE
    };

    unsigned int maxSeqLen;                   // Maximum length of the sequence
    unsigned int maxResListLen;               // Maximum length of the results list
    std::atomic<int>  state{IDLE};            // State of the shared memory
    std::atomic<bool> serverExit{false};      // Has server exited
    unsigned int queryOffset;                 // Offset to the query data
    unsigned int queryLen;                    // Length of the query sequence
    unsigned int resultsOffset;               // Offset to the results data
    unsigned int resultLen;                   // Length of the result list
    unsigned int profileOffset;               // Offset to the profile data

    // Get pointers to the query, results, and profile data sections
    int8_t* getQueryPtr() { return reinterpret_cast<int8_t*>(this) + queryOffset; }
    Marv::Result* getResultsPtr() { return reinterpret_cast<Marv::Result*>(reinterpret_cast<char*>(this) + resultsOffset); }
    int8_t* getProfilePtr() { return reinterpret_cast<int8_t*>(this) + profileOffset; }

    // Calculate the total size needed for the shared memory
    static size_t calculateSize(unsigned int maxSeqLen, unsigned int maxResListLen) {
        return sizeof(GPUSharedMemory) +
               sizeof(char) * maxSeqLen +                              // Size for query data
               sizeof(Marv::Result) * maxResListLen +  // Size for results data
               sizeof(int8_t) * 21 * maxSeqLen;        // Size for profile data
    }

    static std::string getShmHash(const std::string& db) {
        std::string dbpath = FileUtil::getRealPathFromSymLink(PrefilteringIndexReader::dbPathWithoutIndex(db));
        char* visibleDevices = getenv("CUDA_VISIBLE_DEVICES");
        if (visibleDevices) {
            dbpath.append(visibleDevices);
        }
        size_t hash = Util::hash(dbpath.c_str(), dbpath.length());
        return SSTR(hash);
    }

    // Allocate and initialize shared memory
    static GPUSharedMemory* alloc(const std::string& name, unsigned int maxSeqLen, unsigned int maxResListLen) {
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
    static void dealloc(GPUSharedMemory* layout, const std::string& name) {
        if (layout) {
            size_t shm_size = calculateSize(layout->maxSeqLen, layout->maxResListLen);
            if (munmap(layout, shm_size) == -1) {
                Debug(Debug::ERROR) << "Error unmapping shared memory\n";
            }
            if (shm_unlink(name.c_str()) == -1) {
                Debug(Debug::ERROR) << "Error unlinking shared memory\n";
            }
        }
    }

    static void unmap(GPUSharedMemory* layout) {
        if (layout) {
            size_t shm_size = calculateSize(layout->maxSeqLen, layout->maxResListLen);
            if (munmap(layout, shm_size) == -1) {
                Debug(Debug::ERROR) << "Error unmapping shared memory\n";
            }
        }
    }

    // Function to open and map existing shared memory and automatically determine sizes
    static GPUSharedMemory* openSharedMemory(const std::string& name) {
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

};

#endif
