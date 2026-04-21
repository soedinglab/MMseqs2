#ifndef GPUUTIL_H
#define GPUUTIL_H

#include <atomic>
#include <string>
#include <cstddef>
#include <sched.h>
#ifdef USE_GPU_SEM
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <semaphore.h>
#include <fcntl.h>
#endif
#include "marv.h"

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
               sizeof(char) * maxSeqLen +              // Size for query data
               sizeof(Marv::Result) * maxResListLen +  // Size for results data
               sizeof(int8_t) * 21 * maxSeqLen;        // Size for profile data
    }

    static std::string getShmHash(const std::string& db);

    // Allocate and initialize shared memory
    static GPUSharedMemory* alloc(const std::string& name, unsigned int maxSeqLen, unsigned int maxResListLen);

    // Deallocate shared memory
    static void dealloc(GPUSharedMemory* layout, const std::string& name);

    static void unmap(GPUSharedMemory* layout);

    // Function to open and map existing shared memory and automatically determine sizes
    static GPUSharedMemory* openSharedMemory(const std::string& name);

};

struct GPUSharedMemorySem {
#ifdef USE_GPU_SEM
    sem_t* sem;
    std::string shmName;

    GPUSharedMemorySem() : sem(SEM_FAILED) {}
    void create(const std::string& name) {
        shmName = name;
        std::string semName = "/" + name + "_sem";
        sem_unlink(semName.c_str());
        sem = sem_open(semName.c_str(), O_CREAT, 0660, 0);
        if (sem == SEM_FAILED) {
            perror(("sem_open(create) " + semName).c_str());
            exit(EXIT_FAILURE);
        }
    }
    void open(const std::string& name) {
        shmName = name;
        std::string semName = "/" + name + "_sem";
        sem = sem_open(semName.c_str(), 0);
        if (sem == SEM_FAILED) {
            perror(("sem_open " + semName).c_str());
            exit(EXIT_FAILURE);
        }
    }
    void wait()  { while (sem_wait(sem) == -1 && errno == EINTR) {} }
    void post()  { if (sem != SEM_FAILED) sem_post(sem); }
    void close() { if (sem != SEM_FAILED) { sem_close(sem); sem = SEM_FAILED; } }
    void destroy() { close(); sem_unlink(("/" + shmName + "_sem").c_str()); }
#else
    GPUSharedMemorySem() {}
    void create(const std::string&) {}
    void open(const std::string&) {}
    void wait()    { sched_yield(); }
    void post()    {}
    void close()   {}
    void destroy() {}
#endif
};

#endif
