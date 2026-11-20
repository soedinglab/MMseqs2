#ifndef GPUUTIL_H
#define GPUUTIL_H

#include <atomic>
#include <string>
#include <cstddef>
#include "marv.h"

// Ungapped GPU scan; gapped also rescores the top hits with Smith-Waterman end positions.
// gop/gex -11/-1 are the penalties the former Marv GAPLESS_SMITH_WATERMAN mode used.
inline Marv::Stats marvSearch(Marv& marv, bool gapped, const char* query, int queryLength,
                              int8_t* profile, int topN, Marv::Result* results) {
    if (gapped) {
        return marv.search_scoreOnly_then_scoreEndpos_for_tops(results, query, queryLength, profile, topN, -11, -1,
                                                               Marv::GAPLESS, true, true,
                                                               Marv::SMITH_WATERMAN_ENDPOS, true, true);
    }
    return marv.search_scoreOnly(results, query, queryLength, profile, topN, Marv::GAPLESS, 0, 0, true, true);
}

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

#endif
