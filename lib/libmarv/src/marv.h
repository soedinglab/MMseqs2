#ifndef MARV_H
#define MARV_H

#include <string>
#include <cstdint>
#include <vector>

class Marv {
public:
    enum AlignmentType {
        GAPLESS,
        SMITH_WATERMAN,
        GAPLESS_ENDPOS,
        SMITH_WATERMAN_ENDPOS,
        GAPLESS_THEN_SMITH_WATERMAN,
        FOLDSEEK_3DI_AA12
    };

    struct MemoryConfig {
        size_t maxGpuMem = SIZE_MAX; // per device: database cache and transfer buffers
        size_t maxBatchBytes = 128ull * 1024 * 1024;
        size_t maxBatchSequences = 10'000'000;
        size_t maxCachedQueryBytes = 512ull * 1024 * 1024; // score-only workspace kept between searches, 0 disables
    };
    Marv(size_t dbEntries, int alphabetSize);
    Marv(size_t dbEntries, int alphabetSize, const MemoryConfig& memoryConfig);
    void releaseSearchWorkspace();
    ~Marv();

    static std::vector<int> getDeviceIds();
    void* loadDb(char* data, size_t* offset, int32_t* length, size_t dbByteSize);
    void* loadDb(char* data, size_t dbByteSize, void* otherdb);
    void setDb(void* dbhandle);
    void setDbWithAllocation(void* dbhandle, const std::string& allocationinfo);
    std::string getDbMemoryHandle();

    void printInfo();
    void prefetch();

    /*
        Databases that are not in the padded GPU format may store letters in another encoding (e.g. ASCII).
        table has 256 entries mapping a stored byte to its alphabet code; the kernels translate while scanning.
        Only available on the Metal and SIMD backends.
    */
    void setLetterTranslation(const unsigned char* table);

    struct Stats {
        size_t results;
        double seconds;
        double gcups;
    };

    struct Result {
        unsigned int id;
        int score;
        int qEndPos;
        int dbEndPos;

        Result(unsigned int id, int score, int qEndPos, int dbEndPos) :
            id(id), score(score), qEndPos(qEndPos), dbEndPos(dbEndPos) {};
    };

    /*
        Score-only database scan for batch of queries.
        AlignmentType must be GAPLESS or SMITH_WATERMAN. 
        gop / gex are gap penalty for SMITH_WATERMAN. gap length i has score gop + (i-1) * gex.

        return topNSize highest scoring matches per query.
        resulting end positions will be -1,-1
    */
    std::vector<Marv::Stats> search_scoreOnly(
        Marv::Result* results,
        const std::vector<const char*>& queries,
        const std::vector<int>& queryLengths,
        const std::vector<int8_t*>& queryPssms,
        int topNSize, 
        Marv::AlignmentType alignmentType, 
        int gop, 
        int gex,
        bool int8IsAllowed, // allow alignments in 8-bit precision. automatic fallback to greater precision if 8-bit is not available
        bool handleOverflows, // recompute overflows; use 32-bit retries when 16-bit cannot guarantee exact scores
        std::vector<std::int32_t>* subsetIdsPtr_ = nullptr //. dont search full database, but only the db sequences with subsetids. nullptr if unused
    );

    Marv::Stats search_scoreOnly(
        Marv::Result* results,
        const char* query,
        int queryLength,
        int8_t* queryPssm,
        int topNSize, 
        Marv::AlignmentType alignmentType, 
        int gop, 
        int gex,
        bool int8IsAllowed,
        bool handleOverflows,
        std::vector<std::int32_t>* subsetIdsPtr_ = nullptr //nullptr if unused
    );

    /*
        Score + endposition database scan for batch of queries. 
        (
            Note: Internally, end positions will be computed for all alignments, not only for tops.
            Use search_scoreOnly_then_scoreEndpos_for_tops to restrict end position computation to top results
        )
        AlignmentType must be GAPLESS_ENDPOS or SMITH_WATERMAN_ENDPOS. 
        gop / gex are gap penalty for SMITH_WATERMAN. gap length i has score gop + (i-1) * gex.

        return topNSize highest scoring matches per query.
    */
    std::vector<Marv::Stats> search_scoreEndpos(
        Marv::Result* results,
        const std::vector<const char*>& queries,
        const std::vector<int>& queryLengths,
        const std::vector<int8_t*>& queryPssms,
        int topNSize, 
        Marv::AlignmentType alignmentType, 
        int gop, 
        int gex,
        bool int8IsAllowed, // allow alignments in 8-bit precision. automatic fallback to greater precision if 8-bit is not available
        bool handleOverflows, // recompute overflows (currently this means retry 8-bit alignments with 16-bit)
        std::vector<std::int32_t>* subsetIdsPtr_ = nullptr //. dont search full database, but only the db sequences with subsetids. nullptr if unused
    );

    Marv::Stats search_scoreEndpos(
        Marv::Result* results,
        const char* query,
        int queryLength,
        int8_t* queryPssm,
        int topNSize, 
        Marv::AlignmentType alignmentType, 
        int gop, 
        int gex,
        bool int8IsAllowed,
        bool handleOverflows,
        std::vector<std::int32_t>* subsetIdsPtr_ = nullptr //nullptr if unused
    );

    /*
        Score-only database scan for batch of queries, followed by end position computation for the top scoring alignments.
        The alignment score of the end-position pass will be the final score.
        scoreonly AlignmentType must be GAPLESS or SMITH_WATERMAN. 
        endpos AlignmentType must be GAPLESS_ENDPOS or SMITH_WATERMAN_ENDPOS. 
        gop / gex are gap penalty for SMITH_WATERMAN. gap length i has score gop + (i-1) * gex.

        return topNSize highest scoring matches per query.
    */
    std::vector<Marv::Stats> search_scoreOnly_then_scoreEndpos_for_tops(
        Marv::Result* results,
        const std::vector<const char*>& queries,
        const std::vector<int>& queryLengths,
        const std::vector<int8_t*>& queryPssms,
        int topNSize, 
        int gop, 
        int gex,
        Marv::AlignmentType alignmentType_scoreOnly, 
        bool int8IsAllowed_scoreOnly,
        bool handleOverflows_scoreOnly,
        Marv::AlignmentType alignmentType_endpos, 
        bool int8IsAllowed_endpos, 
        bool handleOverflows_endpos
    );

    Marv::Stats search_scoreOnly_then_scoreEndpos_for_tops(
        Marv::Result* results,
        const char* query,
        int queryLength,
        int8_t* queryPssm,
        int topNSize, 
        int gop, 
        int gex,
        Marv::AlignmentType alignmentType_scoreOnly, 
        bool int8IsAllowed_scoreOnly,
        bool handleOverflows_scoreOnly,
        Marv::AlignmentType alignmentType_endpos, 
        bool int8IsAllowed_endpos, 
        bool handleOverflows_endpos
    );

    /*
        Score-only gapless database scan for batch of queries using 3di, followed diagonal rescore of top results using aa12.

        return topNSize highest scoring matches per query.
    */
    std::vector<Marv::Stats> search_scoreOnly_foldseek(
        Marv::Result* results,
        const std::vector<int>& queryLengths,
        const std::vector<int8_t*>& queryPssms_3di,
        const std::vector<int8_t*>& queryPssms_aa12,
        int topNSize, 
        bool int8IsAllowed_scoreOnly,
        bool handleOverflows_scoreOnly,
        char maskingLetter3di,
        int maskingThreshold // in the database, 3di values in runs of length >= maskingThreshold will be replaced by maskingLetter3di
    );

    Marv::Stats search_scoreOnly_foldseek(
        Marv::Result* results,
        int queryLength,
        int8_t* queryPssm_3di,
        int8_t* queryPssm_aa12,
        int topNSize, 
        bool int8IsAllowed,
        bool handleOverflows,
        char maskingLetter3di,
        int maskingThreshold
    );



private:
    size_t dbEntries;
    int alphabetSize;

    // void* db;
    void* dbmanager;
    void* gpuDatabaseVoidPtr;
    void* searchWorkspaceVoidPtr;
};

#endif
