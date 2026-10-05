#ifndef MARV_H
#define MARV_H

#include <vector>

class Marv {
public:
    enum AlignmentType {
        GAPLESS,
        SMITH_WATERMAN,
        GAPLESS_SMITH_WATERMAN
    };

    Marv(size_t dbEntries, int alphabetSize, int maxSeqLength, size_t maxSeqs, AlignmentType alignmentType = AlignmentType::GAPLESS);
    ~Marv();

    static std::vector<int> getDeviceIds();
    void* loadDb(char* data, size_t* offset, int32_t* length, size_t dbByteSize);
    void* loadDb(char* data, size_t dbByteSize, void* otherdb);
    void setDb(void* dbhandle);
    void setDbWithAllocation(void* dbhandle, const std::string& allocationinfo);
    std::string getDbMemoryHandle();

    void printInfo();
    void prefetch();

    void startTimer();
    void stopTimer();

    struct Stats {
        size_t results;
        int numOverflows;
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

    //sequence must be encoded
    Stats scan(const char* sequence,  size_t sequenceLength, int8_t* pssm, Result* results);

    // Optional asynchronous 2-phase form of scan(), letting a caller overlap
    // building the next query's profile with this query's GPU execution.
    // submitScan() dispatches GPU work without blocking; it must be followed
    // by exactly one collectScan() (which blocks until that work completes)
    // before submitScan() may be called again. Only the Metal/MPS backend
    // implements these; the CUDA/HIP backend does not define them, so
    // callers must guard use behind #ifdef HAVE_MPS.
    void submitScan(const char* sequence, size_t sequenceLength, int8_t* pssm);
    Stats collectScan(Result* results);

private:
    size_t dbEntries;
    int alphabetSize;

    void* cudasw;
    // void* db;
    void* dbmanager;
    AlignmentType alignmentType;
};

#endif
