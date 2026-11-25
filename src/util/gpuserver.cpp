#include "Parameters.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "MemoryMapped.h"
#include "IndexReader.h"
#include "SubstitutionMatrix.h"
#include "NucleotideMatrix.h"

// #define HAVE_CUDA 1
#ifdef HAVE_CUDA
#include "GpuUtil.h"
#endif

#include <random>
#include <fcntl.h>
#include <sys/mman.h>
#include <signal.h>

volatile sig_atomic_t keepRunning = 1;
void intHandler(int) {
    keepRunning = 0;
}

int gpuserver(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
#ifdef HAVE_METAL
    Debug(Debug::ERROR) << "gpuserver is not supported with Metal\n";
    EXIT(EXIT_FAILURE);
#elif defined(HAVE_CUDA)
    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader dbrIdx(par.db1, par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0 );
    DBReader<DBKeyType>* dbr = dbrIdx.sequenceReader;

    const bool isGpuDb = DBReader<DBKeyType>::getExtendedDbtype(dbr->getDbtype()) & Parameters::DBTYPE_EXTENDED_GPU;
    if (isGpuDb == false) {
        Debug(Debug::ERROR) << "Database " << FileUtil::baseName(par.db1) << " is not a valid GPU database\n"
                            << "Please call: makepaddedseqdb " << FileUtil::baseName(par.db1) << " " << FileUtil::baseName(par.db1) << "_pad\n";
        EXIT(EXIT_FAILURE);
    }

    const bool packedDb = DBReader<DBKeyType>::getExtendedDbtype(dbr->getDbtype()) & Parameters::DBTYPE_EXTENDED_AUX_SEQ;
    if (packedDb && par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED) {
        Debug(Debug::ERROR) << "--prefilter-mode " << par.prefMode << " is not supported for packed two-alphabet databases\n";
        EXIT(EXIT_FAILURE);
    }
    const unsigned int maxResults = packedDb
        ? (unsigned int) std::max(par.maxResListLen, std::min(dbr->getSize(), (size_t) Parameters::AUX_RESCORE_TOP_N))
        : (unsigned int) par.maxResListLen;

    std::vector<size_t> offsets;
    offsets.reserve(dbr->getSize() + 1);

    std::vector<int32_t> lengths;
    lengths.reserve(dbr->getSize());
    for(size_t id = 0; id < dbr->getSize(); id++){
        offsets.emplace_back(dbr->getIndex()[id].offset);
        lengths.emplace_back(dbr->getIndex()[id].length - 2);
    }
    offsets.emplace_back(offsets.back() + lengths.back());

    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(dbrIdx.sequenceReader->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    }

    const bool gapped = par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED;
    Marv marv(dbr->getSize(), subMat->alphabetSize);
    void* h1 = marv.loadDb(
            dbr->getDataForFile(0), offsets.data(), lengths.data(), dbr->getDataSizeForFile(0)
    );
    marv.setDb(h1);
    marv.prefetch();

    struct sigaction act;
    memset(&act, 0, sizeof(act));
    act.sa_handler = intHandler;

    // Set up the handler for SIGINT and SIGTERM
    sigaction(SIGINT, &act, NULL);
    sigaction(SIGTERM, &act, NULL);

    std::string shmFile = GPUSharedMemory::getShmHash(par.db1);
    GPUSharedMemory* layout = GPUSharedMemory::alloc(shmFile, par.maxSeqLen, maxResults);
    layout->packedDb = packedDb;
    Debug(Debug::WARNING) << shmFile << "\n";
    while (keepRunning) {
        if (layout->state.load(std::memory_order_acquire) == GPUSharedMemory::READY) {
            std::atomic_thread_fence(std::memory_order_acquire);

            Marv::Stats stats;
            if (packedDb) {
                stats = marv.search_scoreOnly_foldseek(
                    layout->getResultsPtr(), layout->queryLen, layout->getProfilePtr(), layout->getProfileAuxPtr(),
                    layout->topNSize, true, true, layout->maskLetter, layout->maskThreshold);
            } else {
                stats = marvSearch(marv, gapped, reinterpret_cast<const char *>(layout->getQueryPtr()), layout->queryLen,
                                   layout->getProfilePtr(), par.maxResListLen, layout->getResultsPtr());
            }
            layout->resultLen = stats.results;

            std::atomic_thread_fence(std::memory_order_release);
            // Debug(Debug::ERROR) << "switch to done\n";
            layout->state.store(GPUSharedMemory::DONE, std::memory_order_release);
        } else {
            std::this_thread::yield();
        }
    }

    // server shutdown
    // Debug(Debug::ERROR) << "switch to exit\n";
    layout->serverExit.store(true, std::memory_order_release);
    std::atomic_thread_fence(std::memory_order_release);
    GPUSharedMemory::dealloc(layout, shmFile);
#endif
    return EXIT_SUCCESS;

}
