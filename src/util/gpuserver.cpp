#include "Parameters.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "MemoryMapped.h"
#include "IndexReader.h"
#include "SubstitutionMatrix.h"
#include "NucleotideMatrix.h"

#ifdef HAVE_CUDA
#include "GpuUtil.h"
#include "marv.h"
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
#ifdef HAVE_CUDA
    bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader dbrIdx(par.db1, par.threads, IndexReader::SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0 );
    DBReader<unsigned int>* dbr = dbrIdx.sequenceReader;

    std::vector<size_t> offsets;
    offsets.reserve(dbr->getSize() + 1);

    std::vector<int32_t> lengths;
    lengths.reserve(dbr->getSize());
    for(size_t id = 0; id < dbr->getSize(); id++){
        offsets.emplace_back(dbr->getIndex()[id].offset);
        lengths.emplace_back(dbr->getIndex()[id].length - 2);
    }
    offsets.emplace_back(offsets.back() + lengths.back());
    int32_t maxTargetLength = lengths.back();

    std::string dbrName  = par.db1;
    std::string dbrRelPath = FileUtil::getRealPathFromSymLink(PrefilteringIndexReader::dbPathWithoutIndex(dbrName));
    size_t hash = Util::hash(dbrRelPath.c_str(), dbrRelPath.length());

    std::string shmFile = SSTR(hash);
    GPUSharedMemory* layout = GPUSharedMemory::alloc(shmFile, par.maxSeqLen , par.maxResListLen); // Adjust sizes as necessary

    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(dbrIdx.sequenceReader->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    }

    Marv::AlignmentType type =  (par.prefMode == Parameters::PREF_MODE_UNGAPPED_AND_GAPPED) ?
                                Marv::AlignmentType::GAPLESS_SMITH_WATERMAN : Marv::AlignmentType::GAPLESS;
    Marv marv(dbr->getSize(), subMat->alphabetSize, maxTargetLength, par.maxResListLen, type);
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
    Debug(Debug::WARNING) << shmFile << "\n";
    while (keepRunning) {
        while (layout->serverReady.load(std::memory_order_acquire) == 0 || layout->clientReady.load(std::memory_order_acquire) == 0) {
            std::this_thread::yield();
            if (keepRunning == false) {
                break;
            }
        }
        if (keepRunning == false) {
            break;
        }
        Marv::Stats stats = marv.scan(reinterpret_cast<const char *>(layout->getQueryPtr()), layout->queryLen, layout->getProfilePtr(), layout->getResultsPtr());
        layout->resultLen = stats.results;
        layout->serverReady.store(UINT_MAX, std::memory_order_release);
        while (layout->clientReady.load(std::memory_order_acquire) != 0) {
            std::this_thread::yield();  // Wait for client to finish
        }
        layout->resetServerAndClientReady();
    }

    GPUSharedMemory::dealloc(layout, shmFile);
#endif
    return EXIT_SUCCESS;

}
