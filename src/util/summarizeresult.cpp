#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Matcher.h"

#ifdef OPENMP
#include <omp.h>
#endif

int summarizeresult(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    MMseqsMPI::init(argc, argv);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

#ifdef HAVE_MPI
    size_t dbFrom = 0;
    size_t dbSize = 0;
    reader.decomposeDomainByAminoAcid(MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db2, par.db2Index, MMseqsMPI::rank);
    const char* outData = tmpOutput.first.c_str();
    const char* outIndex = tmpOutput.second.c_str();
    const bool merge = true;
#else
    size_t dbFrom = 0;
    size_t dbSize = reader.getSize();
    const char* outData = par.db2.c_str();
    const char* outIndex = par.db2Index.c_str();
    const bool merge = false;
#endif

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, dbSize), (size_t)1);
#endif

    DBWriter writer(outData, outIndex, localThreads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    Debug::Progress progress(dbSize);
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        char buffer[1024 + 32768*4];
        std::vector<bool> covered(par.maxSeqLen + 1, false);

#pragma omp for schedule(dynamic, 10)
        for (size_t i = dbFrom; i < dbFrom + dbSize; ++i) {
            progress.updateProgress();
            char *data = reader.getData(i, thread_idx);

            bool readFirst = false;
            writer.writeStart(thread_idx);
            while (*data != '\0') {
                Matcher::result_t domain = Matcher::parseAlignmentRecord(data, true);
                data = Util::skipLine(data);

                if (readFirst == false) {
                    covered.reserve(domain.qLen);
                    std::fill_n(covered.begin(), domain.qLen, false);
                    readFirst = true;
                }

                if (domain.qStartPos > static_cast<int>(domain.qLen) || domain.qEndPos > static_cast<int>(domain.qLen)) {
                    Debug(Debug::WARNING) << "Query alignment start or end is greater than query length! Skipping line.\n";
                    continue;
                }
                if (domain.dbcov < par.covThr) {
                    continue;
                }

                size_t counter = 0;
                for (int j = std::min(domain.qStartPos, domain.qEndPos); j < std::max(domain.qStartPos, domain.qEndPos); ++j) {
                    counter += covered[j] ? 1 : 0;
                }
                const float percentageOverlap = static_cast<float>(counter) / static_cast<float>(std::max(domain.qStartPos, domain.qEndPos) - std::min(domain.qStartPos, domain.qEndPos) + 1);
                if (percentageOverlap <= par.overlap) {
                    for (int j = std::min(domain.qStartPos, domain.qEndPos); j < std::max(domain.qStartPos, domain.qEndPos); ++j) {
                        covered[j] = true;
                    }
                    size_t len = Matcher::resultToBuffer(buffer, domain, par.addBacktrace, false);
                    writer.writeAdd(buffer, len, thread_idx);
                }
            }
            writer.writeEnd(reader.getDbKey(i), thread_idx);
        }
    }
    writer.close(merge);

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    if (MMseqsMPI::isMaster()) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (int i = 0; i < MMseqsMPI::numProc; ++i) {
            splitFiles.push_back(Util::createTmpFileNames(par.db2, par.db2Index, i));
        }
        DBWriter::mergeResults(par.db2, par.db2Index, splitFiles);
    }
#endif

    return EXIT_SUCCESS;
}

