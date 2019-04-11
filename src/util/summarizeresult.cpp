#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "MathUtil.h"
#include "Matcher.h"

#ifdef OPENMP
#include <omp.h>
#endif

static inline float getOverlap(const std::vector<bool> &covered, unsigned int qStart, unsigned int qEnd) {
    size_t counter = 0;
    for (size_t i = qStart; i < qEnd; ++i) {
        counter += covered[i] ? 1 : 0;
    }
    return static_cast<float>(counter) / static_cast<float>(qEnd - qStart + 1);
}

int doSummarize(Parameters &par, DBReader<unsigned int> &resultReader,
                const std::pair<std::string, std::string> &resultdb,
                const size_t dbFrom, const size_t dbSize, bool merge) {
#ifdef OPENMP
    unsigned int totalThreads = par.threads;
#else
    unsigned int totalThreads = 1;
#endif

    unsigned int localThreads = totalThreads;
    if (resultReader.getSize() <= totalThreads) {
        localThreads = resultReader.getSize();
    }

    DBWriter writer(resultdb.first.c_str(), resultdb.second.c_str(), localThreads, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();
    Debug::Progress progress(dbSize);

#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[32768];

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);

        std::string annotation;
        annotation.reserve(1024*1024);

#pragma omp for schedule(dynamic, 100)
        for (size_t i = dbFrom; i < dbFrom + dbSize; ++i) {
            progress.updateProgress();

            unsigned int id = resultReader.getDbKey(i);
            char *tabData = resultReader.getData(i,  thread_idx);
            Matcher::readAlignmentResults(alnResults, tabData);
            if (alnResults.size() == 0) {
                Debug(Debug::WARNING) << "Can not map any alingment results for entry " << id << "!\n";
                continue;
            }

            std::vector<bool> covered(alnResults[0].qLen, false);
            for (size_t i = 0; i < alnResults.size(); i++) {
                Matcher::result_t domain = alnResults[i];
                if (domain.qStartPos > static_cast<int>(domain.qLen) || domain.qEndPos > static_cast<int>(domain.qLen)) {
                    Debug(Debug::WARNING) << "Query alignment start or end is greater than query length! Skipping line.\n";
                    continue;
                }
                if (domain.qStartPos > domain.qEndPos) {
                    Debug(Debug::WARNING) << "Query alignment end is greater than start! Skipping line.\n";
                    continue;
                }
                float percentageOverlap = getOverlap(covered, domain.qStartPos, domain.qEndPos);
                if (domain.dbStartPos > domain.dbEndPos) {
                    Debug(Debug::WARNING) << "Target alignment end is greater than start! Skipping line.\n";
                    continue;
                }
                if (domain.dbStartPos > static_cast<int>(domain.dbLen) || domain.dbEndPos > static_cast<int>(domain.dbLen)) {
                    Debug(Debug::WARNING) << "Target alignment start or end is greater than target length! Skipping line.\n";
                    continue;
                }
                float targetCov = MathUtil::getCoverage(domain.dbStartPos, domain.dbEndPos, domain.dbLen);
                if (percentageOverlap <= par.overlap && targetCov > par.covThr && domain.eval < par.evalThr) {
                    for (int j = domain.qStartPos; j < domain.qEndPos; ++j) {
                        covered[j] = true;
                    }
                    size_t len = Matcher::resultToBuffer(buffer, domain, par.addBacktrace);
                    annotation.append(buffer, len);
                }
            }

            writer.writeData(annotation.c_str(), annotation.length(), id, thread_idx);

            alnResults.clear();
            annotation.clear();
        }
    }
    writer.close(merge);

    return EXIT_SUCCESS;
}

int doSummarize(Parameters &par, const unsigned int mpiRank, const unsigned int mpiNumProc) {
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(reader.getDataSize(), reader.getSeqLens(), reader.getSize(),
                                     mpiRank, mpiNumProc, &dbFrom, &dbSize);
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db2, par.db2Index, mpiRank);

    int status = doSummarize(par, reader, tmpOutput, dbFrom, dbSize, true);

    reader.close();

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // master reduces results
    if (mpiRank == 0) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (unsigned int proc = 0; proc < mpiNumProc; ++proc) {
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(par.db2, par.db2Index, proc);
            splitFiles.push_back(std::make_pair(tmpFile.first, tmpFile.second));
        }
        DBWriter::mergeResults(par.db2, par.db2Index, splitFiles);
    }
    return status;
}

int doSummarize(Parameters &par) {
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t resultSize = reader.getSize();
    int status = doSummarize(par, reader, std::make_pair(par.db2, par.db2Index), 0, resultSize, false);
    reader.close();
    return status;
}

int summarizeresult(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    MMseqsMPI::init(argc, argv);

#ifdef HAVE_MPI
    int status = doSummarize(par, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int status = doSummarize(par);
#endif

    return status;
}
