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

std::vector<Matcher::result_t> mapDomains(const std::vector<Matcher::result_t> &input, float overlap, float minCoverage,
                                          double eValThreshold) {
    std::vector<Matcher::result_t> result;
    if (input.empty()) {
        return result;
    }

    std::vector<bool> covered(input[0].qLen, false);
    for (size_t i = 0; i < input.size(); i++) {
        Matcher::result_t domain = input[i];
        if (domain.qStartPos > domain.qLen || domain.qEndPos > domain.qLen) {
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
        if (domain.dbStartPos > domain.dbLen || domain.dbEndPos > domain.dbLen) {
            Debug(Debug::WARNING) << "Target alignment start or end is greater than target length! Skipping line.\n";
            continue;
        }
        float targetCov = MathUtil::getCoverage(domain.dbStartPos, domain.dbEndPos, domain.dbLen);
        if (percentageOverlap <= overlap && targetCov > minCoverage && domain.eval < eValThreshold) {
            for (int j = domain.qStartPos; j < domain.qEndPos; ++j) {
                covered[j] = true;
            }
            result.push_back(domain);
        }
    }
    return result;
}

int doSummarize(Parameters &par, DBReader<unsigned int> &blastTabReader,
                const std::pair<std::string, std::string> &resultdb,
                const size_t dbFrom, const size_t dbSize) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBWriter writer(resultdb.first.c_str(), resultdb.second.c_str(), static_cast<unsigned int>(par.threads));
    writer.open();
    Debug(Debug::INFO) << "Start writing to file " << resultdb.first << "\n";
#pragma omp parallel for schedule(dynamic, 100)
    for (size_t i = dbFrom; i < dbFrom + dbSize; ++i) {
        unsigned int thread_idx = 0;
        char buffer[1024+32768];
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        unsigned int id = blastTabReader.getDbKey(i);
        char *tabData = blastTabReader.getData(i);
        const std::vector<Matcher::result_t> entries = Matcher::readAlignmentResults(tabData);
        if (entries.size() == 0) {
            Debug(Debug::WARNING) << "Could not map any entries for entry " << id << "!\n";
            continue;
        }

        std::vector<Matcher::result_t> result = mapDomains(entries, par.overlap, par.cov, par.evalThr);
        if (result.size() == 0) {
            Debug(Debug::WARNING) << "Could not map any domains for entry " << id << "!\n";
            continue;
        }
        std::string annotation;
        annotation.reserve(1024*1024);
        for (size_t j = 0; j < result.size(); j++) {
            size_t len = Matcher::resultToBuffer(buffer, result[j], par.addBacktrace);
            annotation.append(buffer, len);
        }
        writer.writeData(annotation.c_str(), annotation.length(), id, thread_idx);
    }
    writer.close();

    return EXIT_SUCCESS;
}

int doSummarize(Parameters &par, const unsigned int mpiRank, const unsigned int mpiNumProc) {
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(reader.getAminoAcidDBSize(), reader.getSeqLens(), reader.getSize(),
                                     mpiRank, mpiNumProc, &dbFrom, &dbSize);
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db2, par.db2Index, mpiRank);

    int status = doSummarize(par, reader, tmpOutput, dbFrom, dbSize);

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
        DBWriter::mergeResults(par.db2, par.db2 + ".index", splitFiles);
    }
    return status;
}

int doSummarize(Parameters &par) {
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t resultSize = reader.getSize();
    int status = doSummarize(par, reader, std::make_pair(par.db2, par.db2Index), 0, resultSize);
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
