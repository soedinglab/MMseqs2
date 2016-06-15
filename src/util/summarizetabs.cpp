#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "CompressedA3M.h"
#include "Alignment.h"
#include "MathUtil.h"
#include "Domain.h"

#include <fstream>
#include <iomanip>

#ifdef OPENMP
#include <omp.h>
#endif

std::vector<Domain> mapDomains(const std::vector<Domain> &input, float overlap, float minCoverage,
                               double eValThreshold) {
    std::vector<Domain> result;
    if(input.size() == 0) {
        return result;
    }

    std::vector<bool> covered(input[0].qLength, false);
    for (std::vector<Domain>::const_iterator it = input.begin(); it != input.end(); ++it) {
        float percentageOverlap = MathUtil::getOverlap(covered, (*it).qStart, (*it).qEnd);
        float targetCov = MathUtil::getCoverage((*it).tStart, (*it).tEnd, (*it).tLength);
        if (percentageOverlap <= overlap && targetCov > minCoverage && (*it).eValue < eValThreshold) {
            // FIXME possible bug with -1
            for (unsigned int i = ((*it).qStart - 1); i < ((*it).qEnd - 1); ++i) {
                covered[i] = true;
            }
            result.push_back(*it);
        }
    }

    return result;
}

std::map<std::string, unsigned int> readLength(const std::string &file) {
    std::fstream mappingStream(file);
    if (mappingStream.fail()) {
        Debug(Debug::ERROR) << "File " << file << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::map<std::string, unsigned int> mapping;
    std::string line;
    while (std::getline(mappingStream, line)) {
        std::vector<std::string> split = Util::split(line, "\t");
        unsigned int length = static_cast<unsigned int>(strtoul(split[1].c_str(), NULL, 10));
        mapping.emplace(split[0], length);
    }

    return mapping;
}

std::vector<Domain> getEntries(char *data, size_t length, const std::map<std::string, unsigned int> &lengths) {
    std::vector<Domain> result;

    std::string line;
    std::istringstream iss(std::string(data, length));
    while (std::getline(iss, line)) {
        size_t offset = 1;
        std::vector<std::string> fields = Util::split(line.c_str(), "\t");

        unsigned int qStart = static_cast<unsigned int>(strtoul(fields[offset + 6].c_str(), NULL, 10)) - 1;
        unsigned int qEnd = static_cast<unsigned int>(strtoul(fields[offset + 7].c_str(), NULL, 10)) - 1;

        std::map<std::string, unsigned int>::const_iterator it = lengths.lower_bound(fields[0]);
        unsigned int qLength;
        if(it != lengths.end()) {
            qLength = (*it).second;
        } else {
            Debug(Debug::WARNING) << "Missing query length! Skipping line.\n";
            continue;
        }


        unsigned int tStart = static_cast<unsigned int>(strtoul(fields[offset + 8].c_str(), NULL, 10)) - 1;
        unsigned int tEnd = static_cast<unsigned int>(strtoul(fields[offset + 9].c_str(), NULL, 10)) - 1;

        it = lengths.lower_bound(fields[offset + 1]);
        unsigned int tLength;
        if(it != lengths.end()) {
            tLength = (*it).second;
        } else {
            Debug(Debug::WARNING) << "Missing target length! Skipping line.\n";
            continue;
        }

        double eValue = strtod(fields[offset + 10].c_str(), NULL);

        result.emplace_back(fields[0], qStart, qEnd, qLength, fields[offset + 1], tStart, tEnd, tLength, eValue);
    }

    std::stable_sort(result.begin(), result.end());

    return result;
}

int doAnnotate(Parameters &par, DBReader<unsigned int> &blastTabReader,
               const std::pair<std::string, std::string>& resultdb,
               const size_t dbFrom, const size_t dbSize) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBWriter writer(resultdb.first.c_str(), resultdb.second.c_str(), static_cast<unsigned int>(par.threads));
    writer.open();

    std::map<std::string, unsigned int> lengths = readLength(par.db3);

    Debug(Debug::INFO) << "Start writing to file " << par.db4 << "\n";

#pragma omp parallel for schedule(dynamic, 100)
    for (size_t i = dbFrom; i < dbSize; ++i) {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        unsigned int id = blastTabReader.getDbKey(i);

        char* tabData = blastTabReader.getData(i);
        size_t tabLength = blastTabReader.getSeqLens(i) - 1;
        const std::vector<Domain> entries = getEntries(tabData, tabLength, lengths);
        if (entries.size() == 0) {
            Debug(Debug::WARNING) << "Could not map any entries for entry " << id << "!\n";
            continue;
        }

        std::vector<Domain> result = mapDomains(entries, par.overlap, par.cov, par.evalThr);
        if (result.size() == 0) {
            Debug(Debug::WARNING) << "Could not map any domains for entry " << id << "!\n";
            continue;
        }

        std::ostringstream oss;
        oss << std::setprecision(std::numeric_limits<float>::digits10);

        for (std::vector<Domain>::const_iterator j = result.begin(); j != result.end(); ++j) {
            (*j).writeResult(oss);
            oss << "\n";
        }

        std::string annotation = oss.str();
        writer.write(annotation.c_str(), annotation.length(), SSTR(id).c_str(), thread_idx);
    }

    writer.close();

    return EXIT_SUCCESS;
}

int doAnnotate(Parameters &par, const unsigned int mpiRank, const unsigned int mpiNumProc) {
    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::NOSORT);

    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(reader.getAminoAcidDBSize(), reader.getSeqLens(), reader.getSize(),
                                     mpiRank, mpiNumProc, &dbFrom, &dbSize);
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db3, par.db3Index, mpiRank);

    int status = doAnnotate(par, reader, tmpOutput, dbFrom, dbSize);

    reader.close();

#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // master reduces results
    if(mpiRank == 0) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for(unsigned int proc = 0; proc < mpiNumProc; ++proc){
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(par.db3, par.db3Index, proc);
            splitFiles.push_back(std::make_pair(tmpFile.first,  tmpFile.first + ".index"));
        }
        Alignment::mergeAndRemoveTmpDatabases(par.db3, par.db3 + ".index", splitFiles);
    }

    return status;
}

int doAnnotate(Parameters &par) {
    size_t resultSize;

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::NOSORT);
    resultSize = reader.getSize();

    int status = doAnnotate(par, reader, std::make_pair(par.db3, par.db3Index), 0, resultSize);

    reader.close();

    return status;
}

int summarizetabs(int argc, const char **argv) {
    MMseqsMPI::init(argc, argv);

    std::string usage("Extract annotations from HHblits blasttab results.\n");
    usage.append("Written by Milot Mirdita (milot@mirdita.de) & Martin Steinegger <martin.steinegger@mpibpc.mpg.de>\n");
    usage.append("USAGE: <blastTabDB> <lengthFile> <outDB>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.summarizetabs, 3);

#ifdef HAVE_MPI
    int status = doAnnotate(par, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int status = doAnnotate(par);
#endif

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return status;
}
