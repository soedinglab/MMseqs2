#include "Align2ClustCheckpoint.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"

#include <cstdio>
#include <cstring>
#include <sstream>

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

static std::string computeFingerprintHex(const Parameters &par, size_t dbSize, size_t endRange,
                                          int mode, size_t chunkSize) {
    // Every field below can change the *result* of clustering and therefore must be
    // part of the fingerprint. Rank count, thread count, and the runner used to
    // launch the job are intentionally excluded so the fingerprint (and thus the
    // checkpoint directory) stays identical across topology changes.
    std::ostringstream key;
    key << "db1=" << par.db1 << ";db1size=" << dbSize
        << ";db2=" << par.db2 << ";endRange=" << endRange
        << ";mode=" << mode << ";chunkSize=" << chunkSize
        << ";covMode=" << par.covMode << ";covThr=" << par.covThr
        << ";seqIdThr=" << par.seqIdThr << ";seqIdMode=" << par.seqIdMode
        << ";evalThr=" << par.evalThr << ";alnLenThr=" << par.alnLenThr
        << ";alignmentMode=" << par.alignmentMode
        << ";compBiasCorrection=" << par.compBiasCorrection
        << ";compBiasCorrectionScale=" << par.compBiasCorrectionScale
        << ";gapOpen=" << par.gapOpen.values.aminoacid()
        << ";gapExtend=" << par.gapExtend.values.aminoacid()
        << ";zdrop=" << par.zdrop
        << ";scoringMatrixFile=" << par.scoringMatrixFile.values.aminoacid()
        << ";addBacktrace=" << par.addBacktrace
        << ";includeAlignFiles=" << par.includeAlignFiles
        << ";filterCluDBFile=" << par.filterCluDBFile
        << ";filterSeqDBFile=" << par.filterSeqDBFile;

    const std::string keyStr = key.str();
    const size_t hash = Util::hash<const char>(keyStr.data(), keyStr.size());

    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%016zx", hash);
    return std::string(buffer);
}

static void ensureDirExists(const std::string &dir) {
    if (FileUtil::directoryExists(dir.c_str())) {
        return;
    }
    if (FileUtil::makeDir(dir.c_str()) == false && FileUtil::directoryExists(dir.c_str()) == false) {
        Debug(Debug::ERROR) << "Could not create directory " << dir << "\n";
        EXIT(EXIT_FAILURE);
    }
}

Align2ClustCheckpoint::Align2ClustCheckpoint(const Parameters &par, const std::string &baseOutputPath,
                                             size_t dbSize, size_t endRange, int mode, size_t chunkSize)
    : parentDir_(baseOutputPath + "_chunks"),
      fingerprint_(computeFingerprintHex(par, dbSize, endRange, mode, chunkSize)),
      chunkDir_(parentDir_ + "/" + fingerprint_) {
}

void Align2ClustCheckpoint::ensureChunkDirExists() const {
    ensureDirExists(parentDir_);
    ensureDirExists(chunkDir_);
}

std::string Align2ClustCheckpoint::dataPath(size_t chunkIndex) const {
    return chunkDir_ + "/chunk_" + SSTR(chunkIndex) + ".bin";
}

std::string Align2ClustCheckpoint::donePath(size_t chunkIndex) const {
    return chunkDir_ + "/chunk_" + SSTR(chunkIndex) + ".done";
}

bool Align2ClustCheckpoint::isChunkDone(size_t chunkIndex) const {
    return FileUtil::fileExists(donePath(chunkIndex).c_str());
}

void Align2ClustCheckpoint::writeChunk(size_t chunkIndex, const std::vector<ClusterResult> &results) const {
    const std::string finalPath = dataPath(chunkIndex);
    const std::string tmpPath = finalPath + ".tmp";

    FILE *fp = fopen(tmpPath.c_str(), "wb");
    if (fp == NULL) {
        Debug(Debug::ERROR) << "Could not open " << tmpPath << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    const uint64_t resultCount = static_cast<uint64_t>(results.size());
    fwrite(&resultCount, sizeof(uint64_t), 1, fp);
    for (const ClusterResult &result : results) {
        const uint64_t sequenceIdx = static_cast<uint64_t>(result.sequenceIdx);
        const uint64_t representativeId = static_cast<uint64_t>(result.representativeId);
        const uint64_t prefSize = static_cast<uint64_t>(result.prefSize);
        const uint64_t memberCount = static_cast<uint64_t>(result.memberIds.size());
        fwrite(&sequenceIdx, sizeof(uint64_t), 1, fp);
        fwrite(&representativeId, sizeof(uint64_t), 1, fp);
        fwrite(&prefSize, sizeof(uint64_t), 1, fp);
        fwrite(&memberCount, sizeof(uint64_t), 1, fp);
        for (DBLocalId memberId : result.memberIds) {
            const uint64_t member = static_cast<uint64_t>(memberId);
            fwrite(&member, sizeof(uint64_t), 1, fp);
        }
    }
    fflush(fp);
    fsync(fileno(fp));
    fclose(fp);

    if (rename(tmpPath.c_str(), finalPath.c_str()) != 0) {
        Debug(Debug::ERROR) << "Could not rename " << tmpPath << " to " << finalPath << "\n";
        EXIT(EXIT_FAILURE);
    }

    // The empty ".done" marker is created only after the data file has been fully
    // written and atomically renamed into place, so isChunkDone() never observes a
    // partially written chunk.
    FILE *doneFp = fopen(donePath(chunkIndex).c_str(), "wb");
    if (doneFp == NULL) {
        Debug(Debug::ERROR) << "Could not create " << donePath(chunkIndex) << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(doneFp);
}

std::vector<ClusterResult> Align2ClustCheckpoint::readChunk(size_t chunkIndex) const {
    std::vector<ClusterResult> results;

    const std::string path = dataPath(chunkIndex);
    FILE *fp = fopen(path.c_str(), "rb");
    if (fp == NULL) {
        Debug(Debug::ERROR) << "Could not open " << path << " for reading\n";
        EXIT(EXIT_FAILURE);
    }

    uint64_t resultCount = 0;
    if (fread(&resultCount, sizeof(uint64_t), 1, fp) != 1) {
        Debug(Debug::ERROR) << "Could not read chunk header from " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    results.reserve(resultCount);
    for (uint64_t i = 0; i < resultCount; ++i) {
        uint64_t sequenceIdx = 0, representativeId = 0, prefSize = 0, memberCount = 0;
        if (fread(&sequenceIdx, sizeof(uint64_t), 1, fp) != 1 ||
            fread(&representativeId, sizeof(uint64_t), 1, fp) != 1 ||
            fread(&prefSize, sizeof(uint64_t), 1, fp) != 1 ||
            fread(&memberCount, sizeof(uint64_t), 1, fp) != 1) {
            Debug(Debug::ERROR) << "Truncated chunk file " << path << "\n";
            EXIT(EXIT_FAILURE);
        }

        ClusterResult result;
        result.sequenceIdx = static_cast<size_t>(sequenceIdx);
        result.representativeId = static_cast<DBLocalId>(representativeId);
        result.prefSize = static_cast<size_t>(prefSize);
        result.memberIds.resize(memberCount);
        for (uint64_t m = 0; m < memberCount; ++m) {
            uint64_t member = 0;
            if (fread(&member, sizeof(uint64_t), 1, fp) != 1) {
                Debug(Debug::ERROR) << "Truncated chunk file " << path << "\n";
                EXIT(EXIT_FAILURE);
            }
            result.memberIds[m] = static_cast<DBLocalId>(member);
        }
        results.push_back(std::move(result));
    }
    fclose(fp);
    return results;
}

void Align2ClustCheckpoint::cleanup() const {
    DIR *dir = opendir(chunkDir_.c_str());
    if (dir != NULL) {
        struct dirent *entry;
        while ((entry = readdir(dir)) != NULL) {
            const std::string name = entry->d_name;
            if (name == "." || name == "..") {
                continue;
            }
            FileUtil::remove((chunkDir_ + "/" + name).c_str());
        }
        closedir(dir);
    }
    rmdir(chunkDir_.c_str());
    // Only removes parentDir_ if this was the last fingerprint subdirectory in it;
    // harmless no-op (ENOTEMPTY) otherwise.
    rmdir(parentDir_.c_str());
}
