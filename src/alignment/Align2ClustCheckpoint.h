#ifndef MMSEQS_ALIGN2CLUST_CHECKPOINT_H
#define MMSEQS_ALIGN2CLUST_CHECKPOINT_H

#include "Align2ClustReducer.h"

#include <cstddef>
#include <string>
#include <vector>

class Parameters;

// Manages on-disk chunk checkpoint state for a (potentially MPI-distributed)
// align2clust run.
//
// The checkpoint directory is named after a fingerprint computed from every
// parameter that affects the *result* of clustering (thresholds, alignment
// mode, input database identity/size, chunk size, ...). It deliberately
// excludes anything about *how* the work is distributed (rank count, thread
// count, MPI runner): those never change the fingerprint, so the exact same
// checkpoint directory -- and therefore the exact same per-chunk completion
// state -- is reused whether a run is resumed on the same topology, a
// different number of nodes, or a single node. A run interrupted on N nodes
// can always be resumed on M nodes (including M == 1) as long as the
// semantic clustering parameters are unchanged.
//
// Known limitation: input database identity is fingerprinted by path and
// sequence count, not file content checksum. Replacing db1/db2 in place with
// different content of the same size, between an interrupted run and its
// resume, will not be detected. Documented as an out-of-scope hardening item.
class Align2ClustCheckpoint {
public:
    // baseOutputPath is typically par.db3 (the align2clust cluster-result DB
    // path); the checkpoint directory is created as a sibling of it.
    Align2ClustCheckpoint(const Parameters &par, const std::string &baseOutputPath,
                          size_t dbSize, size_t endRange, int mode, size_t chunkSize);

    const std::string &getChunkDir() const { return chunkDir_; }
    const std::string &getFingerprint() const { return fingerprint_; }

    // Ensures the checkpoint directory (and its parent) exist on disk.
    void ensureChunkDirExists() const;

    bool isChunkDone(size_t chunkIndex) const;

    // Serializes results to disk and marks the chunk done. The data file is
    // written to a temporary name and rename()'d into place (atomic on a
    // POSIX filesystem) before the ".done" marker is created, so a chunk is
    // only ever observed as done once its data is fully and correctly on
    // disk; a crash mid-write leaves no ".done" marker and the chunk is
    // regenerated on the next run.
    void writeChunk(size_t chunkIndex, const std::vector<ClusterResult> &results) const;

    std::vector<ClusterResult> readChunk(size_t chunkIndex) const;

    // Removes the entire checkpoint directory (and its parent, if left
    // empty). Call only after every chunk has been consumed and the final
    // clustering result has been written successfully -- an interrupted run
    // must keep its checkpoint directory so it can be resumed.
    void cleanup() const;

private:
    std::string dataPath(size_t chunkIndex) const;
    std::string donePath(size_t chunkIndex) const;

    std::string parentDir_;
    std::string fingerprint_;
    std::string chunkDir_;
};

#endif
