#include "Lin8DbReader.h"
#include "Parameters.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include "Timer.h"

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

static const size_t SHARD_BUFFER = 1u << 20;
static const uint64_t REPS_PER_PROGRESS = 1u << 20;

// the cluster index is text, one "key\toffset\tlength" line per representative, keys ascending
static bool nextRepKey(FILE *in, const std::string &path, uint64_t &key) {
    char line[256];
    if (fgets(line, sizeof(line), in) == NULL) {
        if (ferror(in) != 0) {
            Debug(Debug::ERROR) << "Cannot read " << path << "\n";
            EXIT(EXIT_FAILURE);
        }
        return false;
    }
    char *end = NULL;
    const uint64_t read = strtoull(line, &end, 10);
    if (end == line || *end != '\t' || strchr(line, '\n') == NULL) {
        Debug(Debug::ERROR) << path << " holds \"" << line
                            << "\", which is not a key, a tab and a line of its own\n";
        EXIT(EXIT_FAILURE);
    }
    key = read;
    return true;
}

static void writeAll(FILE *out, const std::string &path, const char *from, size_t length) {
    if (fwrite(from, 1, length, out) != length) {
        Debug(Debug::ERROR) << "Cannot write " << length << " byte to " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
}

int lin8createrepseqfasta(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    RunDbReader reader(par.db1, true);
    reader.open();

    const size_t splits = par.fastaSplits > 1 ? (size_t) par.fastaSplits : 1;
    std::vector<std::string> path(splits);
    std::vector<std::string> tmp(splits);
    std::vector<FILE *> out(splits, NULL);
    for (size_t split = 0; split < splits; split++) {
        if (splits == 1) {
            path[split] = par.db3;
        } else {
            char suffix[32];
            snprintf(suffix, sizeof(suffix), ".%05zu", split);
            path[split] = par.db3 + suffix;
        }
        tmp[split] = path[split] + ".tmp";
        out[split] = FileUtil::openAndDelete(tmp[split].c_str(), "w");
        setvbuf(out[split], NULL, _IOFBF, SHARD_BUFFER);
    }

    Timer timer;
    FILE *index = FileUtil::openFileOrDie(par.db2Index.c_str(), "r", true);
    RunDbReader::HeaderStream headers(reader);
    RunDbReader::Cursor cursor;
    Debug::Progress progress;
    const char *begin = NULL;
    size_t length = 0;
    uint64_t rank = 0;
    uint64_t reps = 0;
    uint64_t rep = 0;
    bool haveRep = nextRepKey(index, par.db2Index, rep);
    while (haveRep) {
        if (rep >= reader.getSize()) {
            Debug(Debug::ERROR) << "The clustering names rank " << rep << ", past the database\n";
            EXIT(EXIT_FAILURE);
        }
        while (rank <= rep) {
            if (headers.next(begin, length) == false) {
                Debug(Debug::ERROR) << "The headers hold " << rank << " entries and the clustering "
                                    << "names rank " << rep << "\n";
                EXIT(EXIT_FAILURE);
            }
            rank++;
        }
        const std::string name = Util::parseFastaHeader(begin);
        const uint32_t seqLen = reader.getSeqLen(rep, cursor);
        const char *residues = reader.getData(rep, cursor);
        FILE *shard = out[reps % splits];
        const std::string &shardPath = tmp[reps % splits];
        writeAll(shard, shardPath, ">", 1);
        writeAll(shard, shardPath, name.c_str(), name.size());
        writeAll(shard, shardPath, "\n", 1);
        writeAll(shard, shardPath, residues, seqLen);
        writeAll(shard, shardPath, "\n", 1);
        reps++;
        if (reps % REPS_PER_PROGRESS == 0) {
            progress.updateProgress();
        }
        uint64_t next = 0;
        haveRep = nextRepKey(index, par.db2Index, next);
        if (haveRep && next <= rep) {
            Debug(Debug::ERROR) << par.db2Index << " goes back from key " << rep << " to " << next
                                << ", so the representatives are not in rank order\n";
            EXIT(EXIT_FAILURE);
        }
        rep = next;
    }
    if (fclose(index) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << par.db2Index << "\n";
        EXIT(EXIT_FAILURE);
    }
    for (size_t split = 0; split < splits; split++) {
        if (fclose(out[split]) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << tmp[split] << "\n";
            EXIT(EXIT_FAILURE);
        }
        FileUtil::publishAtomically(tmp[split], path[split]);
    }

    reader.close();
    Debug(Debug::INFO) << "Wrote " << reps << " representatives in " << timer.lap() << "\n";
    return EXIT_SUCCESS;
}
