#include "Lin8Db.h"
#include "Lin8DbReader.h"
#include "Parameters.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"
#include "Timer.h"

#include <algorithm>
#include <cstdio>
#include <vector>

static void readPipelineShape(const std::string &path, unsigned int &nodes, size_t &repRankBlocks, uint64_t &ranks) {
    FILE *in = fopen(path.c_str(), "r");
    if (in == NULL) {
        Debug(Debug::ERROR) << "Cannot open " << path << ". Run lin8-pref first\n";
        EXIT(EXIT_FAILURE);
    }
    nodes = 0;
    repRankBlocks = 0;
    size_t seen = 0;
    const bool read = fscanf(in, "nodes\t%u\nrepRankBlocks\t%zu\nranks\t%zu", &nodes, &repRankBlocks, &seen) == 3;
    fclose(in);
    ranks = seen;
    if (read == false || nodes == 0 || repRankBlocks == 0) {
        Debug(Debug::ERROR) << path << " does not name a set of aligned repRankBlocks\n";
        EXIT(EXIT_FAILURE);
    }
}

static void readNodeRepRankBlock(const std::string &prefix, unsigned int node, size_t repRankBlock,
                                 size_t budget, std::vector<PairRecord> &rows) {
    rows.clear();
    const std::string path = prefix + "." + SSTR(node) + "." + SSTR(repRankBlock);
    FILE *in = NULL;
    // the marker follows the rename, so a miss here is attribute cache lag and never permanent
    for (unsigned int waited = 0; (in = fopen(path.c_str(), "r")) == NULL; waited++) {
        if (waited >= 600) {
            Debug(Debug::ERROR) << "Cannot open " << path << ", which its marker promises\n";
            EXIT(EXIT_FAILURE);
        }
        sleep(1);
    }
    const size_t bytes = FileUtil::getFileSize(path);
    requireMemory("Representative rank block " + SSTR(repRankBlock),
                 bytes / PairRecord::DISK_BYTES * sizeof(PairRecord), budget, "raise --pair-splits");
    rows.reserve(bytes / PairRecord::DISK_BYTES);
    std::vector<PairRecord> buffer(1u << 16);
    size_t read = 0;
    while ((read = readRecords(buffer.data(), buffer.size(), in)) > 0) {
        rows.insert(rows.end(), buffer.begin(), buffer.begin() + read);
    }
    if (ferror(in) != 0) {
        Debug(Debug::ERROR) << "Cannot read " << path << "\n";
        EXIT(EXIT_FAILURE);
    }
    fclose(in);
}

int lin8align2clustmulti(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    FileUtil::fixRlimitNoFile();
    unsigned int alignNodes = 0;
    size_t repRankBlocks = 0;
    uint64_t ranks = 0;
    readPipelineShape(par.db2, alignNodes, repRankBlocks, ranks);
    if (par.lin8RepRankBlock < 0 || (size_t) par.lin8RepRankBlock >= repRankBlocks) {
        Debug(Debug::ERROR) << "--repRankBlock must name one of the " << repRankBlocks << " repRankBlocks\n";
        EXIT(EXIT_FAILURE);
    }
    const size_t firstRepRankBlock = (size_t) par.lin8RepRankBlock;
    const size_t lastRepRankBlock =
        par.lin8RepRankBlockCount <= 0
            ? repRankBlocks
            : std::min(firstRepRankBlock + (size_t) par.lin8RepRankBlockCount, repRankBlocks);
    ClusterAssignmentBitmap assignedCluster;
    assignedCluster.open(par.db3 + ".cluster_assigned", ranks);
    assignedCluster.catchUpTo(par.db3, firstRepRankBlock);

    const size_t budget = Util::computeMemory(par.splitMemoryLimit);
    Timer timer;
    uint64_t clusters = 0;
    uint64_t assigned = 0;
    std::vector<PairRecord> rows;
    std::vector<PairRecord> outBuffer;
    std::vector<uint64_t> members;
    const bool wantText = par.includeAlignFiles;
    std::string text;
    Debug::Progress progress(repRankBlocks);

    for (size_t repRankBlock = firstRepRankBlock; repRankBlock < lastRepRankBlock; repRankBlock++) {
        const std::string outPath = par.db3 + ".0." + SSTR(repRankBlock);
        if (FileUtil::fileExists(outPath.c_str())) {
            assignedCluster.catchUpTo(par.db3, repRankBlock + 1);
            progress.updateProgress(repRankBlock);
            continue;
        }
        const std::string outTmp = outPath + ".tmp" + uniqueTmpSuffix();
        FILE *out = FileUtil::openAndDelete(outTmp.c_str(), "w");
        const std::string textPath = alnTextPath(par.db3, 0, repRankBlock);
        const std::string textTmp = textPath + ".tmp" + uniqueTmpSuffix();
        FILE *textOut = NULL;
        if (wantText) {
            textOut = FileUtil::openAndDelete(textTmp.c_str(), "w");
            setvbuf(textOut, NULL, _IOFBF, 1u << 20);
        }
        uint64_t lastRep = 0;
        for (unsigned int chunk = 0; chunk < alignNodes; chunk++) {
        waitNodeDone(par.db1 + "." + SSTR(repRankBlock), chunk, 3600);
        AlnTextReader *textIn = wantText ? new AlnTextReader(par.db1, chunk, repRankBlock, true) : NULL;
        outBuffer.clear();
        readNodeRepRankBlock(par.db1, chunk, repRankBlock,
                         budget > assignedCluster.bytesHeld() ? budget - assignedCluster.bytesHeld() : 0, rows);
        size_t at = 0;
        while (at < rows.size()) {
            const uint64_t rep = rows[at].rep();
            if (rep < lastRep) {
                Debug(Debug::ERROR) << "Representative rank block " << repRankBlock << " goes back from representative "
                                    << lastRep << " to " << rep
                                    << ". The aligning pass did not cut it in rank order\n";
                EXIT(EXIT_FAILURE);
            }
            lastRep = rep;
            size_t end = at;
            while (end < rows.size() && rows[end].rep() == rep) {
                end++;
            }
            members.clear();
            for (size_t i = at; i < end; i++) {
                members.push_back(rows[i].member());
            }
            const size_t before = outBuffer.size();
            const bool made = assignCluster(rep, members.data(), members.size(), assignedCluster,
                                            outBuffer, assigned);
            clusters += made ? 1 : 0;
            // the text lines are one per pair in the same order, so the accepted ones fall out of a walk
            if (wantText) {
                text.clear();
                if (made) {
                    text.push_back(ALN_TEXT_REP_MARK);
                    text.push_back('\t');
                    text.append(SSTR(rep));
                    text.push_back('\n');
                }
                size_t j = before;
                for (size_t i = at; i < end; i++) {
                    char *line = NULL;
                    size_t length = 0;
                    if (textIn->next(line, length) == false) {
                        Debug(Debug::ERROR) << textIn->name() << " holds fewer lines than the "
                                            << "aligned pairs of repRankBlock " << repRankBlock << "\n";
                        EXIT(EXIT_FAILURE);
                    }
                    if (j < outBuffer.size() && rows[i].member() == outBuffer[j].member()) {
                        text.append(line, length);
                        j++;
                    }
                }
                if (text.empty() == false
                    && fwrite(text.c_str(), 1, text.size(), textOut) != text.size()) {
                    Debug(Debug::ERROR) << "Cannot write " << textTmp << "\n";
                    EXIT(EXIT_FAILURE);
                }
            }
            at = end;
            if (outBuffer.size() >= (1u << 16)) {
                if (writeRecords(outBuffer.data(), outBuffer.size(), out)
                    != outBuffer.size()) {
                    Debug(Debug::ERROR) << "Cannot write " << outTmp << "\n";
                    EXIT(EXIT_FAILURE);
                }
                outBuffer.clear();
            }
        }
        if (outBuffer.empty() == false
            && writeRecords(outBuffer.data(), outBuffer.size(), out) != outBuffer.size()) {
            Debug(Debug::ERROR) << "Cannot write " << outTmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (wantText) {
            char *line = NULL;
            size_t length = 0;
            if (textIn->next(line, length)) {
                Debug(Debug::ERROR) << textIn->name() << " holds more lines than the aligned pairs "
                                    << "of repRankBlock " << repRankBlock << "\n";
                EXIT(EXIT_FAILURE);
            }
            delete textIn;
        }
        }
        if (fclose(out) != 0) {
            Debug(Debug::ERROR) << "Cannot close " << outTmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        // published first, because the pair file is what a rerun takes as the block being decided
        if (wantText) {
            if (fclose(textOut) != 0) {
                Debug(Debug::ERROR) << "Cannot close " << textTmp << "\n";
                EXIT(EXIT_FAILURE);
            }
            FileUtil::publishAtomically(textTmp, textPath);
        }
        FileUtil::publishAtomically(outTmp, outPath);
        publishProgress(par.db3 + ".progress", repRankBlock + 1);
        if (par.removeTmpFiles) {
            removeConsumedBuckets(par.db2, alignNodes, repRankBlock, repRankBlock + 1, 1);
        }
        progress.updateProgress(repRankBlock);
    }
    assignedCluster.save(lastRepRankBlock);

    const std::string shapeTmp = par.db3 + "." + SSTR(firstRepRankBlock) + ".shape.tmp";
    FILE *shape = FileUtil::openAndDelete(shapeTmp.c_str(), "w");
    fprintf(shape, "repRankBlocks\t%zu\nranks\t%zu\n", repRankBlocks, (size_t) ranks);
    if (fclose(shape) != 0) {
        Debug(Debug::ERROR) << "Cannot close " << shapeTmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::publishAtomically(shapeTmp, par.db3);

    Debug(Debug::INFO) << "Made " << clusters << " clusters holding " << (clusters + assigned)
                       << " sequences in " << timer.lap() << "\n";
    return EXIT_SUCCESS;
}
