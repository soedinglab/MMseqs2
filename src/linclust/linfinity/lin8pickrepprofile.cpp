#include "Lin8DbReader.h"
#include "Parameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Util.h"
#include "Timer.h"
#include "Matcher.h"
#include "MultipleAlignment.h"
#include "PSSMCalculator.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "itoa.h"

#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

static const size_t ARENA_BYTES = 64u << 20;

static void appendDecimalKey(std::string &into, uint64_t key) {
    char buffer[32];
    char *end = Itoa::u64toa_sse2(key, buffer);
    into.append(buffer, end - buffer - 1);
    into.push_back('\n');
}

int lin8pickrepprofile(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    size_t blocks = 0;
    while (FileUtil::fileExists(alnTextPath(par.db2, 0, blocks).c_str())) {
        blocks++;
    }
    if (blocks == 0) {
        Debug(Debug::ERROR) << "No alignments beside " << par.db2
                            << ". Run the clustering with --include-align-files 1\n";
        EXIT(EXIT_FAILURE);
    }

    Lin8DbReader reader(par.db1);
    reader.open();
    const unsigned int threads = std::max<unsigned int>(1, par.threads);
    reader.openBatch(threads, ARENA_BYTES, Util::computeMemory(par.splitMemoryLimit),
                     Lin8DbReader::READ_ONCE);

    DBReader<DBKeyType> clusters(par.db3.c_str(), par.db3Index.c_str(), threads,
                                 DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    clusters.open(DBReader<DBKeyType>::NOSORT);
    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), threads, par.compressed,
                    Parameters::DBTYPE_CLUSTER_RES);
    writer.open();

    // the -0.2 match state adjustment result2profile trims its alignments with
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, -0.2f);
    const size_t maxSetSize = clusters.maxCount('\n') + 1;
    const size_t maxSeqLen = std::max<size_t>(reader.getIndex().getMaxSeqLen(), 1);
    Debug(Debug::INFO) << "Scoring the members of " << clusters.getSize()
                       << " clusters by profile PSSM, minimum coverage: " << par.covThr << "\n";

    Timer timer;
    uint64_t rewritten = 0;
    uint64_t switched = 0;
    Debug::Progress progress(blocks);
#pragma omp parallel for schedule(static, 1) num_threads(threads) reduction(+ : rewritten, switched)
    for (unsigned int part = 0; part < threads; part++) {
        MultipleAlignment aligner(maxSeqLen, &subMat);
        PSSMCalculator calculator(&subMat, maxSeqLen, maxSetSize, par.pcmode, par.pca, par.pcb
#ifdef GAP_POS_SCORING
                                  , par.gapOpen.values.aminoacid()
                                  , par.gapPseudoCount
#endif
        );
        Sequence centerSequence(maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false,
                                par.compBiasCorrection);
        Sequence edgeSequence(maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, false);
        std::vector<uint64_t> members;
        std::vector<Matcher::result_t> results;
        std::vector<size_t> order;
        std::vector<uint64_t> sorted;
        std::vector<Matcher::result_t> alnResults;
        std::vector<std::vector<unsigned char> > seqSet;
        std::string entry;

        for (size_t block = blocks * part / threads; block < blocks * (part + 1) / threads; block++) {
            AlnTextReader in(par.db2, 0, block, true);
            char *begin = NULL;
            size_t length = 0;
            bool have = in.next(begin, length);
            while (have) {
                if (begin[0] != ALN_TEXT_REP_MARK) {
                    Debug(Debug::ERROR) << in.name() << " holds a line outside a representative's "
                                        << "group\n";
                    EXIT(EXIT_FAILURE);
                }
                const uint64_t rep = strtoull(begin + 2, NULL, 10);
                members.clear();
                results.clear();
                while ((have = in.next(begin, length)) && begin[0] != ALN_TEXT_REP_MARK) {
                    const uint64_t member = strtoull(begin, NULL, 10);
                    if (member == rep) {
                        continue;
                    }
                    members.push_back(member);
                    results.push_back(Matcher::parseAlignmentRecord(begin));
                }

                // the batch reader wants ascending ranks, the rows have to stay in file order
                order.resize(members.size());
                for (size_t i = 0; i < members.size(); i++) {
                    order[i] = i;
                }
                std::sort(order.begin(), order.end(),
                          [&members](size_t a, size_t b) { return members[a] < members[b]; });
                sorted.clear();
                alnResults.clear();
                for (size_t i = 0; i < order.size(); i++) {
                    sorted.push_back(members[order[i]]);
                    alnResults.push_back(results[order[i]]);
                }

                const uint32_t repLen = reader.getSeqLen(rep);
                seqSet.clear();
                size_t from = 0;
                size_t got = sorted.empty()
                                 ? 0
                                 : reader.startBatch(rep, sorted.data(), sorted.size(), part, 0);
                Lin8DbReader::Cursor cursor;
                while (from < sorted.size()) {
                    reader.awaitBatch(part, 0);
                    if (from == 0) {
                        centerSequence.mapSequence(0, 0, reader.batchQueryAt(part, 0), repLen);
                    }
                    for (size_t k = 0; k < got; k++) {
                        const uint32_t len = reader.getSeqLen(sorted[from + k], cursor);
                        const Matcher::result_t &row = alnResults[from + k];
                        if (row.dbEndPos >= (int) len || row.qEndPos >= (int) repLen) {
                            Debug(Debug::ERROR) << in.name() << " aligns rank " << sorted[from + k]
                                                << " past the end of it or of rank " << rep << "\n";
                            EXIT(EXIT_FAILURE);
                        }
                        edgeSequence.mapSequence(0, 0, reader.batchAt(part, 0, k), len);
                        seqSet.push_back(std::vector<unsigned char>(
                            edgeSequence.numSequence, edgeSequence.numSequence + edgeSequence.L));
                    }
                    from += got;
                    if (from < sorted.size()) {
                        got = reader.startBatch(rep, sorted.data() + from, sorted.size() - from,
                                                part, 0);
                    }
                }
                if (sorted.empty()) {
                    centerSequence.mapSequence(0, 0, reader.getData(rep), repLen);
                }

                uint64_t newRep = rep;
                if (seqSet.empty() == false) {
                    MultipleAlignment::MSAResult res =
                        aligner.computeMSA(&centerSequence, seqSet, alnResults, true);
                    PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(
                        res.setSize, res.centerLength, (const char **) res.msaSequence,
#ifdef GAP_POS_SCORING
                        alnResults,
#endif
                        par.wg, 0.0);
                    double bestScore = -DBL_MAX;
                    bool anyPassed = false;
                    for (size_t row = 0; row < res.setSize; row++) {
                        const char *msaRow = res.msaSequence[row];
                        const uint64_t rank = (row == 0) ? rep : sorted[row - 1];
                        double score = 0.0;
                        size_t aligned = 0;
                        for (size_t pos = 0; pos < res.centerLength; pos++) {
                            const unsigned char state = (unsigned char) msaRow[pos];
                            if (state < Sequence::PROFILE_AA_SIZE) {
                                aligned++;
                                score += (int) pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + state];
                            }
                        }
                        const float coverage = res.centerLength > 0
                                                   ? (float) aligned / (float) res.centerLength
                                                   : 0.0f;
                        if (coverage < par.covThr) {
                            continue;
                        }
                        // the highest scoring row wins, and the lowest rank breaks a tie
                        if (anyPassed == false || score > bestScore
                            || (score == bestScore && rank < newRep)) {
                            bestScore = score;
                            newRep = rank;
                            anyPassed = true;
                        }
                    }
                    MultipleAlignment::deleteMSA(&res);
                }

                const size_t id = clusters.getId((DBKeyType) rep);
                if (id == DB_ENTRY_NOT_FOUND) {
                    Debug(Debug::ERROR) << "Representative " << rep << " of " << in.name()
                                        << " is not a cluster of " << par.db3 << "\n";
                    EXIT(EXIT_FAILURE);
                }
                entry.clear();
                appendDecimalKey(entry, newRep);
                char *data = clusters.getData(id, part);
                while (data != NULL && *data != '\0') {
                    const uint64_t member = strtoull(data, NULL, 10);
                    if (member != newRep) {
                        appendDecimalKey(entry, member);
                    }
                    data = Util::skipLine(data);
                }
                writer.writeData(entry.c_str(), entry.length(), (DBKeyType) newRep, part);
                rewritten++;
                switched += newRep != rep ? 1 : 0;
            }
            progress.updateProgress();
        }
    }
    // the keys are the new representatives, so they no longer ascend
    writer.close(false, true);
    if (rewritten != clusters.getSize()) {
        Debug(Debug::ERROR) << "The alignments name " << rewritten << " clusters and " << par.db3
                            << " holds " << clusters.getSize() << "\n";
        EXIT(EXIT_FAILURE);
    }
    clusters.close();
    reader.close();
    Debug(Debug::INFO) << "Switched the representative of " << switched << " of " << rewritten
                       << " clusters in " << timer.lap() << "\n";
    return EXIT_SUCCESS;
}
