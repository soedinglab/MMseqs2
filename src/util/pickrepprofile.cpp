#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "SubstitutionMatrix.h"
#include "Matcher.h"
#include "MultipleAlignment.h"
#include "PSSMCalculator.h"
#include "EvalueComputation.h"
#include "Sequence.h"

#include <climits>
#include <cfloat>
#include <cmath>

#ifdef OPENMP
#include <omp.h>
#endif

// Pick a new representative by scoring each observed member against the cluster profile
// (PSSM log-odds summed over its aligned match columns, the score a profile search uses)
// and taking the highest. MSA and profile are built like result2profile and are gap-free
// on the center, so each member is scored in one pass over its MSA row. Members below the
// coverage threshold (-c) are excluded so fragments cannot win on a conserved core alone.

int pickrepprofile(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<DBKeyType> seqReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
    seqReader.open(DBReader<DBKeyType>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        seqReader.readMmapedDataInMemory();
    }

    DBReader<DBKeyType> resultReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<DBKeyType>::USE_DATA | DBReader<DBKeyType>::USE_INDEX);
    resultReader.open(DBReader<DBKeyType>::LINEAR_ACCCESS);

    DBWriter resultWriter(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    resultWriter.open();

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, resultReader.getSize()), (size_t)1);
#endif

    // + 1 for the center sequence (see result2profile)
    const size_t maxSetSize = resultReader.maxCount('\n') + 1;
    const unsigned int maxSequenceLength = seqReader.getMaxSeqLen();

    // adjust score of each match state by -0.2 to trim alignment (identical to result2profile)
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, -0.2f);
    EvalueComputation evalueComputation(seqReader.getAminoAcidDBSize(), &subMat, par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid());

    if (seqReader.getDbtype() == -1) {
        Debug(Debug::ERROR) << "Please recreate your database or add a .dbtype file to your sequence database\n";
        return EXIT_FAILURE;
    }
    if (Parameters::isEqualDbtype(seqReader.getDbtype(), Parameters::DBTYPE_HMM_PROFILE)) {
        Debug(Debug::ERROR) << "The sequence database must not be a profile database\n";
        return EXIT_FAILURE;
    }

    Debug(Debug::INFO) << "Scoring observed members by profile PSSM, minimum coverage: " << par.covThr << "\n";

    // set if any alignment record lacks a backtrace and had to be recomputed
    bool missingBacktrace = false;

    Debug::Progress progress(resultReader.getSize());
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        Matcher matcher(seqReader.getDbtype(), maxSequenceLength, &subMat, &evalueComputation,
                        par.compBiasCorrection, par.compBiasCorrectionScale,
                        par.gapOpen.values.aminoacid(), par.gapExtend.values.aminoacid(), 0.0, par.zdrop);
        MultipleAlignment aligner(maxSequenceLength, &subMat);
        PSSMCalculator calculator(
            &subMat, maxSequenceLength, maxSetSize, par.pcmode, par.pca, par.pcb
#ifdef GAP_POS_SCORING
            , par.gapOpen.values.aminoacid()
            , par.gapPseudoCount
#endif
        );
        Sequence centerSequence(maxSequenceLength, seqReader.getDbtype(), &subMat, 0, false, par.compBiasCorrection);
        Sequence edgeSequence(maxSequenceLength, seqReader.getDbtype(), &subMat, 0, false, false);

        char dbKey[255];
        char buffer[1024];

        std::vector<Matcher::result_t> alnResults;
        alnResults.reserve(300);
        std::vector<std::vector<unsigned char>> seqSet;
        seqSet.reserve(300);

        std::string result;
        result.reserve(64);

#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < resultReader.getSize(); id++) {
            progress.updateProgress();

            const DBKeyType queryKey = resultReader.getDbKey(id);
            char *data = resultReader.getData(id, thread_idx);
            if (*data == '\0') {
                // no alignments for this cluster: keep the representative so that no
                // cluster is lost when the mapping is used to rewrite the cluster DB
                result.clear();
                result.append(SSTR(queryKey));
                result.append("\t0\t1.0000\t1\n");
                resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx);
                continue;
            }

            const size_t queryId = seqReader.getId(queryKey);
            if (queryId == DB_ENTRY_NOT_FOUND) {
                Debug(Debug::WARNING) << "Invalid representative sequence " << queryKey << "\n";
                continue;
            }
            centerSequence.mapSequence(queryId, queryKey, seqReader.getData(queryId, thread_idx), seqReader.getSeqLen(queryId));

            bool isQueryInit = false;
            while (*data != '\0') {
                Util::parseKey(data, dbKey);
                const DBKeyType key = Util::fast_atoi<DBKeyType>(dbKey);
                // the representative is the center (row 0); skip its self-hit line
                if (key == queryKey) {
                    data = Util::skipLine(data);
                    continue;
                }

                const size_t edgeId = seqReader.getId(key);
                if (edgeId == DB_ENTRY_NOT_FOUND) {
                    Debug(Debug::ERROR) << "Sequence " << key << " does not exist in the sequence database\n";
                    EXIT(EXIT_FAILURE);
                }
                edgeSequence.mapSequence(edgeId, key, seqReader.getData(edgeId, thread_idx), seqReader.getSeqLen(edgeId));
                seqSet.emplace_back(std::vector<unsigned char>(edgeSequence.numSequence, edgeSequence.numSequence + edgeSequence.L));

                const char *entry[255];
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns > Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    alnResults.emplace_back(Matcher::parseAlignmentRecord(data));
                } else {
                    // Reuse is impossible without a backtrace: recompute the cheap
                    // member-vs-center alignment (never profile-vs-member).
                    missingBacktrace = true;
                    if (isQueryInit == false) {
                        matcher.initQuery(&centerSequence);
                        isQueryInit = true;
                    }
                    alnResults.emplace_back(matcher.getSWResult(&edgeSequence, INT_MAX, false, 0, 0.0, FLT_MAX, Matcher::SCORE_COV_SEQID, 0, false));
                }
                data = Util::skipLine(data);
            }

            // singleton cluster: the representative stays the representative
            if (seqSet.empty()) {
                result.clear();
                result.append(SSTR(queryKey));
                result.append("\t0\t1.0000\t1\n");
                resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx);
                alnResults.clear();
                continue;
            }

            MultipleAlignment::MSAResult res = aligner.computeMSA(&centerSequence, seqSet, alnResults, true);

            PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(res.setSize, res.centerLength,
                                                                            (const char **) res.msaSequence,
#ifdef GAP_POS_SCORING
                                                                            alnResults,
#endif
                                                                            par.wg, 0.0);

            // score every observed member row against the profile
            DBKeyType bestKey = queryKey;
            double bestScore = -DBL_MAX;
            float bestCov = 0.0f;
            int bestOrigScore = INT_MAX;
            bool bestIsRep = true;
            bool anyPassed = false;

            // deterministic fallback if nothing meets the coverage threshold
            DBKeyType repKey = queryKey;
            double repScore = -DBL_MAX;
            float repCov = 0.0f;

            for (size_t row = 0; row < res.setSize; row++) {
                const char *msaRow = res.msaSequence[row];
                const DBKeyType candKey = (row == 0) ? queryKey : alnResults[row - 1].dbKey;
                const int origScore = (row == 0) ? INT_MAX : alnResults[row - 1].score;

                double score = 0.0;
                size_t aligned = 0;
                for (size_t pos = 0; pos < res.centerLength; pos++) {
                    const unsigned char state = (unsigned char) msaRow[pos];
                    if (state < Sequence::PROFILE_AA_SIZE) {
                        aligned++;
                        score += (int) pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + state];
                    }
                }
                const float coverage = (res.centerLength > 0) ? (float) aligned / (float) res.centerLength : 0.0f;

                if (row == 0) {
                    repKey = candKey;
                    repScore = score;
                    repCov = coverage;
                }

                if (coverage < par.covThr) {
                    continue;
                }

                // tie-break (deterministic): score, coverage, original alignment score,
                // prefer the old representative, then smallest key
                bool better = false;
                if (anyPassed == false) {
                    better = true;
                } else if (score != bestScore) {
                    better = score > bestScore;
                } else if (coverage != bestCov) {
                    better = coverage > bestCov;
                } else if (origScore != bestOrigScore) {
                    better = origScore > bestOrigScore;
                } else if ((row == 0) != bestIsRep) {
                    better = (row == 0);
                } else {
                    better = candKey < bestKey;
                }
                if (better) {
                    bestKey = candKey;
                    bestScore = score;
                    bestCov = coverage;
                    bestOrigScore = origScore;
                    bestIsRep = (row == 0);
                    anyPassed = true;
                }
            }

            if (anyPassed == false) {
                // no candidate reached the coverage threshold: keep the old representative
                bestKey = repKey;
                bestScore = repScore;
                bestCov = repCov;
            }

            result.clear();
            result.append(SSTR(bestKey));
            result.push_back('\t');
            snprintf(buffer, sizeof(buffer), "%d\t%.4f\t%zu\n", (int) bestScore, bestCov, res.setSize);
            result.append(buffer);
            resultWriter.writeData(result.c_str(), result.length(), queryKey, thread_idx);

            MultipleAlignment::deleteMSA(&res);
            seqSet.clear();
            alnResults.clear();
        }
    }
    resultWriter.close();
    resultReader.close();
    seqReader.close();

    if (missingBacktrace) {
        Debug(Debug::WARNING) << "Some alignment records had no backtrace and were recomputed on the fly.\n"
                              << "This is slower and can happen if the alignment step was run without '-a'.\n"
                              << "Re-run the alignment with '-a' to store backtraces if this was unintended.\n";
    }

    return EXIT_SUCCESS;
}
