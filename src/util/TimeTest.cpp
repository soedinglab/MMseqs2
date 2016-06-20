#include "QueryTemplateMatcherLocal.h"
#include "DBReader.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "Prefiltering.h"

#include <sys/time.h>

#ifdef OPENMP
#include <omp.h>
#endif

void doTimeTest(const std::string &targetDB,
                const std::string &targetDBIndex,
                const std::string &scoringMatrixFile,
                size_t maxSeqLen,
                const std::string &logFile,
                int threads) {
    Sequence **seqs = new Sequence *[threads];

    size_t QUERY_SET_SIZE = 50000;

    std::ofstream logFileStream;
    logFileStream.open(logFile);

    QueryTemplateMatcher **matchers = new QueryTemplateMatcher *[threads];
    DBReader<unsigned int> tdbr(targetDB.c_str(), targetDBIndex.c_str());
    tdbr.open(DBReader<unsigned int>::NOSORT);

    size_t targetSeqLenSum = 0;
    for (size_t i = 0; i < tdbr.getSize(); i++)
        targetSeqLenSum += tdbr.getSeqLens(i);

    // generate a small random sequence set for testing
    size_t querySetSize = tdbr.getSize();
    if (querySetSize > QUERY_SET_SIZE)
        querySetSize = QUERY_SET_SIZE;

    size_t *querySeqs = new size_t[querySetSize];
    srand(1);
    size_t querySeqLenSum = 0;
    for (size_t i = 0; i < querySetSize; i++) {
        querySeqs[i] = rand() % tdbr.getSize();
        querySeqLenSum += tdbr.getSeqLens(querySeqs[i]);
    }

    short kmerThrPerPosMin = 1;
    short kmerThrPerPosMax = 25;

    SubstitutionMatrix m(scoringMatrixFile.c_str(), 8.0, -0.2);

    // adjust k-mer list length threshold
    for (size_t alphabetSize = 21; alphabetSize <= 21; alphabetSize += 4) {
        Debug(Debug::INFO) << "Target database: " << targetDB << "(size=" << tdbr.getSize() << ")\n";

        BaseMatrix *subMat;
        if (alphabetSize < 21)
            subMat = new ReducedMatrix(m.probMatrix, m.subMatrixPseudoCounts, alphabetSize, 8.0);
        else
            subMat = &m;

        ExtendedSubstitutionMatrix *_2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, alphabetSize);
        ExtendedSubstitutionMatrix *_3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, alphabetSize);


        for (int isSpaced = 1; isSpaced < 2; isSpaced++) {
            for (int isLocal = 1; isLocal < 2; isLocal++) {
                for (unsigned int kmerSize = 6; kmerSize <= 7; kmerSize++) {
#pragma omp parallel for schedule(static)
                    for (int i = 0; i < threads; i++) {
                        int thread_idx = 0;
#ifdef OPENMP
                        thread_idx = omp_get_thread_num();
#endif
                        // TODO: isSpaced is probably wrong
                        seqs[thread_idx] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa,
                                                        Sequence::AMINO_ACIDS, kmerSize, isSpaced, false);
                    }

                    short kmerThrMin = (short) ((float) (kmerThrPerPosMin * kmerSize) *
                                                (pow(((float) alphabetSize / 21.0), 2.0)));
                    int kmerThrMax = kmerThrPerPosMax * kmerSize;

                    Debug(Debug::INFO) << "a = " << alphabetSize << ",  k = " << kmerSize << "\n";
                    IndexTable *indexTable = Prefiltering::generateIndexTable(&tdbr, seqs[0], subMat, alphabetSize, kmerSize, 0,
                                                                              tdbr.getSize(), isLocal, 0);

                    short decr = 1;
                    if (kmerSize == 6 || kmerSize == 7)
                        decr = 2;

                    Debug(Debug::INFO) << "Omitting runs with too short running time...\n";
                    for (int kmerThr = kmerThrMax; kmerThr >= kmerThrMin; kmerThr -= decr) {
                        size_t dbMatchesSum = 0;
                        size_t doubleMatches = 0;
                        double kmersPerPos = 0.0;

                        // determine k-mer match probability and k-mer list length for kmerThr
#pragma omp parallel for schedule(static)
                        for (int i = 0; i < threads; i++) {
                            int thread_idx = 0;
#ifdef OPENMP
                            thread_idx = omp_get_thread_num();
#endif
                            // set a current k-mer list length threshold and a high prefitlering threshold (we don't need the prefiltering results in this test run)
                            if (isLocal == 1) {
                                matchers[thread_idx] = new QueryTemplateMatcherLocal(subMat, indexTable,
                                                                                     tdbr.getSeqLens(), kmerThr, 1.0,
                                                                                     kmerSize,
                                                                                     seqs[0]->getEffectiveKmerSize(),
                                                                                     tdbr.getSize(), maxSeqLen, 300);
                            }
                            matchers[thread_idx]->setSubstitutionMatrix(_3merSubMatrix->scoreMatrix,
                                                                        _2merSubMatrix->scoreMatrix);
                        }

                        struct timeval start, end;
                        gettimeofday(&start, NULL);
#pragma omp parallel for schedule(dynamic, 10) reduction (+: dbMatchesSum, kmersPerPos, doubleMatches)
                        for (size_t i = 0; i < querySetSize; i++) {
                            size_t id = querySeqs[i];

                            int thread_idx = 0;
#ifdef OPENMP
                            thread_idx = omp_get_thread_num();
#endif
                            char *seqData = tdbr.getData(id);
                            seqs[thread_idx]->mapSequence(id, tdbr.getDbKey(id), seqData);

                            matchers[thread_idx]->matchQuery(seqs[thread_idx], UINT_MAX);

                            kmersPerPos += matchers[thread_idx]->getStatistics()->kmersPerPos;
                            dbMatchesSum += matchers[thread_idx]->getStatistics()->dbMatches;
                            doubleMatches += matchers[thread_idx]->getStatistics()->doubleMatches;
                        }
                        gettimeofday(&end, NULL);
                        ssize_t sec = end.tv_sec - start.tv_sec;

                        // too short running time is not recorded
                        if (sec <= 2) {
                            //Debug(Debug::WARNING) << "Time <= 2 sec, not counting this step.\n\n";
                            continue;
                        }

                        Debug(Debug::INFO) << "spaced = " << isSpaced << ", local = " << isLocal << ", k = " << kmerSize << ", a = " << alphabetSize << "\n";
                        Debug(Debug::INFO) << "k-mer threshold = " << kmerThr << "\n";
                        Debug(Debug::INFO) << "Time: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

                        kmersPerPos /= (double) querySetSize;
                        double kmerMatchProb;
                        if (isLocal == 1) {
                            kmerMatchProb =
                                    (((double) doubleMatches) / ((double) (querySeqLenSum * targetSeqLenSum))) / 256;
                        } else {
                            kmerMatchProb = ((double) dbMatchesSum) / ((double) (querySeqLenSum * targetSeqLenSum));
                        }

                        Debug(Debug::INFO) << "kmerPerPos: " << kmersPerPos << "\n";
                        Debug(Debug::INFO) << "k-mer match probability: " << kmerMatchProb << "\n";
                        Debug(Debug::INFO) << "dbMatch: " << dbMatchesSum << "\n";
                        Debug(Debug::INFO) << "doubleMatches: " << doubleMatches << "\n";
                        Debug(Debug::INFO) << "querySeqLenSum: " << querySeqLenSum << "\n";
                        Debug(Debug::INFO) << "targetSeqLenSum: " << targetSeqLenSum << "\n\n";

                        logFileStream << kmersPerPos << "\t" << kmerMatchProb << "\t" << dbMatchesSum << "\t" <<
                        doubleMatches << "\t" << kmerSize << "\t" << alphabetSize << "\t" << isSpaced << "\t" <<
                        isLocal << "\t" << sec << "\n";

                        // running time for the next step will be too long
                        if (sec >= 1200) {
                            Debug(Debug::WARNING) << "Time >= 300 sec, going to the next parameter combination.\n";
                            break;
                        }

                        for (int j = 0; j < threads; j++) {
                            delete matchers[j];
                        }
                    }
                    delete indexTable;
                }
            }
        }

        for (int i = 0; i < threads; i++) {
            delete seqs[i];
        }

        if (alphabetSize < 21) {
            delete subMat;
        }

        delete _2merSubMatrix;
        delete _3merSubMatrix;
    }

    tdbr.close();
    logFileStream.close();
    delete[] querySeqs;
    delete[] matchers;
    delete[] seqs;
}

int timetest(int argc, const char *argv[]) {
    std::string usage("\nRuns time benchmark on a database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: <dbBaseIn> <logOut>");

    Parameters par;
    std::vector<MMseqsParameter> options;
    options.push_back(par.PARAM_V);
    options.push_back(par.PARAM_SUB_MAT);
    options.push_back(par.PARAM_MAX_SEQ_LEN);
    options.push_back(par.PARAM_THREADS);

    par.parseParameters(argc, argv, usage, options, 2);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    doTimeTest(par.db1,
               par.db1Index,
               par.scoringMatrixFile, par.maxSeqLen,
               par.db2,
               par.threads);

    return EXIT_SUCCESS;
}
