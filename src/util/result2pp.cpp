#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"
#include "Matcher.h"
#include "SubstitutionMatrix.h"

#ifdef OPENMP
#include <omp.h>
#endif

float computeNeff(float neffA, float maxNeffA, float neffB, float maxNeffB, float avgNewNeff) {
    float w = (neffA + neffB) / (maxNeffA + maxNeffB);
    //std::cout<<"COmputing new neff:"<<neffA<<","<<maxNeffA<<","<<neffB<<","<<maxNeffB<<","<<avgNewNeff<<","<<w<<","<<std::endl;
    return avgNewNeff + 1 - exp(log(avgNewNeff) * (1 - w));
}

int result2pp(int argc, const char **argv, const Command& command) {
    MMseqsMPI::init(argc, argv);

    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, 0, 0);
    par.evalProfile = (par.evalThr < par.evalProfile) ? par.evalThr : par.evalProfile;
    par.printParameters(command.cmd, argc, argv, *command.params);

    DBReader<unsigned int> qDbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    qDbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *tDbr = &qDbr;
    bool sameDatabase = true;
    if (par.db1.compare(par.db2) != 0) {
        sameDatabase = false;
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        tDbr->open(DBReader<unsigned int>::NOSORT);
    }

    DBReader<unsigned int> resultReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t dbFrom = 0;
    size_t dbSize = resultReader.getSize();
    std::string outDb = par.db4;
    std::string outIndex = par.db4Index;
#ifdef HAVE_MPI
    resultReader.decomposeDomainByAminoAcid(MMseqsMPI::rank, MMseqsMPI::numProc, &dbFrom, &dbSize);
    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << (dbFrom + dbSize) << "\n";
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(par.db4, par.db4Index, MMseqsMPI::rank);
    outDb = tmpOutput.first;
    outIndex = tmpOutput.second;
#endif

    DBWriter resultWriter(outDb.c_str(), outIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_HMM_PROFILE);
    resultWriter.open();

    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0f, 0.0f);

    Debug::Progress progress(dbSize - dbFrom);
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

        Sequence queryProfile(par.maxSeqLen, qDbr.getDbtype(), &subMat, 0, false, par.compBiasCorrection, false);
        Sequence targetProfile(par.maxSeqLen, tDbr->getDbtype(), &subMat, 0, false, par.compBiasCorrection, false);
        float *outProfile = new float[(par.maxSeqLen + 1) * Sequence::PROFILE_AA_SIZE];
        float *neffM = new float[par.maxSeqLen + 1];

        std::string result;
        result.reserve((par.maxSeqLen + 1) * Sequence::PROFILE_READIN_SIZE);

        char dbKey[255 + 1];
        const char *entry[255];

#pragma omp for schedule(dynamic, 10)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            progress.updateProgress();

            unsigned int queryKey = resultReader.getDbKey(id);
            size_t queryId = qDbr.getId(queryKey);
            char *queryData = qDbr.getData(queryId, thread_idx);

            char *results = resultReader.getData(id, thread_idx);
            if (*results == '\0') {
                resultWriter.writeData(queryData, qDbr.getEntryLen(queryId) - 1, queryKey, thread_idx);
                continue;
            }

            queryProfile.mapSequence(queryId, queryKey, queryData, qDbr.getSeqLen(queryId));
            const float *qProfile = queryProfile.getProfile();

            /*const size_t profile_row_size = queryProfile.profile_row_size;
            // init outProfile with query Probs
            for(int l = 0; l < queryProfile.L; l++) {
                for (size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    outProfile[l*Sequence::PROFILE_AA_SIZE + aa_num] = qProfile[l * profile_row_size + aa_num];
                }
            }
            */
            memset(outProfile, 0, queryProfile.L * Sequence::PROFILE_AA_SIZE * sizeof(float));
            float maxNeffQ = 0;
            for (int pos = 0; pos < queryProfile.L; pos++) {
                maxNeffQ = std::max(maxNeffQ, queryProfile.neffM[pos]);
                neffM[pos] = queryProfile.neffM[pos];
            }

            int minqStartPos = INT_MAX;
            int maxqEndPos = 0;
            bool didMerge = false;
            while (*results != '\0') {
                Util::parseKey(results, dbKey);
                const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                if (columns <= Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    Debug(Debug::ERROR) << "Alignment must contain the alignment information. Compute the alignment with option -a.\n";
                    EXIT(EXIT_FAILURE);
                }

                // just add sequences if eval < thr. and if key is not the same as the query in case of sameDatabase
                double evalue = strtod(entry[3], NULL);
                if (evalue <= par.evalProfile && (key != queryKey || sameDatabase == false)) {
                    didMerge = true;
                    const Matcher::result_t res = Matcher::parseAlignmentRecord(results);
                    const size_t edgeId = tDbr->getId(key);
                    targetProfile.mapSequence(edgeId, key, tDbr->getData(edgeId, thread_idx), tDbr->getSeqLen(edgeId));
                    const float *tProfile = targetProfile.getProfile();
                    size_t qPos = res.qStartPos;
                    minqStartPos = std::min(minqStartPos, res.qStartPos);
                    maxqEndPos = std::max(maxqEndPos, res.qEndPos);
                    size_t tPos = res.dbStartPos;
                    size_t aliLength = 0;
                    float avgEntropy = 0.0f;
                    float maxNeffT = 0;
                    for (int pos = 0; pos < targetProfile.L; pos++) {
                        maxNeffT = std::max(maxNeffT, targetProfile.neffM[pos]);
                    }

                    for (size_t btPos = 0; btPos < res.backtrace.size(); btPos++) {
                        aliLength++;
                        char letter = res.backtrace[btPos];
                        for (size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                            // TODO do all alignment states contribute?
//                            if (letter == 'M') {
                            float qProb = qProfile[qPos * Sequence::PROFILE_AA_SIZE + aa_num] * queryProfile.neffM[qPos];
                            float tProb = tProfile[tPos * Sequence::PROFILE_AA_SIZE + aa_num] * targetProfile.neffM[tPos];
                            float mixedProb = qProb + tProb;
                            outProfile[qPos * Sequence::PROFILE_AA_SIZE + aa_num] += mixedProb;
                            mixedProb /= queryProfile.neffM[qPos] + targetProfile.neffM[tPos];
                            avgEntropy += (mixedProb > 0.0) ? (-mixedProb * log(mixedProb)) : 0.0;
//                            } else {
//                                outProfile[qPos * Sequence::PROFILE_AA_SIZE + aa_num] = qProfile[qPos * Sequence::PROFILE_AA_SIZE + aa_num];
//                            }
                        }

                        if (letter == 'M') {
                            qPos++;
                            tPos++;
                        } else if (letter == 'I') {
                            qPos++;
                        } else if (letter == 'D') {
                            tPos++;
                        }
                    }

                    // Normalize probability
                    for (int l = res.qStartPos; l < res.qEndPos; l++) {
                        MathUtil::NormalizeTo1(&outProfile[l * Sequence::PROFILE_AA_SIZE], Sequence::PROFILE_AA_SIZE);
                    }
                    avgEntropy /= aliLength;
                    float avgNewNeff = exp(avgEntropy);

                    // update the Neff of the merge between the target prof and the query prof
                    qPos = res.qStartPos;
                    tPos = res.dbStartPos;
                    for (size_t btPos = 0; btPos < res.backtrace.size(); btPos++) {
                        char letter = res.backtrace[btPos];
                        neffM[qPos] = computeNeff(queryProfile.neffM[qPos], maxNeffQ, targetProfile.neffM[tPos], maxNeffT, avgNewNeff);

                        if (letter == 'M') {
                            qPos++;
                            tPos++;
                        } else if (letter == 'I') {
                            ++qPos;
                        } else if (letter == 'D') {
                            ++tPos;
                        }
                    }
                }
                results = Util::skipLine(results);
            }

            if (didMerge == false) {
                resultWriter.writeData(queryData, qDbr.getEntryLen(queryId) - 1, queryKey, thread_idx);
                continue;
            }

            for (int l = 0; l < minqStartPos; l++) {
                for (size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    outProfile[l*Sequence::PROFILE_AA_SIZE + aa_num] = qProfile[l * Sequence::PROFILE_AA_SIZE + aa_num];
                }
            }
            for (int l = maxqEndPos; l < queryProfile.L; l++) {
                for (size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num] = qProfile[l * Sequence::PROFILE_AA_SIZE + aa_num];
                }
            }
            /*
            float maxNewNeff = 0.0;
            float avgEntropy = 0.0;
            for(int l = 0; l < queryProfile.L; l++) {
                maxNewNeff = std::max(maxNewNeff,neffM[l]);
                for(size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    float mixedProb = outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num];
                    avgEntropy += -mixedProb * log(mixedProb);
                }
            }
            avgEntropy /= queryProfile.L;

            float avgNewNeff = exp(avgEntropy);

            // update the Neff of the merges of all target prof and the query prof
            for(int l = 0; l < queryProfile.L; l++) {
                //for(size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    neffM[l] = computeNeff(queryProfile.neffM[l], maxNeffQ, neffM[l],maxNewNeff,avgNewNeff);
                //}
            } */
//            size_t pos = 0;
            std::string consensus(queryProfile.L, 'X');
            for (int l = 0; l < queryProfile.L; l++) {
                float maxProb = std::numeric_limits<float>::lowest();
                for (size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    result.push_back(Sequence::scoreMask(outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num]));
                    //std::cout<< outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num]<<"\t";
                    if (outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num] > maxProb) {
                        consensus[l] = aa_num;
                        maxProb = outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num];
                    }
                }
                //std::cout<<std::endl;

                // write query, consensus sequence and neffM
                result.push_back(queryProfile.numSequence[l]);
                result.push_back(consensus[l]);
                unsigned char neff = MathUtil::convertNeffToChar(neffM[l]);
                result.push_back(neff);
            }
            //std::cout<<"Query length:"<<queryProfile.L<<", res length:"<<result.size()<<", should be:"<<queryProfile.L*23<<std::endl;
            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx);
            result.clear();
        }
        delete[] outProfile;
        delete[] neffM;
    }
    resultWriter.close(true);
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // master reduces results
    if (MMseqsMPI::isMaster()) {
        std::vector<std::pair<std::string, std::string>> splitFiles;
        for (int procs = 0; procs < MMseqsMPI::numProc; procs++) {
            splitFiles.emplace_back(Util::createTmpFileNames(par.db4, par.db4Index, procs));

        }
        DBWriter::mergeResults(par.db4, par.db4Index, splitFiles);
    }
#endif
    resultReader.close();
    if (!sameDatabase) {
        tDbr->close();
        delete tDbr;
    }
    qDbr.close();

    return EXIT_SUCCESS;
}
