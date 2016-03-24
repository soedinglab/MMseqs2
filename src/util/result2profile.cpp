// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include <string>
#include <vector>
#include <sstream>

#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Log.h"
#include "Util.h"
#include "Debug.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

enum {
    MSA = 0,
    PSSM
};

size_t findMaxSetSize(DBReader<unsigned int> *reader) {
    //Find the max set size
    size_t maxSetSize = 0;
    for (size_t i = 0; i < reader->getSize(); i++) {
        char *data = reader->getData(i);
        size_t entry = 0;
        size_t position = 0;
        while (data[position] != '\0') {
            if (data[position] == '\n') {
                entry++;
            }
            position++;
        }
        maxSetSize = std::max(maxSetSize, entry);
    }
    return maxSetSize + 1;
}

MultipleAlignment::MSAResult computeAlignment(MultipleAlignment &aligner, Sequence *centerSequence,
                                              std::vector<Sequence *> seqSet,
                                              std::vector<Matcher::result_t> alnResults,
                                              bool allowDeletion, bool sameDatabase) {
    if (alnResults.size() > 0) {
        std::vector<Matcher::result_t> alnWithoutIdentity;
        for (size_t i = 0; i < alnResults.size(); i++) {
            if (alnResults[i].dbKey == centerSequence->getDbKey() && sameDatabase == true) { ;
            } else {
                alnWithoutIdentity.push_back(alnResults[i]);
            }
        }
        return aligner.computeMSA(centerSequence, seqSet, alnWithoutIdentity, !allowDeletion);
    } else {
        return aligner.computeMSA(centerSequence, seqSet, !allowDeletion);
    }
}

int result2outputmode(Parameters &par, int mode) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> *qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str());
    qDbr->open(DBReader<unsigned int>::NOSORT);

    std::string headerName(par.db1);
    headerName.append("_h");

    std::string headerIndexName(par.db1);
    headerIndexName.append("_h.index");

    DBReader<unsigned int> *queryHeaderReader = new DBReader<unsigned int>(headerName.c_str(), headerIndexName.c_str());
    queryHeaderReader->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *tDbr = qDbr;
    DBReader<unsigned int> *tempateHeaderReader = queryHeaderReader;

    size_t maxSequenceLength = 0;
    unsigned int *lengths = qDbr->getSeqLens();
    for (size_t i = 0; i < qDbr->getSize(); i++) {
        maxSequenceLength = std::max((size_t) lengths[i], maxSequenceLength);
    }

    bool sameDatabase = true;
    if (par.db1.compare(par.db2) != 0) {
        sameDatabase = false;
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
        tDbr->open(DBReader<unsigned int>::NOSORT);

        lengths = tDbr->getSeqLens();
        for (size_t i = 0; i < tDbr->getSize(); i++) {
            maxSequenceLength = std::max((size_t) lengths[i], maxSequenceLength);
        }

        headerName = par.db2;
        headerName.append("_h");

        headerIndexName = par.db2;
        headerIndexName.append("_h.index");

        tempateHeaderReader = new DBReader<unsigned int>(headerName.c_str(), headerIndexName.c_str());
        tempateHeaderReader->open(DBReader<unsigned int>::NOSORT);
    }

    DBReader<unsigned int> *resultReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str());
    resultReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);

//    FileUtil::errorIfFileExist(par.db4.c_str());
//    FileUtil::errorIfFileExist(par.db4Index.c_str());

    DBWriter profileWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBWriter::BINARY_MODE);
    profileWriter.open();
    DBWriter * concensusWriter;
    if(mode == PSSM) {
        concensusWriter = new DBWriter(std::string(par.db4 + "_consensus").c_str(),
                                       std::string(par.db4 + "_consensus.index").c_str(),
                                       par.threads, DBWriter::ASCII_MODE);
        concensusWriter->open();
    }
    size_t maxSetSize = findMaxSetSize(resultReader);

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);
    Debug(Debug::INFO) << "Start computing " << (!mode ? "MSAs" : "profiles") << ".\n";
#pragma omp parallel
    {
        Matcher matcher(maxSequenceLength, &subMat, tDbr->getAminoAcidDBSize(), tDbr->getSize(),
                        par.compBiasCorrection);

        MultipleAlignment aligner(maxSequenceLength, maxSetSize, &subMat, &matcher);
        PSSMCalculator calculator(&subMat, maxSequenceLength, par.pca, par.pcb);
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat);
        int sequenceType = Sequence::AMINO_ACIDS;
        if (par.profile == true) {
            sequenceType = Sequence::HMM_PROFILE;
        }
        Sequence *centerSequence = new Sequence(maxSequenceLength, subMat.aa2int, subMat.int2aa,
                                                sequenceType, 0, false, par.compBiasCorrection);

#pragma omp for schedule(static)
        for (size_t id = 0; id < resultReader->getSize(); id++) {
            Log::printProgress(id);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char *results = resultReader->getData(id);
            char dbKey[255 + 1];

            // check if already aligned results exists
            std::vector<Matcher::result_t> alnResults;
            char *entry[255];
            size_t columns = Util::getWordsOfLine(results, entry, 255);
            if (columns == Matcher::ALN_RES_WITH_BT_COL_CNT) {
                alnResults = Matcher::readAlignmentResults(results);
            }

            unsigned int queryKey = resultReader->getDbKey(id);
            char *seqData = qDbr->getDataByDBKey(queryKey);
            centerSequence->mapSequence(0, queryKey, seqData);

            std::vector<Sequence *> seqSet;
            while (*results != '\0') {
                Util::parseKey(results, dbKey);
                unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                // just add sequences if its not the same as the query in case of sameDatabase
                if (key != queryKey || sameDatabase == false) {
                    size_t edgeId = tDbr->getId(key);
                    char *dbSeqData = tDbr->getData(edgeId);
                    Sequence *edgeSequence = new Sequence(tDbr->getSeqLens(edgeId), subMat.aa2int, subMat.int2aa,
                                                          Sequence::AMINO_ACIDS, 0, false, false);
                    edgeSequence->mapSequence(0, key, dbSeqData);
                    seqSet.push_back(edgeSequence);
                }
                results = Util::skipLine(results);
            }

            MultipleAlignment::MSAResult res = computeAlignment(aligner, centerSequence, seqSet,
                                                                alnResults, par.allowDeletion, sameDatabase);
            //MultipleAlignment::print(res, &subMat);
            MsaFilter::MsaFilterResult filterRes = filter.filter((const char **) res.msaSequence, res.setSize,
                                                                 res.centerLength, static_cast<int>(par.cov * 100),
                                                                 static_cast<int>(par.qid * 100), par.qsc,
                                                                 static_cast<int>(par.filterMaxSeqId * 100), par.Ndiff);

            char *data;
            size_t dataSize;
            std::stringstream msa;
            std::string result;
            switch (mode) {
                case MSA: {

                    for (size_t i = 0; i < filterRes.setSize; i++) {
                        unsigned int key;
                        char *data;
                        if (i == 0) {
                            key = queryKey;
                            data = queryHeaderReader->getDataByDBKey(key);
                        } else {
                            key = seqSet[i - 1]->getDbKey();
                            data = tempateHeaderReader->getDataByDBKey(key);
                        }
                        if (par.addInternalId) {
                            msa << "#" << key << "\n";
                        }
                        msa << ">" << data;
                        for (size_t pos = 0; pos < res.centerLength; pos++) {
                            char aa = filterRes.filteredMsaSequence[i][pos];
                            msa << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int) aa] : '-');
                        }
                        msa << "\n";
                    }
                    result = msa.str();
                    data = (char *) result.c_str();
                    dataSize = result.length();
                }
                    break;
                case PSSM: {
//                    std::cout << centerSequence->getDbKey() << " " << res.setSize << std::endl;;
//                    for (size_t i = 0; i < res.setSize; i++) {
//                        for (size_t pos = 0; pos < res.msaSequenceLength; pos++) {
//                            char aa = res.msaSequence[i][pos];
//                            std::cout << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int) aa] : '-');
//                        }
//                        std::cout << std::endl;
//                    }
//
//                    if (par.pruneSequences) {
//                        filter.pruneAlignment((char **) res.msaSequence, res.setSize, res.centerLength);
//                    }
//
//                    std::cout << centerSequence->getDbKey() << " " << res.setSize << " " << filterRes.setSize <<
//                    std::endl;
//                    for (size_t i = 0; i < filterRes.setSize; i++) {
//                        for (size_t pos = 0; pos < res.msaSequenceLength; pos++) {
//                            char aa = filterRes.filteredMsaSequence[i][pos];
//                            std::cout << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int) aa] : '-');
//                        }
//                        std::cout << std::endl;
//                    }

/*                  std::cout << centerSequence->getDbKey() << " " << res.setSize << " " << filterRes.setSize << std::endl;
		    for (size_t i = 0; i < filterRes.setSize; i++) {
                       for(size_t pos = 0; pos < res.msaSequenceLength; pos++){
                              char aa = filterRes.filteredMsaSequence[i][pos];
                              std::cout << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int)aa] : '-');
                       }
                       std::cout << std::endl;
                    }
*/		    
		            data = (char *) calculator.computePSSMFromMSA(filterRes.setSize, res.centerLength,
                                                                  filterRes.filteredMsaSequence, par.wg);
                    dataSize = res.centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);

                    std::string consensusStr;
                    for(size_t i = 0; i < res.centerLength; i++){
                        int maxScore = INT_MIN;
                        int maxAA = 0;
                        for(size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++){
                            if(static_cast<int>(data[i * Sequence::PROFILE_AA_SIZE + aa ]) > maxScore){
                                maxScore = data[i * Sequence::PROFILE_AA_SIZE + aa ];
                                maxAA = aa;
                            }
                        }
                        consensusStr.push_back(subMat.int2aa[maxAA]);
                    }
                    consensusStr.push_back('\n');
                    concensusWriter->write(consensusStr.c_str(), consensusStr.length(), SSTR(queryKey).c_str(), thread_idx);

                    for (size_t i = 0; i < dataSize; i++) {
                        // Avoid a null byte result
                        data[i] = data[i] ^ 0x80;
                    }

                    //pssm.printProfile(res.centerLength);
                    //calculator.printPSSM(res.centerLength);
                }
                    break;
                default:
                    EXIT(EXIT_FAILURE);
            }

            profileWriter.write(data, dataSize, SSTR(queryKey).c_str(), thread_idx);

            MultipleAlignment::deleteMSA(&res);
            for (std::vector<Sequence *>::iterator it = seqSet.begin(); it != seqSet.end(); ++it) {
                Sequence *seq = *it;
                delete seq;
            }
        }
        delete centerSequence;
    }

    // cleanup
    profileWriter.close();
    if(mode == PSSM){
        concensusWriter->close();
        delete concensusWriter;
    }
    resultReader->close();
    delete resultReader;

    if (!sameDatabase) {
        tempateHeaderReader->close();
        tDbr->close();
        delete tempateHeaderReader;
        delete tDbr;
    }

    queryHeaderReader->close();
    qDbr->close();
    delete queryHeaderReader;
    delete qDbr;

    Debug(Debug::INFO) << "\nDone.\n";

    return EXIT_SUCCESS;
}


int result2profile(int argc, const char **argv) {
    std::string usage("Calculates profiles from a clustering.\n");
    usage.append("USAGE: <queryDB> <targetDB> <resultDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>\n");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.result2profile, 4);

    // never allow deletions
    par.allowDeletion = false;

    return result2outputmode(par, PSSM);
}

int result2msa(int argc, const char **argv) {
    std::string usage("Calculates MSAs from a clustering.\n");
    usage.append("USAGE: <queryDB> <targetDB> <resultDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@mpibpc.mpg.de>\n");
    usage.append("\t & Milot Mirdita <milot@mirdita.de>");

    Parameters par;
    par.parseParameters(argc, argv, usage, par.result2msa, 4);

    return result2outputmode(par, MSA);
}
