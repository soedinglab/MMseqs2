// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include <string>
#include <vector>
#include <sstream>
#include <MsaFilter.h>

#include "Matcher.h"
#include "SubstitutionMatrix.h"
#include "Parameters.h"
#include "Sequence.h"
#include "MultipleAlignment.h"
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


size_t findMaxSetSize(DBReader<unsigned int>* reader) {
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

int result2outputmode(Parameters par, int mode) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int>*qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str());
    qDbr->open(DBReader<unsigned int>::NOSORT);

    std::string headerName(par.db1);
    headerName.append("_h");

    std::string headerIndexName(par.db1);
    headerIndexName.append("_h.index");

    DBReader<unsigned int>* queryHeaderReader = new DBReader<unsigned int>(headerName.c_str(), headerIndexName.c_str());
    queryHeaderReader->open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int>*tDbr = qDbr;
    DBReader<unsigned int>* tempateHeaderReader = queryHeaderReader;

    bool sameDatabase = true;
    if (par.db1.compare(par.db2) != 0) {
        sameDatabase = false;
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str());
        tDbr->open(DBReader<unsigned int>::NOSORT);

        headerName = par.db2;
        headerName.append("_h");

        headerIndexName = par.db2;
        headerIndexName.append("_h.index");

        tempateHeaderReader = new DBReader<unsigned int>(headerName.c_str(), headerIndexName.c_str());
        tempateHeaderReader->open(DBReader<unsigned int>::NOSORT);
    }

    size_t maxSequenceLength = 0;
    unsigned int* lengths = qDbr->getSeqLens();
    for (size_t i = 0; i < qDbr->getSize(); i++) {
        maxSequenceLength = std::max((size_t) lengths[i], maxSequenceLength);
    }
    lengths = tDbr->getSeqLens();
    for (size_t i = 0; i < tDbr->getSize(); i++) {
        maxSequenceLength = std::max((size_t) lengths[i], maxSequenceLength);
    }

    DBReader<unsigned int>* clusterReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str());
    clusterReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);

    FileUtil::errorIfFileExist(par.db4.c_str());
    FileUtil::errorIfFileExist(par.db4Index.c_str());

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBWriter::BINARY_MODE);
    writer.open();

    size_t maxSetSize = findMaxSetSize(clusterReader);

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.2f);
    Debug(Debug::INFO) << "Start computing " << (!mode ? "MSAs" : "profiles") << ".\n";
#pragma omp parallel
    {
        Matcher matcher(maxSequenceLength, &subMat, tDbr->getAminoAcidDBSize(), tDbr->getSize(),
                        par.compBiasCorrection);

        MultipleAlignment aligner(maxSequenceLength, maxSetSize, &subMat, &matcher);
        PSSMCalculator calculator(&subMat, maxSequenceLength);
        MsaFilter filter(maxSequenceLength, maxSetSize, &subMat);
        Sequence  *centerSequence;
        if(par.profile == true){
            centerSequence = new Sequence(maxSequenceLength, subMat.aa2int, subMat.int2aa, Sequence::HMM_PROFILE, 0, false);
        }else{
            centerSequence = new Sequence(maxSequenceLength, subMat.aa2int, subMat.int2aa, Sequence::AMINO_ACIDS, 0, false);
        }
        Sequence **sequences = new Sequence *[maxSetSize];
        for (size_t i = 0; i < maxSetSize; i++) {
            sequences[i] = new Sequence(maxSequenceLength, subMat.aa2int, subMat.int2aa, Sequence::AMINO_ACIDS, 0,
                                        false);
        }

#pragma omp for schedule(static)
        for (size_t id = 0; id < clusterReader->getSize(); id++) {
            Log::printProgress(id);

            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif

            char *clusters = clusterReader->getData(id);
            unsigned int queryId = clusterReader->getDbKey(id);

            char *seqData = qDbr->getDataByDBKey(queryId);
            centerSequence->mapSequence(0, queryId, seqData);
            std::vector<Sequence *> seqSet;
            size_t position = 0;

            char dbKey[255 + 1];
            while (*clusters != '\0') {
                Util::parseKey(clusters, dbKey);
                unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                if (key != queryId || !sameDatabase) {
                    char *dbSeqData = tDbr->getDataByDBKey(key);
                    sequences[position]->mapSequence(0, key, dbSeqData);
                    seqSet.push_back(sequences[position]);
                    position++;
                }
                clusters = Util::skipLine(clusters);
            }

            MultipleAlignment::MSAResult res = aligner.computeMSA(centerSequence, seqSet, !par.allowDeletion);

            std::stringstream msa;
            std::string result;

            char *data;
            size_t dataSize;
            switch (mode) {
                case MSA:
                {
                    for (size_t i = 0; i < res.setSize; i++) {
                        unsigned int key;
                        char* data;
                        if(i == 0) {
                            key = centerSequence->getDbKey();
                            data = queryHeaderReader->getDataByDBKey(key);
                        } else {
                            key =  sequences[i - 1]->getDbKey();
                            data = tempateHeaderReader->getDataByDBKey(key);
                        }
                        if(par.addInternalId){
                            msa << "#" << key  << "\n";
                        }
                        msa << ">" << data;
                        for(size_t pos = 0; pos < res.msaSequenceLength; pos++){
                            char aa = res.msaSequence[i][pos];
                            msa << ((aa < MultipleAlignment::NAA) ? subMat.int2aa[(int)aa] : '-');
                        }
                        msa << "\n";
                    }
                    result = msa.str();
                    data = (char *) result.c_str();
                    dataSize = result.length();
                }
                    break;
                case PSSM:
                {
                    MsaFilter::MsaFilterResult filterRes = filter.filter(res.msaSequence, res.setSize, res.centerLength, 0, 0, -20.0, 90, 100);
                    data = (char *) calculator.computePSSMFromMSA(filterRes.setSize, res.centerLength,
                                                                  filterRes.filteredMsaSequence, false);
                    dataSize = res.centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);

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

            writer.write(data, dataSize, SSTR(queryId).c_str(), thread_idx);
        }

        for (size_t i = 0; i < maxSetSize; i++) {
            delete sequences[i];
        }
        delete[] sequences;
        delete centerSequence;
    }

    // cleanup
    writer.close();

    clusterReader->close();
    delete clusterReader;

    if (!sameDatabase) {
        tempateHeaderReader->close();
        tDbr->close();
        delete tempateHeaderReader;
        delete tDbr;
    }

    queryHeaderReader->close();
    delete queryHeaderReader;
    qDbr->close();
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
