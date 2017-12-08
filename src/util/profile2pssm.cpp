#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "itoa.h"

#ifdef OPENMP
#include <omp.h>
#endif

int profile2pssm(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2,  true, 0, MMseqsParameter::COMMAND_PROFILE);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> profileReader(par.db1.c_str(), par.db1Index.c_str());
    profileReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    writer.open();

    const bool isDbOutput = par.dbOut;

    size_t entries = profileReader.getSize();

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, 0.0);
    Debug(Debug::INFO) << "Start converting profiles.\n";
#pragma omp parallel
    {
        Sequence seq(par.maxSeqLen, Sequence::HMM_PROFILE, &subMat, 0, false, par.compBiasCorrection);

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        char* buffer = new char[64];
        std::string result;
        result.reserve(10 * 1024 * 1024);

#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < entries; ++i) {
            Debug::printProgress(i);
            result.clear();

            unsigned int key = profileReader.getDbKey(i);
            seq.mapSequence(i, key, profileReader.getData(i));

            if (isDbOutput == false) {
                result.append("Query profile of sequence ");
                Itoa::i32toa_sse2(key, buffer);
                result.append(buffer);
                result.push_back('\n');
            }

            result.append("Pos\tCns");
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                result.push_back('\t');
                result.push_back(subMat.int2aa[aa]);
            }
            result.push_back('\n');

            for (int j = 0; j < seq.L; ++j) {
                Itoa::i32toa_sse2(j, buffer);
                result.append(buffer);
                result.push_back('\t');
                result.push_back(subMat.int2aa[seq.int_consensus_sequence[j]]);
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                    result.push_back('\t');
                    Itoa::i32toa_sse2(seq.profile_for_alignment[aa * seq.L + j], buffer);
                    result.append(buffer);
                }
                result.push_back('\n');
            }
            writer.writeData(result.c_str(), result.length(), key, thread_idx, isDbOutput);
        }
        delete[] buffer;
    }
    writer.close();
    if (isDbOutput == false) {
        remove(par.db2Index.c_str());
    }

    Debug(Debug::INFO) << "\nDone.\n";

    profileReader.close();

    return EXIT_SUCCESS;
}
