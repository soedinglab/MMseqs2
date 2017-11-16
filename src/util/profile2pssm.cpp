#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"

#ifdef OPENMP
#include <omp.h>
#endif

int profile2pssm(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2,  true, false, MMseqsParameter::COMMAND_PROFILE);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> profileReader(par.db1.c_str(), par.db1Index.c_str());
    profileReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> consReader((par.db1 + "_consensus").c_str(),
                                      (par.db1 + "_consensus.index").c_str());
    consReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    writer.open();

    const bool isDbOutput = par.dbOut;

    size_t entries = profileReader.getSize();

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, 0.0);
#define APPEND(x) do { std::string __x__ = (x); writer.writeAdd(__x__.c_str(), __x__.length(), thread_idx); } while (false);
    Debug(Debug::INFO) << "Start converting profiles.\n";
#pragma omp parallel
    {
        Sequence seq(par.maxSeqLen, subMat.aa2int, subMat.int2aa, Sequence::HMM_PROFILE, 0, false, par.compBiasCorrection);

#pragma omp for schedule(dynamic, 100)
        for (size_t i = 0; i < entries; ++i) {
            Debug::printProgress(i);

            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

            unsigned int key = profileReader.getDbKey(i);
            char *data = profileReader.getData(i);
            seq.mapProfile(data);


            writer.writeStart(thread_idx);
            if (isDbOutput == false) {
                APPEND("Query profile of sequence ");
                APPEND(SSTR(key));
                APPEND("\n");
            }
            APPEND("Pos\tCns\t");
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                APPEND(SSTR(subMat.int2aa[aa]));
                if (aa < Sequence::PROFILE_AA_SIZE - 1) {
                    APPEND("\t");
                } else {
                    APPEND("\n");
                }
            }

            for (int j = 0; j < seq.L; j++) {
                APPEND(SSTR(j));
                APPEND("\t");
                APPEND(SSTR(subMat.int2aa[seq.int_sequence[j]]));
                APPEND("\t");
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                    APPEND(SSTR((int)seq.profile_for_alignment[aa * seq.L + j]));
                    if (aa < Sequence::PROFILE_AA_SIZE - 1) {
                        APPEND("\t");
                    } else {
                        APPEND("\n");
                    }
                }
            }

            writer.writeEnd(key, thread_idx, isDbOutput);
        }
    }
#undef APPEND
    writer.close();

    if (isDbOutput == false) {
        remove(par.db2Index.c_str());
    }

    Debug(Debug::INFO) << "\nDone.\n";

    consReader.close();
    profileReader.close();

    return EXIT_SUCCESS;
}

