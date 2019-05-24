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

    DBReader<unsigned int> profileReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    profileReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, shouldCompress, dbType);
    writer.open();

    const bool isDbOutput = par.dbOut;

    size_t entries = profileReader.getSize();
    Debug::Progress progress(entries);

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, 0.0);
    Debug(Debug::INFO) << "Start converting profiles.\n";
#pragma omp parallel
    {
        Sequence seq(par.maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, &subMat, 0, false, par.compBiasCorrection, false);

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        char* buffer = new char[64];
        std::string result;
        result.reserve(10 * 1024 * 1024);

#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < entries; ++i) {
            progress.updateProgress();
            result.clear();

            unsigned int key = profileReader.getDbKey(i);
            seq.mapSequence(i, key, profileReader.getData(i, thread_idx));

            if (isDbOutput == false) {
                result.append("Query profile of sequence ");
                Itoa::u32toa_sse2(key, buffer);
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
                    //probas: result.append(std::to_string(seq.profile[j * Sequence::PROFILE_AA_SIZE + aa]));

                }
                result.push_back('\n');
            }
            writer.writeData(result.c_str(), result.length(), key, thread_idx, isDbOutput);
        }
        delete[] buffer;
    }
    writer.close(isDbOutput == false);
    if (isDbOutput == false) {
        remove(par.db2Index.c_str());
    }



    profileReader.close();

    return EXIT_SUCCESS;
}
