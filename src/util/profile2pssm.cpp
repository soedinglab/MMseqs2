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
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_PROFILE);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    const bool isDbOutput = par.dbOut;
    const bool shouldCompress = isDbOutput == true && par.compressed == true;
    const int dbType = isDbOutput == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, shouldCompress, dbType);
    writer.open();

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0f, 0.0);

    size_t entries = reader.getSize();
    Debug::Progress progress(entries);
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

        Sequence seq(par.maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, &subMat, 0, false, par.compBiasCorrection, false);

        char buffer[256];
        std::string result;
        result.reserve(10 * 1024 * 1024);

#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < entries; ++i) {
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            seq.mapSequence(i, key, reader.getData(i, thread_idx), reader.getSeqLen(i));


            if (isDbOutput == false) {
                result.append("Query profile of sequence ");
                Itoa::u32toa_sse2(key, buffer);
                result.append(buffer);
                result.push_back('\n');
            }

            result.append("Pos\tCns");
            for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                result.push_back('\t');
                result.push_back(subMat.num2aa[aa]);
            }
            result.push_back('\n');

            for (int j = 0; j < seq.L; ++j) {
                Itoa::i32toa_sse2(j, buffer);
                result.append(buffer);
                result.push_back('\t');
                result.push_back(subMat.num2aa[seq.numConsensusSequence[j]]);
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
                    result.push_back('\t');
                    Itoa::i32toa_sse2(seq.profile_for_alignment[aa * seq.L + j], buffer);
                    result.append(buffer);
                    //probas: result.append(std::to_string(seq.profile[j * Sequence::PROFILE_AA_SIZE + aa]));

                }
                result.push_back('\n');
            }
            writer.writeData(result.c_str(), result.length(), key, thread_idx, isDbOutput);
            result.clear();
        }
    }
    writer.close(isDbOutput == false);
    if (isDbOutput == false) {
        remove(par.db2Index.c_str());
    }
    reader.close();

    return EXIT_SUCCESS;
}