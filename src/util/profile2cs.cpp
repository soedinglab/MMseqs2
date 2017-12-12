#include <ProfileStates.h>
#include <PSSMCalculator.h>
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

int profile2cs(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2,  true, 0, MMseqsParameter::COMMAND_PROFILE);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> profileReader(par.db1.c_str(), par.db1Index.c_str());
    profileReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    writer.open();

    size_t entries = profileReader.getSize();

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, 0.0);
    ProfileStates ps(subMat.pBack);

    Debug(Debug::INFO) << "Start converting profiles.\n";
#pragma omp parallel
    {
        Sequence seq(par.maxSeqLen, Sequence::HMM_PROFILE, &subMat, 0, false, false);

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char* buffer = new char[64];
        std::string result;
        result.reserve(par.maxSeqLen);
#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < entries; ++i) {
            Debug::printProgress(i);
            result.clear();

            unsigned int key = profileReader.getDbKey(i);
            seq.mapSequence(i, key, profileReader.getData(i));
            ps.discretize(seq.getProfile(), seq.L, result);
            result.push_back('\0'); // needed to avoid seq. len calculation problems (aa sequence has \n\0, PS \0\0)
            for(size_t i = 0; i < result.size() - 1; i++){ // do not overwrite last \0
//                std::cout << (int)result[i] << std::endl;
                char val = result[i];
                val += 1; // avoid null byte (needed for read in)
                result[i] = val;
            }
            writer.writeData(result.c_str(), result.length(), key, thread_idx);
        }
        delete[] buffer;
    }
    writer.close(DBReader<unsigned int>::DBTYPE_PROFILE_STATE);
    Debug(Debug::INFO) << "\nDone.\n";
    profileReader.close();
    return EXIT_SUCCESS;
}
