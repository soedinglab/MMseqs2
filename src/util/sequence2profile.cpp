// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte


#include "CSProfile.h"
#include "MathUtil.h"
#include "DBReader.h"
#include "Parameters.h"
#include "DBWriter.h"

#include <string>


#ifdef OPENMP
#include <omp.h>
#endif

int sequence2profile(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 8.0, 0.0);

    DBReader<unsigned int> sequenceDb(par.db1.c_str(), par.db1Index.c_str());
    sequenceDb.open(DBReader<unsigned int>::NOSORT);
    DBWriter resultDbw(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    resultDbw.open();
#pragma omp parallel
    {
        Sequence seq(par.maxSeqLen, Sequence::AMINO_ACIDS, &subMat, 0, false, false);
        CSProfile ps(par.maxSeqLen);
        char * data = new char[par.maxSeqLen * Sequence::PROFILE_AA_SIZE];
#pragma omp for schedule(static)
        for (size_t id = 0; id < sequenceDb.getSize(); id++) {
            Debug::printProgress(id);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char *seqData     = sequenceDb.getData(id);
            unsigned int queryKey = sequenceDb.getDbKey(id);
            seq.mapSequence(id, queryKey, seqData);
            float * prob =  ps.computeProfile(&seq, par.neff, par.tau);
            // a bit hacky. I reuse the pssm data as storage
            size_t dataSize = seq.L * Sequence::PROFILE_READIN_SIZE * sizeof(char);
            for (size_t i = 0; i < dataSize; i++) {
                // Avoid a null byte result
                if((i+1) % Sequence::PROFILE_READIN_SIZE == 0){
                    data[i] =  MathUtil::convertNeffToChar(prob[i]); // avoid \0 by setting 0 values to MINIFLOAT_MIN
                }else{
                    // Set to MINIFLOAT_MAX to have a 0 value. MAX is 1.96875 (probs range is 0.0-1.0)
                    // We map this value to 0 in Sequence::mapProfile/mapProfileState
                    char val = MathUtil::convertFloatToChar(prob[i]);
                    data[i] = (val == 0) ? MINIFLOAT_MAX : val; // avoid \0 by setting 0 values to MINIFLOAT_MIN
                }
            }
            resultDbw.writeData(data, dataSize, queryKey, thread_idx);
        }
        delete [] data;
    }
    sequenceDb.close();
    resultDbw.close();
    return 0;
}
