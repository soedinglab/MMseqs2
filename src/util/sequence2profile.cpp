// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte


#include "CSProfile.h"
#include "MathUtil.h"
#include "DBReader.h"
#include "DBWriter.h"

#include <string>


#ifdef OPENMP
#include <omp.h>
#endif

int sequence2profile(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 8.0, 0.0);

    DBReader<unsigned int> sequenceDb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    sequenceDb.open(DBReader<unsigned int>::NOSORT);
    DBWriter resultDbw(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_HMM_PROFILE);
    resultDbw.open();
#pragma omp parallel
    {
        Sequence seq(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, false);
        CSProfile ps(par.maxSeqLen);
        char * data = new char[par.maxSeqLen * Sequence::PROFILE_READIN_SIZE];
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < sequenceDb.getSize(); id++) {
            Debug::printProgress(id);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            char *seqData     = sequenceDb.getData(id, thread_idx);
            unsigned int queryKey = sequenceDb.getDbKey(id);
            seq.mapSequence(id, queryKey, seqData);
            float * prob =  ps.computeProfile(&seq, par.neff, par.tau);

            size_t idx = 0;
            for (int i = 0; i < seq.L; i++) {
                for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE;aa++)
                {
                    data[idx++] = Sequence::scoreMask(prob[i*Sequence::PROFILE_READIN_SIZE + aa]);
                    //std::cout<<"\t"<<(int)data[idx-1];
                }
                //std::cout<<std::endl;
                data[idx++] = static_cast<unsigned char>(seq.int_sequence[i]); // query
                data[idx++] = static_cast<unsigned char>(seq.int_sequence[i]); // consensus
                data[idx++] = MathUtil::convertNeffToChar(prob[i*Sequence::PROFILE_READIN_SIZE + Sequence::PROFILE_AA_SIZE]);
            }
            resultDbw.writeData(data, idx, queryKey, thread_idx);
        }
        delete [] data;
    }
    sequenceDb.close();
    resultDbw.close();
    return 0;
}
