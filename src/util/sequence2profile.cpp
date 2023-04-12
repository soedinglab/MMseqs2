// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte


#include "CSProfile.h"
#include "MathUtil.h"
#include "DBReader.h"
#include "Parameters.h"
#include "DBWriter.h"

#include <string>
#include <PSSMMasker.h>


#ifdef OPENMP
#include <omp.h>
#endif

int sequence2profile(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);

    DBReader<unsigned int> sequenceDb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    sequenceDb.open(DBReader<unsigned int>::NOSORT);

    int type = Parameters::DBTYPE_HMM_PROFILE;
    if (par.pcmode == Parameters::PCMODE_CONTEXT_SPECIFIC) {
        type = DBReader<unsigned int>::setExtendedDbtype(type, Parameters::DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS);
    }
    DBWriter resultDbw(par.db2.c_str(), par.db2Index.c_str(), par.threads,  par.compressed, type);
    resultDbw.open();
    Debug::Progress progress( sequenceDb.getSize());

#pragma omp parallel
    {
        Sequence seq(par.maxSeqLen, sequenceDb.getDbtype(), &subMat, 0, false, false);
        CSProfile ps(par.maxSeqLen);
        ProbabilityMatrix probMatrix(subMat);
        PSSMMasker masker(sequenceDb.getMaxSeqLen(), probMatrix, subMat);
        char * pssm = (char * )mem_align(16, Sequence::PROFILE_AA_SIZE * sequenceDb.getMaxSeqLen() * sizeof(char));
        float * Neff_M = new float[sequenceDb.getMaxSeqLen()];
        std::fill(Neff_M, Neff_M + sequenceDb.getMaxSeqLen(), 1.0f);

        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string result;
        result.reserve(sequenceDb.getMaxSeqLen() * Sequence::PROFILE_READIN_SIZE);
#pragma omp for schedule(static)
        for (size_t id = 0; id < sequenceDb.getSize(); id++) {
            progress.updateProgress();
            char *seqData     = sequenceDb.getData(id, thread_idx);
            unsigned int queryKey = sequenceDb.getDbKey(id);
            unsigned int seqLen = sequenceDb.getSeqLen(id);

            seq.mapSequence(id, queryKey, seqData, seqLen);
            float * profile = ps.computeSequenceCs(seq.numSequence, seq.L, par.tau);
            PSSMCalculator::computeLogPSSM(&subMat, pssm, profile, 8.0,  seq.L, 0.0);
#ifdef GAP_POS_SCORING
            PSSMCalculator::Profile pssmRes(pssm, profile, Neff_M, NULL, NULL, seq.numSequence);
#else
            PSSMCalculator::Profile pssmRes(pssm, profile, Neff_M, seq.numSequence);
#endif
            if (par.maskProfile == true) {
                masker.mask(seq, par.maskProb, pssmRes);
            }
            pssmRes.toBuffer(seq, subMat, result);

            resultDbw.writeData(result.c_str(), result.size(), queryKey, thread_idx);
            result.clear();
        }
        free(pssm);
        delete [] Neff_M;
    }
    sequenceDb.close();
    resultDbw.close();
    return 0;
}
