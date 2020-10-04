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

int profile2seq(int argc, const char **argv, const Command &command, bool consensus) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_PROFILE);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
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
        Sequence seq(par.maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, &subMat, 0, false, false, false);
        std::string result;
        result.reserve(par.maxSeqLen);
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < entries; ++i) {
            progress.updateProgress();
            seq.mapProfile(reader.getData(i, thread_idx), reader.getSeqLen(i));
            unsigned char * sequence = consensus ? seq.numConsensusSequence : seq.numSequence;
            for (int aa = 0; aa < seq.L; aa++) {
                result.append(1, subMat.num2aa[sequence[aa]]);
            }
            result.append(1, '\n');
            writer.writeData(result.c_str(), result.length(), reader.getDbKey(i), thread_idx);
            result.clear();
        }
    }
    writer.close(true);
    reader.close();
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_ANCILLARY);

    return EXIT_SUCCESS;
}

int profile2consensus(int argc, const char **argv, const Command &command) {
    return profile2seq(argc, argv, command, true);
}

int profile2repseq(int argc, const char **argv, const Command &command) {
    return profile2seq(argc, argv, command, false);
}
