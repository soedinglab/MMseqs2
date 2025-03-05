#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Masker.h"

#ifdef OPENMP
#include <omp.h>
#endif

int masksequence(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::NOSORT);

    BaseMatrix *subMat;
    if (Parameters::isEqualDbtype(reader.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    } else {
        // keep score bias at 0.0 (improved ROC)
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    }

    // need to prune low scoring k-mers through masking

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();
#pragma omp parallel
    {

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        Masker masker(*subMat);
        unsigned char *charSequence = new unsigned char[reader.getMaxSeqLen() + 1];
        Sequence seq(reader.getMaxSeqLen(), reader.getDbtype(), subMat,  0, false, false);
#pragma omp for schedule(dynamic, 1)
        for (size_t id = 0; id < reader.getSize(); ++id) {
            seq.mapSequence(id, reader.getDbKey(id), reader.getData(id, thread_idx), reader.getSeqLen(id));
            masker.maskSequence(seq, par.maskMode, par.maskProb, par.maskLowerCaseMode, par.maskNrepeats);
            memcpy(charSequence, seq.getSeqData(), seq.L * sizeof(char));
            masker.applySoftmasking(charSequence, seq.numSequence, seq.L);
            charSequence[seq.L] = '\n';
            writer.writeData((const char *)charSequence, seq.L + 1,  seq.getDbKey(), thread_idx);
        }
        delete[] charSequence;
    }
    writer.close(true);
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_ANCILLARY);
    reader.close();

    delete subMat;
    return EXIT_SUCCESS;
}
