#include "DBReader.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "Prefiltering.h"
#include "Parameters.h"

int createindex(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 1);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str());
    dbr.open(DBReader<unsigned int>::NOSORT);

    BaseMatrix *subMat = Prefiltering::getSubstitutionMatrix(par.scoringMatrixFile, par.alphabetSize, 8.0f, false);

    int kmerSize = par.kmerSize;
    int split = par.split;
    int splitMode = Parameters::TARGET_DB_SPLIT;
    Prefiltering::setupSplit(dbr, subMat->alphabetSize, par.threads, false, &kmerSize, &split, &splitMode);

    PrefilteringIndexReader::createIndexFile(par.db1, &dbr, subMat, par.maxSeqLen,
                                             par.spacedKmer, par.compBiasCorrection, split,
                                             subMat->alphabetSize, kmerSize, par.diagonalScoring,
                                             par.maskMode, par.threads);

    delete subMat;
    dbr.close();

    return EXIT_SUCCESS;
}
