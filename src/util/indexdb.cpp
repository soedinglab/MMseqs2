#include "DBReader.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "Prefiltering.h"
#include "Parameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

void setCreateIndexDefaults(Parameters *p) {
    p->sensitivity = 5;
}

int indexdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setCreateIndexDefaults(&par);
    par.parseParameters(argc, argv, command, 2);

    if (par.split > 1) {
        Debug(Debug::ERROR) << "Creating a split index is not supported anymore.\n";
        Debug(Debug::ERROR) << "Please use the prefilter without a precomputed index if you do not have enough memory.";
        EXIT(EXIT_FAILURE);
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str());
    dbr.open(DBReader<unsigned int>::NOSORT);
    BaseMatrix *subMat = Prefiltering::getSubstitutionMatrix(par.scoringMatrixFile, par.alphabetSize, 8.0f, false, false);

    int kmerSize = par.kmerSize;
    int split = 1;
    int splitMode = Parameters::TARGET_DB_SPLIT;

    size_t memoryLimit;
    if (par.splitMemoryLimit > 0) {
        memoryLimit = static_cast<size_t>(par.splitMemoryLimit) * 1024;
    } else {
        memoryLimit = static_cast<size_t>(Util::getTotalSystemMemory() * 0.9);
    }
    Prefiltering::setupSplit(dbr, subMat->alphabetSize, par.threads, false, par.maxResListLen, memoryLimit,
                             &kmerSize, &split, &splitMode);

    bool kScoreSet = false;
    for (size_t i = 0; i < par.indexdb.size(); i++) {
        if (par.indexdb[i].uniqid == par.PARAM_K_SCORE.uniqid && par.indexdb[i].wasSet) {
            kScoreSet = true;
        }
    }

    if (dbr.getDbtype() != DBReader<unsigned int>::DBTYPE_PROFILE && kScoreSet == false) {
        par.kmerScore = 0;
    }

    // query seq type is actually unknown here, but if we pass HMM_PROFILE then its +20 k-score
    int kmerThr = Prefiltering::getKmerThreshold(par.sensitivity, Sequence::AMINO_ACIDS, par.kmerScore, kmerSize);

    DBReader<unsigned int> *hdbr = NULL;
    if (par.includeHeader == true) {
        std::string hdr_filename(par.db1);
        hdr_filename.append("_h");

        std::string hdr_index_filename(par.db1);
        hdr_index_filename.append("_h.index");

        hdbr = new DBReader<unsigned int>(hdr_filename.c_str(), hdr_index_filename.c_str());
        hdbr->open(DBReader<unsigned int>::NOSORT);
    }

    PrefilteringIndexReader::createIndexFile(par.db2, &dbr, hdbr, subMat, par.maxSeqLen,
                                             par.spacedKmer, par.compBiasCorrection,
                                             subMat->alphabetSize, kmerSize, par.diagonalScoring,
                                             par.maskMode, dbr.getDbtype(), kmerThr, par.threads);

    if (hdbr != NULL) {
        hdbr->close();
        delete hdbr;
    }

    delete subMat;
    dbr.close();

    return EXIT_SUCCESS;
}
