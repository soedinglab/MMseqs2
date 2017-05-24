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

int createindex(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setCreateIndexDefaults(&par);
    par.parseParameters(argc, argv, command, 1);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str());
    dbr.open(DBReader<unsigned int>::NOSORT);

    BaseMatrix *subMat = Prefiltering::getSubstitutionMatrix(par.scoringMatrixFile, par.alphabetSize, 8.0f, false);

    int kmerSize = par.kmerSize;
    int split = 1;
    int splitMode = Parameters::TARGET_DB_SPLIT;

    Prefiltering::setupSplit(dbr, subMat->alphabetSize, par.threads, false, par.maxResListLen, &kmerSize, &split, &splitMode);

    bool kScoreSet = false;
    for (size_t i = 0; i < par.createindex.size(); i++) {
        if (par.createindex[i].uniqid == par.PARAM_K_SCORE.uniqid && par.createindex[i].wasSet) {
            kScoreSet = true;
        }
    }

    if (par.targetSeqType != Sequence::HMM_PROFILE && kScoreSet == false) {
        par.kmerScore = 0;
    }

    int kmerThr = Prefiltering::getKmerThreshold(par.sensitivity, par.querySeqType, par.kmerScore, kmerSize);

    DBReader<unsigned int> *hdbr = NULL;
    if (par.includeHeader == true) {
        std::string hdr_filename(par.db1);
        hdr_filename.append("_h");

        std::string hdr_index_filename(par.db1);
        hdr_index_filename.append("_h.index");

        hdbr = new DBReader<unsigned int>(hdr_filename.c_str(), hdr_index_filename.c_str());
        hdbr->open(DBReader<unsigned int>::NOSORT);
    }

    PrefilteringIndexReader::createIndexFile(par.db1, &dbr, hdbr, subMat, par.maxSeqLen,
                                             par.spacedKmer, par.compBiasCorrection,
                                             subMat->alphabetSize, kmerSize, par.diagonalScoring,
                                             par.maskMode, par.targetSeqType, kmerThr, par.threads);

    if (hdbr != NULL) {
        hdbr->close();
        delete hdbr;
    }

    delete subMat;
    dbr.close();

    return EXIT_SUCCESS;
}
