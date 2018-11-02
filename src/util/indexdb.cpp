#include "DBReader.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "Prefiltering.h"
#include "Parameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

void setCreateIndexDefaults(Parameters *p) {
    p->sensitivity = 5.7;
}

int indexdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setCreateIndexDefaults(&par);
    par.overrideParameterDescription((Command &) command, par.PARAM_MASK_RESIDUES.uniqid, "0: w/o low complexity masking, 1: with low complexity masking, 2: add both masked and unmasked sequences to index", "^[0-2]{1}", par.PARAM_MASK_RESIDUES.category);
    par.parseParameters(argc, argv, command, 2);

    if (par.split > 1) {
        Debug(Debug::ERROR) << "Creating a split index is not supported anymore.\n";
        Debug(Debug::ERROR) << "Please use the prefilter without a precomputed index if you do not have enough memory.";
        EXIT(EXIT_FAILURE);
    }

    const bool sameDB = (par.db1 == par.db2);
    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str());
    dbr.open(DBReader<unsigned int>::NOSORT);
    BaseMatrix *subMat = Prefiltering::getSubstitutionMatrix(par.scoringMatrixFile, par.alphabetSize, 8.0f, false);

    int kmerSize = par.kmerSize;
    int split = 1;
    int splitMode = Parameters::TARGET_DB_SPLIT;

    size_t memoryLimit;
    if (par.splitMemoryLimit > 0) {
        memoryLimit = static_cast<size_t>(par.splitMemoryLimit) * 1024;
    } else {
        memoryLimit = static_cast<size_t>(Util::getTotalSystemMemory() * 0.9);
    }
    Prefiltering::setupSplit(dbr, subMat->alphabetSize, dbr.getDbtype(), par.threads, false, par.maxResListLen, memoryLimit, &kmerSize, &split, &splitMode);

    bool kScoreSet = false;
    for (size_t i = 0; i < par.indexdb.size(); i++) {
        if (par.indexdb[i].uniqid == par.PARAM_K_SCORE.uniqid && par.indexdb[i].wasSet) {
            kScoreSet = true;
        }
    }

    const bool isProfileSearch = dbr.getDbtype() == Sequence::HMM_PROFILE;
    if (isProfileSearch && kScoreSet == false) {
        par.kmerScore = 0;
    }

    // TODO: investigate if it makes sense to mask the profile consensus sequence
    if (isProfileSearch) {
        par.maskMode = 0;
    }

    // query seq type is actually unknown here, but if we pass HMM_PROFILE then its +20 k-score
    const int kmerThr = Prefiltering::getKmerThreshold(par.sensitivity, isProfileSearch, par.kmerScore, kmerSize);

    DBReader<unsigned int> *hdbr1 = NULL;
    DBReader<unsigned int> *hdbr2 = NULL;
    if (par.includeHeader == true) {
        hdbr1 = new DBReader<unsigned int>(par.hdr1.c_str(), par.hdr1Index.c_str());
        hdbr1->open(DBReader<unsigned int>::NOSORT);
        if (sameDB == false) {
            hdbr2 = new DBReader<unsigned int>(par.hdr2.c_str(), par.hdr2Index.c_str());
            hdbr2->open(DBReader<unsigned int>::NOSORT);
        }
    }

    PrefilteringIndexReader::createIndexFile(par.db2, &dbr, hdbr1, hdbr2, subMat, par.maxSeqLen,
                                             par.spacedKmer, par.compBiasCorrection, subMat->alphabetSize,
                                             kmerSize, par.maskMode, kmerThr);

    if (hdbr2 != NULL) {
        hdbr2->close();
        delete hdbr2;
    }

    if (hdbr1 != NULL) {
        hdbr1->close();
        delete hdbr1;
    }

    delete subMat;
    dbr.close();

    return EXIT_SUCCESS;
}
