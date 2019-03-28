#include "DBReader.h"
#include "Util.h"
#include "FileUtil.h"
#include "PrefilteringIndexReader.h"
#include "Prefiltering.h"
#include "Parameters.h"

#ifdef OPENMP
#include <omp.h>
#endif

void setIndexDbDefaults(Parameters *p) {
    p->sensitivity = 5.7;
}

bool isIndexCompatible(DBReader<unsigned int>& index, const Parameters& par, const int dbtype) {
    PrefilteringIndexData meta = PrefilteringIndexReader::getMetadata(&index);
    if (meta.compBiasCorr != par.compBiasCorrection)
        return false;
    if (meta.maxSeqLength != static_cast<int>(par.maxSeqLen))
        return false;
    if (meta.seqType != dbtype)
        return false;
    if (meta.alphabetSize != par.alphabetSize)
        return false;
    if (meta.kmerSize != par.kmerSize)
        return false;
    if (meta.mask != par.maskMode)
        return false;
    if (meta.kmerThr != par.kmerScore)
        return false;
    if (meta.spacedKmer != par.spacedKmer)
        return false;
    if (par.seedScoringMatrixFile != PrefilteringIndexReader::getSubstitutionMatrixName(&index))
        return false;
    if (par.spacedKmerPattern != PrefilteringIndexReader::getSpacedPattern(&index))
        return false;
    if (meta.headers2 == 1 && (par.db1 != par.db2))
        return true;
    return true;
}

int indexdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setIndexDbDefaults(&par);
    par.parseParameters(argc, argv, command, 2);

    if (par.split > 1) {
        Debug(Debug::ERROR) << "Creating a split index is not supported currently.\n";
        Debug(Debug::ERROR) << "Please use the prefilter without a precomputed index.\n";
        EXIT(EXIT_FAILURE);
    }

    const bool sameDB = (par.db1 == par.db2);
    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *dbr2 = NULL;

    if (sameDB == false) {
        dbr2 = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        dbr2->open(DBReader<unsigned int>::NOSORT);
    }

    bool db2IsNucl = dbr2 != NULL && Parameters::isEqualDbtype(dbr2->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    bool db1IsNucl = Parameters::isEqualDbtype(dbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);

    BaseMatrix *seedSubMat = Prefiltering::getSubstitutionMatrix(par.seedScoringMatrixFile, par.alphabetSize, 8.0f, false, (db1IsNucl && db2IsNucl));

    int split = 1;
    int splitMode = Parameters::TARGET_DB_SPLIT;

    size_t memoryLimit;
    if (par.splitMemoryLimit > 0) {
        memoryLimit = static_cast<size_t>(par.splitMemoryLimit) * 1024;
    } else {
        memoryLimit = static_cast<size_t>(Util::getTotalSystemMemory() * 0.9);
    }
    Prefiltering::setupSplit(dbr, seedSubMat->alphabetSize - 1, dbr.getDbtype(), par.threads, false, par.maxResListLen, memoryLimit, &par.kmerSize, &split, &splitMode);

    bool kScoreSet = false;
    for (size_t i = 0; i < par.indexdb.size(); i++) {
        if (par.indexdb[i]->uniqid == par.PARAM_K_SCORE.uniqid && par.indexdb[i]->wasSet) {
            kScoreSet = true;
        }
    }

    const bool isProfileSearch = (Parameters::isEqualDbtype(dbr.getDbtype(), Parameters::DBTYPE_HMM_PROFILE));
    if (isProfileSearch && kScoreSet == false) {
        par.kmerScore = 0;
    }

    // TODO: investigate if it makes sense to mask the profile consensus sequence
    if (isProfileSearch) {
        par.maskMode = 0;
    }

    // query seq type is actually unknown here, but if we pass DBTYPE_HMM_PROFILE then its +20 k-score
    par.kmerScore = Prefiltering::getKmerThreshold(par.sensitivity, isProfileSearch, par.kmerScore, par.kmerSize);

    std::string indexDB = PrefilteringIndexReader::indexName(par.db2);
    if (par.checkCompatible && FileUtil::fileExists(indexDB.c_str())) {
        Debug(Debug::INFO) << "Check index " << indexDB << "\n";
        DBReader<unsigned int> index(indexDB.c_str(), (indexDB + ".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        index.open(DBReader<unsigned int>::NOSORT);
        if (PrefilteringIndexReader::checkIfIndexFile(&index) && isIndexCompatible(index, par, dbr.getDbtype())) {
            Debug(Debug::INFO) << "Index is already up to date and compatible. Force recreation with --check-compatibility 0 parameter.\n";
            return EXIT_SUCCESS;
        } else {
            Debug(Debug::WARNING) << "Index is incompatible and will be recreated.\n";
        }
    }

    DBReader<unsigned int> *hdbr1 = NULL;
    DBReader<unsigned int> *hdbr2 = NULL;
    {
        hdbr1 = new DBReader<unsigned int>(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        hdbr1->open(DBReader<unsigned int>::NOSORT);
        if (sameDB == false) {
            hdbr2 = new DBReader<unsigned int>(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
            hdbr2->open(DBReader<unsigned int>::NOSORT);
        }
    }

    PrefilteringIndexReader::createIndexFile(indexDB, &dbr, dbr2, hdbr1, hdbr2, seedSubMat, par.maxSeqLen,
                                             par.spacedKmer, par.spacedKmerPattern, par.compBiasCorrection,
                                             seedSubMat->alphabetSize, par.kmerSize, par.maskMode, par.maskLowerCaseMode,
                                             par.kmerScore);

    if (hdbr2 != NULL) {
        hdbr2->close();
        delete hdbr2;
    }

    if (hdbr1 != NULL) {
        hdbr1->close();
        delete hdbr1;
    }

    if (dbr2 != NULL) {
        dbr2->close();
        delete dbr2;
    }

    delete seedSubMat;
    dbr.close();

    return EXIT_SUCCESS;
}
