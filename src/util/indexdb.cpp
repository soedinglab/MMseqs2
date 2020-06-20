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

std::string findIncompatibleParameter(DBReader<unsigned int>& index, const Parameters& par, const int dbtype) {
    PrefilteringIndexData meta = PrefilteringIndexReader::getMetadata(&index);
    if (meta.compBiasCorr != par.compBiasCorrection)
        return "compBiasCorrection";
    if (meta.maxSeqLength != static_cast<int>(par.maxSeqLen))
        return "maxSeqLen";
    if (meta.seqType != dbtype)
        return "seqType";
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_NUCLEOTIDES) == false && par.searchType != Parameters::SEARCH_TYPE_NUCLEOTIDES && meta.alphabetSize != par.alphabetSize.aminoacids)
        return "alphabetSize";
    if (meta.kmerSize != par.kmerSize)
        return "kmerSize";
    if (meta.mask != par.maskMode)
        return "maskMode";
    if (meta.kmerThr != par.kmerScore)
        return "kmerScore";
    if (meta.spacedKmer != par.spacedKmer)
        return "spacedKmer";
    if (par.seedScoringMatrixFile != PrefilteringIndexReader::getSubstitutionMatrixName(&index))
        return "seedScoringMatrixFile";
    if (par.spacedKmerPattern != PrefilteringIndexReader::getSpacedPattern(&index))
        return "spacedKmerPattern";
    return "";
}

int indexdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    setIndexDbDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    const bool sameDB = (par.db1 == par.db2);

    std::string path = FileUtil::getRealPathFromSymLink(par.db2dbtype);
    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> *dbr2 = NULL;
    if (sameDB == false) {
        dbr2 = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        dbr2->open(DBReader<unsigned int>::NOSORT);
    }

    const bool db1IsNucl = Parameters::isEqualDbtype(dbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    const bool db2IsNucl = dbr2 != NULL && Parameters::isEqualDbtype(dbr2->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    BaseMatrix *seedSubMat = Prefiltering::getSubstitutionMatrix(par.seedScoringMatrixFile, par.alphabetSize, 8.0f, false, (db1IsNucl && db2IsNucl));

    // memoryLimit in bytes
    size_t memoryLimit=Util::computeMemory(par.splitMemoryLimit);

    int splitMode = Parameters::TARGET_DB_SPLIT;
    par.maxResListLen = std::min(dbr.getSize(), par.maxResListLen);
    Prefiltering::setupSplit(dbr, seedSubMat->alphabetSize - 1, dbr.getDbtype(), par.threads, false, memoryLimit, 1, par.maxResListLen, par.kmerSize, par.split, splitMode);

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

    int status = EXIT_SUCCESS;
    bool recreate = true;
    std::string indexDbType = indexDB + ".dbtype";
    if (par.checkCompatible > 0 && FileUtil::fileExists(indexDbType.c_str())) {
        Debug(Debug::INFO) << "Check index " << indexDB << "\n";
        DBReader<unsigned int> index(indexDB.c_str(), (indexDB + ".index").c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        index.open(DBReader<unsigned int>::NOSORT);

        if (Parameters::isEqualDbtype(dbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES) && par.searchType == Parameters::SEARCH_TYPE_NUCLEOTIDES && par.PARAM_ALPH_SIZE.wasSet) {
            Debug(Debug::WARNING) << "Alphabet size is not taken into account for compatibility check in nucleotide search.\n";
        }

        std::string check;
        const bool compatible = PrefilteringIndexReader::checkIfIndexFile(&index) && (check = findIncompatibleParameter(index, par, dbr.getDbtype())) == "";
        index.close();
        if (compatible) {
            Debug(Debug::INFO) << "Index is up to date and compatible. Force recreation with --check-compatibility 0 parameter.\n";
            recreate = false;
        } else {
            if (par.checkCompatible == 2) {
                Debug(Debug::ERROR) << "Index is incompatible. Incompatible parameter: " << check << "\n";
                recreate = false;
                status = EXIT_FAILURE;
            } else {
                Debug(Debug::WARNING) << "Index is incompatible and will be recreated. Incompatible parameter: " << check << "\n";
                recreate = true;
            }
        }
    }

    if (recreate) {
        DBReader<unsigned int> hdbr1(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
        hdbr1.open(DBReader<unsigned int>::NOSORT);

        DBReader<unsigned int> *hdbr2 = NULL;
        if (sameDB == false) {
            hdbr2 = new DBReader<unsigned int>(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
            hdbr2->open(DBReader<unsigned int>::NOSORT);
        }

        DBReader<unsigned int>::removeDb(indexDB);
        PrefilteringIndexReader::createIndexFile(indexDB, &dbr, dbr2, &hdbr1, hdbr2, seedSubMat, par.maxSeqLen,
                                                 par.spacedKmer, par.spacedKmerPattern, par.compBiasCorrection,
                                                 seedSubMat->alphabetSize, par.kmerSize, par.maskMode, par.maskLowerCaseMode,
                                                 par.kmerScore, par.split);

        if (hdbr2 != NULL) {
            hdbr2->close();
            delete hdbr2;
        }

        hdbr1.close();
    }

    if (dbr2 != NULL) {
        dbr2->close();
        delete dbr2;
    }

    delete seedSubMat;
    dbr.close();

    return status;
}
