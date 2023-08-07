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
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_NUCLEOTIDES) == false && par.searchType != Parameters::SEARCH_TYPE_NUCLEOTIDES && meta.alphabetSize != par.alphabetSize.values.aminoacid())
        return "alphabetSize";
    if (meta.kmerSize != par.kmerSize)
        return "kmerSize";
    if (meta.mask != par.maskMode)
        return "maskMode";
    if (Parameters::isEqualDbtype(dbtype, Parameters::DBTYPE_HMM_PROFILE)) {
        if (meta.kmerThr != par.kmerScore.values.profile())
            return "kmerScore";
    } else {
        if (meta.kmerThr != par.kmerScore.values.sequence())
            return "kmerScore";
    }
    if (meta.spacedKmer != par.spacedKmer)
        return "spacedKmer";
    if (BaseMatrix::unserializeName(par.seedScoringMatrixFile.values.aminoacid().c_str()) != PrefilteringIndexReader::getSubstitutionMatrixName(&index) &&
        BaseMatrix::unserializeName(par.seedScoringMatrixFile.values.nucleotide().c_str()) != PrefilteringIndexReader::getSubstitutionMatrixName(&index))
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
    std::string alnDbtypeFile = par.db1 + "_aln.dbtype";
    std::string alnFile = par.db1 + "_aln";
    std::string alnIndexFile = par.db1 + "_aln.index";
    if(FileUtil::fileExists(alnDbtypeFile.c_str()) == false){
        alnDbtypeFile =  par.db1 + "_clu.dbtype";
        alnFile = par.db1 + "_clu";
        alnIndexFile = par.db1 + "_clu.index";
    }


    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    dbr.open(DBReader<unsigned int>::NOSORT);

    // remove par.indexDbsuffix from db1
    std::string seqDb = par.db1 + "_seq";
    std::string seqDbIndex = par.db1 + "_seq.index";
    std::string seqDbtypeFile = par.db1 + "_seq.dbtype";
    if (par.indexDbsuffix != "") {
        std::string::size_type pos = par.db1.find(par.indexDbsuffix);
        if (pos != std::string::npos) {
            par.db1 = par.db1.substr(0, pos);
        }
        seqDb         = par.db1 + "_seq" + par.indexDbsuffix;
        seqDbIndex    = par.db1 + "_seq" + par.indexDbsuffix + ".index";
        seqDbtypeFile = par.db1 + "_seq" + par.indexDbsuffix + ".dbtype";
    }
    const bool ppDB = FileUtil::fileExists(alnDbtypeFile.c_str()) && FileUtil::fileExists(seqDbtypeFile.c_str());

    std::string db2 = ppDB ? seqDb : par.db2;
    std::string db2Index = ppDB ? seqDbIndex : par.db2Index;

    std::string hdr1 = ppDB ? seqDb + "_h" : par.hdr1;
    std::string hdr1Index = ppDB ?  seqDb + "_h.index" : par.hdr1Index;

    DBReader<unsigned int> *dbr2 = NULL;
    if ((sameDB == false) || ppDB) {
        dbr2 = new DBReader<unsigned int>(db2.c_str(), db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        dbr2->open(DBReader<unsigned int>::NOSORT);
    }

    const bool db1IsNucl = Parameters::isEqualDbtype(dbr.getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    const bool db2IsNucl = dbr2 != NULL && Parameters::isEqualDbtype(dbr2->getDbtype(), Parameters::DBTYPE_NUCLEOTIDES);
    BaseMatrix *seedSubMat = Prefiltering::getSubstitutionMatrix(par.seedScoringMatrixFile, par.alphabetSize, 8.0f, false, (db1IsNucl && db2IsNucl));

    // memoryLimit in bytes
    size_t memoryLimit = Util::computeMemory(par.splitMemoryLimit);

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
        par.kmerScore.values = 0;
    }

    const bool contextPseudoCnts = DBReader<unsigned int>::getExtendedDbtype(dbr.getDbtype()) & Parameters::DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS;

    // TODO: investigate if it makes sense to mask the profile consensus sequence
    if (isProfileSearch) {
        par.maskMode = 0;
    }

    // query seq type is actually unknown here, but if we pass DBTYPE_HMM_PROFILE then its +20 k-score
    int kmerScore = Prefiltering::getKmerThreshold(par.sensitivity, isProfileSearch, contextPseudoCnts, par.kmerScore.values, par.kmerSize);

    const std::string& baseDB = ppDB ? par.db1 + par.indexDbsuffix : par.db2;
    std::string indexDB = PrefilteringIndexReader::indexName(baseDB);

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

    const bool noHeaders = (par.indexSubset & Parameters::INDEX_SUBSET_NO_HEADERS) != 0;
    if (recreate) {
        DBReader<unsigned int> *hdbr1 = NULL;
        if (noHeaders == false) {
            hdbr1 = new DBReader<unsigned int>(hdr1.c_str(), hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
            hdbr1->open(DBReader<unsigned int>::NOSORT);
        }

        DBReader<unsigned int> *hdbr2 = NULL;
        if (sameDB == false && ppDB == false && noHeaders == false) {
            hdbr2 = new DBReader<unsigned int>(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
            hdbr2->open(DBReader<unsigned int>::NOSORT);
        }

        DBReader<unsigned int> *alndbr = NULL;
        const bool noAlignment = (par.indexSubset & Parameters::INDEX_SUBSET_NO_ALIGNMENT) != 0;
        if (ppDB == true && noAlignment == false) {
            alndbr = new DBReader<unsigned int>(alnFile.c_str(), alnIndexFile.c_str(),
                                                par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
            alndbr->open(DBReader<unsigned int>::NOSORT);
        }

        DBReader<unsigned int>::removeDb(indexDB);
        PrefilteringIndexReader::createIndexFile(indexDB, &dbr, dbr2, hdbr1, hdbr2, alndbr, seedSubMat, par.maxSeqLen,
                                                 par.spacedKmer, par.spacedKmerPattern, par.compBiasCorrection,
                                                 seedSubMat->alphabetSize, par.kmerSize, par.maskMode, par.maskLowerCaseMode,
                                                 par.maskProb, kmerScore, par.targetSearchMode, par.split, par.indexSubset);

        if (alndbr != NULL) {
            alndbr->close();
            delete alndbr;
        }

        if (hdbr2 != NULL) {
            hdbr2->close();
            delete hdbr2;
        }

        if (hdbr1 != NULL) {
            hdbr1->close();
            delete hdbr1;
        }
    }

    if (dbr2 != NULL) {
        dbr2->close();
        delete dbr2;
    }

    delete seedSubMat;
    dbr.close();

    return status;
}
