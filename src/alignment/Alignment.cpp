#include "Alignment.h"
#include "BlastScoreUtils.h"
#include "Util.h"
#include "Debug.h"

#include <list>
#include <iomanip>
#include <memory>

#ifdef OPENMP
#include <omp.h>
#endif


Alignment::Alignment(std::string& querySeqDB, std::string& querySeqDBIndex,
                     std::string& targetSeqDB, std::string& targetSeqDBIndex,
                     std::string& prefDB, std::string& prefDBIndex,
                     std::string& outDB, std::string& outDBIndex,
                     Parameters &par){

    this->covThr = par.covThr;
    this->evalThr = par.evalThr;
    this->seqIdThr = par.seqIdThr;
    this->fragmentMerge = par.fragmentMerge;
    this->addBacktrace = par.addBacktrace;
    this->realign = par.realign;
    if(addBacktrace == true){
        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }
    if(realign == true){
        par.alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_ONLY;
        if(addBacktrace == false){
            Debug(Debug::ERROR) << "Realign is just useful in combination with -a.\n";
            EXIT(EXIT_FAILURE);
        }
    }
    this->swMode = Matcher::SCORE_ONLY; // fast
    switch (par.alignmentMode){
        case Parameters::ALIGNMENT_MODE_FAST_AUTO:
            if(this->covThr > 0.0 && this->seqIdThr == 0.0) {
                this->swMode = Matcher::SCORE_COV; // fast
            } else if(this->covThr > 0.0 && this->seqIdThr > 0.0) { // if seq id is needed
                this->swMode = Matcher::SCORE_COV_SEQID; // slowest
            }
            break;
        case Parameters::ALIGNMENT_MODE_SCORE_COV:
            this->swMode = Matcher::SCORE_COV; // fast
            break;
        case Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID:
            this->swMode = Matcher::SCORE_COV_SEQID; // fast
            break;
    }
    // print out mode and check for errors
    switch(swMode){
        case Matcher::SCORE_ONLY:
            Debug(Debug::WARNING) << "Compute score only.\n";
            if(fragmentMerge == true){
                Debug(Debug::ERROR) << "Fragment merge does not work with Score only mode. Set --alignment-mode to 2 or 3.\n";
                EXIT(EXIT_FAILURE);
            }
            break;
        case Matcher::SCORE_COV:
            Debug(Debug::WARNING) << "Compute score and coverage.\n";
            break;
        case Matcher::SCORE_COV_SEQID:
            Debug(Debug::WARNING) << "Compute score, coverage and sequence id.\n";
            break;
        default:
            Debug(Debug::ERROR) << "Wronge swMode mode.\n";
            EXIT(EXIT_FAILURE);
            break;
    }

    if (par.querySeqType == Sequence::AMINO_ACIDS || par.querySeqType == Sequence::HMM_PROFILE){
        // keep score bias to 0.0 (improved ROC)
        this->m = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    }else{
        this->m = new NucleotideMatrix();
    }

    threads = par.threads;
#ifdef OPENMP
    Debug(Debug::INFO) << "Using " << threads << " threads.\n";
#endif

    qSeqs = new Sequence*[threads];
    dbSeqs = new Sequence*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        qSeqs[i]  = new Sequence(par.maxSeqLen, m->aa2int, m->int2aa, par.querySeqType,  0, false, par.compBiasCorrection);
        dbSeqs[i] = new Sequence(par.maxSeqLen, m->aa2int, m->int2aa, par.targetSeqType, 0, false, par.compBiasCorrection);
    }

    // open the sequence, prefiltering and output databases
    qseqdbr = new DBReader<unsigned int>(querySeqDB.c_str(), querySeqDBIndex.c_str());
    qseqdbr->open(DBReader<unsigned int>::NOSORT);
    qseqdbr->readMmapedDataInMemory();
    sameQTDB = (querySeqDB.compare(targetSeqDB) == 0);
    if(sameQTDB == true) {
        this->tseqdbr = qseqdbr;
    }else{
        tseqdbr = new DBReader<unsigned int>(targetSeqDB.c_str(), targetSeqDBIndex.c_str());
        tseqdbr->open(DBReader<unsigned int>::NOSORT);
        tseqdbr->readMmapedDataInMemory();
    }

    includeIdentity = par.includeIdentity;
    prefdbr = new DBReader<unsigned int>(prefDB.c_str(), prefDBIndex.c_str());
    prefdbr->open(DBReader<unsigned int>::LINEAR_ACCCESS);

    this->outDB = outDB;
    this->outDBIndex = outDBIndex;

    matchers = new Matcher*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++) {
        matchers[i] = new Matcher(par.maxSeqLen, this->m, tseqdbr->getAminoAcidDBSize(), tseqdbr->getSize(),
                                  par.compBiasCorrection);
    }
    if(this->realign == true){
        this->realign_m = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, -0.2);
        this->realigner = new Matcher*[threads];
# pragma omp parallel for schedule(static)
        for (int i = 0; i < threads; i++) {
            realigner[i] = new Matcher(par.maxSeqLen, this->realign_m, tseqdbr->getAminoAcidDBSize(),
                                       tseqdbr->getSize(), par.compBiasCorrection);
        }
    }
    dbKeys = new unsigned int[threads];
}

Alignment::~Alignment(){
    if(this->realign == true) {
        for (int i = 0; i < threads; i++) {
            delete realigner[i];
        }
        delete [] realigner;
        delete realign_m;
    }
    for (int i = 0; i < threads; i++){
        delete qSeqs[i];
        delete dbSeqs[i];
        delete matchers[i];
    }
    delete[] qSeqs;
    delete[] dbSeqs;
    delete[] matchers;
    delete[] dbKeys;
    delete m;
    delete qseqdbr;
    if(sameQTDB == false){
        delete tseqdbr;
    }
    delete prefdbr;
}

void Alignment::run (const unsigned int mpiRank, const unsigned int mpiNumProc,
                     const unsigned int maxAlnNum, const unsigned int maxRejected) {

    size_t dbFrom = 0;
    size_t dbSize = 0;
    Util::decomposeDomainByAminoAcid(qseqdbr->getAminoAcidDBSize(), qseqdbr->getSeqLens(), qseqdbr->getSize(),
                                     mpiRank, mpiNumProc, &dbFrom, &dbSize);
    Debug(Debug::WARNING) << "Compute split from " << dbFrom << " to " << dbFrom+dbSize << "\n";
    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outDB, outDBIndex, mpiRank);
    run(tmpOutput.first.c_str(), tmpOutput.second.c_str(), dbFrom, dbSize, maxAlnNum, maxRejected);

    // close reader to reduce memory
    this->closeReader();
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(mpiRank == 0){ // master reduces results
        std::vector<std::pair<std::string, std::string> > splitFiles;
        for(unsigned int procs = 0; procs < mpiNumProc; procs++){
            splitFiles.push_back(Util::createTmpFileNames(outDB, outDBIndex, procs));
        }
        // merge output ffindex databases
        this->mergeAndRemoveTmpDatabases(outDB, outDBIndex, splitFiles);
    }
}


void Alignment::closeReader(){
    qseqdbr->close();
    if(sameQTDB == false) {
        tseqdbr->close();
    }
    prefdbr->close();
}

void Alignment::run(const unsigned int maxAlnNum, const unsigned int maxRejected){
    run(outDB.c_str(), outDBIndex.c_str(), 0, prefdbr->getSize(), maxAlnNum, maxRejected);
    this->closeReader();
}

void Alignment::run (const char * outDB, const char * outDBIndex,
                     const size_t dbFrom, const size_t dbSize,
                     const unsigned int maxAlnNum, const unsigned int maxRejected){

    size_t alignmentsNum = 0;
    size_t totalPassedNum = 0;
    DBWriter dbw(outDB, outDBIndex, threads);
    dbw.open();
    const size_t flushSize = 1000000;
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbSize)/static_cast<double>(flushSize)));
    for(size_t i = 0; i < iterations; i++) {
        size_t start = dbFrom + (i * flushSize);
        size_t bucketSize = std::min(dbSize - (i * flushSize), flushSize);
# pragma omp parallel for schedule(dynamic, 100) reduction (+: alignmentsNum, totalPassedNum)
        for (size_t id = start; id < (start + bucketSize); id++){
            Debug::printProgress(id);
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            // get the prefiltering list
            char* data = prefdbr->getData(id);
            unsigned int queryDbKey = prefdbr->getDbKey(id);
            // map the query sequence
            char* querySeqData = qseqdbr->getDataByDBKey(queryDbKey);
            if (querySeqData == NULL){
# pragma omp critical
                {
                    Debug(Debug::ERROR) << "ERROR: Query sequence " << queryDbKey
                    << " is required in the prefiltering, but is not contained in the query sequence database!\n" <<
                    "Please check your database.\n";
                    EXIT(1);
                }
            }

            qSeqs[thread_idx]->mapSequence(id, queryDbKey, querySeqData);
            matchers[thread_idx]->initQuery(qSeqs[thread_idx]);
            // parse the prefiltering list and calculate a Smith-Waterman alignment for each sequence in the list
            std::vector<Matcher::result_t> swResults;
            size_t passedNum = 0;
            unsigned int rejected = 0;
            while (*data != '\0' && passedNum < maxAlnNum && rejected < maxRejected){
                // DB key of the db sequence
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                dbKeys[thread_idx] = dbKey;
                // sequence are identical if qID == dbID  (needed to cluster really short sequences)
                const bool isIdentity = (queryDbKey == dbKey && (includeIdentity || sameQTDB))? true : false;
                // map the database sequence
                char* dbSeqData = tseqdbr->getDataByDBKey(dbKeys[thread_idx]);
                if (dbSeqData == NULL){
# pragma omp critical
                    {
                        Debug(Debug::ERROR) << "ERROR: Sequence " << dbKeys[thread_idx]
                        << " is required in the prefiltering, but is not contained in the target sequence database!\n" <<
                        "Please check your database.\n";
                        EXIT(EXIT_FAILURE);
                    }
                }
                //char *maskedDbSeq = seg[thread_idx]->maskseq(dbSeqData);
                dbSeqs[thread_idx]->mapSequence(-1, dbKeys[thread_idx], dbSeqData);
                // check if the sequences could pass the coverage threshold
                if(fragmentMerge == false){
                    if ( (((float) qSeqs[thread_idx]->L) / ((float) dbSeqs[thread_idx]->L) < covThr) ||
                         (((float) dbSeqs[thread_idx]->L) / ((float) qSeqs[thread_idx]->L) < covThr) ) {
                        rejected++;
                        data = Util::skipLine(data);
                        continue;
                    }
                }
                // calculate Smith-Waterman alignment
                Matcher::result_t res = matchers[thread_idx]->getSWResult(dbSeqs[thread_idx], tseqdbr->getSize(), evalThr, this->swMode);
                alignmentsNum++;
                //set coverage and seqid if identity
                if (isIdentity){
                    res.qcov=1.0f;
                    res.dbcov=1.0f;
                    res.seqId=1.0f;
                }

                // check first if it is identity
                if (isIdentity
                    ||
                    // general accaptance criteria
                    (res.eval <= evalThr && res.seqId >= seqIdThr && res.qcov >= covThr && res.dbcov >= covThr)
                    ||
                    // check for fragment
                    (( swMode == Matcher::SCORE_COV_SEQID || swMode == Matcher::SCORE_COV) && fragmentMerge == true && res.dbcov >= 0.95 && res.seqId >= 0.9 ))
                {
                    swResults.push_back(res);
                    passedNum++;
                    totalPassedNum++;
                    rejected = 0;
                }
                else{
                    rejected++;
                }
                data = Util::skipLine(data);
            }
            // write the results
            std::sort(swResults.begin(), swResults.end(), Matcher::compareHits);
            if(realign == true){
                realigner[thread_idx]->initQuery(qSeqs[thread_idx]);
                for (size_t i = 0; i < swResults.size(); i++) {
                    char* dbSeqData = tseqdbr->getDataByDBKey(swResults[i].dbKey);
                    dbSeqs[thread_idx]->mapSequence(-1, swResults[i].dbKey, dbSeqData);
                    Matcher::result_t res = realigner[thread_idx]->getSWResult(dbSeqs[thread_idx],
                                                                               tseqdbr->getSize(), 0.0,
                                                                               Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID);
                    swResults[i].backtrace = res.backtrace;
                    swResults[i].qStartPos = res.qStartPos;
                    swResults[i].qEndPos = res.qEndPos;
                    swResults[i].dbStartPos = res.dbStartPos;
                    swResults[i].dbEndPos = res.dbEndPos;
                }
            }
            // put the contents of the swResults list into ffindex DB
            std::stringstream swResultsString;
            for (size_t i = 0; i < swResults.size(); i++){
                swResultsString << Matcher::resultToString(swResults[i], addBacktrace);
            }

            std::string swResultString = swResultsString.str();
            const char* swResultsData = swResultString.c_str();
            dbw.writeData(swResultsData, swResultString.length(), SSTR(qSeqs[thread_idx]->getDbKey()).c_str(),
                          thread_idx);
            swResults.clear();
//        prefdbr->unmapDataById(id);
        }
        prefdbr->remapData();
    }
    dbw.close();
    Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "All sequences processed.\n\n";
    Debug(Debug::INFO) << alignmentsNum << " alignments calculated.\n";
    Debug(Debug::INFO) << totalPassedNum << " sequence pairs passed the thresholds (" << ((float)totalPassedNum/(float)alignmentsNum) << " of overall calculated).\n";
    size_t hits = totalPassedNum / dbSize;
    size_t hits_rest = totalPassedNum % dbSize;
    float hits_f = ((float) hits) + ((float)hits_rest) / (float) dbSize;
    Debug(Debug::INFO) << hits_f << " hits per query sequence.\n";
}

void Alignment::mergeAndRemoveTmpDatabases(const std::string& out,const  std::string& outIndex,
                                           const std::vector<std::pair<std::string, std::string >>& files) {
    const char ** datafilesNames = new const char*[files.size()];
    const char ** indexFilesNames= new const char*[files.size()];
    for(size_t i = 0; i < files.size(); i++){
        datafilesNames[i] = files[i].first.c_str();
        indexFilesNames[i] = files[i].second.c_str();
    }
    DBWriter::mergeResults(out.c_str(), outIndex.c_str(), datafilesNames, indexFilesNames, files.size());
    delete [] datafilesNames;
    delete [] indexFilesNames;
}
