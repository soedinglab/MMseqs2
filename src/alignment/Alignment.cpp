#include "Alignment.h"


Alignment::Alignment(std::string querySeqDB, std::string querySeqDBIndex, 
        std::string targetSeqDB, std::string targetSeqDBIndex,
        std::string prefDB, std::string prefDBIndex, 
        std::string outDB, std::string outDBIndex,
        Parameters par){

    BUFFER_SIZE = 10000000;

    this->covThr = par.covThr;
    this->evalThr = par.evalThr;
    this->seqIdThr = par.seqIdThr;
    if(this->covThr == 0.0 && this->seqIdThr == 0.0){
        Debug(Debug::WARNING) << "Compute score only.\n";
        this->mode = Matcher::SCORE_ONLY; //fastest
    } else if(this->covThr > 0.0 && this->seqIdThr == 0.0) {
        Debug(Debug::WARNING) << "Compute score and coverage.\n";
        this->mode = Matcher::SCORE_COV; // fast
    } else { // if seq id is needed
        Debug(Debug::WARNING) << "Compute score, coverage and sequence id.\n";
        this->mode = Matcher::SCORE_COV_SEQID; // slowest
    }
    if (par.querySeqType == Sequence::AMINO_ACIDS || par.querySeqType == Sequence::HMM_PROFILE)
        this->m = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0);
    else
        this->m = new NucleotideMatrix();

    threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
    Debug(Debug::INFO) << "Using " << threads << " threads.\n";
#endif

    qSeqs = new Sequence*[threads];
    dbSeqs = new Sequence*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        qSeqs[i]  = new Sequence(par.maxSeqLen, m, par.querySeqType, 0, false);
        dbSeqs[i] = new Sequence(par.maxSeqLen, m, par.targetSeqType, 0, false);
    }



    // open the sequence, prefiltering and output databases
    qseqdbr = new DBReader(querySeqDB.c_str(), querySeqDBIndex.c_str());
    qseqdbr->open(DBReader::NOSORT);

    tseqdbr = new DBReader(targetSeqDB.c_str(), targetSeqDBIndex.c_str());
    tseqdbr->open(DBReader::NOSORT);

    prefdbr = new DBReader(prefDB.c_str(), prefDBIndex.c_str());
    prefdbr->open(DBReader::NOSORT);

    dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), threads);
    dbw->open();

    matchers = new Matcher*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++) {
        matchers[i] = new Matcher(par.maxSeqLen, this->m);
    }

    dbKeys = new char*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        dbKeys[i] = new char[100];

    outBuffers = new char*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        outBuffers[i] = new char[BUFFER_SIZE];

}

Alignment::~Alignment(){
    for (int i = 0; i < threads; i++){
        delete qSeqs[i];
        delete dbSeqs[i];
        delete matchers[i];
        delete[] dbKeys[i];
        delete[] outBuffers[i];
    }

    delete[] qSeqs;
    delete[] dbSeqs;
    delete[] matchers;
    delete[] dbKeys;
    delete[] outBuffers;

    delete m;

    delete qseqdbr;
    delete tseqdbr;
    delete prefdbr;
    delete dbw;
}

void Alignment::run (const unsigned int maxAlnNum, const unsigned int maxRejected){

    size_t alignmentsNum = 0;
    size_t totalPassedNum = 0;

# pragma omp parallel for schedule(dynamic, 10) reduction (+: alignmentsNum, totalPassedNum)
    for (size_t id = 0; id < prefdbr->getSize(); id++){
        Log::printProgress(id);

        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        // get the prefiltering list
        char* prefList = prefdbr->getData(id);
        char* queryDbKey = prefdbr->getDbKey(id);

        // map the query sequence
        char* querySeqData = qseqdbr->getDataByDBKey(queryDbKey);
        qSeqs[thread_idx]->mapSequence(id, queryDbKey, querySeqData);
        matchers[thread_idx]->initQuery(qSeqs[thread_idx]);
        // parse the prefiltering list and calculate a Smith-Waterman alignment for each sequence in the list 
        std::list<Matcher::result_t>* swResults = new std::list<Matcher::result_t>();
        std::stringstream lineSs (prefList);
        std::string val;
        size_t passedNum = 0;
        unsigned int rejected = 0;
        while (std::getline(lineSs, val, '\t') && passedNum < maxAlnNum && rejected < maxRejected){
            // DB key of the db sequence
            for (unsigned int j = 0; j < val.length(); j++)
                dbKeys[thread_idx][j] = val.at(j);
            dbKeys[thread_idx][val.length()] = '\0';
            // prefiltering score
            std::getline(lineSs, val, '\t');
            //float prefScore = atof(val.c_str());
            // prefiltering e-value
            std::getline(lineSs, val, '\n');
            //float prefEval = atof(val.c_str());

            // map the database sequence
            char* dbSeqData = tseqdbr->getDataByDBKey(dbKeys[thread_idx]);
            if (dbSeqData == NULL){
# pragma omp critical
                {
                    Debug(Debug::ERROR) << "ERROR: Sequence " << dbKeys[thread_idx]
                            << " is required in the prefiltering, but is not contained in the input sequence database!\n" <<
                               "Please check your database.\n";
                    EXIT(1);
                }
            }
            dbSeqs[thread_idx]->mapSequence(-1, dbKeys[thread_idx], dbSeqData);

            // check if the sequences could pass the coverage threshold 
            if ( (((float) qSeqs[thread_idx]->L) / ((float) dbSeqs[thread_idx]->L) < covThr) ||
                 (((float) dbSeqs[thread_idx]->L) / ((float) qSeqs[thread_idx]->L) < covThr) ) {
                rejected++;
                continue;
            }

            // calculate Smith-Waterman alignment
            Matcher::result_t res = matchers[thread_idx]->getSWResult(dbSeqs[thread_idx], tseqdbr->getAminoAcidDBSize(), evalThr, this->mode);
            alignmentsNum++;

            if ((res.eval <= evalThr || (mode != Matcher::SCORE_ONLY && res.seqId == 1.0)) &&
                (res.qcov >= covThr && res.dbcov >= covThr) && //TODO is qcov not enough? maybe two thresholds?
                (res.seqId > seqIdThr)) {
                swResults->push_back(res);
                passedNum++;
                totalPassedNum++;
            }
            else{
                rejected++;
            }
        }
        // write the results
        swResults->sort(Matcher::compareHits);
        std::list<Matcher::result_t>::iterator it;
        std::stringstream swResultsSs;

        // put the contents of the swResults list into ffindex DB
        for (it = swResults->begin(); it != swResults->end(); ++it){
                swResultsSs << it->dbKey << "\t";
                swResultsSs << it->score << "\t"; //TODO fix for formats
                swResultsSs << std::fixed << std::setprecision(3) << it->qcov << "\t";
                swResultsSs << it->dbcov << "\t";
                swResultsSs << it->seqId << "\t";
                swResultsSs << std::scientific << it->eval << "\n";
       }
        std::string swResultsString = swResultsSs.str();
        const char* swResultsStringData = swResultsString.c_str();
        if (BUFFER_SIZE <= swResultsString.length()){
            Debug(Debug::ERROR) << "Output buffer size < result size! ("
                                << BUFFER_SIZE << " <= " << swResultsString.length()
                                << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
            EXIT(1);
        }
        memcpy(outBuffers[thread_idx], swResultsStringData, swResultsString.length()*sizeof(char));
        dbw->write(outBuffers[thread_idx], swResultsString.length(), qSeqs[thread_idx]->getDbKey(), thread_idx);

        delete swResults;

    }
    Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "All sequences processed.\n\n";
    Debug(Debug::INFO) << alignmentsNum << " alignments calculated.\n";
    Debug(Debug::INFO) << totalPassedNum << " sequence pairs passed the thresholds (" << ((float)totalPassedNum/(float)alignmentsNum) << " of overall calculated).\n";
    size_t hits = totalPassedNum / prefdbr->getSize();
    size_t hits_rest = totalPassedNum % prefdbr->getSize();
    float hits_f = ((float) hits) + ((float)hits_rest) / (float) prefdbr->getSize();
    Debug(Debug::INFO) << hits_f << " hits per query sequence.\n";

    qseqdbr->close();
    tseqdbr->close();
    prefdbr->close();
    dbw->close();

}

