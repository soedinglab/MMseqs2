#include "Alignment.h"

Alignment::Alignment(std::string seqDB, std::string prefDB, std::string outDB, std::string matrixFile, double evalThr, double covThr, int maxSeqLen, int seqType){

    BUFFER_SIZE = 10000000;

    this->covThr = covThr;

    this->evalThr = evalThr;

    std::cout << "Initialising data structures...\n";
    std::string seqDBIndex = seqDB + ".index";
    std::string prefDBIndex = prefDB + ".index";
    std::string outDBIndex = outDB + ".index";

    BaseMatrix* m;
    if (seqType == Sequence::AMINO_ACIDS)
        m = new SubstitutionMatrix(matrixFile.c_str(), 2.0);
    else
        m = new NucleotideMatrix();

    threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
    std::cout << "Using " << threads << " threads.\n";
#endif

    qSeqs = new Sequence*[threads];
    dbSeqs = new Sequence*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        qSeqs[i] = new Sequence(maxSeqLen, m->aa2int, m->int2aa, seqType);
        dbSeqs[i] = new Sequence(maxSeqLen, m->aa2int, m->int2aa, seqType);
    }

    matchers = new Matcher*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        matchers[i] = new Matcher(m, maxSeqLen);

    // open the sequence, prefiltering and output databases
    seqdbr = new DBReader(seqDB.c_str(), seqDBIndex.c_str());
    seqdbr->open(DBReader::NOSORT);

    prefdbr = new DBReader(prefDB.c_str(), prefDBIndex.c_str());
    prefdbr->open(DBReader::NOSORT);

    dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), threads);
    dbw->open();

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
    seqdbr->close();
    prefdbr->close();
    dbw->close();

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
}

void Alignment::run (int maxAlnNum){

    int alignmentsNum = 0;
    int passedNum = 0;
# pragma omp parallel for schedule(dynamic, 10) reduction (+: alignmentsNum, passedNum)
    for (unsigned int id = 0; id < prefdbr->getSize(); id++){
        if (id % 1000000 == 0 && id > 0){
            std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";
            fflush(stdout);
        }
        else if (id % 10000 == 0 && id > 0) {
            std::cout << ".";
            fflush(stdout);
        }

        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        // get the prefiltering list
        char* prefList = prefdbr->getData(id);
        char* queryDbKey = prefdbr->getDbKey(id);
        std::stringstream lineSs (prefList);
        std::string val;

        // map the query sequence
        char* querySeqData = seqdbr->getDataByDBKey(queryDbKey);
        qSeqs[thread_idx]->mapSequence(id, queryDbKey, querySeqData);

        // parse the prefiltering list and calculate a Smith-Waterman alignment for each sequence in the list 
        std::list<Matcher::result_t>* swResults = new std::list<Matcher::result_t>();
        int cnt = 0;
        while (std::getline(lineSs, val, '\t') && cnt < maxAlnNum){
            // DB key of the db sequence
            for (unsigned int j = 0; j < val.length(); j++)
                dbKeys[thread_idx][j] = val.at(j);
            dbKeys[thread_idx][val.length()] = '\0';
            // prefiltering score
            std::getline(lineSs, val, '\t');
            float prefScore = atof(val.c_str());
            // prefiltering e-value
            std::getline(lineSs, val, '\n');
            float prefEval = atof(val.c_str());

            // map the database sequence
            char* dbSeqData = seqdbr->getDataByDBKey(dbKeys[thread_idx]);
            if (dbSeqData == NULL){
# pragma omp critical
                {
                    std::cerr << "ERROR: Sequence " << dbKeys[thread_idx] << " required in the prefiltering not contained in the input sequence database! Please check your database.\n";
                    exit(1);
                }
            }
            dbSeqs[thread_idx]->mapSequence(-1, dbKeys[thread_idx], dbSeqData);

            // check if the sequences could pass the coverage threshold 
            if ( (((float) qSeqs[thread_idx]->L) / ((float) dbSeqs[thread_idx]->L) < covThr) ||
                    (((float) dbSeqs[thread_idx]->L) / ((float) qSeqs[thread_idx]->L) < covThr) )
                continue;

            // calculate Smith-Waterman alignment
            Matcher::result_t res = matchers[thread_idx]->getSWResult(qSeqs[thread_idx], dbSeqs[thread_idx], seqdbr->getSize());
            swResults->push_back(res);

            cnt++;
        }

        // write the results
        swResults->sort(Matcher::compareHits);
        std::list<Matcher::result_t>::iterator it;
        std::stringstream swResultsSs;

        // put the contents of the swResults list into ffindex DB
        for (it = swResults->begin(); it != swResults->end(); ++it){
            if (it->eval <= evalThr && it->qcov >= covThr && it->dbcov >= covThr){
                swResultsSs << it->dbKey << "\t" << it->score << "\t" << it->qcov << "\t" << it->dbcov << "\t" << it->eval << "\n";
                passedNum++;
            }
            alignmentsNum++;
        }
        std::string swResultsString = swResultsSs.str();
        const char* swResultsStringData = swResultsString.c_str();
        if (BUFFER_SIZE <= swResultsString.length()){
            std::cerr << "Output buffer size < result size! (" << BUFFER_SIZE << " <= " << swResultsString.length() << ")\nIncrease buffer size or reconsider your parameters - output buffer is already huge ;-)\n";
            exit(1);
        }
        memcpy(outBuffers[thread_idx], swResultsStringData, swResultsString.length()*sizeof(char));
        dbw->write(outBuffers[thread_idx], swResultsString.length(), qSeqs[thread_idx]->getDbKey(), thread_idx);

        delete swResults;

    }
    std::cout << "\n";
    std::cout << "All sequences processed.\n\n";
    std::cout << alignmentsNum << " alignments calculated.\n";
    std::cout << passedNum << " sequences passed the thresholds (" << ((float)passedNum/(float)alignmentsNum) << " of overall calculated).\n";
    std::cout << ((float)passedNum/(float)prefdbr->getSize()) << " hits per query sequence.\n";

}

