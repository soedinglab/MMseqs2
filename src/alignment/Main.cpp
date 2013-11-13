#ifdef OPENMP
#include <omp.h>
#endif

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include <string>
#include <list>

#include "../commons/DBReader.h"
#include "../commons/DBWriter.h"
#include "../commons/NucleotideMatrix.h"
#include "../commons/SubstitutionMatrix.h"

#include "Matcher.h"

void printUsage(){

    std::string usage("\nCalculates Smith-Waterman alignment scores.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: kClust2_pref ffindexSequenceDBBase ffindexPrefilteringDBBase ffindexOutDBBase [opts]\n"
            "-m\t[file]\tAmino acid substitution matrix file.\n"
            "-e\t[float]\tMaximum e-value (default=0.01).\n"
            "-c\t[float]\tMinimum alignment coverage (default=0.8).\n"
            "-s\t[int]\tMaximum sequence length (default=50000).\n"
            "-r\t[int]\tMaximum result alignment number per query sequence (default=10).\n"
            "-n\t\tNucleotide sequences input.\n");
    std::cout << usage;
}

void parseArgs(int argc, char** argv, std::string* seqDB, std::string* prefDB, std::string* matrixFile, std::string* outDB, double* evalThr, double* covThr, int* maxSeqLen, int* maxAlnNum, int* nucl){
    if (argc < 4){
        printUsage();
        exit(EXIT_FAILURE);
    }

    seqDB->assign(argv[1]);
    prefDB->assign(argv[2]);
    outDB->assign(argv[3]);
    int i = 4;
    while (i < argc){
        if (strcmp(argv[i], "-m") == 0){
            if (*nucl == 1){
                std::cerr << "No scoring matrix is allowed for nucleotide sequences.\n";
                exit(EXIT_FAILURE);
            }
            if (++i < argc){
                matrixFile->assign(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -m\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-e") == 0){
            if (++i < argc){
                *evalThr = atof(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -e\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-c") == 0){
            if (++i < argc){
                *covThr = atof(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -c\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-s") == 0){
            if (++i < argc){
                *maxSeqLen = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -s\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-r") == 0){
            if (++i < argc){
                *maxAlnNum = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -a\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-n") == 0){
            if (strcmp(matrixFile->c_str(), "") != 0){
                std::cerr << "No scoring matrix is allowed for nucleotide sequences.\n";
                exit(EXIT_FAILURE);
            }
            *nucl = 1;
            i++;
        }
        else {
            printUsage();
            std::cerr << "Wrong argument: " << argv[i] << "\n";
            exit(EXIT_FAILURE);
        }
    }
}

bool compareHits (Matcher::result_t first, Matcher::result_t second){
    if (first.score > second.score)
        return true;
    return false;
}

int main(int argc, char **argv){
    std::string seqDB = "";
    std::string prefDB = ""; 
    std::string matrixFile = ""; 
    std::string outDB = "";

    double evalThr = 0.001;
    double covThr = 0.8;
    int maxSeqLen = 50000;
    int maxAlnNum = 10;
    size_t BUFFER_SIZE = 10000000;
    int nucleotides = 0;

    parseArgs(argc, argv, &seqDB, &prefDB, &matrixFile, &outDB, &evalThr, &covThr, &maxSeqLen, &maxAlnNum, &nucleotides);

    std::cout << "Calculating Smith-Waterman alignments with following parameters:\nmax. evalue:\t" << evalThr 
        << "\nmin. sequence coverage:\t" << covThr 
        << "\nmax. sequence length:\t" << maxSeqLen 
        << "\nmax. alignment number per query:\t" << maxAlnNum << "\n\n";

    std::cout << "Initialising data structures...\n";
    std::string seqDBIndex = seqDB + ".index";
    std::string prefDBIndex = prefDB + ".index";
    std::string outDBIndex = outDB + ".index";

    BaseMatrix* m;
    if (nucleotides == 0)
        m = new SubstitutionMatrix(matrixFile.c_str(), 2.0);
    else
        m = new NucleotideMatrix();

    int threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
#endif

    // 2 Sequence objects for each thread: one for the query, one for the DB sequence
    Sequence** qSeqs = new Sequence*[threads];
    Sequence** dbSeqs = new Sequence*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        qSeqs[i] = new Sequence(maxSeqLen, m->aa2int, m->int2aa, nucleotides);
        dbSeqs[i] = new Sequence(maxSeqLen, m->aa2int, m->int2aa, nucleotides);
    }

    Matcher** matchers = new Matcher*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        matchers[i] = new Matcher(m, maxSeqLen);

    // open the sequence, prefiltering and output databases
    DBReader* seqdbr = new DBReader(seqDB.c_str(), seqDBIndex.c_str());
    seqdbr->open();

    DBReader* prefdbr = new DBReader(prefDB.c_str(), prefDBIndex.c_str());
    prefdbr->open();

    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), threads);
    dbw->open();

    // buffers for the database keys (needed during the processing of the prefilterings lists)
    char** dbKeys = new char*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        dbKeys[i] = new char[100];

    // output buffers
    char** outBuffers = new char*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++)
        outBuffers[i] = new char[BUFFER_SIZE];

    std::cout << "Starting Smith-Waterman scores calculation.\n";

    int alignmentsNum = 0;
    int passedNum = 0;
# pragma omp parallel for schedule(dynamic, 10) reduction (+: alignmentsNum, passedNum)
    for (unsigned int id = 0; id < prefdbr->getSize(); id++){
        if (id % 1000000 == 0 && id > 0)
            std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";

        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        // get the prefiltering list
        char* prefList = prefdbr->getData(id);
        char* queryDbKey = prefdbr->getDbKey(id);
        std::cout << "query: " << queryDbKey << "\n";
        std::stringstream lineSs (prefList);
        std::string val;

        // map the query sequence
        char* querySeqData = seqdbr->getDataByDBKey(queryDbKey);
        qSeqs[thread_idx]->mapSequence(id, queryDbKey, querySeqData);

        // parse the prefiltering list and calculate a Smith-Waterman alignment for each sequence in the list 
        std::list<Matcher::result_t>* swResults = new std::list<Matcher::result_t>();
        int cnt = 0;
        while (std::getline(lineSs, val, '\t') && cnt < 10){
            // DB key of the db sequence
            for (unsigned int j = 0; j < val.length(); j++)
                dbKeys[thread_idx][j] = val.at(j);
            dbKeys[thread_idx][val.length()] = '\0';
            std::cout << "\t" << dbKeys[thread_idx] << "\n";
//            if (strcmp (dbKeys[thread_idx], "UniRef100_E9IKE5") != 0)
//                continue;
            // prefiltering score
            std::getline(lineSs, val, '\t');
            float prefScore = atof(val.c_str());
            // prefiltering e-value
            std::getline(lineSs, val, '\n');
            float prefEval = atof(val.c_str());

            // map the database sequence
            char* dbSeqData = seqdbr->getDataByDBKey(dbKeys[thread_idx]);
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
        swResults->sort(compareHits);
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
    std::cout << "All sequences processed.\n\n";
    std::cout << alignmentsNum << " alignments calculated.\n";
    std::cout << passedNum << " sequences passed the thresholds (" << ((float)passedNum/(float)alignmentsNum) << " of overall calculated).\n";
    std::cout << ((float)passedNum/(float)prefdbr->getSize()) << " hits per query sequence.\n";

    seqdbr->close();
    prefdbr->close();
    dbw->close();

    return 0;
}


