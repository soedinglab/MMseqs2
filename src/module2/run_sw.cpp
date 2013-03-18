#include "run_sw.h"

void printUsage(){

    std::string usage("\nCalculates Smith-Waterman scores with SSE2 Smith Waterman between all sequence pairs contained in prefiltering score lists.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: run_sw dbBase prefilteringScoresBase scoringMatrixFile outBase [opts]\n" 
            "-e\t[double]\tE-value threshold.\n" 
            "-c\t[double]\tMinimum coverage of the longest sequence.\n");
    std::cout << usage;
}

void parseArgs(int argc, char** argv, std::string* ffindexSWBase, std::string* ffindexPrefBase, std::string* mFile, std::string* ffindexDbBase, double* evalThr, double* covThr){
    if (argc < 5){
        printUsage();
        exit(EXIT_FAILURE);
    }

    ffindexDbBase->assign(argv[1]);
    ffindexPrefBase->assign(argv[2]);
    mFile->assign(argv[3]);
    ffindexSWBase->assign(argv[4]);
    int i = 5;
    while (i < argc){
        if (strcmp(argv[i], "-e") == 0){
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
        else {
            printUsage();
            std::cerr << "Wrong argument: " << argv[i] << "\n";
            exit(EXIT_FAILURE);
        }
    }
}

bool compareHits (hit_t first, hit_t second){
    if (first.score > second.score)
        return true;
    return false;
}

void initFFIndexRead(std::string ffindexBase, char** data, ffindex_index_t** index, size_t* data_size){
    std::string dataFileName = ffindexBase;
    std::string indexFileName = ffindexBase + ".index";

    FILE* dataFile = fopen(dataFileName.c_str(), "r");
    FILE* indexFile = fopen(indexFileName.c_str(), "r");

    if( indexFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", indexFileName.c_str());  exit(EXIT_FAILURE); }
    if( dataFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", dataFileName.c_str());  exit(EXIT_FAILURE); }

    *data = ffindex_mmap_data(dataFile, data_size);

    *index = ffindex_index_parse(indexFile, 0);
    if(index == NULL)
    {   
        fferror_print(__FILE__, __LINE__, "ffindex_index_parse", (char*)indexFileName.c_str());
        exit(EXIT_FAILURE);
    }
}

void initFFIndexWrite(std::string ffindexBase, FILE** dataFile, FILE** indexFile){
    std::string dataFileName = ffindexBase;
    std::string indexFileName = ffindexBase + ".index";

    struct stat st;
    if(stat(dataFileName.c_str(), &st) == 0) { errno = EEXIST; perror(dataFileName.c_str()); exit(EXIT_FAILURE); }
    if(stat(indexFileName.c_str(), &st) == 0) { errno = EEXIST; perror(indexFileName.c_str()); exit(EXIT_FAILURE); }

    *dataFile = fopen(dataFileName.c_str(), "w");
    *indexFile = fopen(indexFileName.c_str(), "w");

    if( *indexFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", indexFileName.c_str());  exit(EXIT_FAILURE); }
    if( *dataFile == NULL) { fferror_print(__FILE__, __LINE__, "run_sw", dataFileName.c_str());  exit(EXIT_FAILURE); }
}

sequence_t seq2sequence_t(char* seq, size_t len, int* aa2int){
    // count spaces
    int spaces = 0;
    for (int i = 0; i < len; i++){
        if (seq[i] == '\n' || seq[i] == '\0'){
            spaces++;
        }
    }
    // convert input into a unsigned char array (spaces and \0 omitted)
    unsigned char* seq_uc = new unsigned char[len-spaces];
    int pos = 0;
    for (int i = 0; i < len; i++){
        if (seq[i] == '\n' || seq[i] == '\0')
            continue;
        seq_uc[pos] = (unsigned char)aa2int[seq[i]];
        pos++;
    }
    int newLen = len - spaces;

    sequence_t ret = {newLen, seq_uc};
    return ret;
}

std::list<hit_t>* getSWScoresForSequence(char* querySeqData, int querySeqLen, char* prefList, char* dbData, ffindex_index_t* dbIndex, void* workspace, int* aa2int, int** scMatrix, unsigned char bias, int dbSize)
{
    // convert the sequence into the numeric encoding (in the unsigned char range) and get the query profile
    sequence_t querySeq = seq2sequence_t(querySeqData, querySeqLen, aa2int);
    unsigned char* queryProfByte = getQueryProfileByte(querySeq.sequence, querySeq.length, scMatrix, bias);
    unsigned short* queryProfWord = getQueryProfileWord(querySeq.sequence, querySeq.length, scMatrix);

    // parse the prefiltering scores list
    char id[20];
    char scoreString[20];
    double score;
    int j = 0;
    std::list<hit_t>* hitList = new std::list<hit_t>();

    while(sscanf(prefList + j, "%s %s\n", id, scoreString) != EOF){
        score = atof(scoreString);

        // get the database sequence data
        ffindex_entry_t* dbSeqEntry = ffindex_bsearch_get_entry(dbIndex, id);
        char* dbSeqData = ffindex_get_data_by_entry(dbData, dbSeqEntry);
        sequence_t dbSeq = seq2sequence_t(dbSeqData, dbSeqEntry->length, aa2int);

        void* Hmatrix = memalign(16,(querySeq.length + 7)/8 * 8 * 8 * dbSeq.length * 16);   // 2GB für 36805*36805 (Q3ASY8_CHLCH)
        void* Ematrix = memalign(16,(querySeq.length + 7)/8 * 8 * 8 * dbSeq.length * 16);   
        void* Fmatrix = memalign(16,(querySeq.length + 7)/8 * 8 * 8 * dbSeq.length * 16);

        std::cout << "\t" << std::string(id) << "\n";
        // run SW alignment for the query sequence and the prefiltering list sequence
        //      int sw_score = smith_waterman_sse2_byte(querySeq.sequence, queryProfByte, querySeq.length, dbSeq.sequence, dbSeq.length, bias, GAP_OPEN, GAP_EXTEND, workspace);
        //      if (sw_score >= 255)
        aln_t aln = {0, 0, 0, 0};
        std::cout << "calculating score\n";
        int sw_score = smith_waterman_sse2_word(querySeq.sequence, queryProfWord, querySeq.length, dbSeq.sequence, dbSeq.length, GAP_OPEN, GAP_EXTEND, workspace, Hmatrix, Ematrix, Fmatrix, &aln.qEndPos, &aln.dbEndPos);
        std::cout << "score: " << sw_score << "\n";
        traceback_word((short*) Hmatrix, (short*) Ematrix, (short*) Fmatrix, querySeq.sequence, queryProfWord, querySeq.length, dbSeq.sequence, dbSeq.length, aln.qEndPos, aln.dbEndPos, GAP_OPEN, GAP_EXTEND, &aln.qStartPos, &aln.dbStartPos);

        float qcov = (aln.qEndPos - aln.qStartPos)/(float)querySeq.length;
        float dbcov = (aln.dbEndPos - aln.dbStartPos)/(float)dbSeq.length;
        double evalue = (double)(dbSize * querySeq.length * dbSeq.length) * fpow2((double)-sw_score/BIT_FACTOR);

        /*std::cout << "score: " << sw_score << "\n";
          std::cout << "Query coverage: " << qcov << "\n";
          std::cout << "DB coverage: " << dbcov << "\n";
          std::cout << "E-value: " << evalue << "\n\n";
          */
        delete[] dbSeq.sequence;

        if (sw_score > SW_SCORE_THR){
            hit_t h = {std::string(id), sw_score, qcov, dbcov, evalue};
            hitList->push_back(h);
        }

        j += strlen(id) + strlen(scoreString) + 2;
        free(Hmatrix);
        free(Ematrix);
        free(Fmatrix);

    }
    hitList->sort(compareHits);

    delete[] querySeq.sequence;
    free(queryProfByte);
    free(queryProfWord);

    return hitList;
}

void runSWParallel(std::string ffindexPrefBase, std::string ffindexDbBase, std::string ffindexSWBase, int* aa2int, int** scMatrix, double evalThr, double covThr){

    // set scoring bias
    unsigned char bias = getMatrixMinValue(scMatrix);

    // open prefiltering scores ffindex 
    char* prefData;
    ffindex_index_t* prefIndex;

    size_t prefDBSize;
    initFFIndexRead(ffindexPrefBase, &prefData, &prefIndex, &prefDBSize);

    // open sequence database ffindex
    char* dbData;
    ffindex_index_t* dbIndex;

    size_t dbSize;
    initFFIndexRead(ffindexDbBase, &dbData, &dbIndex, &dbSize);

    // allocate aligned memory for each thread 
    int threads = 1;
#ifdef OPENMP
    threads = omp_get_max_threads();
#endif
    void ** workspace = new void*[threads];
# pragma omp parallel for schedule(static)
    for (int i = 0; i < threads; i++){
        // the longest sequence in UniProt is 36805 amino acids
        // 1.2 MB per thread
        void * workspace_memory  = (void *)memalign(16,2*16*40000+256);
        workspace[i] = (void *) ((((size_t) workspace_memory) + 255) & (~0xff));
    }

    // generate output ffindex databases for each thread
    FILE* outDataFile[threads];
    FILE* outIndexFile[threads];
    for (int i = 0; i < threads; i++){
        std::stringstream ss;
        ss << i;
        std::string outBase = ffindexSWBase + "." + ss.str();
        initFFIndexWrite(outBase, &outDataFile[i], &outIndexFile[i]);
    }

    size_t* offset = new size_t[threads];
    for (int i = 0; i < threads; i++)
        offset[i] = 0;

    char** outBuffer = new char*[threads];
    for (int i = 0; i < threads; i++)
        outBuffer[i] = new char[1000];

    // for loop over all the sequences in the prefiltering
    # pragma omp parallel for schedule(static)
    for(int i = 0; i < prefIndex->n_entries; i++){
        int thr_idx = 0;
#ifdef OPENMP
        thr_idx = omp_get_thread_num();
#endif

        // get prefiltering results for the query sequence
        ffindex_entry_t* prefEntry = ffindex_get_entry_by_index(prefIndex, i);
        if(prefEntry == NULL) { perror(prefEntry->name); continue; }
        
        //std::cout << prefEntry->name << " -> ";
        // get the query sequence
        char* prefList = ffindex_get_data_by_entry(prefData, prefEntry);
        std::cout << prefEntry->name << "\n";
        ffindex_entry_t* querySeqEntry = ffindex_bsearch_get_entry(dbIndex, prefEntry->name);
        char* querySeqData =  ffindex_get_data_by_entry(dbData, querySeqEntry);

        // calculate SW scores for all sequences in the prefiltering list
        std::list<hit_t>* hitList = getSWScoresForSequence(querySeqData, querySeqEntry->length, prefList, dbData, dbIndex, workspace[thr_idx], aa2int, scMatrix, bias, dbSize);

        // write the hits into ffindex
        std::list<hit_t>::iterator it;
        std::stringstream hitListStringStream;

        // put the contents of hitList into a string
        for (it = hitList->begin(); it != hitList->end(); ++it){
            if (it->eval <= evalThr && it->qcov >= covThr && it->dbcov >= covThr){
                hitListStringStream << it->name << " " << it->score << " " << it->qcov << " " << it->dbcov << " " << it->eval << "\n";
            }
        }

        std::string hitListString = hitListStringStream.str();
        if (hitListString.size() > 0){
            const char* hitListStringData = hitListString.c_str();
            memcpy(outBuffer[thr_idx], hitListStringData, hitListString.length()*sizeof(char));

            ffindex_insert_memory(outDataFile[thr_idx], outIndexFile[thr_idx], &offset[thr_idx], outBuffer[thr_idx], hitListString.length(), prefEntry->name);
        }
        
        delete hitList;
    } // end parallel for
    
    // merge ffindexes from each thread into one ffindex
    char* merge_command  = (char*) malloc(FILENAME_MAX * 5);
    for(int j = 0; j < threads; j++)
    {
        fclose(outDataFile[j]);
        fclose(outIndexFile[j]);
        std::string indexSuffix("index");
        snprintf(merge_command, FILENAME_MAX, "ffindex_build -as %s %s.%s -d %s.%d -i %s.%d.%s",
                ffindexSWBase.c_str(), ffindexSWBase.c_str(), indexSuffix.c_str(), ffindexSWBase.c_str(), j, ffindexSWBase.c_str(), j, indexSuffix.c_str());
        system(merge_command);
    }

}

int main(int argc, char **argv)
{

    std::string ffindexSWBase = "";
    std::string ffindexPrefBase = "";
    std::string mFile = "";
    std::string ffindexDbBase = "";
    double evalThr = 0.01;
    double covThr = 0.8;

    parseArgs(argc, argv, &ffindexSWBase, &ffindexPrefBase, &mFile, &ffindexDbBase, &evalThr, &covThr);

    int* aa2int = get_aa2int(get_int2aa());
    int** scMatrix = getBiasedScoringMatrix(mFile, aa2int);

    runSWParallel(ffindexPrefBase, ffindexDbBase, ffindexSWBase, aa2int, scMatrix, evalThr, covThr);

}
