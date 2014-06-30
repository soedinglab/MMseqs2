#include "WorkflowFunctions.h"

ffindex_index_t* openIndex(const char* indexFileName){
    // count the number of entries in the clustering
    char line [1000];
    int cnt = 0;
    std::ifstream index_file(indexFileName);
    if (index_file.is_open()) {
        while ( index_file.getline (line, 1000) ){
            cnt++;
        }
        index_file.close();
    }
    else{
        std::cerr << "Could not open ffindex index file " << indexFileName << "\n";
        exit(EXIT_FAILURE);
    }
    // open clustering ffindex
    FILE* indexFile = fopen(indexFileName, "r");
    if( indexFile == NULL) { fferror_print(__FILE__, __LINE__, "DBReader", indexFileName);  exit(EXIT_FAILURE); }

    ffindex_index_t* index = ffindex_index_parse(indexFile, cnt);
    return index;
}

std::string runStep(std::string inDBData, std::string inDBWorkingIndex, std::string targetDBData, std::string targetDBIndex, std::string tmpDir,
        std::string scoringMatrixFile, int maxSeqLen, int seqType,
        int kmerSize, int alphabetSize, size_t maxResListLen, int split, int skip, bool aaBiasCorrection, float zscoreThr, float sensitivity,
        double evalThr, double covThr, int maxRejects,
        int step_num, int restart, bool search){

    std::cout << "------------------------------------------------------------\n";
    std::cout << "GENERAL PARAMETERS:\n";
    std::cout << "max. sequence length: " << maxSeqLen << "\n";
    std::cout << "max. result list length: " << maxResListLen << "\n";
    std::cout << "PREFILTERING PARAMETERS:\n";
    std::cout << "k-mer size: " << kmerSize << "\n";
    std::cout << "z-score threshold: " << zscoreThr << "\n";
    std::cout << "sensitivity: " << sensitivity << "\n";
    std::cout << "ALIGNMENT PARAMETERS:\n";
    std::cout << "e-value threshold: " << evalThr << "\n";
    std::cout << "coverage threshold: " << covThr << "\n";
    std::cout << "Maximum number of rejects: " << maxRejects << "\n";
    std::cout << "------------------------------------------------------------\n\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // prefiltering step
    std::stringstream ss;
    ss << step_num;
    std::string prefDB_step = tmpDir + "/db_pref_step" + ss.str();

    int sec;
    if (restart <= 1){
        std::cout << "\n##### Step " << ss.str() << ": PREFILTERING #####\n\n";
        Prefiltering* pref = new Prefiltering (inDBData, inDBWorkingIndex,
                inDBData, inDBWorkingIndex,
                prefDB_step, prefDB_step+ ".index",
                scoringMatrixFile, sensitivity, kmerSize, alphabetSize, zscoreThr, maxSeqLen, seqType, aaBiasCorrection, split, skip);
        std::cout << "Starting prefiltering scores calculation.\n";
        pref->run(maxResListLen);
        delete pref;

        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        std::cout << "\nTime for step " << step_num << " prefiltering: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);
    }

    // alignment step
    std::string alnDB_step = tmpDir + "/db_aln_step" + ss.str();

    if (restart <= 2){
        std::cout << "\n##### Step " << ss.str() << ": ALIGNMENT #####\n\n";
        Alignment* aln = new Alignment(inDBData, inDBWorkingIndex,
                inDBData, inDBWorkingIndex,
                prefDB_step, prefDB_step + ".index",
                alnDB_step, alnDB_step + ".index",
                scoringMatrixFile, evalThr, covThr, maxSeqLen, seqType);
        std::cout << "Starting alignments calculation.\n";
        aln->run(maxResListLen, maxRejects);
        delete aln;

        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        std::cout << "\nTime for step " << step_num << " alignment: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);
    }

    if (search)
        return alnDB_step;

    // clustering step
    std::string cluDB_step = tmpDir + "/db_clu_step" + ss.str();

    if (restart <= 3){
        std::cout << "\n##### Step " << ss.str() << ": CLUSTERING #####\n\n";
        Clustering* clu = new Clustering(inDBData, inDBWorkingIndex,
                alnDB_step, alnDB_step + ".index",
                cluDB_step, cluDB_step + ".index",
                0.0, 0, maxResListLen);
        clu->run(Clustering::SET_COVER);
        delete clu;

        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        std::cout << "\nTime for step " << step_num << " clustering: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);
    }

    return cluDB_step;
}

void copy(std::string inFile, std::string outFile){
    // buffer
    std::ifstream inf(inFile.c_str(), std::ios::binary);
    std::ofstream outf(outFile.c_str(), std::ios::binary);
    if (inf.is_open()) {
        outf << inf.rdbuf();
        inf.close();
        outf.close();
    }
    else{
        std::cerr << "Could not open file " << inFile << "\n";
        exit(EXIT_FAILURE);
    }
}

float getZscoreForSensitivity (float sensitivity){
    float zscoreThr = 50.0;

    if (1.0 <= sensitivity && sensitivity < 2.0)
        zscoreThr = 500.0;
    else if (2.0 <= sensitivity && sensitivity < 3.0)
        zscoreThr = 200.0;
    else if (3.0 <= sensitivity && sensitivity < 4.0)
        zscoreThr = 100.0;
    else if (4.0 <= sensitivity && sensitivity < 5.0)
        zscoreThr = 50.0;
    else if (5.0 <= sensitivity && sensitivity < 6.0)
        zscoreThr = 40.0;
    else if (6.0 <= sensitivity && sensitivity < 7.0)
        zscoreThr = 30.0;
    else if (7.0 <= sensitivity && sensitivity < 8.0)
        zscoreThr = 20.0;
    else if (8.0 <= sensitivity && sensitivity <= 9.0)
        zscoreThr = 10.0;

    return zscoreThr;
}

