#include <iostream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <list>

#include "WorkflowFunctions.h"

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

void printUsage(){

    std::string usage("\nCalculates the clustering of the sequences in the input database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de))\n\n");
    usage.append("USAGE: mmseqs_clustering [sequenceDBBase] [outDBBase] [tmpDir] [opts]\n"
            "GENERAL OPTIONS:\n"
            "--cascaded      \t\tStart the cascaded instead of simple clustering workflow.\n"
            "-s              \t[float]\tTarget sensitivity in the range [2:9] (default=4).\n"
            "--max-seqs      \tMaximum result sequences per query (default=300)\n"
            "--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n"
//            "--restart          \t[int]\tRestart the clustering workflow starting with alignment or clustering.\n"
//            "                \t     \tThe value is in the range [1:3]: 1: restart from prefiltering  2: from alignment; 3: from clustering.\n"
//            "CASCADED CLUSTERING OPTIONS:\n"
            "--sub-mat       \t[file]\tAmino acid substitution matrix file.\n"
/*            "--step          \t[int]\tRestart the step of the cascaded clustering. For values in [1:3], the resprective step number, 4 is only the database merging.\n"
            "\nRESTART OPTIONS NOTE:\n"
            "                \t     \tIt is assumed that all valid intermediate results exist.\n"
            "                \t     \tValid intermediate results are taken from the tmp directory specified by the user.\n"
            "                \t     \tFor the cascaded clustering, --restart and --step options can be combined.\n"*/
            );
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexInDBBase, std::string* ffindexOutDBBase, std::string* tmpDir, std::string* scoringMatrixFile, size_t* maxSeqLen, bool* cascaded, float* sens, size_t* maxResListLen, int* restart, int* step){
    if (argc < 4){
        printUsage();
        exit(EXIT_FAILURE);
    }

    ffindexInDBBase->assign(argv[1]);
    ffindexOutDBBase->assign(argv[2]);
    tmpDir->assign(argv[3]);

    int i = 4;
    while (i < argc){

        if (strcmp(argv[i], "--restart") == 0){
            if (++i < argc){
                *restart = atoi(argv[i]);
                i++;
            }
            else{
                printUsage();
                std::cerr << "No value provided for " << argv[i] << "\n";
                exit(EXIT_FAILURE);

            }   
        }
        else if (strcmp(argv[i], "--step") == 0){
            if (++i < argc){
                *step = atoi(argv[i]);
                i++;
            }
            else{
                printUsage();
                std::cerr << "No value provided for " << argv[i] << "\n";
                exit(EXIT_FAILURE);

            }   
        }
        else if (strcmp(argv[i], "-m") == 0){
            if (++i < argc){
                scoringMatrixFile->assign(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--sub-mat") == 0){
            if (++i < argc){
                scoringMatrixFile->assign(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--max-seq-len") == 0){
            if (++i < argc){
                *maxSeqLen = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for " << argv[i] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "--cascaded") == 0){
            *cascaded = true;
            i++;
        }
        else if (strcmp(argv[i], "-s") == 0){
            if (++i < argc){
                *sens = atof(argv[i]);
                if (*sens < 2.0 || *sens > 9.0){
                    Debug(Debug::ERROR) << "Please choose sensitivity in the range [2:9].\n";
                    exit(EXIT_FAILURE);
                }
                i++;
            }
        }
        else if (strcmp(argv[i], "--max-seqs") == 0){
            if (++i < argc){
                *maxResListLen = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                exit(EXIT_FAILURE);
            }
        }
        else {
            printUsage();
            std::cerr << "Wrong argument: " << argv[i] << "\n";
            exit(EXIT_FAILURE);
        }
    }
    if (strcmp (scoringMatrixFile->c_str(), "") == 0){
        printUsage();
        std::cerr << "\nPlease provide a scoring matrix file. You can find scoring matrix files in $INSTALLDIR/data/.\n";
        exit(EXIT_FAILURE);
    }
    if ((!*cascaded) && *step > 1){
        printUsage();
        std::cout << "Only cascaded clustering has clustering steps. Please set clustering mode to cascaded or leave the step parameter out.\n";
        exit(EXIT_FAILURE);
    }
}

void extractNewIndex(std::string seqDBIndex, std::string cluDBIndex, std::string newIndexFileName){

    ffindex_index_t* seq_index = openIndex(seqDBIndex.c_str());

    ffindex_index_t* clu_index = openIndex(cluDBIndex.c_str());

    FILE* new_index_file = fopen(newIndexFileName.c_str(), "w");

    for (unsigned int i = 0; i < clu_index->n_entries; i++){
        // get the key in the clustering index
        char* dbKey = &(ffindex_get_entry_by_index(clu_index, i)->name[0]);
        // get the entry from the sequence index
        ffindex_entry_t* e = ffindex_get_entry_by_name(seq_index, dbKey);
        // write the sequence index entry into the new index
        fprintf(new_index_file, "%s\t%zd\t%zd\n", e->name, e->offset, e->length);
    }

    // close all the files
    fclose(seq_index->file);
    fclose(clu_index->file);
    fclose(new_index_file);

}

void mergeClusteringResults(std::string seqDB, std::string outDB, std::list<std::string> cluSteps){

    // open the sequence database
    // it will serve as the reference for sequence indexes
    std::string seqDBIndex = seqDB + ".index";
    DBReader* dbr = new DBReader(seqDB.c_str(), seqDBIndex.c_str());
    dbr->open(DBReader::NOSORT);

    // init the structure for cluster merging
    std::list<int>** mergedClustering = new std::list<int>*[dbr->getSize()];
    for (unsigned int i = 0; i < dbr->getSize(); i++){
        mergedClustering[i] = new std::list<int>();
    }

    // read the clustering from the first clustering step
    std::string firstCluStep = cluSteps.front();
    cluSteps.pop_front();
    std::string firstCluStepIndex = firstCluStep + ".index";
    DBReader* cluStepDbr = new DBReader(firstCluStep.c_str(), firstCluStepIndex.c_str());
    cluStepDbr->open(DBReader::NOSORT);

    for (unsigned int i = 0; i < cluStepDbr->getSize(); i++){
        int cluId = dbr->getId(cluStepDbr->getDbKey(i));
        std::stringstream lineSs (cluStepDbr->getData(i));
        std::string val;
        // go through the sequences in the cluster and add them to the initial clustering
        while (std::getline(lineSs, val)){
            int seqId = dbr->getId(val.c_str()); 
            mergedClustering[cluId]->push_back(seqId);
        }
    }
    cluStepDbr->close();

    std::cout << "Clustering step 1...\n";

    // merge later clustering steps into the initial clustering step
    int cnt = 2;
    while(!cluSteps.empty()){
        // open the next clustering database
        std::string cluStep = cluSteps.front();
        std::string cluStepIndex = cluStep + ".index";
        cluSteps.pop_front();

        cluStepDbr = new DBReader(cluStep.c_str(), cluStepIndex.c_str());
        cluStepDbr->open(DBReader::NOSORT);

        // go through the clusters and merge them into the clusters from the previous clustering step
        for (unsigned int i = 0; i < cluStepDbr->getSize(); i++){
            int cluId = dbr->getId(cluStepDbr->getDbKey(i));
            char* cluData = cluStepDbr->getData(i);
            std::stringstream lineSs(cluData);
            std::string val;
            // go through the sequences in the cluster and add them and their clusters to the cluster of cluId
            // afterwards, delete the added cluster from the clustering
            while (std::getline(lineSs, val, '\n')){
                int seqId = dbr->getId(val.c_str());
                mergedClustering[cluId]->splice(mergedClustering[cluId]->end(), *mergedClustering[seqId]);
            }
        }
        cluStepDbr->close();

        std::cout << "Clustering step " << cnt << "...\n";
        cnt++;
    }

    std::cout << "Writing the results...\n";

    std::string outDBIndex = outDB + ".index";
    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str());
    dbw->open();

    size_t BUFFER_SIZE = 1000000;
    char* outBuffer = new char[BUFFER_SIZE];
    // go through all sequences in the database
    for (unsigned int i = 0; i < dbr->getSize(); i++){

        // no cluster for this representative
        if (mergedClustering[i]->size() == 0)
            continue;

        // representative
        char* dbKey = dbr->getDbKey(i);

        std::stringstream res;
        for(std::list<int>::iterator it = mergedClustering[i]->begin(); it != mergedClustering[i]->end(); ++it){
            res << dbr->getDbKey(*it) << "\n";
        }

        std::string cluResultsOutString = res.str();
        const char* cluResultsOutData = cluResultsOutString.c_str();
        if (BUFFER_SIZE < strlen(cluResultsOutData)){
            std::cerr << "Tried to process the clustering list for the query " << dbKey << " , the length of the list = " << mergedClustering[i]->size() << "\n";
            std::cerr << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " << cluResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters -> output buffer is already huge ;-)\n";
            continue;
        }
        memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length()*sizeof(char));
        dbw->write(outBuffer, cluResultsOutString.length(), dbKey);
    }
    dbw->close();
    delete[] outBuffer;

    // delete the clustering data structure
    for (unsigned int i = 0; i < dbr->getSize(); i++){
        delete mergedClustering[i];
    }
    delete[] mergedClustering;
    std::cout << "...done.\n";

}

void runClustering(float sensitivity, size_t maxSeqLen, int seqType, 
        int kmerSize, int alphabetSize, size_t maxResListLen, int split, int skip, bool aaBiasCorrection, 
        double evalThr, double covThr, 
        std::string inDB, std::string outDB, std::string scoringMatrixFile, std::string tmpDir, int restart){

    std::string inDBIndex = inDB + ".index";

    float zscoreThr = getZscoreForSensitivity(sensitivity);

    std::cout << "\nRunning the clustering with sensitivity " << sensitivity << "\n";
    std::cout << "Z-score threshold: " << zscoreThr << "\n";

    std::string cluDB = runStep(inDB, inDBIndex, inDB, inDBIndex, tmpDir, scoringMatrixFile, maxSeqLen, seqType, kmerSize, alphabetSize, maxResListLen, split, skip, aaBiasCorrection, zscoreThr, sensitivity, evalThr, covThr, 10, 1, restart, false);

    std::string cluDBIndex = cluDB + ".index";
    std::string outDBIndex = outDB + ".index";

    // copy the clustering databases to the right location
    copy(cluDBIndex, outDBIndex);
    copy(cluDB, outDB); 

}

void runSearch(float sensitivity, size_t maxSeqLen, int seqType,
        int kmerSize, int alphabetSize, size_t maxResListLen, int split, int skip, bool aaBiasCorrection,
        double evalThr, double covThr,
        std::string inDB, std::string targetDB, std::string outDB, std::string scoringMatrixFile, std::string tmpDir, int restart){

    std::string inDBIndex = inDB + ".index";
    std::string targetDBIndex = targetDB + ".index";
    float zscoreThr = getZscoreForSensitivity(sensitivity);

    std::string alnDB = runStep(inDB, inDBIndex, targetDB, targetDBIndex, tmpDir, scoringMatrixFile, maxSeqLen, seqType, kmerSize, alphabetSize, maxResListLen, split, skip, aaBiasCorrection, zscoreThr, sensitivity, evalThr, covThr, INT_MAX, 1, restart, true);

    std::string alnDBIndex = alnDB + ".index";
    std::string outDBIndex = outDB + ".index";

    // copy the clustering databases to the right location
    copy(alnDBIndex, outDBIndex);
    copy(alnDB, outDB);
}

void runCascadedClustering(float targetSensitivity, size_t maxSeqLen, int seqType,
        int kmerSize, int alphabetSize, size_t maxResListLen, int split, int skip, bool aaBiasCorrection,
        double evalThr, double covThr,
        std::string inDB, std::string outDB, std::string scoringMatrixFile, std::string tmpDir, int restart, int step){

    std::cout << "\nRunning cascaded clustering for the database " << inDB << "\n";
    std::cout << "Target sensitivity: " << targetSensitivity << "\n";

    std::string inDBIndex = inDB + ".index";

    // copy index to a new location, it will be overwritten
    std::string inDBWorkingIndex = tmpDir + "/input_seqs.index";

    // save the original sequences to the working location
    // index will be overwritten in each clustering step
    ffindex_index_t* index_orig = openIndex(inDBIndex.c_str());
    FILE* index_cpy_file = fopen(inDBWorkingIndex.c_str(), "w");
    ffindex_write(index_orig, index_cpy_file);
    fclose(index_cpy_file);

    std::list<std::string> cluSteps;
    std::string cluDB = "";

    std::cout << "\n\n";
    std::cout << "------------------------------- Step 1 ----------------------------------------------\n";

    float sens = 1.0;
    float zscoreThr = getZscoreForSensitivity(sens);

    int local_restart;
    if (step > 1)
        // skip step
        local_restart = 4;
    else if (step == 1)
        // set the starting point of the clustering step
        local_restart = restart;
    else
        local_restart = 1;
    cluDB = runStep(inDB, inDBWorkingIndex, inDB, inDBWorkingIndex, tmpDir, scoringMatrixFile, maxSeqLen, seqType, kmerSize, alphabetSize, 50, split, skip, aaBiasCorrection, zscoreThr, sens, evalThr, covThr, 10, 1, local_restart, false);
    cluSteps.push_back(cluDB);

    std::cout << "\n##### Updating databases #####\n\n";
    // extract new index containing only sequences remained after the clustering and overwrite old index with the extracted index
    std::string cluDBIndex = cluDB + ".index";
    std::string newIndexFileName =  tmpDir + "/seqs_tmp.index";
    extractNewIndex(inDBWorkingIndex, cluDBIndex, newIndexFileName);
    copy(newIndexFileName, inDBWorkingIndex);

    std::cout << "------------------------------- Step 2 ----------------------------------------------\n";
    sens = targetSensitivity/2.0;
    zscoreThr = getZscoreForSensitivity(sens);
    
    if (step > 2)
        // skip step
        local_restart = 4;
    else if (step == 2)
        // set the starting point of the clustering step
        local_restart = restart;
    else
        local_restart = 1;

    cluDB = runStep(inDB, inDBWorkingIndex, inDB, inDBWorkingIndex, tmpDir, scoringMatrixFile, maxSeqLen, seqType, kmerSize, alphabetSize, maxResListLen, split, skip, aaBiasCorrection, zscoreThr, sens, evalThr, covThr, 10, 2, local_restart, false);
    cluSteps.push_back(cluDB);

    std::cout << "\n##### Updating databases #####\n\n";
    // extract new index containing only sequences remained after the clustering and overwrite old index with the extracted index
    cluDBIndex = cluDB + ".index";
    extractNewIndex(inDBWorkingIndex, cluDBIndex, newIndexFileName);
    copy(newIndexFileName, inDBWorkingIndex);

    std::cout << "------------------------------- Step 3 ----------------------------------------------\n";
    sens = targetSensitivity;
    zscoreThr = getZscoreForSensitivity(sens);
    
    if (step > 3)
        // skip step
        local_restart = 4;
    else if (step == 3)
        // set the starting point of the clustering step
        local_restart = restart;
    else
        local_restart = 1;

    cluDB = runStep(inDB, inDBWorkingIndex, inDB, inDBWorkingIndex, tmpDir, scoringMatrixFile, maxSeqLen, seqType, kmerSize, alphabetSize, maxResListLen, split, skip, aaBiasCorrection, zscoreThr, sens, evalThr, covThr, INT_MAX, 3, local_restart, false);
    cluSteps.push_back(cluDB);

    std::cout << "--------------------------- Merging databases ---------------------------------------\n";
    mergeClusteringResults(inDB, outDB, cluSteps);

}

int main (int argc, const char * argv[]){

    // general parameters
    bool cascaded = false;
    size_t maxSeqLen = 50000;
    int seqType = Sequence::AMINO_ACIDS;
    float targetSens = 4.0;
    int restart = 0;
    int step = 1;
    size_t maxResListLen = 300;

    // parameter for the prefiltering
    int kmerSize = 6;
    int alphabetSize = 21;
    int split = 0;
    int skip = 0;
    bool aaBiasCorrection = true;

    // parameters for the alignment
    double evalThr = 0.001;
    double covThr = 0.8;

    std::string inDB = "";
    std::string outDB = "";
    std::string tmpDir = "";

    // get the path of the scoring matrix
    char* mmdir = getenv ("MMDIR");
    if (mmdir == 0){
        std::cerr << "Please set the environment variable $MMDIR to your MMSEQS installation directory.\n";
        exit(1);
    }
    std::string scoringMatrixFile(mmdir);
    scoringMatrixFile = scoringMatrixFile + "/data/blosum62.out";

    parseArgs(argc, argv, &inDB, &outDB, &tmpDir, &scoringMatrixFile, &maxSeqLen, &cascaded, &targetSens, &maxResListLen, &restart, &step);

    if (cascaded)
        runCascadedClustering(targetSens, maxSeqLen, seqType,
                kmerSize, alphabetSize, maxResListLen, split, skip, aaBiasCorrection, 
                evalThr, covThr,
                inDB, outDB, scoringMatrixFile, tmpDir, restart, step);
    else
        runClustering(targetSens, maxSeqLen, seqType,
                kmerSize, alphabetSize, maxResListLen, split, skip, aaBiasCorrection,
                evalThr, covThr,
                inDB, outDB, scoringMatrixFile, tmpDir, restart);

}
