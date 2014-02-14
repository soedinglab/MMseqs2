#include <iostream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <list>

#include "../prefiltering/Prefiltering.h"
#include "../alignment/Alignment.h"
#include "../clustering/Clustering.h"

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

void printUsage(){

    std::string usage("\nCalculates similarity scores between all sequences in the query database and all sequences in the target database using cascaded clustering algorithm.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de))\n\n");
    usage.append("USAGE: cascaded_clustering ffindexInDBBase ffindexOutDBBase tmpDir [opts]\n"
            "-m              \t[file]\tAmino acid substitution matrix file.\n"
            "--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n");
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexInDBBase, std::string* ffindexOutDBBase, std::string* tmpDir, std::string* scoringMatrixFile, size_t* maxSeqLen){
    if (argc < 4){
        printUsage();
        exit(EXIT_FAILURE);
    }

    ffindexInDBBase->assign(argv[1]);
    ffindexOutDBBase->assign(argv[2]);
    tmpDir->assign(argv[3]);

    int i = 4;
    while (i < argc){
        if (strcmp(argv[i], "-m") == 0){
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
}

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

void extractNewIndex(std::string seqDBIndex, std::string cluDBIndex, std::string newIndexFileName){

    ffindex_index_t* seq_index = openIndex(seqDBIndex.c_str());

    ffindex_index_t* clu_index = openIndex(cluDBIndex.c_str());

    FILE* new_index_file = fopen(newIndexFileName.c_str(), "w");

    for (int i = 0; i < clu_index->n_entries; i++){
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

std::string runStep(std::string inDBData, std::string inDBWorkingIndex, std::string tmpDir, 
        std::string scoringMatrixFile, int maxSeqLen, int seqType, 
        int kmerSize, int alphabetSize, size_t maxResListLen, int skip, bool aaBiasCorrection, float zscoreThr, float sensitivity, 
        double evalThr, double covThr, int maxAlnNum, int step_num){

    struct timeval start, end;
    gettimeofday(&start, NULL);

    // prefiltering step
    std::stringstream ss;
    ss << step_num;
    std::cout << "\n##### Step " << ss.str() << ": PREFILTERING #####\n\n";
    std::string prefDB_step = tmpDir + "/db_pref_step" + ss.str();
    Prefiltering* pref = new Prefiltering (inDBData, inDBWorkingIndex, 
            inDBData, inDBWorkingIndex, 
            prefDB_step, prefDB_step+ ".index", 
            scoringMatrixFile, sensitivity, kmerSize, alphabetSize, zscoreThr, maxSeqLen, seqType, aaBiasCorrection, skip);
    std::cout << "Starting prefiltering scores calculation.\n";
    pref->run(maxResListLen);
    delete pref;

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for step " << step_num << " prefiltering: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    // alignment step
    std::cout << "\n##### Step " << ss.str() << ": ALIGNMENT #####\n\n"; 
    std::string alnDB_step = tmpDir + "/db_aln_step" + ss.str();
    Alignment* aln = new Alignment(inDBData, inDBWorkingIndex, 
            inDBData, inDBWorkingIndex,
            prefDB_step, prefDB_step + ".index", 
            alnDB_step, alnDB_step + ".index", 
            scoringMatrixFile, evalThr, covThr, maxSeqLen, seqType);
    std::cout << "Starting alignments calculation.\n";
    aln->run(maxResListLen);
    delete aln;

    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for step " << step_num << " alignment: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    // clustering step
    std::cout << "\n##### Step " << ss.str() << ": CLUSTERING #####\n\n";
    std::string cluDB_step = tmpDir + "/db_clu_step" + ss.str();
    Clustering* clu = new Clustering(inDBData, inDBWorkingIndex,
            alnDB_step, alnDB_step + ".index",
            cluDB_step, cluDB_step + ".index",
            0.0);
    clu->run(Clustering::SET_COVER);
    delete clu;

    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for step " << step_num << " clustering: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    std::cout << "\n##### Updating databases #####\n\n";
    // extract new index containing only sequences remained after the clustering and overwrite old index with the extracted index
    std::string cluDBIndex = cluDB_step + ".index";
    std::string newIndexFileName =  tmpDir + "/db_seq_step" + ss.str();
    extractNewIndex(inDBWorkingIndex, cluDBIndex, inDBWorkingIndex);

/*    // Overwrite old index with the new index
    // (this is not the original index!)
    ffindex_index_t* seq_index = openIndex(newIndexFileName.c_str());

    FILE* index_file = fopen(inDBWorkingIndex.c_str(), "w");

    for (int i = 0; i < seq_index->n_entries; i++){
        // get the entry from the sequence index
        ffindex_entry_t* e = ffindex_get_entry_by_index(seq_index, i);
        // write the sequence index entry into the new index
        fprintf(index_file, "%s\t%zd\t%zd\n", e->name, e->offset, e->length);
    }

    // close all the files
    fclose(seq_index->file);
    fclose(index_file);*/

    return cluDB_step;
}

void mergeClusteringResults(std::string seqDB, std::string outDB, std::list<std::string> cluSteps){

    // open the sequence database
    // it will serve as the reference for sequence indexes
    std::string seqDBIndex = seqDB + ".index";
    DBReader* dbr = new DBReader(seqDB.c_str(), seqDBIndex.c_str());
    dbr->open(DBReader::NOSORT);

    // init the structure for cluster merging
    std::list<int>** mergedClustering = new std::list<int>*[dbr->getSize()];
    for (int i = 0; i < dbr->getSize(); i++){
        mergedClustering[i] = new std::list<int>();
    }

    // read the clustering from the first clustering step
    std::string firstCluStep = cluSteps.front();
    cluSteps.pop_front();
    std::string firstCluStepIndex = firstCluStep + ".index";
    DBReader* cluStepDbr = new DBReader(firstCluStep.c_str(), firstCluStepIndex.c_str());
    cluStepDbr->open(DBReader::NOSORT);

    for (int i = 0; i < cluStepDbr->getSize(); i++){
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
/*    // print out the merged clustering
    for (int i = 0; i < dbr->getSize(); i++){
        std::cout << dbr->getDbKey(i) << ": ";
        for(std::list<int>::iterator it = mergedClustering[i]->begin(); it != mergedClustering[i]->end(); ++it){
            std::cout << dbr->getDbKey(*it) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";*/

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
        for (int i = 0; i < cluStepDbr->getSize(); i++){
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
/*        // print out the merged clustering
        for (int i = 0; i < dbr->getSize(); i++){
            std::cout << dbr->getDbKey(i) << ": ";
            for(std::list<int>::iterator it = mergedClustering[i]->begin(); it != mergedClustering[i]->end(); ++it){
                std::cout << dbr->getDbKey(*it) << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";*/
        cnt++;
    }

    std::cout << "Writing the results...\n";
    
    std::string outDBIndex = outDB + ".index";
    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str());
    dbw->open();

    size_t BUFFER_SIZE = 1000000;
    char* outBuffer = new char[BUFFER_SIZE];
    // go through all sequences in the database
    for (int i = 0; i < dbr->getSize(); i++){

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
    for (int i = 0; i < dbr->getSize(); i++){
        delete mergedClustering[i];
    }
    delete[] mergedClustering;
    std::cout << "...done.\n";

}

int main (int argc, const char * argv[]){

    // general parameters
    size_t maxSeqLen = 50000;
    int seqType = Sequence::AMINO_ACIDS;
    
    // parameter for the prefiltering
    int kmerSize = 6;
    int alphabetSize = 21;
    size_t maxResListLen = 100;
    int skip = 0;
    bool aaBiasCorrection = false;
    float zscoreThr = 300.0f;

    // parameters for the alignment
    double evalThr = 0.001;
    double covThr = 0.8;
    int maxAlnNum = 10;

    std::string inDB = "";
    std::string outDB = "";
    std::string scoringMatrixFile = "";
    std::string tmpDir = "";

    parseArgs(argc, argv, &inDB, &outDB, &tmpDir, &scoringMatrixFile, &maxSeqLen);

    std::string inDBIndex = inDB + ".index";

    // copy index to a new location, it will be overwritten
    std::string inDBWorkingIndex = tmpDir + "/input_seqs.index";

    // save the original sequences to the working location
    // index will be overwritten in each clustering step
    ffindex_index_t* index_orig = openIndex(inDBIndex.c_str());
    FILE* index_cpy_file = fopen(inDBWorkingIndex.c_str(), "w");
    ffindex_write(index_orig, index_cpy_file);
    fclose(index_cpy_file);

    // save the data
/*    int BUF_SIZE = 100000;
    char buf[BUF_SIZE];
    std::ifstream data_file;
    data_file.open(inDB.c_str());
 
    std::ofstream data_cpy_file;
    data_cpy_file.open(inDBWorking.c_str(), std::ios::out);

   if (data_file.is_open()) {
        while (data_file.getline (buf, BUF_SIZE, '\0'))
            data_cpy_file << buf << '\0';
   }

   data_file.close();
   data_cpy_file.close();*/

   std::list<std::string> cluSteps;
   std::string cluDB = "";

   std::cout << "\n\n";
   std::cout << "------------------------------- Step 1 ----------------------------------------------\n";
   cluDB = runStep(inDB, inDBWorkingIndex, tmpDir, scoringMatrixFile, maxSeqLen, seqType, kmerSize, alphabetSize, maxResListLen, skip, aaBiasCorrection, zscoreThr, 2.0, evalThr, covThr, maxAlnNum, 1);
   cluSteps.push_back(cluDB);

   std::cout << "------------------------------- Step 2 ----------------------------------------------\n";
   cluDB = runStep(inDB, inDBWorkingIndex, tmpDir, scoringMatrixFile, maxSeqLen, seqType, kmerSize, alphabetSize, maxResListLen, skip, aaBiasCorrection, 100.0, 5.0, evalThr, covThr, maxAlnNum, 2);
   cluSteps.push_back(cluDB);

   std::cout << "------------------------------- Step 3 ----------------------------------------------\n";
   cluDB = runStep(inDB, inDBWorkingIndex, tmpDir, scoringMatrixFile, maxSeqLen, seqType, kmerSize, alphabetSize, maxResListLen, skip, aaBiasCorrection, 50.0, 7.2, evalThr, covThr, maxAlnNum, 3);
   cluSteps.push_back(cluDB);

   std::cout << "--------------------------- Merging databases ---------------------------------------\n";
   mergeClusteringResults(inDB, outDB, cluSteps);
}
