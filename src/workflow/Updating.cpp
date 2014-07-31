#include <iostream>
#include <string>
#include <time.h>
#include <sys/time.h>

#include "../prefiltering/Prefiltering.h"
#include "../alignment/Alignment.h"
#include "../clustering/Clustering.h"
#include "WorkflowFunctions.h"

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

struct clu_entry_t {
    unsigned int id;
    clu_entry_t* next;
};

struct cluster_t {
    clu_entry_t* first;
    clu_entry_t* last;
    int clu_size;
};


int oldDBSize;
int newDBSize;

int deletedSeqs;
int sharedSeqs;
int newSeqs;

int seqsWithMatches;
int seqsWithoutMatches;
int newClus;

void printUsage(){

    std::string usage("\nUpdates the existing clustering of the previous database version with new sequences from the current version of the same database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: mmseqs_update <oldDB> <newDB> <oldDB_clustering> <outDB> <tmpDir> [opts]\n"
            "--sub-mat       \t[file]\tAmino acid substitution matrix file.\n"
            "--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n");
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexLastSeqDBBase, std::string* ffindexCurrentSeqDBBase, std::string* ffindexCluDBBase, std::string* ffindexOutDBBase, 
        std::string* tmpDir, std::string* scoringMatrixFile, size_t* maxSeqLen){
    if (argc < 4){
        printUsage();
        exit(EXIT_FAILURE);
    }

    ffindexLastSeqDBBase->assign(argv[1]);
    ffindexCurrentSeqDBBase->assign(argv[2]);
    ffindexCluDBBase->assign(argv[3]);
    ffindexOutDBBase->assign(argv[4]);
    tmpDir->assign(argv[5]);

    int i = 6;
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

void writeIndexes(std::string A_indexFile, std::string B_indexFile, std::string oldDBIndex, std::string newDBIndex){

    FILE* A_index_file = fopen(A_indexFile.c_str(), "w");
    FILE* B_index_file = fopen(B_indexFile.c_str(), "w");

    ffindex_index_t* index_old = openIndex(oldDBIndex.c_str());
    ffindex_index_t* index_new = openIndex(newDBIndex.c_str());

    // positions in the databases
    unsigned int i = 0;
    unsigned int j = 0;

    int deleted_cnt = 0;
    int new_cnt = 0;
    int shared_cnt = 0;
    while (i < index_old->n_entries && j < index_new->n_entries){
        ffindex_entry_t* e_i = ffindex_get_entry_by_index(index_old, i);
        ffindex_entry_t* e_j = ffindex_get_entry_by_index(index_new, j);
        int cmp = strcmp(&(e_i->name[0]), &(e_j->name[0]));
        if (cmp == 0){
            // this sequence is in both databases
            fprintf(A_index_file, "%s\t%zd\t%zd\n", e_j->name, e_j->offset, e_j->length);
            shared_cnt++;
            i++;
            j++;
        }
        else if (cmp < 0){
            // sequence was deleted from the old database
            deleted_cnt++;
            i++;
        }
        else{
            // this sequence is new
            fprintf(B_index_file, "%s\t%zd\t%zd\n", e_j->name, e_j->offset, e_j->length);
            new_cnt++;
            j++;
        }
    }
    while (i < index_old->n_entries){
        deleted_cnt++;
        i++;
    }
    // add the rest of the new database to the new sequences
    while (j < index_new->n_entries){
        ffindex_entry_t* e_j = ffindex_get_entry_by_index(index_new, j);
        fprintf(B_index_file, "%s\t%zd\t%zd\n", e_j->name, e_j->offset, e_j->length);
        new_cnt++;
        j++;
    }

    // set the global count variables
    oldDBSize = index_old->n_entries;
    newDBSize = index_new->n_entries;
    deletedSeqs = deleted_cnt;
    sharedSeqs = shared_cnt;
    newSeqs = new_cnt;

    fclose(A_index_file);
    fclose(B_index_file);
}


std::string runScoresCalculation(std::string queryDB, std::string queryDBIndex,
        std::string targetDB, std::string targetDBIndex,
        std::string tmpDir,
        std::string scoringMatrixFile, int maxSeqLen, int seqType,
        int kmerSize, int alphabetSize, size_t maxResListLen, int split, int skip, bool aaBiasCorrection, float zscoreThr, float sensitivity,
        double evalThr, double covThr, int maxAlnNum, std::string dbName, std::list<std::string>* tmpFiles){

    struct timeval start, end;
    gettimeofday(&start, NULL);

    // prefiltering step
    std::cout << "\n----------------------------- Prefiltering ------------------------\n";
    std::string prefDB = tmpDir + "/db_pref_" + dbName;
    std::string prefDBIndex = prefDB + ".index";
    tmpFiles->push_back(prefDB);
    tmpFiles->push_back(prefDBIndex);
    Prefiltering* pref = new Prefiltering (queryDB, queryDBIndex,
            targetDB, targetDBIndex,
            prefDB, prefDBIndex,
            scoringMatrixFile, sensitivity, kmerSize, alphabetSize, zscoreThr, maxSeqLen, seqType, aaBiasCorrection, split, skip);
    std::cout << "Starting prefiltering scores calculation.\n";
    pref->run(maxResListLen);
    delete pref;

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for the prefiltering: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    // alignment step
    std::cout << "------------------------------ Alignment --------------------------\n";
    std::string alnDB = tmpDir + "/db_aln_" + dbName;
    std::string alnDBIndex = alnDB + ".index";
    tmpFiles->push_back(alnDB);
    tmpFiles->push_back(alnDBIndex);
    Alignment* aln = new Alignment(queryDB, queryDBIndex,
            targetDB, targetDBIndex,
            prefDB, prefDBIndex,
            alnDB, alnDBIndex,
            scoringMatrixFile, evalThr, covThr, maxSeqLen, seqType);
    std::cout << "Starting alignments calculation.\n";
    aln->run(maxResListLen, 10);
    delete aln;

    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for alignments: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    return alnDB;
}


int readClustering(DBReader* currSeqDbr, std::string cluDB, unsigned int* id2rep, char** rep2cluName, cluster_t* clusters){

    DBReader* cluDbr = new DBReader(cluDB.c_str(), (cluDB + ".index").c_str());
    cluDbr->open(DBReader::NOSORT);

    int ret = cluDbr->getSize();

    for (unsigned int i = 0; i < cluDbr->getSize(); i++){
        id2rep[i] = UINT_MAX;
    }

    char* buf = new char[1000000];
    for (unsigned int i = 0; i < cluDbr->getSize(); i++){
        unsigned int repId = UINT_MAX;
        
        // parse the cluster
        char* cluData = cluDbr->getData(i);
        strcpy(buf, cluData);

        // first cluster member
        char* cluMemDbKey = strtok(buf, "\n");
        clu_entry_t* prev = 0;
        clu_entry_t* curr = 0;

        // store cluster members
        while (cluMemDbKey != 0){
            unsigned int cluMemId = currSeqDbr->getId(cluMemDbKey);
            // this cluster member is contained in the newest DB version
            if (cluMemId != UINT_MAX){
                // define a cluster representative if the representative is not set yet
                if (repId == UINT_MAX){
                    repId = cluMemId;
                    // remember the name of the cluster
                    strcpy(rep2cluName[repId], cluDbr->getDbKey(i));
                }
                id2rep[cluMemId] = repId;
                // create a cluster member entry
                // ATTENTION: consider counting cluster members first and allocate the memory at one piece (faster, no memory fragmentation)
                // if the program becomes too slow and/or the memory consumption is too high
                curr = new clu_entry_t;
                curr->id = cluMemId;
                curr->next = 0;
                if (prev != 0){
                    prev->next = curr;
                }
                prev = curr;
                // update the current clustering entry
                clusters[repId].clu_size++;
                if (clusters[repId].first == 0){
                    clusters[repId].first = curr;
                }
                clusters[repId].last = curr;
            }
            cluMemDbKey = strtok(NULL, "\n");
        }
    }
    cluDbr->close();
    return ret;
}

void appendToClustering(DBReader* currSeqDbr, std::string BIndexFile, std::string BA_base, unsigned int* id2rep, cluster_t* clusters, std::string Brest_indexFile){

    DBReader* BADbr = new DBReader(BA_base.c_str(), (BA_base + ".index").c_str());
    BADbr->open(DBReader::NOSORT);

    ffindex_index_t* Bindex = openIndex(BIndexFile.c_str());

    FILE* Brest_index_file = fopen(Brest_indexFile.c_str(), "w");

    seqsWithMatches = 0;
    seqsWithoutMatches = 0;
    char* buf = new char[1000000];
    for (unsigned int i = 0; i < BADbr->getSize(); i++){
        char* qKey = BADbr->getDbKey(i);
        unsigned int qId = currSeqDbr->getId(qKey);

        // find out which cluster the sequence belongs to
        char* alnData = BADbr->getData(i);
        strcpy(buf, alnData);

        char* tKey = strtok(buf, "\t");
        if (tKey != 0){
            unsigned int tId = currSeqDbr->getId(tKey);
            if (tId == UINT_MAX){
                std::cerr << "ERROR: DB key " << tKey << " is in the B->A alignment lists, but not in the new database!\n";
                exit(EXIT_FAILURE);
            }
            // find out the representative sequence of the cluster of the hit
            int repId = id2rep[tId];

            if(repId == -1){
                std::cout << "ERROR: database sequence " << tKey << " is not in the clustering!\n";
                std::cout << "Query from B: " << qKey << " matched this sequence.\n";
                continue;
            }

            // append new member to the cluster
            clu_entry_t* curr = new clu_entry_t; 
            curr->id = qId;
            curr->next = 0;

            clusters[repId].last->next = curr;
            clusters[repId].last = curr;
            clusters[repId].clu_size++;

            seqsWithMatches++;
        }
        else{
            ffindex_entry_t* e = ffindex_get_entry_by_name(Bindex, qKey);
            fprintf(Brest_index_file, "%s\t%zd\t%zd\n", e->name, e->offset, e->length);

            seqsWithoutMatches++;
        }
    }

    BADbr->close();
    fclose(Brest_index_file);
}

void writeResults(cluster_t* clusters, char** rep2cluName, DBReader* seqDbr, int seqDbSize, std::string outDB){

    DBWriter* dbw = new DBWriter(outDB.c_str(), (outDB + ".index").c_str());
    dbw->open();

    size_t BUFFER_SIZE = 1000000;
    char* outBuffer = new char[BUFFER_SIZE];

    for (int i = 0; i < seqDbSize; i++){
        // check if sequence i is a representative
        if (clusters[i].clu_size == 0)
            continue;

        // get the cluster name
        char* cluName = rep2cluName[i];
        std::stringstream res;
        clu_entry_t* e = clusters[i].first;
        while (e != 0){
            res << seqDbr->getDbKey(e->id) << "\n";
            e = e->next;
        }
        std::string cluResultsOutString = res.str();
        const char* cluResultsOutData = cluResultsOutString.c_str();
        if (BUFFER_SIZE < strlen(cluResultsOutData)){
            std::cerr << "Tried to process the clustering list for the cluster " << cluName << " , length of the list = " << clusters[i].clu_size << "\n";
            std::cerr << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " << cluResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters -> output buffer is already huge ;-)\n";
            continue;
        }
        memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length()*sizeof(char));
        dbw->write(outBuffer, cluResultsOutString.length(), cluName);
    }

    dbw->close();
}

int main (int argc, const char * argv[]){

    struct timeval start, end;
    gettimeofday(&start, NULL);

    // general parameters
    size_t maxSeqLen = 50000;
    int seqType = Sequence::AMINO_ACIDS;

    // parameter for the prefiltering
    int kmerSize = 6;
    int alphabetSize = 21;
    size_t maxResListLen = 100;
    int split = 0;
    int skip = 0;
    bool aaBiasCorrection = true;
    float zscoreThr = 50.0f;
    float sensitivity = 4.0;

    // parameters for the alignment
    double evalThr = 0.001;
    double covThr = 0.8;
    int maxAlnNum = 10;

    std::string lastSeqDB = "";
    std::string currentSeqDB = "";
    std::string cluDB = ""; 
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

    parseArgs(argc, argv, &lastSeqDB, &currentSeqDB, &cluDB, &outDB, &tmpDir, &scoringMatrixFile, &maxSeqLen);

    std::string lastSeqDBIndex = lastSeqDB + ".index";
    std::string currentSeqDBIndex = currentSeqDB + ".index";
    std::string cluDBIndex = cluDB + ".index";
    std::string outDBIndex = outDB + ".index";

    std::list<std::string>* tmpFiles = new std::list<std::string>();
    std::string AIndex = tmpDir + "/A.index";
    std::string BIndex = tmpDir + "/B.index";
    tmpFiles->push_back(AIndex);
    tmpFiles->push_back(BIndex);

    std::string Brest_indexFile = tmpDir + "/Brest.index";
    tmpFiles->push_back(Brest_indexFile);
    
    std::string BB_clu = tmpDir + "/BB_clu";
    std::string BB_clu_index = BB_clu + ".index";
    tmpFiles->push_back(BB_clu);
    tmpFiles->push_back(BB_clu_index);
    
    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    std::cout << "///////                   Init                             /////////////\n";
    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    // extract three indexes:
    // - A: last database version without deleted sequences
    // - B: sequences which are new in the database
    writeIndexes(AIndex, BIndex, lastSeqDBIndex, currentSeqDBIndex);


    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    std::cout << "///////            Calculating B->A scores                 /////////////\n";
    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    // calculate score for the updating
    // B->A scores
    std::string BA_base = runScoresCalculation(currentSeqDB, BIndex,
            currentSeqDB, AIndex,
            tmpDir,
            scoringMatrixFile, maxSeqLen, seqType,
            kmerSize, alphabetSize, maxResListLen, split, skip, aaBiasCorrection, zscoreThr, sensitivity,
            evalThr, covThr, maxAlnNum, "BA", tmpFiles);

    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    std::cout << "///////      Adding sequences to existing clusters         /////////////\n";
    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    // update the clustering
    DBReader* currSeqDbr = new DBReader(currentSeqDB.c_str(), currentSeqDBIndex.c_str());
    currSeqDbr->open(DBReader::NOSORT);

    // data structures for the clustering
    int seqDBSize = currSeqDbr->getSize();
    unsigned int* id2rep = new unsigned int[seqDBSize];
    char** rep2cluName = new char*[seqDBSize];
    for (int i = 0; i < seqDBSize; i++)
        rep2cluName[i] = new char[FFINDEX_MAX_ENTRY_NAME_LENTH];
    cluster_t* clusters = new cluster_t[seqDBSize];
    for (int i = 0; i < seqDBSize; i++){
        clusters[i].clu_size = 0;
        clusters[i].first = 0;
        clusters[i].last = 0;
    }

    std::cout << "Read the existing clustering...\n";
    // Read the existing clustering
    readClustering(currSeqDbr, cluDB, id2rep, rep2cluName, clusters);

    std::cout << "Append new sequences to the existing clustering...\n";
    // append sequences from the new database to the existing clustering based on the B->A alignment scores
    // write sequences without a match to a separate index (they will be clustered separately)
    appendToClustering(currSeqDbr, BIndex, BA_base, id2rep, clusters, Brest_indexFile);

    if (seqsWithoutMatches > 0){
        std::cout << "////////////////////////////////////////////////////////////////////////\n";
        std::cout << "///////            Calculating B->B scores                 /////////////\n";
        std::cout << "////////////////////////////////////////////////////////////////////////\n";
        // B->B scores
        std::string BB_base = runScoresCalculation(currentSeqDB, Brest_indexFile, 
                currentSeqDB, Brest_indexFile,
                tmpDir,
                scoringMatrixFile, maxSeqLen, seqType,
                kmerSize, alphabetSize, maxResListLen, split, skip, aaBiasCorrection, zscoreThr, sensitivity,
                evalThr, covThr, maxAlnNum, "BB", tmpFiles);

        std::cout << "////////////////////////////////////////////////////////////////////////\n";
        std::cout << "///////             Appending new clusters                 /////////////\n";
        std::cout << "////////////////////////////////////////////////////////////////////////\n";
        std::cout << "Cluster new sequences without a match to the existing clusters...\n";
        // cluster sequences without a match to the existing clusters separately
        // use the index generated in the previous step
        Clustering* clu = new Clustering(currentSeqDB, currentSeqDBIndex,
                BB_base, BB_base + ".index",
                BB_clu, BB_clu_index,
                0.0, 0, maxResListLen);
        clu->run(Clustering::SET_COVER); 

        std::cout << "Append generated clusters to the complete clustering...\n";
        // append B->B clusters to the clustering
        newClus = readClustering(currSeqDbr, BB_clu, id2rep, rep2cluName, clusters);
    }

    // write new clustering
    std::cout << "Write clustering results...\n";
    writeResults(clusters, rep2cluName, currSeqDbr, seqDBSize, outDB);
    std::cout << "done.\n";

    currSeqDbr->close();

    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    std::cout << "///////                   Statistics                            ////////\n";
    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    std::cout << "\nPrevios database version: " << oldDBSize << " entries.\n";
    std::cout << "New database vesion     : " << newDBSize << " entries.\n";
    std::cout << deletedSeqs << " entries were deleted,\n";
    std::cout << newSeqs << " entries are new,\n";
    std::cout << sharedSeqs << " entries are shared.\n\n";

    std::cout << seqsWithMatches << " new sequences had matches to the previous database version.\n";
    std::cout << "Remaining " << seqsWithoutMatches << " were grouped into " << newClus << " new clusters.\n";
 
    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for updating: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";

    deleteTmpFiles(tmpFiles);
    delete tmpFiles;

}
