#include <iostream>
#include <string>
#include <time.h>
#include <sys/time.h>

#include "../prefiltering/Prefiltering.h"
#include "../alignment/Alignment.h"
#include "../clustering/Clustering.h"

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

struct clu_entry_t {
    int id;
    clu_entry_t* next;
};

struct cluster_t {
    clu_entry_t* first;
    clu_entry_t* last;
    int clu_size;
};

void printUsage(){

    std::string usage("\nUpdates the existing clustering of the previous database version with new sequences from the current version of the same database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de))\n\n");
    usage.append("USAGE: update ffindexOldSeqDBBase ffindexCurrentSeqDBBase ffindexcluDBBase ffindexOutDBBase tmpDir [opts]\n"
            "-m              \t[file]\tAmino acid substitution matrix file.\n"
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
    std::cout << "-------\n";
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

    std::cout << "Previos database version: " << index_old->n_entries << " entries.\n";
    std::cout << "New database vesion     : " << index_new->n_entries << " entries.\n";
    std::cout << deleted_cnt << " entries were deleted,\n";
    std::cout << new_cnt << " entries are new,\n";
    std::cout << shared_cnt << " entries are shared.\n";
    
    fclose(A_index_file);
    fclose(B_index_file);
}


std::string runScoresCalculation(std::string queryDB, std::string queryDBIndex,
        std::string targetDB, std::string targetDBIndex,
        std::string tmpDir,
        std::string scoringMatrixFile, int maxSeqLen, int seqType,
        int kmerSize, int alphabetSize, size_t maxResListLen, int skip, bool aaBiasCorrection, float zscoreThr, float sensitivity,
        double evalThr, double covThr, int maxAlnNum, std::string dbName){

    struct timeval start, end;
    gettimeofday(&start, NULL);

    // prefiltering step
    std::cout << "\n----------------------------- Prefiltering ------------------------\n";
    std::string prefDB = tmpDir + "/db_pref_" + dbName;
    Prefiltering* pref = new Prefiltering (queryDB, queryDBIndex,
            targetDB, targetDBIndex,
            prefDB, prefDB + ".index",
            scoringMatrixFile, sensitivity, kmerSize, alphabetSize, zscoreThr, maxSeqLen, seqType, aaBiasCorrection, skip);
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
    Alignment* aln = new Alignment(queryDB, queryDBIndex,
            targetDB, targetDBIndex,
            prefDB, prefDB + ".index",
            alnDB, alnDB + ".index",
            scoringMatrixFile, evalThr, covThr, maxSeqLen, seqType);
    std::cout << "Starting alignments calculation.\n";
    aln->run(maxResListLen);
    delete aln;

    gettimeofday(&end, NULL);
    sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for alignments: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    gettimeofday(&start, NULL);

    return alnDB;
}


void readClustering(DBReader* currSeqDbr, std::string cluDB, int* id2clu, cluster_t* clusters){

    DBReader* cluDbr = new DBReader(cluDB.c_str(), (cluDB + ".index").c_str());
    cluDbr->open(DBReader::NOSORT);

    char* buf = new char[1000000];
    for (unsigned int i = 0; i < cluDbr->getSize(); i++){
        char* repDbKey = cluDbr->getDbKey(i);
        unsigned int repId = currSeqDbr->getId(repDbKey);

        // the representative is not in the newest DB version
        if (repId == UINT_MAX)
            continue;
        
        // the representative is contained in the newest DB version
        // parse the cluster
        id2clu[repId] = repId;

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
                id2clu[cluMemId] = repId;
                // create a cluster member entry
                // ATTENTION: consider counting cluster members first and allocate the memory at one piece (faster, no memory fragmentation)
                // if the program becomes too slow or the memory consumption is too high
                curr = new clu_entry_t;
                curr->id = cluMemId;
                curr->next = 0;
                if (prev != 0){
                    prev->next = curr;
                }
                prev = curr;
                // update clustering entry
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
}

int appendToClustering(DBReader* currSeqDbr, std::string BIndexFile, std::string BA_base, int* id2clu, cluster_t* clusters, std::string Brest_indexFile){

    DBReader* BADbr = new DBReader(BA_base.c_str(), (BA_base + ".index").c_str());
    BADbr->open(DBReader::NOSORT);

    ffindex_index_t* Bindex = openIndex(BIndexFile.c_str());

    FILE* Brest_index_file = fopen(Brest_indexFile.c_str(), "w");

    int cnt = 0;
    char* buf = new char[1000000];
    for (unsigned int i = 0; i < BADbr->getSize(); i++){
        char* qKey = BADbr->getDbKey(i);
        std::cout << "Query: " << qKey << "\n";
        unsigned int qId = currSeqDbr->getId(qKey);

        // find out which cluster the sequence belongs to
        char* alnData = BADbr->getData(i);
        strcpy(buf, alnData);

        char* tKey = strtok(buf, "\t");
        if (tKey != 0){
            std::cout << " -> " << tKey << "\n";
            unsigned int tId = currSeqDbr->getId(tKey);
            if (tId == UINT_MAX){
                std::cerr << "ERROR: DB key " << tKey << " contained in the old database not found in the new database!\n";
                exit(1);
            }
            // find out the representative sequence of the cluster of the hit
            unsigned int repId = id2clu[tId];
            // append new member to the cluster
            clu_entry_t* curr = new clu_entry_t;
            curr->id = qId;
            curr->next = 0;

            clusters[repId].last->next = curr;
            clusters[repId].last = curr;
            clusters[repId].clu_size++;
        }
        else{
            std::cout << " without a match\n";
            ffindex_entry_t* e = ffindex_get_entry_by_name(Bindex, qKey);
            fprintf(Brest_index_file, "%s\t%zd\t%zd\n", e->name, e->offset, e->length);
            cnt++;
        }
    }

    fclose(Brest_index_file);
    BADbr->close();

    return cnt;
}

void writeResults(cluster_t* clusters, DBReader* seqDbr, int seqDbSize, std::string outDB){

    DBWriter* dbw = new DBWriter(outDB.c_str(), (outDB + ".index").c_str());
    dbw->open();

    size_t BUFFER_SIZE = 1000000;
    char* outBuffer = new char[BUFFER_SIZE];

    for (int i = 0; i < seqDbSize; i++){
        // check if sequence i is a representative
        if (clusters[i].clu_size == 0)
            continue;

        char* repDbKey = seqDbr->getDbKey(i);
        std::stringstream res;
        clu_entry_t* e = clusters[i].first;
        while (e != 0){
            res << seqDbr->getDbKey(e->id) << "\n";
            e = e->next;
        }
        std::string cluResultsOutString = res.str();
        const char* cluResultsOutData = cluResultsOutString.c_str();
        if (BUFFER_SIZE < strlen(cluResultsOutData)){
            std::cerr << "Tried to process the clustering list for the query " << repDbKey << " , length of the list = " << clusters[i].clu_size << "\n";
            std::cerr << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " << cluResultsOutString.length() << ")\nIncrease buffer size or reconsider your parameters -> output buffer is already huge ;-)\n";
            continue;
        }
        memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length()*sizeof(char));
        dbw->write(outBuffer, cluResultsOutString.length(), repDbKey);
    }

    dbw->close();
}

void updateClustering(std::string cluDB_base, std::string BA_base, std::string BB_base, std::string currSeqDB_base, std::string outDB, std::string tmpDir){
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
    int skip = 0;
    bool aaBiasCorrection = false;
    float zscoreThr = 300.0f;
    float sensitivity = 7.2;

    // parameters for the alignment
    double evalThr = 0.001;
    double covThr = 0.8;
    int maxAlnNum = 10;

    std::string lastSeqDB = "";
    std::string currentSeqDB = "";
    std::string cluDB = ""; 
    std::string outDB = "";
    std::string scoringMatrixFile = "";
    std::string tmpDir = "";

    parseArgs(argc, argv, &lastSeqDB, &currentSeqDB, &cluDB, &outDB, &tmpDir, &scoringMatrixFile, &maxSeqLen);

    std::string lastSeqDBIndex = lastSeqDB + ".index";
    std::string currentSeqDBIndex = currentSeqDB + ".index";
    std::string cluDBIndex = cluDB + ".index";
    std::string outDBIndex = outDB + ".index";

    std::string AIndex = tmpDir + "/A.index";
    std::string BIndex = tmpDir + "/B.index";

    std::string Brest_indexFile = tmpDir + "/Brest.index";
    
    std::string BB_clu = tmpDir + "/BB_clu";
    
    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    std::cout << "///////                   Init                             /////////////\n";
    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    // extract three indexes:
    // - A: last database version without deleted sequences
    // - B: sequences which are new in the database
    // - deleted: sequences from the last database version which are deleted in the current database version
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
            kmerSize, alphabetSize, maxResListLen, skip, aaBiasCorrection, zscoreThr, sensitivity,
            evalThr, covThr, maxAlnNum, "BA");


    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    std::cout << "///////      Adding sequences to existing clusters         /////////////\n";
    std::cout << "////////////////////////////////////////////////////////////////////////\n";
    // update the clustering
    DBReader* currSeqDbr = new DBReader(currentSeqDB.c_str(), currentSeqDBIndex.c_str());
    currSeqDbr->open(DBReader::NOSORT);

    // data structures for the clustering
    int seqDBSize = currSeqDbr->getSize();
    int* id2clu = new int[seqDBSize];
    cluster_t* clusters = new cluster_t[seqDBSize];
    for (int i = 0; i < seqDBSize; i++){
        clusters[i].clu_size = 0;
        clusters[i].first = 0;
        clusters[i].last = 0;
    }

    std::cout << "Read the existing clustering...\n";
    // Read the existing clustering
    readClustering(currSeqDbr, cluDB, id2clu, clusters);

    std::cout << "Append new sequences to the existing clustering...\n";
    // append sequences from the new database to the existing clustering based on the B->A alignment scores
    // write sequences without a match to a separate index (they will be clustered separately)
    int cnt = appendToClustering(currSeqDbr, BIndex, BA_base, id2clu, clusters, Brest_indexFile);

    if (cnt > 0){
        std::cout << "////////////////////////////////////////////////////////////////////////\n";
        std::cout << "///////            Calculating B->B scores                 /////////////\n";
        std::cout << "////////////////////////////////////////////////////////////////////////\n";
        // B->B scores
        std::string BB_base = runScoresCalculation(currentSeqDB, Brest_indexFile, 
                currentSeqDB, Brest_indexFile,
                tmpDir,
                scoringMatrixFile, maxSeqLen, seqType,
                kmerSize, alphabetSize, maxResListLen, skip, aaBiasCorrection, zscoreThr, sensitivity,
                evalThr, covThr, maxAlnNum, "BB");

        std::cout << "////////////////////////////////////////////////////////////////////////\n";
        std::cout << "///////             Appending new clusters                 /////////////\n";
        std::cout << "////////////////////////////////////////////////////////////////////////\n";
        std::cout << "Cluster new sequences without a match to the existing clusters...\n";
        // cluster sequences without a match to the existing clusters separately
        // use the index generated in the previous step
        Clustering* clu = new Clustering(currentSeqDB, currentSeqDBIndex,
                BB_base, BB_base + ".index",
                BB_clu, BB_clu + ".index",
                0.0);
        clu->run(Clustering::SET_COVER); 

        std::cout << "Append generated clusters to the complete clustering...\n";
        // append B->B clusters to the clustering
        readClustering(currSeqDbr, BB_clu, id2clu, clusters);
    }

    std::cout << "Write clustering results...\n";
    // write new clustering
    writeResults(clusters, currSeqDbr, seqDBSize, outDB);

    currSeqDbr->close();

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for updating: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";

}
