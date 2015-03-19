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
#include "Util.h"

void extractNewIndex(std::string seqDBIndex, std::string cluDBIndex, std::string newIndexFileName){

    ffindex_index_t* seq_index = Util::openIndex(seqDBIndex.c_str());

    ffindex_index_t* clu_index = Util::openIndex(cluDBIndex.c_str());

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
    // it has the size of all possible cluster (sequence amount)
    std::list<int>** mergedClustering = new std::list<int>*[dbr->getSize()];
    std::cout << "List amount "<< dbr->getSize() << std::endl;
    for (size_t i = 0; i < dbr->getSize(); i++){
        mergedClustering[i] = new std::list<int>();
    }

    // read the clustering from the first clustering step
    std::string firstCluStep = cluSteps.front();
    cluSteps.pop_front();
    std::string firstCluStepIndex = firstCluStep + ".index";
    DBReader* cluStepDbr = new DBReader(firstCluStep.c_str(), firstCluStepIndex.c_str());
    cluStepDbr->open(DBReader::NOSORT);

    for (size_t i = 0; i < cluStepDbr->getSize(); i++){
        size_t cluId = dbr->getId(cluStepDbr->getDbKey(i));
        std::stringstream lineSs (cluStepDbr->getData(i));
        std::string val;
        // go through the sequences in the cluster and add them to the initial clustering
        while (std::getline(lineSs, val)){
            size_t seqId = dbr->getId(val.c_str());
            mergedClustering[cluId]->push_back(seqId);
        }
    }
    cluStepDbr->close();
    delete cluStepDbr;
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
        for (size_t i = 0; i < cluStepDbr->getSize(); i++){
            size_t cluId = dbr->getId(cluStepDbr->getDbKey(i));
            char* cluData = cluStepDbr->getData(i);
            std::stringstream lineSs(cluData);
            std::string val;
            // go through the sequences in the cluster and add them and their clusters to the cluster of cluId
            // afterwards, delete the added cluster from the clustering
            while (std::getline(lineSs, val, '\n')){
                size_t seqId = dbr->getId(val.c_str());
                if(seqId != cluId) // to avoid copies of the same cluster list
                    mergedClustering[cluId]->splice(mergedClustering[cluId]->end(), *mergedClustering[seqId]);
            }

        }
        cluStepDbr->close();
        delete cluStepDbr;
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
    for (size_t i = 0; i < dbr->getSize(); i++){

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
    delete dbw;
    delete[] outBuffer;

    // delete the clustering data structure
    for (unsigned int i = 0; i < dbr->getSize(); i++){
        delete mergedClustering[i];
    }
    delete[] mergedClustering;
    std::cout << "...done.\n";

}

void runClustering(Parameters par,
        std::string inDB, std::string outDB, std::string tmpDir, int restart){

    std::string inDBIndex = inDB + ".index";
    float targetSens = par.sensitivity;
    float zscoreThr = getZscoreForSensitivity(targetSens);

    std::cout << "\nRunning the clustering with sensitivity " << par.sensitivity << "\n";
    std::cout << "Z-score threshold: " << par.zscoreThr << "\n";

    std::list<std::string>* tmpFiles = new std::list<std::string>();
    int oldMaxRejected = par.maxRejected;
    par.maxRejected = 10;
    par.zscoreThr = zscoreThr;
    std::string cluDB = runStep(inDB, inDBIndex, inDB, inDBIndex, tmpDir,
                                par, 1, restart, false, tmpFiles);
    par.maxRejected = oldMaxRejected;
    std::string cluDBIndex = cluDB + ".index";
    std::string outDBIndex = outDB + ".index";

    // copy the clustering databases to the right location
    copy(cluDBIndex, outDBIndex);
    copy(cluDB, outDB);

    if(!par.keepTempFiles) {
        deleteTmpFiles(tmpFiles);
    }
    delete tmpFiles;
}

void runCascadedClustering(Parameters par,
        std::string inDB, std::string outDB, std::string tmpDir, int restart, int step){

    std::cout << "\nRunning cascaded clustering for the database " << inDB << "\n";
    std::cout << "Target sensitivity: " << par.sensitivity << "\n";

    std::list<std::string>* tmpFiles = new std::list<std::string>();

    std::string inDBIndex = inDB + ".index";

    // copy index to a new location, it will be overwritten
    std::string inDBWorkingIndex = tmpDir + "/input_seqs.index";
    tmpFiles->push_back(inDBWorkingIndex);

    // save the original sequences to the working location
    // index will be overwritten in each clustering step
    ffindex_index_t* index_orig = Util::openIndex(inDBIndex.c_str());
    FILE* index_cpy_file = fopen(inDBWorkingIndex.c_str(), "w");
    ffindex_write(index_orig, index_cpy_file);
    fclose(index_cpy_file);

    std::list<std::string> cluSteps;
    std::string cluDB = "";

    std::cout << "\n\n";
    std::cout << "------------------------------- Step 1 ----------------------------------------------\n";

    float targetSensitivity = par.sensitivity;

    int local_restart;
    if (step > 1) // skip step
        local_restart = 4;
    else if (step == 1) // set the starting point of the clustering step
        local_restart = restart;
    else
        local_restart = 1;
    // set parameter for first step
    par.sensitivity = 1; // 1 is lowest sens
    par.zscoreThr = getZscoreForSensitivity( par.sensitivity );
    int oldMaxRejected = par.maxRejected;
    par.maxRejected = 10;
    int oldMaxResListLen= par.maxResListLen;
    par.maxResListLen = 50;
//    int oldSkip = par.skip;
//    par.skip = 2;
    cluDB = runStep(inDB, inDBWorkingIndex, inDB, inDBWorkingIndex, tmpDir,
                    par, 1, local_restart, false, tmpFiles);
    par.maxRejected   = oldMaxRejected;
    par.maxResListLen = oldMaxResListLen;
//    par.skip = oldSkip;

    cluSteps.push_back(cluDB);

    std::cout << "\n##### Updating databases #####\n\n";
    // extract new index containing only sequences remained after the clustering
    // and overwrite old index with the extracted index
    std::string cluDBIndex = cluDB + ".index";
    std::string newIndexFileName =  tmpDir + "/seqs_tmp.index";
    tmpFiles->push_back(newIndexFileName);
    // copy representing cluster sequences.
    // They are the input for the next iteration
    extractNewIndex(inDBWorkingIndex, cluDBIndex, newIndexFileName);
    copy(newIndexFileName, inDBWorkingIndex);

    std::cout << "------------------------------- Step 2 ----------------------------------------------\n";

    if (step > 2)
        // skip step
        local_restart = 4;
    else if (step == 2)
        // set the starting point of the clustering step
        local_restart = restart;
    else
        local_restart = 1;

    par.sensitivity = targetSensitivity / 2.0;
    par.zscoreThr =  getZscoreForSensitivity( par.sensitivity );
    oldMaxResListLen = par.maxResListLen;
    par.maxResListLen = 100;
    oldMaxRejected = par.maxRejected;
    par.maxRejected = 10;
    cluDB = runStep(inDB, inDBWorkingIndex, inDB, inDBWorkingIndex, tmpDir,
            par, 2, local_restart, false, tmpFiles);
    par.maxRejected   = oldMaxRejected;
    par.maxResListLen = oldMaxResListLen;
    cluSteps.push_back(cluDB);

    std::cout << "\n##### Updating databases #####\n\n";
    // extract new index containing only sequences remained after the clustering and overwrite old index with the extracted index
    cluDBIndex = cluDB + ".index";
    // copy cluster index into new index
    extractNewIndex(inDBWorkingIndex, cluDBIndex, newIndexFileName);
    copy(newIndexFileName, inDBWorkingIndex);

    std::cout << "------------------------------- Step 3 ----------------------------------------------\n";
    if (step > 3)
        // skip step
        local_restart = 4;
    else if (step == 3)
        // set the starting point of the clustering step
        local_restart = restart;
    else
        local_restart = 1;
    par.sensitivity = targetSensitivity;
    par.zscoreThr = getZscoreForSensitivity( par.sensitivity );
    cluDB = runStep(inDB, inDBWorkingIndex, inDB, inDBWorkingIndex, tmpDir,
                    par, 3, local_restart, false, tmpFiles);
    cluSteps.push_back(cluDB);

    std::cout << "--------------------------- Merging databases ---------------------------------------\n";
    mergeClusteringResults(inDB, outDB, cluSteps);

    if(!par.keepTempFiles) {
        deleteTmpFiles(tmpFiles);
    }
    delete tmpFiles;
}

void setWorkflowDefaults(Parameters* p) {
    p->localSearch = false;;
    p->spacedKmer = false;
}

int clusteringworkflow (int argc, const char * argv[]){

    std::string usage("\nCalculates the clustering of the sequences in the input database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: mmseqs_cluster <sequenceDB> <outDB> <tmpDir> [opts]\n");
                 //            "--restart          \t[int]\tRestart the clustering workflow starting with alignment or clustering.\n"
                 //            "                \t     \tThe value is in the range [1:3]: 1: restart from prefiltering  2: from alignment; 3: from clustering.\n"
                 //            "CASCADED CLUSTERING OPTIONS:\n"
    /*            "--step          \t[int]\tRestart the step of the cascaded clustering. For values in [1:3], the resprective step number, 4 is only the database merging.\n"
     "\nRESTART OPTIONS NOTE:\n"
     "                \t     \tIt is assumed that all valid intermediate results exist.\n"
     "                \t     \tValid intermediate results are taken from the tmp directory specified by the user.\n"
     "                \t     \tFor the cascaded clustering, --restart and --step options can be combined.\n"*/
    std::vector<MMseqsParameter> perfPar = {
        Parameters::PARAM_S,
        Parameters::PARAM_SUB_MAT,
        Parameters::PARAM_MAX_SEQS,
        Parameters::PARAM_CASCADED,
        Parameters::PARAM_MAX_SEQ_LEN,
        Parameters::PARAM_RESTART,
        Parameters::PARAM_STEP,
        Parameters::PARAM_V,
        Parameters::PARAM_SEARCH_MODE,
        Parameters::PARAM_SPACED_KMER_MODE,
        Parameters::PARAM_KEEP_TEMP_FILES
    };

    Parameters par;
    setWorkflowDefaults(&par);
    par.parseParameters(argc, argv, usage, perfPar, 3);
    par.localSearch = par.localSearch;
    par.spacedKmer = par.spacedKmer;
    par.covThr = 0.8;
    par.evalThr = 0.001;
    Debug::setDebugLevel(par.verbosity);
    int restart = par.restart;
    int step = par.step;
    
    if (par.cascaded)
        runCascadedClustering(par, par.db1, par.db2, par.db3, restart, step);
    else
        runClustering(par, par.db1, par.db2, par.db3, restart);

    return 0;
}
