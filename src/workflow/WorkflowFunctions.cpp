#include "WorkflowFunctions.h"



std::string runAlignment(std::string inDBData, std::string inDBWorkingIndex, std::string targetDBData, std::string perfResult,
                         std::string targetDBIndex, std::string tmpDir,
                         Parameters par, int step_num,
                         std::list<std::string>* tmpFiles){
    // alignment step
    struct timeval start, end;
    gettimeofday(&start, NULL);
    std::stringstream ss;
    ss << step_num;
    std::string alnDB_step = tmpDir + "/db_aln_step" + ss.str();
    std::string alnDB_step_index = alnDB_step + ".index";
    tmpFiles->push_back(alnDB_step);
    tmpFiles->push_back(alnDB_step_index);

    std::cout << "\n##### Step " << ss.str() << ": ALIGNMENT #####\n\n";
    Alignment* aln = new Alignment(inDBData, inDBWorkingIndex,
                                       targetDBData, targetDBIndex,
                                       perfResult, (perfResult + ".index"),
                                       alnDB_step, alnDB_step_index,
                                       par);
    std::cout << "Starting alignments calculation.\n";
    aln->run(par.maxResListLen, par.maxRejected);
    delete aln;

    gettimeofday(&end, NULL);
    int sec = end.tv_sec - start.tv_sec;
    std::cout << "\nTime for step " << step_num << " alignment: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    return alnDB_step;
}


std::string runPrefilter(std::string inDBData, std::string inDBWorkingIndex, std::string targetDBData,
                    std::string targetDBIndex, std::string tmpDir,
                    Parameters par, int step_num,
                    std::list<std::string>* tmpFiles){
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // prefiltering step
    std::stringstream ss;
    ss << step_num;
    std::string prefDB_step = tmpDir + "/db_pref_step" + ss.str();
    std::string prefDB_step_index = prefDB_step+ ".index";
    tmpFiles->push_back(prefDB_step);
    tmpFiles->push_back(prefDB_step_index);

    std::cout << "\n##### Step " << ss.str() << ": PREFILTERING #####\n\n";
    Prefiltering* pref = new Prefiltering (inDBData, inDBWorkingIndex,
                                           targetDBData, targetDBIndex,
                                           prefDB_step, prefDB_step_index,
                                           par);
    std::cout << "Starting prefiltering scores calculation.\n";
    pref->run();
    delete pref;

    gettimeofday(&end, NULL);
    int sec= end.tv_sec - start.tv_sec;
    std::cout << "\nTime for step " << step_num << " prefiltering: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
    return prefDB_step;
}

std::string runStep(std::string inDBData, std::string inDBWorkingIndex, std::string targetDBData,
                    std::string targetDBIndex, std::string tmpDir,
                    Parameters par, int step_num, int restart,
                    bool search, std::list<std::string>* tmpFiles){

    std::cout << "------------------------------------------------------------\n";
    std::cout << "GENERAL PARAMETERS:\n";
    std::cout << "max. sequence length: " << par.maxSeqLen << "\n";
    std::cout << "max. result list length: " << par.maxResListLen << "\n";
    std::cout << "PREFILTERING PARAMETERS:\n";
    std::cout << "k-mer size: " << par.kmerSize << "\n";
    std::cout << "z-score threshold: " << par.zscoreThr << "\n";
    std::cout << "sensitivity: " << par.sensitivity << "\n";
    std::cout << "ALIGNMENT PARAMETERS:\n";
    std::cout << "e-value threshold: " << par.evalThr << "\n";
    std::cout << "coverage threshold: " << par.covThr << "\n";
    std::cout << "Maximum number of rejects: " << par.maxRejected << "\n";
    std::cout << "------------------------------------------------------------\n\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);

    // prefiltering step
    std::stringstream ss;
    ss << step_num;
    std::string prefDB_step = tmpDir + "/db_pref_step" + ss.str();
    std::string prefDB_step_index = prefDB_step+ ".index";
    tmpFiles->push_back(prefDB_step);
    tmpFiles->push_back(prefDB_step_index);

    int sec;
    if (restart <= 1){
        std::cout << "\n##### Step " << ss.str() << ": PREFILTERING #####\n\n";
        Prefiltering* pref = new Prefiltering (inDBData, inDBWorkingIndex,
                targetDBData, targetDBIndex,
                prefDB_step, prefDB_step_index,
                par);
        std::cout << "Starting prefiltering scores calculation.\n";
        pref->run();
        delete pref;

        gettimeofday(&end, NULL);
        sec = end.tv_sec - start.tv_sec;
        std::cout << "\nTime for step " << step_num << " prefiltering: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n\n";
        gettimeofday(&start, NULL);
    }

    // alignment step
    std::string alnDB_step = tmpDir + "/db_aln_step" + ss.str();
    std::string alnDB_step_index = alnDB_step + ".index";
    tmpFiles->push_back(alnDB_step);
    tmpFiles->push_back(alnDB_step_index);

    if (restart <= 2){
        std::cout << "\n##### Step " << ss.str() << ": ALIGNMENT #####\n\n";
        Alignment* aln = new Alignment(inDBData, inDBWorkingIndex,
                targetDBData, targetDBIndex,
                prefDB_step, prefDB_step_index,
                alnDB_step, alnDB_step_index,
                par);
        std::cout << "Starting alignments calculation.\n";
        aln->run(par.maxResListLen, par.maxRejected);
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
    std::string cluDB_step_index = cluDB_step + ".index";
    tmpFiles->push_back(cluDB_step);
    tmpFiles->push_back(cluDB_step_index);

    if (restart <= 3){
        std::cout << "\n##### Step " << ss.str() << ": CLUSTERING #####\n\n";
        Clustering* clu = new Clustering(inDBData, inDBWorkingIndex,
                alnDB_step, alnDB_step_index,
                cluDB_step, cluDB_step_index,
                0, par.maxIteration,par.similarityScoreType, par.threads);
        clu->run(Parameters::SET_COVER);
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
        EXIT(EXIT_FAILURE);
    }
}

float getZscoreForSensitivity(int sensitivity){
    switch (sensitivity){
        case 1:
            return 400.0;
        case 2:
            return 200.0;
        case 3:
            return 100.0;
        case 4:
            return 50.0;
        case 5:
            return 40.0;
        case 6:
            return 30.0;
        case 7:
            return 20.0;
        case 8:
        case 9:
        case 10:
            return 10.0;
    }
    return 50.0;
}

void deleteTmpFiles(std::list<std::string>* tmpFiles){
    for (std::list<std::string>::iterator it = tmpFiles->begin(); it != tmpFiles->end(); it++){
        std::cout << "Deleting " << *it << "\n";
        if( remove((*it).c_str()) != 0 )
            std::cerr << "Error deleting file " << *it << "\n";
    }
}

