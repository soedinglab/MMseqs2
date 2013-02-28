#include <iostream>

#include "DBReader.h"
#include "DBWriter.h"
#include "Sequence.h"
#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "QueryTemplateMatcher.h"


void printUsage(){

    std::string usage("\nCalculates similarity scores between all sequences in the query database and all sequences in the target database.\n");
    usage.append("Written by Maria Hauser (mhauser@genzentrum.lmu.de)\n\n");
    usage.append("USAGE: kClust2_pref ffindexQueryDBBase ffindexTargetDBBase scoringMatrixFile ffindexOutDBBase [opts]\n"
            "-t\t[double]\tSimilarity threshold (Prefiltering threshold).\n");
    std::cout << usage;
}

void parseArgs(int argc, const char** argv, std::string* ffindexQueryDBBase, std::string* ffindexTargetDBBase, std::string* scoringMatrixFile, std::string* ffindexOutDBBase, short* prefThr){
    if (argc < 5){
        printUsage();
        exit(EXIT_FAILURE);
    }

    ffindexQueryDBBase->assign(argv[1]);
    ffindexTargetDBBase->assign(argv[2]);
    scoringMatrixFile->assign(argv[3]);
    ffindexOutDBBase->assign(argv[4]);
    int i = 5;
    while (i < argc){
        if (strcmp(argv[i], "-t") == 0){
            if (++i < argc){
                *prefThr = atof(argv[i]);
                i++;
            }
            else {
                printUsage();
                std::cerr << "No value provided for -t\n";
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

int main (int argc, const char * argv[])
{
    std::string queryDB = "";
    std::string queryDBIndex = "";
    std::string targetDB = "";
    std::string targetDBIndex = "";
    std::string outDB = "";
    std::string outDBIndex = "";
    std::string scoringMatrixFile = "";
    
    short prefThr = 6000;
    int kmerSize =  6;
    short kmerThr = 27;
    int maxSeqLen = 40000;

    parseArgs(argc, argv, &queryDB, &targetDB, &scoringMatrixFile, &outDB, &prefThr);

    queryDBIndex = queryDB + ".index";
    targetDBIndex = targetDB + ".index";
    outDBIndex = outDB + ".index";

    std::cout << "Init data structures:\n";
    DBReader* qdbr = new DBReader(queryDB.c_str(), queryDBIndex.c_str());
    qdbr->open();

    DBReader* tdbr = new DBReader(targetDB.c_str(), targetDBIndex.c_str());
    tdbr->open();

    DBWriter* dbw = new DBWriter(outDB.c_str(), outDBIndex.c_str(), 1);
    dbw->open();

    // init the substitution matrices
    SubstitutionMatrix* subMat = new SubstitutionMatrix (scoringMatrixFile.c_str());
//    ReducedMatrix* subMat = new ReducedMatrix(sMat->probMatrix, 13);
//    BaseMatrix::print(subMat->subMatrix, subMat->int2aa, subMat->alphabetSize);
    
    int targetDBSize = tdbr->getSize();
    Sequence* seq = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa);
    
    std::cout << "Index table fill...\n";
    // fill and init the index table
    IndexTable* indexTable = new IndexTable(subMat->alphabetSize, kmerSize);
    for (int id = 0; id < targetDBSize; id++){
        if (id > 0 && id % 100000 == 0){
            std::cout << id << " ";
            indexTable->checkSizeAndCapacity();
            std::cout << "\n";
        }
        seq->id = id;
        char* seqData = tdbr->getData(id);
        std::string str(seqData);
        seq->mapSequence(seqData);
        indexTable->addSequence(seq);
    }
    std::cout << "Index table init...\n";
    indexTable->init();

    std::cout << "Substitution matrices...\n";
    ExtendedSubstitutionMatrix* _2merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 2, subMat->alphabetSize);
    ExtendedSubstitutionMatrix* _3merSubMatrix = new ExtendedSubstitutionMatrix(subMat->subMatrix, 3, subMat->alphabetSize);
    QueryTemplateMatcher* matcher = new QueryTemplateMatcher(_2merSubMatrix, _3merSubMatrix, indexTable, kmerThr, prefThr, kmerSize, targetDBSize, subMat->alphabetSize); 

    // calculate prefiltering scores for each sequence in the database
    std::cout << "Calculating prefiltering scores!\n";
    std::list<hit_t>* prefResults;
    int queryDBSize = qdbr->getSize();
    char* outBuffer = new char[1000000];

    double kmersPerPos = 0.0;
    int dbMatches = 0;
    int resSize = 0;
    for (int id = 0; id < queryDBSize; id++){
        char* seqData = qdbr->getData(id);
        seq->id = id;
        seq->mapSequence(seqData);
        std::cout << "Sequence: " << qdbr->getDbKey(id) << " (length: " << seq->L << ")\n";
//        seq->print();

        prefResults = matcher->matchQuery(seq);

        std::list<hit_t>::iterator iter;
        std::stringstream prefResultsOut;

        for (iter = prefResults->begin(); iter != prefResults->end(); iter++){
//            std::cout << tdbr->getDbKey(iter->seqId);
//            std::cout << "\tscore: " << iter->prefScore << "\n";
            prefResultsOut << tdbr->getDbKey(iter->seqId) << "\t" << iter->prefScore << "\n";
        }

        std::string prefResultsOutString = prefResultsOut.str();
        const char* prefResultsOutData = prefResultsOutString.c_str();
        memcpy(outBuffer, prefResultsOutData, prefResultsOutString.length()*sizeof(char));
        dbw->write(outBuffer, prefResultsOutString.length(), qdbr->getDbKey(id), 0);

        kmersPerPos += seq->stats->kmersPerPos;
        dbMatches += seq->stats->dbMatches;
        resSize += prefResults->size();
    }
    kmersPerPos /= queryDBSize;
    double dbMatchesPerSeq =  (double)dbMatches/(double)queryDBSize;
    double prefPassedPerSeq = (double)resSize/(double)queryDBSize;
    std::cout << kmersPerPos << " k-mers per position.\n";
    std::cout << dbMatchesPerSeq << " DB matches per sequence.\n";
    std::cout << prefPassedPerSeq << " sequences passed prefiltering per query sequence.\n";

    qdbr->close();
    if (queryDB.compare(targetDB) != 0)
        tdbr->close();
    dbw->close();

    return 0;
}
