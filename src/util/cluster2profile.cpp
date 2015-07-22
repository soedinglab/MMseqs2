// Written by Martin Steinegger martin.steinegger@campus.lmu.de
//
// Converts PSSM or HHM to MMseqs profile format.
// MMseqs just stores the position specific score in 1 byte
//

#include "cluster2profile.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <BaseMatrix.h>
#include <SubstitutionMatrix.h>
#include <Parameters.h>
#include <Sequence.h>
#include <MultipleAlignment.h>
#include <PSSMCalculator.h>
extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include "DBReader.h"
#include "Util.h"
#include "Debug.h"
#define MAX_FILENAME_LIST_FILES 4096


int runResult2Profile(std::string resultDb, std::string queryDb, std::string targetDb,
                      std::string outDb, std::string subMatPath){
    int err = EXIT_SUCCESS;

    struct stat st;

    size_t offset_sequence = 0;
    size_t entry_num = 0;

    DBReader dbr_cluster_data(resultDb.c_str(), (resultDb + ".index").c_str());
    dbr_cluster_data.open(DBReader::NOSORT);
    DBReader * qdbr = NULL;
    DBReader * tdbr = NULL;
    bool sameDatabase = true;
    qdbr = new DBReader(queryDb.c_str(), (queryDb + ".index").c_str());
    qdbr->open(DBReader::NOSORT);
    tdbr = qdbr;
    if(queryDb.compare(targetDb) != 0){
        sameDatabase = false;
        tdbr = new DBReader(targetDb.c_str(), (targetDb + ".index").c_str());
        tdbr->open(DBReader::NOSORT);
    }

    if(stat(outDb.c_str(), &st) == 0) { errno = EEXIST; perror(outDb.c_str()); return EXIT_FAILURE; }
    FILE * data_file  = fopen(outDb.c_str(), "wb"); // binary file
    if( data_file == NULL) { perror(outDb.c_str()); return EXIT_FAILURE; }
    std::string outDbIndex = (outDb+".index");
    if(stat(outDbIndex.c_str(), &st) == 0) { errno = EEXIST; perror(outDbIndex.c_str()); return EXIT_FAILURE; }
    FILE * index_file = fopen(outDbIndex.c_str(), "w+");
    if(index_file == NULL) { perror(outDbIndex.c_str()); return EXIT_FAILURE; }
    // find longest sequence
    size_t maxSeqLen = 0;
    for(size_t i = 0; i < qdbr->getSize(); i++) {
        maxSeqLen = std::max((size_t) qdbr->getSeqLens()[i], maxSeqLen);
    }
    for(size_t i = 0; i < tdbr->getSize(); i++) {
        maxSeqLen = std::max((size_t) tdbr->getSeqLens()[i], maxSeqLen);
    }
    //Find the max set size
    size_t maxSetSize = 0;
    for(size_t i = 0; i < dbr_cluster_data.getSize(); i++) {
        char * data = dbr_cluster_data.getData(i);
        size_t currClusterEntry = 0;
        size_t pos = 0;
        while(data[pos] != '\0'){
            if(data[pos] == '\n'){
                currClusterEntry++;
            }
            pos++;
        }
        maxSetSize = std::max(maxSetSize, currClusterEntry);
    }
    maxSetSize += 1;
    SubstitutionMatrix * subMat = new SubstitutionMatrix(subMatPath.c_str(), 2.0);

    Debug(Debug::WARNING) << "Start computing profiles.\n";
    MultipleAlignment msaAligner(maxSeqLen, maxSetSize, subMat);
    PSSMCalculator pssmCalculator(subMat, maxSeqLen);
    Sequence centerSeq(maxSeqLen, subMat->aa2int, subMat->int2aa, Sequence::AMINO_ACIDS, 0, false);
    Sequence ** dbSeqs = new Sequence *[maxSetSize];
    for(size_t i = 0; i < maxSetSize; i++) {
        dbSeqs[i] = new Sequence(maxSeqLen, subMat->aa2int, subMat->int2aa, Sequence::AMINO_ACIDS, 0, false);
    }

    for(size_t i = 0; i < dbr_cluster_data.getSize(); i++){
        char dbKey[255+1];
        char * data = dbr_cluster_data.getData(i);
        char * queryId   = dbr_cluster_data.getDbKey(i);
        char * seqData = qdbr->getDataByDBKey(queryId);
        centerSeq.mapSequence(0, queryId, seqData);
        std::vector<Sequence *> seqSet;
        size_t pos = 0;

        while (*data != '\0' ){
            Util::parseKey(data, dbKey);
            if(strcmp(dbKey, queryId) != 0){
                char * dbSeqData = tdbr->getDataByDBKey(dbKey);
                dbSeqs[pos]->mapSequence(0, dbKey, dbSeqData);
                seqSet.push_back(dbSeqs[pos]);
                pos++;
            }
            data = Util::skipLine(data);
        }
        //parseHMM(data, profileBuffer, &elementSize, id, subMat);
        MultipleAlignment::MSAResult res = msaAligner.computeMSA(&centerSeq, seqSet, true);
        //MultipleAlignment::print(res);

        char * pssmData = (char*) pssmCalculator.computePSSMFromMSA(res.setSize, res.centerLength, res.msaSequence);
        size_t pssmDataSize = res.centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);
        for(size_t i = 0; i < pssmDataSize; i++){
            pssmData[i] = pssmData[i] ^ 0x80;
        }
        //pssm.printProfile(res.centerLength);
        //pssmCalculator.printPSSM(res.centerLength);

        ffindex_insert_memory(data_file,  index_file,     &offset_sequence,
                              (char *) pssmData,  pssmDataSize , queryId);

        seqSet.clear();
    }
    // clean memeory
    for(size_t i = 0; i < maxSetSize; i++){
        delete dbSeqs[i];
    }
    delete [] dbSeqs;
    delete subMat;

    entry_num = dbr_cluster_data.getSize();
    // close reader
    dbr_cluster_data.close();
    qdbr->close();
    delete qdbr;
    if(sameDatabase == false) {
        tdbr->close();
        delete tdbr;
    }

    fclose(data_file);

    /* Sort the index entries and write back */
    fclose(index_file);
	index_file = fopen(outDbIndex.c_str(), "r+");
    ffindex_index_t* index = ffindex_index_parse(index_file, entry_num);
    if(index == NULL)
    {
        perror("ffindex_index_parse failed");
        EXIT(EXIT_FAILURE);
    }
    fclose(index_file);
    ffindex_sort_index_file(index);
    index_file = fopen(outDbIndex.c_str(), "w");
    if(index_file == NULL) { perror(outDbIndex.c_str()); return EXIT_FAILURE; }
    err += ffindex_write(index, index_file);
    Debug(Debug::WARNING) << "Done.\n";
    return err;
}


int cluster2profile(int argn,const char **argv)
{

    std::string usage;
    usage.append("Converts a ffindex profile database to ffindex.\n");
    usage.append("USAGE: <clusteredDB> <queryDB> <targetDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@campus.lmu.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.createprofiledb, 4);

    return runResult2Profile(par.db1, par.db2, par.db3,
                             par.db4, par.scoringMatrixFile);
}
