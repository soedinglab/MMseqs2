// Written by Martin Steinegger martin.steinegger@campus.lmu.de
//
// Converts PSSM or HHM to MMseqs profile format.
// MMseqs just stores the position specific score in 1 byte
//

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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


int cluster2profile(int argn,const char **argv)
{
    int err = EXIT_SUCCESS;

    std::string usage;
    usage.append("Converts a ffindex profile database to ffindex.\n");
    usage.append("USAGE: <clusteredDB> <sequenceDB> <ffindexDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@campus.lmu.de>.\n");
    std::vector<MMseqsParameter> profile_ffindex_par = {
            Parameters::PARAM_SUB_MAT,
            Parameters::PARAM_V};
    Parameters par;
    par.parseParameters(argn, argv, usage, profile_ffindex_par, 3);

    struct stat st;


    size_t offset_sequence = 0;
    size_t entry_num = 0;

    DBReader dbr_cluster_data(par.db1.c_str(), par.db1Index.c_str());
    dbr_cluster_data.open(DBReader::NOSORT);


    DBReader dbr_sequence_data(par.db2.c_str(), par.db2Index.c_str());
    dbr_sequence_data.open(DBReader::NOSORT);

    if(stat(par.db3.c_str(), &st) == 0) { errno = EEXIST; perror(par.db3.c_str()); return EXIT_FAILURE; }
    FILE * data_file  = fopen(par.db3.c_str(), "wb"); // binary file
    if( data_file == NULL) { perror(par.db3.c_str()); return EXIT_FAILURE; }

    if(stat(par.db3Index.c_str(), &st) == 0) { errno = EEXIST; perror(par.db3Index.c_str()); return EXIT_FAILURE; }
    FILE * index_file = fopen(par.db3Index.c_str(), "w+");
    if(index_file == NULL) { perror(par.db3Index.c_str()); return EXIT_FAILURE; }


    size_t maxSeqLen = 0;
    for(size_t i = 0; i < dbr_sequence_data.getSize(); i++) {
        maxSeqLen = std::max((size_t) dbr_sequence_data.getSeqLens()[i], maxSeqLen);
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
    SubstitutionMatrix * subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0);

    Debug(Debug::WARNING) << "Start computing profiles.\n";
    MultipleAlignment msaAligner(maxSeqLen, 10, subMat);
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
        char * seqData = dbr_sequence_data.getDataByDBKey(queryId);
        centerSeq.mapSequence(0, queryId, seqData);
        std::vector<Sequence *> seqSet;
        size_t pos = 0;

        while (*data != '\0' ){
            Util::parseKey(data, dbKey);
            if(strcmp(dbKey, queryId) != 0){
                char * dbSeqData = dbr_sequence_data.getDataByDBKey(dbKey);
                std::cout << dbKey << " " << dbSeqData << std::endl;
                dbSeqs[pos]->mapSequence(0, dbKey, dbSeqData);
                seqSet.push_back(dbSeqs[pos]);
                pos++;
            }
            data = Util::skipLine(data);

        }
        //parseHMM(data, profileBuffer, &elementSize, id, subMat);
        MultipleAlignment::MSAResult res = msaAligner.computeMSA(&centerSeq, seqSet, true);
        MultipleAlignment::print(res);

        const char * pssmData = pssmCalculator.computePSSMFromMSA(res.setSize, res.centerLength, res.msaSequence);
        size_t pssmDataSize = res.centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);

        //pssm.printProfile(res.centerLength);
        //pssm.printPSSM(res.centerLength);

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
    dbr_sequence_data.close();

    fclose(data_file);

    /* Sort the index entries and write back */
    rewind(index_file);
    ffindex_index_t* index = ffindex_index_parse(index_file, entry_num);
    if(index == NULL)
    {
        perror("ffindex_index_parse failed");
        exit(EXIT_FAILURE);
    }
    fclose(index_file);
    ffindex_sort_index_file(index);
    index_file = fopen(par.db3Index.c_str(), "w");
    if(index_file == NULL) { perror(par.db3Index.c_str()); return EXIT_FAILURE; }
    err += ffindex_write(index, index_file);
    Debug(Debug::WARNING) << "Done.\n";

    return err;
}

