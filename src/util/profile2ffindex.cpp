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

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include "DBReader.h"
#include "Util.h"
#include "Debug.h"
#define MAX_FILENAME_LIST_FILES 4096



void parsePSSM(char *data, char *profileBuffer, size_t *size, char *id, BaseMatrix * subMat){
    // go to readin position
    for(size_t i = 0; i < 2; i++){
        data = Util::skipLine(data);
    }
    // read aa_index
    char * words[22];
    Util::getWordsOfLine(data, words, 22);
    // A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    // header line
    int aa_index[20];
    for(size_t i = 0; i < 20; i++) {
        aa_index[i] = subMat->aa2int[(int)words[i][0]];
    }
    data = Util::skipLine(data);
    size_t curr_pos = 0;
    while (data[0] != '\n') {
        Util::getWordsOfLine(data, words, 22);
        for (size_t i = 0; i < 20; i++) {
            size_t writePos = curr_pos + aa_index[i];
            profileBuffer[writePos] = atoi(words[2+i]);
            // shifted score by -128 to avoid \0
            profileBuffer[writePos] = (profileBuffer[writePos] ^ 0x80);
        }
        data = Util::skipLine(data);
        curr_pos += 20;
    }
    *size = curr_pos;
}


// pasre HMM format
void parseHMM(char *data, char *profileBuffer, size_t *size, char *id, BaseMatrix * subMat){
    size_t l = 0;
    // find beging of profile information
    while( data[0] != '#') {
        data = Util::skipLine(data);
    }
    // go to readin position
    for(int i = 0; i < 5; i++)
        data = Util::skipLine(data);
    //ammino acids are ordered in HMM
    char * words[22];
    size_t curr_pos = 0;
    while (data[0] != '/' &&  data[1] != '/'){
        Util::getWordsOfLine(data, words, 22);
        for(size_t aa_num = 0; aa_num < 20; aa_num++) {
            // * entry: 0.0 probability
            if (words[aa_num+2][0] == '*'){
                Debug(Debug::ERROR) << "ERROR: 0 PROBABILITY FOR " << id << ".hhm AT " << l << "," << aa_num <<"\n";
                profileBuffer[curr_pos] = (char) -1;
            }
                // 0 entry: 1.0 probability
            else if (words[aa_num+2][0] == '0'){// integer number entry: 0.0 < probability < 1.0
                float score = Util::flog2(1.0f / subMat->getBackgroundProb(aa_num)) * subMat->getBitFactor();
                profileBuffer[curr_pos] = (char) floor (score + 0.5);
            } else {
                int entry = Util::fast_atoi(words[aa_num+2]);
                const float p = Util::fpow2( -(entry/1000.0f)); // back scaling from hhm
                const float backProb  = subMat->getBackgroundProb(aa_num);
                const float bitFactor = subMat->getBitFactor(); //TODO solve somehow this?!?

                double score = Util::flog2(p / backProb) * bitFactor;
                profileBuffer[curr_pos] = (char) floor (score + 0.5); // rounding
//                std::cout << aa_num << " " << subMat->int2aa[aa_num] << " " << profile_score[pos_in_profile] << " " << score << " " << entry << " " << p << " " << backProb << " " << bitFactor << std::endl;
            }
            // shifted score by -128 to avoid \0
            profileBuffer[curr_pos] = (profileBuffer[curr_pos] ^ 0x80);
            if(profileBuffer[curr_pos] == 0){
                Debug(Debug::ERROR) << "ERROR: 0 PSSM score is to big at id: " << id << ".hhm, pos: " << curr_pos << ", score:" <<
                (char)(profileBuffer[curr_pos] ^ 0x80) << "\n";
                EXIT(EXIT_FAILURE);
            }
            curr_pos++;
        }
        // go to next entry start
        for(int i = 0; i < 3; i++) // skip transitions
            data = Util::skipLine(data);
    }
    // return size of buffer
    *size = curr_pos;
}


int createprofiledb(int argn,const char **argv)
{
    int err = EXIT_SUCCESS;

    std::string usage;
    usage.append("Converts a ffindex profile database to ffindex.\n");
    usage.append("USAGE: <ffindexProfileDB>  <ffindexDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@campus.lmu.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.createprofiledb, 2);

    struct stat st;
    const char* data_filename = par.db2.c_str();
    if(stat(data_filename, &st) == 0) { errno = EEXIST; perror(data_filename); return EXIT_FAILURE; }
    FILE * data_file  = fopen(data_filename, "wb"); // binary file
    if( data_file == NULL) { perror(data_filename); return EXIT_FAILURE; }

    const char* index_filename = par.db2Index.c_str();
    if(stat(index_filename, &st) == 0) { errno = EEXIST; perror(index_filename); return EXIT_FAILURE; }
    FILE * index_file = fopen(index_filename, "w+");
    if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }

    size_t offset_sequence = 0;
    size_t entry_num = 0;

    DBReader dbr_data(par.db1.c_str(), std::string(par.db1  + ".index").c_str());
    dbr_data.open(DBReader::NOSORT);

    size_t maxElementSize = 0;
    for(size_t i = 0; i < dbr_data.getSize(); i++) {
        maxElementSize = std::max((size_t) dbr_data.getSeqLens()[i], maxElementSize);
    }
    BaseMatrix * subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);

    char * profileBuffer= new char[maxElementSize * Sequence::PROFILE_AA_SIZE];
    Debug(Debug::WARNING) << "Start converting profile to mmseqs profile.\n";
    for(size_t i = 0; i < dbr_data.getSize(); i++){
        char * data = dbr_data.getData(i);
        char * id   = dbr_data.getDbKey(i);
        size_t elementSize = 0;
        if(par.profileMode == Parameters::PROFILE_MODE_HMM){
            parseHMM(data, profileBuffer, &elementSize, id, subMat);
        } else if(par.profileMode == Parameters::PROFILE_MODE_HMM) {
            parsePSSM(data, profileBuffer, &elementSize, id, subMat);
        } else{
            Debug(Debug::ERROR) << "Wrong profile mode.\n";
            EXIT(EXIT_FAILURE);
        }

        ffindex_insert_memory(data_file,  index_file,     &offset_sequence,
                              profileBuffer,  elementSize , id);
    }
    entry_num = dbr_data.getSize();

    delete [] profileBuffer;
    delete subMat;
    dbr_data.close();
    fclose(data_file);

    /* Sort the index entries and write back */
    fclose(index_file);
    index_file = fopen(index_filename, "r+");
    ffindex_index_t* index = ffindex_index_parse(index_file, entry_num);
    if(index == NULL)
    {
        perror("ffindex_index_parse failed");
        EXIT(EXIT_FAILURE);
    }
    fclose(index_file);
    ffindex_sort_index_file(index);
    index_file = fopen(par.db2Index.c_str(), "w");
    if(index_file == NULL) { perror(par.db2Index.c_str()); return EXIT_FAILURE; }
    err += ffindex_write(index, index_file);
    Debug(Debug::WARNING) << "Done.\n";

    return err;
}

