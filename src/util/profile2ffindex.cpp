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

int openFile(const char *filename, FILE **pFILE);
int sortIndex(char *index_filename, size_t entry_num);

void parsePSSM(char *data, char *profileBuffer, size_t *size, const char *id, BaseMatrix * subMat){
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
void parseHMM(char *data, std::string *header, char *profileBuffer, size_t *size, const char *id, BaseMatrix *subMat) {
    size_t l = 0;
    // find name tag
    while(data[0] !='N' && data[1] != 'A' && data[2] != 'M' && data[3] != 'E') {
        data = Util::skipLine(data);
    }
    // parse NAME entry
    const char * startData = data;
    data = Util::skipLine(data);
    const char * endData = data;
    header->append(startData + 6, endData - (startData + 6));
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
    std::string index_filename_str(par.db2);
    index_filename_str.append(".index");
    char *index_filename = (char *) index_filename_str.c_str();
    std::string data_filename_hdr_str(par.db2);
    data_filename_hdr_str.append("_h");
    char *data_filename_hdr  = (char *)data_filename_hdr_str.c_str() ;
    std::string index_filename_hdr_str(par.db2);
    index_filename_hdr_str.append("_h.index");
    char *index_filename_hdr = (char *)index_filename_hdr_str.c_str() ;
    FILE *data_file, *index_file, *fasta_file, *data_file_hdr, *index_file_hdr;
    struct stat st;

    openFile(par.db2.c_str(), &data_file);
    openFile(index_filename, &index_file);
    openFile(data_filename_hdr, &data_file_hdr);
    openFile(index_filename_hdr, &index_file_hdr);

    size_t offset_profile = 0;
    size_t offset_header = 0;

    size_t entry_num = 0;

    DBReader<unsigned int> dbr_data(par.db1.c_str(), std::string(par.db1  + ".index").c_str());
    dbr_data.open(DBReader<unsigned int>::NOSORT);

    size_t maxElementSize = 0;
    for(size_t i = 0; i < dbr_data.getSize(); i++) {
        maxElementSize = std::max((size_t) dbr_data.getSeqLens()[i], maxElementSize);
    }
    BaseMatrix * subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);

    char * profileBuffer= new char[maxElementSize * Sequence::PROFILE_AA_SIZE];
    Debug(Debug::WARNING) << "Start converting profile to mmseqs profile.\n";
    for(size_t i = 0; i < dbr_data.getSize(); i++){
        char * data = dbr_data.getData(i);
        std::string idStr = SSTR(i);
        size_t elementSize = 0;
        std::string header;
        if(par.profileMode == Parameters::PROFILE_MODE_HMM){
            parseHMM(data, &header, profileBuffer, &elementSize, idStr.c_str(), subMat);
        } else if(par.profileMode == Parameters::PROFILE_MODE_HMM) {
            parsePSSM(data, profileBuffer, &elementSize, idStr.c_str(), subMat);
        } else{
            Debug(Debug::ERROR) << "Wrong profile mode.\n";
            EXIT(EXIT_FAILURE);
        }
        ffindex_insert_memory(data_file_hdr,  index_file_hdr,     &offset_header,
                              (char*)header.c_str(),  header.length() , (char *)idStr.c_str());
        ffindex_insert_memory(data_file,  index_file,     &offset_profile,
                              profileBuffer,  elementSize , (char *)idStr.c_str());
    }
    entry_num = dbr_data.getSize();

    delete [] profileBuffer;
    delete subMat;
    dbr_data.close();
    fclose(data_file);
    fclose(data_file_hdr);

    /* Sort the index entries and write back */
    fclose(index_file);
    sortIndex(index_filename, entry_num);
    fclose(index_file_hdr);
    sortIndex(index_filename_hdr, entry_num);
    Debug(Debug::WARNING) << "Done.\n";

}

int sortIndex(char *index_filename, size_t entry_num) {
    FILE * index_file = fopen(index_filename, "r+");
    ffindex_index_t* index = ffindex_index_parse(index_file, entry_num);
    if(index == NULL)
    {
        perror("ffindex_index_parse failed");
        EXIT(EXIT_FAILURE);
    }
    fclose(index_file);
    ffindex_sort_index_file(index);
    index_file = fopen(index_filename, "w");
    if(index_file == NULL) { perror(index_filename); return EXIT_FAILURE; }
    return ffindex_write(index, index_file);
}

int openFile(const char *filename, FILE **pFILE) {
    struct stat st;
    if(stat(filename, &st) == 0) { errno = EEXIST; perror(filename); return EXIT_FAILURE; }
    *pFILE = fopen(filename, "w+");
    if(filename == NULL) { perror(filename); return EXIT_FAILURE; }
}

