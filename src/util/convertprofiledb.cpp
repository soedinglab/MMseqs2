// Written by Martin Steinegger martin.steinegger@mpibpc.mpg.de
//
// Converts PSSM or HHM to MMseqs profile format.
// MMseqs just stores the position specific score in 1 byte
//
#include <unistd.h>
#include <limits.h>
#include <stdlib.h>
#include "SubstitutionMatrix.h"
#include "Parameters.h"
#include "Sequence.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"

void parsePSSM(char *data, std::string * sequence, char *profileBuffer, size_t *size, BaseMatrix *subMat) {
    // go to read in position
    for (size_t i = 0; i < 2; i++) {
        data = Util::skipLine(data);
    }

    // read aa_index
    char *words[22];
    Util::getWordsOfLine(data, words, 22);

    // A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    // header line
    int aa_index[20];
    for (size_t i = 0; i < 20; i++) {
        aa_index[i] = subMat->aa2int[(int) words[i][0]];
    }

    data = Util::skipLine(data);
    size_t curr_pos = 0;
    while (data[0] != '\n') {
        Util::getWordsOfLine(data, words, 22);
        sequence->push_back(words[1][0]);
        for (size_t i = 0; i < 20; i++) {
            size_t writePos = curr_pos + aa_index[i];
            profileBuffer[writePos] = atoi(words[2 + i]);
            // shifted score by -128 to avoid \0
            profileBuffer[writePos] = (profileBuffer[writePos] ^ 0x80);
        }

        data = Util::skipLine(data);
        curr_pos += 20;
    }
    sequence->push_back('\n');

    *size = curr_pos;
}

/////////////////////////////////////////////////////////////////////////////////////
void parseHMMer(char *data, std::string *sequence, std::string *header, char *profileBuffer, size_t *size, const char *id, BaseMatrix *subMat) {
    size_t l = 0;
    // find name tag
    while (data[0] != 'N' || data[1] != 'A' || data[2] != 'M' || data[3] != 'E') {
        data = Util::skipLine(data);
    }
    float pb[20];
    pb[0] = 0.0787945;             /* A */
    pb[1] = 0.0151600;             /* C */
    pb[2] = 0.0535222;             /* D */
    pb[3] = 0.0668298;             /* E */
    pb[4] = 0.0397062;             /* F */
    pb[5] = 0.0695071;             /* G */
    pb[6] = 0.0229198;             /* H */
    pb[7] = 0.0590092;             /* I */
    pb[8] = 0.0594422;             /* K */
    pb[9] = 0.0963728;             /* L */
    pb[10]= 0.0237718;             /* M */
    pb[11]= 0.0414386;             /* N */
    pb[12]= 0.0482904;             /* P */
    pb[13]= 0.0395639;             /* Q */
    pb[14]= 0.0540978;             /* R */
    pb[15]= 0.0683364;             /* S */
    pb[16]= 0.0540687;             /* T */
    pb[17]= 0.0673417;             /* V */
    pb[18]= 0.0114135;             /* W */
    pb[19]= 0.0304133;             /* Y */

    // parse NAME entry
    const char *startData = data;
    data = Util::skipLine(data);
    const char *endData = data;
    header->append(startData + 6, endData - (startData + 6));



    // find beginning of profile information
    while (data[0] != 'H' || data[1] != 'M' || data[2] != 'M') {
        data = Util::skipLine(data);
    }
    // go to readin position
    for (int i = 0; i < 5; i++)
        data = Util::skipLine(data);

    //ammino acids are ordered in HMM
    char *words[23];
    size_t curr_pos = 0;
    while (data[0] != '/' && data[1] != '/') {
        Util::getWordsOfLine(data, words, 23);
        char aa = std::toupper(words[22][0]);
        sequence->push_back(aa);
        for (size_t aa_num = 0; aa_num < 20; aa_num++) {
            // entry: 0.0 probability
            if (words[aa_num + 1][0] == '*') {
                Debug(Debug::ERROR) << "ERROR: 0 PROBABILITY FOR " << id << ".hhm AT " << l << "," << aa_num << "\n";
                profileBuffer[curr_pos] = (char) -1;
            }
                // 0 entry: 1.0 probability
            else if (words[aa_num + 1][0] == '0') {
                // integer number entry: 0.0 < probability < 1.0
                float score = MathUtil::flog2(1.0f / subMat->getBackgroundProb(aa_num)) * subMat->getBitFactor();
                profileBuffer[curr_pos] = (char) floor(score + 0.5);
            } else {
                float entry = strtof(words[aa_num + 1],NULL);
                // back scaling from hhm
                // fprintf(fp, " %*.5f", fieldwidth, -logf(p)
                const float p =  (float) exp(-1.0 * entry);
                const float backProb = pb[aa_num];
                //TODO solve somehow this?!?
                const float bitFactor = subMat->getBitFactor();

                double score = MathUtil::flog2(p / backProb) * bitFactor;
                // rounding
                profileBuffer[curr_pos]  = static_cast<char>((score < 0.0) ? score - 0.5 : score + 0.5);
//                Debug(Debug::INFO) << aa_num << " " << subMat->int2aa[aa_num] << " " << profile_score[pos_in_profile] << " " << score << " " << entry << " " << p << " " << backProb << " " << bitFactor << std::endl;
            }
            // shifted score by -128 to avoid \0
            profileBuffer[curr_pos] = (profileBuffer[curr_pos] ^ 0x80);
            if (profileBuffer[curr_pos] == 0) {
                Debug(Debug::ERROR) << "ERROR: 0 PSSM score is too large at id: " << id << ".hhm, pos: " << curr_pos <<
                ", score:" <<
                (char) (profileBuffer[curr_pos] ^ 0x80) << "\n";
                EXIT(EXIT_FAILURE);
            }
            curr_pos++;
        }

        // go to next entry start and skip transitions
        for (int i = 0; i < 3; i++)
            data = Util::skipLine(data);
    }

    sequence->push_back('\n');

    // return size of buffer
    *size = curr_pos;
}



void parseHMM(char *data, std::string *sequence, std::string *header, char *profileBuffer, size_t *size, const char *id, BaseMatrix *subMat) {
    size_t l = 0;
    // find name tag
    while (data[0] != 'N' || data[1] != 'A' || data[2] != 'M' || data[3] != 'E') {
        data = Util::skipLine(data);
    }

    // parse NAME entry
    const char *startData = data;
    data = Util::skipLine(data);
    const char *endData = data;
    header->append(startData + 6, endData - (startData + 6));

    // >Consensus
    while (data[0] != '>') {
        data = Util::skipLine(data);
    }
    // skip over Cons. header
    data = Util::skipLine(data);
    // find first line after >Consensus that starts with a >
    while (data[0] != '>' ) {
        data = Util::skipLine(data);
    }
    data = Util::skipLine(data);
    char * seqStartPos = data;
    // copy sequence
    while (data[0] != '>' && data[0] != '#'  ) {
        data = Util::skipLine(data);
    }
    char * seqEndPos = data;
    size_t len = (seqEndPos - seqStartPos);
    for(size_t i = 0; i < len; i++){
        if(seqStartPos[i] != '\n')
            sequence->push_back(seqStartPos[i]);
    }
    sequence->push_back('\n');

    // find beginning of profile information
    while (data[0] != '#') {
        data = Util::skipLine(data);
    }

    // go to readin position
    for (int i = 0; i < 5; i++)
        data = Util::skipLine(data);

    //ammino acids are ordered in HMM
    char *words[22];
    size_t curr_pos = 0;
    while (data[0] != '/' && data[1] != '/') {
        Util::getWordsOfLine(data, words, 22);
        for (size_t aa_num = 0; aa_num < 20; aa_num++) {
            // entry: 0.0 probability
            if (words[aa_num + 2][0] == '*') {
                Debug(Debug::ERROR) << "ERROR: 0 PROBABILITY FOR " << id << ".hhm AT " << l << "," << aa_num << "\n";
                profileBuffer[curr_pos] = (char) -1;
            }
            // 0 entry: 1.0 probability
            else if (words[aa_num + 2][0] == '0') {
                // integer number entry: 0.0 < probability < 1.0
                float score = MathUtil::flog2(1.0f / subMat->getBackgroundProb(aa_num)) * subMat->getBitFactor();
                profileBuffer[curr_pos] = (char) floor(score + 0.5);
            } else {
                int entry = Util::fast_atoi(words[aa_num + 2]);
                // back scaling from hhm
                const float p = MathUtil::fpow2(-(entry / 1000.0f));
                const float backProb = subMat->getBackgroundProb(aa_num);
                //TODO solve somehow this?!?
                const float bitFactor = subMat->getBitFactor();

                double score = MathUtil::flog2(p / backProb) * bitFactor;
                // rounding
                profileBuffer[curr_pos]  = static_cast<char>((score < 0.0) ? score - 0.5 : score + 0.5);
//                Debug(Debug::INFO) << aa_num << " " << subMat->int2aa[aa_num] << " " << profile_score[pos_in_profile] << " " << score << " " << entry << " " << p << " " << backProb << " " << bitFactor << std::endl;
            }
            // shifted score by -128 to avoid \0
            profileBuffer[curr_pos] = (profileBuffer[curr_pos] ^ 0x80);
            if (profileBuffer[curr_pos] == 0) {
                Debug(Debug::ERROR) << "ERROR: 0 PSSM score is too large at id: " << id << ".hhm, pos: " << curr_pos <<
                ", score:" <<
                (char) (profileBuffer[curr_pos] ^ 0x80) << "\n";
                EXIT(EXIT_FAILURE);
            }
            curr_pos++;
        }

        // go to next entry start and skip transitions
        for (int i = 0; i < 3; i++)
            data = Util::skipLine(data);
    }

    // return size of buffer
    *size = curr_pos;
}

int convertprofiledb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    switch (par.profileMode) {
        case Parameters::PROFILE_MODE_HMM:
        case Parameters::PROFILE_MODE_PSSM:
        case Parameters::PROFILE_MODE_HMM3:
            break;
        default:
            Debug(Debug::ERROR) << "Unsupported profile mode " << par.profileMode << "\n";
            return EXIT_FAILURE;
    }

    DBReader<std::string> dataIn(par.db1.c_str(), par.db1Index.c_str());
    dataIn.open(DBReader<std::string>::NOSORT);

    DBWriter dataOut(par.db2.c_str(), par.db2Index.c_str());
    dataOut.open();

    DBWriter seqOut(std::string(par.db2 +"_seq").c_str(),std::string(par.db2 +"_seq.index").c_str());
    seqOut.open();

    std::string headerFileName(par.db2);
    headerFileName.append("_h");

    std::string headerIndexFileName(par.db2);
    headerIndexFileName.append("_h.index");

    DBWriter headerOut(headerFileName.c_str(), headerIndexFileName.c_str());
    headerOut.open();

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);

    unsigned int *lengths = dataIn.getSeqLens();
    unsigned int maxElementSize = 0;
    for (size_t i = 0; i < dataIn.getSize(); i++) {
        maxElementSize = std::max(lengths[i], maxElementSize);
    }
    std::string sequence;
    std::string header;

    Debug(Debug::INFO) << "Start converting profile to MMseqs profile.\n";
    char *profileBuffer = new char[maxElementSize * Sequence::PROFILE_AA_SIZE];
    for (size_t i = 0; i < dataIn.getSize(); i++) {
        char *data = dataIn.getData(i);
        std::string idStr = SSTR(i);
        size_t elementSize = 0;
        if (par.profileMode == Parameters::PROFILE_MODE_HMM) {
            parseHMM(data, &sequence, &header, profileBuffer, &elementSize, idStr.c_str(), &subMat);
        } else if (par.profileMode == Parameters::PROFILE_MODE_HMM3) {
            parseHMMer(data, &sequence, &header, profileBuffer, &elementSize, idStr.c_str(), &subMat);
        } else if (par.profileMode == Parameters::PROFILE_MODE_PSSM) {
            parsePSSM(data, &sequence, profileBuffer, &elementSize, &subMat);
            header.append(dataIn.getDbKey(i));
            header.append(" \n");
        }
        seqOut.writeData(sequence.c_str(), sequence.size(), (char *) idStr.c_str());
        dataOut.writeData(profileBuffer, elementSize, (char *) idStr.c_str());
        headerOut.writeData((char *) header.c_str(), header.length(), (char *) idStr.c_str());
        sequence.clear();
        header.clear();
    }
    delete[] profileBuffer;
    char *absHeaderFileName = realpath(headerFileName.c_str(), NULL);
    symlink(absHeaderFileName, std::string(par.db2 +"_seq_h").c_str());
    free(absHeaderFileName);
    char *absHeaderIndexFileName = realpath(headerIndexFileName.c_str(), NULL);
    symlink(absHeaderIndexFileName, std::string(par.db2 +"_seq_h.index").c_str());
    free(absHeaderIndexFileName);
    headerOut.close();
    dataOut.close();
    seqOut.close();
    dataIn.close();

    Debug(Debug::INFO) << "Done.\n";

    return EXIT_SUCCESS;
}
