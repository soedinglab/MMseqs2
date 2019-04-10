// Written by Martin Steinegger martin.steinegger@mpibpc.mpg.de
//
// Converts PSSM or HHM to MMseqs profile format.
// MMseqs just stores the position specific score in 1 byte

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
#include "FileUtil.h"

#ifdef OPENMP
#include "omp.h"
#endif


void parseHMM(char *data, std::string *sequence, std::string *header, char *profileBuffer, size_t *size, unsigned int id, BaseMatrix *subMat) {
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
    while (strncmp(">Consensus", data, 10) != 0) {
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
    const char *words[22];
    float probs[20];
    int seq_pos = 0;
    size_t curr_pos = 0;
    while (data[0] != '/' && data[1] != '/') {
        Util::getWordsOfLine(data, words, 22);
        for (size_t aa_num = 0; aa_num < 20; aa_num++) {
            // entry: 0.0 probability
            if (words[aa_num + 2][0] == '*') {
                probs[aa_num] = 0.0;
            }
            // 0 entry: 1.0 probability
            else if (words[aa_num + 2][0] == '0') {
                // integer number entry: 0.0 < probability < 1.0
                probs[aa_num] = 1.0;
            } else {
                int entry = Util::fast_atoi<int>(words[aa_num + 2]);
                // back scaling from hhm
                const float p = MathUtil::fpow2(-(entry / 1000.0f));
                probs[aa_num] = p;
//                const float backProb = subMat->getBackgroundProb(aa_num);
//                float score = MathUtil::flog2(p / backProb) * Sequence::PROFILE_SCALING;
//                float truncPssmVal =  std::min(score, 127.0f);
//                truncPssmVal       =  std::max(-128.0f, truncPssmVal);
                // rounding
//                profileBuffer[curr_pos]  = static_cast<char>((truncPssmVal < 0.0) ? truncPssmVal - 0.5 : truncPssmVal + 0.5);
//                Debug(Debug::INFO) << aa_num << " " << subMat->int2aa[aa_num] << " " << profile_score[pos_in_profile] << " " << score << " " << entry << " " << p << " " << backProb << " " << bitFactor << std::endl;
            }
            // shifted score by -128 to avoid \0
            profileBuffer[curr_pos] = Sequence::scoreMask(probs[aa_num]);

            if (profileBuffer[curr_pos] == 0) {
                Debug(Debug::ERROR) << "PSSM score of 0 is too large at id: " << id << ".hhm, pos: " << curr_pos <<
                ", score:" <<
                (char) (profileBuffer[curr_pos] ^ 0x80) << "\n";
                EXIT(EXIT_FAILURE);
            }
            curr_pos++;
        }

        float maxw = 0.0;
        int maxa = 21;
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; ++aa) {
            float prob = probs[aa];
            const float backProb = subMat->getBackgroundProb(aa);
            if (prob - backProb > maxw) {
                maxw = prob - backProb;
                maxa = aa;
            }
        }
        // write query, consensus and neff
        profileBuffer[curr_pos] = static_cast<char>(subMat->aa2int[(int)sequence->at(seq_pos)]);
        curr_pos++;
        profileBuffer[curr_pos] = maxa;
        curr_pos++;
        Util::getWordsOfLine(data, words, 22);
        int entry = Util::fast_atoi<int>(words[7]); // NEFF value
        const float neff = static_cast<float>(entry) / 1000.0f;
        profileBuffer[curr_pos] = MathUtil::convertNeffToChar(neff);
        curr_pos++;
        seq_pos++;
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

    DBReader<std::string> dataIn(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    dataIn.open(DBReader<std::string>::NOSORT);

    DBWriter dataOut(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_HMM_PROFILE);
    dataOut.open();

    DBWriter seqOut(std::string(par.db2 +"_seq").c_str(),std::string(par.db2 +"_seq.index").c_str(), par.threads, par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    seqOut.open();

    DBWriter headerOut(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headerOut.open();

    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);

    unsigned int *lengths = dataIn.getSeqLens();
    unsigned int maxElementSize = 0;
    for (size_t i = 0; i < dataIn.getSize(); i++) {
        maxElementSize = std::max(lengths[i], maxElementSize);
    }

    Debug(Debug::INFO) << "Start converting profiles to MMseqs2 profiles.\n";
    #pragma omp parallel
    {
        char *profileBuffer = new char[maxElementSize * Sequence::PROFILE_READIN_SIZE];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif

#pragma omp for  schedule(dynamic, 100)
        for (size_t i = 0; i < dataIn.getSize(); i++) {
            char *data = dataIn.getData(i, thread_idx);

            std::string sequence;
            std::string header;
            size_t elementSize = 0;

            parseHMM(data, &sequence, &header, profileBuffer, &elementSize, i, &subMat);

            seqOut.writeData(sequence.c_str(), sequence.size(), i, thread_idx);
            dataOut.writeData(profileBuffer, elementSize, i, thread_idx);
            headerOut.writeData((char *) header.c_str(), header.length(), i, thread_idx);
        }
        delete[] profileBuffer;
    }
    headerOut.close(true);
    dataOut.close(true);
    seqOut.close(true);

    std::string base = FileUtil::baseName(par.db2 + "_seq_h");
    FileUtil::symlinkAlias(par.hdr2, base);
    FileUtil::symlinkAlias(par.hdr2Index, base + ".index");

    dataIn.close();



    return EXIT_SUCCESS;
}
