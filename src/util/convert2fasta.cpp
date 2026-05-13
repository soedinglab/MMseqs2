/*
 * convert2fasta
 * written by Milot Mirdita <milot@mirdita.de>
 */

#include <cstring>
#include <cstdio>
#include <string>
#include <vector>

#include "Parameters.h"
#include "DBReader.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"

const char headerStart[] = {'>'};
const char newline[] = {'\n'};

int convert2fasta(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> db(par.db1.c_str(), par.db1Index.c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    db.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> db_header(par.hdr1.c_str(), par.hdr1Index.c_str(), 1, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    db_header.open(DBReader<unsigned int>::NOSORT);

    FILE* fastaFP = fopen(par.db2.c_str(), "w");
    if(fastaFP == NULL) {
        perror(par.db2.c_str());
        EXIT(EXIT_FAILURE);
    }

    // Check for combined primary+aux format
    const Sequence::SeqAuxInfo* auxInfo = Sequence::getAuxInfo(db.getDbtype());
    FILE* auxFastaFP = NULL;
    std::string auxPath;
    SubstitutionMatrix* subMat = NULL;
    if (auxInfo != NULL) {
        const std::string ext = ".fasta";
        if (par.db2.length() >= ext.length() &&
            par.db2.substr(par.db2.length() - ext.length()) == ext) {
            auxPath = par.db2.substr(0, par.db2.length() - ext.length()) + ".aux.fasta";
        } else {
            auxPath = par.db2 + ".aux";
        }
        auxFastaFP = fopen(auxPath.c_str(), "w");
        if (auxFastaFP == NULL) {
            perror(auxPath.c_str());
            fclose(fastaFP);
            EXIT(EXIT_FAILURE);
        }
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, 0.0);
    }

    DBReader<unsigned int>* from = &db;
    if(par.useHeaderFile) {
        from = &db_header;
    }

    Debug(Debug::INFO) << "Start writing file to " << par.db2 << "\n";
    std::vector<char> decodedBuf;
    for(size_t i = 0; i < from->getSize(); i++){
        unsigned int key = from->getDbKey(i);
        unsigned int headerKey = db_header.getId(key);
        const char* headerData = db_header.getData(headerKey, 0);
        const size_t headerLen = db_header.getEntryLen(headerKey);

        fwrite(headerStart, sizeof(char), 1, fastaFP);
        fwrite(headerData, sizeof(char), headerLen - 2, fastaFP);
        fwrite(newline, sizeof(char), 1, fastaFP);

        unsigned int bodyKey = db.getId(key);
        const char* bodyData = db.getData(bodyKey, 0);
        const size_t bodyLen = db.getEntryLen(bodyKey);
        size_t seqLen = bodyLen - 2;

        if (auxInfo != NULL) {
            if (decodedBuf.size() < seqLen) {
                decodedBuf.resize(seqLen);
            }
            for (size_t j = 0; j < seqLen; j++) {
                decodedBuf[j] = subMat->num2aa[auxInfo->primaryRemap[(unsigned char)bodyData[j]]];
            }
            fwrite(decodedBuf.data(), sizeof(char), seqLen, fastaFP);
            fwrite(newline, sizeof(char), 1, fastaFP);

            fwrite(headerStart, sizeof(char), 1, auxFastaFP);
            fwrite(headerData, sizeof(char), headerLen - 2, auxFastaFP);
            fwrite(newline, sizeof(char), 1, auxFastaFP);
            for (size_t j = 0; j < seqLen; j++) {
                decodedBuf[j] = subMat->num2aa[auxInfo->auxRemap[(unsigned char)bodyData[j]]];
            }
            fwrite(decodedBuf.data(), sizeof(char), seqLen, auxFastaFP);
            fwrite(newline, sizeof(char), 1, auxFastaFP);
        } else {
            fwrite(bodyData, sizeof(char), seqLen, fastaFP);
            fwrite(newline, sizeof(char), 1, fastaFP);
        }
    }
    if (fclose(fastaFP) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << par.db2 << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (auxFastaFP != NULL && fclose(auxFastaFP) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << auxPath << "\n";
        EXIT(EXIT_FAILURE);
    }
    delete subMat;
    db_header.close();
    db.close();

    return EXIT_SUCCESS;
}
