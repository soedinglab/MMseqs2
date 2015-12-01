// Written by Martin Steinegger martin.steinegger@campus.lmu.de
//
// Computes PSSM from clustering or alignment result
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
#include <DBWriter.h>
#include <Log.h>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include "DBReader.h"
#include "Util.h"
#include "Debug.h"

int runResult2Profile(std::string queryDb, std::string targetDb, std::string resultDb, std::string outDb,
                      std::string subMatPath, int cpu, bool aaBiasCorrection) {
    int err = EXIT_SUCCESS;


    bool sameDatabase = true;
    DBReader * qdbr  = new DBReader(queryDb.c_str(), (queryDb + ".index").c_str());
    qdbr->open(DBReader::NOSORT);
    DBReader * tdbr = NULL;
    tdbr = qdbr;
    if(queryDb.compare(targetDb) != 0){
        sameDatabase = false;
        tdbr = new DBReader(targetDb.c_str(), (targetDb + ".index").c_str());
        tdbr->open(DBReader::NOSORT);
    }

    DBReader resultdbr(resultDb.c_str(), (resultDb + ".index").c_str());
    resultdbr.open(DBReader::NOSORT);

    DBWriter::errorIfFileExist(outDb.c_str());
    DBWriter::errorIfFileExist((std::string(outDb)+".index").c_str());
    DBWriter dbWriter(outDb.c_str(), (std::string(outDb)+".index").c_str(), cpu, DBWriter::BINARY_MODE);
    dbWriter.open();

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
    for(size_t i = 0; i < resultdbr.getSize(); i++) {
        char * data = resultdbr.getData(i);
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
    SubstitutionMatrix subMat(subMatPath.c_str(), 2.0, -0.2);
    Debug(Debug::WARNING) << "Start computing profiles.\n";
#pragma omp parallel
    {
        Matcher aligner(maxSeqLen, &subMat, tdbr->getAminoAcidDBSize(), tdbr->getSize(), aaBiasCorrection);
        MultipleAlignment msaAligner(maxSeqLen, maxSetSize, &subMat, &aligner);
        PSSMCalculator pssmCalculator(&subMat, maxSeqLen);
        Sequence centerSeq(maxSeqLen, subMat.aa2int, subMat.int2aa, Sequence::AMINO_ACIDS, 0, false);
        Sequence **dbSeqs = new Sequence *[maxSetSize];
        for (size_t i = 0; i < maxSetSize; i++) {
            dbSeqs[i] = new Sequence(maxSeqLen, subMat.aa2int, subMat.int2aa, Sequence::AMINO_ACIDS, 0, false);
        }
#pragma omp for schedule(dynamic, 1000)
        for (size_t id = 0; id < resultdbr.getSize(); id++) {
            Log::printProgress(id);

            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif

            char dbKey[255 + 1];
            char *data = resultdbr.getData(id);
            std::string queryId = resultdbr.getDbKey(id).c_str();
            char *seqData = qdbr->getDataByDBKey(queryId.c_str());
            centerSeq.mapSequence(0, (char*)queryId.c_str(), seqData);
            std::vector<Sequence *> seqSet;
            size_t pos = 0;

            while (*data != '\0') {
                Util::parseKey(data, dbKey);
                if ( (strcmp(dbKey, queryId.c_str()) == 0 && sameDatabase) == false ) {
                    char *dbSeqData = tdbr->getDataByDBKey(dbKey);
                    dbSeqs[pos]->mapSequence(0, dbKey, dbSeqData);
                    seqSet.push_back(dbSeqs[pos]);
                    pos++;
                }
                data = Util::skipLine(data);
            }
            //parseHMM(data, profileBuffer, &elementSize, id, subMat);
            MultipleAlignment::MSAResult res = msaAligner.computeMSA(&centerSeq, seqSet, true);
            //MultipleAlignment::print(res);

            char *pssmData = (char *) pssmCalculator.computePSSMFromMSA(res.setSize, res.centerLength, res.msaSequence);
            size_t pssmDataSize = res.centerLength * Sequence::PROFILE_AA_SIZE * sizeof(char);
            for (size_t i = 0; i < pssmDataSize; i++) {
                pssmData[i] = pssmData[i] ^ 0x80;
            }
            //pssm.printProfile(res.centerLength);
            //pssmCalculator.printPSSM(res.centerLength);

            dbWriter.write(pssmData, pssmDataSize, (char*)queryId.c_str(), thread_idx);
            seqSet.clear();
        }
        // clean memeory
        for (size_t i = 0; i < maxSetSize; i++) {
            delete dbSeqs[i];
        }
        delete[] dbSeqs;
    }
    // close reader
    resultdbr.close();
    qdbr->close();
    delete qdbr;
    if(sameDatabase == false) {
        tdbr->close();
        delete tdbr;
    }

    dbWriter.close();
    Debug(Debug::WARNING) << "\nDone.\n";
    return err;
}


int result2profile(int argn, const char **argv)
{
    std::string usage;
    usage.append("Converts a ffindex profile database to ffindex.\n");
    usage.append("USAGE: <queryDB> <targetDB> <resultDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Martin Steinegger <martin.steinegger@campus.lmu.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.createprofiledb, 4);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    return runResult2Profile(par.db1,
                             par.db2,
                             par.db3,
                             par.db4,
                             par.scoringMatrixFile,
                             par.threads,
                             par.compBiasCorrection);
}