#include <cstddef>
#include <stdio.h>
#include <Matcher.h>
#include <SubstitutionMatrix.h>
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Prefiltering.h"
#include "Util.h"
typedef std::map<unsigned int, std::string *>::iterator SwapIt;

void readAllKeysIntoMap(std::string datafile,
                        std::map<unsigned int, std::string *> &map) {
    char dbKey[255 + 1];
    std::ifstream resultFile(datafile);
    std::string line;
    while(std::getline(resultFile, line)){
        size_t i;
        for(i = 0; i < line.size() - 1; i++ ){
            if(line[i] != '\0') {
                break;
            }
        }
        Util::parseKey((char*)line.c_str() + i, dbKey);
        unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);

        if(dbKey[0]=='\0')
            continue;
        SwapIt it = map.find(key);
        if (it == map.end()) {
            map[key] = new std::string();
        }
    }
    resultFile.close();
}

void processSplit(DBReader<unsigned int> &dbr,
                  FILE * dataFile, std::map<unsigned int, std::string *> &map,
                  size_t startIndex, size_t domainSize) {
    std::string resultData;
    for (size_t i = startIndex; i < (startIndex + domainSize); i++) {
        char dbKey[255 + 1];
        std::string outerKey = SSTR(dbr.getDbKey(i));
        int c1;
        while(( c1=fgetc(dataFile)) != EOF && c1 != (int) '\0'){
            char file1Char = (char) c1;
            resultData.push_back(file1Char);
        }
        resultData.push_back('\0');
        char * data = (char *) resultData.c_str();
        while (*data != '\0') {
            // extract key from results (ids must be always at the first position)
            Util::parseKey(data, dbKey);
            unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
            std::string *entry = NULL;
            entry = map[key];
            // write data to map
            entry->append(outerKey);
            // next db key
            char *endPosOfId = data + Util::skipNoneWhitespace(data);
            data = Util::skipLine(data);
            entry->append(endPosOfId, data);
        }
        resultData.clear();
    }
}

int swapresults (int argc, const char * argv[]){
    std::string usage("Swaps results of ffindex database. A -> A, B, C to A->A, B->A, C->A \n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de).\n\n");
    usage.append("USAGE: <queryDb> <targetDb> <ffindexDB> <fastaDB> [ffindexHeaderDB]\n");


    Parameters par;
    par.parseParameters(argc, argv, usage, par.swapresults, 4);
    size_t splitSize = par.split;
    Debug(Debug::INFO) << "FFindex input file is " << par.db1 << "\n";
    std::pair<std::string, std::string> queryDb = Util::databaseNames(par.db1);
    std::pair<std::string, std::string> targetDb = Util::databaseNames(par.db2);
    std::pair<std::string, std::string> resultDb = Util::databaseNames(par.db3);
    DBReader<unsigned int> qdbr(queryDb.first.c_str(), queryDb.second.c_str());
    qdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> tdbr(targetDb.first.c_str(), targetDb.second.c_str());
    tdbr.open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> rdbr(resultDb.first.c_str(), resultDb.second.c_str());
    rdbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Debug(Debug::INFO) << "Start to swap results. Write to " << par.db2 << ".\n";
    size_t entries_num = 0;
    std::vector<std::pair<std::string, std::string> > filesToDelete;
    std::map<unsigned int, std::string *> swapMap;
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
    double * kmnByLen = new double[par.maxSeqLen];

    BlastScoreUtils::BlastStat stats = BlastScoreUtils::getAltschulStatsForMatrix(subMat.getMatrixName(), Matcher::GAP_OPEN, Matcher::GAP_EXTEND);
    for(int len = 0; len < par.maxSeqLen; len++){
        kmnByLen[len] = BlastScoreUtils::computeKmn(len, stats.K, stats.lambda, stats.alpha, stats.beta, qdbr.getAminoAcidDBSize(), qdbr.getSize());
    }


    // read all keys
    readAllKeysIntoMap(resultDb.first, swapMap);
    FILE * dataFile = fopen(resultDb.first.c_str(), "r");
    for (size_t split = 0; split < splitSize; split++) {
        // create and open db write
        // create splite file name
        std::string splitName(par.db4);
        splitName.append("_").append(SSTR(split));
        std::pair<std::string, std::string> splitNames = Util::databaseNames(splitName);
        filesToDelete.push_back(std::pair<std::string, std::string>(splitNames.first, splitNames.second));

        DBWriter splitWrite(splitNames.first.c_str(), splitNames.second.c_str(), 1);
        splitWrite.open();

        size_t startIndex = 0;
        size_t domainSize = 0;
        Util::decomposeDomain(rdbr.getSize(), split, splitSize, &startIndex, &domainSize);
        Debug(Debug::INFO) << "Process split " << split << " from " << startIndex << " to " << (startIndex + domainSize) << " ... ";

        processSplit(rdbr, dataFile, swapMap, startIndex, domainSize);
        // Write sorted results and delete memory
        char * wordCnt[255];
        for (SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
            std::string resultData;
            resultData = std::string(iterator->second->c_str());
            if(iterator->second->size() > 1){
                size_t cols = Util::getWordsOfLine((char*)iterator->second->c_str(), wordCnt, 254);
                if(Matcher::ALN_RES_WITH_OUT_BT_COL_CNT >= cols){
                    std::vector<Matcher::result_t> alnRes = Matcher::readAlignmentResults((char*)iterator->second->c_str());
                    std::stringstream swResultsSs;
                    for(size_t i = 0; i < alnRes.size(); i++){
                        Matcher::result_t res = alnRes[i];
                        //TODO I can not use score here because its bit score and not RAW
                        res.eval = BlastScoreUtils::computeEvalue(res.score, kmnByLen[res.dbLen], stats.lambda);
                        size_t qstart = res.qStartPos;
                        size_t qend    = res.qEndPos;
                        size_t qLen    = res.qLen;
                        res.qStartPos  = res.dbStartPos;
                        res.qEndPos    = res.dbEndPos;
                        res.qLen       = res.dbLen;
                        res.dbStartPos = qstart;
                        res.dbEndPos   = qend;
                        res.dbLen      = qLen;
                    }
                    std::sort(alnRes.begin(), alnRes.end(), Matcher::compareHits);
                    for(size_t i = 0; i < alnRes.size(); i++){
                        swResultsSs << Matcher::resultToString(alnRes[i], (Matcher::ALN_RES_WITH_BT_COL_CNT==cols));
                    }
                    resultData = swResultsSs.str();
                }
            }
            splitWrite.write((char *) resultData.c_str(), iterator->second->size(),
                             (char *) SSTR(iterator->first).c_str(), 0);
            entries_num++;
            // remove just the value (string *) not the keys
            // the keys are needed for the merging step later
            delete iterator->second;
            iterator->second = new std::string();
        }
        Debug(Debug::INFO) << "Done.\n";
        splitWrite.close();
    }
    qdbr.close();
    tdbr.close();
    rdbr.close();
    fclose(dataFile);
    delete [] kmnByLen;
    for (SwapIt iterator = swapMap.begin(); iterator != swapMap.end(); iterator++) {
        delete iterator->second;
    }
    // merge output of all swap splits
    Prefiltering::mergeOutput(par.db4, std::string(par.db4 + ".index"), filesToDelete);
    for (size_t i = 0; i < filesToDelete.size(); i++) {
        remove(filesToDelete[i].first.c_str());
        remove(filesToDelete[i].second.c_str());
    }

    return 0;
}
