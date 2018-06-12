#include "DistanceCalculator.h"
#include "Util.h"
#include "Parameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "QueryMatcher.h"

#include <string>
#include <vector>
#include <sys/time.h>
#include <NucleotideMatrix.h>
#include <ReducedMatrix.h>

#ifdef OPENMP
#include <omp.h>
#endif

void setAlignByKmerDefaults(Parameters *p) {
    p->kmerSize = 4;
    p->alphabetSize = 11;
}


float parseKmertoSeqidLib(std::string scoreFile, double targetSeqid, double targetCov, double targetPrecision){
    std::stringstream in(scoreFile);
    std::string line;
    // find closest lower seq. id in a grid of size 5
    int intTargetSeqid = static_cast<int>((targetSeqid+0.0001)*100);
    int seqIdRest = (intTargetSeqid % 5);
    targetSeqid = static_cast<float>(intTargetSeqid- seqIdRest)/100;
    // find closest lower cov. id in a grid of size 10
    targetCov = static_cast<float>(static_cast<int>((targetCov + 0.0001)* 10))/10;
    while (std::getline(in, line)) {
        std::vector<std::string> values = Util::split(line, " ");
        float cov = strtod(values[0].c_str(),NULL);
        float seqid = strtod(values[1].c_str(),NULL);
        float scorePerCol = strtod(values[2].c_str(),NULL);
        float precision = strtod(values[3].c_str(),NULL);
        if(MathUtil::AreSame(cov, targetCov)
           && MathUtil::AreSame(seqid, targetSeqid)
           && precision >= targetPrecision){
            return scorePerCol;
        }
    }
    Debug(Debug::WARNING) << "Could not find any score per column for "
            "cov="<< targetCov << " seq.id=" << targetSeqid << "\n"
                                  "No hit will be filtered.\n";

    return 0;
}


int alignbykmer(int argc, const char **argv, const Command &command) {
    Debug(Debug::INFO) << "Rescore diagonals.\n";
    struct timeval start, end;
    gettimeofday(&start, NULL);
    Parameters &par = Parameters::getInstance();
    setAlignByKmerDefaults(&par);
    par.parseParameters(argc, argv, command, 4);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> *qdbr = NULL;

    Debug(Debug::INFO) << "Query  file: " << par.db1 << "\n";
    qdbr = new DBReader<unsigned int>(par.db1.c_str(), (par.db1 + ".index").c_str());
    qdbr->open(DBReader<unsigned int>::NOSORT);
    int querySeqType  =  qdbr->getDbtype();

    DBReader<unsigned int> *tdbr = NULL;
    BaseMatrix *subMat;
    if (querySeqType == Sequence::NUCLEOTIDES) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.c_str(), 1.0, 0.0);
    }else {
        if (par.alphabetSize == 21) {
            subMat = new SubstitutionMatrix(par.scoringMatrixFile.c_str(), 2.0, 0.0);
        } else {
            SubstitutionMatrix sMat(par.scoringMatrixFile.c_str(), 2.0, 0.0);
            subMat = new ReducedMatrix(sMat.probMatrix, sMat.subMatrixPseudoCounts, par.alphabetSize, 2.0);
        }
    }

    qdbr->readMmapedDataInMemory();
    Debug(Debug::INFO) << "Target  file: " << par.db2 << "\n";
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), (par.db2 + ".index").c_str());
        tdbr->open(DBReader<unsigned int>::NOSORT);
        tdbr->readMmapedDataInMemory();
    }
    float scorePerColThr = 0.0;
    if(par.filterHits){
        if(par.rescoreMode == Parameters::RESCORE_MODE_HAMMING){
            Debug(Debug::WARNING) << "HAMMING distance can not be used to filter hits. Using --rescore-mode 1\n";
            par.rescoreMode = Parameters::RESCORE_MODE_SUBSTITUTION;
        }

        std::string libraryString = "";
        scorePerColThr = parseKmertoSeqidLib(libraryString, par.seqIdThr, par.covThr, 0.99);
    }
    Debug(Debug::WARNING) << "Prefilter database: " << par.db3 << "\n";
    DBReader<unsigned int> dbr_res(par.db3.c_str(), std::string(par.db3 + ".index").c_str());
    dbr_res.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    Debug(Debug::INFO) << "Result database: " << par.db4 << "\n";
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads);
    resultWriter.open();

    struct KmerPos{
        unsigned int ij;
        unsigned short i;
        unsigned short j;
        KmerPos(unsigned int ij, unsigned short i, unsigned short j):
                ij(ij), i(i), j(j) {}
        KmerPos(){};
        // need for sorting the results
        static bool compareKmerPos (const KmerPos &first, const KmerPos &second){
            //return (first.eval < second.eval);
            if(first.i < second.i )
                return true;
            if(second.i < first.i )
                return false;
            if(first.ij < second.ij )
                return true;
            if(second.ij < first.ij )
                return false;
            if(first.j <  second.j )
                return true;
            if(second.j < first.j )
                return false;


            return false;
        }
    };


    size_t totalMemory = Util::getTotalSystemMemory();
    size_t flushSize = 100000000;
    if(totalMemory > dbr_res.getDataSize()){
        flushSize = dbr_res.getSize();
    }
    size_t iterations = static_cast<int>(ceil(static_cast<double>(dbr_res.getSize()) / static_cast<double>(flushSize)));
    for (size_t i = 0; i < iterations; i++) {
        size_t start = (i * flushSize);
        size_t bucketSize = std::min(dbr_res.getSize() - (i * flushSize), flushSize);
#pragma omp parallel
        {
            Sequence query(par.maxSeqLen, querySeqType, subMat, par.kmerSize, true, false);
            Sequence target(par.maxSeqLen, querySeqType, subMat, par.kmerSize, true, false);
            size_t lookupSize = MathUtil::ipow<size_t>(par.alphabetSize, par.kmerSize);
            unsigned short * queryPosLookup = new unsigned short[lookupSize];
            memset(queryPosLookup, 255, lookupSize * sizeof(unsigned short) );
            unsigned short * diagonalCount = new unsigned short[lookupSize];
            memset(diagonalCount, 0, lookupSize * sizeof(unsigned short) );
            Indexer idxer(subMat->alphabetSize, par.kmerSize);
            std::vector<KmerPos> kmerPosVec;
#pragma omp for schedule(dynamic, 1)
            for (size_t id = start; id < (start + bucketSize); id++) {
                Debug::printProgress(id);
                char buffer[1024 + 32768];
                std::string prefResultsOutString;
                prefResultsOutString.reserve(1000000);
                unsigned int thread_idx = 0;
#ifdef OPENMP
                thread_idx = (unsigned int) omp_get_thread_num();
#endif
                char *data = dbr_res.getData(id);
                unsigned int queryId = qdbr->getId(dbr_res.getDbKey(id));
                char *querySeq = qdbr->getData(queryId);
                query.mapSequence(id, queryId, querySeq);

                while(query.hasNextKmer()) {
                    const int *kmer = query.nextKmer();
                    unsigned short kmerIdx = idxer.int2index(kmer);
                    unsigned short pos = query.getCurrentPosition();
                    if(queryPosLookup[kmerIdx] == USHRT_MAX){
                        queryPosLookup[kmerIdx] = pos;
                    }
                }


                resultWriter.writeStart(thread_idx);

                while (*data != '\0') {
                    // DB key of the db sequence
                    char dbKeyBuffer[255 + 1];
                    Util::parseKey(data, dbKeyBuffer);
                    const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                    unsigned int targetId = tdbr->getId(dbKey);
                    char *targetSeq = tdbr->getData(targetId);
                    const bool isIdentity = (queryId == targetId && (par.includeIdentity || sameDB)) ? true : false;
                    target.mapSequence(targetId, dbKey, targetSeq);
                    while(target.hasNextKmer()) {
                        const int *kmer = target.nextKmer();
                        unsigned short kmerIdx = idxer.int2index(kmer);
                        if(queryPosLookup[kmerIdx]!=USHRT_MAX){
                            unsigned short pos_j = target.getCurrentPosition();
                            unsigned short pos_i = queryPosLookup[kmerIdx];
                            unsigned short ij = static_cast<unsigned short>(pos_i) -
                                                static_cast<unsigned short>(pos_j);
                            kmerPosVec.emplace_back(ij, pos_i, pos_j);
                        }
                    }

//                    for(size_t kmerIdx = 0; kmerIdx < kmerPosVec.size(); kmerIdx++) {
//                        diagonalCount[kmerPosVec[kmerIdx].ij]++;
//                    }
//                    size_t writePos = 0;
//                    for(size_t kmerIdx = 0; kmerIdx < kmerPosVec.size(); kmerIdx++) {
//                        if(diagonalCount[kmerPosVec[kmerIdx].ij] > 1){
//                            kmerPosVec[writePos].ij = kmerPosVec[kmerIdx].ij;
//                            kmerPosVec[writePos].i = kmerPosVec[kmerIdx].i;
//                            kmerPosVec[writePos].j = kmerPosVec[kmerIdx].j;
//                            writePos++;
//                        }else{
//                            printf("%d\n", kmerPosVec[kmerIdx].ij);
//                        }
//                    }

                    std::sort(kmerPosVec.begin(), kmerPosVec.end(), KmerPos::compareKmerPos);
                    unsigned int kmerCnt = 0;
                    unsigned short min_i = USHRT_MAX;
                    unsigned short max_i = 0;
                    unsigned short min_j = USHRT_MAX;
                    unsigned short max_j = 0;

                    if(kmerPosVec.size() > 1){
                        unsigned int prevDiagonal = UINT_MAX;
                        unsigned short prev_i = 0;
                        unsigned short prev_j = 0;


                        for(size_t kmerIdx = 0; kmerIdx < kmerPosVec.size(); kmerIdx++){
                            unsigned int currDiagonal = static_cast<unsigned int>(kmerPosVec[kmerIdx].i)- static_cast<unsigned int>(kmerPosVec[kmerIdx].j);
                            unsigned short curr_i = kmerPosVec[kmerIdx].i;
                            unsigned short curr_j = kmerPosVec[kmerIdx].j;
                            unsigned int nextDiagonal = UINT_MAX;
                            if(kmerIdx <  kmerPosVec.size()-1 ) {
                                nextDiagonal= static_cast<unsigned int>(kmerPosVec[kmerIdx+1].i)- static_cast<unsigned int>(kmerPosVec[kmerIdx+1].j);
                            }
                            // skip noise
                            if(nextDiagonal==prevDiagonal && currDiagonal != prevDiagonal){
                                continue;
                            }


                            if((nextDiagonal == currDiagonal || prevDiagonal == currDiagonal)
                               && prev_i <= curr_i && prev_j <= curr_j  ){
//                std::cout << "ij: " << kmerPosVec[kmerIdx].ij << " i: " << kmerPosVec[kmerIdx].i
//                      << " j:" << kmerPosVec[kmerIdx].j << " Diag: " <<  kmerPosVec[kmerIdx].i - kmerPosVec[kmerIdx].j << std::endl;
                                min_i = std::min(min_i, curr_i);
                                max_i = std::max(max_i, curr_i);
                                min_j = std::min(min_j, curr_j);
                                max_j = std::max(max_j, curr_j);
                                kmerCnt++;
                            }
                            if(prev_i!=prev_j){
                                prevDiagonal = currDiagonal;
                            }
                            prev_i=curr_i;
                            prev_j=curr_j;
                        }
                    }

                    float queryCov = SmithWaterman::computeCov(min_i, max_i,  query.L);
                    float targetCov = SmithWaterman::computeCov(min_j, max_j, target.L);
                    int alnLen = std::min(max_i - min_i, max_j - min_j);
                    if(min_i == USHRT_MAX || min_j == USHRT_MAX){
                        queryCov = 0.0f;
                        targetCov = 0.0f;
                        alnLen = 0;
                    }

                    const float seqId = 1.0;
                    const float evalue = 0.0;
                    // query/target cov mode
                    const bool hasCov = Util::hasCoverage(par.covThr, par.covMode, queryCov, targetCov);
                    // --min-seq-id
                    const bool hasSeqId = seqId >= (par.seqIdThr - std::numeric_limits<float>::epsilon());
                    const bool hasEvalue = (evalue <= par.evalThr);
                    // --filter-hits
                    const bool hasToFilter = (par.filterHits == true  && kmerCnt >= scorePerColThr);
                    if (isIdentity || hasToFilter || (hasCov && hasSeqId && hasEvalue)) {
                        Matcher::result_t result = Matcher::result_t(dbKey, kmerCnt, queryCov, targetCov, 1.0, 0.0,
                                                                     alnLen,
                                                                     min_i, max_i, query.L, min_j, max_j,
                                                                     target.L, std::string(""));
                        size_t len = Matcher::resultToBuffer(buffer, result, true, false);
                        resultWriter.writeAdd(buffer, len, thread_idx);
                    }
                    kmerPosVec.clear();
                    data = Util::skipLine(data);
                }
                memset(queryPosLookup, 255, lookupSize * sizeof(unsigned short) );
                memset(diagonalCount, 0, lookupSize * sizeof(unsigned short) );
                resultWriter.writeEnd(qdbr->getDbKey(queryId), thread_idx, true);
            }
        }
        dbr_res.remapData();
    }
    Debug(Debug::INFO) << "Done." << "\n";
    dbr_res.close();
    resultWriter.close();
    qdbr->close();
    delete qdbr;
    delete subMat;
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }
    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::INFO) << "Time for diagonal calculation: " << (sec / 3600) << " h "
                       << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";
    return EXIT_SUCCESS;
}


