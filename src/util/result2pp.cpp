// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include "Alignment.h"
#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBWriter.h"
#include "DBReader.h"
#include "HeaderSummarizer.h"
#include "CompressedA3M.h"
#include "Debug.h"
#include "Util.h"
#include "ProfileStates.h"
#include "MathUtil.h"
#include "SubstitutionMatrix.h"
#include <string>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif


float computeNeff(float neffA, float maxNeffA, float neffB, float maxNeffB, float avgNewNeff) {
    
    float w = (neffA + neffB) / (maxNeffA + maxNeffB);
    //std::cout<<"COmputing new neff:"<<neffA<<","<<maxNeffA<<","<<neffB<<","<<maxNeffB<<","<<avgNewNeff<<","<<w<<","<<std::endl;
    return avgNewNeff + 1 - exp(log(avgNewNeff) * (1-w));
     
}


int computeProfileProfile(Parameters &par,const std::string &outpath,
                          const size_t dbFrom, const size_t dbSize) {
    DBReader<unsigned int> *qDbr = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qDbr->open(DBReader<unsigned int>::NOSORT);
    DBReader<unsigned int> *tDbr = qDbr;
    bool sameDatabase = true;
    if (par.db1.compare(par.db2) != 0) {
        sameDatabase = false;
        tDbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tDbr->open(DBReader<unsigned int>::NOSORT);
    }

    DBReader<unsigned int> *resultReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter resultWriter(outpath.c_str(), (outpath + ".index").c_str(), par.threads, par.compressed, Parameters::DBTYPE_HMM_PROFILE);
    resultWriter.open();
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0f, 0.0f);

//#pragma omp parallel
    {
        Sequence queryProfile(par.maxSeqLen, qDbr->getDbtype(), &subMat, 0, false,
                              par.compBiasCorrection,false);
        Sequence targetProfile(par.maxSeqLen, tDbr->getDbtype(), &subMat, 0, false,
                               par.compBiasCorrection,false);
        float * outProfile = new float[(par.maxSeqLen + 1) * Sequence::PROFILE_AA_SIZE];
        float * neffM = new float[par.maxSeqLen + 1];
        std::string result;
        result.reserve((par.maxSeqLen + 1) * Sequence::PROFILE_READIN_SIZE);

        const char *entry[255];
        Debug::Progress progress(dbSize-dbFrom);

#pragma omp for schedule(dynamic, 10)
        for (size_t id = dbFrom; id < (dbFrom + dbSize); id++) {
            progress.updateProgress();
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = (unsigned int) omp_get_thread_num();
#endif
            char *results = resultReader->getData(id, thread_idx);
            char dbKey[255 + 1];
            // Get the sequence from the queryDB
            unsigned int queryKey = resultReader->getDbKey(id);
            char *queryData = qDbr->getDataByDBKey(queryKey, thread_idx);
            queryProfile.mapSequence(id, queryKey, queryData, resultReader->getSeqLen(queryKey));
            const float * qProfile =  queryProfile.getProfile();
	    /*
            const size_t profile_row_size = queryProfile.profile_row_size;
            // init outProfile with query Probs
            for(int l = 0; l < queryProfile.L; l++) {
                for (size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    outProfile[l*Sequence::PROFILE_AA_SIZE + aa_num] = qProfile[l * profile_row_size + aa_num];
                }
            }
            */

            float maxNeffQ = 0;
            for (int pos = 0; pos<queryProfile.L; pos++)
            {
                maxNeffQ = std::max(maxNeffQ,queryProfile.neffM[pos]);
                neffM[pos] = queryProfile.neffM[pos];
            }
            
            
            
            memset(outProfile, 0, queryProfile.L * Sequence::PROFILE_AA_SIZE * sizeof(float));
            while (*results != '\0') {
                Util::parseKey(results, dbKey);
                const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                double evalue = 0.0;
                const size_t columns = Util::getWordsOfLine(results, entry, 255);
                // its an aln result
                if (columns > Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    evalue = strtod(entry[3], NULL);
                }else{
                    Debug(Debug::ERROR) << "Alignment must contain the alignment information. Compute the alignment with option -a.\n";
                    EXIT(EXIT_FAILURE);
                }
                
                // just add sequences if eval < thr. and if key is not the same as the query in case of sameDatabase
                if (evalue <= par.evalProfile && (key != queryKey || sameDatabase == false)) {
                    Matcher::result_t res = Matcher::parseAlignmentRecord(results);
                    const size_t edgeId = tDbr->getId(key);
                    char *dbSeqData = tDbr->getData(edgeId, thread_idx);
                    targetProfile.mapSequence(0, key, dbSeqData, tDbr->getSeqLen(edgeId));
                    const float * tProfile = targetProfile.getProfile();
                    size_t qPos = res.qStartPos;
                    size_t tPos = res.dbStartPos;
                    size_t aliLength = 0;
                    float avgEntropy = 0.0f;
                    float maxNeffT = 0;
                    for (int pos = 0; pos<targetProfile.L; pos++)
                    {
                        maxNeffT = std::max(maxNeffT,targetProfile.neffM[pos]);
                    } 
                    
                    for(size_t btPos = 0; btPos < res.backtrace.size(); btPos++){
                        aliLength++;
                        char letter = res.backtrace[btPos];
//                        std::cout << letter;

//                        float qNeff = queryProfile.neffM[qPos];
//                        float tNeff = targetProfile.neffM[tPos];
                        
                        for(size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                            //TODO do all alignment states contribute?
                            if(letter == 'M') {
                                float qProb = qProfile[qPos * Sequence::PROFILE_AA_SIZE + aa_num] * queryProfile.neffM[qPos];
                                float tProb = tProfile[tPos * Sequence::PROFILE_AA_SIZE + aa_num] * targetProfile.neffM[tPos];
                                float mixedProb = qProb + tProb;
                                outProfile[qPos * Sequence::PROFILE_AA_SIZE + aa_num] += mixedProb;
                                mixedProb /= queryProfile.neffM[qPos] + targetProfile.neffM[tPos];
                                avgEntropy += (mixedProb>0.0) ? -mixedProb * log(mixedProb):0.0;

                            }
                        }
      
              
                        if(letter == 'M'){
                            qPos++;
                            tPos++;
                        } else if (letter == 'I') {
                            ++qPos;
                            //backtrace.append("I");
                        } else if (letter == 'D'){
                            ++tPos;
                            //backtrace.append("D");
                        }

                    }
                    
                    

                   
                    // Normalize probability
                    for(int l = res.qStartPos; l < res.qEndPos; l++) {
                        MathUtil::NormalizeTo1(&outProfile[l * Sequence::PROFILE_AA_SIZE], Sequence::PROFILE_AA_SIZE);
                    }
                    
                    avgEntropy /= aliLength;
                    float avgNewNeff = exp(avgEntropy);
                    
                    // update the Neff of the merge between the target prof and the query prof
                    qPos = res.qStartPos;
                    tPos = res.dbStartPos;
                    for(size_t btPos = 0; btPos < res.backtrace.size(); btPos++){
                        char letter = res.backtrace[btPos];

//                        float qNeff = queryProfile.neffM[qPos];
//                        float tNeff = targetProfile.neffM[tPos];
                        
                        neffM[qPos] = computeNeff(queryProfile.neffM[qPos], maxNeffQ, targetProfile.neffM[tPos],maxNeffT,avgNewNeff);
                        
                        if(letter == 'M'){
                            qPos++;
                            tPos++;
                        } else if (letter == 'I') {
                            ++qPos;
                            //backtrace.append("I");
                        } else if (letter == 'D'){
                            ++tPos;
                            //backtrace.append("D");
                        }

                    }
                   
                }
                results = Util::skipLine(results);
            }
            /*
            float maxNewNeff = 0.0;
            float avgEntropy = 0.0;
            for(int l = 0; l < queryProfile.L; l++) {
                maxNewNeff = std::max(maxNewNeff,neffM[l]);
                for(size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    float mixedProb = outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num];
                    avgEntropy += -mixedProb * log(mixedProb);
                }
            }
            avgEntropy /= queryProfile.L;
            
            float avgNewNeff = exp(avgEntropy);
            
            // update the Neff of the merges of all target prof and the query prof
            for(int l = 0; l < queryProfile.L; l++) {
                //for(size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    neffM[l] = computeNeff(queryProfile.neffM[l], maxNeffQ, neffM[l],maxNewNeff,avgNewNeff);
                //}
            }
                    */
//            size_t pos = 0;
            std::string consensus(queryProfile.L,'X');
            result.clear();
            for(int l = 0; l < queryProfile.L; l++) {
                float maxProb = 0.0f;
                for(size_t aa_num = 0; aa_num < Sequence::PROFILE_AA_SIZE; aa_num++) {
                    result.push_back(Sequence::scoreMask(outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num]));
                    //std::cout<< outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num]<<"\t";
                    if (outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num] > maxProb)
                    {
                        consensus[l] = aa_num;
                        maxProb = outProfile[l * Sequence::PROFILE_AA_SIZE + aa_num];
                    }
                }
                //std::cout<<std::endl;
                
                // write query, consensus sequence and neffM
                result.push_back(queryProfile.numSequence[l]);
                result.push_back(consensus[l]);
                unsigned char neff = MathUtil::convertNeffToChar(neffM[l]);
                result.push_back(neff);
            }
            //std::cout<<"Query length:"<<queryProfile.L<<", res length:"<<result.size()<<", should be:"<<queryProfile.L*23<<std::endl;
            resultWriter.writeData(result.c_str(), result.size(), queryKey, thread_idx);
        }
        delete [] outProfile;
        delete [] neffM;
    }

    // cleanup
    resultWriter.close(true);

    resultReader->close();
    delete resultReader;

    if (!sameDatabase) {
        tDbr->close();
        delete tDbr;
    }

    qDbr->close();
    delete qDbr;


    return EXIT_SUCCESS;
}

int computeProfileProfile(Parameters &par) {


    DBReader<unsigned int> *resultReader = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader->open(DBReader<unsigned int>::NOSORT);
    size_t resultSize = resultReader->getSize();
    resultReader->close();
    delete resultReader;

    std::string outname = par.db4;

    return computeProfileProfile(par, outname, 0, resultSize);

}

int computeProfileProfile(Parameters &par,const unsigned int mpiRank, const unsigned int mpiNumProc) {
    DBReader<unsigned int> *qDbr = new DBReader<unsigned int>(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qDbr->open(DBReader<unsigned int>::NOSORT);

    size_t dbFrom = 0;
    size_t dbSize = 0;
    qDbr->decomposeDomainByAminoAcid(mpiRank, mpiNumProc, &dbFrom, &dbSize);
    qDbr->close();
    delete qDbr;


    Debug(Debug::INFO) << "Compute split from " << dbFrom << " to " << dbFrom+dbSize << "\n";

    std::string outname = par.db4;

    std::pair<std::string, std::string> tmpOutput = Util::createTmpFileNames(outname, "", mpiRank);
    int status = computeProfileProfile(par, tmpOutput.first, dbFrom, dbSize );


#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // master reduces results
    if(mpiRank == 0) {
        std::vector<std::pair<std::string, std::string> > splitFiles;
        for(unsigned int procs = 0; procs < mpiNumProc; procs++){
            std::pair<std::string, std::string> tmpFile = Util::createTmpFileNames(outname, "", procs);
            splitFiles.push_back(std::make_pair(tmpFile.first ,  tmpFile.first + ".index"));

        }
        DBWriter::mergeResults(outname , outname + ".index", splitFiles);
    }

    return status;
}

int result2pp(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    par.evalProfile = (par.evalThr < par.evalProfile) ? par.evalThr : par.evalProfile;
    std::vector<MMseqsParameter*>* params = command.params;
    par.printParameters(command.cmd, argc, argv, *params);
    MMseqsMPI::init(argc, argv);

    // never allow deletions
    par.allowDeletion = false;
    Debug(Debug::INFO) << "Compute profile.\n";

#ifdef HAVE_MPI
    int retCode = computeProfileProfile(par, MMseqsMPI::rank, MMseqsMPI::numProc);
#else
    int retCode = computeProfileProfile(par);
#endif

    return retCode;
}
