#include <ProfileStates.h>
#include <PSSMCalculator.h>
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Util.h"
#include "Debug.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "itoa.h"

#ifdef OPENMP
#include <omp.h>
#endif

int profile2cs(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.alphabetSize = 8;
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_PROFILE);

    DBReader<unsigned int> profileReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    profileReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    int alphabetSize[] = {219, 255};
    for(size_t i = 0; i < 2; i++){
        std::string dbName = par.db2;
        std::string dbIndex = par.db2;
        if(i>=1){
            dbName +=   "." + SSTR(alphabetSize[i]);
            dbIndex +=  "." + SSTR(alphabetSize[i]);
        }
        dbIndex += ".index";
        DBWriter writer(dbName.c_str(), dbIndex.c_str(), par.threads, par.compressed, Parameters::DBTYPE_PROFILE_STATE_SEQ);
        writer.open();
        size_t alphSize = alphabetSize[i];
        size_t entries = profileReader.getSize();

        SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0f, 0.0);
        Debug::Progress progress(entries);

        Debug(Debug::INFO) << "Start converting profiles.\n";
#pragma omp parallel
        {
            Sequence seq(par.maxSeqLen, Parameters::DBTYPE_HMM_PROFILE, &subMat, 0, false, false, true);
            unsigned int thread_idx = 0;
#ifdef OPENMP
            thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
            std::string result;
            result.reserve(par.maxSeqLen + 1);
            ProfileStates ps(alphabetSize[i], subMat.pBack);

#pragma omp for schedule(dynamic, 1000)
            for (size_t i = 0; i < entries; ++i) {
                progress.updateProgress();
                result.clear();

                unsigned int key = profileReader.getDbKey(i);
                seq.mapSequence(i, key, profileReader.getData(i, thread_idx), profileReader.getSeqLen(i));
                if(alphSize == 219){
                    ps.discretizeCs219(seq.getProfile(), seq.L, result);
                }else {
                    ps.discretize(seq.getProfile(), seq.L, result);
                }
                
                //std::cout<<result.size()<<" vs "<<seq.L<<std::endl;
                
                // DEBUG: in case of pure state library, check when the wrong pure state has been chosen
                /*for (size_t k = 0; k<result.size();k++)
                {
                    if (subMat.subMatrix[ProfileStates::hh2mmseqsAAorder((int)result[k])][seq.consensus_sequence[k]]<0)
                    {
                        std::cout<<"Pos: "<<k<<", "<<subMat.num2aa[ProfileStates::hh2mmseqsAAorder((int)result[k])]<<"-"<<subMat.num2aa[seq.consensus_sequence[k]]<<"("<<subMat.subMatrix[ProfileStates::hh2mmseqsAAorder((int)result[k])][seq.consensus_sequence[k]]<<") \n";
                        for (size_t aa = 0; aa<20;aa++)
                        {
                            std::cout<<subMat.num2aa[aa]<<"\t";
                        }
                        std::cout<<"\n";
                        for (size_t aa = 0; aa<20;aa++)
                        {
                            std::cout<<(float)seq.profile[k * Sequence::PROFILE_AA_SIZE + aa]<<"\t";
                        }
                        std::cout<<"\n";std::cout<<"\n";
                        
                    }
                }
                std::cout<<std::endl;*/
                
                    
                result.push_back('\0'); // needed to avoid seq. len calculation problems (aa sequence has \n\0, PS \0\0)
                for(size_t i = 0; i < result.size() - 1; i++){ // do not overwrite last \0
//                std::cout << (int)result[i] << std::endl;
                    char val = result[i];
                    val += 1; // avoid null byte (needed for read in)
                    result[i] = val;
                }
                writer.writeData(result.c_str(), result.length(), key, thread_idx);
            }
        }
        writer.close();
    }

    profileReader.close();
    return EXIT_SUCCESS;
}
