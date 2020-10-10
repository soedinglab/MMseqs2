
#include "ProfileStates.h"
#include "Util.h"
#include "MathUtil.h"
#include <sstream>
#include <simd/simd.h>
//#include <LibraryPure.lib.h>
//#include <LibraryPureMMorder.lib.h>
//#include <Library_Training2_run17.lib.h>
//#include <LibraryExpOpt.lib.h>
//#include <LibraryExpOpt7_10_polished.lib.h>
//#include <LibraryMix.lib.h>
#include <ExpOpt3_8_polished.cs32.lib.h>
#include <libPolished_8.lib.h>
#include <Library255_may17.lib.h>
#include <cs219.lib.h>

//#include <libPure_blosum62_255.lib.h>
//#include <libPure_blosum62_32.lib.h>
//#include <libPureMMorder20_blosum62_255.lib.h>
//#include <LibraryExpOpt3_8_polished2.lib.h>

// **********************************************
// ********** ProfileStates *********************
// **********************************************
ProfileStates::ProfileStates( int pAlphSize, double * pBack)
{
    //std::string libraryString((const char *)LibraryExpOpt_lib, LibraryExpOpt_lib_len);
    //std::string libraryString((const char *) LibraryPureMMorder_lib, LibraryPureMMorder_lib_len);
    //std::string libraryString((const char *) LibraryExpOpt7_10_polished_lib, LibraryExpOpt7_10_polished_lib_len);
    //std::string libraryString((const char *) Library_Training2_run17_lib, Library_Training2_run17_lib_len);
    std::string libraryString;
    switch (pAlphSize){
        case 8:
            //libraryString=std::string((const char *)libPure_blosum62_32_lib, libPure_blosum62_32_lib_len);
            libraryString=std::string((const char *)libPolished_8_lib, libPolished_8_lib_len);
            break;
        case 32:
            //libraryString=std::string((const char *)libPure_blosum62_32_lib, libPure_blosum62_32_lib_len);
            libraryString=std::string((const char *)ExpOpt3_8_polished_cs32_lib, ExpOpt3_8_polished_cs32_lib_len);
            break;
        case 219:
            //libraryString=std::string((const char *)libPure_blosum62_255_lib, libPure_blosum62_255_lib_len);
            //libraryString=std::string((const char *)libPureMMorder20_blosum62_255_lib, libPureMMorder20_blosum62_255_lib_len);
            libraryString=std::string((const char *)cs219_lib, cs219_lib_len);
            break;
        case 255:
            //libraryString=std::string((const char *)libPure_blosum62_255_lib, libPure_blosum62_255_lib_len);
            //libraryString=std::string((const char *)libPureMMorder20_blosum62_255_lib, libPureMMorder20_blosum62_255_lib_len);
            libraryString=std::string((const char *)Library255_may17_lib, Library255_may17_lib_len);
            break;
        default:
            Debug(Debug::ERROR) << "Could not load library for alphabet size " << alphSize << "\n";
            EXIT(EXIT_FAILURE);
    }
    //std::string libraryString((const char *)LibraryExpOpt3_8_polished2_lib, LibraryExpOpt3_8_polished2_lib_len);
    //std::string libraryString((const char *)LibraryMix_lib, LibraryMix_lib_len);

    background = new float[20];
    for (size_t k = 0; k< 20 ; k++){
        background[k] = (float) pBack[k];
    }
    alphSize = 0;
    read(libraryString);
}

ProfileStates::~ProfileStates()
{
    for (size_t k = 0; k< alphSize ; k++){
        free(normalizedProfiles[k]);
    }
    delete [] normalizedProfiles;
    for (size_t k = 0; k< alphSize; k++){
        free(profiles[k]);
        free(discProfScores[k]);
    }
    delete [] discProfScores;
    delete [] profiles;
    delete [] background;
    free(prior);
}

int ProfileStates::readProfile(std::stringstream &in, float * profile,  float * normalizedProfile, float &prior) {
    // Parse and check header information
    if (!reader.StreamStartsWith(in, "ContextProfile"))
    {
        Debug(Debug::WARNING) << "Stream does not start with class id 'ContextProfile'!\n";
        return -1;
    }

    std::string line;
    line = reader.getline(in);
    if (strstr(line.c_str(), "NAME")) {
        std::string name = reader.ReadString(line.c_str(), "NAME", "Unable to parse context profile 'NAME'!");
        names.push_back(name);
        line = reader.getline(in);
    }
    else{
        names.push_back("0"); // default name
    }

    prior = reader.ReadDouble(line.c_str(), "PRIOR", "Unable to parse context profile 'PRIOR'!");
    line = reader.getline(in);
    if (strstr(line.c_str(), "COLOR")) {
        std::string coldef;
        coldef = reader.ReadString(line.c_str(), "COLOR", "Unable to parse context profile 'COLOR'!");
        colors.push_back(Color(coldef));
        line = reader.getline(in);
    } else {
        colors.push_back(Color(0.0,0.0,0.0));
    }
    reader.ReadBool(line.c_str(), "ISLOG", "Unable to parse context profile 'ISLOG'!");
    line = reader.getline(in);
    reader.ReadInt(line.c_str(), "LENG", "Unable to parse context profile 'LENG'!");
    line = reader.getline(in);
    int nalph = reader.ReadInt(line.c_str(), "ALPH", "Unable to parse context profile 'ALPH'!");

    //if (is_log) prior = log(prior);
    //assert(len == 1); // no context

    if (nalph != 20)
    {
        Debug(Debug::WARNING) << "Alphabet size of serialized context profile should be " << 20 << " but is actually "<< nalph <<"!\n";
        return -1;
    }

    line = reader.getline(in);
    if (strstr(line.c_str(), "PROBS"))
    {
        line = reader.getline(in);
        //float *profile = (float *)mem_align(16, 20);
        char *pos = (char *) line.c_str();

        // Skip the first fied that is equal to 1
        pos += Util::skipWhitespace(pos);
        pos += Util::skipNoneWhitespace(pos);
        pos += Util::skipWhitespace(pos);
        float s = 0.0;

        // Store the probabilities
        for (int k = 0 ; k < nalph ; k++)
        {
            float score = std::strtod(pos, NULL);
            float prob = MathUtil::fpow2(-score/kScale);
            profile[ProfileStates::hh2mmseqsAAorder(k)] = prob;// /background[k];
            s+=prob;
            char* oldPos = pos;
            pos += Util::skipNoneWhitespace(pos);
            pos += Util::skipWhitespace(pos);

            if (pos == oldPos) // we reached the end of the line
            {
                Debug(Debug::WARNING) << "Context profile should have " << nalph << " columns but actually has " << k+1 <<"!\n";
                return -1;
            }
        }

        float norm = sqrt(ScalarProd20(profile,profile));
        for (int k = 0 ; k < nalph ; k++)
        {
            normalizedProfile[k] = profile[k] / norm;
        }
    } else {
        Debug(Debug::WARNING) << "Cannot find the probabilities of the state in the file !\n";
        return -1;
    }

    if (!reader.StreamStartsWith(in, "//"))
    {
        Debug(Debug::WARNING) << "Expected end of profile description!\n";
        return -1;
    }
    return 0;
}

int ProfileStates::read(std::string libraryData) {

    std::stringstream in(libraryData);
    // Parse and check header information
    if (!reader.StreamStartsWith(in, "ContextLibrary"))
    {
        Debug(Debug::WARNING) << "LibraryData does not start with ContextLibrary" << "!\n";
        return -1;
    }

    std::string line;

    if ((line = reader.getline(in)) != "")
        alphSize = reader.ReadInt(line.c_str(), "SIZE", "Unable to parse context library 'SIZE'!");
    if ((line = reader.getline(in)) != "")
        reader.ReadInt(line.c_str(), "LENG", "Unable to parse context library 'LENG'!");



    profiles = new float*[alphSize];
    normalizedProfiles = new float*[alphSize];
    prior = (float*) mem_align(ALIGN_FLOAT, alphSize * sizeof(float));

    // Read profiles
    size_t k;
    float zPrior = 0.0;
    for (k = 0; k < alphSize && in.good(); ++k)
    {
        profiles[k]           = (float *)mem_align(ALIGN_FLOAT, 20 * sizeof(float));
        normalizedProfiles[k] = (float *)mem_align(ALIGN_FLOAT, 20 * sizeof(float));
        readProfile(in, profiles[k], normalizedProfiles[k],prior[k]);
        zPrior += prior[k];
	
    }

    if (zPrior == 0.0)// In case prior is not def in the lib, approximate it by
    {		      // projection on the bck proba
	for (k = 0; k < alphSize && in.good(); ++k)
    	{
       		for (size_t a = 0 ; a < 20 ; a++)
           		prior[k] += profiles[k][a] * background[a];
        	zPrior += prior[k];
	}
    }

    if (k != alphSize)
    {
        Debug(Debug::WARNING) << "Serialized context library should have "<< alphSize << " profiles but actually has " << (unsigned int)k << "\n";
        return -1;
    }


    for (k = 0; k < alphSize; ++k)
    {
        prior[k] /= zPrior;
        //DEBUG: std::cout<<"Prior["<<k<<"] = "<<prior[k]<<std::endl;
    }
    
 /*   zPrior = 0;
    for (k = 0; k < 20; ++k)
    {
        zPrior += prior[k];
    }
    

    
    for (k = 20; k < alphSize; ++k)
    {
        prior[k] = 0;
    }*/


    discProfScores = new float*[alphSize];
    for (k = 0; k< alphSize ; k++)
    {
        unsigned int ceilAlphSize = MathUtil::ceilIntDivision(static_cast<unsigned int>(alphSize), static_cast<unsigned int >(VECSIZE_FLOAT));
        discProfScores[k] = (float*) mem_align(ALIGN_FLOAT, sizeof(float)*VECSIZE_FLOAT*ceilAlphSize);
        memset(discProfScores[k], 0,ceilAlphSize*VECSIZE_FLOAT* sizeof(float) );
        
        
        for (size_t l = 0; l< alphSize ; l++)
        {
            discProfScores[k][l] = score(k,l);
            // DEBUG: print the disc sub mat
            //std::cout<<discProfScores[k][l]<<"\t";
        }
        //std::cout<<"\n";
    }


    return 0;
}

float ProfileStates::entropy(float *distribution)
{
    float entropy = 0.0;
    for (size_t a = 0 ; a < 20 ; a++)
    {
        entropy -= distribution[a] * MathUtil::flog2(distribution[a]);
    }
    return entropy;
}

float ProfileStates::getScoreNormalization() {

    

    float *maxScore;
    maxScore = new float[alphSize];

    for (size_t k = 0; k<alphSize; k++)
    {
        maxScore[k] = FLT_MIN;
        for (size_t a = 0; a < Sequence::PROFILE_AA_SIZE; a++)
        {
            maxScore[k] = std::max(maxScore[k], profiles[k][a] / background[a]);
        }
    }
    float exp = 0.0;
    for (size_t k = 0; k<alphSize; k++) {
        exp += MathUtil::flog2(maxScore[k]) * prior[k];
    }
    exp /= entropy(background);

    exp = 1.0 + (exp - 1.0)/2;

    Debug(Debug::INFO) << "Score normalization : " << 1.0/exp << "\n";

    exp = 1.0;
    return exp;
}

void ProfileStates::discretize(const float* sequence, size_t length, std::string &result)
{
    
    char closestState = 0;
    float curDiffScore;
    float* profileCol;
    float* repScore = (float*)mem_align(ALIGN_FLOAT, 256*sizeof(float));
    memset(repScore, 0, sizeof(float)*256);
    for (size_t i = 0 ; i<length ; i++)
    {
        float minDiffScore = FLT_MAX;
        profileCol = (float*)
                &sequence[i * Sequence::PROFILE_AA_SIZE];

        float maxScore = -FLT_MIN;
//        size_t maxScoreIndex = 0;
	
	/*for (size_t a = 0 ; a < 20 ; a++)
        {
                std::cout<<profileCol[a]<<"\t";
        }*/

//        for(size_t pos = 0; pos < 20; pos++){
//            printf("%.3f\t", profileCol[pos]);
//        }
//        std::cout << std::endl;
        // S(profile, c_k)
        for (size_t k=0;k<alphSize;k++)
        {
            repScore[k] = score(profileCol, profiles[k]);

//            for(size_t pos = 0; pos < 20; pos++){
//                printf("%.3f\t", profiles[k][pos]);
//            }
//            std::cout << repScore[k] << std::endl;
            if (repScore[k]>maxScore)
            {
                maxScore = repScore[k];
//                maxScoreIndex = k;
            }
        }

        // Find the k that minimizes sum_l prior_l*(S(profile, c_l) - S(c_k,c_l))^2
        for (size_t k=0;k<alphSize;k++)
        {
            curDiffScore = 0.0;
            simd_float curDiffScoreSimd=simdf32_setzero(0);
            unsigned int ceilAlphSize = MathUtil::ceilIntDivision(static_cast<unsigned int >(alphSize), static_cast<unsigned int >(VECSIZE_FLOAT));
//            for (size_t l=0;l<alphSize;l++)
//            {
//                float diff = repScore[l] - discProfScores[k][l];
//                tmp += diff*diff;
//            }


            for (size_t l=0;l<ceilAlphSize*VECSIZE_FLOAT;l+=VECSIZE_FLOAT)
            {
                simd_float priorLSimd = simdf32_load(&prior[l]);
                simd_float repScoreLSimd = simdf32_load(&repScore[l]);
                simd_float discProfScoresLSimd = simdf32_load(&discProfScores[k][l]);
                simd_float diff = simdf32_sub(repScoreLSimd, discProfScoresLSimd);
                simd_float diffSquared = simdf32_mul(diff, diff);
                simd_float postDiff = simdf32_mul(priorLSimd, diffSquared);
                curDiffScoreSimd = simdf32_add(curDiffScoreSimd, postDiff);
            }
            float * curDiffScoreSimdFlt = (float*) &curDiffScoreSimd;
            //std::cout<<k<<":"<<*curDiffScoreSimdFlt<<"\n";
            /*if (alphSize==255)
		{

		    for (size_t y = 0; y<21;y++)
			{
			std::cout<<k<<" "<<prior[y]<<" "<<repScore[y]<<" "<<discProfScores[k][y]<<" "<<std::endl;
			}
		}*/
            for (size_t l=0;l<VECSIZE_FLOAT;l++){
                curDiffScore+= curDiffScoreSimdFlt[l];
            }

            if (curDiffScore < minDiffScore)
            {
                minDiffScore = curDiffScore;
                closestState = k;
            }
        }
	//std::cout<<"Pos "<<i<<", closest state: "<<(int)closestState<<", max score index: "<<maxScoreIndex<<"\n";
        result.push_back(closestState);//(maxScoreIndex);//
    }
    free(repScore);
}



void ProfileStates::discretizeCs219(const float* sequence, size_t length, std::string &result)
{
//    float curDiffScore;
    float* profileCol;
    float* repScore = (float*)mem_align(ALIGN_FLOAT, 256*sizeof(float));
    memset(repScore, 0, sizeof(float)*256);
    for (size_t i = 0 ; i<length ; i++)
    {
//        float minDiffScore = FLT_MAX;
        profileCol = (float*)
                &sequence[i * Sequence::PROFILE_AA_SIZE];
        // Calculate posterior probabilities given sequence window around 'i'
        double max = -FLT_MAX;
        size_t k_max = 0;
        for (size_t k = 0; k < alphSize; ++k) {
            repScore[k] = prior[k] * score(profiles[k], profileCol);
            k_max = (repScore[k] > max ) ? k : k_max;
            max   = (repScore[k] > max ) ? repScore[k] : max;
        }
        result.push_back(k_max);
    }
    free(repScore);
}



float ProfileStates::score(float* profileCol, size_t state)
{
    float *stateProfile = profiles[state];
    return score(profileCol,stateProfile);
}

float ProfileStates::score(size_t stateA, size_t stateB)
{
    float *stateProfileA = profiles[stateA];
    float *stateProfileB = profiles[stateB];
    return score(stateProfileA,stateProfileB);
}

float ProfileStates::score(float* profileA, float* profileB)
{
    return score(profileA,background,profileB);
}

float ProfileStates::distance(float* profileA, float* profileB)
{
    return score(profileA,profileA) + score(profileB,profileB) -2*score(profileA,profileB);
}


