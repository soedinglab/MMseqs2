
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
//#include <LibraryExpOpt3_8_polished2.lib.h>

// **********************************************
// ********** ProfileStates *********************
// **********************************************
ProfileStates::ProfileStates( double * pBack)
{
    //std::string libraryString((const char *)LibraryExpOpt_lib, LibraryExpOpt_lib_len);
    //std::string libraryString((const char *) LibraryPureMMorder_lib, LibraryPureMMorder_lib_len);
    //std::string libraryString((const char *) LibraryExpOpt7_10_polished_lib, LibraryExpOpt7_10_polished_lib_len);
    //std::string libraryString((const char *) Library_Training2_run17_lib, Library_Training2_run17_lib_len);
    std::string libraryString((const char *)ExpOpt3_8_polished_cs32_lib, ExpOpt3_8_polished_cs32_lib_len);
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
    delete prior;
}

int ProfileStates::readProfile(std::stringstream &in, float * profile,  float * normalizedProfile) {
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

    reader.ReadDouble(line.c_str(), "PRIOR", "Unable to parse context profile 'PRIOR'!");
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
        Debug(Debug::WARNING) << "Alphabet size of serialized context profile should be " << 20 << " but is acutally "<< nalph <<"!\n";
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
        for (size_t k = 0 ; k < nalph ; k++)
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
        for (size_t k = 0 ; k < nalph ; k++)
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
    prior = new float[alphSize];

    // Read profiles
    size_t k;
    float zPrior = 0.0;
    for (k = 0; k < alphSize && in.good(); ++k)
    {
        profiles[k]           = (float *)mem_align(16, 20 * sizeof(float));
        normalizedProfiles[k] = (float *)mem_align(16, 20 * sizeof(float));
        readProfile(in, profiles[k], normalizedProfiles[k]);
        prior[k] = 0.0;
        for (size_t a = 0 ; a < 20 ; a++)
            prior[k] += profiles[k][a] * background[a];
        zPrior += prior[k];
    }

    if (k != alphSize)
    {
        Debug(Debug::WARNING) << "Serialized context library should have "<< alphSize << " profiles but actually has " << (unsigned int)k << "\n";
        return -1;
    }

    for (k = 0; k < alphSize; ++k)
        prior[k] /= zPrior;


    discProfScores = new float*[alphSize];
    for (k = 0; k< alphSize ; k++)
    {
        discProfScores[k] = new float[alphSize];
        for (size_t l = 0; l< alphSize ; l++)
            discProfScores[k][l] = score(k,l);
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

    //std::cout << "Score normalization : " << 1.0/exp << std::endl;

    return 1.0/exp;
}

void ProfileStates::discretize(const float* sequence, size_t length, std::string &result)
{
    float minDiffScore;
    char closestState;
    float curDiffScore;
    float* profileCol;
    float* repScore;
    repScore = new float[alphSize];

    for (size_t i = 0 ; i<length ; i++)
    {
        profileCol = (float *)sequence + i*Sequence::PROFILE_AA_SIZE;


        // S(profile, c_k)
        for (size_t k=0;k<alphSize;k++)
        {
            repScore[k] = score(profileCol,profiles[k]);
        }

        // FInd the k that minimizes sum_l (S(profile, c_l) - S(c_k,c_l))^2
        for (size_t k=0;k<alphSize;k++)
        {
            curDiffScore = 0.0;
            for (size_t l=0;l<alphSize;l++)
            {
                float diff = repScore[l] - discProfScores[k][l];
                curDiffScore += diff*diff;
            }

            if (!k||(curDiffScore < minDiffScore))
            {
                minDiffScore = curDiffScore;
                closestState = k;
            }
        }

        result.push_back(closestState);
    }
    delete [] repScore;
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


