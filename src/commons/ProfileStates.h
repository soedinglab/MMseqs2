#ifndef MMSEQS_PROFILESTATES
#define MMSEQS_PROFILESTATES

#include <cstdlib>

#include "Util.h"
#include "Debug.h"
#include "SubstitutionMatrix.h"
#include "MathUtil.h"
#include "LibraryReader.h"
#include "Sequence.h"

#define kScale 1000 // Scaling factor for the profile scores in the library file

struct Color {
    Color(double r = 1.0, double g = 1.0, double b = 1.0)
            : red(r), green(g), blue(b) {}

    Color(std::string coldef) {
        std::vector<std::string> tokens;
        tokens = LibraryReader::tokenize(coldef.c_str(), ',');
        red   = atof(tokens[0].c_str());
        green = atof(tokens[1].c_str());
        blue  = atof(tokens[2].c_str());
    }
    bool operator< (const Color& rhs) const {
        if (red != rhs.red)     return (red < rhs.red);
        if (green != rhs.green) return (green < rhs.green);
        return (blue < rhs.blue);
    }
    double red, green, blue;
};


class ProfileStates {
public:

    ProfileStates(int alphSize, double * pBack);
    ~ProfileStates();
    int read(std::string libraryData);
    int readProfile(std::stringstream &in, float * profile, float * normalizedProfile, float &prior);
    void discretize(const float* sequence, size_t length, std::string &result);
	void discretizeCs219(const float* sequence, size_t length, std::string &result);
    float getScoreNormalization();

    float* getProfile(size_t state)
    {
        return profiles[state];
    }

    float score(float* profile, size_t state);
    float score(size_t stateA, size_t stateB);

    // Score with local AA bias correction
    float score(float* profileCol, float* avgProfCol, size_t state)
    {
        return score(profileCol,avgProfCol,profiles[state]);
    }

    // Score with local AA bias correction
    float score(float* profileColA, float* avgProfColA, float* profileColB)
    {
        if (0){ // Correlation score
		float result = 0.0;
		float avgA = 0.0;
		float avgB = 0.0;
		float varA = 0.0;
		float varB = 0.0;
		float upper = 0.0;

		for (size_t aa = 0 ; aa < Sequence::PROFILE_AA_SIZE; aa++){
		    avgA += profileColA[aa];
		    avgB += profileColB[aa];

		    /*// Test of scoring scheme
		    float la = MathUtil::flog2(profileColA[aa]);
		    float lna = MathUtil::flog2(1-profileColA[aa]);
		    float lb = MathUtil::flog2(profileColB[aa]);
		    float lnb = MathUtil::flog2(1-profileColB[aa]);
		    result += (profileColB[aa] - avgProfColA[aa])*(la-lna) + (profileColA[aa] - avgProfColA[aa])*(lb-lnb);
			*/	
		    //result += profileColB[aa] * profileColA[aa] / avgProfColA[aa];
		}
		avgA /= Sequence::PROFILE_AA_SIZE;
		avgB /= Sequence::PROFILE_AA_SIZE;
		for (size_t aa = 0 ; aa < Sequence::PROFILE_AA_SIZE; aa++){
			varA += (profileColA[aa]-avgA) * (profileColA[aa]-avgA);
			varB += (profileColB[aa]-avgB) * (profileColB[aa]-avgB);
			upper += (profileColA[aa]-avgA) * (profileColB[aa]-avgB);
		}
		//result  = MathUtil::flog2(result);
		result = upper/sqrt(varA)/sqrt(varB);
		return result/2;
	} else { // HHBlits score
	        float result = 0.0;
       		for (size_t aa = 0 ; aa < Sequence::PROFILE_AA_SIZE; aa++){
			result += profileColB[aa] * profileColA[aa] / avgProfColA[aa];
		}
		result  = MathUtil::flog2(result);
		return result;
	}
    }
    float distance(float* profileA, float* profileB);

    size_t getAlphSize() {return alphSize;};


    static size_t hh2mmseqsAAorder(size_t k){
        int order[] = {0 , 14 , 11 , 2 , 1 , 13 , 3 , 5 , 6 , 7 , 9 , 8 , 10 , 4 , 12 , 15 , 16 , 18 , 19 , 17};
        return order[k];
    }

    float* prior;
    
private:
    LibraryReader reader;
    float entropy (float*);
    std::vector<Color> colors;
    std::vector<std::string> names;
    float* background;
    float score(float* profileA, float* profileB);
    size_t alphSize;
    float ** profiles;
    float ** normalizedProfiles;
    float ** discProfScores;
};

#endif
