#ifndef MMSEQS_PROFILESTATES
#define MMSEQS_PROFILESTATES

#include <cstdlib>
#include "Util.h"
#include "Debug.h"
#include "SubstitutionMatrix.h"

#define kScale 1000 // Scaling factor for the profile scores in the library file
#define KB 1024
#define HALF_LOCAL_WINDOW 3



/*
 * Example of usage
 * 
 *     
 * 
    profileStates ps("../data/Library1.lib","/home/clovis/Software/MMseqs2/data/blosum62.out");
    size_t L = 10;
    float profile[20*L];
    float avgProfile[20*L];
    size_t disc[L];
    // build a profile sequence
    for (size_t i =0 ; i<20*L/2;i+=20)
    {
        for(size_t aa = 0 ; aa < 20 ; aa++)
        {
            profile[i+aa] = 0.0;
            avgProfile[i+aa] = 0.0;
            if (aa == 2)
                profile[i+aa] = 1.0;
        }
    }
    for (size_t i =20*L/2 ; i<20*L;i+=20)
    {
        for(size_t aa = 0 ; aa < 20 ; aa++)
        {
            avgProfile[i+aa] = 0.0;
            profile[i+aa] = 0.0;
            if (aa == 4)
                profile[i+aa] = 1.0;
        }
    }
    
    ps.slidingAverage(profile,L,avgProfile);
    ps.discretize(profile,avgProfile,L,disc);
    for (size_t i = 0 ; i< L ; i++)
        std::cout<<i<<","<<disc[i]<<std::endl;
        
 * 
 * 
 * */

class libraryReader {
public:
    bool StreamStartsWith(FILE* fp, const char* id);
    int ReadInt(const char* line,const char* label,const char* errmsg);
    double ReadDouble(const char* line,const char* label,const char* errmsg);
    std::string ReadString(const char* line,const char* label,const char* errmsg);
    bool ReadBool(const char* line,const char* label,const char* errmsg);
    const char* strscn(const char* str) ;
    int chomp(char* str);
    static std::vector<std::string> tokenize(const char* str, char sep);
                                  
    char* fgetline(char* str, int maxlen, FILE* file);
};

struct Color {
    Color(double r = 1.0, double g = 1.0, double b = 1.0)
	    : red(r), green(g), blue(b) {}

    Color(std::string coldef) {
        std::vector<std::string> tokens;
        tokens = libraryReader::tokenize(coldef.c_str(), ',');
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


class profileStates {
public:
    
    profileStates() {};
    profileStates(std::string filename, std::string substitutionMatrixFile);
    ~profileStates();
    
    int read(std::string filename);
    int readProfile(FILE* fin);
    
   
    size_t discretize(float * profile);
    void discretize(float * sequence, float* avgSequence, size_t length, size_t *result);
    float *getProfile(size_t state);
     
    float score(float* profile, size_t state);
    float score(size_t stateA, size_t stateB);
    // Score a sequence with local AA bias correction
    float score(float* profile, float* avgProf, size_t state);
    float score(float* profileA, float* avgProfA, float* profileB);
    
    
    void slidingAverage(float* sequence, size_t length, float* avgSequence);


  
    size_t getAlphSize() {return alphSize;};
    
private:
    libraryReader reader;
    std::vector<float*> profiles;
    std::vector<float*> normalizedProfiles;
    std::vector<Color> colors;
    std::vector<std::string> names;
    
    SubstitutionMatrix* subMat;
    float* background;
    
    
    float distance(float* profileA, float* profileB);
    float score(float* profileA, float* profileB);
    float innerProd(float* profileA, float* profileB);
    
    size_t alphSize;
};

#endif