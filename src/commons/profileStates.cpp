
#include "profileStates.h"
#include "Util.h"
#include "MathUtil.h"
#include <sstream>


// Returns pointer to first non-white-space character in str OR to NULL if none
// found
const char* libraryReader::strscn(const char* str) {
  if (!str) return NULL;
  const char* ptr = str;
  while (*ptr != '\0' && isspace(*ptr)) ++ptr;
  return (*ptr == '\0') ? NULL : ptr;
}

// Returns true iff next non-blank line in file 'fp' contains string 'id'.
bool libraryReader::StreamStartsWith(FILE* fp, const char* id) {
  char buffer[KB];
  while (fgetline(buffer, KB, fp)) if (strscn(buffer)) break;
  return strstr(buffer, id) == buffer;
}



// Parse serialization record and return integer value following label 'str' in
// line read from file pointer 'fp'.
int libraryReader::ReadInt(const char* line,
                   const char* label,
                   const char* errmsg = NULL) {
  int rv = 0;
  if (strstr(line, label)) {
    const char* ptr = line + strlen(label);
    rv = atoi(ptr);
  } else if (errmsg) {
        Debug(Debug::WARNING) << errmsg;
  }
  return rv;
}

// Parse serialization record and return double value following label 'label' in
// line read from file pointer 'fp'.
double libraryReader::ReadDouble(const char* line,
                         const char* label,
                         const char* errmsg = NULL) {
  double rv = 0;
  if (strstr(line, label)) {
    rv = atof(line + strlen(label));
  } else if (errmsg) {
    Debug(Debug::WARNING) << errmsg;
  }
  return rv;
}


// Parse serialization record and return string following label 'label' in
// line read from file pointer 'fp'.
std::string  libraryReader::ReadString(const char* line,
                              const char* label,
                              const char* errmsg = NULL) {
  std::string rv;
  if (strstr(line, label)) {
    const char* ptr = strscn(line + strlen(label));
    rv = ptr;
  } else if (errmsg) {
    Debug(Debug::WARNING) << errmsg;
  }
  return rv;
}


// Removes the newline and other control characters at the end of a string (if
//  present) and returns the new length of the string (-1 if str is NULL)
int libraryReader::chomp(char* str) {
  if (!str) return -1;
  int l = 0;
  for (l = strlen(str) - 1; l >= 0 && str[l] < 32; --l)
    /* do nothing */;
  str[++l] = '\0';
  return l;
}

char* libraryReader::fgetline(char* str, int maxlen, FILE* file) {
  if (!fgets(str, maxlen, file)) return NULL;
  if (chomp(str) + 1 >= maxlen)  // if line is cut after maxlen characters...
    while (fgetc(file) != '\n')  // ... read in rest of line
      /* do nothing */;
  return(str);
}



// Parse serialization record and return bool value following label 'str' in
// line read from file pointer 'fp'.
bool libraryReader::ReadBool(const char* line,
                     const char* label,
                     const char* errmsg = NULL) {
  bool rv = false;
  if (strstr(line, label)) {
    const char* ptr = line + strlen(label);
    if (strchr(ptr, 'T') != NULL || strchr(ptr, '1') != NULL)
      rv = true;
    else if (strchr(ptr, 'F') != NULL || strchr(ptr, '0') != NULL)
      rv = false;
    else if (errmsg)
        Debug(Debug::WARNING) << errmsg;
  } else if (errmsg) {
    Debug(Debug::WARNING) << errmsg;
  }
  return rv;
}


std::vector<std::string> libraryReader::tokenize(const char* str, char sep = ' ') {
    
    std::stringstream tok;
    std::vector<std::string> tokens;
    size_t pos=0;
    
    while (str[pos] != '\0')
    {
        if (str[pos] == sep)
        {
            tokens.push_back(tok.str());
            tok.str("");
        }
        
        // go to the next non-empty field
        while (str[pos] != '\0' && str[pos] == sep)
            pos++;
            
        tok << str[pos];
        pos++;
    }
    tokens.push_back(tok.str()); 
    return tokens;
}


// **********************************************
// ********** profileStates *********************
// **********************************************

profileStates::profileStates(std::string filename, std::string substitutionMatrixFile)
{
    subMat = new SubstitutionMatrix(substitutionMatrixFile.c_str(), 8.0, -0.2);
    background = new float[20];
    
    for (size_t k = 0; k< 20 ; k++)
        background[k] = sqrt((float) subMat->pBack[k]);
        
    alphSize = 0;
    read(filename);
}

int profileStates::readProfile(FILE* fin) {
    // Parse and check header information
    if (!reader.StreamStartsWith(fin, "ContextProfile"))
    {
        Debug(Debug::WARNING) << "Stream does not start with class id 'ContextProfile'!\n";
        return -1;
    }
        
    
    char buffer[KB];
    reader.fgetline(buffer, KB, fin);
    if (strstr(buffer, "NAME")) {
        std::string name = reader.ReadString(buffer, "NAME", "Unable to parse context profile 'NAME'!");
        names.push_back(name);
        reader.fgetline(buffer, KB, fin);
    }
    else{
        names.push_back("0"); // default name
    }
    
    
    reader.ReadDouble(buffer, "PRIOR", "Unable to parse context profile 'PRIOR'!");
    reader.fgetline(buffer, KB, fin);
    if (strstr(buffer, "COLOR")) {
        std::string coldef;
        coldef = reader.ReadString(buffer, "COLOR", "Unable to parse context profile 'COLOR'!");
        colors.push_back(Color(coldef));
        reader.fgetline(buffer, KB, fin);
    } else {
        colors.push_back(Color(0.0,0.0,0.0));
    }
    



    reader.ReadBool(buffer, "ISLOG", "Unable to parse context profile 'ISLOG'!");
    reader.fgetline(buffer, KB, fin);
    reader.ReadInt(buffer, "LENG", "Unable to parse context profile 'LENG'!");
    reader.fgetline(buffer, KB, fin);
    int nalph = reader.ReadInt(buffer, "ALPH", "Unable to parse context profile 'ALPH'!");
    
    //if (is_log) prior = log(prior);
    //assert(len == 1); // no context
    
    if (nalph != 20)
    {
        Debug(Debug::WARNING) << "Alphabet size of serialized context profile should be " << 20 << " but is acutally "<< nalph <<"!\n";
        return -1;      
    }

    reader.fgetline(buffer, KB, fin);
    if (strstr(buffer, "PROBS"))
    {
        reader.fgetline(buffer, KB, fin);
        float *profile;
        profile = new float[20];
        char *pos = buffer;
        
        // Skip the first fied that is equal to 1
        pos += Util::skipWhitespace(pos);
        pos += Util::skipNoneWhitespace(pos);
        pos += Util::skipWhitespace(pos);
        float s = 0.0;
        
        // Store the probabilities
        for (size_t k = 0 ; k < nalph ; k++)
        {
            float score = strtof(pos, NULL);
            
            profile[k] = MathUtil::fpow2(-score/kScale);// /background[k];
            s+=profile[k];
            char* oldPos = pos;
            pos += Util::skipNoneWhitespace(pos);
            pos += Util::skipWhitespace(pos);
            
            if (pos == oldPos) // we reached the end of the line
            {
                Debug(Debug::WARNING) << "Context profile should have " << nalph << " columns but actually has " << k+1 <<"!\n";
                return -1;      
            }
        }
        profiles.push_back(profile);
        
        float norm = sqrt(innerProd(profile,profile));
        float *normalizedProfile;
        normalizedProfile = new float[20]; // q/||q||
        for (size_t k = 0 ; k < nalph ; k++)
        {
            normalizedProfile[k] = profile[k]/norm;
        }
        normalizedProfiles.push_back(normalizedProfile);
        
    } else {
        Debug(Debug::WARNING) << "Cannot find the probabilities of the state in the file !\n";
        return -1;   
    }

    if (!reader.StreamStartsWith(fin, "//"))
    {
        Debug(Debug::WARNING) << "Expected end of profile description!\n";
        return -1;
    }
    
    return 0;


}



int profileStates::read(std::string stateFileName) {
    FILE* fin;
    fin = fopen(stateFileName.c_str(),"r");
    
    
    
    // Parse and check header information
    if (!reader.StreamStartsWith(fin, "ContextLibrary"))
    {
        Debug(Debug::WARNING) << "File " << stateFileName << " does not start with ContextLibrary" << "!\n";
        fclose(fin);
        return -1;
    }

    char buffer[KB];
    
    if (reader.fgetline(buffer, KB, fin))
        alphSize = reader.ReadInt(buffer, "SIZE", "Unable to parse context library 'SIZE'!");
    if (reader.fgetline(buffer, KB, fin))
        reader.ReadInt(buffer, "LENG", "Unable to parse context library 'LENG'!");

    
    // Read profiles
    size_t k;
    for (k = 0; k < alphSize && !feof(fin); ++k)
    {
        readProfile(fin);
    }
    
    if (k != alphSize)
    {
        Debug(Debug::WARNING) << "Serialized context library should have "<< alphSize << " profiles but actually has " << (unsigned int)k << "\n";
        fclose(fin);
        return -1;
    }
    
    fclose(fin);
    return 0;
}

float profileStates::score(float* profileA, float* profileB)
{
    return MathUtil::flog2(innerProd(profileA,profileB));
        
}

float profileStates::innerProd(float* profileA, float* profileB)
{
    float s = 0.0;
    
    for (size_t k = 0 ; k<20;k++)
        s += profileA[k]*profileB[k]; // profiles already normalized by the background
    return s;
        
}




size_t profileStates::discretize(float* profile)
{
    // There is no need to compute everything in the distance
    // we only seek the state s that maximize p.s/||s|| (where everything
    // is background-normalized)
    
    float maxScore = 0.0;
    size_t closestState = 0;
    float curScore;
    /*
    float profileBck[20];
    for (size_t i=0;i<20;i++)
    {
        profileBck[i] = profile[i] / background[i];
    }*/
    
    for (size_t k=0;k<alphSize;k++)
    {
        curScore = score(profile,normalizedProfiles[k]);
        if (curScore > maxScore)
        {
            maxScore = curScore;
            closestState = k;
        }
    }
    return closestState;
}


void profileStates::discretize(float* sequence, float* avgSequence, size_t length, size_t* result)
{
    float maxScore;
    size_t closestState;
    float curScore;
    float* profileCol;
    float* avgProfileCol;

    for (size_t i = 0 ; i<length ; i++)
    {
        profileCol = sequence + i*20;
        avgProfileCol = avgSequence + i*20;
        
        for (size_t k=0;k<alphSize;k++)
        {
            curScore = score(profileCol,avgProfileCol,normalizedProfiles[k]);
            
            if (!k||curScore > maxScore)
            {
                maxScore = curScore;
                closestState = k;
            }
        }
        result[i] = closestState;
    }
}


void profileStates::slidingAverage(float* sequence, size_t length, float* avgSequence)
{
    float windowSize = 2*HALF_LOCAL_WINDOW + 1;
    
    if (!length)
    {
        return;
    }
    
    for (size_t p = 0 ; p<20 ; p++)
    {
        avgSequence[p] = background[p] * HALF_LOCAL_WINDOW;
    }
    for (size_t p = 0 ; p<(HALF_LOCAL_WINDOW + 1)*20 ; p++)
    {
        avgSequence[p%20] += sequence[p];
    }
    
    for (size_t aa = 0 ; aa<20 ; aa++)
    {
        avgSequence[aa] /= windowSize;
    }
    
    for (size_t i = 1 ; i<length ; i++)
    {
        if (i>HALF_LOCAL_WINDOW)
            for (size_t aa = 0 ; aa<20 ; aa++)
            {
                avgSequence[i*20 + aa] = avgSequence[(i-1)*20 + aa] - sequence[(i - HALF_LOCAL_WINDOW - 1)*20 + aa] / windowSize;
                if (avgSequence[i*20 + aa]<0)
                    avgSequence[i*20 + aa] = 0.0;
            }
        else {
            for (size_t aa = 0 ; aa<20 ; aa++)
            {
                avgSequence[i*20 + aa] = avgSequence[(i-1)*20 + aa] - background[aa] / windowSize;
                if (avgSequence[i*20 + aa]<0)
                    avgSequence[i*20 + aa] = 0.0;
            }
            
        }
            
        if (i<length - HALF_LOCAL_WINDOW)
        {
            for (size_t aa = 0 ; aa<20 ; aa++)
            {
                avgSequence[i*20 + aa] += sequence[(i + HALF_LOCAL_WINDOW) *20 + aa ] / windowSize;
            }
        }
        else {
            for (size_t aa = 0 ; aa<20 ; aa++)
            {
                avgSequence[i*20 + aa] += background[aa] / windowSize;
            }
        }
        
    }
}

float profileStates::score(float* profileCol, float* avgProfCol, size_t state)
{
    return score(profileCol,avgProfCol,profiles[state]);
}



float profileStates::score(float* profileColA, float* avgProfColA, float* profileColB)
{
    float result = 0.0;
    for (size_t aa = 0 ; aa < 20 ; aa++)
    {
        result += profileColA[aa] * profileColB[aa] / avgProfColA[aa];
    }
    return MathUtil::flog2(result);
}


float profileStates::score(float* profileCol, size_t state)
{
    float result = 0.0;
    float *stateProfile = profiles[state];
    
    for (size_t aa = 0 ; aa < 20 ; aa++)
        result += profileCol[aa] * stateProfile[aa]; // profileCol already background corrected
        
    return MathUtil::flog2(result);
}


float profileStates::score(size_t stateA, size_t stateB)
{
    float result = 0.0;
    float *stateProfileA = profiles[stateA];
    float *stateProfileB = profiles[stateB];
    
    for (size_t aa = 0 ; aa < 20 ; aa++)
        result += stateProfileA[aa] * stateProfileB[aa] / background[aa]; 
        
    return MathUtil::flog2(result);
}

float profileStates::distance(float* profileA, float* profileB)
{
    float ip = innerProd(profileA,profileB);
    return MathUtil::flog2(innerProd(profileA,profileA) * innerProd(profileB,profileB) / (ip*ip));
}

float* profileStates::getProfile(size_t state)
{
    return profiles[state];
}


profileStates::~profileStates()
{
    for (size_t k = 0; k<profiles.size() ; k++)
        delete profiles[k];
    delete background;
    delete subMat;
}

