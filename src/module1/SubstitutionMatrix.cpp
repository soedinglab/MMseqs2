#include "SubstitutionMatrix.h"

SubstitutionMatrix::SubstitutionMatrix(char* scoringMatrixFileName_, size_t alphabetSize_):
scoringMatrixFileName(scoringMatrixFileName_),
ALPHABET_SIZE(alphabetSize_)
{
    // init [amino acid <-> int] mappings
    int2aa = new char[ALPHABET_SIZE];
    int2aa[0] = 'A';
    int2aa[1] = 'R';
    int2aa[2] = 'N';
    int2aa[3] = 'D';
    int2aa[4] = 'C';
    int2aa[5] = 'Q';
    int2aa[6] = 'E';
    int2aa[7] = 'G';
    int2aa[8] = 'H';
    int2aa[9] = 'I';
    int2aa[10] = 'L';
    int2aa[11] = 'K';
    int2aa[12] = 'M';
    int2aa[13] = 'F';
    int2aa[14] = 'P';
    int2aa[15] = 'S';
    int2aa[16] = 'T';
    int2aa[17] = 'W';
    int2aa[18] = 'Y';
    int2aa[19] = 'V';
    
    aa2int = new int['Z'+1];
    for (int i = 0; i <= 'Z'; ++i) aa2int[i]=-1;
    for (int i = 0; i < ALPHABET_SIZE; ++i){
        aa2int[int2aa[i]] = i;
    }
    
    
    pBackground = new double[ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; i++)
        pBackground[i] = 0.0;
    
    // read amino acid substitution matrix from file
    std::string fileName(scoringMatrixFileName);
    if (fileName.substr(fileName.length()-4, 4).compare(".mat") == 0)
        readScoringMatrix();
    else
        readBiasedScoringMatrix(2.0, 0.0);
    
}

SubstitutionMatrix::~SubstitutionMatrix(){
    
    delete[] aa2int;
    delete[] int2aa;
    for (int i = 0; i < ALPHABET_SIZE; ++i)
        delete[] scMatrix[i];
    delete[] scMatrix;
    delete[] pBackground;
}

std::vector<std::string>& SubstitutionMatrix::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    
    std::stringstream ss(s);
    std::string item;
    
    while(std::getline(ss, item, delim)) {
        if (item.compare("") != 0)
            elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> SubstitutionMatrix::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

void SubstitutionMatrix::readScoringMatrix(){
    
    scMatrix = new short* [ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; i++){
        scMatrix[i] = new short[ALPHABET_SIZE];
    }
    
    std::ifstream in(scoringMatrixFileName);
    if(in.fail()){
        std::cerr << "Cannot read " << scoringMatrixFileName << "\n";
        exit(1);
    }
    std::string line;
    getline( in, line );
    while (line.at(0) == '#')
        getline (in, line);
    trim(line);
    std::vector<std::string> aas = split(line, ' ');
    
    char* aa_entries = new char [aas.size()];
    int pos = 0;
    for (int i = 0; i < aas.size(); i++){
        aa_entries[pos] = aas[i].at(0);
        pos++;
    }
    
    int row = 0;
    while(in.good()){
        getline (in, line);
        if (!isalpha(line[0]) || line[0] == 'B' || line[0] == 'Z' || line[0] == 'J' || line[0] == 'X')
            continue;
        
        std::vector<std::string> scores = split(line, ' ');
        for (int i = 1; i < ALPHABET_SIZE+1; i++){
            scMatrix[aa2int[aa_entries[row]]][aa2int[aa_entries[i-1]]] = (short)atoi(scores[i].c_str());
        }
        row++;
    }
    
    in.close();
    
    delete[] aa_entries;
}

void SubstitutionMatrix::readBiasedScoringMatrix(double bitFactor, double scoringBias){
    
    scMatrix = new short* [ALPHABET_SIZE];
    probMatrix = new double * [ALPHABET_SIZE];
    for (int i = 0; i < ALPHABET_SIZE; i++){
        scMatrix[i] = new short[ALPHABET_SIZE];
        probMatrix[i] = new double [ALPHABET_SIZE];
    }
    
    for( size_t i=0; i<ALPHABET_SIZE; ++i ) pBackground[i] = 0.0f;
    
    std::ifstream in(scoringMatrixFileName);
    if( in.fail() ) {
        std::cerr << "Cannot read " << scoringMatrixFileName << "\n";
        exit(1);
    }
    int row    = 0;
    int column = 0;
    std::string line;
    bool capture = false;
    while( in.good() ){
        getline( in, line );
        if( line.length()>11 && line.substr(0, 11)!="Frequencies" && !capture )
            continue;
        if( line.length()>11 && line.substr(0, 11)=="Frequencies"){
            capture=true;
            continue;
        }
        if( row==ALPHABET_SIZE ) break;
        std::stringstream stream(line);
        std::string h;
        stream >> h;
        if( h=="" ) continue;
        if (!isalpha(h.at(0))){
            column = 0;
            stream.clear();
            stream.str(line);
            float f;
            while( stream >> f ){
                probMatrix[row][column] = f;
                probMatrix[column][row] = f;
                ++column;
            }
            ++row;
        }
    }
    in.close();
    
    double sum=0.0;
    for(size_t i=0; i<ALPHABET_SIZE; ++i)
        for(size_t j=0; j<ALPHABET_SIZE; ++j){
            if( i==j ) pBackground[i] += probMatrix[i][j];
            else       pBackground[i] += (probMatrix[i][j]/2.0f);
            if( j<=i ) sum += probMatrix[i][j];
        }
    
    const double _2sum = 2.0*sum;
    double pbsum = 0.0;
    for(size_t i=0; i<ALPHABET_SIZE; ++i){
        pbsum += pBackground[i];
        for(size_t j=0; j<ALPHABET_SIZE; ++j)
            if( i==j ) probMatrix[i][j] = probMatrix[i][j] / sum;
            else       probMatrix[i][j] = probMatrix[i][j] / _2sum;
    }
    
    for(size_t i=0; i<ALPHABET_SIZE; ++i)pBackground[i] /= sum;
    
    for(size_t i=0; i<ALPHABET_SIZE; ++i){
        for(size_t j=0; j<ALPHABET_SIZE; ++j){
            scMatrix[i][j] = (short)floor (bitFactor * _log2(probMatrix[i][j]/(pBackground[i]*pBackground[j])) + scoringBias + 0.5);
        }
    }
}
