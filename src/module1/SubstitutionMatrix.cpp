#include "SubstitutionMatrix.h"

SubstitutionMatrix::SubstitutionMatrix(char* scoringMatrixFileName)
{
    this->scoringMatrixFileName = scoringMatrixFileName;
    // read amino acid substitution matrix from file
    std::string fileName(scoringMatrixFileName);
    if (fileName.substr(fileName.length()-4, 4).compare(".out") == 0)
        readProbMatrix();
    else{
        std::cerr << "Invalid format of the substitution matrix input file! Only .out files are accepted.\n";
        exit(1);
    }

    generateSubMatrix(this->probMatrix, this->subMatrix, this->alphabetSize, 2.0, 0.0);
}

SubstitutionMatrix::~SubstitutionMatrix(){
}

void SubstitutionMatrix::readProbMatrix(){
    
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
        if( row==alphabetSize ) break;
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
    for(size_t i=0; i<alphabetSize; ++i)
        for(size_t j=0; j<alphabetSize; ++j){
            if( i==j ) pBack[i] += probMatrix[i][j];
            else       pBack[i] += (probMatrix[i][j]/2.0f);
            if( j<=i ) sum += probMatrix[i][j];
        }
    
    const double _2sum = 2.0*sum;
    double pbsum = 0.0;
    for(size_t i=0; i<alphabetSize; ++i){
        pbsum += pBack[i];
        for(size_t j=0; j<alphabetSize; ++j)
            if( i==j ) probMatrix[i][j] = probMatrix[i][j] / sum;
            else       probMatrix[i][j] = probMatrix[i][j] / _2sum;
    }
    
    for(size_t i=0; i<alphabetSize; ++i)pBack[i] /= sum;
    
}

