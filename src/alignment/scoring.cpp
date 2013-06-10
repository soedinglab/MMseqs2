#include "scoring.h"

char* get_int2aa(){
    char* int2aa = new char[AMINOACID_DIM];
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

    return int2aa;
}

int* get_aa2int(char* int2aa){

    int* aa2int = new int['Z'+1];
    for (int i = 0; i <= 'Z'; ++i) aa2int[i]=-1;
    for (int i = 0; i < AMINOACID_DIM; ++i){
        aa2int[int2aa[i]] = i;
    }

    return aa2int;
}

int** getScoringMatrix(std::string mFile, int* aa2int){
    
    int** sc_matrix = new int* [AMINOACID_DIM];
    for (int i = 0; i < AMINOACID_DIM; i++){
        sc_matrix[i] = new int[AMINOACID_DIM];
    }
    
    std::ifstream in(mFile.c_str());
    if(in.fail()){
        std::cerr << "Cannot read " << mFile << "\n";
        exit(1);
    }
    std::string line;
    getline( in, line );
    while (line.at(0) == '#')
        getline (in, line);
    boost::trim(line);
    std::vector<std::string> aas;
    boost::split_regex(aas, line, boost::regex(" +"));

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

        std::vector<std::string> scores;
        boost::split_regex(scores, line, boost::regex(" +"));
        for (int i = 1; i < AMINOACID_DIM+1; i++){
            sc_matrix[aa2int[aa_entries[row]]][aa2int[aa_entries[i-1]]] = atoi(scores[i].c_str());
        }
        row++;
    }

    in.close();

    delete[] aa_entries;

    return sc_matrix;
}

int** getBiasedScoringMatrix(std::string mFile, int* aa2int){

    int** sc_matrix = new int* [AMINOACID_DIM];
    for (int i = 0; i < AMINOACID_DIM; i++){
        sc_matrix[i] = new int[AMINOACID_DIM];
    }
 
    float original[AMINOACID_DIM][AMINOACID_DIM];
    float p_background[AMINOACID_DIM];
    for( size_t i=0; i<AMINOACID_DIM; ++i ) p_background[i] = 0.0f;
    
    std::ifstream in(mFile.c_str());
    if( in.fail() ) {
        std::cerr << "Cannot read " << mFile << "\n";
        exit(1);
    }
    int c      = 0;
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
        if( row==AMINOACID_DIM ) break;
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
                original[row][column] = f;
                original[column][row] = f;
                ++column;
            }
            ++row;
        }
    }
    in.close();

    float sum=0.0f;
    for(size_t i=0; i<AMINOACID_DIM; ++i)
        for(size_t j=0; j<AMINOACID_DIM; ++j){
            if( i==j ) p_background[i] += original[i][j];
            else       p_background[i] += (original[i][j]/2.0f);
            if( j<=i ) sum += original[i][j]; 
        }

    const float _2sum = 2.0*sum;    
    float pbsum = 0.0f;
    for(size_t i=0; i<AMINOACID_DIM; ++i){
        pbsum += p_background[i];
        for(size_t j=0; j<AMINOACID_DIM; ++j)
            if( i==j ) original[i][j] = original[i][j] / sum;
            else       original[i][j] = original[i][j] / _2sum;
    }

    for(size_t i=0; i<AMINOACID_DIM; ++i)p_background[i] /= sum;

    for(size_t i=0; i<AMINOACID_DIM; ++i){
        for(size_t j=0; j<AMINOACID_DIM; ++j){
            sc_matrix[i][j] = floor (BIT_FACTOR * _log2(original[i][j]/(p_background[i]*p_background[j])) + SCORING_BIAS + 0.5);
/*            if (sc_matrix[i][j] >= 0)
                std::cout << " ";
            std::cout << sc_matrix[i][j] << " ";*/
        }
//        std::cout << "\n";
    }

    return sc_matrix;
}

unsigned char getMatrixMinValue(int** sc_matrix){
    int min = 0;
    for (int i = 0; i < AMINOACID_DIM; i++){
        for (int j = 0; j < AMINOACID_DIM; j++) {
            if (sc_matrix[i][j] < min)
                min = sc_matrix[i][j];
        }
    }

    min = (-1)*min;
    unsigned char ret = (unsigned char) min;
    return ret;
}

unsigned char* getQueryProfileByte(unsigned char* query, int qlen, int** sc_matrix, unsigned char bias){

    unsigned char* qProfByte = (unsigned char*)memalign(16,AMINOACID_DIM*(qlen+15)*sizeof(unsigned char));   // allocate memory for the query profile
    int segLen = (qlen+15) / 16;

    char* int2aa = get_int2aa();

    int a,h,i,j,k;
    for (a=0; a < AMINOACID_DIM; ++a)
    {
        h = a * segLen * 16;
        for (i=0; i < segLen; ++i)
        {
            j = i;
            for (k = 0; k < 16; ++k)
            {
                if (j >= qlen){
                    qProfByte[h] = bias;
//                    std::cout << (int)bias << " ";
                }
                else {
                    qProfByte[h] = (unsigned char) (sc_matrix[query[j]][a]+(int)bias);
//                    std::cout << int2aa[(int)query[j]] << "_" << (sc_matrix[query[j]][a]+(int)bias) << " ";
                }
                ++h;
                j+=segLen;
            }
//            std::cout << "\n";
        }
    }
    return qProfByte;
}

unsigned short* getQueryProfileWord(unsigned char* query, int qlen, int** sc_matrix){

    unsigned short* qProfWord = (unsigned short*)memalign(16,AMINOACID_DIM*(qlen+7)*sizeof(unsigned short));   // allocate memory for the query profile
    int segLen = (qlen+7) / 8;

    int a,h,i,j,k;
    for (a=0; a < AMINOACID_DIM; ++a)
    {   
        h = a * segLen * 8;
        for (i=0; i < segLen; ++i)
        {   
            j = i;
            for (k = 0; k < 8; ++k)
            {   
                if (j >= qlen)
                    qProfWord[h] = 0;
                else {
                    qProfWord[h] = (unsigned short) (sc_matrix[query[j]][a]);
                }
                ++h; 
                j+=segLen;
            }
        }
    }
    return qProfWord;
}
