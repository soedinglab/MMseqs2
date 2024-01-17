#include "SubstitutionMatrix.h"
#include "Util.h"
#include "Debug.h"
#include "lambda_calculator.h"


#include <cstring>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <climits>

SubstitutionMatrix::SubstitutionMatrix(const char *filename, float bitFactor, float scoreBias) : bitFactor(bitFactor) {
    std::pair<std::string, std::string> parsedMatrix = BaseMatrix::unserialize(filename);
    if(parsedMatrix.second != "") {
        // the filename can contain the substituion matrix
        // SUBMATNAME.out:DATA
        // this is used for index databases
        matrixName = parsedMatrix.first;
        matrixData = parsedMatrix.second;
    } else {
        // read amino acid substitution matrix from file
        std::string fileName(parsedMatrix.first.c_str());
        matrixName = Util::base_name(fileName, "/\\");
        if (fileName.length() < 4 || fileName.substr(fileName.length() - 4, 4).compare(".out") != 0) {
            Debug(Debug::ERROR) << "Invalid format of the substitution matrix input file! Only .out files are accepted.\n";
            EXIT(EXIT_FAILURE);
        }
        std::ifstream in(fileName);
        if (in.fail()) {
            Debug(Debug::ERROR) << "Cannot read " << filename << "\n";
            EXIT(EXIT_FAILURE);
        }
        matrixData = std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        in.close();
    }

    std::pair<int, bool> alphSizeAndX = setAaMappingDetectAlphSize(matrixData);
    alphabetSize = alphSizeAndX.first;
    if(alphabetSize == -1){
        Debug(Debug::ERROR) << "Could not estimate alphabet size.\n";
        EXIT(EXIT_FAILURE);
    }
    initMatrixMemory(alphabetSize);
    readProbMatrix(matrixData, alphSizeAndX.second);
    if(mappingHasAminoAcidLetters()){
        setupLetterMapping();
    }else {
        for (int letter = 0; letter < UCHAR_MAX; letter++) {
            char upperLetter = toupper(static_cast<char>(letter));
            aa2num[letter] = (aa2num[static_cast<unsigned char>(upperLetter)] == UCHAR_MAX)
                             ? alphabetSize-1 : aa2num[static_cast<int>(upperLetter)];
        }
    }

    //print(probMatrix, num2aa, alphabetSize);
    generateSubMatrix(probMatrix, subMatrixPseudoCounts, subMatrix, alphabetSize, true, bitFactor, scoreBias);
}


bool SubstitutionMatrix::estimateLambdaAndBackground(const double **scoreMatrix,
                                                     int alphabetSize, double *pBack, double &lambda) {
    // We need to pass the parameters as 1-based pointers, hence the +1s and -1s.
    std::vector<double> cells(alphabetSize * (alphabetSize + 1));
    std::vector<const double *> pointers(alphabetSize + 1);

    for (int i = 0; i < alphabetSize; ++i) {
        pointers[i + 1] = &cells[i * alphabetSize];
        for (int j = 0; j < alphabetSize; ++j) {
            cells[i * alphabetSize + j + 1] = scoreMatrix[i][j];
        }
    }

    std::vector<double> letterProbs1(alphabetSize, 0);
    std::vector<double> letterProbs2(alphabetSize, 0);

    lambda = calculate_lambda(&pointers[0], alphabetSize,
                              &letterProbs1[0] - 1,
                              &letterProbs2[0] - 1);

    for (int i = 0; i < alphabetSize; i++) {
        pBack[i] = letterProbs1[i];
    }

    if (lambda < 0)
        return false; //bad
    else
        return true; //good
}


void SubstitutionMatrix::calcLocalAaBiasCorrection(const BaseMatrix *m,
                                                   const unsigned char *int_sequence,
                                                   const int N,
                                                   float *compositionBias,
                                                   float scale) {
    const int windowSize = 40;
    for (int i = 0; i < N; i++) {
        const int minPos = std::max(0, (i - windowSize / 2));
        const int maxPos = std::min(N, (i + windowSize / 2));
        const int windowLength = maxPos - minPos;

        // negative score for the amino acids in the neighborhood of i
        int sumSubScores = 0;
        short *subMat = m->subMatrix[int_sequence[i]];
        for (int j = minPos; j < maxPos; j++) {
            sumSubScores += subMat[int_sequence[j]];
        }

        // remove own amino acid
        sumSubScores -= subMat[int_sequence[i]];
        float deltaS_i = (float) sumSubScores;
        // negative avg.
        deltaS_i /= -1.0 * static_cast<float>(windowLength);
        // positive score for the background score distribution for i
        for (int a = 0; a < m->alphabetSize; a++) {
            deltaS_i += m->pBack[a] * static_cast<float>(subMat[a]);
        }
        compositionBias[i] = scale * deltaS_i;
//        std::cout << i << " " << compositionBias[i] << std::endl;
    }
}


void SubstitutionMatrix::calcProfileProfileLocalAaBiasCorrection(short *profileScores,
                                                                 const size_t profileAASize,
                                                                 const int N, size_t alphabetSize) {

    const int windowSize = 40;

    float * pnul  = new float[alphabetSize];
    float * aaSum = new float[alphabetSize];

    memset(pnul, 0, sizeof(float) * alphabetSize);


    for (int pos = 0; pos < N; pos++) {
        const short * subMat = profileScores + (pos * profileAASize);
        for(size_t aa = 0; aa < alphabetSize; aa++) {
            pnul[aa] += subMat[aa]  ;
        }
    }
    for(size_t aa = 0; aa < alphabetSize; aa++)
        pnul[aa] /= N;
    for (int i = 0; i < N; i++){
        const int minPos = std::max(0, (i - windowSize/2));
        const int maxPos = std::min(N, (i + windowSize/2));
        const int windowLength = maxPos - minPos;
        // negative score for the amino acids in the neighborhood of i
        memset(aaSum, 0, sizeof(float) * alphabetSize);

        for (int j = minPos; j < maxPos; j++){
            const short * subMat = profileScores + (j * profileAASize);
            if( i == j )
                continue;
            for(size_t aa = 0; aa < alphabetSize; aa++){
                aaSum[aa] += subMat[aa] - pnul[aa];
            }
        }
        for(size_t aa = 0; aa < alphabetSize; aa++) {
            profileScores[i*profileAASize + aa] = static_cast<int>((profileScores + (i * profileAASize))[aa] - aaSum[aa]/windowLength);
        }
    }
    delete [] aaSum;
    delete [] pnul;
}

void SubstitutionMatrix::calcProfileProfileLocalAaBiasCorrectionAln(int8_t *profileScores,
                                                                    unsigned int N, size_t alphabetSize, BaseMatrix *subMat) {

    const int windowSize = 40;

    float * pnul = new float[N]; // expected score of the prof ag. a random (blosum bg dist) seq
    memset(pnul, 0, sizeof(float) * N);
    float * aaSum = new float[alphabetSize];

    ProfileStates ps(alphabetSize,subMat->pBack);

    for (unsigned int pos = 0; pos < N; pos++) {
        for(size_t aa = 0; aa < alphabetSize; aa++) {
            pnul[pos] += *(profileScores + pos + N*aa) * ps.prior[aa];
        }
    }

    for (unsigned int i = 0; i < N; i++){
        const int minPos = std::max(0, ((int)i - windowSize/2));
        const unsigned int maxPos = std::min(N, (i + windowSize/2));
        const int windowLength = maxPos - minPos;
        // negative score for the amino acids in the neighborhood of i
        memset(aaSum, 0, sizeof(float) * alphabetSize);

        for (unsigned int j = minPos; j < maxPos; j++){
            if (i == j) {
                continue;
            }
            for (size_t aa = 0; aa < alphabetSize; aa++) {
                aaSum[aa] += *(profileScores + aa*N + j) - pnul[j];
            }
        }
        for (size_t aa = 0; aa < alphabetSize; aa++) {
            profileScores[i + aa*N] = static_cast<int8_t>(*(profileScores + i + N*aa) - aaSum[aa]/windowLength);
        }
    }
    delete[] aaSum;
    delete[] pnul;
}



/* Compute aa correction
   => p(a) =  ( \prod_{i=1}^L pi(a) )^(1/L)
   => p(a) = 2^[ (1/L) * \log2 ( \prod_{i=1}^L pi(a) )
   => p(a) = 2^[ (1/L) * \sum_{i=1}^L  \log2 pi(a) ]
   => p(a) = 2^[ (1/L) * \sum_{i=1}^L  \log2 ( pi(a) / f(a) )  + log2 f(a) ]
   => p(a) = f(a) * 2^[ (1/L) * \sum_{i=1}^L  S(pi, a) ]
 */

void SubstitutionMatrix::calcGlobalAaBiasCorrection(const BaseMatrix *m,
                                                    char *profileScores,
                                                    float *pNullBuffer,
                                                    const size_t profileAASize,
                                                    const int N) {
    memset(pNullBuffer, 0, sizeof(float) * N);
    const int windowSize = 40;
    for (int pos = 0; pos < N; pos++) {
        const char * subMat = profileScores + (pos * profileAASize);
        for(size_t aa = 0; aa < 20; aa++) {
            pNullBuffer[pos] += m->pBack[aa] * static_cast<float>(subMat[aa]);
        }
    }
//    for(size_t aa = 0; aa < 20; aa++)
//        pNullBuffer[aa] /= N;
    for (int i = 0; i < N; i++) {
        const int minPos = std::max(0, (i - windowSize / 2));
        const int maxPos = std::min(N, (i + windowSize / 2));
        const int windowLength = maxPos - minPos;
        // negative score for the amino acids in the neighborhood of i
        float aaSum[20];
        memset(aaSum, 0, sizeof(float) * 20);

        for (int j = minPos; j < maxPos; j++) {
            const char *subMat = profileScores + (j * profileAASize);
            if (i == j) {
                continue;
            }
            for (size_t aa = 0; aa < 20; aa++) {
                aaSum[aa] += subMat[aa] - pNullBuffer[j];
            }
        }
        for (size_t aa = 0; aa < 20; aa++) {
            profileScores[i * profileAASize + aa] = static_cast<int>((profileScores + (i * profileAASize))[aa] -
                                                                     aaSum[aa] / windowLength);
//            avg += static_cast<int>((profileScores + (i * profileAASize))[aa] -  aaSum[aa]/windowLength);
        }
    }
}


SubstitutionMatrix::~SubstitutionMatrix() {}

bool SubstitutionMatrix::mappingHasAminoAcidLetters(){
    std::string lettersToCheck = "ATGCDEFHIKLMNPQRSVWYX";
    size_t cnt = 0;
    for(size_t i = 0; i < lettersToCheck.size(); i++){
        cnt += (aa2num[static_cast<int>(lettersToCheck[i])] != UCHAR_MAX);
    }
    return (cnt == lettersToCheck.size());
}

void SubstitutionMatrix::setupLetterMapping(){
    for(int letter = 0; letter < UCHAR_MAX; letter++){
        char upperLetter = toupper(static_cast<char>(letter));
        switch(upperLetter){
            case 'A':
            case 'T':
            case 'G':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
            case 'H':
            case 'I':
            case 'K':
            case 'L':
            case 'M':
            case 'N':
            case 'P':
            case 'Q':
            case 'R':
            case 'S':
            case 'V':
            case 'W':
            case 'Y':
            case 'X':
                this->aa2num[static_cast<int>(letter)] = this->aa2num[static_cast<int>(upperLetter)];
                break;
            case 'J':
                this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'L'];
                break;
            case 'U':
            case 'O':
                this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'X'];
                break;
            case 'Z': this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'E']; break;
            case 'B': this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'D']; break;
            default:
                this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'X'];
                break;
        }
    }
}

int SubstitutionMatrix::parseAlphabet(char *word, char *num2aa, int *aa2num) {
    char *charReader = word;
    int minAAInt = INT_MAX;
    // find amino acid with minimal int value
    while (isalpha(*charReader)) {
        const char aa = *charReader;
        const int intAA = aa2num[static_cast<int>(aa)];
        minAAInt = std::max(minAAInt, intAA);
        charReader++;
    }
    if(minAAInt==-1){

    }
    char minAAChar = num2aa[minAAInt];
    // do alphabet reduction
    charReader = word;
    while (isalpha(*charReader)) {
        const char aa = *charReader;
        const int intAA = aa2num[static_cast<int>(aa)];
        aa2num[static_cast<int>(aa)] = minAAInt;
        num2aa[intAA] = minAAChar;
        charReader++;
    }
    return minAAInt;
}

void SubstitutionMatrix::readProbMatrix(const std::string &matrixData, const bool containsX) {
    std::stringstream in(matrixData);
    std::string line;
    bool probMatrixStart = false;
    const char *words[256];
    bool hasLambda = false;
    bool hasBackground = false;
    while (in.good()) {
        getline(in, line);
        size_t wordCnt = Util::getWordsOfLine((char *) line.c_str(), words, 256);
        // skip comments
        if (line[0] == '#') {
            if (line.find("# Background (precomputed optional):") == 0) {
                for (size_t i = 4; i < wordCnt; i++) {
                    double f = strtod(words[i], NULL);
                    pBack[i-4] = f;
                }
                hasBackground = true;
            }
            if (line.find("# Lambda     (precomputed optional):") == 0) {
                double f = strtod(words[4], NULL);
                lambda = f;
                hasLambda = true;
            }
            continue;
        }

        if (wordCnt > 1 && probMatrixStart == false) {
            probMatrixStart = true;
            continue;
        }
        if (wordCnt > 1 && probMatrixStart == true) {
            if (isalpha(words[0][0]) == false) {
                Debug(Debug::ERROR) << "First element in probability line must be an alphabet letter.\n";
                EXIT(EXIT_FAILURE);
            }
            int aa = static_cast<int>(aa2num[toupper(words[0][0])]);
            for (int i = 0; i < alphabetSize; i++) {
                double f = strtod(words[i + 1], NULL);
                probMatrix[aa][i] = f; // divided by 2 because we scale bit/2 ot bit
            }
        }
    }
    bool xIsPositive = false;
    if( containsX == true ){
        for (int j = 0; j < alphabetSize; j++) {
            int xIndex = static_cast<int>(aa2num[static_cast<int>('X')]);
            if ((probMatrix[xIndex][j] > 0) || (probMatrix[j][xIndex] > 0)) {
                xIsPositive = true;
                break;
            }
        }
    }


    if (containsX == false) {
        Debug(Debug::ERROR) << "Please add X to your substitution matrix.\n";
        EXIT(EXIT_FAILURE);
    }

    if(hasLambda == false || hasBackground == false){
        if (estimateLambdaAndBackground(const_cast<const double **>(probMatrix), alphabetSize - ((xIsPositive) ? 0 : 1),
                                        pBack, lambda) == false) {
            Debug(Debug::ERROR) << "Computing inverse of substitution matrix failed\n";
            EXIT(EXIT_FAILURE);
        }
        pBack[static_cast<int>(aa2num[static_cast<int>('X')])]=ANY_BACK;
    }
    if(xIsPositive == false){
        for (int i = 0; i < alphabetSize - 1; i++) {
            pBack[i] = pBack[i] * (1.0 - pBack[static_cast<int>(aa2num[static_cast<int>('X')])]);
        }
    }
    // Reconstruct Probability Sab=(Pab/Pa*Pb) -> Pab = exp(Sab) * Pa * Pb
    for (int i = 0; i < alphabetSize; i++) {
        //smat[i] = smatData+((subMat.alphabetSize-1)*i);
        for (int j = 0; j < alphabetSize; j++) {
            probMatrix[i][j] = std::exp(lambda * probMatrix[i][j]) * pBack[i] * pBack[j];
        }
    }
}

std::pair<int, bool> SubstitutionMatrix::setAaMappingDetectAlphSize(std::string &matrixData){
    std::stringstream in(matrixData);
    std::string line;
    const char *words[256];
    int alphabetSize = 0;
    bool containsX;
    while (in.good()) {
        getline(in, line);
        size_t wordCnt = Util::getWordsOfLine((char *) line.c_str(), words, 256);

        if (line[0] == '#') {
            continue;
        }
        if (wordCnt > 1) {
            for (size_t i = 0; i < wordCnt; i++) {
                if (isalpha(words[i][0]) == false) {
                    Debug(Debug::ERROR) << "Probability matrix must start with alphabet header.\n";
                    EXIT(EXIT_FAILURE);
                }
                int aa = toupper(words[i][0]);
                aa2num[aa] = static_cast<unsigned char>(i);
                num2aa[i] = aa;
                if (aa == 'X') {
                    containsX = true;
                }
//                column_aa[i] = parseAlphabet(words[i], num2aa, aa2num);
            }
            alphabetSize = wordCnt;
            return std::make_pair(alphabetSize, containsX);
        }
    }
    return std::make_pair(-1, false);
}




