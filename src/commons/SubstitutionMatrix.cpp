#include "SubstitutionMatrix.h"
#include "Util.h"
#include "Debug.h"
#include "MathUtil.h"
#include "lambda_calculator.h"
#include <cstring>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <climits>
#include "blosum62.out.h"
#include "nucleotide.out.h"



SubstitutionMatrix::SubstitutionMatrix(const char *scoringMatrixFileName_,
                                       float bitFactor, float scoreBias = -0.2) :
        scoringMatrixFileName(scoringMatrixFileName_) {
    setupLetterMapping();
    if(strcmp(scoringMatrixFileName,"blosum62.out") != 0 && strcmp(scoringMatrixFileName,"nucleotide.out")!=0 ) {
        // read amino acid substitution matrix from file
        std::string fileName(scoringMatrixFileName);
        matrixName = Util::base_name(fileName, "/\\");
        matrixName = Util::remove_extension(matrixName);
        if (fileName.substr(fileName.length() - 4, 4).compare(".out") == 0){
            std::ifstream in(fileName);
            if (in.fail()) {
                Debug(Debug::ERROR) << "Cannot read " << scoringMatrixFileName << "\n";
                EXIT(EXIT_FAILURE);
            }
            std::string str((std::istreambuf_iterator<char>(in)),
                            std::istreambuf_iterator<char>());
            int alphabetSize = readProbMatrix(str);
            if (alphabetSize < this->alphabetSize - 1) {
                this->alphabetSize = alphabetSize;
            }
            in.close();
        }
        else {
            Debug(Debug::ERROR) << "Invalid format of the substitution matrix input file! Only .out files are accepted.\n";
            EXIT(EXIT_FAILURE);
        }
    } else if(strcmp(scoringMatrixFileName,"nucleotide.out") == 0){
        std::string submat((const char *)nucleotide_out,nucleotide_out_len);
        matrixName = "nucleotide.out";
        int alphabetSize = readProbMatrix(submat);
        if (alphabetSize < this->alphabetSize - 1) {
            this->alphabetSize = alphabetSize;
        }
    } else{
        std::string submat((const char *)blosum62_out,blosum62_out_len);
        matrixName = "blosum62.out";
        int alphabetSize = readProbMatrix(submat);
        if (alphabetSize < this->alphabetSize - 1) {
            this->alphabetSize = alphabetSize;
        }
    }

    //print(probMatrix, int2aa, alphabetSize);
    generateSubMatrix(this->probMatrix, this->subMatrixPseudoCounts,
                      this->subMatrix, this->subMatrix2Bit,
                      this->alphabetSize, true, bitFactor, scoreBias);
    this->bitFactor = bitFactor;
}


bool SubstitutionMatrix::estimateLambdaAndBackground(const double **scoreMatrix,
                                                     int alphabetSize, double *pBack, double &lambda) {
    // We need to pass the parameters as 1-based pointers, hence the +1s
    // and -1s.

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
                                                   const int *int_sequence,
                                                   const int N,
                                                   float *compositionBias) {
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
        compositionBias[i] = deltaS_i;
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
                                                             int N, size_t alphabetSize, BaseMatrix *subMat) {

    const int windowSize = 40;

    float * pnul = new float[N]; // expected score of the prof ag. a random (blosum bg dist) seq
    memset(pnul, 0, sizeof(float) * N);
    float * aaSum = new float[alphabetSize];

    ProfileStates ps(alphabetSize,subMat->pBack);

    for (int pos = 0; pos < N; pos++) {
        for(size_t aa = 0; aa < alphabetSize; aa++) {
            pnul[pos] += *(profileScores + pos + N*aa) * ps.prior[aa];
        }
    }

    for (int i = 0; i < N; i++){
        const int minPos = std::max(0, (i - windowSize/2));
        const int maxPos = std::min(N, (i + windowSize/2));
        const int windowLength = maxPos - minPos;
        // negative score for the amino acids in the neighborhood of i
        memset(aaSum, 0, sizeof(float) * alphabetSize);

        for (int j = minPos; j < maxPos; j++){
            if( i == j )
                continue;
            for(size_t aa = 0; aa < alphabetSize; aa++){
                aaSum[aa] += *(profileScores + aa*N + j) - pnul[j];
            }
        }
        for(size_t aa = 0; aa < alphabetSize; aa++) {
            profileScores[i + aa*N] = static_cast<int8_t>(*(profileScores + i + N*aa) - aaSum[aa]/windowLength);
        }
    }
    delete [] aaSum;
    delete [] pnul;
}



/* Compute aa correction
   => p(a) =  ( \prod_{i=1}^L pi(a) )^(1/L)
   => p(a) = 2^[ (1/L) * \log2 ( \prod_{i=1}^L pi(a) )
   => p(a) = 2^[ (1/L) * \sum_{i=1}^L  \log2 pi(a) ]
   => p(a) = 2^[ (1/L) * \sum_{i=1}^L  \log2 ( pi(a) / f(a) )  + log2 f(a) ]
   => p(a) = f(a) * 2^[ (1/L) * \sum_{i=1}^L  S(pi, a) ]
 */

void SubstitutionMatrix::calcGlobalAaBiasCorrection(const BaseMatrix *m,
                                                    short *profileScores,
                                                    float *pNullBuffer,
                                                    const size_t profileAASize,
                                                    const int N) {
    memset(pNullBuffer, 0, sizeof(float) * N);
    const int windowSize = 40;
    for (int pos = 0; pos < N; pos++) {
        const short * subMat = profileScores + (pos * profileAASize);
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
            const short *subMat = profileScores + (j * profileAASize);
            if (i == j) {
                continue;
            }
            for (size_t aa = 0; aa < 20; aa++) {
                aaSum[aa] += subMat[aa] - pNullBuffer[j];
            }
        }
        for (size_t aa = 0; aa < 20; aa++) {
//            printf("%d\t%d\t%2.3f\t%d\n", i, (profileScores + (i * profileAASize))[aa],
//                   aaSum[aa]/windowLength,
//                   static_cast<int>((profileScores + (i * profileAASize))[aa] -  aaSum[aa]/windowLength) );
            //std::cout << i << "\t" << (profileScores + (i * profileAASize))[aa] << "\t" <<  aaSum[aa]/windowLength << "\t" <<  (profileScores + (i * profileAASize))[aa] -  aaSum[aa]/windowLength << std::endl;
            profileScores[i * profileAASize + aa] = static_cast<int>((profileScores + (i * profileAASize))[aa] -
                                                                     aaSum[aa] / windowLength);
//            avg += static_cast<int>((profileScores + (i * profileAASize))[aa] -  aaSum[aa]/windowLength);
        }
    }
//    std::cout << "avg=" << avg/(N*20) << std::endl;
}


SubstitutionMatrix::~SubstitutionMatrix() {
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
                    this->aa2int[static_cast<int>(letter)] = this->aa2int[static_cast<int>(upperLetter)];
                    break;
                case 'J':
                    this->aa2int[static_cast<int>(letter)] = this->aa2int[(int)'L'];
                    break;
                case 'U':
                case 'O':
                    this->aa2int[static_cast<int>(letter)] = this->aa2int[(int)'X'];
                    break;
                case 'Z': this->aa2int[static_cast<int>(letter)] = this->aa2int[(int)'E']; break;
                case 'B': this->aa2int[static_cast<int>(letter)] = this->aa2int[(int)'D']; break;
                default:
                    this->aa2int[static_cast<int>(letter)] = this->aa2int[(int)'X'];
                    break;
            }
        }
}

int SubstitutionMatrix::parseAlphabet(char *word, char *int2aa, int *aa2int) {
    char *charReader = word;
    int minAAInt = INT_MAX;
    // find amino acid with minimal int value
    while (isalpha(*charReader)) {
        const char aa = *charReader;
        const int intAA = aa2int[static_cast<int>(aa)];
        minAAInt = std::min(minAAInt, intAA);
        charReader++;
    }
    char minAAChar = int2aa[minAAInt];
    // do alphbet reduction
    charReader = word;
    while (isalpha(*charReader)) {
        const char aa = *charReader;
        const int intAA = aa2int[static_cast<int>(aa)];
        aa2int[static_cast<int>(aa)] = minAAInt;
        int2aa[intAA] = minAAChar;
        charReader++;
    }
    return minAAInt;
}

int SubstitutionMatrix::readProbMatrix(const std::string &matrixData) {
    std::stringstream in(matrixData);
    std::string line;
    bool probMatrixStart = false;

    char *words[256];
    int column_aa[32];
    int column_aa_sorted[32];
    int alphabetSize = 0;
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
            for (size_t i = 0; i < wordCnt; i++) {
                if (isalpha(words[i][0]) == false) {
                    Debug(Debug::ERROR) << "Probability matrix must start with alphabet header.\n";
                    EXIT(EXIT_FAILURE);
                }
                column_aa[i] = parseAlphabet(words[i], int2aa, aa2int);
            }
            alphabetSize = wordCnt;
            if (alphabetSize > 32) {
                Debug(Debug::ERROR) << "Only alphabets with up to 32 letters are allowed.\n";
                EXIT(EXIT_FAILURE);
            }
            memcpy(column_aa_sorted, column_aa, sizeof(int) * alphabetSize);
            std::sort(column_aa_sorted, column_aa_sorted + alphabetSize);
            int column_old_aa[32];
            memcpy(column_old_aa, column_aa, sizeof(int) * alphabetSize);
            std::map<int, int> mapping;
            for (int i = 0; i < alphabetSize; i++) {
                for (int j = 0; j < alphabetSize; j++) {
                    if (column_aa_sorted[i] == column_aa[j]) {
                        const char repAA = int2aa[column_aa[j]];
                        for (size_t z = 'A'; z < 'Z'; z++) {
                            aa2int[z] = (aa2int[z] == column_aa_sorted[i]) ? i : aa2int[z];
                        }
                        int2aa[i] = repAA;
                        column_aa[j] = i;
                    }
                }
            }
            probMatrixStart = true;
            continue;
        }
        if (wordCnt > 1 && probMatrixStart == true) {
            if (isalpha(words[0][0]) == false) {
                Debug(Debug::ERROR) << "First element in probability line must be an alphabet letter.\n";
                EXIT(EXIT_FAILURE);
            }
            int aa = parseAlphabet(words[0], int2aa, aa2int);
            for (int i = 0; i < alphabetSize; i++) {
                double f = strtod(words[i + 1], NULL);
                probMatrix[aa][column_aa[i]] = f; // divided by 2 because we scale bit/2 ot bit
            }
        }
    }
    bool containsX = false;
    bool xIsPositive = false;
    for (int i = 0; i < alphabetSize; i++) {
        if (column_aa[i] == aa2int[(int)'X']) {
            containsX = true;
            for (int j = 0; j < alphabetSize; j++) {
                int xIndex = aa2int[(int)'X'];
                if ((probMatrix[xIndex][j] > 0) || (probMatrix[j][xIndex] > 0)) {
                    xIsPositive = true;
                    break;
                }
            }
            break;
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
        pBack[aa2int[(int)'X']]=ANY_BACK;
    }
    if(xIsPositive == false){
        for (int i = 0; i < alphabetSize - 1; i++) {
            pBack[i] = pBack[i] * (1.0 - pBack[aa2int[(int)'X']]);
        }
    }
        // Reconstruct Probability Sab=(Pab/Pa*Pb) -> Pab = exp(Sab) * Pa * Pb
    for (int i = 0; i < alphabetSize; i++) {
        //smat[i] = smatData+((subMat.alphabetSize-1)*i);
        for (int j = 0; j < alphabetSize; j++) {
            probMatrix[i][j] = std::exp(lambda * probMatrix[i][j]) * pBack[i] * pBack[j];
        }
    }

    return alphabetSize;
}



