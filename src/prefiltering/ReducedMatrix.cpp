#include <cmath>
#include "ReducedMatrix.h"
#include "Util.h"

ReducedMatrix::ReducedMatrix(double **probMatrix, float ** rMatrix,
                             unsigned char* aa2num, char* num2aa,
                             size_t orgAlphabetSize,
                             size_t reducedAlphabetSize, float bitFactor){
    if(reducedAlphabetSize >= orgAlphabetSize) {
        Debug(Debug::ERROR) << "Reduced alphabet has to be smaller than the original one!";
        EXIT(EXIT_FAILURE);
    }
    initMatrixMemory(orgAlphabetSize);
    // swap the matrix and alphabet mappings
    this->origAlphabetSize = orgAlphabetSize;
    this->orig_aa2num = new unsigned char[UCHAR_MAX];
    memcpy(orig_aa2num, aa2num, sizeof(unsigned char) * UCHAR_MAX);
    this->orig_num2aa = new char[orgAlphabetSize];
    memcpy(orig_num2aa, num2aa, sizeof(char) * orgAlphabetSize);

    for(size_t i = 0; i < this->origAlphabetSize; i++) {
        for (size_t j = 0; j < this->origAlphabetSize; j++) {
            this->probMatrix[i][j] = probMatrix[i][j];
        }
    }
    // initialize new matrices and alphabet mappings
    this->alphabetSize = reducedAlphabetSize;
    for (size_t i = 0; i < UCHAR_MAX; ++i) { this->aa2num[i] = orig_aa2num[i]; };
    for (size_t i = 0; i < origAlphabetSize; ++i){
        this->num2aa[i] = orig_num2aa[i];
        reducedAlphabet.push_back(this->num2aa[i]);
    }

    double ** subMatrix_tmp=new double*[origAlphabetSize-1];
    double ** probMatrix_new=new double*[origAlphabetSize-1];
    for(size_t i = 0; i<origAlphabetSize-1;i++){
        subMatrix_tmp[i]=new double[origAlphabetSize-1];
        probMatrix_new[i]=new double[origAlphabetSize-1];
    }

    generateSubMatrix(this->probMatrix, subMatrix_tmp, rMatrix,  origAlphabetSize-1, false);

//    double info = calculateMutualInformation(probMatrix, subMatrix_tmp, origAlphabetSize-1);
//    Debug(Debug::INFO) << "20 " << info << "\n";
    //print(subMatrix, origAlphabetSize -1,  )

    size_t reduce_steps = origAlphabetSize - reducedAlphabetSize;


    for(size_t step = 0; step < reduce_steps; step++){
        // Ensuring every element is 0.
        for(size_t j = 0; j < this->origAlphabetSize-1; j++)
        {
            for(size_t i = 0; i < this->origAlphabetSize-1; i++)
            {
                probMatrix_new[i][j] = 0;
            }
        }
        //This is where the function to couple the two bases is called.
        std::pair<int,int> reduce_bases=coupleWithBestInfo(this->probMatrix, probMatrix_new, rMatrix,  origAlphabetSize-1-step);

        int reduced_index=reduce_bases.first;
        int lost_index=reduce_bases.second;

        char reduced_aa= reducedAlphabet.at(reduced_index);
        char lost_aa   = reducedAlphabet.at(lost_index);

        // Debug(Debug::INFO)  << lost_aa  << " -> " << reduced_aa << "\n";
        reducedAlphabet.erase(reducedAlphabet.begin() + lost_index);

        int reduced_int= this->orig_aa2num[(int)reduced_aa];
        int lost_int   = this->aa2num[static_cast<int>(lost_aa)];

        for (size_t i = 0; i < this->origAlphabetSize; i++) {
            if(this->num2aa[i]==lost_aa){
                this->num2aa[i]=reduced_aa;
            }
        }
        for (int i =0; i < UCHAR_MAX; i++) {
            if (this->aa2num[i]==lost_int) {
                this->aa2num[i] = (int) reduced_int;
            }
        }
        copyMatrix(probMatrix_new, this->probMatrix, origAlphabetSize-1);
    }

    // map big index to new small index
    Debug(Debug::INFO) << "Reduced amino acid alphabet: ";
    unsigned char* aa2num_new = new unsigned char[UCHAR_MAX+1];
    for (int i = 0; i <= UCHAR_MAX; ++i) {
        aa2num_new[i] = UCHAR_MAX;
    }
    char* num2aa_new = new char[origAlphabetSize];
    for(size_t i = 0; i<reducedAlphabet.size(); i++){
        const char representative_aa = reducedAlphabet.at(i);
        Debug(Debug::INFO) << "(" << representative_aa;
        for(size_t j =0; j < UCHAR_MAX; j++){
            if(this->aa2num[static_cast<int>(j)] == this->aa2num[static_cast<int>(representative_aa)]){
                if(j>=65 && j <=90 && static_cast<char>(j) != representative_aa && representative_aa != 'X'){ // only upper case letters
                    Debug(Debug::INFO) << " " << static_cast<char>(j);
                }
                aa2num_new[j] = i;
            }
        }
        Debug(Debug::INFO) << ") ";
        num2aa_new[i] = representative_aa;
    }
    Debug(Debug::INFO) << "\n";

    // compute background
    computeBackground(probMatrix_new, pBack, alphabetSize, true);

    // compute X background
    for (int i = 0; i < alphabetSize - 1; i++) {
        pBack[i] = pBack[i] * (1.0 - pBack[aa2num[(int)'X']]);
    }

    double * origpBack=new double[origAlphabetSize];
    computeBackground(probMatrix, origpBack, origAlphabetSize, true);
    // copy old X state
    for (int i = 0; i < this->alphabetSize; i++) {
        int oldIndex = aa2num[(int)num2aa_new[i]];
        double Pab = probMatrix[oldIndex][origAlphabetSize-1] / ( origpBack[oldIndex] * origpBack[origAlphabetSize-1]);
        probMatrix_new[alphabetSize-1][i] = Pab * pBack[i] * pBack[alphabetSize-1];
        probMatrix_new[i][alphabetSize-1] = Pab * pBack[alphabetSize-1] * pBack[i];
    }
    delete [] origpBack;
    generateSubMatrix(probMatrix_new, rMatrix, this->subMatrix, alphabetSize, true, bitFactor, 0.0);


    delete[] this->num2aa;
    delete[] this->aa2num;

    this->num2aa = num2aa_new;
    this->aa2num = aa2num_new;


    setupLetterMapping();
    for (size_t i = 0; i < origAlphabetSize-1; i++) {
        delete[] probMatrix_new[i];
        delete[] subMatrix_tmp[i];
    }
    delete[] subMatrix_tmp;
    delete[] probMatrix_new;
}

ReducedMatrix::~ReducedMatrix(){
    delete[] orig_num2aa;
    delete[] orig_aa2num;
}

void ReducedMatrix::copyMatrix(double ** input,double ** output, size_t size){
    for (size_t i=0; i< size; i++){
        for (size_t j=0; j< size; j++){
            output[i][j] = input[i][j];
        }
    }
}

double ReducedMatrix::calculateMutualInformation(double ** pMatrix, double ** subMatrix, size_t size){
    double mutualInfo = 0;

    for (size_t i=0; i< size; i++){
        for (size_t j=0; j< size; j++){
            mutualInfo += pMatrix[i][j]*subMatrix[i][j];
        }
    }
    return mutualInfo;
}

void ReducedMatrix::coupleBases(double ** input, double ** output, size_t size, size_t base1, size_t base2){

    double ** temp=new double *[this->origAlphabetSize-1];
    //To ensure every element of temp is set to 0.
    for(size_t i = 0; i < this->origAlphabetSize-1; i++)
    {
        temp[i]=new double[this->origAlphabetSize-1];

        for(size_t j = 0; j < this->origAlphabetSize-1; j++)
        {
            temp[i][j] = 0;
        }
    }

    //Add corresponding columns first.
    addTwoColumns(input, temp, size, base1, base2);

    //Add the corresponding rows.
    addTwoRows(temp, output, size, base1, base2);

    for (size_t i = 0; i < this->origAlphabetSize-1; i++)
    {
        delete [] temp[i];
    }
    delete [] temp;
}

std::pair<size_t,size_t> ReducedMatrix::coupleWithBestInfo(double ** pinput, double ** pMatrix, float ** rMatrix, size_t size){
    double bestInfo = 0;
    size_t besti = 0, bestj = 0;


    // To store the temporary substitution matrix after each coupling.
    double ** tempsub=new double *[size];
    // To store the temporary probability matrix after each coupling.
    double ** tempp=new double *[size];

    for(size_t i = 0; i < size; i++){
        tempsub[i]=new double [size];
        tempp[i]=new double [size];
    }

    for (size_t i=0; i < size; i++){

        // To store the mutual information of the matrix.
        double temp = 0;

        for (size_t j=i+1; j  < size; j++){
            coupleBases(pinput, tempp, size, i, j);
            // Generate the new substitution matrix after two bases have been coupled.
            generateSubMatrix(tempp, tempsub, rMatrix, size-1, false);
            // Storing mutual information in temp.
            temp = calculateMutualInformation(tempp, tempsub, size-1);
            if (temp > bestInfo) {bestInfo = temp; besti = i; bestj = j;}
//            Debug(Debug::INFO) << " i = " << i << "; j = " << j << " info " << temp << '\n';
        }
    }
    //Debug(Debug::INFO) << (size-1) <<  " " << bestInfo << "\n";
    // Finally coupling the best option.
    coupleBases(pinput, pMatrix, size, besti, bestj);
    for (size_t i = 0; i < size; i++)
    {
        delete[]tempsub[i];
        delete[]tempp[i];
    }
    delete[]tempsub;
    delete[]tempp;
    return std::make_pair(besti,bestj);
}

void ReducedMatrix::addTwoColumns(double ** input, double ** output, size_t size, size_t col1, size_t col2 ){

    for(size_t i = 0; i < size; i++){
        //copy the same data until col2 (excluding)
        for(size_t j = 0; j < col2; j++) {
            output[i][j] = input[i][j];
        }

        //Add col2 to col1.
        output[i][col1] = (input[i][col1] + input[i][col2]);

        //shift the rest of the columns left by 1 place.
        for(size_t j = col2 ; j < size-1; j++) {
            output[i][j] = input[i][j+1];
        }

        // set the last column to 0.0
        for (size_t i = 0; i < size; i++)
            output[i][size-1] = 0.0;
    }
}

void ReducedMatrix::addTwoRows(double ** input, double ** output, size_t size, size_t row1, size_t row2 )
{
    //copy the same data until row2 (excluding)
    for(size_t i = 0; i < row2; i++){
        for(size_t j = 0; j < size; j++)
            output[i][j] = input[i][j];
    }

    //add row2 to row1
    for(size_t j = 0; j < size; j++){
        output[row1][j] = (input[row1][j] + input[row2][j]);
    }

    //shift the rest of rows above by 1 place
    for(size_t i = row2; i < size-1; i++){
        for(size_t j = 0; j < size; j++)
            output[i][j] = input[i+1][j];
    }

    for (size_t j = 0; j < size; j++)
        output[size-1][j] = 0.0;
}
