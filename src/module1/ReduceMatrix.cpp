#include "ReduceMatrix.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

ReduceMatrix::ReduceMatrix(double **probMatrix, int * aa2int,char * int2aa, 
                           size_t alphabet_size, size_t reduced_alphabet_size){
    
    this->alphabet_size = alphabet_size;
    this->reduced_alphabet_size = reduced_alphabet_size;
    this->aa2int = aa2int;
    this->int2aa = int2aa;
    this->reduced_aa2int = new int['Z'+1];
    this->reduced_int2aa = new char[alphabet_size];
    std::vector<char> reduced_int2aa_vec;
    for (size_t i = 0; i <= 'Z'; ++i) reduced_aa2int[i]=aa2int[i];
    for (size_t i = 0; i < alphabet_size; ++i){
        reduced_int2aa[i] = int2aa[i];
        reduced_int2aa_vec.push_back(int2aa[i]);
    }
    
    size_t numRows=alphabet_size;
    size_t numCols=alphabet_size;
    
    
    double ** subMatrix=new double*[numRows];
    double ** probMatrix_new=new double*[numRows];
    this->reduced_Matrix=new short *[numRows];
    for(size_t i = 0; i<numRows;i++){
        subMatrix[i]=new double[numCols];
        probMatrix_new[i]=new double[numCols];
        this->reduced_Matrix[i]=new short[numCols];
    }
    
    
    
    generateSubMatrix(probMatrix, subMatrix, numRows, numCols);
    
    size_t reduce_steps=alphabet_size-reduced_alphabet_size;
    
    
    for(size_t step = 0 ; step < reduce_steps; step++){
        // Ensuring every element is 0.
        for(size_t j = 0; j < this->alphabet_size; j++)
        {
            for(size_t i = 0; i < this->alphabet_size; i++)
            {
                probMatrix_new[i][j] = 0;
            }
        }
        //This is where the function to couple the two bases is called.
        std::pair<size_t,size_t> reduce_bases=coupleWithBestInfo(probMatrix, probMatrix_new, numRows-step, numCols-step);
        
        size_t reduced_index=reduce_bases.first;
        size_t lost_index=reduce_bases.second;
        
        char reduced_aa=reduced_int2aa_vec.at(reduced_index);
        char lost_aa   =reduced_int2aa_vec.at(lost_index);
        
        printf("Reduced aa: %c Lost aa: %c\n",reduced_aa, lost_aa);
        reduced_int2aa_vec.erase(reduced_int2aa_vec.begin()+lost_index);
        
        size_t reduced_int=this->aa2int[reduced_aa];
        size_t lost_int   =this->reduced_aa2int[lost_aa];
        for(size_t i =0; i < this->alphabet_size;i++){
            if(reduced_int2aa[i]==lost_aa){
                this->reduced_int2aa[i]=reduced_aa;
            }
        }
        for(size_t i =0; i < 'Z'; i++){
            if(reduced_aa2int[i]==lost_int){
                this->reduced_aa2int[i]=(int)reduced_int;
            }
        }
        copyMatrix(probMatrix_new,probMatrix, numRows, numCols);
        
    }
    
    // map big index to new small index
    for(size_t i2a = 0; i2a<reduced_int2aa_vec.size();i2a++){
        const char reduced_aa = reduced_int2aa_vec[i2a];
        int big_int= this->aa2int[reduced_aa];
        for(size_t i =0; i < 'Z'; i++){
            if(reduced_aa2int[i]==big_int){
                this->reduced_aa2int[i]=(int)i2a;
            }
        }
    }

    generateBiasedSubMatrix(probMatrix_new, reduced_Matrix,2.0, -0.2, numRows-reduce_steps, numCols-reduce_steps);
    for (size_t i = 0; i < numRows; i++)
    {
        delete[]subMatrix[i];
        delete[]probMatrix_new[i];
    }
    delete[]subMatrix;
    delete[]probMatrix_new;
}





ReduceMatrix::~ReduceMatrix(){
    delete[] this->reduced_aa2int;
    delete[] this->reduced_int2aa;
    for(size_t i = 0; i<alphabet_size;i++){
        delete[] reduced_Matrix[i];
    }
    delete[] reduced_Matrix;

}



void ReduceMatrix::copyMatrix(double ** input,double ** output, size_t numRows,size_t numCols){
    
    for (size_t i=0; i< numRows; i++){
        for (size_t j=0; j< numCols; j++){
            output[i][j] = input[i][j];
        }
    }
}


double ReduceMatrix::calculateMutualInformation(double ** pMatrix, double ** subMatrix, size_t numRows, size_t numCols){
    double mutualInfo = 0;
    
    for (size_t i=0; i< numRows; i++){
        for (size_t j=0; j< numCols; j++){
            mutualInfo += pMatrix[i][j]*subMatrix[i][j];
        }
    }
    
    // This is to incorporate the factor of 3 in the formula and also the fact
    // that the pMatrix being used by us is in 10E-6.
    
    return mutualInfo;
}



void ReduceMatrix::coupleBases(double ** input, double ** output, size_t numRows, size_t numCols, size_t base1, size_t base2){
    
    double ** temp=new double *[this->alphabet_size];
    //To ensure every element of temp is set to 0.
    for(size_t i = 0; i < this->alphabet_size; i++)
    {
        temp[i]=new double[this->alphabet_size];
        
        for(size_t j = 0; j < this->alphabet_size; j++)
        {
            temp[i][j] = 0;
        }
    }
    
    //Add corresponding columns first.
    addTwoColumns(input, temp, numRows, numCols, base1, base2);
    
    //Add the corresponding rows.
    addTwoRows(temp, output, numRows, numCols, base1, base2);
    
    for (size_t i = 0; i < this->alphabet_size; i++)
    {
        delete[]temp[i];
    }
    delete[]temp;
}



std::pair<size_t,size_t> ReduceMatrix::coupleWithBestInfo(double ** pinput, double ** pMatrix, size_t numRows, size_t numCols){
    double bestInfo=0;
    size_t besti = 0, bestj = 0;
    
    
    // To store the temporary substitution matrix after each coupling.
    double ** tempsub=new double *[numRows];
    // To store the temporary probability matrix after each coupling.
    double ** tempp=new double *[numRows];
    
    for(size_t i = 0; i<numRows;i++){
        tempsub[i]=new double [numCols];
        tempp[i]=new double [numCols];
        
    }
    
    for (size_t i=0; i< numRows; i++){
        
        // To store the mutual information of the matrix.
        double temp = 0;
        
        for (size_t j=i+1; j< numRows; j++){
            coupleBases(pinput, tempp, numRows, numCols, i, j);
            
            // Generate the new substitution matrix after two bases have been coupled.
            generateSubMatrix(tempp, tempsub, numRows -1, numCols -1);
            
            // Storing mutual information in temp.
            temp = calculateMutualInformation(tempp, tempsub, numRows-1, numCols-1);
            
            if (temp > bestInfo) {bestInfo = temp; besti = i; bestj = j;}
            
            std::cout << " i = " << i << "; j = " << j << " info " << temp << '\n';
        }
        
    }
    
    
    static int i =0;
    std::cout << "Called: " << i++ << '\n';
    std::cout << "Chosen i and j" << besti << " " << bestj << '\n';
    // Finally coupling the best option.
    coupleBases(pinput, pMatrix, numRows, numCols, besti, bestj);
    for (size_t i = 0; i < numRows; i++)
    {
        delete[]tempsub[i];
        delete[]tempp[i];
    }
    delete[]tempsub;
    delete[]tempp;
    return std::make_pair<size_t,size_t>(besti,bestj);
}





void ReduceMatrix::genProbBaseArray(double ** pmatrix, double * prob, size_t numRows, size_t numCols){
    
    for (size_t i = 0; i < numRows; i++){
        for (size_t j = 0; j < numCols; j++){
            prob[i] += pmatrix[i][j];
        }
    }
    
}

void ReduceMatrix::generateSubMatrix(double ** pmatrix, double ** subMatrix, size_t numRows, size_t numCols){
    
    double prob[numRows];
    genProbBaseArray(pmatrix, prob, numRows, numCols);
    for (size_t i = 0; i < numRows; i++){
        for (size_t j = 0; j < numCols; j++){
            double temp = pmatrix[i][j]/(prob[i]*prob[j]);
            subMatrix[i][j] = _log2 (temp);
        }
    }
    
}


void ReduceMatrix::generateBiasedSubMatrix(double ** pmatrix, short ** subMatrix,
                                           double bitFactor, double scoringBias,
                                           size_t numRows, size_t numCols){
    
    double prob[numRows];
    genProbBaseArray(pmatrix, prob, numRows, numCols);
    for (size_t i = 0; i < numRows; i++){
        for (size_t j = 0; j < numCols; j++){
            subMatrix[i][j] = (short)floor (bitFactor * _log2(pmatrix[i][j]/(prob[i]*prob[j])) + scoringBias + 0.5);
        }
    }
    
}


void ReduceMatrix::addTwoColumns(double ** input, double ** output, size_t numRows , size_t numCols, size_t col1, size_t col2 ){
    
    for(size_t i = 0; i < numRows; i++){
        
        //copy the same data until col2 (excluding)
        for(size_t j = 0; j < col2; j++) {
            output[i][j] = input[i][j];
        }
        
        //Add col2 to col1.
        output[i][col1] = input[i][col1] + input[i][col2];
        
        //shift the rest of the columns left by 1 place.
        for(size_t j = col2 ; j < numCols-1; j++) {
            output[i][j] = input[i][j+1];
        }
    }
}



void ReduceMatrix::addTwoRows(double ** input, double ** output, size_t numRows , 
                              size_t numCols, size_t row1, size_t row2 )
{
    
    //copy the same data until row2 (excluding)
    for(size_t i = 0; i < row2; i++){
        for(size_t j = 0; j < numCols; j++)
            output[i][j] = input[i][j];
    }
    
    //add row2 to row1
    for(size_t j = 0; j < numCols; j++){
        output[row1][j] = input[row1][j] + input[row2][j];
    }
    
    //shift the rest of rows above by 1 place
    for(size_t i = row2; i < numRows-1; i++){
        for(size_t j = 0; j < numCols; j++)
            output[i][j] = input[i+1][j];
    }
    
    
}
