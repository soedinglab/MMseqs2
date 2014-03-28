#ifndef ReducedMatrix_H
#define ReducedMatrix_H
#include <fstream>
#include <string>
#include <vector>

#include "../commons/BaseMatrix.h"
#include "../commons/Debug.h"

class ReducedMatrix : public BaseMatrix {
    public: 
        ReducedMatrix(double **probMatrix, size_t reducedAlphabetSize);
        virtual ~ReducedMatrix();
        /*contains the original matrix before the alphabet reduction*/
        short ** origSubMatrix;
        int*   orig_aa2int;
        char*  orig_int2aa;
        /* size of the original alphabet*/
        size_t origAlphabetSize;

        // base class aa2int and int2aa mappings contain now:
        // aa2int: mapping aa (orig. alphabet) -> int code of the representative amino acid
        // int2aa: mapping int code (orig. alphabet) -> the representative amino acid char
        
        // reducedAlphabet contains only the "representative" amino acids
        std::vector<char>* reducedAlphabet;

    private: 

        /* The function adds two rows of a given m x n input matrix and produces a
         * m-1 x n matrix. row1 and row2 (where row1 < row2) are the rows that are
         * to be added. row1 is replaced by row1 + row2 in the output matrix.
         */
        void addTwoRows(double ** input, double ** output, size_t size, size_t row1, size_t row2 );
        /* The function adds two columns of a given m x n input matrix and produces a
         * m x n-1 matrix. col1 and col2 (where col1 < col2) are the columns that are
         * to be added. col1 is replaced by col1 + col2 in the output matrix
         */
        void addTwoColumns(double ** input, double ** output, size_t size, size_t col1, size_t col2 );

        /* Copy from array to array */
        void copyMatrix(double ** input,double ** output, size_t size);
        /* This function generates the corresponding substitution matrix given the corresponding probability
         * matrix.
         */
        void coupleBases(double ** input, double ** output, size_t size, size_t base1, size_t base2);
        /* This function calculates the mutual information of a given substitution
         * matrix (subMatrix) given that its Probability matrix (pMatrix) is also known.
         * numRows and numCols are the number of rows and columns respectively in the
         * matrices that contain meaningful information.
         *
         * Mutual information measures the information that two random variables X and Y share: 
         * it measures how much knowing one of these variables reduces uncertainty about the other.
         * (Knowing the amino acid X - how much we know about the aligned amino acid?).
         */
        double calculateMutualInformation(double ** pMatrix, double ** subMatrix, size_t size);
        /* This function finds the two best bases to couple such that we loose the minimum amount of information.
         * Returns the amount of mutual information in the best pairing.
         */
        std::pair<size_t,size_t> coupleWithBestInfo(double ** pinput, double ** pMatrix, size_t size);

};


#endif  
