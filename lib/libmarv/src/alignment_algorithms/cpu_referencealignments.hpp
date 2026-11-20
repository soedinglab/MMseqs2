#ifndef LIBMARV_CPU_REFERENCE_ALIGNMENTS_HPP
#define LIBMARV_CPU_REFERENCE_ALIGNMENTS_HPP

#include <vector>

#include "../types.hpp"


#include "../namespace.hpp"
LIBMARV_NAMESPACE_BEGIN

    //sequences must be in to ncbi converted format
    int GaplessFilter_host_protein_converted_blosum62(
        const char* seq1,
        const char* seq2,
        const int length1,
        const int length2
    ) {

        //const int NEGINFINITY = -10000;
        std::vector<int> penalty_H(2*(length2+1));

        int maxi = 0, result;
        for (int index = 0; index <= length2; index++) {
            penalty_H[index] = 0;
        }

        auto BLOSUM = libmarv::BLOSUM62_20::get2D();
        
        //std::cout << "CPU:\n";
        for (int row = 1; row <= length1; row++) {
            char seq1_char = seq1[row-1];
            char seq2_char;

            const int target_row = row & 1;
            const int source_row = !target_row;
            penalty_H[target_row*(length2+1)] = 0; //gap_open + (row-1)*gap_extend;
            for (int col = 1; col <= length2; col++) {
                const int diag = penalty_H[source_row*(length2+1)+col-1];
                seq2_char = seq2[col-1];

                const int residue = BLOSUM[seq1_char][seq2_char];
                result = std::max(0, diag + residue);
                penalty_H[target_row*(length2+1)+col] = result;
                if (result > maxi) maxi = result;
            }

            // for (int col = 1; col <= length2; col++) {
            //     printf("%2d ", penalty_H[target_row*(length2+1)+col]);
            // }
            // printf(", max %2d\n", maxi);
        }

        return maxi;
    }

    //sequences must be in to ncbi converted format
    //<score, <seq1end, seq2end>>, endposition exclusive
    std::pair<int, std::vector<std::pair<int,int>>>  GaplessFilter_host_protein_converted_blosum62_endPos(
        const char* seq1,
        const char* seq2,
        const int length1,
        const int length2
    ) {

        //const int NEGINFINITY = -10000;
        std::vector<int> penalty_H(2*(length2+1));

        int maxi = 0, result;
        for (int index = 0; index <= length2; index++) {
            penalty_H[index] = 0;
        }

        auto BLOSUM = libmarv::BLOSUM62_20::get2D();

        std::vector<std::pair<int,int>> bestEndPositions;
        
        //std::cout << "CPU:\n";
        for (int row = 1; row <= length1; row++) {
            char seq1_char = seq1[row-1];
            char seq2_char;

            const int target_row = row & 1;
            const int source_row = !target_row;
            penalty_H[target_row*(length2+1)] = 0; 
            for (int col = 1; col <= length2; col++) {
                const int diag = penalty_H[source_row*(length2+1)+col-1];
                seq2_char = seq2[col-1];

                const int residue = BLOSUM[seq1_char][seq2_char];
                result = std::max(0, diag + residue);
                penalty_H[target_row*(length2+1)+col] = result;
                if (result >= maxi){
                    if(result > maxi){
                        bestEndPositions.clear();
                    }
                    maxi = result;
                    bestEndPositions.emplace_back(row, col);
                }
            }

            // for (int col = 1; col <= length2; col++) {
            //     printf("%2d ", penalty_H[target_row*(length2+1)+col]);
            // }
            // printf(", max %2d\n", maxi);
        }

        return std::make_pair(maxi, bestEndPositions);
    }



    //sequences must be in to ncbi converted format
        //<score, <seq1end, seq2end>>, endposition exclusive
        std::pair<int, std::vector<std::pair<int,int>>> affine_local_DP_host_protein_blosum62_converted_endPos(
            const char* seq1,
            const char* seq2,
            const int length1,
            const int length2,
            const int gap_open,
            const int gap_extend
        ) {
            const int NEGINFINITY = -10000;
            std::vector<int> penalty_H(2*(length2+1));
            std::vector<int> penalty_F(2*(length2+1));

            // std::cout << "length1 " << length1 << ", length2 " << length2 << "\n";

            // for(int i = 0; i < length1; i++){
            //     std::cout << int(seq1[i]) << " ";
            // }
            // std::cout << "\n";

            // for(int i = 0; i < length2; i++){
            //     std::cout << int(seq2[i]) << " ";
            // }
            // std::cout << "\n";

            int E, F, maxi = 0, result;
            penalty_H[0] = 0;
            penalty_F[0] = NEGINFINITY;
            for (int index = 1; index <= length2; index++) {
                penalty_H[index] = 0;
                penalty_F[index] = NEGINFINITY;
            }

            auto BLOSUM = libmarv::BLOSUM62_20::get2D();

            std::vector<std::pair<int,int>> bestEndPositions;

            for (int row = 1; row <= length1; row++) {
                int seq1_char = seq1[row-1];
                int seq2_char;

                const int target_row = row & 1;
                const int source_row = !target_row;
                penalty_H[target_row*(length2+1)] = 0; //gap_open + (row-1)*gap_extend;
                penalty_F[target_row*(length2+1)] = gap_open + (row-1)*gap_extend;
                E = NEGINFINITY;
                for (int col = 1; col <= length2; col++) {
                    const int diag = penalty_H[source_row*(length2+1)+col-1];
                    const int abve = penalty_H[source_row*(length2+1)+col+0];
                    const int left = penalty_H[target_row*(length2+1)+col-1];
                    seq2_char = seq2[col-1];
                    const int residue = BLOSUM[seq1_char][seq2_char];
                    E = std::max(E+gap_extend, left+gap_open);
                    F = std::max(penalty_F[source_row*(length2+1)+col+0]+gap_extend, abve+gap_open);
                    result = std::max(0, std::max(diag + residue, std::max(E, F)));
                    penalty_H[target_row*(length2+1)+col] = result;
                    if (result >= maxi){
                        if(result > maxi){
                            bestEndPositions.clear();
                        }
                        maxi = result;
                        bestEndPositions.emplace_back(row, col);
                    }
                    penalty_F[target_row*(length2+1)+col] = F;

                    //std::cout << maxi << " ";
                }
                //std::cout << "\n";
            }
            return std::make_pair(maxi, bestEndPositions);
        }

        //sequences must be in to ncbi converted format
        int affine_local_DP_host_protein_blosum62_converted(
            const char* seq1,
            const char* seq2,
            const int length1,
            const int length2,
            const int gap_open,
            const int gap_extend
        ) {
            const int NEGINFINITY = -10000;
            std::vector<int> penalty_H(2*(length2+1));
            std::vector<int> penalty_F(2*(length2+1));

            // std::cout << "length1 " << length1 << ", length2 " << length2 << "\n";

            // for(int i = 0; i < length1; i++){
            //     std::cout << int(seq1[i]) << " ";
            // }
            // std::cout << "\n";

            // for(int i = 0; i < length2; i++){
            //     std::cout << int(seq2[i]) << " ";
            // }
            // std::cout << "\n";

            int E, F, maxi = 0, result;
            penalty_H[0] = 0;
            penalty_F[0] = NEGINFINITY;
            for (int index = 1; index <= length2; index++) {
                penalty_H[index] = 0;
                penalty_F[index] = NEGINFINITY;
            }

            auto BLOSUM = libmarv::BLOSUM62_20::get2D();

            for (int row = 1; row <= length1; row++) {
                int seq1_char = seq1[row-1];
                int seq2_char;

                const int target_row = row & 1;
                const int source_row = !target_row;
                penalty_H[target_row*(length2+1)] = 0; //gap_open + (row-1)*gap_extend;
                penalty_F[target_row*(length2+1)] = gap_open + (row-1)*gap_extend;
                E = NEGINFINITY;
                for (int col = 1; col <= length2; col++) {
                    const int diag = penalty_H[source_row*(length2+1)+col-1];
                    const int abve = penalty_H[source_row*(length2+1)+col+0];
                    const int left = penalty_H[target_row*(length2+1)+col-1];
                    seq2_char = seq2[col-1];
                    const int residue = BLOSUM[seq1_char][seq2_char];
                    E = std::max(E+gap_extend, left+gap_open);
                    F = std::max(penalty_F[source_row*(length2+1)+col+0]+gap_extend, abve+gap_open);
                    result = std::max(0, std::max(diag + residue, std::max(E, F)));
                    penalty_H[target_row*(length2+1)+col] = result;
                    if (result > maxi) maxi = result;
                    penalty_F[target_row*(length2+1)+col] = F;

                    //std::cout << maxi << " ";
                }
                //std::cout << "\n";
            }
            return maxi;
        }

LIBMARV_NAMESPACE_END


#endif