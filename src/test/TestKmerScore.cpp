//
// Created by mad on 10/26/15.
//

#include <SubstitutionMatrix.h>
#include <Sequence.h>

int main (int argc, const char * argv[]) {

    const size_t kmer_size = 6;

    SubstitutionMatrix subMat("/Users/mad/Documents/workspace/mmseqs/data/blosum62.out", 8.0, 0);
    std::cout << "Subustitution matrix:\n";
    SubstitutionMatrix::print(subMat.subMatrix, subMat.int2aa, subMat.alphabetSize);

    char * ref   = "GKILII";
    Sequence refSeq(1000, subMat.aa2int, subMat.int2aa, 0,kmer_size, false);
    refSeq.mapSequence(0, 0, ref);

    char * similar   = "GKVLYL";
    Sequence similarSeq(1000, subMat.aa2int, subMat.int2aa, 0, kmer_size, false);
    similarSeq.mapSequence(0, 1, similar);


    short score = 0;
        for(size_t i = 0; i < kmer_size; i++){
            score += subMat.subMatrix[refSeq.int_sequence[i]][similarSeq.int_sequence[i]];
        }
    std::cout << score << std::endl;

    return 0;
}
