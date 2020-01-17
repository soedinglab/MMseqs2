#include "Sequence.h"
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "BaseMatrix.h"
#include "Timer.h"

#include <iostream>
#include <unistd.h>
#include <string>

const char* binary_name = "test_kmergeneratorperf";

int main (int argc, const char * argv[])
{
    const size_t kmer_size = 7;

    SubstitutionMatrix subMat("../../data/blosum62.out",8.0);
    std::cout << "Subustitution matrix:\n";
    std::cout << "lala matrix:\n";

 //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
//    ReducedMatrix redMat(subMat.probMatrix, 20);
 //   BaseMatrix::print(redMat.subMatrix, redMat.alphabetSize);
    std::cout << "\n";
 
    const int  testSeq[]={1,2,3,1,1,1};
    ExtendedSubstitutionMatrix extMattwo(subMat.subMatrix, 2,subMat.alphabetSize);
    ExtendedSubstitutionMatrix extMatthree(subMat.subMatrix, 3,subMat.alphabetSize);

    Indexer idx(subMat.alphabetSize,kmer_size);
    
    std::cout << "Sequence (id 0):\n";
    char* sequence = (char *) argv[1];
    std::cout << sequence << "\n\n";
 
    Sequence* s = new Sequence (10000, subMat.aa2num, subMat.num2aa,0);
    s->mapSequence(0,"lala",sequence);
    
    printf("Normal alphabet : ");
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("%c\t",subMat.num2aa[i]);
    printf("\nNormal int code: ");
    for(int i = 'A'; i<'Z';i++)
        printf("%d\t",subMat.aa2num[i]);
    
    
    std::cout << "\nInt reduced sequence:\n";
    for (int i = 0; i < s->L; i++)
        std::cout << s->int_sequence[i] << " ";
    std::cout << "\nChar reduced sequence:\n";
    
    KmerGenerator kmerGen(kmer_size,subMat.alphabetSize,10);
    kmerGen.setDivideStrategy(extMatthree.scoreMatrix, extMattwo.scoreMatrix );
    
    int* testKmer = new int[kmer_size];
    Timer timer;
    for (int x = 0; x <10; x++) {
	    s->mapSequence(x, "LALA", sequence);
        int i = 0;
        while (s->hasNextKmer(kmer_size)) {
            const int * curr_pos = s->nextKmer(kmer_size);
            printf("pos: %d kmerpos1: %d\tkmerpos2: %d\n",i++,curr_pos[0],curr_pos[1]);

            unsigned int idx_val=idx.int2index(curr_pos);
            std::cout << "Index:    " <<idx_val << "\n";
    //        std::cout << "MaxScore: " << extMattwo.scoreMatrix[idx_val]->back().first<< "\n";

            ScoreMatrix kmer_list = kmerGen.generateKmerList(curr_pos);

            std::cout << "Similar k-mer list size:" << kmer_list.elementSize << "\n\n";
            std::cout << "Similar " << kmer_size << "-mer list for pos 0:\n";
        }
    }
    std::cout << "Time: " << timer.lap() << "\n";

    return EXIT_SUCCESS;
}

