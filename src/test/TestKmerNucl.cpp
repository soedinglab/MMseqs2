#include "Sequence.h"
#include "Indexer.h"
#include "NucleotideMatrix.h"
#include "BaseMatrix.h"
#include "Parameters.h"
#include "Orf.h"

const char* binary_name = "test_kmernucl";

std::string kmerToSting(size_t idx, int size) {
    char output[32];
    char nuclCode[4] = {'A','C','T','G'};
    for (int i = size-1; i >= 0; i--) {
        output[i] = nuclCode[idx & 3];
        idx = idx >> 2;
    }
    output[size]='\0';
    return output;
}

int main (int, const char**) {
    const size_t kmer_size = 6;

    Parameters& par = Parameters::getInstance();
    par.initMatrices();
    NucleotideMatrix subMat(par.scoringMatrixFile.nucleotides, 2.0, -0.0f);

    Indexer idx((size_t)subMat.alphabetSize, kmer_size);
    std::cout << "Sequence: ";
    const char* sequence = "GATACAGATACAGATACAGATACA";
    std::cout << sequence << "\n";
    Sequence* s = new Sequence(10000, Parameters::DBTYPE_NUCLEOTIDES, &subMat, kmer_size, false, false);
    s->mapSequence(0, 0, sequence, strlen(sequence));

    size_t i = 0;
    while (s->hasNextKmer()) {
        const unsigned char* curr_pos = s->nextKmer();
        printf("Pos1: %zu\n", i++);

        size_t idx_val = idx.int2index(curr_pos);
        std::cout << "Index: " << idx_val << " ";
        idx.printKmer(idx_val, kmer_size, subMat.num2aa);
        std::cout << std::endl;

        std::string kmerStr1;
        uint64_t kmerIdx = 0;
        for(size_t kmerPos = 0; kmerPos < kmer_size; kmerPos++){
            kmerIdx = kmerIdx << 2;
            kmerIdx = kmerIdx | curr_pos[kmerPos];
            kmerStr1.append(1, subMat.num2aa[curr_pos[kmerPos]]);
        }

        std::string revStr1;
        for(int kmerPos = kmer_size - 1; kmerPos >= 0; kmerPos--){
            revStr1.append(1, Orf::complement(kmerStr1[kmerPos]));
        }
        std::string kmerStr2 = kmerToSting(kmerIdx, kmer_size);
        size_t revComp = Util::revComplement(kmerIdx, kmer_size);
        std::string revStr2 = kmerToSting(revComp, kmer_size);

        std::cout << kmerStr1 << " " << kmerStr2 << "\n";
        std::cout << revStr1 << " " << revStr2 << "\n";
        if (kmerStr1 != kmerStr2 || revStr1 != revStr2) {
            std::cout << "fail";
            return EXIT_FAILURE;
        }
    }
    size_t expectedKmers = strlen(sequence) - kmer_size;
    if (i - 1 != expectedKmers) {
        std::cout << "Invalid number of k-mers. Expected " << expectedKmers << " found " << (i - 1) << ".\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

