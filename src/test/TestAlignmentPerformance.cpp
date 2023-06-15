//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>

#include "kseq.h"
#include "Util.h"
#include "Parameters.h"
#include "Sequence.h"
#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrix.h"
#include "StripedSmithWaterman.h"

const char* binary_name = "test_alignmentperformance";

#define MAX_FILENAME_LIST_FILES 4096

KSEQ_INIT(int, read)


std::vector<std::string> readData(std::string fasta_filename){
    std::vector<std::string> retVec;
    kseq_t *seq;
    FILE* fasta_file = fopen(fasta_filename.c_str(), "r");
    if(fasta_file == NULL) {std::cout << "Could not open " << fasta_filename<<std::endl; EXIT(1); }
    seq = kseq_init(fileno(fasta_file));
    size_t entries_num = 0;
    while (kseq_read(seq) >= 0) {
        if (entries_num > 1000)
            break;
        if (seq->seq.l > 500) {
            std::string sequence = seq->seq.s;
            retVec.push_back(sequence);
            entries_num++;
        }
    }
	kseq_destroy(seq);
    fclose(fasta_file);
    return retVec;
}
int main (int, const char**) {
    const size_t kmer_size=6;

    Parameters& par = Parameters::getInstance();
    par.initMatrices();
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0);
    std::cout << "Substitution matrix:\n";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.num2aa,subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    std::cout << "\n";


    std::cout << "Sequence (id 0):\n";
    Sequence* query = new Sequence(10000, 0, &subMat, kmer_size, true, false);
    Sequence* dbSeq = new Sequence(10000, 0, &subMat, kmer_size, true, false);
    //dbSeq->mapSequence(1,"lala2",ref_seq);
    SmithWaterman aligner(15000, subMat.alphabetSize, false, 1.0, Parameters::DBTYPE_AMINO_ACIDS);
    int8_t * tinySubMat = new int8_t[subMat.alphabetSize*subMat.alphabetSize];
    for (int i = 0; i < subMat.alphabetSize; i++) {
        for (int j = 0; j < subMat.alphabetSize; j++) {
            std::cout << ( i*subMat.alphabetSize + j) << " " << subMat.subMatrix[i][j] << " ";

            tinySubMat[i*subMat.alphabetSize + j] = (int8_t)subMat.subMatrix[i][j];
        }
        std::cout << std::endl;
    }

    int gap_open = 10;
    int gap_extend = 1;
    int mode = 0;
    size_t cells = 0;
    std::vector<std::string> sequences = readData("/Users/mad/Documents/databases/rfam/Rfam.fasta");
    for(size_t seq_i = 0; seq_i < sequences.size(); seq_i++){
        query->mapSequence(1,1,sequences[seq_i].c_str(), sequences[seq_i].size());
        aligner.ssw_init(query, tinySubMat, &subMat);

        for(size_t seq_j = 0; seq_j < sequences.size(); seq_j++) {
            dbSeq->mapSequence(2, 2, sequences[seq_j].c_str(),  sequences[seq_j].size());
            int32_t maskLen = query->L / 2;
            EvalueComputation evalueComputation(100000, &subMat, gap_open, gap_extend);
            std::string backtrace;
            s_align alignment = aligner.ssw_align(
                    dbSeq->numSequence,
                    dbSeq->numConsensusSequence,
                    dbSeq->getAlignmentProfile(),
                    dbSeq->L,
                    backtrace,
                    gap_open, gap_extend,
                    0,
                    10000,
                    &evalueComputation,
                    0, 0.0,
                    0.0,
                    maskLen,
                    dbSeq->getId()
            );
            if(mode == 0 ){
                cells += query->L * dbSeq->L;
                std::cout << alignment.qEndPos1 << " " << alignment.dbEndPos1 << "\n";
            } else {
                std::cout << alignment.qStartPos1 << " " << alignment.qEndPos1 << " "
                        << alignment.dbStartPos1 << " " << alignment.dbEndPos1 << "\n";
            }
        }
    }
    std::cerr << "Cells : " << cells << std::endl;
    delete [] tinySubMat;
    delete query;
    delete dbSeq;
    return 0;
}

