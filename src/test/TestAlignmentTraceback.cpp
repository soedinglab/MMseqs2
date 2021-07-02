//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include "Sequence.h"
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"

#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "BaseMatrix.h"
#include "StripedSmithWaterman.h"
#include "Util.h"
#include "Parameters.h"

const char* binary_name = "test_alignmenttraceback";

struct scores{
    short H;
    short F;
    short E;
};

scores workspace[10000*2 + 2];
unsigned char bt[10000*10000];

void sw(
        const unsigned char *db_sequence, const unsigned char *query_sequence,
        const short ** profile_word,
        int32_t query_start, int32_t query_end,
        int32_t target_start, int32_t target_end,
        const short gap_open, const short gap_extend, SubstitutionMatrix &subMat){
    const int M = 2;
    const int F = 1;
    const int E = 0;

    scores * curr_sM_G_D_vec = &workspace[0];
    int query_length = (query_end - query_start);
    int target_length = (target_end - target_start);
    scores * prev_sM_G_D_vec = &workspace[query_length + 1];
    memset(prev_sM_G_D_vec, 0, sizeof(scores) * (query_length + 1));
    short goe = gap_open + gap_extend;
    int pos = 0;
    for (int i = target_start; LIKELY(i < target_end); i++) {
        prev_sM_G_D_vec[query_start-1].H = 0;
        curr_sM_G_D_vec[query_start-1].H = 0;
        curr_sM_G_D_vec[query_start-1].E = 0;
        const short * profile = profile_word[db_sequence[i]];
        //TODO problem if query_start == 0
        for (int j = query_start; LIKELY(j < query_end); j++) {
            curr_sM_G_D_vec[j].E = std::max(curr_sM_G_D_vec[j-1].H - goe, curr_sM_G_D_vec[j-1].E - gap_extend); // j-1
            curr_sM_G_D_vec[j].F = std::max(prev_sM_G_D_vec[j].H   - goe, prev_sM_G_D_vec[j].F - gap_extend);   // i-1
            const short H = prev_sM_G_D_vec[j-1].H + profile[j]; // i - 1, j - 1
            curr_sM_G_D_vec[j].H = std::max(H, curr_sM_G_D_vec[j].E);
            curr_sM_G_D_vec[j].H = std::max(curr_sM_G_D_vec[j].H, curr_sM_G_D_vec[j].F);
            const unsigned char mode1 = (curr_sM_G_D_vec[j].H == H) ? M : E;
            const unsigned char mode2 = (curr_sM_G_D_vec[j].H == curr_sM_G_D_vec[j].F) ? F : E;
            const unsigned char mode = std::max(mode1, mode2);
            bt[pos/4] |= mode << (pos % 4) * 2;
            pos++;
        }
        std::cout << std::endl;

        // swap rows
        scores * tmpPtr = prev_sM_G_D_vec;
        prev_sM_G_D_vec = curr_sM_G_D_vec;
        curr_sM_G_D_vec = tmpPtr;
    }

    // 0x03 00000011
#define get_val(bt, i, j) ( bt[(i * query_length + j)/4] >> (((i * query_length + j) % 4) * 2)  & 0x03 )
//    std::cout << std::endl;
//    std::cout << std::endl;
//  PRINT BT matrix
//    for (int i = 0; LIKELY(i < target_length); i++) {
//        std::cout << i << ": ";
//        for (int j = 0; LIKELY(j < query_length); j++) {
//            std::cout << get_val(bt, i, j) << " ";
//        }
//        std::cout << std::endl;
//    }

    // backtrace
    int i=target_length - 1;
    int j=query_length - 1;
    int step = 0;
    int state = get_val(bt, i, j);
    int matched_cols = 0;
    while (state!=-1)     // while (state!=STOP)  because STOP=0
    {
        step++;
        //std::cout << step<< " " << i << " " << j << " " << state << std::endl;
        switch (state) {
            case M: // current state is MM, previous state is bMM[i][j]
                matched_cols++;
                fprintf(stdout,"%c", subMat.num2aa[db_sequence[target_start+i]]);
                if(query_sequence[query_start + j] == db_sequence[target_start + i]){
                    fprintf(stdout, "|");
                }else{
                    fprintf(stdout, "*");
                }
                fprintf(stdout,"%c", subMat.num2aa[query_sequence[query_start+j]]);
                i--; j--;
                state = (i < 0 || j < 0) ? -1 : get_val(bt, i, j);
                break;
            case E: // current state is GD
                fprintf(stdout, "|  ");

                if (j < 0){
                    state = -1;
                }else{
                    j--;
                    state =  get_val(bt, i, j);
                }
                break;
            case F:
                fprintf(stdout, "  |");

                if (i < 0){
                    state = -1;
                }else{
                    i--;
                    state = get_val(bt, i, j);
                }
                break;
        }
        std::cout << std::endl;
    }

}

int main(int, const char**) {
    const size_t kmer_size=6;

    Parameters& par = Parameters::getInstance();
    par.initMatrices();
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
    std::cout << "Substitution matrix:\n";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.num2aa,subMat.alphabetSize);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    //    ReducedMatrix subMat(subMat.probMatrix, 20);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";
    //static const char ref_seq[40] = {'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T', 'C', 'T', 'G', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'T',
    //						'C', 'A', 'A', 'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'A', 'A', 'A', '\0'};
    //static const char read_seq[16] = {'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G', 'G', 'T', 'A', 'A', 'A', 'T', 'C', '\0'};	// read sequence
//	std::string tim = "APRKFFVGGNWKMNGKRKSLGELIHTLDGAKLSADTEVVCGAPSIYLDFARQKLDAKIGVAAQNCYKVPKGAFTGEISPAMIKDIGAAWVILGH"
//                      "SERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFQETKAIADNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLR"
//			          "GWLKTHVSDAVAVQSRIIYGGSVTGGNCKELASQHDVDGFLVGGASLKPEFVDIINAKH";
    std::string tim1 = "MDDVKIERLKRLNEDVLEDLIEVYMRGYEGLEEYGGEGRDYARDYIKWCWKKAPDGFFVAKVGDRIVGFIVCDRDWYSRYEGKIVGAIHEFVVDKGWQGKGIGKKLLTKCLEFLGKYNDTIELWVGEKNFGAMRLYEKFGFKKVGKSGIWIRMVRRQLS";
    std::string tim2 = "LRSKETFNDMNLPSRHAIAKVVSIEQQLYDNLAYPELLFYQAAHQWPNSQFICRDNNDILAYAMYAPAEKANTLWLMSAAVKPGCQGRGVGTKLLSDSLRSLDEQGVTCVLLSVAPSNAAAISVYQKLGFEVVRKAEHYLKNLREQGLRMTREIIHK";
    std::cout << "Sequence (id 0):\n";
    //const char* sequence = read_seq;
    std::cout << tim1 << "\n\n";
    Sequence* s = new Sequence(10000, 0, &subMat, kmer_size, true, false);
    s->mapSequence(0,0,tim1.c_str(), tim1.size());
    Sequence* dbSeq = new Sequence(10000, 0, &subMat, kmer_size, true, false);
    //dbSeq->mapSequence(1,"lala2",ref_seq);
    dbSeq->mapSequence(1,1,tim2.c_str(), tim2.size());
    SmithWaterman aligner(15000, subMat.alphabetSize, false);
    int8_t * tinySubMat = new int8_t[subMat.alphabetSize*subMat.alphabetSize];
    for (int i = 0; i < subMat.alphabetSize; i++) {
        for (int j = 0; j < subMat.alphabetSize; j++) {
            std::cout << ( i*subMat.alphabetSize + j) << " " << subMat.subMatrix[i][j] << " ";

            tinySubMat[i*subMat.alphabetSize + j] = (int8_t)subMat.subMatrix[i][j];
        }
        std::cout << std::endl;
    }
    short ** profile = new short*[subMat.alphabetSize];
    short * profile_data = new short[10000*subMat.alphabetSize];

    for(int i = 0; i< subMat.alphabetSize; i++) {
        profile[i] = &profile_data[i*s->L];
        for (int j = 0; j < s->L; j++) {
            profile[i][j] = tinySubMat[i*subMat.alphabetSize + s->numSequence[j]];
        }
    }

    sw(dbSeq->numSequence, s->numSequence, (const short ** ) profile, 92, 157, 80, 146, 11, 1, subMat);
    // calcuate stop score

    delete [] tinySubMat;

    delete s;
    delete dbSeq;
    return 0;
}

