//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <smith_waterman_sse2.h>
#include "Sequence.h"
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "BaseMatrix.h"
#include "../alignment/smith_waterman_sse2.h"
int main (int argc, const char * argv[])
{

    const size_t kmer_size=6;

    SubstitutionMatrix subMat("/Users/mad/Documents/workspace/mmseqs/data/blosum62.out", 2.0, 0);
    std::cout << "Subustitution matrix:\n";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.int2aa,subMat.alphabetSize);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    //    ReducedMatrix subMat(subMat.probMatrix, 20);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";
    //static const char ref_seq[40] = {'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T', 'C', 'T', 'G', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'T',
    //						'C', 'A', 'A', 'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'A', 'A', 'A', '\0'};
    //static const char read_seq[16] = {'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G', 'G', 'T', 'A', 'A', 'A', 'T', 'C', '\0'};	// read sequence
    std::string profile = "}|z{|zz\x80{~~z{{{{}\x86z|}{\x84\x82z{||}{|}\x85}|\x82}\x81z{{z\x87~z||z}zz\x83|}|}|zyz}|z{}zz\x7f{\x82~z{{{{}\x86{||z\x85\x86z{}z}z{~|~|}|zz{\x83{}\x84{|||~\x80|||~\x84}||{{\x86}{|{}{||||||||~}}{{}|z{}zz\x82{\x81~z{{{{}\x85{||z\x82\x84z|}z~z{\x83|~\x84\x82}{z{{{zz\x88z|}z~}zzzz{{|~\x80{{z{\x80z\x7f|{||{z|{{||\x7f\x89}|||{{||}||}|}\x83~\x87}{{|z\x85~z{|z\x82z{~{~\x85}|zz{\x82{|}{|\x88{}{|}{\x84}}|{{~||zz~yz\x7f{\x86\x7fzz{{{|~{|}{|{z\x87{y|yz}{{{}|z{z{{zz\x87z|}z\x82~zzzz{{}~\x7f\x81{\x83\x82{||}}||}|\x83|}\x82\x82z{||z{}z{~}\x84~{{|\x84{|\x82{|\x81{{}{{||\x82|||{~\x86}|\x81z{|{{}\x84{|}}\x81\x86|{\x85|||}|~}{|~z\x83\x84z\x82z|}|\x85}\x82}{{|\x84|{|\x84|{||||{\x85|{}|||~\x85|||{\x85{{|{{|||{~||{{}z||z\x83{z|z{|\x88|{||{zz\x84|{|{\x82{}||||||{\x82}\x83z{|{|\x83\x86z|}|~\x86{{}|||}}~\x86}{||}{}|\x82}{|||~}}{|||\x82{}z{\x84{\x82\x84{{{{{|\x7f{|||z{}\x83z~{\x84~{{{{{|\x81{|}{|}z\x84|z\x82{|}{\x83\x83}\x82{{{|z|}z{}z\x7f{||\x85~\x87||{z{|\x86\x85\x82z\x81|{|z{}{}{}|{zz|z\x84}z\x83|z}z{}\x85|\x82}|zzz||zz\x83z{\x7f{\x84~zz{{{|\x84|}|z\x84\x82z{}{~{{}{~\x85}\x82{z{||zz\x83z{~{\x85\x7fzz{{{|~|}{|yz~zz~{\x85~zz{{{|}\x88}||zz}zz\x82{\x84\x7fzz{{{|\x84{|\x83|}}{}|{\x83{|\x83|}}\x84}|z{\x82{|}z\x85|z}z{}{\x84|}|{{{|z||z\x81{z|z{{\x88|{||{zz\x81{{{{\x84{||}\x86|\x83|{}||{{\x83|||{\x84{{|{|}|||\x82\x84|{{}{||z\x84|z}z{}|}\x84\x85}{z{\x86}{|{}{||||||||~}}{{\x82|}|{\x84|{}{|}|}|\x84~|{{|z}~z|}z\x82{|\x84|\x85\x84}|{z{\x81{||z{{{}{{|\x87|\x82}|{z{|{{{}\x83{~|\x81\x88{{|{||}{|\x82|||z||{|{{}\x86||\x85}|z{\x82{\x86~z\x82|z|z{~||{}|{zz}{|}z\x86|z}z{}{\x84|}|z{{\x81{|\x82{{|}~||||}\x84}\x82\x82z{}|{||{{~}\x81}|{|\x84}\x84\x83{||z\x84}z{{z|z{|\x88|{||{yz|{|\x82\x84{}|\x82\x81}|{\x85}||||}\x83{||z|{{|{{{\x88|{}|{z{}z||z\x86|z}z{}{|\x83}|z{z}{||{\x86\x85z|z{}{|{}|z{|{{yz~{{{z{|yy{zz{z\x8c\x7f|{~}z|}{}z{\x87{}}~\x84{z{|z|}z{}z\x7f{||\x85~\x87||{z{|\x85zz\x83zz\x83z~~zzzz{|\x84{}|{|}|{\x87}}|}}{\x85}||\x83{}||zz\x83zz\x83{\x83~zz{{{|\x82|}|{|\x84|z|\x82}}}|{}\x83|}\x84z|}|z{|zz\x80{~~z{{{{}\x86z|\x83{\x84~z||{}{{}\x85\x83|}|{z{{z\x87\x7fz||z|zz~|}|}|zyz||zz}yz\x82{\x85\x7fzz{{{|~{|\x83{\x84\x84z||{}{{}\x83}|}|{z{\x86|{|{}||}||||}\x82~}}{{\x83||\x84|{|~|\x82}{{}|}}\x82{|\x81|z{|{z\x7f{\x82~z{{{|}\x85{|\x83{\x84}z||{|{{}\x85||~\x82|z{\x82z|\x82z|}z\x82{|}|~\x86}|{z{\x82|z{}{{~{\x84\x84z{|{||~{||z{}{{||\x7f\x82}|{~\x87|||{|\x84{}\x82z|}{~{|||\x85\x83}}|{{\x86|||{}|||{|\x83|||~}}z{}z|\x82z\x86{y|z{}{|{}|z{z\x84|z{||{\x82{\x82~{{{{}}\x82{|}|||{{{||\x81}|\x85||\x82\x85}z{{{{|\x87z|}|}}{z\x84|{{|~\x7f|z|\x83{{}{\x7f\x81}|{~\x86||{{{|{~|z\x84}z|z{\x86{}|}}zz{|z\x84\x86z{}z~z{}|~}}|zz{||zz}zz\x85{~~z{{z{}\x84z|}|{||z{\x82|~~{{\x84||}\x85{|}{\x82\x83z||{}{|}|\x83}\x82\x85|z{}{||{\x86\x85z|z{}{|{}|z{|\x83{{|{|{\x81}||{\x86|\x82}|}z{}{{{|\x86{z{z{|{{{|{z\x88|}{|{z\x87{y|yz}{{{}|z{z|z|~z{}z\x84{|}|\x83\x86}|{z{|{|~{{}|~\x82}}|\x88~}||{|\x82|{{|{{\x85|~}|{|{\x83}\x7f{|||zz~yz\x7f{\x86\x7fzz{{{|~{|\x82\x86z{}{z~{\x84~z{{{||}{|\x83|}\x84{{|}}|}||}\x83}}\x83z{{z\x87\x7fz||z|zz~|}|}|zyz|z||z{{z|z{{\x89|{||{yz\x83|||\x84}|||||}|||\x85~||}}{|{z\x87{y|yz}{{{}|z{z|{\x7f}z}~z}z{\x88{}}~}{z{}|{{}z{\x7f{\x82~{\x84{{|}\x85{|}|z{}zz\x84{\x7f~z{{{{}\x85{|}z\x7f\x87z{}z~{{}|\x7f}}}{{{";

    std::cout << "Sequence (id 0):\n";
    //const char* sequence = read_seq;
    const char* sequence = profile.c_str();
    std::cout << sequence << "\n\n";
    Sequence* s = new Sequence(10000, subMat.aa2int, subMat.int2aa, Sequence::HMM_PROFILE, kmer_size, true);
    s->mapSequence(0,"lala",sequence);
    s->printProfile();
    Sequence* dbSeq = new Sequence(10000, subMat.aa2int, subMat.int2aa, Sequence::AMINO_ACIDS, kmer_size, true);
    //dbSeq->mapSequence(1,"lala2",ref_seq);
    const char* sequence2 = "IIRLNHVAVATLQLEKLTSFYRDTLGLQVSEPVPQKEHGVTTVFVDVGNTKFELLLPLGDKSPIANFLEKNKGGGAHHVCLEVDDIEAAVADLKXXGIRMLAEKTRIGAHGKPVMFLHPKDCGGVLVELEQ\n";

    dbSeq->mapSequence(1,"lala2",sequence2);
    SmithWaterman aligner(15000, subMat.alphabetSize);
    int8_t * tinySubMat = new int8_t[subMat.alphabetSize*subMat.alphabetSize];

    aligner.ssw_init(s, s->getAlignmentProfile(), &subMat, subMat.alphabetSize, 2);
    int32_t maskLen = s->L / 2;
    int gap_open = 10;
    int gap_extend = 1;
    s_align * alignment = aligner.ssw_align(dbSeq->int_sequence, dbSeq->L, gap_open, gap_extend, 2, 55, 0, maskLen);
    if(alignment->cigar){
        std::cout << "Cigar" << std::endl;

        int32_t targetPos = alignment->dbStartPos1, queryPos = alignment->qStartPos1;
        for (int32_t c = 0; c < alignment->cigarLen; ++c) {
            char letter = SmithWaterman::cigar_int_to_op(alignment->cigar[c]);
            uint32_t length = SmithWaterman::cigar_int_to_len(alignment->cigar[c]);
            for (uint32_t i = 0; i < length; ++i){
                if (letter == 'M') {
                    fprintf(stdout,"%c",subMat.int2aa[dbSeq->int_sequence[targetPos]]);
                    if (dbSeq->int_sequence[targetPos] == s->int_sequence[queryPos]){
                        fprintf(stdout, "|");
                    }
                    else fprintf(stdout, "*");
                    fprintf(stdout,"%c",subMat.int2aa[s->int_sequence[queryPos]]);
                    ++queryPos;
                    ++targetPos;
                } else {
                    if (letter == 'I'){
                        fprintf(stdout,"%c |",subMat.int2aa[s->int_sequence[queryPos]]);
                        ++queryPos;
                    } else{
                        fprintf(stdout,"| %c",subMat.int2aa[dbSeq->int_sequence[targetPos]]);
                        ++targetPos;
                    }
                }
                std::cout << std::endl;
            }
        }
    }
    std::cout << alignment->score1 << " "<< alignment->qEndPos1 << " " << alignment->qStartPos1  << " "<< alignment->qEndPos1 << " "
            << alignment->dbStartPos1 << " "<< alignment->dbEndPos1 << std::endl;
    delete [] tinySubMat;
    delete [] alignment->cigar;
    delete alignment;
    delete s;
    delete dbSeq;
    return 0;
}

