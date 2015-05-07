//
// Created by lars on 12.04.15.
//

#include "CompareGOTerms.h"

int main(int argc, char **argv)
{

    CompareGOTerms* go=new CompareGOTerms("/home/lars/masterarbeit/data/GO/go-basic_ffindex","/home/lars/masterarbeit/data/GO/go-basic_ffindex.index",
            "/home/lars/masterarbeit/data/uniprot/release-2015_04/uniprot_sprot_go_db", "/home/lars/masterarbeit/data/uniprot/release-2015_04/uniprot_sprot_go_db.index","/home/lars/masterarbeit/data/uniprot/release-2015_04/evaluation");

    go->init();

    //go->all_against_all_comparison();
  //  go->all_against_all_comparison_proteinset();
    go->run_evaluation_mmseqsclustering("/home/lars/masterarbeit/db/sprot/uniprot_sprot_s4_affinity","/home/lars/masterarbeit/db/sprot/uniprot_sprot_s4_affinity.index");
}

