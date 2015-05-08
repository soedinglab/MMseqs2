//
// Created by lars on 12.04.15.
//

#include "CompareGOTerms.h"

int main(int argc, char **argv)
{
    std::string gofolder="/home/lars/masterarbeit/data/GO/db/";
    std::string uniprot_go_folder="/home/lars/masterarbeit/data/uniprot/release-2015_04/uniprot_go/";


    std::string* goCategories= new std::string[3];
    goCategories[0]="_C";
    goCategories[1]="_F";
    goCategories[2]="_P";

    std::string* evidenceCategories= new std::string[2];
    evidenceCategories[0]="";
    evidenceCategories[1]="_EXP";

    for (int j = 0; j <2 ; ++j) {


        for (int i = 0; i < 3; ++i) {
            CompareGOTerms *go = new CompareGOTerms(gofolder + "go-fasta_db" + goCategories[i],
                                                    gofolder + "go-fasta_db" + goCategories[i] + ".index",
                                                    uniprot_go_folder + "uniprot_sprot.dat_go_db" +
                                                    evidenceCategories[j] + goCategories[i],
                                                    uniprot_go_folder + "uniprot_sprot.dat_go_db" +
                                                    evidenceCategories[j] + goCategories[i] + ".index",
                                                    "/home/lars/masterarbeit/data/uniprot/release-2015_04/evaluation");

            go->init();

            //go->all_against_all_comparison();
            //  go->all_against_all_comparison_proteinset();
            go->run_evaluation_mmseqsclustering("/home/lars/masterarbeit/db/sprot/uniprot_sprot_s4_affinity",
                                                "/home/lars/masterarbeit/db/sprot/uniprot_sprot_s4_affinity.index",
                                                "affinity_", evidenceCategories[j] + goCategories[i]);
        }
    }
    }

