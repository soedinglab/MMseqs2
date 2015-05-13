//
// Created by lars on 12.04.15.
//

#include <Util.h>
#include "CompareGOTerms.h"

void printHelp();

int main(int argc, char **argv)
{

    if (argc < 1) {
        printHelp();

    }
   // for (int k = 0; k < 8; ++k) {
    //    Debug(Debug::INFO) << argv[k] <<"\n";
    //}
    if(strcmp(argv[0],"-go")){
        Debug(Debug::INFO) <<"GO-Evaluation" <<"\n";

        if(argc != 7){
            Debug(Debug::INFO) << argc << "\n";
            printHelp();

        }

        std::string gofolder=argv[2];
        std::string uniprot_go_folder=argv[3];
        std::string clustering_file=argv[4];
        std::string prefix=argv[5];
        std::string outputfolder=argv[6];

        //"-go <gofolder> <prot_go_folder> <clustering_file> <prefix> <outputfolder>"
        //std::string gofolder="/home/lars/masterarbeit/data/GO/db/";
        //std::string uniprot_go_folder="/home/lars/masterarbeit/data/uniprot/release-2015_04/uniprot_go/";
        //"/home/lars/masterarbeit/data/uniprot/release-2015_04/evaluation"
          //      "/home/lars/masterarbeit/db/sprot/uniprot_sprot_s4_affinity"

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
                                                        outputfolder);
                go->init();
                //go->all_against_all_comparison();
                //  go->all_against_all_comparison_proteinset();
                go->run_evaluation_mmseqsclustering(clustering_file,
                                                    clustering_file+".index",
                                                    prefix, evidenceCategories[j] + goCategories[i]);
                go->~CompareGOTerms();
            }
        }
    }




    }

void printHelp() {
    std::string usage("\nEvaluation commands\n");
    usage.append("-go <gofolder> <prot_go_folder> <clustering_file> <prefix> <outputfolder>");
    Debug(Debug::INFO) << usage << "\n";
    EXIT(EXIT_FAILURE);
}

