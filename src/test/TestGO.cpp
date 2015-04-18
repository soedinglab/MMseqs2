//
// Created by lars on 12.04.15.
//

#include "CompareGOTerms.h"

int main(int argc, char **argv)
{

    CompareGOTerms* go=new CompareGOTerms("/home/lars/masterarbeit/data/GO/go-basic_ffindex","/home/lars/masterarbeit/data/GO/go-basic_ffindex.index",
            "/home/lars/masterarbeit/data/uniprot/release-2015_04/uniprot_sprot_go_db", "/home/lars/masterarbeit/data/uniprot/release-2015_04/uniprot_sprot_go_db.index");

    go->init();
}

