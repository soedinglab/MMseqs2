#include <iostream>
#include <cassert>
#include "Util.h"
#include "TaxonomyExpression.h"

const char* binary_name = "test_taxexpr";



int main (int, const char**) {
//    TaxonomyExpression expression("(2,2157),4751,33208,33090,(2759,!4751,!33208,!33090)");
    TaxonomyExpression expression1("2");
    std::string path = "/Users/mad/Documents/databases//swissprot/sprot_new";
    NcbiTaxonomy * taxonomy = NcbiTaxonomy::openTaxonomy(path);

    if(expression1.isAncestorOf(*taxonomy, 1117) != -1){
        std::cout << "Found bacteria" << std::endl;
    } else{
        assert(false);
    }

    if(expression1.isAncestorOf(*taxonomy, 33630) == -1){
        std::cout << "Alveolata not ancestor" << std::endl;
    } else{
        assert(false);
    }

    TaxonomyExpression expression2("(2759&!9606)");

    if(expression2.isAncestorOf(*taxonomy, 33630) != -1){
        std::cout << "Found Alveolata" << std::endl;
    } else{
        assert(false);
    }

    if(expression2.isAncestorOf(*taxonomy, 9606) == -1){
        std::cout << "Homo sapiens not ancestor" << std::endl;
    } else{
        assert(false);
    }

    TaxonomyExpression expression3("(2759&!61964),10239");

    if(expression3.isAncestorOf(*taxonomy, 114777) == 1){
        std::cout << "Found Natrialba phage PhiCh1" << std::endl;
    } else{
        assert(false);
    }

    if(expression3.isAncestorOf(*taxonomy, 2759) == 0){
        std::cout << "Found Eukaryota" << std::endl;
    } else{
        assert(false);
    }

    if(expression3.isAncestorOf(*taxonomy, 61964) == -1){
        std::cout << "Enviromental sample in not in" << std::endl;
    } else{
        assert(false);
    }

    TaxonomyExpression expression4("2759,10239");

    if(expression4.isAncestorOf(*taxonomy, 114777) == 1){
        std::cout << "Found Natrialba phage PhiCh1" << std::endl;
    } else{
        assert(false);
    }

    if(expression4.isAncestorOf(*taxonomy, 2759) == 0){
        std::cout << "Found Eukaryota" << std::endl;
    } else{
        assert(false);
    }

    if(expression4.isAncestorOf(*taxonomy, 61964) == 0){
        std::cout << "Found Enviromental sample" << std::endl;
    } else{
        assert(false);
    }


    TaxonomyExpression expression5("!2759");

    if(expression5.isAncestorOf(*taxonomy, 2) == 0){
        std::cout << "Found Bacteria" << std::endl;
    } else{
        assert(false);
    }

    if(expression5.isAncestorOf(*taxonomy, 2759) == -1){
        std::cout << "Eukaryota not in" << std::endl;
    } else{
        assert(false);
    }


    TaxonomyExpression expression6("(2|2759)");

    if(expression6.isAncestorOf(*taxonomy, 2) == 0){
        std::cout << "Found Bacteria" << std::endl;
    } else{
        assert(false);
    }

    if(expression6.isAncestorOf(*taxonomy, 2759) == 0){
        std::cout << "Found Eukaryota" << std::endl;
    } else{
        assert(false);
    }

    if(expression6.isAncestorOf(*taxonomy, 10239) == -1){
        std::cout << "Virus sample not in" << std::endl;
    } else{
        assert(false);
    }

    delete taxonomy;
}

