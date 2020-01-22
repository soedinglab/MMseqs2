#include <iostream>
#include <cassert>
#include <tinyexpr.h>
#include <ExpressionParser.h>
#include "Util.h"
#include "TaxonomyExpression.h"

const char* binary_name = "test_taxexpr";



int main (int, const char**) {
    std::string path = "/Users/mad/Documents/databases//swissprot/sprot_new";
    NcbiTaxonomy *taxonomy = NcbiTaxonomy::openTaxonomy(path);
    TaxonomyExpression parser("!2", *taxonomy);
    if(parser.isAncestor(9606) == true){
        std::cout << "Found human" << std::endl;
    } else{
        assert(false);
    }

    if(parser.isAncestor(1117) == false){
        std::cout << "Alveolata not ancestor" << std::endl;
    } else{
        assert(false);
    }

    TaxonomyExpression expression2("(2759&&!9606)",*taxonomy);

    if(expression2.isAncestor(33630)){
        std::cout << "Found Alveolata" << std::endl;
    } else{
        assert(false);
    }

    if(expression2.isAncestor( 9606) == false){
        std::cout << "Homo sapiens not ancestor" << std::endl;
    } else{
        assert(false);
    }

    TaxonomyExpression expression3("2759||10239",*taxonomy);

    if(expression3.isAncestor(114777)){
        std::cout << "Found Natrialba phage PhiCh1" << std::endl;
    } else{
        assert(false);
    }

    if(expression3.isAncestor(2759)){
        std::cout << "Found Eukaryota" << std::endl;
    } else{
        assert(false);
    }

    if(expression3.isAncestor(61964)){
        std::cout << "Found Enviromental sample" << std::endl;
    } else{
        assert(false);
    }


    TaxonomyExpression expression5("!2759",*taxonomy);

    if(expression5.isAncestor(2)){
        std::cout << "Found Bacteria" << std::endl;
    } else{
        assert(false);
    }

    if(expression5.isAncestor(2759) == false){
        std::cout << "Eukaryota not in" << std::endl;
    } else{
        assert(false);
    }


    TaxonomyExpression expression6("(2||2759)",*taxonomy);

    if(expression6.isAncestor(2)){
        std::cout << "Found Bacteria" << std::endl;
    } else{
        assert(false);
    }

    if(expression6.isAncestor(2759)){
        std::cout << "Found Eukaryota" << std::endl;
    } else{
        assert(false);
    }

    if(expression6.isAncestor( 10239) == false){
        std::cout << "Virus sample not in" << std::endl;
    } else{
        assert(false);
    }

    TaxonomyExpression expression7("(2&&!1117)",*taxonomy);

    if(expression7.isAncestor(57723)){
        std::cout << "Found Acidobacteria" << std::endl;
    } else{
        assert(false);
    }

    if(expression7.isAncestor(1117) == false){
        std::cout << "Cyanobacteria not in" << std::endl;
    } else{
        assert(false);
    }

    if(expression7.isAncestor(9606) == false){
        std::cout << "Human not in" << std::endl;
    } else{
        assert(false);
    }

    delete taxonomy;
}

