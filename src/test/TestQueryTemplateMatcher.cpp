
#include <iostream>
#include "Sequence.h"
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "QueryTemplateMatcher.h"
int main (int argc, const char * argv[])
{
    
    const size_t kmer_size=7;

    SubstitutionMatrix subMat("/cluster/user/maria/kClust2/data/blosum62.out");
    
    QueryTemplateMatcher queryTemplateMatcher();
    
    return 0;
}

