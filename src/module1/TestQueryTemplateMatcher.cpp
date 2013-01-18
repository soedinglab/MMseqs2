//
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
#include "QueryTemplateMatcher.h"
int main (int argc, const char * argv[])
{
    
    const size_t kmer_size=7;
    
    
    SubstitutionMatrix subMat("/cluster/user/maria/kClust2/data/blosum30.out",20);
    
    QueryTemplateMatcher queryTemplateMatcher();
    
    return 0;
}

