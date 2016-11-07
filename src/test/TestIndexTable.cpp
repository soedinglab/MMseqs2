//
// Written by Maria Hauser mhauser@genzentrum.lmu.de
//
// Test class for k-mer generation and index table testing.
//


#include <cstdio>
#include <iostream>
#include <SubstitutionMatrix.h>
#include <IndexTable.h>
#include <Prefiltering.h>



int main (int argc, const char * argv[])
{

    SubstitutionMatrix subMat("blosum62.out", 8.0, -0.2);
    DBReader<unsigned int> dbr(
                               "/Users/mad/Documents/databases/db_small/db_small",
//                               "/Users/mad/Documents/databases/mmseqs_benchmark/benchmarks/protein_search_uniscop/db/mmseqs/db_sw",
                               "/Users/mad/Documents/databases/db_small/db_small.index"
//                               "/Users/mad/Documents/databases/mmseqs_benchmark/benchmarks/protein_search_uniscop/db/mmseqs/db_sw.index"
                               );
    dbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    Sequence* s = new Sequence(32000, subMat.aa2int, subMat.int2aa, Sequence::AMINO_ACIDS, 6, true, false);
    IndexTable t(subMat.alphabetSize, 6);
    Prefiltering::fillDatabase(&dbr, s, &t, &subMat, 0, dbr.getSize(), false, 1 );

    delete s;
    return 0;
}

