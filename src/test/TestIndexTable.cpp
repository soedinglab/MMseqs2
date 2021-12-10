//
// Written by Maria Hauser mhauser@genzentrum.lmu.de
//
// Test class for k-mer generation and index table testing.
//

#include <cstdio>
#include <iostream>
#include "SubstitutionMatrix.h"
#include "IndexTable.h"
#include "IndexBuilder.h"
#include "Parameters.h"

const char* binary_name = "test_indextable";

int main (int, const char**) {
    Parameters &par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2f);
    DBReader<unsigned int> dbr(
                               "/Users/mad/Documents/databases/db_small/db_small",
//                               "/Users/mad/Documents/databases/mmseqs_benchmark/benchmarks/protein_search_uniscop/db/mmseqs/db_sw",
                               "/Users/mad/Documents/databases/db_small/db_small.index"
//                               "/Users/mad/Documents/databases/mmseqs_benchmark/benchmarks/protein_search_uniscop/db/mmseqs/db_sw.index"
                               , 1, 0);
    dbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Sequence *s = new Sequence(32000, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 6, true, false);
    IndexTable t(subMat.alphabetSize, 6, false);
    IndexBuilder::fillDatabase(&t, NULL, NULL, subMat, s, &dbr, 0, dbr.getSize(), 0, 1, 1,0.9);
    t.printStatistics(subMat.num2aa);

    delete s;
    dbr.close();

    return 0;
}

