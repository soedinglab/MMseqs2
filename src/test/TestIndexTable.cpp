//
// Written by Maria Hauser mhauser@genzentrum.lmu.de
//
// Test class for k-mer generation and index table testing.
//

#include <cstdio>
#include <iostream>
#include "SubstitutionMatrix.h"
#include "IndexTable.h"
#include "PrefilteringIndexReader.h"
#include "Parameters.h"

const char* binary_name = "test_indextable";

int main(int argc, const char *argv[]) {
    Parameters &par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 8.0, -0.2f);
    DBReader<unsigned int> dbr(
                               "/Users/mad/Documents/databases/db_small/db_small",
//                               "/Users/mad/Documents/databases/mmseqs_benchmark/benchmarks/protein_search_uniscop/db/mmseqs/db_sw",
                               "/Users/mad/Documents/databases/db_small/db_small.index"
//                               "/Users/mad/Documents/databases/mmseqs_benchmark/benchmarks/protein_search_uniscop/db/mmseqs/db_sw.index"
                               );
    dbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Sequence *s = new Sequence(32000, Sequence::AMINO_ACIDS, &subMat, 6, true, false);
    IndexTable t(subMat.alphabetSize, 6, false);
    PrefilteringIndexReader::fillDatabase(&dbr, s, &t, &subMat, 0, dbr.getSize(), false, 1, NULL, 0, 1);
    t.printStatistics(subMat.int2aa);

    delete s;
    dbr.close();

    return 0;
}

