#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "PrefilteringIndexReader.h"
#include "Prefiltering.h"
#include "Parameters.h"

#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>      // std::ifstream
#include <vector>
#include <utility>      // std::pair

int createindex (int argc, const char **argv, const Command& command)
{
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 1);


#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str());
    dbr.open(DBReader<unsigned int>::NOSORT);

    BaseMatrix* subMat = Prefiltering::getSubstitutionMatrix(par.scoringMatrixFile, par.alphabetSize, 8.0f, false);

    PrefilteringIndexReader::createIndexFile(par.db1, &dbr, subMat, par.maxSeqLen, par.spacedKmer, par.compBiasCorrection, par.split,
                                             subMat->alphabetSize, par.kmerSize, par.diagonalScoring, par.threads);

    // write code
    dbr.close();
        
    delete subMat;
    return 0;
}
