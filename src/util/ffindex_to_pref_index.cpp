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

int createindex (int argc, const char * argv[])
{
    
    std::string usage("\nCreate index for fast readin.\n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE: <ffindexDB> <ffindexOutDB> \n");


    Parameters par;
    par.parseParameters(argc, argv, usage, par.createindex, 2);

    Debug::setDebugLevel(par.verbosity);
    
    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str());
    dbr.open(DBReader<unsigned int>::SORT);

    BaseMatrix* subMat = Prefiltering::getSubstitutionMatrix(par.scoringMatrixFile, par.alphabetSize, 8.0f, false);

    Sequence seq(par.maxSeqLen, subMat->aa2int, subMat->int2aa, Sequence::AMINO_ACIDS, par.kmerSize, par.spacedKmer);

    PrefilteringIndexReader::createIndexFile(par.db2, par.db2Index, &dbr, &seq, par.split, subMat->alphabetSize, par.kmerSize, par.skip, par.spacedKmer, par.searchMode);

    // write code
    dbr.close();
        
    delete subMat;
    return 0;
}
