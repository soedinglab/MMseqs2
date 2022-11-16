#include "Parameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "IndexReader.h"
#include "NucleotideMatrix.h"

#include <climits>

#ifdef OPENMP
#include <omp.h>
#endif

int countkmer(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.verbosity = 1;
    par.kmerSize = 5;
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, true, 0, 0);
    std::vector<std::string> ids = Util::split(par.idList, ",");
    int indexSrcType = IndexReader::SEQUENCES;

    IndexReader reader(par.db1, par.threads, indexSrcType, 0);
    int seqType = reader.sequenceReader->getDbtype();
    BaseMatrix * subMat;
    size_t isNucl=Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES);
    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES)) {
        subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    } else {
        subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);
    }
    size_t maxLen = 0;
    for(size_t i = 0; i < reader.sequenceReader->getSize(); i++){
        maxLen = std::max(maxLen, reader.sequenceReader->getSeqLen(i));
    }
    size_t idxSize = MathUtil::ipow<size_t>(subMat->alphabetSize-1, par.kmerSize);
    unsigned int * kmerCountTable=new unsigned int[idxSize];
    memset(kmerCountTable, 0, sizeof(unsigned int)*idxSize);
#pragma omp parallel
    {
        Indexer idx(subMat->alphabetSize-1, par.kmerSize);
        Sequence s(maxLen, seqType, subMat,
                          par.kmerSize, par.spacedKmer, false);

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader.sequenceReader->getSize(); i++) {
            char *data = reader.sequenceReader->getData(i, 0);
            s.mapSequence(i, 0, data, reader.sequenceReader->getSeqLen(i));
            while (s.hasNextKmer()) {
                const unsigned char *kmer = s.nextKmer();
                if(s.kmerContainsX()){
                    continue;
                }

                size_t kmerIdx = (isNucl) ? Indexer::computeKmerIdx(kmer, par.kmerSize) : idx.int2index(kmer, 0, par.kmerSize);
                __sync_fetch_and_add(&(kmerCountTable[kmerIdx]), 1);
            }
        }
    }
    Indexer idx(subMat->alphabetSize-1, par.kmerSize);
    for(size_t i = 0; i < idxSize; i++){
        std::cout << i << "\t";
        if(isNucl){
            Indexer::printKmer(i, par.kmerSize);
        }else{
            idx.index2int(idx.workspace, i, par.kmerSize);
            for(int k = 0; k < par.kmerSize; k++){
                std::cout << subMat->num2aa[idx.workspace[k]];
            }
        }
        std::cout << "\t" << kmerCountTable[i] << std::endl;
    }
    delete [] kmerCountTable;
    EXIT(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}
