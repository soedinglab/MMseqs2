//
// Created by Martin Steinegger on 2019-01-04.
//

#ifndef MMSEQS_KMERSEARCH_H
#define MMSEQS_KMERSEARCH_H

#include "kmermatcher.h"
#include "KmerIndex.h"

class KmerSearch{

public:
    template  <int TYPE>
    static std::pair<KmerPosition *,size_t > searchInIndex( KmerPosition *kmers, size_t kmersSize, KmerIndex &kmerIndex);

    template  <int TYPE>
    static void writeResult(DBWriter & dbw, KmerPosition *kmers, size_t kmerCount);

    static std::pair<size_t, KmerPosition *> extractKmerAndSort(size_t splitKmerCount, size_t split, size_t splits,
                                                                DBReader<unsigned int> &seqDbr, Parameters &par, BaseMatrix *subMat,
                                                                size_t KMER_SIZE, size_t chooseTopKmer);
};


#endif //MMSEQS_KMERSEARCH_H
