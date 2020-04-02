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
    static std::pair<KmerPosition<short> *,size_t > searchInIndex( KmerPosition<short> *kmers, size_t kmersSize, KmerIndex &kmerIndex, int resultDirection);

    template  <int TYPE>
    static void writeResult(DBWriter & dbw, KmerPosition<short> *kmers, size_t kmerCount);


    struct ExtractKmerAndSortResult{
        ExtractKmerAndSortResult(size_t kmerCount, KmerPosition<short> * kmers, size_t adjustedKmer)
                : kmerCount(kmerCount), kmers(kmers), adjustedKmer(adjustedKmer)  {}
        size_t kmerCount;
        KmerPosition<short> * kmers;
        size_t adjustedKmer;
    };
    static ExtractKmerAndSortResult extractKmerAndSort(size_t splitKmerCount, size_t split, size_t splits,
                                                       DBReader<unsigned int> &seqDbr, Parameters &par, BaseMatrix *subMat);
};


#endif //MMSEQS_KMERSEARCH_H
