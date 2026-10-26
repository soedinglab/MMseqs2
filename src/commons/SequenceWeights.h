//
// Created by annika on 09.11.22.
//

#ifndef MMSEQS_SEQUENCEWEIGHTS_H
#define MMSEQS_SEQUENCEWEIGHTS_H

#include "IndexTypes.h"

class SequenceWeights{
public:
    struct WeightIndexEntry {
        DBKeyType id;
        float weight;

        static bool compareByIdOnly(const WeightIndexEntry &x, const WeightIndexEntry &y) {
            return x.id <= y.id;
        }
    };

    WeightIndexEntry *weightIndex;
    size_t indexSize;

    SequenceWeights(const char* dataFileName);

    ~SequenceWeights();

    float getWeightById(DBKeyType id);
};


#endif //MMSEQS_SEQUENCEWEIGHTS_H
