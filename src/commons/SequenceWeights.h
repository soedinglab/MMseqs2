//
// Created by annika on 09.11.22.
//

#ifndef MMSEQS_SEQUENCEWEIGHTS_H
#define MMSEQS_SEQUENCEWEIGHTS_H
#include "Parameters.h"

class SequenceWeights{
public:
    struct WeightIndexEntry {
        KeyType id;
        float weight;

        static bool compareByIdOnly(const WeightIndexEntry &x, const WeightIndexEntry &y) {
            return x.id <= y.id;
        }
    };

    WeightIndexEntry *weightIndex;
    unsigned int indexSize;

    SequenceWeights(const char* dataFileName);

    ~SequenceWeights();

    float getWeightById(KeyType id);
};


#endif //MMSEQS_SEQUENCEWEIGHTS_H
