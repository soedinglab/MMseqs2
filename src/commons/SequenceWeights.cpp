//
// Created by annika on 09.11.22.
//
#include "SequenceWeights.h"
#include "Util.h"
#include "Debug.h"

#include <algorithm>
#include <fstream>

SequenceWeights::SequenceWeights(const char* dataFileName) {

    //parse file and fill weightIndex
    std::ifstream tsv(dataFileName);
    if (tsv.fail()) {
        Debug(Debug::ERROR) << "File " << dataFileName << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    char keyData[255];
    std::string line;
    this->indexSize = 0;
    unsigned int  pos = 0;
    while(std::getline(tsv, line)) {
        this->indexSize++;
    }

    this->weightIndex = new WeightIndexEntry[this->indexSize];

    tsv.clear();
    tsv.seekg(0);

    while(std::getline(tsv, line)) {
        char *current = (char *) line.c_str();
        Util::parseKey(current, keyData);
        const std::string key(keyData);
        unsigned int keyId = strtoull(key.c_str(), NULL, 10);

        char *restStart = current + key.length();
        restStart = restStart + Util::skipWhitespace(restStart);
        float weight = static_cast<float>(strtod(restStart, NULL));
        this->weightIndex[pos].id = keyId;
        this->weightIndex[pos].weight = weight;
        pos++;
    }
}

float SequenceWeights::getWeightById(unsigned int id) {

    WeightIndexEntry val;
    val.id = id;
    size_t pos = std::upper_bound(weightIndex, weightIndex + indexSize, val, WeightIndexEntry::compareByIdOnly) - weightIndex;
    return (pos < indexSize && weightIndex[pos].id == id ) ? weightIndex[pos].weight : std::numeric_limits<float>::min();
}

SequenceWeights::~SequenceWeights() {
    delete[] weightIndex;
}

