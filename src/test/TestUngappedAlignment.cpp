#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "Parameters.h"
#include "Sequence.h"
#include "SequenceLookup.h"
#include "SubstitutionMatrix.h"
#include "UngappedAlignment.h"

const char* binary_name = "test_ungappedalignment";
DEFAULT_PARAMETER_SINGLETON_INIT

namespace {

const size_t hitCount = 8;
const unsigned int longTargetLength = 32768;

std::string longTarget(size_t matchingPrefixLength) {
    std::string target(longTargetLength, 'X');
    target.replace(0, matchingPrefixLength, matchingPrefixLength, 'A');
    return target;
}

int runScenario(SubstitutionMatrix &subMat, const unsigned char *remap, bool twoLongTargets) {
    const char *mode = remap == NULL ? "direct" : "identity-remap";
    const char *scenario = twoLongTargets ? "two-long" : "single-long";
    const std::string queryText(32, 'A');
    // Unsorted short lengths give distinct scores and exercise sorted writeback.
    const size_t shortLengths[] = {7, 1, 6, 2, 5, 3, 4};
    std::vector<std::string> targets;
    for (size_t i = 0; i < hitCount - 1; ++i) {
        targets.push_back(std::string(shortLengths[i], 'A'));
    }
    targets.push_back(longTarget(queryText.size()));
    if (twoLongTargets) {
        targets[0] = longTarget(3);
        targets[7] = longTarget(6);
    }

    size_t totalLength = 0;
    for (size_t i = 0; i < targets.size(); ++i) {
        totalLength += targets[i].size();
    }
    SequenceLookup lookup(hitCount, totalLength);
    Sequence target(longTargetLength, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, false);
    for (size_t i = 0; i < targets.size(); ++i) {
        target.mapSequence(i, i, targets[i].c_str(), targets[i].size());
        lookup.addSequence(&target);
    }
    Sequence query(longTargetLength, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, false);
    query.mapSequence(0, 0, queryText.c_str(), queryText.size());
    UngappedAlignment matcher(longTargetLength, &subMat, &lookup, remap);
    matcher.createProfile(&query, NULL);

    CounterResult hits[hitCount] = {};
    int referenceScores[hitCount] = {};
    int failures = 0;
    int firstShortScore = -1;
    bool distinctShortScores = false;
    for (size_t i = 0; i < hitCount; ++i) {
        hits[i].id = i;
        hits[i].diagonal = 0;
        hits[i].count = 0;
        referenceScores[i] = matcher.scoreSingelSequenceByCounterResult(hits[i]);
        const bool isLong = targets[i].size() == longTargetLength;
        if (isLong ? (referenceScores[i] <= 0 || (twoLongTargets && referenceScores[i] >= 255))
                   : (referenceScores[i] < 0 || referenceScores[i] >= 255)) {
            std::cerr << scenario << " " << mode << " hit " << hits[i].id
                      << ": invalid fixture reference score " << referenceScores[i] << '\n';
            ++failures;
        }
        if (!isLong) {
            if (firstShortScore == -1) {
                firstShortScore = referenceScores[i];
            } else if (referenceScores[i] != firstShortScore) {
                distinctShortScores = true;
            }
        }
    }
    if (!distinctShortScores) {
        std::cerr << scenario << " " << mode << ": short reference scores must differ\n";
        ++failures;
    }
    if (twoLongTargets && referenceScores[0] == referenceScores[7]) {
        std::cerr << scenario << " " << mode << ": long reference scores must differ\n";
        ++failures;
    }
    if (failures != 0) {
        return failures;
    }

    // Eight same-diagonal hits fill one eight-lane bin or two four-lane bins.
    for (size_t i = 0; i < hitCount; ++i) {
        hits[i].count = 0;
    }
    matcher.align(hits, hitCount);
    for (size_t i = 0; i < hitCount; ++i) {
        const int expected = std::min(255, referenceScores[i]);
        const int actual = hits[i].count;
        if (actual != expected) {
            std::cerr << scenario << " " << mode << " hit " << hits[i].id
                      << ": expected " << expected << ", actual " << actual << '\n';
            ++failures;
        }
    }
    std::cout << scenario << " " << mode << ": "
              << (failures == 0 ? "PASS" : "FAIL") << '\n';
    return failures;
}

} // namespace

int main(int, const char**) {
#ifdef AVX2
    std::cout << "Compiled mode: AVX2 defined (eight lanes)\n";
#else
    std::cout << "Compiled mode: AVX2 undefined (four lanes)\n";
#endif
    Parameters &par = Parameters::getInstance();
    par.initMatrices();
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2);
    unsigned char identityRemap[256];
    for (size_t i = 0; i < sizeof(identityRemap); ++i) {
        identityRemap[i] = static_cast<unsigned char>(i);
    }

    int failures = 0;
    failures += runScenario(subMat, NULL, false);
    failures += runScenario(subMat, identityRemap, false);
    failures += runScenario(subMat, NULL, true);
    failures += runScenario(subMat, identityRemap, true);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
