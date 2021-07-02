#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrix.h"
#include "ScoreMatrix.h"
#include "Parameters.h"
#include "Debug.h"

const char* binary_name = "test_scorematrixserialization";

int main (int, const char**) {
    Parameters& par = Parameters::getInstance();
    par.initMatrices();
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 8.0, 0);
    ScoreMatrix extMattwo = ExtendedSubstitutionMatrix::calcScoreMatrix(subMat, 2);

    Debug(Debug::INFO) << extMattwo.elementSize << " " << extMattwo.rowSize  << " "
                       << extMattwo.score[0]    << " " << extMattwo.index[0] << "\n";

    char* serialized = ScoreMatrix::serialize(extMattwo);
    ExtendedSubstitutionMatrix::freeScoreMatrix(extMattwo);

    ScoreMatrix unserialized = ScoreMatrix::unserialize(serialized, subMat.alphabetSize, 2);
    Debug(Debug::INFO) << unserialized.elementSize << " " << unserialized.rowSize  << " "
                       << unserialized.score[0]    << " " << unserialized.index[0] << "\n";
    free(serialized);

    return EXIT_SUCCESS;
}
