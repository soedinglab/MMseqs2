
#include "DistanceCalculator.h"


void DistanceCalculator::prepareGlobalAliParam(const BaseMatrix &subMat)
{
    globalAliMu = 0;
    globalAliSigma = 0;

    for (size_t i = 0; i<subMat.alphabetSize - 1;i++)
    {
        
        for (size_t j = 0; j<subMat.alphabetSize - 1;j++)
        {
            globalAliMu += subMat.pBack[i] * subMat.pBack[j] * log(subMat.probMatrix[i][j] / (subMat.pBack[i] * subMat.pBack[j]) );
        }
    }
    
    for (size_t i = 0; i<subMat.alphabetSize - 1;i++)
    {
        
        for (size_t j = 0; j<subMat.alphabetSize - 1;j++)
        {
            double distToMean = (log(subMat.probMatrix[i][j] / (subMat.pBack[i] * subMat.pBack[j]) ) - globalAliMu);
            globalAliSigma += subMat.pBack[i] * subMat.pBack[j] * distToMean*distToMean;
        }
    }
    globalAliSigma = sqrt(globalAliSigma);
    
}


double DistanceCalculator::getPvalGlobalAli(float score)
{
    
    return 0.5 - 0.5*erf((score-globalAliMu)/(sqrt(2.0)*globalAliSigma));
}