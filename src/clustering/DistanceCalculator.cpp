
#include "DistanceCalculator.h"


void DistanceCalculator::prepareGlobalAliParam(const BaseMatrix &subMat)
{
    globalAliMu = 0;
    globalAliSigma = 0;

    for (size_t i = 0; i<subMat.alphabetSize - 1;i++)
    {
        
        for (size_t j = 0; j<subMat.alphabetSize - 1;j++)
        {
            globalAliMu += subMat.pBack[i] * subMat.pBack[j] * subMat.subMatrix[i][j];
        }
    }
    
    for (size_t i = 0; i<subMat.alphabetSize - 1;i++)
    {
        
        for (size_t j = 0; j<subMat.alphabetSize - 1;j++)
        {
            double distToMean = (subMat.subMatrix[i][j] - globalAliMu);
            globalAliSigma += subMat.pBack[i] * subMat.pBack[j] * distToMean*distToMean;
        }
    }
    globalAliSigma = sqrt(globalAliSigma);
    
}


double DistanceCalculator::getPvalGlobalAli(float score,size_t len)
{
    return 0.5 - 0.5*erf((score/len-globalAliMu)/(sqrt(2.0/sqrt((float)len))*globalAliSigma));
}