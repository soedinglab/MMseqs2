#include <cfloat>
#include <sstream>
#include "Aggregation.h"
#include "Util.h"
#include "Parameters.h"
#include <algorithm>
#include <tgmath.h>
#include <iostream>
#include <fstream>
#include "Debug.h"

#ifdef OPENMP
#include <omp.h>
#endif

#define EVAL_COLUMN 3





struct compareByStart {
    bool operator()(const std::pair <long,long> &lhs,
                    const std::pair <long,long> &rhs) const {
        return (lhs.first < rhs.first);
    }
};

// TODO : Get rid of stringStreams
Aggregation::Aggregation(std::string argInputDBname, std::string argOutputDBname, unsigned int argNbrThread,
                         size_t argTargetSetColumn) : inputDBname(std::move(argInputDBname)), outputDBname(std::move(argOutputDBname)),
                                                    nbrThread(argNbrThread), targetSetColumn(argTargetSetColumn){}

// build a map with the value in [target column] field as a key and the rest of the line, cut in fields, as values
bool Aggregation::buildMap(std::stringstream& dataStream, std::map<unsigned int, std::vector<std::vector<std::string> > > &dataToAggregate){
    dataToAggregate.clear() ;
    std::string singleLine;
    while(std::getline(dataStream, singleLine)){
        if(!singleLine.empty()) {
            std::vector<std::string> splittedLine = Util::split(singleLine, "\t");
            dataToAggregate[stoul(splittedLine[this->targetSetColumn])].push_back(splittedLine);
        }
    }
    return true;
}

// General workflow for search Result processing
void Aggregation::runAggregate(){

    std::string inputDBIndex = inputDBname + ".index";
    DBReader<unsigned int>* inputDB = new DBReader<unsigned int> (inputDBname.c_str(), inputDBIndex.c_str()) ;
    inputDB->open(DBReader<unsigned int>::NOSORT);

    std::string outputDBIndex = outputDBname + ".index";
    DBWriter* outputDB = new DBWriter (outputDBname.c_str(), outputDBIndex.c_str(), nbrThread) ;
    outputDB->open();

#ifdef OPENMP
    omp_set_num_threads(this->nbrThread);
#endif
    std::map<unsigned int, std::vector<std::vector<std::string> > > dataToMerge ;
    unsigned int thread_idx = 0;
#pragma omp parallel for private(thread_idx,dataToMerge)
    for (size_t i = 0; i < inputDB->getSize(); i++) {
        //std::cout << "Debut \n" ;
#ifdef OPENMP
        thread_idx = (unsigned int)omp_get_thread_num();
#endif
        unsigned int key = inputDB->getDbKey(i);
        //std::cout << "key : " << key << "\n" ;
        std::string geneResults(inputDB->getDataByDBKey(key));
        std::stringstream geneResultsStream(geneResults);
        //std::cout << A << "\n" ;


        if (this->buildMap(geneResultsStream,dataToMerge)){
            std::string outputBuffer;
            outputBuffer.reserve(1048576) ;
            for( std::map<unsigned int, std::vector<std::vector<std::string>>>::iterator it=dataToMerge.begin() ; it!= dataToMerge.end() ; ++it) {
                aggregParams params = {.querySetKey = key, .targetSetKey = it->first};
                outputBuffer.append(this->aggregateEntry(it->second, &params));
                outputBuffer.append("\n");
            }

            outputDB->writeData(outputBuffer.c_str(), outputBuffer.length(), key, thread_idx);


        } else {
            Debug(Debug::ERROR) << "buildMap failed \n";
            EXIT(EXIT_FAILURE);
        }
    }
    outputDB->close();
    delete outputDB;
}


BestHitAggregator::BestHitAggregator(std::string argInputDBname, std::string argOutputDBname,std::string argTargetSetSizeName,
                                     size_t argTargetColumn, unsigned int argNbrThread,  size_t argScoreColumn, bool argSimpleBestHitMode) :
                                     Aggregation(std::move(argInputDBname), std::move(argOutputDBname), argNbrThread,
                                     argTargetColumn), scoreColumn(argScoreColumn), targetSetSizeName(std::move(argTargetSetSizeName)), simpleBestHitMode(argSimpleBestHitMode){
    std::string sizeDBIndex = targetSetSizeName + ".index";
    this->targetSetSizeDB = new DBReader<unsigned int> (targetSetSizeName.c_str(), sizeDBIndex.c_str()) ;
    this->targetSetSizeDB->open(DBReader<unsigned int>::NOSORT);
}

// Only Keep the best hits of a protein against each Target Set
std::string BestHitAggregator::aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* vparams) {

    std::string buffer;
    double bestScore = 0;//DBL_MAX;
    size_t maxId = 0;
    size_t bestEvalId = 0;
    double secondBestScore = 0;
    double correctedPval = 0;

    double evalDouble = 0 ;
    double bestEval = DBL_MAX ;
    // Look for the lowest p-value and retain only this line
    // dataToAggregate = [nbrTargetGene][Field of result]
    for(size_t i =0 ; i < dataToAggregate.size() ; i++){
        std::istringstream score(dataToAggregate[i][this->scoreColumn]) ;
        double scoreInDouble;
        score >> scoreInDouble;
        std::istringstream eVal(dataToAggregate[i][EVAL_COLUMN]) ;
        eVal >> evalDouble ;

        //std::cout << doubleEval << "\n" ;
        if (scoreInDouble > bestScore) {
            secondBestScore = bestScore;
            bestScore = scoreInDouble;
            maxId = i;
        } else if (scoreInDouble > secondBestScore) {
            secondBestScore = scoreInDouble;
        }
        if (bestEval>evalDouble) {
            bestEval = evalDouble;
            bestEvalId = i;
        }

    }

    unsigned int nbrGenes = (unsigned int)std::stoul(this->targetSetSizeDB->getDataByDBKey(
            static_cast<unsigned int>(std::stoul(dataToAggregate[maxId][targetSetColumn]))));


    if(dataToAggregate.size() < 2) {
        secondBestScore = 2.0*log((nbrGenes+1)/2)/log(2.0);//(2.00/(nbrGenes+1.00)) ; // (2.00*nbrGenes/(nbrGenes+1.00)) ;
    }

    if (simpleBestHitMode) {
        maxId = bestEvalId;
        correctedPval = bestEval/nbrGenes;
    } else { // ortholog mode
        correctedPval = pow(2.0, secondBestScore / 2 - bestScore / 2);//exp(log(secondBestEval / nbrGenes) - log(secondBestEval / nbrGenes));
    }


    char tmpBuf[15];
    sprintf(tmpBuf,"%.3E",correctedPval);
    dataToAggregate[maxId][this->scoreColumn] = std::string(tmpBuf);

    // Aggregate the full line in one string
    for (std::vector<std::string>::iterator it = dataToAggregate[maxId].begin();
         it != dataToAggregate[maxId].end(); ++it) {
        buffer.append(*it + "\t");
    }
    return buffer;
}


PvalAggregator::PvalAggregator(std::string argInputDBname, std::string argOutputDBname, unsigned int arg_nbrThread,
                               std::string argQuerySetSizeDBname,std::string argTargetSetSizeDBname, size_t argTargetSetColumn, float argAlpha, size_t argScoreColumn) :
        Aggregation(std::move(argInputDBname), std::move(argOutputDBname), arg_nbrThread, argTargetSetColumn), scoreColumn(argScoreColumn),
        querySetSizeDBname(std::move(argQuerySetSizeDBname)),targetSetSizeDBname(std::move(argTargetSetSizeDBname)), alpha(argAlpha){

    std::string sizeDBIndex = querySetSizeDBname + ".index";
    this->querySetSizeDB = new DBReader<unsigned int> (querySetSizeDBname.c_str(), sizeDBIndex.c_str()) ;
    this->querySetSizeDB->open(DBReader<unsigned int>::NOSORT);

    sizeDBIndex = targetSetSizeDBname + ".index";
    this->targetSetSizeDB = new DBReader<unsigned int> (targetSetSizeDBname.c_str(), sizeDBIndex.c_str()) ;
    this->targetSetSizeDB->open(DBReader<unsigned int>::NOSORT);
}

/*double kthOrderProba(size_t k, size_t N, double p) {
    return (double)0.5*std::erfc(((double)k+0.5-N*p)/sqrt(N*p*(1.0-p)*2.0)) ;
}

double BinCoeff (int M, int k) {
    double result;
    result = exp(lgamma(M+1)-lgamma(M-k+1)-lgamma(k+1)) ;
    return result;
}*/
double LBinCoeff (int M, int k) {
    double result;
    result = lgamma(M+1)-lgamma(M-k+1)-lgamma(k+1);
    return result;
}

double factorial (size_t i) {
    double fac = 1;
    for (int j = 1; j < i+1; j++) {
        fac *= j;
    }
    return fac;
}

//Get all result of a single Query Set VS a Single Target Set and return the multiple-match p-value for it
std::string PvalAggregator::aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* params){
    std::stringstream Buffer ;
    unsigned int M = (unsigned int)std::stoul(querySetSizeDB->getDataByDBKey(params->querySetKey)) ;
    unsigned nbrTargetGene = (unsigned int)std::stoul(targetSetSizeDB->getDataByDBKey(params->targetSetKey)) ;
    double dP0 = this->alpha/nbrTargetGene;//
    double updatedPval;
    double r = 0;
    double pValInDouble;
    size_t k =0;
    std::vector<double> pvals;
    for (auto &i : dataToAggregate) {
        pValInDouble = std::strtod(i[scoreColumn].c_str(), nullptr);
        if (pValInDouble < dP0){
            k++;
            r-=log(pValInDouble/dP0);
        }
    }


    double leftSum = 0;
    double rightSum = 0;
    for(size_t i =0 ; i < k ; i++) { // LeftSum
        rightSum = 0.0;
        for(size_t j = i+1 ; j < k+1 ; j++) { // RightSum
            rightSum += exp(LBinCoeff(M, j) + j*log(dP0) + (M - j)*log(1.0 - dP0));//BinCoeff(M, j) * pow(P0, j) * pow((1.0 - P0), (M - j));
        }
        leftSum += (pow(r, i) / factorial(i))*rightSum;
    }

    std::vector<double> ppvals;
    double I=0;
    if(r==0){I=1;}
    updatedPval = (1.0-pow((1.0-dP0), M))*I + exp(-r)*leftSum ;
    double updatedEval = updatedPval*M;

    if (std::isinf(r)) {
        Debug(Debug::WARNING) << "R is infinity !\n";
        updatedEval = 0;
    }

    if(updatedEval == 0) {
        //Debug(Debug::WARNING) << "Eval is 0 !\n";
        for (auto &i : dataToAggregate) {
            pValInDouble = std::strtod(i[scoreColumn].c_str(), nullptr);
        }
    }

    Buffer << params->targetSetKey << "\t" << updatedEval ;//* targetSetSizeDB->getSize() ;
    return Buffer.str();
}

GetHitDistance::GetHitDistance(std::string arg_inputDBname, std::string arg_outputDBname, std::string argTargetSetSizeDB, std::string argQuerySetSizeDB,std::string argTargetGenomesDB,
                               unsigned int arg_nbrThread, float argAlpha, size_t arg_targetColumn)
        : Aggregation(std::move(arg_inputDBname), std::move(arg_outputDBname), arg_nbrThread, arg_targetColumn), targetSetSizeDBname(std::move(argTargetSetSizeDB)),
          querySetSizeDBname(std::move(argQuerySetSizeDB)), alpha(argAlpha){

    std::string nbrGeneQueryDBIndex = querySetSizeDBname + ".index";
    this->querySetSizeDB = new DBReader<unsigned int> (querySetSizeDBname.c_str(), nbrGeneQueryDBIndex.c_str());
    this->querySetSizeDB->open(DBReader<unsigned int>::NOSORT);

    std::string nbrGeneDBIndex = targetSetSizeDBname + ".index";
    this->targetSetSizeDB = new DBReader<unsigned int> (targetSetSizeDBname.c_str(), nbrGeneDBIndex.c_str());
    this->targetSetSizeDB->open(DBReader<unsigned int>::NOSORT);

    std::string targetGenomnesDBIndex = argTargetGenomesDB + ".index";
    this->targetSetGenomes = new DBReader<unsigned int> (argTargetGenomesDB.c_str(), targetGenomnesDBIndex.c_str()) ;
    this->targetSetGenomes->open(DBReader<unsigned int>::USE_INDEX);
}

GetHitDistance::~GetHitDistance() {

    this->targetSetGenomes->close();
    this->targetSetSizeDB->close();
    this->querySetSizeDB->close();

    delete this->targetSetGenomes;
    delete this->targetSetSizeDB;
    delete this->querySetSizeDB;
}
/*struct GenePosition {
    size_t start;
    size_t end;
};
*/

double geomMedianProba(size_t i, size_t K, double rate) {
    return pow((1.0-pow(1.0-rate,i))* pow(1.0-rate,(i+1)),(K/2));
}

double normalizedSpreadPval(size_t median,size_t K, double rate) {
    double s = 0.0;

    if (median >= 1.0/rate)
        return 1.0;


    for (size_t i = 0; i<=median; i++) {
        s+=geomMedianProba(i, K, rate);
    }

    double sNorm = s;
    for (size_t i = median+1; i<=1.0/rate; i++) {
        sNorm += geomMedianProba(i, K, rate);
    }

    // would need *exp(LBinCoeff(K,K/2)); but does not matter since we normalize
    return s/sNorm;
}
double spreadPval(long median, double rate, size_t K)
{
    if (median < 1)
        median = 1.0;

    if (median > 100000)
        return 1.0;

    return normalizedSpreadPval(median,K,rate);
}


// return the median of the distance between genes of dataToAggregate
std::string GetHitDistance::aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate,
                                                 aggregParams *params) {
    long nbrTargetGene = std::stol(targetSetSizeDB->getDataByDBKey(params->targetSetKey));
    double P0 = this->alpha/nbrTargetGene;
    int i = 0;
    unsigned int indexOfInterSpaceToReturn =0;
    std::vector<std::pair<long,long> > genesPositions;
    std::vector<long> interGeneSpaces;
    long currentInterGenePosition = 0;
    std::stringstream buffer;
    double Pval;
    double meanEval = 0;
    std::string positionsStr;
    std::string genesID;
    std::string eVals;
    unsigned int nbrGoodEvals = 0;
    std::vector<unsigned long> checkGenesID; //DELETE

    for (auto &it : dataToAggregate) {
        Pval=std::strtod(it[3].c_str(), nullptr);
        if(Pval < P0) {
            checkGenesID.emplace_back(std::strtoul(it[0].c_str(), nullptr, 0));
            meanEval+=log10(std::strtod(it[3].c_str(), nullptr));
            eVals += it[3]+",";
            if(std::strtod(it[3].c_str(), nullptr)<1e-10){nbrGoodEvals++;}
            unsigned long start = static_cast<unsigned long>(std::stol(it[8]));
            unsigned long stop = static_cast<unsigned long>(std::stol(it[10]));
            genesID += it[0]+",";
            genesPositions.emplace_back(std::make_pair(start, stop));
            positionsStr += std::to_string(start) + "," + std::to_string(stop) + ",";
            i++;
        }
    }
    std::sort(genesPositions.begin(), genesPositions.end(), compareByStart());

    double genomeSize = (targetSetGenomes->getSeqLens(targetSetGenomes->getId(params->targetSetKey)) - 2);
    double rate = ((double) i) / genomeSize;


    if(i>1) {
        for (size_t pos = 0; pos < i-1; pos++) {
            currentInterGenePosition = genesPositions[pos+1].first - genesPositions[pos].second;
            interGeneSpaces.push_back(currentInterGenePosition);
        }
        std::sort(begin(interGeneSpaces), end(interGeneSpaces));
        //if odd number
        if(interGeneSpaces.size() == 1) {
            indexOfInterSpaceToReturn = 0;
        }
        else if(interGeneSpaces.size() % 2 && interGeneSpaces.size() !=1) {
            indexOfInterSpaceToReturn = static_cast<unsigned int>((interGeneSpaces.size() + 1) / 2);
        }
        else{
            indexOfInterSpaceToReturn = static_cast<unsigned int>(interGeneSpaces.size() / 2);
        }
        // T    HitDist     nbrGeneT        nbrHits
        buffer << params->targetSetKey << "\t" << interGeneSpaces[indexOfInterSpaceToReturn] << "\t" << std::stol(targetSetSizeDB->getDataByDBKey(params->targetSetKey))
               << "\t" << i << "\t" << meanEval/i << "\t" << std::stol(querySetSizeDB->getDataByDBKey(params->querySetKey))  << "\t"
               << spreadPval(interGeneSpaces[indexOfInterSpaceToReturn],rate,i) << "\t" << genomeSize << "\t" << nbrGoodEvals;
    }
    else {
        buffer << params->targetSetKey << "\t0" << "\t" << std::stol(targetSetSizeDB->getDataByDBKey(params->targetSetKey)) << "\t" << i
               << "\t" << meanEval/i << "\t" << std::stol(querySetSizeDB->getDataByDBKey(params->querySetKey)) << "\t1.0" << "\t"
               << genomeSize << "\t" << nbrGoodEvals; //NaMM
    }
    buffer<<"\t" <<positionsStr << "\t" << genesID << "\t" <<eVals;

    std::string bufferString= buffer.str();
    return bufferString;
}
