#include <cfloat>
#include <sstream>
#include "Aggregation.h"
#include "Util.h"
#include "Parameters.h"
#include <iostream>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#include <tgmath.h>

#endif

Aggregation::Aggregation(std::string arg_inputDBname, std::string arg_outputDBname, unsigned int arg_nbrThread,
                         size_t arg_targetColumn) : inputDBname(std::move(arg_inputDBname)), outputDBname(std::move(arg_outputDBname)),
                                                    nbrThread(arg_nbrThread), targetColumn(arg_targetColumn){}


bool Aggregation::buildMap(std::stringstream& dataStream, std::map<unsigned int, std::vector<std::vector<std::string> > > &dataToAggregate){
    dataToAggregate.clear() ;
    std::string singleLine;
    while(std::getline(dataStream, singleLine)){
        if(!singleLine.empty()) {
            std::vector<std::string> splittedLine = Util::split(singleLine, "\t");
            dataToAggregate[stoul(splittedLine[this->targetColumn])].push_back(splittedLine);
        }
    }
    return true;
}


void Aggregation::runAggregate(){

    std::string inputDBIndex = inputDBname + ".index";
    auto inputDB = new DBReader<unsigned int> (inputDBname.c_str(), inputDBIndex.c_str()) ;
    inputDB->open(DBReader<unsigned int>::NOSORT);

    std::string outputDBIndex = outputDBname + ".index";
    auto outputDB = new DBWriter (outputDBname.c_str(), outputDBIndex.c_str(), nbrThread) ;
    outputDB->open();

#ifdef OPENMP
    omp_set_num_threads(1);
#endif

    unsigned int thread_idx = 0;
#pragma omp parallel for private(thread_idx)
    for (size_t i = 0; i < inputDB->getSize(); i++) {
        //std::cout << "Debut \n" ;
#ifdef OPENMP
        thread_idx = (unsigned int)omp_get_thread_num();
#endif
        unsigned int key = inputDB->getDbKey(i);
        //std::cout << "key : " << key << "\n" ;
        std::stringstream geneResultsStream(inputDB->getDataByDBKey(key));
        std::string A = geneResultsStream.str() ;
        //std::cout << A << "\n" ;
        std::map<unsigned int, std::vector<std::vector<std::string> > > dataToMerge;


        if (this->buildMap(geneResultsStream,dataToMerge)){
            std::string outputBuffer;
            for( auto &it : dataToMerge) {
                //std::cout << "For auto" << "\n" ;
                aggregParams params = {.querySetKey = key, .targetSetKey = it.first};
                outputBuffer.append(this->aggregateEntry(it.second, &params));
                outputBuffer.append("\n");
            }
            outputDB->writeData(outputBuffer.c_str(), outputBuffer.length(), key, thread_idx);
        } else {
            std::cout << "buildMap Failed \n";
        }
    }
    outputDB->close();
    delete outputDB;
}


bestHitAggregator::bestHitAggregator(std::string arg_inputDBname, std::string arg_outputDBname,std::string arg_querySetSizeName,
                                     unsigned int arg_nbrThread, size_t arg_targetColumn, size_t arg_scoreColumn,
                                     DBReader<unsigned int>* querySetSize) :
                                     Aggregation(std::move(arg_inputDBname), std::move(arg_outputDBname), arg_nbrThread,
                                     arg_targetColumn), pValColumn(arg_scoreColumn), querySetSizeName(std::move(arg_querySetSizeName)){
    std::string sizeDBIndex = querySetSizeName + ".index";
    this->querySetSize = new DBReader<unsigned int> (querySetSizeName.c_str(), sizeDBIndex.c_str()) ;
    this->querySetSize->open(DBReader<unsigned int>::NOSORT);
}

// Process a single Search Entry (1 gene) and only keep the best Hit
std::string bestHitAggregator::aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* vparams) {



    std::string buffer ;
    auto bestEval = DBL_MAX;
    size_t maxId = 0;
    double secondBestEval = 0 ;
    double correctedPval = 0 ;

    // Look for the lowest e-value and retain only this line
    // dataToAggregate = [nbrTargetGene][Field of result]
    for(size_t i =0 ; i < dataToAggregate.size() ; i++){
        std::istringstream score(dataToAggregate[i][this->pValColumn]) ;
        double pValInDouble ;
        score >> pValInDouble ;
        if (pValInDouble < bestEval) {
            secondBestEval = bestEval ;
            bestEval = pValInDouble;
            maxId = i;
        } else if (pValInDouble < secondBestEval) {
            secondBestEval = pValInDouble;
        }
    }


    //Add parameters to choose between the 2 possibilities

    unsigned int nbrGenes = (unsigned int)std::stoul(this->querySetSize->getDataByDBKey(
            static_cast<unsigned int>(std::stoul(dataToAggregate[maxId][10])))) ;
    if(dataToAggregate.size() < 2) {
        secondBestEval = (2.00*nbrGenes/(nbrGenes+1.00)) ;
    }

    correctedPval = exp(log(bestEval / nbrGenes) - log(secondBestEval / nbrGenes));
    //std::cout << "maxVal : " << bestPval << " \t secondMaxVal :" << secondBestPval << "\t Number of entries ; " << dataToAggregate.size() << " corrected:"<<correctedPval << "\n" ;
    if (correctedPval < 1e-300) {
        char tmpBuf[15];
        sprintf(tmpBuf,"%.3E",correctedPval);
        std::cout << "Corrected Pval : " << correctedPval << "\n" << "char CorrectedPval : " << tmpBuf << "\n" ;
        double Test = 0 ;
        Test = std::atof(tmpBuf) ;
        Test = std::strtod(tmpBuf, nullptr) ;
        std::cout << " atof : " << Test << "\n" ;
    }
    char tmpBuf[15];
    sprintf(tmpBuf,"%.3E",correctedPval);
    dataToAggregate[maxId][this->pValColumn] = std::string(tmpBuf);



    // Aggregate the full line in one string
    for(auto &it : dataToAggregate[maxId]){
        buffer.append(it + "\t");
    }
    return buffer ;
}

pvalAggregator::pvalAggregator(std::string arg_inputDBname, std::string arg_outputDBname, unsigned int arg_nbrThread,
                               std::string arg_querySetSizeDBname, std::string arg_targetSetSizeDBname,
                               DBReader<unsigned int>* querySetSizeDB, DBReader<unsigned int>* targetSetSizeDB,
                               size_t arg_targetColumn, size_t arg_scoreColumn) :
        Aggregation(std::move(arg_inputDBname), std::move(arg_outputDBname), arg_nbrThread, arg_targetColumn), pValColumn(arg_scoreColumn),
        querySetSizeDBname(std::move(arg_querySetSizeDBname)), targetSetSizeDBname(std::move(arg_targetSetSizeDBname)){

    std::string sizeDBIndex = querySetSizeDBname + ".index";
    this->querySetSizeDB = new DBReader<unsigned int> (querySetSizeDBname.c_str(), sizeDBIndex.c_str()) ;
    this->querySetSizeDB->open(DBReader<unsigned int>::NOSORT);/*

    sizeDBIndex = targetSetSizeDBname + ".index";
    this->targetSetSizeDB = new DBReader<unsigned int> (targetSetSizeDBname.c_str(), sizeDBIndex.c_str()) ;
    this->targetSetSizeDB->open(DBReader<unsigned int>::NOSORT);*/

    /*targetGlobalSize = 0;
    for (size_t i = 0; i<this->targetSetSizeDB->getSize(); i++) {
        size_t curGenomeSize = std::stoul(this->targetSetSizeDB->getData(i));
        targetGlobalSize += curGenomeSize;
    }*/
}

/*double kthOrderProba(size_t k, size_t N, double p) {
    return (double)0.5*std::erfc(((double)k+0.5-N*p)/sqrt(N*p*(1.0-p)*2.0)) ;
}*/

double BinCoeff (int M, int k) {
    double result;
    result = exp(lgamma(M+1)-lgamma(M-k+1)-lgamma(k+1)) ;
    return result;
}

double factorial (size_t i) {
    double fac = 1;
    for (int j = 1; j < i+1; j++) {
        fac *= j;
    }
    return fac;
}
std::string pvalAggregator::aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* params){
    std::stringstream Buffer ;
    unsigned int M = (unsigned int)std::stoul(querySetSizeDB->getDataByDBKey(params->querySetKey)) ;
    //unsigned int nbrTargetSetElements = (unsigned int)std::stoul(targetSetSizeDB->getDataByDBKey(params->targetSetKey)) ;
    //unsigned int crossLinks = nbrTargetSetElements*nbrQuerySetElements;
    double P0 = 0.001 ;
    double updatedPval;
    double r = 0 ;
    double pValInDouble ;
    size_t K =0 ;
    std::vector<double> pvals;
    for (auto &i : dataToAggregate) {
        pValInDouble = std::strtod(i[pValColumn].c_str(), nullptr) ;
        //std::cout << "i[pValColumn] : " << pValInDouble << "\n" ;
        if (pValInDouble < P0){// pvals.push_back(std::stod(i[pValColumn]) / targetGlobalSize);
            K++ ;
            std::cout << " pVal : " << pValInDouble << "\n" ;
            r-=log(pValInDouble/P0);
        }
    }

    //size_t  K = pvals.size();
    //std::cout << " nbr Pval < P0 : " << K << "\n" ;
    //std::sort(pvals.begin(),pvals.end());
    double leftSum = 0 ;
    double rightSum = 0 ;
    for(size_t i =0 ; i < K ; i++) { // LeftSum
        for(size_t j = i+1 ; j < K+1 ; j++) { // RightSum
            rightSum += BinCoeff(M, j) * pow(P0, j) * pow((1.0 - P0), (M - j));
            //std::cout << rightSum << "\n" ;
        }
        leftSum += (pow(r, i) / factorial(i))*rightSum;
        std::cout << i << "^" << r << "/" << factorial(i) << "*" << rightSum << "=" <<  (pow(r, i) / factorial(i))*rightSum << "\n" ;
    }

    //std::cout  << "LeftSum : " << leftSum << "\tRightSum : " << rightSum << "\n" ;

    std::vector<double> ppvals;
    /*for (size_t i=0; i<pvals.size(); i++) {
        //ppvals.push_back(kthOrderProba(i,crossLinks,pvals[i]));
    }*/
    //double pMin = *(std::min_element(ppvals.begin(),ppvals.end()));
    double I=0 ;
    if(r==0){I=1;}
    updatedPval = (1.0-pow((1.0-P0), M))*I + exp(-r)*leftSum ;
    std::cout << "rest equal : " << pow((1.0-P0), K)*I << "\n" ;
    std::cout << "updatedPval = " << "(1-(1-"<<P0<<")^"<< M <<")*"<<I<<"+e"<<-r<<"*"<<leftSum << " = " << updatedPval <<"\n" ;
    std::cout << ".\n.\n" ;
    //updatedPval = kthOrderProba(0,pvals.size(),pMin);
    //std::cout << updatedPval << "\n" ;
    Buffer << params->targetSetKey << "\t" << updatedPval ;//* targetSetSizeDB->getSize() ;
    return Buffer.str() ;
}

