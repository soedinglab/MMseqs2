#include "Aggregation.h"
#include "Util.h"
#include "Parameters.h"
#include "Debug.h"

#include <algorithm>
#include <tgmath.h>
#include <cfloat>
#include <sstream>

#ifdef OPENMP
#include <omp.h>
#endif



struct compareByStart {
    bool operator()(const std::pair <long,long> &lhs,
                    const std::pair <long,long> &rhs) const {
        return (lhs.first < rhs.first);
    }
};

// TO DO : Get rid of stringStreams
Aggregation::Aggregation(std::string arg_inputDBname, std::string arg_outputDBname, unsigned int arg_nbrThread,
                         size_t arg_targetColumn) : inputDBname(std::move(arg_inputDBname)), outputDBname(std::move(arg_outputDBname)),
                                                    nbrThread(arg_nbrThread), targetColumn(arg_targetColumn){}

// build a map with the value in [target column] field as a key and the rest of the line, cut in fields, as values
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

    unsigned int thread_idx = 0;
#pragma omp parallel for private(thread_idx)
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
        std::map<unsigned int, std::vector<std::vector<std::string> > > dataToMerge ;

        if (this->buildMap(geneResultsStream,dataToMerge)){
            std::string outputBuffer;
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


BestHitAggregator::BestHitAggregator(std::string arg_inputDBname, std::string arg_outputDBname,std::string arg_targetSetSizeName,
                                     size_t arg_targetColumn, unsigned int arg_nbrThread,  size_t arg_scoreColumn) :
                                     Aggregation(std::move(arg_inputDBname), std::move(arg_outputDBname), arg_nbrThread,
                                     arg_targetColumn), pValColumn(arg_scoreColumn), targetSetSizeName(std::move(arg_targetSetSizeName)){
    std::string sizeDBIndex = targetSetSizeName + ".index";
    this->targetSetSizeDB = new DBReader<unsigned int> (targetSetSizeName.c_str(), sizeDBIndex.c_str()) ;
    this->targetSetSizeDB->open(DBReader<unsigned int>::NOSORT);
}

// Only Keep the best hits of a protein against each Target Set
std::string BestHitAggregator::aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, aggregParams* vparams) {

    std::string buffer ;
    double bestEval = DBL_MAX;
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
    unsigned int nbrGenes = (unsigned int)std::stoul(this->targetSetSizeDB->getDataByDBKey(
            static_cast<unsigned int>(std::stoul(dataToAggregate[maxId][targetColumn])))) ;

    if(dataToAggregate.size() < 2) {
        secondBestEval = (2.00*nbrGenes/(nbrGenes+1.00)) ;
    }

    correctedPval = exp(log(bestEval / nbrGenes) - log(secondBestEval / nbrGenes));
    //std::cout << "maxVal : " << bestEval << " \t secondMaxVal :" << secondBestEval << "\t Number of entries ; "
              //<< dataToAggregate.size() << " corrected:"<<correctedPval << "nbr Genes : " << "\n" ;

    char tmpBuf[15];
    sprintf(tmpBuf,"%.3E",correctedPval);
    dataToAggregate[maxId][this->pValColumn] = std::string(tmpBuf);

    // Aggregate the full line in one string
    for( std::vector<std::string>::iterator it=dataToAggregate[maxId].begin() ; it!= dataToAggregate[maxId].end() ; ++it){
        buffer.append(*it + "\t");
    }
    return buffer ;
}


PvalAggregator::PvalAggregator(std::string arg_inputDBname, std::string arg_outputDBname, unsigned int arg_nbrThread,
                               std::string arg_querySetSizeDBname, size_t arg_targetColumn, size_t arg_scoreColumn) :
        Aggregation(std::move(arg_inputDBname), std::move(arg_outputDBname), arg_nbrThread, arg_targetColumn), pValColumn(arg_scoreColumn),
        querySetSizeDBname(std::move(arg_querySetSizeDBname)){

    std::string sizeDBIndex = querySetSizeDBname + ".index";
    this->querySetSizeDB = new DBReader<unsigned int> (querySetSizeDBname.c_str(), sizeDBIndex.c_str()) ;
    this->querySetSizeDB->open(DBReader<unsigned int>::NOSORT);
}

/*double kthOrderProba(size_t k, size_t N, double p) {
    return (double)0.5*std::erfc(((double)k+0.5-N*p)/sqrt(N*p*(1.0-p)*2.0)) ;
}*/

double BinCoeff (int M, int k) {
    double result;
    result = exp(lgamma(M+1)-lgamma(M-k+1)-lgamma(k+1)) ;
    return result;
}
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
    double P0 = 0.001 ;
    double updatedPval;
    double r = 0 ;
    double pValInDouble ;
    size_t k =0 ;
    std::vector<double> pvals;
    for (auto &i : dataToAggregate) {
        pValInDouble = std::strtod(i[pValColumn].c_str(), nullptr) ;
        if (pValInDouble < P0){
            k++ ;
            r-=log(pValInDouble/P0);
        }
    }

    double leftSum = 0 ;
    double rightSum = 0 ;
    for(size_t i =0 ; i < k ; i++) { // LeftSum
        for(size_t j = i+1 ; j < k+1 ; j++) { // RightSum
            rightSum += exp(LBinCoeff(M, j) + j*log(P0) + (M - j)*log(1.0 - P0));//BinCoeff(M, j) * pow(P0, j) * pow((1.0 - P0), (M - j));
        }
        leftSum += (pow(r, i) / factorial(i))*rightSum;
    }

    std::vector<double> ppvals;
    double I=0 ;
    if(r==0){I=1;}
    updatedPval = (1.0-pow((1.0-P0), M))*I + exp(-r)*leftSum ;
    Buffer << params->targetSetKey << "\t" << updatedPval ;//* targetSetSizeDB->getSize() ;
    return Buffer.str() ;
}

ClusteringAggregator::ClusteringAggregator(std::string arg_inputDBname, std::string arg_outputDBname,
                                           unsigned int arg_nbrThread, size_t arg_targetColumn)
        : Aggregation(std::move(arg_inputDBname), std::move(arg_outputDBname), arg_nbrThread, arg_targetColumn) {}

// return the median of the distance between genes of dataToAggregate
std::string ClusteringAggregator::aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate,
                                                 aggregParams *params) {
    int i = 0 ;
    unsigned int indexOfInterSpaceToReturn =0 ;
    std::vector<std::pair<long,long> > genesPositions;
    std::vector<long> interGeneSpaces;
    long currentInterGenePosition = 0;
    std::stringstream buffer ;
    for (auto &it : dataToAggregate) {
        genesPositions.emplace_back(std::make_pair(std::stol(it[8]),std::stol(it[10]))) ;
        i++ ;
    };
    std::sort(begin(genesPositions), end(genesPositions), compareByStart()) ;
    if(i>2) {
        for (size_t pos = 0; pos < i-1; pos++) {
            currentInterGenePosition = genesPositions[pos+1].second - genesPositions[pos].first;
            interGeneSpaces.push_back(currentInterGenePosition);
        }
        std::sort(begin(interGeneSpaces), end(interGeneSpaces)) ;
        //if odd number
        if(interGeneSpaces.size() % 2) {
            indexOfInterSpaceToReturn = static_cast<unsigned int>((interGeneSpaces.size() + 1) / 2);
        }
        else{
            indexOfInterSpaceToReturn = static_cast<unsigned int>(interGeneSpaces.size() / 2);
        }
        if(interGeneSpaces[indexOfInterSpaceToReturn] > 1e9) {
            std::cout << indexOfInterSpaceToReturn << " : " << interGeneSpaces[indexOfInterSpaceToReturn] << "\n";
            for (auto it : interGeneSpaces) { std::cout << it << "\n"; }
        }
        buffer << params->targetSetKey << "\t" << interGeneSpaces[indexOfInterSpaceToReturn] << "\t" << dataToAggregate.size()  ;
    }
    else {
        buffer << params->targetSetKey << "\t0" << "\t" << dataToAggregate.size(); //NaMM
    }
    std::string bufferString= buffer.str() ;
    return bufferString ;
}













































