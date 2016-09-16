// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>
#include <unordered_map>
#include <cmath>

#include "Alignment.h"
#include "MsaFilter.h"
#include "Parameters.h"
#include "PSSMCalculator.h"
#include "DBReader.h"
#include "DBConcat.h"
#include "HeaderSummarizer.h"
#include "result2stats.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif



statsComputer::statsComputer(Parameters &par)//:par(par)
{
    this->par = &par;
    if (par.stat == STAT_LINECOUNT_STR)
    {
        stat = STAT_LINECOUNT;
    } else if (par.stat == STAT_MEAN_STR)
    {
        stat = STAT_MEAN;
    } else if (par.stat == STAT_DOOLITTLE_STR)
    {
        doolittleValues['a'] = 6.3;
        doolittleValues['r'] = 0.0;
        doolittleValues['n'] = 1.0;
        doolittleValues['d'] = 1.0;
        doolittleValues['c'] = 7.0;
        doolittleValues['q'] = 1.0;
        doolittleValues['e'] = 1.0;
        doolittleValues['g'] = 4.1;
        doolittleValues['h'] = 1.3;
        doolittleValues['i'] = 9.0;
        doolittleValues['l'] = 5.2;
        doolittleValues['k'] = 0.6;
        doolittleValues['m'] = 6.4;
        doolittleValues['f'] = 7.2;
        doolittleValues['p'] = 2.9;
        doolittleValues['s'] = 3.6;
        doolittleValues['t'] = 3.8;
        doolittleValues['w'] = 3.6;
        doolittleValues['y'] = 3.2;
        doolittleValues['v'] = 8.7;
        doolittleValues['x'] = 0.0;
        doolittleValues['0'] = 0.0; // N-ter
        doolittleValues['1'] = 0.0; // C-ter

        stat = STAT_DOOLITTLE;
    } else if (par.stat == STAT_CHARGES_STR)
    {
        std::unordered_map<char,float> pKs,chargeSign;
        
        pH = 7.0;
        
        // pKs values:
        // Bjellqvist Dawson EMBOSS Lehninger Murray Rodwell Sillero Solomon Stryer
        pKs['c'] = 9.00;//    8.3    8.5      8.18   8.33    8.33     9.0     8.3    8.5
        pKs['d'] = 4.05;//    3.9    3.9      3.65   3.68    3.86     4.0     3.9    4.4
        pKs['e'] = 4.45;//    4.3    4.1      4.25   4.25    4.25     4.5     4.3    4.4
        pKs['h'] = 5.98;//    6.0    6.5      6.00   6.00    6.00     6.4     6.0    6.5
        pKs['k'] = 10.00;//   10.5   10.8     10.53  11.50   11.50    10.4    10.5   10.0
        pKs['r'] = 12.00;//   12.0   12.5     12.48  11.50   11.50    12.0    12.5   12.0
        pKs['y'] = 10.00;//   10.1   10.1     10.07  10.07   10.70    10.0    10.1   10.0
        pKs['1'] = 3.55;//    3.2    3.6      2.34   2.15    3.10     3.2     3.2    3.2 // C ter
        pKs['0'] = 7.50;//    8.2    8.6      9.69   9.52    8.00     8.2     8.2    8.2 // N ter
        
        chargeSign['c'] = -1.0;
        chargeSign['d'] = -1.0;
        chargeSign['e'] = -1.0;
        chargeSign['y'] = -1.0;
        chargeSign['h'] = 1.0;
        chargeSign['k'] = 1.0;
        chargeSign['r'] = 1.0;
        chargeSign['1'] = -1.0; // C ter
        chargeSign['0'] = 1.0; // N ter

        for (std::unordered_map<char,float>::iterator k = pKs.begin(); k != pKs.end();k++)
                chargeValues[k->first] = chargeSign[k->first] /(1+pow(10,(chargeSign[k->first] * (pH - pKs[k->first]))));

        stat = STAT_CHARGES;
    }  else if (par.stat == STAT_SEQLEN_STR)
    {
        stat = STAT_SEQLEN;
    } else {
        stat = STAT_UNKNOWN;
        Debug(Debug::WARNING) << "Unrecognized statistics: " << par.stat << "\n";
    }
    

}


int statsComputer::run(){
    
    
    if (stat != STAT_UNKNOWN)
    {
        resultReader = new DBReader<unsigned int>(par->db3.c_str(), par->db3Index.c_str());
        resultReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);

        queryReader = new DBReader<unsigned int>(par->db1.c_str(), par->db1Index.c_str());
        queryReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);
        
        targetReader = new DBReader<unsigned int>(par->db2.c_str(), par->db2Index.c_str());
        targetReader->open(DBReader<unsigned int>::LINEAR_ACCCESS);

        statWriter = new DBWriter(par->db4.c_str(), par->db4Index.c_str(), par->threads, DBWriter::BINARY_MODE);
        statWriter->open();
    }
    
    switch (stat)
    {
        case STAT_LINECOUNT:
            return countNumberOfLines();
        case STAT_MEAN:
            return meanValue();
        case STAT_DOOLITTLE:
            return sequenceWise(&statsComputer::doolittle);
        case STAT_CHARGES:
            return sequenceWise(&statsComputer::charges);
        case STAT_SEQLEN:
            return sequenceWise(&statsComputer::strlen);
        default:
            return 0;
        
    }
}

statsComputer::~statsComputer()
{
    if (stat != STAT_UNKNOWN)
    {
        statWriter->close();
        queryReader->close();
        targetReader->close();
        resultReader->close();
        delete statWriter;
        delete queryReader;
        delete targetReader;
        delete resultReader;
    }

}

int statsComputer::countNumberOfLines()
{
    #pragma omp parallel for schedule(static)
    for (size_t id = 0; id < resultReader->getSize(); id++) {
            Debug::printProgress(id);
            unsigned int thread_idx = 0;
            #ifdef OPENMP
            thread_idx = omp_get_thread_num();
            #endif

            unsigned int lineCount(0);
            std::string lineCountString;
            
            char *results = resultReader->getData(id);
            while(*results!='\0')
            {
                if (*results == '\n')
                    lineCount++;
                results++;
            }
            
            lineCountString = std::to_string(lineCount) + '\n';

        statWriter->writeData(lineCountString.c_str(), lineCountString.length(),
                              SSTR(resultReader->getDbKey(id)).c_str(),
                              thread_idx);
    }
    return 0;
    
}


int statsComputer::meanValue()
{
    #pragma omp parallel for schedule(static)
    for (size_t id = 0; id < resultReader->getSize(); id++) {
            Debug::printProgress(id);
            unsigned int thread_idx = 0;
            #ifdef OPENMP
            thread_idx = omp_get_thread_num();
            #endif
            float meanVal(0.0);
            std::string meanValString;
            
            char *results = resultReader->getData(id);
            size_t dataLength = resultReader->getSeqLens(id);
            size_t nbSeq = 0;
            while(*results!='\0')
            {
                meanVal += atof(results);
                nbSeq++;
                results = Util::skipLine(results);
            }
            
            if (!nbSeq)
                nbSeq = 1; // to have a mean value equal to 0
                
                
            meanValString = std::to_string(meanVal/nbSeq) + '\n';

            statWriter->writeData(meanValString.c_str(), meanValString.length(), SSTR(resultReader->getDbKey(id)).c_str(),
                              thread_idx);
    }
    return 0;
    
}

float statsComputer::strlen(char *seq) {
    return (float)std::strlen(seq);
}

float statsComputer::doolittle(char *seq) {
    return averageValueOnAminoAcids(doolittleValues,seq);
}


float statsComputer::charges(char *seq) {
    return averageValueOnAminoAcids(chargeValues,seq);
}

float statsComputer::averageValueOnAminoAcids(std::unordered_map<char,float> values,char* seq)
{
    char *seqPointer = seq;
    float ret = values['0'] + values['1']; // C ter and N ter values
    std::unordered_map<char,float>::iterator k;
    
    while (*seqPointer != '\0' && *seqPointer != '\n') {
        if ((k = values.find(tolower(*seqPointer))) != values.end())
            ret += k->second;
            
        seqPointer++;
    }
    
    size_t seqLen = seqPointer - seq;
    
    if (!seqLen)
        seqLen = 1;
        
    return ret/seqLen;
}

int statsComputer::sequenceWise(float (statsComputer::*statFunction)(char*))
{
    #pragma omp parallel for schedule(static)
    for (size_t id = 0; id < resultReader->getSize(); id++) {
            Debug::printProgress(id);
            unsigned int thread_idx = 0;
            #ifdef OPENMP
            thread_idx = omp_get_thread_num();
            #endif
            char dbKey[255 + 1];

            
            std::string statString;
            char *results = resultReader->getData(id);
            while(*results!='\0') // for every hit
            {
                Util::parseKey(results, dbKey);
                const unsigned int key = (unsigned int) strtoul(dbKey, NULL, 10);
                
                const size_t edgeId = targetReader->getId(key);
                char *dbSeqData = targetReader->getData(edgeId);
                float stat = (this->*statFunction)(dbSeqData);
                
                statString += std::to_string(stat) + '\n';
                
                results = Util::skipLine(results);
            }


            statWriter->writeData(statString.c_str(), statString.length(), SSTR(resultReader->getDbKey(id)).c_str(),
                              thread_idx);
    }
    
    return 0;
}

int result2stats(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    //MMseqsMPI::init(argc, argv);

    struct timeval start, end;
    gettimeofday(&start, NULL);
    int retCode;

    statsComputer computeStats(par);

    retCode = computeStats.run();

    gettimeofday(&end, NULL);
    time_t sec = end.tv_sec - start.tv_sec;
    Debug(Debug::WARNING) << "Time for processing: " << (sec / 3600) << " h " << (sec % 3600 / 60) << " m " << (sec % 60) << "s\n";

    /*
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    */

    return retCode;
}
