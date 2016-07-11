// Computes either a PSSM or a MSA from clustering or alignment result
// For PSSMs: MMseqs just stores the position specific score in 1 byte

#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>

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
#include "Log.h"

#ifdef OPENMP
#include <omp.h>
#endif



statsComputer::statsComputer(Parameters &par)//:par(par)
{
    this->par = &par;
    if (par.stat == STAT_LINECOUNT_STR)
    {
        stat = STAT_LINECOUNT;
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
        
        statWriter = new DBWriter(par->db4.c_str(), par->db4Index.c_str(), par->threads, DBWriter::BINARY_MODE);
        statWriter->open();
    }
    
    switch (stat)
    {
        case STAT_LINECOUNT:
            return countNumberOfLines();
        default:
            return 0;
        
    }
}

statsComputer::~statsComputer()
{
    if (stat != STAT_UNKNOWN)
    {
        statWriter->close();
        resultReader->close();
        delete statWriter;
        delete resultReader;
    }
}

int statsComputer::countNumberOfLines()
{
    for (size_t id = 0; id < resultReader->getSize(); id++) {
            Log::printProgress(id);
            unsigned int thread_idx = 0;
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
            
            statWriter->write(lineCountString.c_str(), lineCountString.length(), SSTR(resultReader->getDbKey(id)).c_str(),
                                           thread_idx);
    }
    return 0;
    
}


int result2stats(int argc, const char **argv) {
    //MMseqsMPI::init(argc, argv);

    std::string usage("Compute statistics from a result.\n");
    usage.append("USAGE: <queryDB> <targetDB> <resultDB> <outDB>\n");
    usage.append("\nDesigned and implemented by Clovis Galiez <clovis.galiez@mpibpc.mpg.de>\n");
    
    
    Parameters par;
    par.parseParameters(argc, argv, usage, par.result2stats, 4);

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
