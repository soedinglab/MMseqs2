#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include <vector>
#include "Parameters.h"
#include <utility>      // std::pair



int mergeffindex (int argc, const char * argv[])
{
    std::string usage("\nMerge multiple ffindex files based on similar id into one file. \n");
    usage.append("Written by Martin Steinegger (martin.steinegger@mpibpc.mpg.de).\n\n");
    usage.append("USAGE: <queryDB> <outDB> <ffindexDB1> ... <ffindexDBn>\n");
    Parameters par;
    par.parseParameters(argc, argv, usage, par.onlyverbosity, 4, true, true);

    std::vector<std::pair<std::string, std::string>> filenames;
    for(int i = 2; i < argc; i++){
        filenames.push_back( std::make_pair<std::string, std::string>(std::string(argv[i]), std::string(argv[i])+".index" ));
    }
    DBReader<unsigned int> qdbr(par.db1.c_str(), std::string(par.db1+".index").c_str());
    qdbr.open(DBReader<unsigned int>::NOSORT);
    DBWriter writer(par.db2.c_str(), std::string( par.db2 +".index").c_str());
    writer.open();
    writer.mergeFiles(qdbr, filenames);
    writer.close();
    qdbr.close();
    return 0;
}
