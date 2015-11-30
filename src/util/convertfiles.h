//
// Created by lars on 02.06.15.
//

#ifndef MMSEQS_CONVERTFILES_H
#define MMSEQS_CONVERTFILES_H
#include <Util.h>
#include <DBReader.h>

class convertfiles{

public:

    convertfiles(std::string sequencedb);
    void getAlignmentscoresForCluster(std::string clusteringfile, std::string alignmentfile, std::string outputfile);

    void convertDomainFileToFFindex(std::string domainscorefile, std::string domainIdentifierFile, std::string outputfile);

    void getDomainScoresForCluster(std::string clusteringfile, std::string alignmentfile, std::string outputfolder, std::string prefix, bool allagainstall,bool randomized);

    void convertFfindexToTsv(std::string clusteringfile, std::string prefix, std::string outputfolder);

private:
    DBReader* targetdb_header;
    std::string getProteinNameForID(const char * dbKey);

    };


#endif //MMSEQS_CONVERTFILES_H
