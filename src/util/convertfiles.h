#ifndef MMSEQS_CONVERTFILES_H
#define MMSEQS_CONVERTFILES_H

#include <DBReader.h>

class convertfiles {
public:
    explicit convertfiles(std::string sequencedb, bool use_header);

    void getAlignmentscoresForCluster(std::string clusteringfile, std::string alignmentfile, std::string outputfile);

    void convertDomainFileToFFindex(std::string domainscorefile, std::string domainIdentifierFile,
                                    std::string outputfile);

    void getDomainScoresForCluster(std::string clusteringfile, std::string alignmentfile, std::string outputfolder,
                                   std::string prefix, bool allagainstall, bool randomized);

    void convertFfindexToTsv(std::string clusteringfile, std::string prefix, std::string outputfolder);

private:
    bool use_header;
    DBReader<unsigned int> *targetdb_header;
    std::string getProteinNameForID(unsigned int dbKey);
};


#endif //MMSEQS_CONVERTFILES_H
