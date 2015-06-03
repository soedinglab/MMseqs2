//
// Created by lars on 02.06.15.
//

#ifndef MMSEQS_CONVERTFILES_H
#define MMSEQS_CONVERTFILES_H
#include <Util.h>

class convertfiles{

public:
    void getAlignmentscoresForCluster(std::string clusteringfile, std::string alignmentfile, std::string outputfile);

    void convertDomainFileToFFindex(std::string domainscorefile, std::string domainIdentifierFile, std::string outputfile);
};


#endif //MMSEQS_CONVERTFILES_H
