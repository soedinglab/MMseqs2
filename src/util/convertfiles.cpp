//
// Created by lars on 02.06.15.
//

#include <Debug.h>
#include <DBReader.h>
#include <bits/stl_set.h>
#include <sstream>
#include <bits/stl_list.h>
#include <DBWriter.h>
#include "convertfiles.h"

void convertfiles::convertFfindexToTsv(std::string clusteringfile,std::string suffix,
                                                std::string outputfile) {
    Debug(Debug::INFO) <<"transfering: "<<clusteringfile <<" to "<<outputfile<<suffix<<"cluster.tsv";

    std::string cluster_ffindex_indexfile=clusteringfile+".index";

    DBReader* cluster_ffindex_reader = new DBReader(clusteringfile.c_str(), cluster_ffindex_indexfile.c_str());
    cluster_ffindex_reader->open(DBReader::SORT);

       std::ofstream outfile_stream;
    outfile_stream.open(outputfile+suffix+"cluster.tsv");

    outfile_stream<<"algorithm\tclusterid\tid2\n";
    for (int i = 0; i < cluster_ffindex_reader->getSize(); ++i) {


        char *representative=cluster_ffindex_reader->getDbKey(i);
        char *data = cluster_ffindex_reader->getData(i);
        char *idbuffer = new char[255 + 1];
        while (*data != '\0') {
            Util::parseKey(data, idbuffer);
            outfile_stream<<suffix<<"\t"<<representative<<"\t"<<idbuffer<<"\n";
            data = Util::skipLine(data);
        }
        outfile_stream.flush();
    }
    cluster_ffindex_reader->close();
    cluster_ffindex_reader->~DBReader();
    outfile_stream.close();

}


void convertfiles::getAlignmentscoresForCluster(std::string clusteringfile, std::string alignmentfile,
                                                std::string outputfile) {
    //Debug(Debug::INFO) <<clusteringfile <<alignmentfile<<outputfile;

    std::string cluster_ffindex_indexfile=clusteringfile+".index";
    std::string alignment_ffindex_indexfile=alignmentfile+".index";
    DBReader* cluster_ffindex_reader = new DBReader(clusteringfile.c_str(), cluster_ffindex_indexfile.c_str());
    cluster_ffindex_reader->open(DBReader::SORT);
    DBReader* alignment_ffindex_reader = new DBReader(alignmentfile.c_str(), alignment_ffindex_indexfile.c_str());
    alignment_ffindex_reader->open(DBReader::SORT);

    std::ofstream outfile_stream;
    outfile_stream.open(outputfile);

    outfile_stream<<"clusterid\tid2\taliscore\tqcov\tdbcov\tseqId\teval\n";
    for (int i = 0; i < cluster_ffindex_reader->getSize(); ++i) {


        char *representative=cluster_ffindex_reader->getDbKey(i);
        char *data = cluster_ffindex_reader->getData(i);
        char *idbuffer = new char[255 + 1];
        char *linebuffer=new char[255+1];


        std::set<std::string> clusterset;




        while (*data != '\0') {
            Util::parseKey(data, idbuffer);
                clusterset.insert(std::string(idbuffer));
            data = Util::skipLine(data);
        }

        char *data_alignment = alignment_ffindex_reader->getDataByDBKey(representative);
        if (data_alignment== NULL) {
          //  Debug(Debug::INFO) <<representative<<"\n";
            continue;
        }
        while (*data_alignment != '\0') {
            Util::parseKey(data_alignment, idbuffer);
            //Debug(Debug::INFO) <<idbuffer;
            if(clusterset.find(idbuffer)!= clusterset.end()){
                outfile_stream<<representative<<"\t"<<Util::getLine(data_alignment,linebuffer)<<"\n";
            }
            data_alignment = Util::skipLine(data_alignment);
        }
        outfile_stream.flush();

    }




    alignment_ffindex_reader->close();
    alignment_ffindex_reader->~DBReader();
    cluster_ffindex_reader->close();
    cluster_ffindex_reader->~DBReader();
    outfile_stream.close();

}

void convertfiles::convertDomainFileToFFindex(std::string domainscorefile, std::string domainIdentifierFile,
                                              std::string outputfile) {
 /*
  * Read in Domain identifiers
  * Line number (1 based index) is the index of the identifier in the domainfile
  */
    std::ifstream domainIdentifierFile_stream;
    domainIdentifierFile_stream.open(domainIdentifierFile);
    std::list<std::string> domainIdentifiers;
    std::string line;
    while (std::getline(domainIdentifierFile_stream, line))
    {

        domainIdentifiers.push_back(std::string(line));

    }
    int number=1;
    std::string *domainIdentifier=new std::string[domainIdentifiers.size()+1];
    for(std::string id:domainIdentifiers){
        domainIdentifier[number] =id;
       // Debug(Debug::INFO)<<id<<"\n";
        number++;
    }
    /*
     * read in scorefile , replace identifier and transfer to ffindex format
     */
    std::ifstream domainscorefile_stream;
    domainscorefile_stream.open(domainscorefile);
    size_t BUFFER_SIZE = 1000000;
    char* outBuffer = new char[BUFFER_SIZE];
    std::string outputfileindex =outputfile+".index";
    DBWriter *dbw = new DBWriter(outputfile.c_str(), outputfileindex.c_str());
    dbw->open();
    std::stringstream res;
    int lastid=-1;
    char * clusterid=new char[20];
    int a, b;
    double c;
    while (std::getline(domainscorefile_stream, line))
    {

        std::istringstream iss(line);

        if (!(iss >> a >> b>> c)) {
            Debug(Debug::INFO)<<line<<"\n";
            break; }
        res<<domainIdentifier[b]<<"\t"<<c<<"\n";

       // domainIdentifier[a]<<"\t"<<domainIdentifier[b]<<"\t"<<res.str().length()<<"\n";

        if(lastid !=-1 && a !=lastid) {
            std::string cluResultsOutString = res.str();
            const char* cluResultsOutData = cluResultsOutString.c_str();
            if (BUFFER_SIZE < strlen(cluResultsOutData)){
                Debug(Debug::ERROR) << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " << cluResultsOutString.length()
                << ")\nIncrease bclusterset.erase(idbuffer);uffer size or reconsider your parameters -> output buffer is already huge ;-)\n";
                continue;
            }


            memcpy(clusterid, domainIdentifier[lastid].c_str(), domainIdentifier[lastid].length() * sizeof(char));
            clusterid[domainIdentifier[lastid].length()]='\0';
            memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length() * sizeof(char));
            dbw->write(outBuffer, cluResultsOutString.length(), clusterid);
            //clear res
            res.str("");
            res.clear();
        }
            lastid=a;
    }
    std::string cluResultsOutString = res.str();
    const char* cluResultsOutData = cluResultsOutString.c_str();
    if (BUFFER_SIZE < strlen(cluResultsOutData)){
        Debug(Debug::ERROR) << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " << cluResultsOutString.length()
        << ")\nIncrease buffer size or reconsider your parameters -> output buffer is already huge ;-)\n";

    }

    //save last entry
    memcpy(clusterid, domainIdentifier[a].c_str(), domainIdentifier[a].length() * sizeof(char));
    clusterid[domainIdentifier[a].length()]='\0';
    memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length() * sizeof(char));
    dbw->write(outBuffer, cluResultsOutString.length(), clusterid);
    //clear res
    res.str("");
    res.clear();

    delete[] outBuffer;


    dbw->close();
    dbw->~DBWriter();
    domainIdentifierFile_stream.close();

}

void convertfiles::getDomainScoresForCluster(std::string clusteringfile, std::string alignmentfile,
                                                std::string outputfolder, std::string prefix, bool allagainstall) {
   // Debug(Debug::INFO) <<clusteringfile <<alignmentfile<< outputfolder;

    std::string cluster_ffindex_indexfile=clusteringfile+".index";
    std::string alignment_ffindex_indexfile=alignmentfile+".index";
    DBReader* cluster_ffindex_reader = new DBReader(clusteringfile.c_str(), cluster_ffindex_indexfile.c_str());
    cluster_ffindex_reader->open(DBReader::SORT);
    DBReader* alignment_ffindex_reader = new DBReader(alignmentfile.c_str(), alignment_ffindex_indexfile.c_str());
    alignment_ffindex_reader->open(DBReader::SORT);

    std::ofstream outfile_stream;
    outfile_stream.open(outputfolder +"/"+ prefix +"domainscore.tsv");

    outfile_stream<<"algorithm\tclusterid\tid2\tdomain_score\n";
    for (int i = 0; i < cluster_ffindex_reader->getSize(); ++i) {


        char *representative=cluster_ffindex_reader->getDbKey(i);
        char *data = cluster_ffindex_reader->getData(i);
        char *idbuffer = new char[255 + 1];
        char *linebuffer=new char[255+1];


        std::set<std::string> clusterset;
        std::set<std::string> clusterset2;



        while (*data != '\0') {
            Util::parseKey(data, idbuffer);
            clusterset.insert(std::string(idbuffer));
            data = Util::skipLine(data);
        }

        char *data_alignment = alignment_ffindex_reader->getDataByDBKey(representative);
        if (data_alignment== NULL) {
            //  Debug(Debug::INFO) <<representative<<"\n";
            continue;
        }
        Debug(Debug::INFO) <<representative<<"\n";
        while (*data_alignment != '\0') {
            Util::parseKey(data_alignment, idbuffer);
           // Debug(Debug::INFO) <<idbuffer;
            if(clusterset.find(idbuffer)!= clusterset.end()){
                clusterset2.insert(std::string(idbuffer));
                outfile_stream<< prefix <<"\t"<<representative<<"\t"<<Util::getLine(data_alignment,linebuffer)<<"\n";
            }
            data_alignment = Util::skipLine(data_alignment);
        }
        outfile_stream.flush();



        for(std::string id :clusterset){
            strncpy(idbuffer, id.c_str(), id.length());
            idbuffer[id.length()] = '\0';
            if(strcmp(representative,id.c_str())==0||clusterset2.find(id)!= clusterset2.end()){

            }else{
                if(alignment_ffindex_reader->getDataByDBKey(idbuffer)!=NULL){
                    outfile_stream<< prefix <<"\t"<<representative<<"\t"<<id<<"\t"<<"0"<<"\n";
                }
            }

        }


        outfile_stream.flush();

    }




    alignment_ffindex_reader->close();
    alignment_ffindex_reader->~DBReader();
    cluster_ffindex_reader->close();
    cluster_ffindex_reader->~DBReader();
    outfile_stream.close();

}