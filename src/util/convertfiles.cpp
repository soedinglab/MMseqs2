//
// Created by lars on 02.06.15.
//

#include <Debug.h>
#include <DBReader.h>
#include <bits/stl_set.h>
#include "convertfiles.h"


void convertfiles::getAlignmentscoresForCluster(std::string clusteringfile, std::string alignmentfile,
                                                std::string outputfile) {
    Debug(Debug::INFO) <<clusteringfile <<alignmentfile<<outputfile;

    std::string cluster_ffindex_indexfile=clusteringfile+".index";
    std::string alignment_ffindex_indexfile=alignmentfile+".index";
    DBReader* cluster_ffindex_reader = new DBReader(clusteringfile.c_str(), cluster_ffindex_indexfile.c_str());
    cluster_ffindex_reader->open(DBReader::SORT);
    DBReader* alignment_ffindex_reader = new DBReader(alignmentfile.c_str(), alignment_ffindex_indexfile.c_str());
    alignment_ffindex_reader->open(DBReader::SORT);

    std::ofstream outfile_stream;
    outfile_stream.open(outputfile);

    outfile_stream<<"clusterid\tid2\taliscore\tqcov\tdbcov\tseqId\teval";
    for (int i = 0; i < cluster_ffindex_reader->getSize(); ++i) {


        char *representative=cluster_ffindex_reader->getDbKey(i);
        char *data = cluster_ffindex_reader->getData(i);
        char *idbuffer = new char[255 + 1];
        char *linebuffer=new char[255+1];
        double sumofscore =0;
        double minscore=1;
        double maxscore=0;

        std::set<std::string> clusterset;




        while (*data != '\0') {
            Util::parseKey(data, idbuffer);
                clusterset.insert(std::string(idbuffer));
            data = Util::skipLine(data);
        }
        char *data_alignment = alignment_ffindex_reader->getDataByDBKey(representative);
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