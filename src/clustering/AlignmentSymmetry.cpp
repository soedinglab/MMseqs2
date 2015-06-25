//
// Created by lars on 10.06.15.
//

#include "AlignmentSymmetry.h"
#include <Debug.h>
#include "AffinityClustering.h"

AlignmentSymmetry::AlignmentSymmetry(DBReader * seqDbr, DBReader * alnDbr, float seqIdThr, float coverage){

    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->seqIdThr=seqIdThr;
    this->coverage=coverage;
    this->dbSize=alnDbr->getSize();

}

std::list<set *>  AlignmentSymmetry::execute() {
    unsigned int dbcutnumber=5;
    unsigned int stepsize=dbSize/dbcutnumber+1;
    unsigned int missingnumber=0;
    unsigned int connections=0;

    for (int step = 0; step < dbcutnumber; ++step) {
        std::set<int>*resultsets=new std::set<int>[stepsize];
        std::set<std::string>*missinglines=new std::set<std::string>[stepsize];
        unsigned int start= step *stepsize;
        unsigned int end=std::min((step+1) *stepsize,dbSize);
      //  Debug(Debug::INFO)<<start<<"\t"<<end<<"\n";
//initialise datastructures
        char *idbuffer1 = new char[255 + 1];
        char *linebuffer = new char[255 + 1];
        memset(resultsets,0,sizeof(resultsets));
        memset(missinglines,0,sizeof(missinglines));
//read in datachunk
        for (int i = start; i < end; i++) {
        char *data = alnDbr->getData(i);
        while (*data != '\0') {
            Util::parseKey(data, idbuffer1);
            resultsets[i-start].insert(alnDbr->getId(idbuffer1));
            data = Util::skipLine(data);
            }
        }

        Debug(Debug::INFO)<<"reading data finished:\t"<<step<<"\n";

        //determine missing connections
        for (int j = 0; j < dbSize ; ++j) {
            char *data = alnDbr->getData(j);
            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);
                int targetid=alnDbr->getId(idbuffer1);
                if(targetid>=start && targetid<end){
                    if(resultsets[targetid-start].find(j) ==resultsets[targetid-start].end()){
                        missingnumber++;
                        Util::getLine(data,linebuffer);
                    //    missinglines[targetid-start].insert(std::string(linebuffer));
                    }else{
                        connections++;
                    }

                }
                data = Util::skipLine(data);
            }
        }
        //print
        for (int i = start; i < end; i++) {
            char *data = alnDbr->getData(i);
            while (*data != '\0') {
                Util::parseKey(data, idbuffer1);

                data = Util::skipLine(data);
            }
        }
        delete[]resultsets;
        delete[]missinglines;
        Debug(Debug::INFO)<<step<<"\t"<<missingnumber<<"\t"<<connections<<"\n";
    }
    Debug(Debug::INFO)<<missingnumber<<"\n";
}



