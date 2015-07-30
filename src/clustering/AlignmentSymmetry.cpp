//
// Created by lars on 10.06.15.
//

#include "AlignmentSymmetry.h"
#include <Debug.h>
#include "AffinityClustering.h"

AlignmentSymmetry::AlignmentSymmetry(DBReader * seqDbr, DBReader * alnDbr, DBWriter* alnWr, float seqIdThr, float coverage){

    this->seqDbr=seqDbr;
    this->alnDbr=alnDbr;
    this->seqIdThr=seqIdThr;
    this->coverage=coverage;
    this->dbSize=alnDbr->getSize();
    this->alnWr=alnWr;

}

void AlignmentSymmetry::execute() {
    char *similarity = new char[255+1];
    unsigned int dbcutnumber=5;
    unsigned int stepsize=dbSize/dbcutnumber+1;
    unsigned int missingnumber=0;
    unsigned int connections=0;

    size_t BUFFER_SIZE = 1000000;
    char* outBuffer = new char[BUFFER_SIZE];

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
                        std::string stringbuffer=std::string(linebuffer);
                        missinglines[targetid-start].insert(alnDbr->getDbKey(j)+stringbuffer.substr(stringbuffer.find("\t",0)));
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
            std::string cluResultsOutString = std::string("");
            while (*data != '\0') {
		 Util::parseKey(data, idbuffer1);
                Util::getLine(data, linebuffer);
                cluResultsOutString=cluResultsOutString+linebuffer+"\n";
                data = Util::skipLine(data);
            }

                  // Debug(Debug::INFO)<<data;
                if(!missinglines[i-start].empty()){
                    for(std::string line :missinglines[i-start]){
                      //  Debug(Debug::INFO)<<line<<"\n";
                        cluResultsOutString=cluResultsOutString+line+"\n";
                    }
                }

            const char* cluResultsOutData = cluResultsOutString.c_str();
            if (BUFFER_SIZE < strlen(cluResultsOutData)){
                Debug(Debug::ERROR) << "Tried to process the alignment list for the query " << alnDbr->getDbKey(i)
                << " , length of the list = " << "\n";
                Debug(Debug::ERROR) << "Output buffer size < clustering result size! (" << BUFFER_SIZE << " < " << cluResultsOutString.length()
                << ")\nIncrease buffer size or reconsider your parameters -> output buffer is already huge ;-)\n";
                continue;
            }
            memcpy(outBuffer, cluResultsOutData, cluResultsOutString.length()*sizeof(char));
            alnWr->write(outBuffer, cluResultsOutString.length(), alnDbr->getDbKey(i));

              //  data = Util::skipLine(data);
          //  }
        }
        delete[]resultsets;
        delete[]missinglines;
        Debug(Debug::INFO)<<step<<"\t"<<missingnumber<<"\t"<<connections<<"\n";
    }
    Debug(Debug::INFO)<<missingnumber<<"\n";
}



