#include <unistd.h>
#include <string>
#include <limits.h>
#include <stdlib.h>

#include "Orf.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "TranslateNucl.h"

#ifdef OPENMP
#include <omp.h>
#endif

int translatenucs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    std::string in_header_filename = std::string(par.db1 + "_h");
    std::string in_header_index_filename = std::string(par.db1 + "_h.index");
    std::string out_header_filename  = std::string(par.db2 + "_h");
    std::string out_header_index_filename = std::string(par.db2 + "_h.index");

    // set links to header
    Debug(Debug::INFO) << "Set sym link from " << in_header_filename << " to " << out_header_filename << "\n";
    char *abs_in_header_filename = realpath(in_header_filename.c_str(), NULL);
    symlink(abs_in_header_filename, out_header_filename.c_str());
    free(abs_in_header_filename);
    char *abs_in_header_index_filename = realpath(in_header_index_filename.c_str(), NULL);
    Debug(Debug::INFO) << "Set sym link from " << in_header_index_filename << " to " << out_header_index_filename << "\n";
    symlink(abs_in_header_index_filename, out_header_index_filename.c_str());
    free(abs_in_header_index_filename);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::NOSORT);

    bool addOrfStop = par.addOrfStop;
    DBReader<unsigned int> * header;
    if(addOrfStop == true) {
        header=new DBReader<unsigned int>((par.db1 + "_h").c_str(), (par.db1 + "_h.index").c_str());
        header->open(DBReader<unsigned int>::NOSORT);
    }
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    writer.open();
    size_t entries = reader.getSize();
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
#pragma omp parallel
    {
        char* aa = new char[par.maxSeqLen/3 + 3 + 1];
#pragma omp  for schedule(static)
        for (size_t i = 0; i < entries; ++i) {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            unsigned int key = reader.getDbKey(i);
            char* data = reader.getData(i);
            bool addStopAtStart = false;
            bool addStopAtEnd = false;
            if(addOrfStop == true){
                char* headData = header->getData(i);
                char * entry[255];
                size_t columns = Util::getWordsOfLine(headData, entry, 255);
                size_t col;
                bool found = false;
                for(col = 0; col < columns; col++){
                    if(entry[col][0] == '[' && entry[col][1] == 'O' && entry[col][2] == 'r' && entry[col][3] == 'f' & entry[col][4] == ':'){
                        found=true;
                        break;
                    }
                }
                if(found==false ){
                    Debug(Debug::ERROR) << "Could not find Orf information in header.\n";
                    EXIT(EXIT_FAILURE);
                }
                Orf::SequenceLocation loc;
                int strand;
                int retCode = sscanf(entry[col], "[Orf: %zu, %zu, %d, %d, %d]", &loc.from, &loc.to, &strand, &loc.hasIncompleteStart, &loc.hasIncompleteEnd);
                if(retCode < 5){
                    Debug(Debug::ERROR) << "Could not parse Orf " << entry[col] << ".\n";
                    EXIT(EXIT_FAILURE);
                }
                loc.strand =  static_cast<Orf::Strand>(strand);
                addStopAtStart=!(loc.hasIncompleteStart);
                addStopAtEnd=!(loc.hasIncompleteEnd);
            }

            //190344_chr1_1_129837240_129837389_3126_JFOA01000125.1 Prochlorococcus sp. scB245a_521M10 contig_244, whole genome shotgun sequence  [Orf: 1, 202, -1, 1, 0]
            // ignore null char at the end
            // needs to be int in order to be able to check
            int length = reader.getSeqLens(i) - 1;

            if((data[length] != '\n' && length % 3 != 0) && (data[length - 1] == '\n' && (length - 1) % 3 != 0)) {
                Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is not divisible by three. Adjust length to (lenght=" <<  length - (length % 3) << ").\n";
                length = length - (length % 3);
            }

            if(length < 3)  {
                Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is too short. Skipping entry.\n";
                continue;
            }
//        std::cout << data << std::endl;
            char * writeAA;
            if(addStopAtStart){
                aa[0]='*';
                writeAA = aa + 1;
            }else{
                writeAA = aa;
            }
            translateNucl.translate(writeAA, data, length);

            if(addStopAtEnd){
                writeAA[length/3] = '*';
                writeAA[length/3+1] = '\n';
            }else{
                writeAA[length/3] = '\n';
            }

//        std::cout << aa << std::endl;
            writer.writeData(aa, (length / 3) + 1 + addStopAtStart + addStopAtEnd, key, thread_idx);
        }
        delete[] aa;
    }
    writer.close(DBReader<unsigned int>::DBTYPE_AA);
    if(addOrfStop == true) {
        header->close();
    }
    reader.close();

    return EXIT_SUCCESS;
}


