#include <unistd.h>
#include <string>
#include <limits.h>
#include <stdlib.h>
#include <SubstitutionMatrix.h>
#include <Sequence.h>

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

int translateaa(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads);
    writer.open();
    size_t entries = reader.getSize();
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0f, -0.0f);

    char lookupAA[21][3];
    const char nucLookup[4] = {'A','C', 'G','T'};
    char data[3];
    char writeAA[1];
    for(size_t i = 0; i < 20; i++){
        bool found = false;
        for(size_t nuc1 = 0; nuc1 < 4 && found == false; nuc1++){
            for(size_t nuc2 = 0; nuc2 < 4 && found == false; nuc2++){
                for(size_t nuc3 = 0; nuc3 < 4 && found == false; nuc3++) {
                    data[0] = nucLookup[nuc1];
                    data[1] = nucLookup[nuc2];
                    data[2] = nucLookup[nuc3];
                    translateNucl.translate(writeAA, data, 3);
                    if(writeAA[0] == subMat.int2aa[i]){
                        lookupAA[i][0]=data[0];
                        lookupAA[i][1]=data[1];
                        lookupAA[i][2]=data[2];
                        found = true;
                    }
                }
            }
        }
    }

#pragma omp parallel
    {
        char* aa = new char[par.maxSeqLen/3 + 3 + 1];
        std::string nucSeq;
        nucSeq.reserve(10000);
        Sequence centerSequence(par.maxSeqLen, Sequence::AMINO_ACIDS, &subMat, 0, false, par.compBiasCorrection);

#pragma omp for schedule(dynamic, 5)
        for (size_t i = 0; i < entries; ++i) {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif
            unsigned int key = reader.getDbKey(i);
            char* data = reader.getData(i);
            centerSequence.mapSequence(0, key, data);

            //190344_chr1_1_129837240_129837389_3126_JFOA01000125.1 Prochlorococcus sp. scB245a_521M10 contig_244, whole genome shotgun sequence  [Orf: 1, 202, -1, 1, 0]
            // ignore null char at the end
            for(size_t pos = 0; pos < centerSequence.L; pos++){
                nucSeq.append(lookupAA[centerSequence.int_sequence[pos]],3);
            }

            nucSeq.push_back('\n');
//        std::cout << aa << std::endl;
            writer.writeData(nucSeq.c_str(), nucSeq.size(), key, thread_idx);
            nucSeq.clear();
        }
        delete[] aa;
    }
    writer.close(DBReader<unsigned int>::DBTYPE_AA);

    reader.close();

    return EXIT_SUCCESS;
}


