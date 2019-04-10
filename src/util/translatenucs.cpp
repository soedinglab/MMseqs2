#include <string>
#include <limits.h>

#include "Orf.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "TranslateNucl.h"
#include "FileUtil.h"

#ifdef OPENMP
#include <omp.h>
#endif

int translatenucs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    bool addOrfStop = par.addOrfStop;
    DBReader<unsigned int> *header = NULL;
    if (addOrfStop == true) {
        header = new DBReader<unsigned int>(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        header->open(DBReader<unsigned int>::NOSORT);
    }

    size_t entries = reader.getSize();

#ifdef OPENMP
    unsigned int totalThreads = par.threads;
#else
    unsigned int totalThreads = 1;
#endif

    unsigned int localThreads = totalThreads;
    if (entries <= totalThreads) {
        localThreads = entries;
    }
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), localThreads, par.compressed, Parameters::DBTYPE_AMINO_ACIDS);
    writer.open();
    Debug::Progress progress(entries);

    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));
#pragma omp parallel num_threads(localThreads)
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        const char * entry[255];
        char* aa = new char[par.maxSeqLen + 3 + 1];
#pragma omp for schedule(dynamic, 5)
        for (size_t i = 0; i < entries; ++i) {
            progress.updateProgress();

            bool addStopAtStart = false;
            bool addStopAtEnd = false;

            unsigned int key = reader.getDbKey(i);
            char* data = reader.getData(i, thread_idx);
            if (addOrfStop == true) {
                char* headData = header->getDataByDBKey(key, thread_idx);
                size_t columns = Util::getWordsOfLine(headData, entry, 255);
                size_t col;
                bool found = false;
                for (col = 0; col < columns; col++) {
                    if (entry[col][0] == '[' && entry[col][1] == 'O' && entry[col][2] == 'r' && entry[col][3] == 'f' && entry[col][4] == ':') {
                        found=true;
                        break;
                    }
                }
                if (found == false) {
                    Debug(Debug::ERROR) << "Could not find Orf information in header.\n";
                    EXIT(EXIT_FAILURE);
                }

                Orf::SequenceLocation loc = Orf::parseOrfHeader(entry[col]);

                addStopAtStart=!(loc.hasIncompleteStart);
                addStopAtEnd=!(loc.hasIncompleteEnd);
            }

            //190344_chr1_1_129837240_129837389_3126_JFOA01000125.1 Prochlorococcus sp. scB245a_521M10 contig_244, whole genome shotgun sequence  [Orf: 1, 202, -1, 1, 0]
            // ignore null char at the end
            // needs to be int in order to be able to check
            size_t length = reader.getSeqLens(i) - 1;
            if ((data[length] != '\n' && length % 3 != 0) && (data[length - 1] == '\n' && (length - 1) % 3 != 0)) {
                Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is not divisible by three. Adjust length to (lenght=" <<  length - (length % 3) << ").\n";
                length = length - (length % 3);
            }

            if (length < 3) {
                Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is too short. Skipping entry.\n";
                continue;
            }

            if (length > (3 * par.maxSeqLen)) {
                Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is too long. Trimming entry.\n";
                length = (3 * par.maxSeqLen);
            }

            char * writeAA;
            if (addStopAtStart) {
                aa[0]='*';
                writeAA = aa + 1;
            } else {
                writeAA = aa;
            }
            translateNucl.translate(writeAA, data, length);

            if (addStopAtEnd && writeAA[(length/3)-1]!='*') {
                writeAA[length/3] = '*';
                writeAA[length/3+1] = '\n';
            } else {
                addStopAtEnd =false;
                writeAA[length/3] = '\n';
            }

            writer.writeData(aa, (length / 3) + 1 + addStopAtStart + addStopAtEnd, key, thread_idx);
        }
        delete[] aa;
    }
    writer.close(true);

    FileUtil::symlinkAbs(par.hdr1, par.hdr2);
    FileUtil::symlinkAbs(par.hdr1Index, par.hdr2Index);
    DBWriter::writeDbtypeFile(par.hdr2.c_str(), Parameters::DBTYPE_GENERIC_DB, par.compressed);

    if (addOrfStop == true) {
        header->close();
    }
    reader.close();

    return EXIT_SUCCESS;
}



