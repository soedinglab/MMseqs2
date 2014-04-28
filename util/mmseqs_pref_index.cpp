#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Prefiltering.h"

#include <vector>
#include <utility>      // std::pair

void printUsage(){
    std::string usage("\nMerge multiple ffindex files based on simular id into one file. \n");
    usage.append("Written by Martin Steinegger (Martin.Steinegger@campus.lmu.de) & Maria Hauser (mhauser@genzentrum.lmu.de).\n\n");
    usage.append("USAGE: ffindex_database_merge ffindexQueryDB ffindexOutputDB ffindexOutDB ffindexFILES*\n");
    usage.append("-s\t[int]\tcreate splitted index for n nodes\n");
    usage.append("-k              \t[int]\tk-mer size in the range [4:7] (default=6).\n");
    usage.append("-a              \t[int]\tAmino acid alphabet size (default=21).\n");
    usage.append("--max-seq-len   \t[int]\tMaximum sequence length (default=50000).\n");
    usage.append("--skip          \t[int]\tNumber of skipped k-mers during the index table generation.\n");


    Debug(Debug::ERROR) << usage;
}

void parseArgs(int argc, const char** argv,
               std::string* ffindexSeqDB,
               std::string* ffindexOutDB,
               int *kmerSize, int *alphabetSize,
               int *maxSeqLen, int * skip, int * splitt){
    if (argc < 2){
        printUsage();
        exit(EXIT_FAILURE);
    }
    ffindexSeqDB->assign(argv[1]);
    ffindexOutDB->assign(argv[2]);

    
    int i = 3;
    while (i < argc){
        if (strcmp(argv[i], "-k") == 0){
            if (++i < argc){
                *splitt = atoi(argv[i]);
                i++;
            }
            else {
                printUsage();
                Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
                EXIT(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "-k") == 0){
        if (++i < argc){
            *kmerSize = atoi(argv[i]);
            if (*kmerSize < 4 || *kmerSize > 7){
                Debug(Debug::ERROR) << "Please choose k in the range [4:7].\n";
                EXIT(EXIT_FAILURE);
            }
            i++;
        }
        else {
            printUsage();
            Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    else if (strcmp(argv[i], "-a") == 0){
        if (++i < argc){
            *alphabetSize = atoi(argv[i]);
            i++;
        }
        else {
            printUsage();
            Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
    else if (strcmp(argv[i], "--max-seq-len") == 0){
        if (++i < argc){
            *maxSeqLen = atoi(argv[i]);
            i++;
        }
        else {
            printUsage();
            Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
            EXIT(EXIT_FAILURE);
        }
    }        else if (strcmp(argv[i], "--skip") == 0){
        if (++i < argc){
            *skip = atoi(argv[i]);
            i++;
        }
        else {
            printUsage();
            Debug(Debug::ERROR) << "No value provided for " << argv[i-1] << "\n";
            EXIT(EXIT_FAILURE);
        }
    } else {
            Debug(Debug::ERROR) << "WRONG PARAMETER";
            printUsage();
            exit(EXIT_FAILURE);
            break;
        }
    }
}



int main (int argc, const char * argv[])
{
    
    std::string seqDB = "";
    std::string outDB = "";

    int splitt = 1;
    int kmerSize =  6;
    int alphabetSize = 21;
    int maxSeqLen = 50000;
    int dbFrom = 0;
    int skip = 0;
    std::string scoringMatrixFile = "/Users/mad/Documents/workspace/mmseqs/data/blosum62.out";
    parseArgs(argc, argv, &seqDB, &outDB, &kmerSize, &alphabetSize, &maxSeqLen, &skip, &splitt);
    DBReader dbr(seqDB.c_str(), std::string(seqDB+".index").c_str());
    dbr.open(DBReader::NOSORT);

    
    BaseMatrix* subMat = Prefiltering::getSubstitutionMatrix(scoringMatrixFile, alphabetSize, 8.0f);

    size_t dbSize = dbr.getSize();

    Sequence seq(maxSeqLen, subMat->aa2int, subMat->int2aa, Sequence::AMINO_ACIDS, subMat);
    IndexTable * indexTable = Prefiltering::getIndexTable(&dbr, &seq, alphabetSize, kmerSize, dbFrom, dbFrom + dbSize , skip);
    DBWriter writer(outDB.c_str(), std::string( outDB +".index").c_str(),DBWriter::BINARY_MODE);
    writer.open();
    char * entries = (char *) indexTable->getEntries(); // pointer to first datapoint
    writer.write(entries, indexTable->tableEntriesNum * sizeof(int), "MMSEQSDATA", 0);
    
    char * sizes = (char *) indexTable->getSizes();
    writer.write(sizes, indexTable->tableSize * sizeof(int), "MMSEQSIZES", 0);

    int64_t metadata[] = {kmerSize, alphabetSize, skip, indexTable->tableSize, indexTable->tableEntriesNum};
    char * metadataptr = (char *) &metadata;
    writer.write(metadataptr, 5 * sizeof(int64_t), "MMSEQSMETA", 0);
    
    // write code
    writer.close();
    dbr.close();
    delete indexTable;
    delete subMat;
    return 0;
}
