#include <unistd.h>
#include <string>
#include <limits.h>
#include <stdlib.h>

#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"

#include <seqan/translation/translation_tables.h>
#include <seqan/translation/translation.h>

#include "Util.h"

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

    DBReader<std::string> reader(par.db1.c_str(), par.db1Index.c_str());
    reader.open(DBReader<std::string>::NOSORT);
    
    DBWriter writer(par.db2.c_str(), par.db2Index.c_str());
    writer.open();
    
    size_t entries = reader.getSize();

    for (size_t i = 0; i < entries; ++i) {
        std::string key = reader.getDbKey(i);
        char* data = reader.getData(i);
        
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


        
        char* aa = new char[length/3 + 1];

        switch(par.translationTable) {
            case seqan::GeneticCodeSpec::CANONICAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::CANONICAL>());
                break;
            case seqan::GeneticCodeSpec::VERT_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::VERT_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::YEAST_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::YEAST_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::MOLD_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::MOLD_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::INVERT_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::INVERT_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::CILIATE:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::CILIATE>());
                break;
            case seqan::GeneticCodeSpec::FLATWORM_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::FLATWORM_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::EUPLOTID:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::EUPLOTID>());
                break;
            case seqan::GeneticCodeSpec::PROKARYOTE:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::PROKARYOTE>());
                break;
            case seqan::GeneticCodeSpec::ALT_YEAST:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::ALT_YEAST>());
                break;
            case seqan::GeneticCodeSpec::ASCIDIAN_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::ASCIDIAN_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::ALT_FLATWORM_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::ALT_FLATWORM_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::BLEPHARISMA:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::BLEPHARISMA>());
                break;
            case seqan::GeneticCodeSpec::CHLOROPHYCEAN_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::CHLOROPHYCEAN_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::TREMATODE_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::TREMATODE_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::SCENEDESMUS_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::SCENEDESMUS_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::THRAUSTOCHYTRIUM_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::THRAUSTOCHYTRIUM_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::PTEROBRANCHIA_MITOCHONDRIAL:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::PTEROBRANCHIA_MITOCHONDRIAL>());
                break;
            case seqan::GeneticCodeSpec::GRACILIBACTERIA:
                seqan::translateString(aa, data, length,
                                       seqan::GeneticCode<seqan::GeneticCodeSpec::GRACILIBACTERIA>());
                break;
            default:
                Debug(Debug::ERROR) << "Invalid translation table selected!\n";
                return EXIT_FAILURE;
        }

        aa[length/3] = '\n';

        writer.writeData(aa, (length / 3) + 1, (char *) key.c_str());
        delete[] aa;
    }

    writer.close();
    reader.close();

    return EXIT_SUCCESS;
}

