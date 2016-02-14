#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>

#include <string>

#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"

#include <seqan/translation/translation_tables.h>
#include <seqan/translation/translation.h>

#include "Util.h"

#define MAX_FILENAME_LIST_FILES 4096

int translatenucleotide(int argn, const char **argv)
{
    std::string usage;
    usage.append("Translate nucleotide sequences into aminoacid sequences in a FFindex database.\n");
    usage.append("USAGE: <ffindexInDB>  <ffindexOutDB>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.translateNucleotide, 2);

    const char* in_filename = par.db1.c_str();
    const char* in_index_filename = par.db1Index.c_str();
    
    const char *out_filename  = par.db2.c_str();
    const char *out_index_filename = par.db2Index.c_str();
    
    DBReader<std::string> reader(in_filename, in_index_filename);
    reader.open(DBReader<std::string>::NOSORT);
    
    DBWriter writer(out_filename, out_index_filename);
    writer.open();
    
    size_t entries = reader.getSize();

    for (size_t i = 0; i < entries; ++i) {
        std::string key = reader.getDbKey(i);
        char* data = reader.getData(i);
        
        // ignore null char at the end
        unsigned int length = reader.getSeqLens(i) - 1;
        
        if(length < 3)  {
            Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is too short. Skipping entry.\n";
            continue;
        }

        if((data[length] != '\n' && length % 3 != 0) && (data[length - 1] == '\n' && (length - 1) % 3 != 0)) {
            Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is not divisible by three. Skipping entry.\n";
            continue;
        }
        
        char aa[length/3 + 1];

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
        
        writer.write(aa, (length / 3) + 1, (char*)key.c_str());
    }
    
    writer.close();
    reader.close();
    
    return EXIT_SUCCESS;
}

