#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>

#include <string>

#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"

extern "C" {
#include "ffutil.h"
}

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
    std::vector<MMseqsParameter> orf_par = {
        Parameters::PARAM_V
    };
    
    Parameters par;
    par.parseParameters(argn, argv, usage, orf_par, 2);
    
    Debug::setDebugLevel(par.verbosity);
    
    const char* in_filename = par.db1.c_str();
    const char* in_index_filename = par.db1Index.c_str();
    
    const char *out_filename  = par.db2.c_str();
    const char *out_index_filename = par.db2Index.c_str();
    
    DBReader reader(in_filename, in_index_filename);
    reader.open(DBReader::NOSORT);
    
    DBWriter writer(out_filename, out_index_filename);
    writer.open();
    
    size_t entries = reader.getSize();
    unsigned int* entries_length = reader.getSeqLens();
    
    for (size_t i = 0; i < entries; ++i) {
        char* key = reader.getDbKey(i);
        char* data = reader.getData(i);
        
        // ignore null char at the end
        unsigned int length = entries_length[i] - 1;
        
        if(length < 3)  {
            Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is too short. Skipping entry.\n";
            continue;
        }

        if((data[length] != '\n' && length % 3 != 0) && (data[length - 1] == '\n' && (length - 1) % 3 != 0)) {
            Debug(Debug::WARNING) << "Nucleotide sequence entry " << key << " length (" << length << ") is not divisible by three. Skipping entry.\n";
            continue;
        }
        
        char aa[length/3 + 1];
        seqan::translateString(aa, data, length,
                               seqan::GeneticCode<seqan::GeneticCodeSpec::CANONICAL>());
        aa[length/3] = '\n';
        
        writer.write(aa, (length / 3) + 1, key);
    }
    
    writer.close();
    reader.close();
    
    return EXIT_SUCCESS;
}

