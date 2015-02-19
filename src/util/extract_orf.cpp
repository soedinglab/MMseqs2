#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>

#include <string>

#include "Parameters.h"
#include "DBWriter.h"

#include "Orf.h"
#include "Util.h"

#include "kseq.h"
#define MAX_FILENAME_LIST_FILES 4096

KSEQ_INIT(int, read)

int extractorf(int argn, const char **argv)
{
    std::string usage;
    usage.append("Extract all open reading frames from a nucleotide fasta file into a ffindex database.\n");
    usage.append("USAGE: <fastaDB>  <ffindexDB>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");
    std::vector<MMseqsParameter> orf_par = {
        Parameters::PARAM_ORF_MIN_LENGTH,
        Parameters::PARAM_ORF_MAX_LENGTH,
        Parameters::PARAM_ORF_MAX_GAP
    };
    
    Parameters par;
    par.parseParameters(argn, (char**)argv, usage, orf_par, 2);
    
    const char* fasta_folder = par.db1.c_str();
    
    const char *fasta_filename = fasta_folder;
    FILE* fasta_file = fopen(fasta_filename, "r");
    if(fasta_file == NULL) { perror(fasta_filename); return EXIT_FAILURE; }
    
    const char *data_filename  = par.db2.c_str();
    const char *index_filename = par.db2Index.c_str();

    std::string data_filename_hdr_str(data_filename);
    data_filename_hdr_str.append("_h");
    const char *data_filename_hdr  = data_filename_hdr_str.c_str();

    std::string index_filename_hdr_str(data_filename);
    index_filename_hdr_str.append("_h.index");
    const char *index_filename_hdr = index_filename_hdr_str.c_str();
    
    DBWriter seq_writer(data_filename, index_filename);
    seq_writer.open();
    DBWriter hdr_writer(data_filename_hdr, index_filename_hdr);
    hdr_writer.open();

    size_t entries_num = 0;

    kseq_t *seq = kseq_init(fileno(fasta_file));

    std::vector<Orf::SequenceLocation> res;
    
    char header_buffer[LINE_MAX];
    char id_buffer[LINE_MAX];
    while (kseq_read(seq) >= 0) {
        std::string id = Util::parseFastaHeader(std::string(seq->name.s));
        
        Orf orf(seq->seq.s);
        orf.FindOrfs(res, par.min_length, par.max_length, par.max_gaps);

        size_t orfs_num = 0;
        for(std::vector<Orf::SequenceLocation>::const_iterator it = res.begin(); it != res.end(); it++) {
            snprintf(id_buffer, LINE_MAX, "%s_%zu", id.c_str(), orfs_num);

            Orf::SequenceLocation loc = (Orf::SequenceLocation)*it;

            if(seq->comment.l)  {
                snprintf(header_buffer, LINE_MAX, "%s %s [Orf: %zu, %zu, %d, %d, %d]\n", seq->name.s, seq->comment.s, loc.from, loc.to, loc.strand, loc.uncertainty_from, loc.uncertainty_to);
            } else {
                snprintf(header_buffer, LINE_MAX, "%s [Orf: %zu, %zu, %d, %d, %d]\n", seq->name.s, loc.from, loc.to, loc.strand, loc.uncertainty_from, loc.uncertainty_to);
            }
            hdr_writer.write(header_buffer, strlen(header_buffer), id_buffer);

            std::unique_ptr<char> sequence(orf.View(loc));
            seq_writer.write(sequence.get(), strlen(sequence.get()), id_buffer);

            orfs_num++;
            entries_num++;
        }
        
        res.clear();
    }

    kseq_destroy(seq);
    hdr_writer.close();
    seq_writer.close();
    fclose(fasta_file);

    return EXIT_SUCCESS;
}

