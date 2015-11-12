#define _GNU_SOURCE 1
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include <unistd.h>
#include <Util.h>
#include <climits>

#include "Parameters.h"
#include "DBWriter.h"

#include "Orf.h"

#include "kseq.h"

KSEQ_INIT(int, read)

int orfFastaToFFindex(
        const char* fasta_filename,
        const char* data_filename, const char* index_filename,
        const char* data_filename_hdr, const char* index_filename_hdr,
        Parameters* par)
{
    FILE* fasta_file = fopen(fasta_filename, "r");
    if (fasta_file == NULL) {
        perror(fasta_filename);
        return EXIT_FAILURE;
    }

    DBWriter seq_writer(data_filename, index_filename);
    seq_writer.open();
    DBWriter hdr_writer(data_filename_hdr, index_filename_hdr);
    hdr_writer.open();

    kseq_t* seq = kseq_init(fileno(fasta_file));

    char header_buffer[LINE_MAX];

    size_t total_num = 0;
    while (kseq_read(seq) >= 0) {
        Orf orf(seq->seq.s);
        std::vector<Orf::SequenceLocation> res;
        orf.FindOrfs(res, par->orfMinLength, par->orfMaxLength, par->orfMaxGaps);

        size_t orf_num = 0;
        for (std::vector<Orf::SequenceLocation>::const_iterator it = res.begin(); it != res.end(); ++it) {
            Orf::SequenceLocation loc = *it;

            std::string id;
            if (par->orfUseNumericIndices) {
                id = SSTR(UINT_MAX - total_num);
            } else {
                id = Util::parseFastaHeader(seq->name.s) + std::string("_") + SSTR(orf_num);
            }

            if (id.length() >= 31) {
                std::cerr << "Id: " << id << " is too long. Maximum of 32 characters are allowed." << std::endl;
                EXIT(EXIT_FAILURE);
            }

            if (par->orfSkipIncomplete && (loc.hasIncompleteStart == true || loc.hasIncompleteEnd == true))
                continue;

            if (seq->comment.l) {
                snprintf(header_buffer, LINE_MAX, "%s %s [Orf: %zu, %zu, %d, %d, %d]\n", seq->name.s, seq->comment.s, loc.from, loc.to, loc.strand, loc.hasIncompleteStart, loc.hasIncompleteEnd);
            }
            else {
                snprintf(header_buffer, LINE_MAX, "%s [Orf: %zu, %zu, %d, %d, %d]\n", seq->name.s, loc.from, loc.to, loc.strand, loc.hasIncompleteStart, loc.hasIncompleteEnd);
            }

            hdr_writer.write(header_buffer, strlen(header_buffer), (char *)id.c_str());

            char* sequence = orf.View(loc);
            size_t length = strlen(sequence);
            seq_writer.write(sequence, length, (char *)id.c_str());
            delete sequence;

            orf_num++;
            total_num++;
        }
    }

    kseq_destroy(seq);
    hdr_writer.close();
    seq_writer.close();

    fclose(fasta_file);

    return EXIT_SUCCESS;
}

int extractorf(int argn, const char** argv)
{
    std::string usage;
    usage.append("Extract all open reading frames from a nucleotide fasta file into a ffindex database.\n");
    usage.append("USAGE: <fastaDB> <ffindexDB>\n");
    usage.append("\nDesigned and implemented by Milot Mirdita <milot@mirdita.de>.\n");

    Parameters par;
    par.parseParameters(argn, argv, usage, par.extractorf, 2);

    const char* input = par.db1.c_str();

    const char* data_filename = par.db2.c_str();
    const char* index_filename = par.db2Index.c_str();

    std::string data_filename_hdr_str(data_filename);
    data_filename_hdr_str.append("_h");
    const char* data_filename_hdr = data_filename_hdr_str.c_str();

    std::string index_filename_hdr_str(data_filename);
    index_filename_hdr_str.append("_h.index");
    const char* index_filename_hdr = index_filename_hdr_str.c_str();

    return orfFastaToFFindex(input,
                             data_filename, index_filename,
                             data_filename_hdr, index_filename_hdr,
                             &par);
}
