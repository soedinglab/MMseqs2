#include "Orf.h"
#include <algorithm>
#include <cstring>

//note: N->N, S->S, W->W, U->A, T->A
static const char* iupac_revcomp_table =
"................................................................"
".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
"................................................................"
"................................................................";

inline char Complement(const char c)
{
    return iupac_revcomp_table[static_cast<unsigned char>(c)];
}

Orf::Orf(const char* sequence) {
    seq_length  = strlen(sequence);
    seq = sequence;
    revcomp = new char[seq_length + 1];
    strncpy(revcomp, seq, seq_length);
    revcomp[seq_length] = '\0';
    
    // compute the reverse complement
    std::reverse(revcomp, revcomp + seq_length);
    for(size_t i = 0; i < seq_length; ++i) {
        revcomp[i] = Complement(revcomp[i]);
    }
}

bool IsGapOrN(const char* codon)
{
    return codon[0] == 'N' || codon[0] == 'n' || Complement(codon[0]) == '.';
}

bool isStart(const char* codon) {
    return  (codon[0] == 'A' && codon[1] == 'U' && codon[2] == 'G')
          ||(codon[0] == 'A' && codon[1] == 'T' && codon[2] == 'G');
}

bool isStop(const char* codon) {
    return  (codon[0] == 'U' && codon[1] == 'A' && codon[2] == 'G')
          ||(codon[0] == 'U' && codon[1] == 'A' && codon[2] == 'A')
          ||(codon[0] == 'U' && codon[1] == 'G' && codon[2] == 'A')
          ||(codon[0] == 'T' && codon[1] == 'A' && codon[2] == 'G')
          ||(codon[0] == 'T' && codon[1] == 'A' && codon[2] == 'A')
          ||(codon[0] == 'T' && codon[1] == 'G' && codon[2] == 'A');
}

static void FindForwardOrfs(const char* sequence, size_t seq_length, std::vector<Orf::Range>& ranges, size_t max_seq_gap) {

    // An open reading frame can beginning in any of the three codon start position
    // Frame 0:  AGA ATT GCC TGA ATA AAA GGA TTA CCT TGA TAG GGT AAA
    // Frame 1: A GAA TTG CCT GAA TAA AAG GAT TAC CTT GAT AGG GTA AA
    // Frame 2: AG AAT TGC CTG AAT AAA AGG ATT ACC TTG ATA GGG TAA A
    const int FRAMES = 3;

    // We want to walk over the memory only once so we calculate which codon we are in
    // and save the values of our state machine in these arrays

    // we also initialize our state machine with being inside an orf
    // this is to handle edge case 1 where we find an end codon but no start codon
    // in this case we just add an orf from the start to the found end codon
    bool inOrf[FRAMES] = {true, true, true};
    bool firstOrfFound[FRAMES] = {false, false, false};

    size_t currentGaps[FRAMES] = {0, 0, 0};
    size_t currentLength[FRAMES] = {0, 0, 0};

    // Offset the values by reading frame
    Orf::SequencePosition from[FRAMES] = {0, 1, 2};
    Orf::SequencePosition to[FRAMES] = {0, 0, 0};

    for (size_t i = 0;  i < seq_length - 2;  i += 3) {
        for(size_t pos = i; pos < i + FRAMES; pos++) {
            const char* current_codon = &sequence[pos];
            size_t cur_frame = pos % FRAMES;

            if(!inOrf[cur_frame] && isStart(current_codon)) {
                inOrf[cur_frame] = true;
                from[cur_frame] = pos;
            }

            if(inOrf[cur_frame]) {
                currentLength[cur_frame]++;
            }

            // ignore orfs with too many gaps or unknown codons
            if(inOrf[cur_frame] && IsGapOrN(current_codon)) {
                currentGaps[cur_frame]++;

                if(currentGaps[cur_frame] >= max_seq_gap) {
                    inOrf[cur_frame] = false;
                    currentGaps[cur_frame] = 0;
                    currentLength[cur_frame] = 0;
                }
            }

            if(inOrf[cur_frame] && isStop(current_codon)) {
                inOrf[cur_frame] = false;
                to[cur_frame] = pos;

                // edge case 1: see above
                if(from[cur_frame] == 0 && firstOrfFound[cur_frame] == false) {
                    firstOrfFound[cur_frame] = true;
                    ranges.emplace_back(from[cur_frame], to[cur_frame]);
                    currentGaps[cur_frame] = 0;
                    currentLength[cur_frame] = 0;
                    continue;
                }

                ranges.emplace_back(from[cur_frame], to[cur_frame]);
                currentGaps[cur_frame] = 0;
                currentLength[cur_frame] = 0;
                continue;
            }
        }
    }

    // edge case 2: we did not find an end codon, add the last orf
    // from the last found start to the end of the sequence
    for(size_t frame = 0; frame < FRAMES; frame++) {
        if(from[frame] > to[frame]) {
            Orf::SequencePosition to = seq_length - ((seq_length - from[frame]) % 3);
            ranges.emplace_back(from[frame], to);
        }
        
        currentGaps[frame] = 0;
        currentLength[frame] = 0;
    }
}

std::unique_ptr<char> Orf::View(SequenceLocation& loc) {
    size_t length = loc.to - loc.from;
    std::unique_ptr<char> buffer(new char[length + 2]);
    
    if(loc.strand == Orf::STRAND_PLUS) {
        strncpy(buffer.get(), seq + loc.from, length);
    } else {
        strncpy(buffer.get(), revcomp + loc.from, length);
    }
    
    buffer.get()[length] = '\n';
    buffer.get()[length + 1] = '\0';
    
    return buffer;
}

//
// find ORFs in a string
void Orf::FindOrfs(std::vector<Orf::SequenceLocation>& results,
                    size_t min_length,
                    size_t max_length,
                    size_t max_seq_gap)
{
    std::vector<Orf::Range> ranges;

    // find ORFs on the forward sequence and report them as-is
    FindForwardOrfs(seq, seq_length, ranges, max_seq_gap);
    for(std::vector<Orf::Range>::const_iterator it = ranges.begin(); it != ranges.end(); ++it) {
        Orf::SequencePosition from = it->from, to = it->to;
        Orf::Uncertainty uncertainty_from = Orf::UNCERTAINTY_UNKOWN;
        Orf::Uncertainty uncertainty_to   = Orf::UNCERTAINTY_UNKOWN;

        if (from < 3) {
            // "beginning" of ORF at beginning of sequence
           uncertainty_from = Orf::UNCERTAINTY_LESS;
        }

        if (to >= seq_length) {
            // "end" of ORF is really end of sequence
            uncertainty_to = Orf::UNCERTAINTY_GREATER;
        }

        size_t length = to - from;
        if(length < min_length || length >= max_length)
            continue;

        results.emplace_back(from, to, uncertainty_from, uncertainty_to, Orf::STRAND_PLUS);
    }
    ranges.clear();

    // find ORFs on the reverse complement
    FindForwardOrfs(revcomp, seq_length, ranges, max_seq_gap);
    for(std::vector<Orf::Range>::const_iterator it = ranges.begin(); it != ranges.end(); ++it) {
        Orf::SequencePosition from =  it->from ;
        Orf::Uncertainty uncertainty_from = Orf::UNCERTAINTY_UNKOWN;
        Orf::Uncertainty uncertainty_to   = Orf::UNCERTAINTY_UNKOWN;
        Orf::SequencePosition to =  it->to;

        if (from < 3) {
            // "end" of ORF is beginning of sequence
           uncertainty_from = Orf::UNCERTAINTY_LESS;
        }

        if (to >= seq_length) {
            // "beginning" of ORF is really end of sequence
            uncertainty_to = Orf::UNCERTAINTY_GREATER;
        }

        size_t length = to - from;
        if(length < min_length || length >= max_length)
            continue;

        results.emplace_back(from, to, uncertainty_from, uncertainty_to, Orf::STRAND_MINUS);
    }
}
