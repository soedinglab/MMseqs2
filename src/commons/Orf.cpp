#include "Orf.h"
#include <cstring>
#include <cassert>

//note: N->N, S->S, W->W, U->A, T->A
static const char* iupacReverseComplementTable =
"................................................................"
".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
"................................................................"
"................................................................";

inline char Complement(const char c)
{
    return iupacReverseComplementTable[static_cast<unsigned char>(c)];
}

Orf::Orf() : sequence(NULL), reverseComplement(NULL) {}

bool Orf::setSequence(const char* seq) {
    cleanup();

    sequenceLength = strlen(seq);

    if(sequenceLength < 3)
        return false;

    sequence = strdup(seq);
    for(size_t i = 0; i < sequenceLength; ++i) {
        sequence[i] = static_cast<char>(toupper(static_cast<int>(seq[i])));
    }

    reverseComplement = strdup(sequence);
    for(size_t i = 0; i < sequenceLength; ++i) {
        reverseComplement[i] = Complement(sequence[sequenceLength - i - 1]);
        if (reverseComplement[i] == '.') {
            return false;
        }
    }

    return true;
}

std::string Orf::View(SequenceLocation& location) {
    assert(location.to > location.from);
    
    size_t length = location.to - location.from;
    if(location.strand == Orf::STRAND_PLUS) {
        return std::string(&sequence[location.from], length);
    } else {
        return std::string(&reverseComplement[location.from], length);
    }
}

//
// find ORFs in a string
void Orf::FindOrfs(std::vector<Orf::SequenceLocation>& results,
                    size_t minLength,
                    size_t maxLength,
                    size_t maxGaps,
                    int forwardFrames,
                    int reverseFrames,
                    int extendMode)
{
    // find ORFs on the forward sequence and report them as-is
    FindForwardOrfs(sequence, sequenceLength, results,
                        minLength, maxLength, maxGaps, forwardFrames, extendMode, STRAND_PLUS);

    // find ORFs on the reverse complement
    FindForwardOrfs(reverseComplement, sequenceLength, results,
                        minLength, maxLength, maxGaps, reverseFrames, extendMode, STRAND_MINUS);
}

inline bool isIncomplete(const char* codon) {
    return codon[0] == 0 || codon[1] == 0 || codon[2] == 0;
}

inline bool isGapOrN(const char *codon) {
    return codon[0] == 'N' || Complement(codon[0]) == '.'
        || codon[1] == 'N' || Complement(codon[1]) == '.'
        || codon[2] == 'N' || Complement(codon[2]) == '.';
}

inline bool isStart(const char* codon) {
    return (codon[0] == 'A' && codon[1] == 'U' && codon[2] == 'G')
        || (codon[0] == 'A' && codon[1] == 'T' && codon[2] == 'G');
}

inline bool isStop(const char* codon) {
    return (codon[0] == 'U' && codon[1] == 'A' && codon[2] == 'G')
        || (codon[0] == 'U' && codon[1] == 'A' && codon[2] == 'A')
        || (codon[0] == 'U' && codon[1] == 'G' && codon[2] == 'A')
        || (codon[0] == 'T' && codon[1] == 'A' && codon[2] == 'G')
        || (codon[0] == 'T' && codon[1] == 'A' && codon[2] == 'A')
        || (codon[0] == 'T' && codon[1] == 'G' && codon[2] == 'A');
}

void FindForwardOrfs(const char* sequence, size_t sequenceLength, std::vector<Orf::SequenceLocation>& ranges,
    size_t minLength, size_t maxLength, size_t maxGaps, int frames, int extendMode, Orf::Strand strand) {
    if (frames == 0)
        return;

    // An open reading frame can beginning in any of the three codon start position
    // Frame 0:  AGA ATT GCC TGA ATA AAA GGA TTA CCT TGA TAG GGT AAA
    // Frame 1: A GAA TTG CCT GAA TAA AAG GAT TAC CTT GAT AGG GTA AA
    // Frame 2: AG AAT TGC CTG AAT AAA AGG ATT ACC TTG ATA GGG TAA A
    const int FRAMES = 3;
    const int frameLookup[FRAMES] = {Orf::FRAME_1, Orf::FRAME_2, Orf::FRAME_3};
    const size_t frameOffset[FRAMES] = {0, 1 , 2};

    // We want to walk over the memory only once so we calculate which codon we are in
    // and save the values of our state machine in these arrays

    // we also initialize our state machine with being inside an orf
    // this is to handle edge case 1 where we find an end codon but no start codon
    // in this case we just add an orf from the start to the found end codon
    bool isInsideOrf[FRAMES]     = {true,  true,  true };
    bool isFirstOrfFound[FRAMES] = {false, false, false};

    size_t countGaps[FRAMES]   = {0, 0, 0};
    size_t countLength[FRAMES] = {0, 0, 0};

    // Offset the start position by reading frame
    size_t from[FRAMES] = {frameOffset[0], frameOffset[1], frameOffset[2]};

    for (size_t i = 0;  i < sequenceLength - (FRAMES - 1);  i += FRAMES) {
        for(size_t position = i; position < i + FRAMES; position++) {
            const char* codon = sequence + position;
            size_t frame = position % FRAMES;

            // skip frames outside of out the frame mask
            if(!(frames & frameLookup[frame])) {
                continue;
            }

            bool nextIncomplete = isIncomplete(codon + FRAMES);

            // if we have the start extend mode the returned orf should return the longest
            // possible orf with possibly multiple start codons
            bool shouldStart = false;
            if((extendMode & Orf::EXTEND_START)) {
                shouldStart = isInsideOrf[frame] == false && isStart(codon);
            } else {
                shouldStart = isStart(codon);
            }

            // do not start a new orf on the last codon
            if(shouldStart && nextIncomplete == false) {
                isInsideOrf[frame] = true;
                from[frame] = position;

                countGaps[frame] = 0;
                countLength[frame] = 0;
            }

            if(isInsideOrf[frame]) {
                countLength[frame]++;

                if(isGapOrN(codon)) {
                    countGaps[frame]++;
                }
            }

            bool thisIncomplete = isIncomplete(codon);
            bool isOrfEnd = isInsideOrf[frame] && isStop(codon);
            if(isOrfEnd || thisIncomplete || nextIncomplete) {
                // possibly bail early if we have an orf shorter than minLength
                // so we can find another longer one
                 if ((extendMode & Orf::EXTEND_END) && countLength[frame] <= minLength) {
                    continue;
                 }

                isInsideOrf[frame] = false;

                // we do not include the stop codon here
                size_t to = position + (nextIncomplete ? 3 : 0);

                // this could happen if the first codon is a stop codon
                if(to == from[frame])
                    continue;

                assert(to > from[frame]);

                // ignore orfs with too many gaps or unknown codons
                // also ignore orfs shorter than the min size and longer than max
                if ((countGaps[frame] > maxGaps)
                || (countLength[frame] > maxLength)
                || (countLength[frame] <= minLength)) {
                    continue;
                }

                // edge case 1: see above
                if(from[frame] == frameOffset[frame] && isFirstOrfFound[frame] == false) {
                    isFirstOrfFound[frame] = true;
                }

                ranges.emplace_back(Orf::SequenceLocation{from[frame], to,
                    isFirstOrfFound[frame], thisIncomplete || nextIncomplete, strand});
            }
        }
    }
}
