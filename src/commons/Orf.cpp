#include "Orf.h"
#include <algorithm>
#include <cstring>
#include <cassert>
#include <cctype>

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

Orf::Orf(const char* seq) {
    sequenceLength = strlen(seq);

    sequence = new char[sequenceLength + 1];
    strncpy(sequence, seq, sequenceLength);
    sequence[sequenceLength + 0] = '\0';

    for(size_t i = 0; i < sequenceLength; ++i) {
        sequence[i] = toupper(sequence[i]);
    }

    reverseComplement = new char[sequenceLength + 1];
    strncpy(reverseComplement, sequence, sequenceLength);
    reverseComplement[sequenceLength + 0] = '\0';

    // compute the reverse complement
    std::reverse(reverseComplement, reverseComplement + sequenceLength);
    for(size_t i = 0; i < sequenceLength; ++i) {
        reverseComplement[i] = Complement(reverseComplement[i]);
    }
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
                    int reverseFrames)
{
    std::vector<Orf::SequenceLocation> resForward;
    std::vector<Orf::SequenceLocation> resBackward;
    #pragma omp parallel sections
    {
        // find ORFs on the forward sequence and report them as-is
        #pragma omp section
        {

            FindForwardOrfs(sequence, sequenceLength, resForward,
                            minLength, maxLength, maxGaps, forwardFrames, STRAND_PLUS);
        }
        
        // find ORFs on the reverse complement
        #pragma omp section
        {

            FindForwardOrfs(reverseComplement, sequenceLength, resBackward,
                            minLength, maxLength, maxGaps, reverseFrames, STRAND_MINUS);
        }
    }

    results.insert(results.end(), resForward.begin(), resForward.end());
    results.insert(results.end(), resBackward.begin(), resBackward.end());
}

inline bool isIncomplete(const char* codon)
{
    return codon[0] == 0 || codon[1] == 0 || codon[2] == 0;
}

inline bool IsGapOrN(const char* codon)
{
    if(isIncomplete(codon))
        return false;

    return codon[0] == 'N' || Complement(codon[0]) == '.'
        || codon[1] == 'N' || Complement(codon[1]) == '.'
        || codon[2] == 'N' || Complement(codon[2]) == '.';
}

inline bool isStart(const char* codon) {
    if(isIncomplete(codon))
        return false;

    return (codon[0] == 'A' && codon[1] == 'U' && codon[2] == 'G')
        || (codon[0] == 'A' && codon[1] == 'T' && codon[2] == 'G');
}

inline bool isStop(const char* codon) {
    if(isIncomplete(codon))
        return false;

    return (codon[0] == 'U' && codon[1] == 'A' && codon[2] == 'G')
        || (codon[0] == 'U' && codon[1] == 'A' && codon[2] == 'A')
        || (codon[0] == 'U' && codon[1] == 'G' && codon[2] == 'A')
        || (codon[0] == 'T' && codon[1] == 'A' && codon[2] == 'G')
        || (codon[0] == 'T' && codon[1] == 'A' && codon[2] == 'A')
        || (codon[0] == 'T' && codon[1] == 'G' && codon[2] == 'A');
}

void FindForwardOrfs(const char* sequence, size_t sequenceLength, std::vector<Orf::SequenceLocation>& ranges,
    size_t minLength, size_t maxLength, size_t maxGaps, int frames, Orf::Strand strand) {
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

            if(!(frames & frameLookup[frame])) {
                continue;
            }
            
            if(!isInsideOrf[frame] && isStart(codon)) {
                isInsideOrf[frame] = true;
                from[frame] = position;

                countGaps[frame] = 0;
                countLength[frame] = 0;
            }

            if(isInsideOrf[frame]) {
                countLength[frame]++;

                if(IsGapOrN(codon)) {
                    countGaps[frame]++;
                }
            }

            bool isStopCodon = isStop(codon);
            bool isLastCodonIncomplete = !isStopCodon && (i + FRAMES >= sequenceLength);
            if(isInsideOrf[frame] && (isStopCodon || isLastCodonIncomplete)) {
                isInsideOrf[frame] = false;
                size_t to = position;

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
                bool hasIncompleteStart = false;
                if(from[frame] == frameOffset[frame] && isFirstOrfFound[frame] == false) {
                    isFirstOrfFound[frame] = true;
                    hasIncompleteStart = true;
                }

                ranges.emplace_back(Orf::SequenceLocation{from[frame], to, hasIncompleteStart, isLastCodonIncomplete, strand});
            }
        }
    }
}
