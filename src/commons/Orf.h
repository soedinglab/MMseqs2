#ifndef ORF_H
#define ORF_H

#include <vector>
#include <string>
#include <cstdlib>

class Orf
{
public:
    enum Strand {
        STRAND_PLUS = 1,
        STRAND_MINUS = -1
    };

    enum Frame {
        FRAME_1 = (1u << 0),
        FRAME_2 = (1u << 1),
        FRAME_3 = (1u << 2)
    };

    struct SequenceLocation {
        size_t from, to;
        bool hasIncompleteStart, hasIncompleteEnd;
        Strand strand;
    };

    Orf();

    bool setSequence(const char* sequence);
    
    ~Orf() {
        cleanup();
    }

    /// Find all ORFs in both orientations that are at least orfMinLength and at most orfMaxLength long.
    /// Report results as SequenceLocations.
    /// seq must be in iupac.
    /// Do not allow more than max_seq_gap consecutive N-or-gap bases in an ORF
    void FindOrfs(std::vector<SequenceLocation>& results,
                  size_t minLength = 1,
                  size_t maxLength = SIZE_MAX,
                  size_t maxGaps = 30,
                  int forwardFrames = FRAME_1 | FRAME_2 | FRAME_3,
                  int reverseFrames = FRAME_1 | FRAME_2 | FRAME_3);

    std::string View(SequenceLocation& location);
    
private:
    size_t sequenceLength;
    char* sequence;
    char* reverseComplement;

    void cleanup() {
        if (sequence) {
            free(sequence);
            sequence = NULL;
        }
        if (reverseComplement) {
            free(reverseComplement);
            reverseComplement = NULL;
        }
    }
};

void FindForwardOrfs(const char* sequence, size_t sequenceLength, std::vector<Orf::SequenceLocation>& ranges,
                     size_t minLength, size_t maxLength, size_t maxGaps, int frames, Orf::Strand strand);

#endif
