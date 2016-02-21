#ifndef ORF_H
#define ORF_H

#include <vector>
#include <string>

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

    enum ExtendMode {
        EXTEND_START = (1u << 0),
        EXTEND_END = (1u << 1),
    };

    struct SequenceLocation {
        size_t from, to;
        bool hasIncompleteStart, hasIncompleteEnd;
        Strand strand;
    };

    Orf();

    bool setSequence(const char* sequence, size_t sequenceLength);
    
    ~Orf() {
        cleanup();
    }

    /// Find all ORFs in both orientations that are at least orfMinLength and at most orfMaxLength long.
    /// Report results as SequenceLocations.
    /// seq must be in iupac.
    /// Do not allow more than max_seq_gap consecutive N-or-gap bases in an ORF
    void findAll(std::vector<SequenceLocation> &result,
                 const size_t minLength = 1,
                 const size_t maxLength = SIZE_MAX,
                 const size_t maxGaps = 30,
                 const unsigned int forwardFrames = FRAME_1 | FRAME_2 | FRAME_3,
                 const unsigned int reverseFrames = FRAME_1 | FRAME_2 | FRAME_3,
                 const unsigned int extendMode = 0);

    static void findForward(const char *sequence, const size_t sequenceLength,
                            std::vector<Orf::SequenceLocation> &result,
                            const size_t minLength, const size_t maxLength, const size_t maxGaps,
                            const unsigned int frames, const unsigned int extendMode, const Strand strand);

    std::string view(const SequenceLocation &location);
    
private:
    size_t sequenceLength;
    char* sequence;
    char* reverseComplement;

    void cleanup();
};

#endif
