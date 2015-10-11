#ifndef ORF_H
#define ORF_H
#pragma once

#include <memory>
#include <vector>

//
// class Orf implements a simple open reading frame search.
//
// This class reports, as a series of seq-locs, all the possible open
// reading frames in a given sequence.  The sequence can be represented
// as either
/// This class provides functions for finding all the ORFs
/// of a specified minimum length in a DNA sequence.
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

    explicit Orf(const char* sequence);
    
    ~Orf() {
        delete[] reverseComplement;
    }

    /// Find all ORFs in both orientations that are at least orfMinLength and at most orfMaxLength long.
    /// Report results as SequenceLocations.
    /// seq must be in iupac.
    /// Do not allow more than max_seq_gap consecutive N-or-gap bases in an ORF
    void FindOrfs(std::vector<SequenceLocation>& results,
                  size_t minLength = 1,
                  size_t maxLength = SIZE_MAX,
                  size_t maxGaps = 30);

    char* View(SequenceLocation& location);
    
private:
    size_t sequenceLength;
    char* sequence;
    char* reverseComplement;
};

void FindForwardOrfs(const char* sequence, size_t sequenceLength, std::vector<Orf::SequenceLocation>& ranges,
                     size_t minLength, size_t maxLength, size_t maxGaps, int frames, Orf::Strand strand);

#endif
