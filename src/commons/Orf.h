#ifndef ORF_H
#define ORF_H
#pragma once

#include <memory>
#include <vector>
#include <cstdint>
#include <cstdlib>

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
    typedef size_t SequencePosition;

    enum Strand {
        STRAND_PLUS = 1,
        STRAND_MINUS = -1
    };

    // matches NCBI ELim values
    enum Uncertainty {
        UNCERTAINTY_UNKOWN = 0,
        UNCERTAINTY_GREATER = 1,
        UNCERTAINTY_LESS = 2,
        UNCERTAINTY_OTHER = 255
    };

    class SequenceLocation {
    public:
        SequencePosition from, to;
        Strand strand;
        Uncertainty uncertainty_from, uncertainty_to;

        SequenceLocation(SequencePosition from, SequencePosition to, Uncertainty uncertainty_from, Uncertainty uncertainty_to, Strand strand)
            : from(from), to(to), strand(strand), uncertainty_from(uncertainty_from), uncertainty_to(uncertainty_to) {}
    };

    class Range {
    public:
        SequencePosition from, to;
        Range(SequencePosition from, SequencePosition to) : from(from), to(to) {}
    };

    static const size_t k_default_max_seq_gap = 30;

    explicit Orf(const char* sequence);
    
    ~Orf() {
        delete[] revcomp;
    }

    /// Find all ORFs in both orientations that are at least orfMinLength and at most orfMaxLength long.
    /// Report results as SequenceLocations.
    /// seq must be in iupac.
    /// Do not allow more than max_seq_gap consecutive N-or-gap bases in an ORF
    void FindOrfs(std::vector<SequenceLocation>& results,
                  size_t min_length = 1,
                  size_t max_length = SIZE_MAX,
                  size_t max_seq_gap = k_default_max_seq_gap);

    std::unique_ptr<char[]> View(SequenceLocation& orf);
    
private:
    size_t seq_length;
    const char* seq;
    char* revcomp;
};

#endif
