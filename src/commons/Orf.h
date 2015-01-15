/*
 *
 * Adapted from NCBI C++ Toolkit:
 *   $Id: orf.hpp 65735 2014-12-23 18:23:27Z astashya $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Josh Cherry
 *
 * File Description:  code for finding and representing ORFs
 *
 */

#include <vector>


//
// class COrf implements a simple open reading frame search.
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
        //TODO: NCBI CRange defines to as to(to+1)??
        Range(SequencePosition from, SequencePosition to) : from(from), to(to) {}
    };

    static const size_t k_default_max_seq_gap = 30;

    /// Find all ORFs in both orientations that
    /// are at least min_length_bp long.
    /// Report results as Seq-locs.
    /// seq must be in iupac.
    /// If allowable_starts is empty (the default), any sense codon can begin
    /// an ORF.  Otherwise, only codons in allowable_starts can do so.
    /// Do not allow more than max_seq_gap consecutive N-or-gap bases in an ORF
    static void FindOrfs(char* seq,
                         std::vector<SequenceLocation*>& results,
                         unsigned int min_length_bp = 3,
                         size_t max_seq_gap = k_default_max_seq_gap);
};
