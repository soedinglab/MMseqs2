/*
 *
 * Adapted from NCBI C++ Toolkit:
 * orf.cpp 65735 2014-12-23 18:23:27Z astashya $
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
 * Authors:  Mike DiCuccio
 */


#include "Orf.h"

#include <algorithm>
#include <vector>

#include <cstdio>
#include <cstring>

char Complement(const char c)
{
    //note: N->N, S->S, W->W, U->A, T->A
    static const char* iupac_revcomp_table =
        "................................................................"
        ".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
        "................................................................"
        "................................................................";
    return iupac_revcomp_table[static_cast<unsigned char>(c)];
}

bool IsGapOrN(const char c)
{
    return c == 'N' || c == 'n' || Complement(c) == '.';
}
/*

    for (size_t i = 2;  i < seq.size() - 2;  i += 3) {
        for (size_t pos = 0;  pos < 3;  pos++) {
            if () {
                stops[(i + pos - 2) % 3].push_back(i + pos - 2);
            }
        }
    }
*/

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

static void FindForwardOrfs(const char* sequence, size_t seq_length, std::vector<Orf::Range>& ranges, unsigned int min_length_bp, size_t max_seq_gap) {

    const int FRAMES = 3;
    // blbala ill do this tomorrow
    bool inOrf[FRAMES] = {true, true, true};

    bool firstOrfFound[FRAMES] = {false, false, false};
    Orf::SequencePosition from[FRAMES] = {0, 1, 2};
    Orf::SequencePosition to[FRAMES] = {0, 0, 0};
    for (size_t i = 0;  i < seq_length - 2;  i += 3) {
        for(size_t pos = i; pos < i + FRAMES; pos++) {
            const char* current_codon = &sequence[pos];
            size_t cur_frame = pos % FRAMES;

            if(!inOrf[cur_frame] && isStart(current_codon)) {
                inOrf[cur_frame] = true;
                from[cur_frame] = pos;
                continue;
            }

            if(inOrf[cur_frame] && isStop(current_codon)) {
                inOrf[cur_frame] = false;
                to[cur_frame] = pos;

                if(from[cur_frame] == 0 && firstOrfFound[cur_frame] == false) {
                    firstOrfFound[cur_frame] = true;
                    ranges.push_back(Orf::Range(from[cur_frame], to[cur_frame]));
                    continue;
                }

                ranges.push_back(Orf::Range(from[cur_frame], to[cur_frame]));

                continue;
            }
        }
    }

    for(size_t frame = 0; frame < FRAMES; frame++) {
        if(from[frame] > to[frame]) {
            Orf::SequencePosition to = seq_length - ((seq_length - from[frame]) % 3);
            ranges.push_back(Orf::Range(from[frame], to));
        }
    }
}

//
// find ORFs in a string
void Orf::FindOrfs(char* seq,
                    std::vector<Orf::SequenceLocation*>& results,
                    unsigned int min_length_bp,
                    size_t max_seq_gap)
{
    size_t seq_length = strlen(seq);
    char *seq_end = strchr(seq, 0);

    std::vector<Orf::Range> ranges;

    // This code might be sped up by a factor of two
    // by use of a state machine that does all six frames
    // in a single pass.

    // find ORFs on the forward sequence and report them as-is
    FindForwardOrfs(seq, seq_length, ranges, min_length_bp, max_seq_gap);
    for(std::vector<Orf::Range>::const_iterator it = ranges.begin(); it != ranges.end(); ++it) {
        Orf::SequencePosition from = it->from, to = it->to;
        Orf::Uncertainty uncertainty_from = Orf::UNCERTAINTY_UNKOWN;
        Orf::Uncertainty uncertainty_to   = Orf::UNCERTAINTY_UNKOWN;
/*
        if (from < 3) {
            // "beginning" of ORF at beginning of sequence
           uncertainty_from = Orf::UNCERTAINTY_LESS;
        }
        if (to + 3 >= seq_length) {
            // "end" of ORF is really end of sequence
            uncertainty_to = Orf::UNCERTAINTY_GREATER;
        } else {
            // ORF was ended by a stop, rather than end of sequence
            to += 3;
        }
*/
        results.push_back(new Orf::SequenceLocation(from, to, uncertainty_from, uncertainty_to, Orf::STRAND_PLUS));
    }

    // find ORFs on the complement and munge the numbers
    ranges.clear();

    // compute the reverse complement;
    std::reverse(seq, seq_end);
    size_t i = 0;
    while(seq[i] != NULL) {
        seq[i] = Complement(seq[i]);

        i++;
    }

    FindForwardOrfs(seq, seq_length, ranges, min_length_bp, max_seq_gap);
    for(std::vector<Orf::Range>::const_iterator it = ranges.begin(); it != ranges.end(); ++it) {
        Orf::SequencePosition from =  it->from ;
        Orf::Uncertainty uncertainty_from = Orf::UNCERTAINTY_UNKOWN;
        Orf::Uncertainty uncertainty_to   = Orf::UNCERTAINTY_UNKOWN;
                Orf::SequencePosition to =  it->to;

/*
        if (from < 3) {
            // "end" of ORF is beginning of sequence
           uncertainty_from = Orf::UNCERTAINTY_LESS;
        } else {
            // ORF was ended by a stop, rather than beginning of sequence
            from -= 3;
        }
        if (to + 3 >= seq_length) {
            // "beginning" of ORF is really end of sequence
            uncertainty_to = Orf::UNCERTAINTY_GREATER;
        }
*/
        results.push_back(new Orf::SequenceLocation(from, to, uncertainty_from, uncertainty_to, Orf::STRAND_MINUS));
    }
}
