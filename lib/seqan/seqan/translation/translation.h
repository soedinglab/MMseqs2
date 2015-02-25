// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Code for Dna(5) to AminoAcid Translation
// ==========================================================================


#ifndef INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_
#define INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_

#include <cassert>
#include <cstddef>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function _ord()
// --------------------------------------------------------------------------


inline size_t
_ord(unsigned char const & c)
{
    return TranslateTableCharToDna5_<>::VALUE[c];
}

// --------------------------------------------------------------------------
// Function _translateTriplet()
// --------------------------------------------------------------------------


template <GeneticCodeSpec CODE_SPEC>
inline char _translateTriplet(size_t const & c1,
                  size_t const & c2,
                  size_t const & c3,
                  GeneticCode<CODE_SPEC> const & /**/)
{
    return (( c1 > 3 ) || ( c2 > 3 ) || ( c3 > 3 ) )
            ? 'X'
            : TranslateTableDnaToAminoAcid_<
                GeneticCode<CODE_SPEC> >::VALUE[c1][c2][c3];
}

// --------------------------------------------------------------------------
// Function _translateString()
// --------------------------------------------------------------------------

template <GeneticCodeSpec CODE_SPEC>
    inline void translateString(char * target,
                 char * source,
                 size_t source_length,
                 GeneticCode<CODE_SPEC> const & /**/)
{
    for (size_t i = 0; i+2 < source_length; i+=3)
    {
        target[i/3] = _translateTriplet(_ord(source[i]),
                                        _ord(source[i+1]),
                                        _ord(source[i+2]),
                                        GeneticCode<CODE_SPEC>());
    }
}

}

#endif // INCLUDE_SEQAN_TRANSLATION_TRANSLATION_H_
