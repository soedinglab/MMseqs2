/* $Id: $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
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
* ===========================================================================*/

/*****************************************************************************

File name: njn_random.cpp

Author: John Spouge

Contents: Additive random number generator

//	Modelled after "Algorithm A" in
//	D.E. Knuth (1981) The art of computer programming, volume 2, page 27.

//	7/26/90 Warren Gish

******************************************************************************/



#include "njn_random.hpp"
#include <cstring>

using namespace Njn;


namespace {

	const size_t r_off = 12;

	long	state [33] = {
	static_cast <long> (0xd53f1852),  static_cast <long> (0xdfc78b83),  static_cast <long> (0x4f256096),  static_cast <long> (0xe643df7),
	static_cast <long> (0x82c359bf),  static_cast <long> (0xc7794dfa),  static_cast <long> (0xd5e9ffaa),  static_cast <long> (0x2c8cb64a),
	static_cast <long> (0x2f07b334),  static_cast <long> (0xad5a7eb5),  static_cast <long> (0x96dc0cde),  static_cast <long> (0x6fc24589),
	static_cast <long> (0xa5853646),  static_cast <long> (0xe71576e2),  static_cast <long> (0xdae30df),   static_cast <long> (0xb09ce711),
	static_cast <long> (0x5e56ef87),  static_cast <long> (0x4b4b0082),  static_cast <long> (0x6f4f340e),  static_cast <long> (0xc5bb17e8),
	static_cast <long> (0xd788d765),  static_cast <long> (0x67498087),  static_cast <long> (0x9d7aba26),  static_cast <long> (0x261351d4),
	static_cast <long> (0x411ee7ea),  static_cast <long> (0x393a263),   static_cast <long> (0x2c5a5835),  static_cast <long> (0xc115fcd8),
	static_cast <long> (0x25e9132c),  static_cast <long> (0xd0c6e906),  static_cast <long> (0xc2bc5b2d),  static_cast <long> (0x6c065c98),
	static_cast <long> (0x6e37bd55)};

   long	*rJ = &state [r_off];
	long	*rK = &state [sizeof state / sizeof *state - 1];

}
void Random::seed (long x)
{
	register size_t i;

	state [0] = x;
   
   // linear congruential initializer
	for (i = 1; i < sizeof state / sizeof *state; ++i) {
		state [i] = 1103515245 * state [i - 1] + 12345;
	}

	rJ = &state [r_off];
	rK = &state [sizeof state / sizeof *state - 1];

	for (i = 0; i < 10 * sizeof state / sizeof *state; ++i) number ();
}

long Random::number () // uniform random x : 0 <= x <= exp2 (31) - 1

{
	register long	r;

	r = *rK;
	r += *rJ--;
	*rK-- = r;

   if (rK < state) {
		rK = &state [sizeof state / sizeof *state - 1];
   } else if (rJ < state) {
			rJ = &state [sizeof state / sizeof *state - 1];
   }
	
   return (r >> 1) &0x7fffffff; // discard the least-random bit
}

