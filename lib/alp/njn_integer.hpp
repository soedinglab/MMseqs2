#ifndef INCLUDED_NJN_INTEGER
#define INCLUDED_NJN_INTEGER

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

File name: njn_integer.hpp

Author: John Spouge

Contents: 

******************************************************************************/

#include "njn_ioutil.hpp"

namespace Njn {
	namespace Integer {


        template <class Int>
        Int minimum (
        Int i, // number being divided
        Int j) // divisor 
        {
            return i < j ? i : j;
        }

        template <class Int>
        Int maximum (
        Int i, // number being divided
        Int j) // divisor 
        {
            return i > j ? i : j;
        }

        template <class Int>
        Int mod ( // Portable modular function %
        Int i, // number being divided
        Int j) // divisor 
        {
            Int abs_j = j >= 0 ? j : -j;

            if (j == 0) 
            {
                Njn::IoUtil::abort ("Nks_Mod : j == 0");   
            }

            if (i < 0) 
            {
                Int k = (i >= 0 ? i : -i) % abs_j;
                return k == 0 ? k : abs_j - k;
            } 
            else 
            {
                return i % abs_j;
            }
        }

        template <class Int>
        Int euclidAlgorithm ( // Euclid's algorithm returning the greatest common divisor (i, j)
        Int i_,
        Int j_)
        {
            Int abs_i = maximum <Int> (i_ >= 0 ? i_ : -i_, j_ >= 0 ? j_ : -j_);
            Int abs_j = minimum <Int> (i_ >= 0 ? i_ : -i_, j_ >= 0 ? j_ : -j_);
            Int remainder = 0;

            while (0 < abs_j) 
            {
                remainder = mod <Int> (abs_i, abs_j);
                abs_i = abs_j;
                abs_j = remainder;
            }

            return abs_i;
        }

        template <class Int, class ConstIntegerPtr>
        Int euclidAlgorithmVector ( // Euclid's algorithm returning the greatest common divisor (i, j)
        ConstIntegerPtr begin_,
        ConstIntegerPtr end_)
        {
            assert (begin_ <= end_);

            if (begin_ == end_) return 0;

            Int gcd = *begin_;

            for (ConstIntegerPtr k = begin_ + 1; k != end_; k++) 
            {
                gcd = euclidAlgorithm (gcd, *k);
            }

            return gcd;
        }

        template <class Int>
        Int minusOnePower ( // (-1)^i
        Int i) // exponent
        {
            return mod <Int> (i, 2) == 0 ? 1 : -1;
        }

        template <class Real, class Int>
        Real integerPower (Real x, Int n) // integer power function
        {
            if (x == 0.0) 
            {
                if (n < 0) 
                {
                    Njn::IoUtil::abort ("Int::integerPower <class Real, class Int> : negative exponent of zero");
                } 
                else if (n == 0) 
                {
                    return 1.0;
                } 
                else 
                {
                    return 0.0;
                }
            }

            Real y = 1.0;
            Int i = 0;

            for (i = n > 0 ? n : -n; i > 0 ; i /= 2) 
            {
                if (i % 2 != 0) 
                {
                    y *= x;
                }

                x *= x;
            }

            if (n < 0) 
            {
                y = 1.0 / y;
            }
    
            return y;
        }

        template <class Real, class UnsignedInteger>
        Real integerPositivePower (Real x, UnsignedInteger n) // integer power function
        {
            if (x == 0.0) 
            {
                if (n == 0) 
                {
                    return 1.0;
                } 
                else 
                {
                    return 0.0;
                }
            }

            Real y = 1.0;
            UnsignedInteger i = 0;

            for (i = n; i > 0 ; i /= 2) 
            {
                if (i % 2 != 0) 
                {
                    y *= x;
                }

                x *= x;
            }
            return y;
        }

        template <class Int, class IntTwo>
        Int intPower (Int x, IntTwo n) // integer power function
        {
            assert (n >= 0);

            if (x == 0) 
            {
                return n == 0 ? 1 : 0;
            }

            Int y = 1;
            IntTwo i = 0;

            for (i = n > 0 ? n : -n; i > 0 ; i /= 2) 
            {
                if (i % 2 != 0)
                {
                    y *= x;
                }

                x *= x;
            }
            return y;
        }

	}
}

#endif 

