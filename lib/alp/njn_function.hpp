#ifndef INCLUDED_NJN_FUNCTION
#define INCLUDED_NJN_FUNCTION

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

File name: njn_function.hpp

Author: John Spouge

Contents: 

******************************************************************************/
#include "njn_doubletype.hpp"


namespace Njn {
	namespace Function {


      template <typename T> inline T bitsToNats (T x_);
      template <typename T> inline T natsToBits (T x_);
      template <typename T> inline T hartleysToNats (T x_);
      template <typename T> inline T natsToHartleys (T x_);
      template <typename T> inline T hartleysToBits (T x_);
      template <typename T> inline T bitsToHartleys (T x_);

      template <typename T> inline T exp2 (T x_);
      template <typename T> inline T log2 (T x_);
      template <typename T> inline T exp10 (T x_);
      template <typename T> inline T log10 (T x_);

      template <typename T> inline T max (T x_, T y_);
      template <typename T> inline T max3 (T a_, T b_, T c_);
      template <typename T> inline T min (T x_, T y_);
      template <typename T> inline T plus (T x_); // x_^+
      template <typename T> inline T minus (T x_); // x_^-
      template <typename T> inline T signum (T x_);
      template <typename T> inline T heaviside (T x_);

      template <typename T> inline T psqrt (T x_); // square-root of x_^+
      template <typename T> inline T square (T x_); // x_ * x_

      template <typename T> inline T bound (T x_, T xlo_, T xhi_); // the closest value to x_ in the interval [xlo_, xhi_]
      template <typename T> inline T probability (T x_); // the closest value to x_ in the interval [0.0, 1.0]
      template <typename T> inline T prob (T x_); // the closest value to x_ in the interval [0.0, 1.0]

      template <typename T> inline bool isProb (T x_); // the closest value to x_ in the interval [0.0, 1.0]


//
// There are no more declarations beyond this point.
//

      template <typename T> T bitsToNats (T x_) {return x_ * DoubleType::LN_2;}
      template <typename T> T natsToBits (T x_) {return x_ / DoubleType::LN_2;}
      template <typename T> T hartleysToNats (T x_) {return x_ * DoubleType::LN_10;}
      template <typename T> T natsToHartleys (T x_) {return x_ / DoubleType::LN_10;}
      template <typename T> T hartleysToBits (T x_) {return hartleysToNats (natsToBits (x_));}
      template <typename T> T bitsToHartleys (T x_) {return bitsToNats (natsToHartleys (x_));}

      template <typename T> T exp2 (T x_) {return exp (DoubleType::LN_2 * (x_));}
      template <typename T> T log2 (T x_) {return log (x_) / DoubleType::LN_2;}
      template <typename T> T exp10 (T x_) {return exp (DoubleType::LN_10 * (x_));}
      template <typename T> T log10 (T x_) {return log (x_) / DoubleType::LN_10;}

      template <typename T> T max (T x_, T y_) {return x_ > y_ ? x_ : y_;}
      template <typename T> T max3 (T a_, T b_, T c_) {return max <T> (a_, max <T> (b_, c_));}
      template <typename T> T min (T x_, T y_) {return x_ < y_ ? x_ : y_;}
      template <typename T> T positive (T x_) {return max <T> (x_, static_cast <T> (0));}
      template <typename T> T negative (T x_) {return max <T> (-x_, static_cast <T> (0));}
      template <typename T> T signum (T x_) {return x_ > 0.0 ? 1.0 : (x_ == 0.0 ? 0.0 : -1.0);}
      template <typename T> T heaviside (T x_) {return x_ < 0.0 ? 0.0 : 1.0;}

      template <typename T> T psqrt (T x_) {return positive <T> (sqrt (x_));}
      template <typename T> T square (T x_) {return x_ * x_;}

      template <typename T> T bound (T x_, T xlo_, T xhi_) {
         if (x_ < xlo_) return xlo_;
         if (x_ < xhi_) return x_;
         return xhi_;
      }

      template <typename T> T probability (T x_) {return bound <T> (x_, 0.0, 1.0);}
      template <typename T> T prob (T x_) {return probability <T> (x_);}

      template <typename T> bool isProb (T x_) {return 0.0 <= x_ && x_ <= 1.0;}/*sls deleted;*/ // the closest value to x_ in the interval [0.0, 1.0]

		}
	}

#endif //! INCLUDED

