#ifndef INCLUDED_NJN_APPROX
#define INCLUDED_NJN_APPROX

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

File name: njn_approx.hpp

Author: John Spouge

Contents: 

******************************************************************************/

#include <float.h>
#include <cassert>
#include <math.h>


namespace Njn { 
	namespace Approx{

      // Approximation

      const float FLT_THRESHOLD =       10.0F; // Rounding threshold
      const float FLT_ROUND  =          FLT_THRESHOLD * FLT_EPSILON; // Rounding error threshold

      const double DBL_THRESHOLD =      100.0; // Rounding threshold
      const double DBL_ROUND  =         DBL_THRESHOLD * DBL_EPSILON; // Rounding error threshold

      template <typename T> inline bool approx (T x_, T y_, T eps_);
      template <typename T> inline bool relApprox (T x_, T y_, T eps_);
      template <typename T> inline bool absRelApprox (T x_, T y_, T tol_, T rtol_);

      // Rounding

      template <typename T> inline bool eq (T x_, T y_, T round_);
      template <typename T> inline bool ge (T x_, T y_, T round_);
      template <typename T> inline bool gt (T x_, T y_, T round_);
      template <typename T> inline bool ne (T x_, T y_, T round_);
      template <typename T> inline bool lt (T x_, T y_, T round_);
      template <typename T> inline bool le (T x_, T y_, T round_);


      // Approximation

      template <typename T> bool approx (T x_, T y_, T eps_) {return fabs (x_ - y_) <= fabs (eps_);}
      template <typename T> bool relApprox (T x_, T y_, T eps_) {return approx <T> (x_, y_, eps_ * y_);}
      template <typename T> bool absRelApprox (T x_, T y_, T tol_, T rtol_) {return approx <T> (x_, y_, tol_) || relApprox <T> (x_, y_, rtol_);}

      // Rounding

      template <typename T> bool eq (T x_, T y_, T round_) {return relApprox <T> (x_, y_, round_);}
      template <typename T> bool ge (T x_, T y_, T round_) {return (x_ - y_) >= -(fabs (y_)) * round_;}
      template <typename T> bool gt (T x_, T y_, T round_) {return (x_ - y_) >   (fabs (y_)) * round_;}
      template <typename T> bool ne (T x_, T y_, T round_) {return ! eq (x_, y_, round_);}
      template <typename T> bool lt (T x_, T y_, T round_) {return ! ge (x_, y_, round_);}
      template <typename T> bool le (T x_, T y_, T round_) {return ! gt (x_, y_, round_);}

		}
	}


#endif //! INCLUDED

