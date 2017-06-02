#ifndef INCLUDED_NJN_UNIFORM
#define INCLUDED_NJN_UNIFORM

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

File name: njn_uniform.hpp

Author: John Spouge

Contents: 

******************************************************************************/



#include "njn_random.hpp"

#include <assert.h>


namespace Njn {

   namespace Uniform { // uniform distribution 

      template <typename T>
      inline T variate (T a_ = static_cast <T> (0.0), T b_= static_cast <T> (1.0)); // returns uniform random variate [a_, b_) 

      template <typename T>
      inline T standardVariate ();

//
// There are no more declarations beyond this point.
//
      template <typename T>
      T variate (T a_, T b_)
      {
         assert (a_ != b_);

         //if (b_ < a_) std::swap <T> (a_, b_); // a_ < b_
		 if (b_ < a_) std::swap (a_, b_); // a_ < b_/*sls deleted <T>*/

         long random = 0;
         while ((random = Njn::Random::number ()) == 0x7fffffff);

         return a_ + static_cast <T> ((b_ - a_) * static_cast <double> (Random::number ()) / static_cast <double> (0x7fffffff)); 
      }
   }
}

#endif //! INCLUDED

