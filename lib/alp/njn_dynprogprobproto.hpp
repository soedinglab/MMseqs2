#ifndef INCLUDED_NJN_DYNPROGPROBPROTO
#define INCLUDED_NJN_DYNPROGPROBPROTO

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

File name: njs_dynprogprobproto.hpp

Author: John Spouge

Contents: 

******************************************************************************/

#include <vector>
#include <cstdlib>

namespace Njn {

    class	DynProgProbProto {

        // DynProgProbProto performs updates for probabilities in a dynamic programming computation.

        // The object behaves as follows:
        //
        // Default:
        //    (1) The initial state of the dynamic programming computation is 0 with probability 1.0.
        //
        // Behavior:
        //    (2) If input_ is the computation's input, it replaces oldValue_ with NewValueFct (oldValue_, input_)
        //    (3) The dynamic programming function can be reset with setValueFct (...).
        //    (4) The probability for the input can be reset with setInput (...).
        //    (5) The probability of input_ = [0, dimInputProb_) is inputProb_ [input_].
        //    (6) getProb (long int i_) returns the probability corresponding to the long int value i_.

        public:

        typedef long int ValueFct (long int oldValue_, size_t input_); 
        // function for producing new value from old value and input state

        virtual ~DynProgProbProto ();

        virtual inline operator bool () const = 0; // ? is the object ready for computation ?

        virtual void clear () = 0; // restarts the computation

        virtual void update () = 0; // updates dynamic prog probs 

        virtual inline double getProb (long int value_) const = 0; // probability value

        virtual inline size_t getStep () const = 0; // current index : starts at 0 

        virtual inline long int getValueLower () const = 0; // present lower long int value in the array
        virtual inline long int getValueUpper () const = 0; // one beyond present upper long int value in the array
    };

}

#endif  

