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

File name: njn_dynprogproblim.cpp

Author: John Spouge

Contents: 

******************************************************************************/

#include <assert.h>
#include "njn_dynprogproblim.hpp"
#include "njn_memutil.hpp"


using namespace Njn;


void DynProgProbLim::setLimits ( 
long int valueBegin_, // range for values is [valueBegin_, valueEnd_).
long int valueEnd_) // range for values is [valueBegin_, valueEnd_).
// assert (valueBegin_ < valueEnd_);
{
    assert (valueBegin_ < valueEnd_);

    // take care of lost probability
    long int value = 0;

    for (value = getValueLower (); value < valueBegin_; value++) 
    {
        d_probLost += getProb (value);
    }

    for (value = valueEnd_; value < getValueUpper (); value++) 
    {
        d_probLost += getProb (value);
    }

    size_t arrayCapacity = static_cast <size_t> (valueEnd_ - valueBegin_);

    if (getArrayCapacity () <= arrayCapacity) 
    {
        reserve (arrayCapacity);
        setValueBegin (valueBegin_);
    } 
    else 
    {
        setValueBegin (valueBegin_);
        reserve (arrayCapacity);
    }
}

void DynProgProbLim::reserve (size_t arrayCapacity_) // new array capacity
// reduces capacity of and copies d_array_p from its begin
{
    if (arrayCapacity_ == getArrayCapacity ()) return;

    if (getArrayCapacity () < arrayCapacity_) 
    {
        DynProgProb::reserve (arrayCapacity_);
        return;
    }

    assert (arrayCapacity_ < getArrayCapacity ());

    double *array = new double [getArrayCapacity ()];

    for (size_t i = 0; i < 2; i++) 
    {
        MemUtil::memCpy (array, getArray () [i], sizeof (double) * arrayCapacity_);
        delete [] lgetArray () [i]; lgetArray () [i] = 0;

        lgetArray () [i] = new double [arrayCapacity_];
        MemUtil::memCpy (lgetArray () [i], array, sizeof (double) * arrayCapacity_);
    }

    lgetArrayCapacity () = arrayCapacity_;

    delete [] array; array = 0;
}

void DynProgProbLim::setValueBegin (long int valueBegin_)
// resets the offset d_valueBegin
// assert (offSet < getArrayCapacity ());
// if (getValueBegin () < valueBegin_) the caller must ensure consistency 
{
    if (valueBegin_ <= getValueBegin ()) 
    {
        DynProgProb::setValueBegin (valueBegin_);
        return;
    }

    assert (getValueBegin () < valueBegin_);
    size_t offSet = static_cast <size_t> (valueBegin_ - getValueBegin ());

    double *array = new double [getArrayCapacity ()];

    for (size_t i = 0; i < 2; i++) 
    {
        MemUtil::memCpy (array, getArray () [i], sizeof (double) * getArrayCapacity ());
        MemUtil::memZero (lgetArray () [i], sizeof (double) * getArrayCapacity ());
        
        if (offSet < getArrayCapacity ()) 
        {
            MemUtil::memCpy (lgetArray () [i], array + offSet, sizeof (double) * (getArrayCapacity () - offSet));
        }
    }

    delete [] array; array = 0;

    lgetValueBegin () = valueBegin_;
}

void DynProgProbLim::update () // updates dynamic prog probs 
    {
    assert (getValueFct ());
    assert (getDimInputProb ());
    assert (getInputProb ());
    assert (0 < getArrayCapacity ());

    long int i = 0;
    size_t j = 0;

    const double *oldArray = 0;
    double *array = 0;
    long int value = 0;
    long int valueLower = 0;
    long int valueUpper = 0;
    /*sls deleted size_t arrayPos = 0;*/
    double prob = 0.0;

    oldArray = getArray () [getStep () % 2];
    array = lgetArray () [(getStep () + 1) % 2];
    valueLower = LONG_MAX;
    valueUpper = LONG_MIN;

    MemUtil::memZero (array, sizeof (double) * getArrayCapacity ());
    for (i = getValueLower (); i < getValueUpper (); i++) 
    {
        if (oldArray [getArrayPos (i)] == 0.0) continue;

        for (j = 0; j < getDimInputProb (); j++) 
        {
            if (getInputProb () [j] == 0.0) continue;

            value = getValueFct () (i, j);
            prob = oldArray [getArrayPos (i)] * getInputProb () [j];

            if (value < getValueBegin () || getValueEnd () <= value) 
            {
                d_probLost += prob;
            } 
            else 
            {
                if (value < valueLower) valueLower = value;
                if (valueUpper < value) valueUpper = value;

                // add the probability
                assert (getValueBegin () <= i);
                assert (i < getValueEnd ());
                array [getArrayPos (value)] += prob;
            }
        }
    }

    lgetValueLower () = valueLower;
    lgetValueUpper () = valueUpper + 1;
    lgetStep ()++;
}

