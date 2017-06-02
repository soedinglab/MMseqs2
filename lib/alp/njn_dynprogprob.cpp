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

File name: njn_dynprogprob.cpp

Author: John Spouge

Contents: 

******************************************************************************/

#include <assert.h>

#include "njn_dynprogprob.hpp"
#include "njn_memutil.hpp"

using namespace Njn;


const size_t DynProgProb::ARRAY_CAPACITY = 256;

void DynProgProb::init (size_t arrayCapacity_) // range for d_array_p [0,1][0...arrayCapacity_ - 1]
{
    if (arrayCapacity_ > 0) 
    {
	for (size_t i = 0; i < 2; i++) 
	{
	    d_array_p [i] = new double [arrayCapacity_];
	    MemUtil::memZero (d_array_p [i], sizeof (double) * arrayCapacity_);
	}
    }

    d_arrayCapacity = arrayCapacity_;
}

void DynProgProb::free2 ()
{
    if (getArrayCapacity () > 0) 
    {
	for (size_t i = 0; i < 2; i++) 
	{
	    delete [] d_array_p [i]; d_array_p [i] = 0;
	}
    }

    d_arrayCapacity = 0;
}

void DynProgProb::clear (
long int valueBegin_, // lower limit for long int values in the array (an offset)
size_t arrayCapacity_) // new array capacity 
// resets the computation
{
    free2 ();
    init (arrayCapacity_);
    d_valueBegin = valueBegin_;
    d_step = 0;
}

void DynProgProb::clear ( // initializes the "probabilities" with non-negative weights
long int valueLower_, // lower long int value corresponding to the "probability" array
long int valueUpper_, // one beyond present upper long int value corresponding to the "probability" array
const double *prob_) // "probabilities" prob [valueLower_, valueUpper_) corresponding to the long ints
{
	assert ((! prob_ && valueLower_ <= 0 && 0 <= valueUpper_) || 
		/*sls added "("*/(prob_ && valueLower_ < valueUpper_)/*sls added ")"*/ );

    if (prob_) 
    {
	for (size_t i = 0; i < static_cast <size_t> (valueUpper_ - valueLower_); i++) 
	{
	    assert (0.0 <= prob_ [i]);
	}

	clear (valueLower_, static_cast <size_t> (valueUpper_ - valueLower_));

	d_valueLower = valueLower_;
	d_valueUpper = valueUpper_;
	MemUtil::memCpy (d_array_p [0], prob_, sizeof (double) * getArrayCapacity ());

	return;
    }

    if (valueLower_ == 0 && valueUpper_ == 0) 
    {
	clear (-static_cast <long int> (ARRAY_CAPACITY / 2) + 1, ARRAY_CAPACITY);
    } 
    else 
    {
	clear (valueLower_, static_cast <size_t> (valueUpper_ - valueLower_));
    }

    d_valueLower = 0;
    d_valueUpper = 1;
    d_array_p [0][getArrayPos (0)] = 1.0;

    return;
}

void DynProgProb::copy (
size_t step_, // current index : starts at 0 
const double *const *array_, // two corresponding arrays of probabilities 
size_t arrayCapacity_, // present capacity of the array
long int valueBegin_, // lower limit for long int values in the array (an offset)
long int valueLower_, // present lower long int value in the array
long int valueUpper_, // one beyond present upper long int value in the array
ValueFct *newStateFct_, // function for updating dynamic programming values
size_t dimInputProb_, 
const double *inputProb_) // array of input states : d_inputProb_p [0...dimStateProb - 1]
{
    if (arrayCapacity_ != getArrayCapacity ()) 
    {
	free2 ();
	init (arrayCapacity_);
    }

    d_step = step_;
    for (size_t i = 0; i < 2; i++) 
    {
	if (getArrayCapacity () > 0) MemUtil::memCpy (d_array_p [i], array_ [i], sizeof (double) * getArrayCapacity ());
    }

    d_valueBegin = valueBegin_;
    d_valueLower = valueLower_;
    d_valueUpper = valueUpper_;

    setValueFct (newStateFct_);
    setInput (dimInputProb_, inputProb_);
}

void DynProgProb::initInput (size_t dimInputProb_) // array of input states : d_inputProb_p [0...dimStateProb - 1]
{
    if (dimInputProb_ > 0) 
    {
	d_inputProb_p = new double [dimInputProb_];
	MemUtil::memZero (d_inputProb_p, sizeof (double) * dimInputProb_);
    }

    d_dimInputProb = dimInputProb_;
}

void DynProgProb::freeInput ()
{
    if (getDimInputProb () > 0) 
    {
	delete [] d_inputProb_p; d_inputProb_p = 0;
    }

    d_dimInputProb = 0;
}

void DynProgProb::setInput (
size_t dimInputProb_, 
const double *inputProb_) // array of input states : d_inputProb_p [0...dimStateProb - 1]
{
    if (dimInputProb_ != getDimInputProb ()) 
    {
	freeInput ();
	initInput (dimInputProb_);
    }

    if (getDimInputProb () > 0) MemUtil::memCpy (d_inputProb_p, inputProb_, sizeof (double) * getDimInputProb ());
}

void DynProgProb::update () // updates dynamic prog probs 
{
    assert (getValueFct ());
    assert (getDimInputProb ());
    assert (getInputProb ());

    const size_t ARRAY_FAC = 2;
    assert (1 < ARRAY_FAC);

    long int i = 0;
    size_t j = 0;

    const double *oldArray = 0;
    double *array = 0;
    long int value = 0;
    long int valueBegin = 0;
    long int valueLower = 0;
    long int valueUpper = 0;
    /*sls deleted size_t arrayPos = 0;*/

    oldArray = d_array_p [d_step % 2];
    array = d_array_p [(d_step + 1) % 2];
    valueLower = LONG_MAX;
    valueUpper = LONG_MIN;

    MemUtil::memZero (array, sizeof (double) * getArrayCapacity ());

    for (i = getValueLower (); i < getValueUpper (); i++) 
    {
	if (oldArray [getArrayPos (i)] == 0.0) continue;

	for (j = 0; j < getDimInputProb (); j++) 
	{
	    if (getInputProb () [j] == 0.0) continue;

	    // adjust the reserve, if necessary	    
	    value = getValueFct () (i, j);
	    while (value < getValueBegin () || getValueEnd () <= value) {
	    valueBegin = getValueBegin ();
	    if (value < getValueBegin ()) valueBegin -= (ARRAY_FAC - 1) * getArrayCapacity ();
	    reserve (ARRAY_FAC * getArrayCapacity ());
	    setValueBegin (valueBegin);
	    oldArray = d_array_p [d_step % 2];
	    array = d_array_p [(d_step + 1) % 2];
	    }

	    if (value < valueLower) valueLower = value;
	    if (valueUpper < value) valueUpper = value;

	    // add the probability
	    assert (getValueBegin () <= i);
	    assert (i < getValueEnd ());
	    array [getArrayPos (value)] += oldArray [getArrayPos (i)] * getInputProb () [j];
	}
    }

    d_valueLower = valueLower;
    d_valueUpper = valueUpper + 1;
    d_step++;
}

void DynProgProb::reserve (size_t arrayCapacity_) // new array capacity
// increases capacity of and copies d_array_p, while updating other variables
{
    assert (getArrayCapacity () < arrayCapacity_);

    double *array = new double [getArrayCapacity ()];

    for (size_t i = 0; i < 2; i++) 
    {
	MemUtil::memCpy (array, d_array_p [i], sizeof (double) * getArrayCapacity ());
	delete [] d_array_p [i]; d_array_p [i] = 0;
	d_array_p [i] = new double [arrayCapacity_];
	MemUtil::memZero (d_array_p [i], sizeof (double) * arrayCapacity_);
	MemUtil::memCpy (d_array_p [i], array, sizeof (double) * getArrayCapacity ());
    }

    d_arrayCapacity = arrayCapacity_;

    delete [] array; array = 0;
}

void DynProgProb::setValueBegin (long int valueBegin_)
// resets the offset d_valueBegin
// assert (valueBegin_ <= getValueBegin ()) : enlarge the array only 
// assert (offSet < getArrayCapacity ());
{
    assert (valueBegin_ <= getValueBegin ());
    size_t offSet = static_cast <size_t> (getValueBegin () - valueBegin_);
    if (offSet == 0) return; // nothing to do

    assert (offSet < getArrayCapacity ());

    double *array = new double [getArrayCapacity ()];

    for (size_t i = 0; i < 2; i++) 
    {
	MemUtil::memCpy (array, d_array_p [i], sizeof (double) * getArrayCapacity ());
	MemUtil::memZero (d_array_p [i], sizeof (double) * getArrayCapacity ());
	MemUtil::memCpy (d_array_p [i] + offSet, array, sizeof (double) * (getArrayCapacity () - offSet));
    }

    delete [] array; array = 0;

    d_valueBegin = valueBegin_;
}

