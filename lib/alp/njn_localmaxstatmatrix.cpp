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

File name: njn_localmaxstatmatrix.cpp

Author: John Spouge

Contents: 

******************************************************************************/

#include "njn_localmaxstatmatrix.hpp"
#include "njn_localmaxstatutil.hpp"
#include "njn_memutil.hpp"


using namespace Njn;

void LocalMaxStatMatrix::init (size_t dimMatrix_, size_t dimMatrix2_)
{
    if (dimMatrix2_ == 0) dimMatrix2_ = dimMatrix_;

    if (dimMatrix_ > 0 && dimMatrix2_ > 0) 
    {
        d_scoreMatrix_p = MemUtil::newMatrix <long int> (dimMatrix_, dimMatrix2_); 
        d_p_p = new double [dimMatrix_]; 
        d_p2_p = new double [dimMatrix2_]; 
    }

    d_dimMatrix = dimMatrix_;
    d_dimMatrix2 = dimMatrix2_;
}

void LocalMaxStatMatrix::free2 ()
{
    if (getDimMatrix () > 0 && getDimMatrix2 () > 0) 
    {
        MemUtil::deleteMatrix <long int> (d_scoreMatrix_p, getDimMatrix (), getDimMatrix2 ());
        d_scoreMatrix_p = 0;
        delete [] d_p_p; d_p_p = 0; 
        delete [] d_p2_p; d_p2_p = 0; 
    }

    d_dimMatrix = 0;
    d_dimMatrix2 = 0;
}

void LocalMaxStatMatrix::copy (
size_t dimMatrix_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *const *scoreMatrix_, // score matrix [0...dimMatrix_)[0...dimMatrix_)
const double *p_, // probability of "letters" p_ [0...dimMatrix_)
const double *p2_, // probability of "letters" p2_ [0...dimMatrix_), the second (j) set of letter-probabilities
size_t dimMatrix2_) // #(distinct letters) in the second alphabet         
{
    if (! p2_) p2_ = p_;
    if (dimMatrix2_ == 0) dimMatrix2_ = dimMatrix_;

    free2 ();
    init (dimMatrix_, dimMatrix2_);

    if (getDimMatrix () == 0) 
    {
        LocalMaxStat::copy (0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        return;
    }

    size_t i = 0;
    size_t j = 0;

    d_dimMatrix = dimMatrix_;
    d_dimMatrix2 = dimMatrix2_;

    for (i = 0; i < getDimMatrix (); i++) 
    {
        memcpy (d_scoreMatrix_p [i], scoreMatrix_ [i], sizeof (long int) * getDimMatrix2 ());
    }

    memcpy (d_p_p, p_, sizeof (double) * getDimMatrix ());
    memcpy (d_p2_p, p2_, sizeof (double) * getDimMatrix2 ());

    size_t dim = 0;
    long int *score = 0;
    double *p = 0;

    double **probMatrix = MemUtil::newMatrix <double> (getDimMatrix (), getDimMatrix2 ());
    
    for (i = 0; i != getDimMatrix (); i++) 
    {
        for (j = 0; j != getDimMatrix2 (); j++) 
        {
            probMatrix [i][j] = p_ [i] * p2_ [j];
        }
    }

    LocalMaxStatUtil::flatten (getDimMatrix (), getScoreMatrix (), probMatrix, &dim, &score, &p, getDimMatrix2 ());
    LocalMaxStat::copy (dim, score, p);

    delete [] p; p = 0;
    delete [] score; score = 0;
    MemUtil::deleteMatrix <double> (probMatrix, getDimMatrix (), getDimMatrix2 ());
}

void LocalMaxStatMatrix::copy (
LocalMaxStat localMaxStat_, // base object 
size_t dimMatrix_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *const *scoreMatrix_, // score matrix [0...dimMatrix_)[0...dimMatrix_)
const double *p_, // probability of "letters" p_ [0...dimMatrix_)
const double *p2_, // probability of "letters" p2_ [0...dimMatrix_), the second (j) set of letter-probabilities
size_t dimMatrix2_) // #(distinct letters) in the second alphabet         
{
    if (! p2_) p2_ = p_;
    if (dimMatrix2_ == 0) dimMatrix2_ = dimMatrix_;

    free2 ();
    init (dimMatrix_);

    size_t i = 0;
    /*sls deleted size_t j = 0;*/

    d_dimMatrix = dimMatrix_;
    d_dimMatrix2 = dimMatrix2_;

    for (i = 0; i < getDimMatrix (); i++) 
    {
        memcpy (d_scoreMatrix_p [i], scoreMatrix_ [i], sizeof (long int) * getDimMatrix2 ());
    }

    memcpy (d_p_p, p_, sizeof (double) * getDimMatrix ());
    memcpy (d_p2_p, p2_, sizeof (double) * getDimMatrix2 ());

    LocalMaxStat::copy (localMaxStat_);
}

