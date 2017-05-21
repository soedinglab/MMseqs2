#ifndef INCLUDED_NJN_LOCALMAXSTATUTIL
#define INCLUDED_NJN_LOCALMAXSTATUTIL

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

File name: njn_localmaxstatutil.hpp

Author: John Spouge

Contents: Random walk parameters

******************************************************************************/

#include "njn_matrix.hpp"


namespace Njn {
	namespace LocalMaxStatUtil {


        const double REL_TOL = 1.0e-6;

        void flatten ( // allocates memory for linear probabilities and scores
        size_t dimension_, // dimension of equilProb_
        const long int *const *scoreMatrix_, // packed scoring matrix [0...dimension_)[0...dimension_)
        const double *const *prob_, // prob_ [0...dimension_)[0...dimension_) : distribution of scores sum to 1.0
        size_t *dim_, // dimension of p_
        long int **score_, // score [0...dim_) in increasing order
        double **p_, // linear p_ [0...dim_) : distribution of scores
        size_t dimension2_ = 0); // dimension2 of equilProb_ : defaults to dimension_
        // asserts (sum (p_) == 1.0);

        double lambda (
        size_t dimMatrix_, // dimension of equilProb_
        const long int *const *scoreMatrix_, // packed scoring matrix [0...dimension_)[0...dimension_)
        const double *q_); // q_ [0...dimension_) : distribution of independent letters

        double mu (
        size_t dimension_, // #(distinct values)          
        const long int *score_, // scores in increasing order
        const double *prob_); // probability of corresponding value  

        double lambda (
        size_t dimension_, // #(distinct values)          
        const long int *score_, // scores in increasing order
        const double *prob_); // probability of corresponding value  
        // assumes logarithmic regime

        double muAssoc (
        size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
        const long int *score_, // scores in increasing order
        const double *prob_, // corresponding probabilities
        double lambda_ = 0.0); // lambda

        double muPowerAssoc (
        size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
        const long int *score_, // scores in increasing order
        const double *prob_, // corresponding probabilities
        double lambda_ = 0.0, // lambda
        long int power_ = 1); // power

        double thetaMin ( // minimizing value for r(theta)
        size_t dimension_, // #(distinct values)          
        const long int *score_, // scores in increasing order
        const double *prob_, // probability of corresponding value 
        double lambda_ = 0.0); // lambda
        // assumes logarithmic regime

        double rMin ( // minimum value of r(theta)
        size_t dimension_, // #(distinct values)          
        const long int *score_, // scores in increasing order
        const double *prob_, // probability of corresponding value 
        double lambda_ = 0.0, // lambda
        double thetaMin_ = 0.0); // argument of rate
        // assumes logarithmic regime

        double r ( // r(theta)
        size_t dimension_, // #(distinct values)          
        const long int *score_, // scores in increasing order
        const double *prob_, // probability of corresponding value 
        double theta_); // argument of rate
        // assumes logarithmic regime

        long int delta ( // theta [minus delta] for ungapped sequence comparison
        size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
        const long int *score_); // scores 

        double thetaMinusDelta ( // theta [minus delta] for ungapped sequence comparison
        double lambda_, // lambda, the exponential rate for the local maximum         
        size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
        const long int *score_); // scores 

        void descendingLadderEpoch (
        size_t dimension_, // #(distinct values)          
        const long int *score_, // values 
        const double *prob_, // probability of corresponding value 
        double *eSumAlpha_ = 0, // expectation (sum [alpha])
        double *eOneMinusExpSumAlpha_ = 0, // expectation [1.0 - exp (sum [alpha])]
        bool isStrict_ = false, // ? is this a strict descending ladder epoch
        double lambda0_ = 0.0, // lambda for flattened distribution (avoid recomputation)
        double mu0_ = 0.0, // mean of flattened distribution (avoid recomputation)
        double muAssoc0_ = 0.0, // mean of associated flattened distribution (avoid recomputation)
        double thetaMin0_ = 0.0, // thetaMin of flattened distribution (avoid recomputation)
        double rMin0_ = 0.0, // rMin of flattened distribution (avoid recomputation)
        double time_ = 0.0, // get time for the dynamic programming computation
        bool *terminated_ = 0); // ? Was the dynamic programming computation terminated prematurely ?

        void descendingLadderEpochRepeat (
        size_t dimension_, // #(distinct values)          
        const long int *score_, // values 
        const double *prob_, // probability of corresponding value 
        double *eSumAlpha_ = 0, // expectation (sum [alpha])
        double *eOneMinusExpSumAlpha_ = 0, // expectation [1.0 - exp (sum [alpha])]
        bool isStrict_ = false, // ? is this a strict descending ladder epoch
        double lambda_ = 0.0, // lambda for repeats : default is lambda0_ below
        size_t endW_ = 0, // maximum w plus 1
        double *pAlphaW_ = 0, // probability {alpha = w} : pAlphaW_ [0, wEnd)
        double *eOneMinusExpSumAlphaW_ = 0, // expectation [1.0 - exp (sum [alpha]); alpha = w] : eOneMinusExpSumAlphaW_ [0, wEnd)
        double lambda0_ = 0.0, // lambda for flattened distribution (avoid recomputation)
        double mu0_ = 0.0, // mean of flattened distribution (avoid recomputation)
        double muAssoc0_ = 0.0, // mean of associated flattened distribution (avoid recomputation)
        double thetaMin0_ = 0.0, // thetaMin of flattened distribution (avoid recomputation)
        double rMin0_ = 0.0, // rMin of flattened distribution (avoid recomputation)
        double time_ = 0.0, // get time for the dynamic programming computation
        bool *terminated_ = 0); // ? Was the dynamic programming computation terminated prematurely ?
        // assumes logarithmic regime

        bool isProbDist (
        size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
        const double *prob_); // corresponding probabilities

        bool isScoreIncreasing (
        size_t dimension_, // #(distinct values)          
        const long int *score_); // scores in increasing order

        bool isLogarithmic (
        size_t dimension_, // #(distinct values)          
        const long int *score_, // scores in increasing order 
        const double *prob_); // probability of corresponding value  

	}
}

#endif  

