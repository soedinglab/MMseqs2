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

File name: njn_localmaxstat.cpp

Author: John Spouge

Contents: 

******************************************************************************/

#include "sls_basic.hpp"
#include "njn_localmaxstat.hpp"
#include "njn_memutil.hpp"
#include "njn_dynprogproblim.hpp"
#include "njn_function.hpp"
#include "njn_integer.hpp"
#include "njn_localmaxstatutil.hpp"

using namespace Njn;

double LocalMaxStat::s_time = 0.0;

void LocalMaxStat::init (size_t dimension_)
{
    if (dimension_ > 0) 
    {
        d_score_p = new long int [dimension_]; 
        d_prob_p = new double [dimension_]; 
    }

    d_dimension = dimension_;
}

void LocalMaxStat::free2 ()
{
    if (getDimension () > 0) 
    {
        delete [] d_score_p; d_score_p = 0; 
        delete [] d_prob_p; d_prob_p = 0; 
    }

    d_dimension = 0;
}

void LocalMaxStat::clear ()
{
    free2 ();
    init (0);

    d_lambda = 0.0;
    d_k = 0.0;
    d_c = 0.0;
    d_thetaMin = 0.0;
    d_rMin = 0.0;
    d_delta = 0; 
    d_thetaMinusDelta = 0.0;
    d_mu = 0.0;
    d_sigma = 0.0;
    d_muAssoc = 0.0;
    d_sigmaAssoc = 0.0;
    d_meanWDLE = 0.0;
    d_terminated = false;
}

void LocalMaxStat::copy (
size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *score_, // scores 
const double *prob_, // probabilities
double lambda_, // lambda for associated random walk
double k_, // k for random walk : exponential prefactor
double c_, // c for random walk : exponential prefactor (global alignment)
double thetaMin_, // theta for minimum expectation (exp (theta * score))
double rMin_, // minimum expectation (exp (theta * score))
long int delta_, // span 
double thetaMinusDelta_, // renewal span parameter
double mu_, // n_step mean for random walk
double sigma_, // n_step standard deviation for random walk
double muAssoc_, // n_step mean for associated random walk (relative entropy)
double sigmaAssoc_, // n_step standard deviation for associated random walk
double meanLength_, // expected renewal length
bool terminated_) // ? Was the dynamic programming computation terminated prematurely ?
{
    free2 ();
    init (dimension_);

    memcpy (d_score_p, score_, sizeof (long int) * getDimension ());
    memcpy (d_prob_p, prob_, sizeof (double) * getDimension ());

    d_lambda = lambda_;
    d_k = k_;
    d_c = c_;
    d_thetaMin = thetaMin_;
    d_rMin = rMin_;
    d_delta = delta_; 
    d_thetaMinusDelta = thetaMinusDelta_;
    d_mu = mu_;
    d_sigma = sigma_;
    d_muAssoc = muAssoc_;
    d_sigmaAssoc = sigmaAssoc_;
    d_meanWDLE = meanLength_;
    d_terminated = terminated_;
}

void LocalMaxStat::copy (
size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *score_, // scores in increasing order
const double *prob_) // corresponding probabilities
{
    if (dimension_ == 0) 
    {
        clear ();
        return;
    }

    if (! LocalMaxStatUtil::isLogarithmic (dimension_, score_, prob_))
    {
        //IoUtil::abort ("LocalMaxStat::copy : ! isLogarithmic");
		throw Sls::error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
    }

    size_t i = 0;
    /*sls deleted size_t j = 0;*/
    /*sls deleted long int iter = 0;*/
    /*sls deleted long int value = 0;*/

    free2 ();
    init (dimension_);

    memcpy (d_score_p, score_, sizeof (long int) * getDimension ());
    memcpy (d_prob_p, prob_, sizeof (double) * getDimension ());

    d_mu = LocalMaxStatUtil::mu (getDimension (), getScore (), getProb ());
    d_sigma = 0.0;

    for (i = 0; i < dimension_; i++) 
    {
        d_sigma += static_cast <double> (score_ [i]) * static_cast <double> (score_ [i]) * prob_ [i];
    }

    d_sigma -= getMu () * getMu ();
    d_sigma = Function::psqrt (getSigma ());

    // calculate lambda

    d_lambda = LocalMaxStatUtil::lambda (getDimension (), getScore (), getProb ());
    d_muAssoc = LocalMaxStatUtil::muAssoc (getDimension (), getScore (), getProb (), getLambda ());
    d_sigmaAssoc = 0.0;

    for (i = 0; i < getDimension (); i++) 
    {
        d_sigmaAssoc += static_cast <double> (getScore () [i]) * static_cast <double> (getScore () [i]) * 
            getProb () [i] * exp (getLambda () * static_cast <double> (getScore () [i]));
    }

    d_sigmaAssoc -= getMuAssoc () * getMuAssoc ();
    d_sigmaAssoc = Function::psqrt (d_sigmaAssoc);

    d_thetaMin = LocalMaxStatUtil::thetaMin (getDimension (), getScore (), getProb (), getLambda ());
    d_rMin = LocalMaxStatUtil::rMin (getDimension (), getScore (), getProb (), getLambda (), getThetaMin ());

    d_delta = LocalMaxStatUtil::delta (getDimension (), getScore ());
    d_thetaMinusDelta = LocalMaxStatUtil::thetaMinusDelta (getLambda (), getDimension (), getScore ());

    dynProgCalc ();
}

double LocalMaxStat::getR (double theta_) const 
{
    return LocalMaxStatUtil::r (d_dimension, d_score_p, d_prob_p, theta_);
}

void LocalMaxStat::dynProgCalc ()
// k for random walk : exponential prefactor 
// expected renewal length for weak ladder epochs
{
    double eSumAlpha_ = 0.0;
    double eOneMinusExpSumAlpha_ = 0.0;
    LocalMaxStatUtil::descendingLadderEpoch (getDimension (), getScore (), getProb (), 
        &eSumAlpha_, &eOneMinusExpSumAlpha_, false, 
        getLambda (), getMu (), getMuAssoc (), getThetaMin (), getRMin (), getTime (), &d_terminated);

    if (getTerminated ()) return;

    // fluctuation sum quantities
    double ratio = eOneMinusExpSumAlpha_ / eSumAlpha_;
    d_k = getMu () * getMu () / getThetaMinusDelta () / getMuAssoc () * ratio * ratio;
    d_meanWDLE = eSumAlpha_ / getMu ();
    d_c = getK () * getMeanWDLE () / eOneMinusExpSumAlpha_;
}

