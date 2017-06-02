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

File name: njn_localmaxstatutil.cpp

Author: John Spouge

Contents: 

******************************************************************************/

#include <assert.h>
#include "njn_localmaxstatutil.hpp"
#include "njn_approx.hpp"
#include "njn_dynprogproblim.hpp"
#include "njn_integer.hpp"
#include "njn_memutil.hpp"
#include "njn_root.hpp"

#include "sls_basic.hpp"

using namespace Njn;


void LocalMaxStatUtil::flatten ( // allocates memory for linear probabilities and scores
size_t dimension_, // dimension of equilProb_
const long int *const *scoreMatrix_, // packed scoring matrix [0...dimension_)[0...dimension2_)
const double *const *prob_, // prob_ [0...dimension_)[0...dimension2_) : distribution of scores sum to 1.0
size_t *dim_, // dimension of p_
long int **score_, // score [0...dim_) in increasing order
double **p_, // linear p_ [0...dim_) : distribution of scores
size_t dimension2_) // dimension2 of equilProb_
{
    if (dimension2_ == 0) dimension2_ = dimension_;

    size_t i = 0;
    size_t j = 0;

    double sum = 0.0;

    for (i = 0; i < dimension_; i++) 
    {
      for (j = 0; j < dimension2_; j++) 
      {
         sum += prob_ [i][j];
      }
    }

    const double FUDGE = 20.0;
    assert (Approx::relApprox (sum, 1.0, FUDGE * REL_TOL));

    long int s = 0;
    long int min = LONG_MAX;
    long int max = LONG_MIN;

    for (i = 0; i < dimension_; i++) 
    {
      for (j = 0; j < dimension2_; j++) 
      {
         if (scoreMatrix_ [i][j] < min) min = scoreMatrix_ [i][j];
         if (max < scoreMatrix_ [i][j]) max = scoreMatrix_ [i][j];
      }
    }

    assert (min <= max);

    size_t dim = static_cast <size_t> (max - min + 1);
    double *p = new double [dim];
    for (i = 0; i < dim; i++) p [i] = 0.0;

    for (i = 0; i < dimension_; i++) 
    {
      for (j = 0; j < dimension2_; j++) 
      {
         p [scoreMatrix_ [i][j] - min] += prob_ [i][j];
      }
    }

    *dim_ = 0;

    for (s = min; s <= max; s++) 
    {
      if (0.0 < p [s - min]) ++*dim_;
    }

    *p_ = new double [*dim_];
    *score_ = new long int [*dim_];
    *dim_ = 0;

    for (s = min; s <= max; s++) 
    {
      if (0.0 < p [s - min]) {
         (*score_) [*dim_] = s;
         (*p_) [*dim_] = p [s - min];
         ++*dim_;
      }
    }

    delete [] p; p = 0;
}

double LocalMaxStatUtil::lambda (
size_t dimension_, // dimension of equilProb_
const long int *const *scoreMatrix_, // packed scoring matrix [0...dimension_)[0...dimension_)
const double *q_) // q_ [0...dimension_) : distribution of independent letters
{
    size_t i = 0;
    size_t j = 0;

   double **prob = MemUtil::newMatrix <double> (dimension_, dimension_);

    for (i = 0; i < dimension_; i++) 
    {
        for (j = 0; j < dimension_; j++) 
        {
            prob [i][j] = q_ [i] * q_ [j];
        }
    }

    size_t dim = 0;
    long int *score = 0;
    double *p = 0;

    flatten (dimension_, scoreMatrix_, prob, &dim, &score, &p);

    for (i = 0; i < dimension_; i++) 
    {
        delete prob [i];
    }

    double lambdaHat = LocalMaxStatUtil::lambda (dim, score, p);

    delete p; p = 0;
    delete score; score = 0;

    return lambdaHat;
}


   size_t n_dimension = 0; // dimension of matrices
   const long int *n_score = 0; // score_ [0...dimension_ - 1]
   const double *n_prob = 0; // prob_ [0...dimension_ - 1]
   long int n_morgue = 0; // score_ [0] - 1
   long int n_entry = 0; // n_entry = 0 : weak descending ladder epoch ; n_entry = -1 : strict descending ladder epoch 

   void n_setParameters (
   size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
   const long int *score_, // scores 
   const double *prob_, // probabilities
   long int entry_ = 0) // entry_ = 0 : weak descending ladder epoch ; entry_ = -1 : strict descending ladder epoch
   {
      n_dimension = dimension_;
      n_score = score_;
      n_prob = prob_;
      n_morgue = score_ [0] - 1;
      n_entry = entry_;
   }

   double n_totalProbAssoc (double x_)
   {
      double sum = 0.0;
      for (size_t i = 0; i < n_dimension; i++) {
         sum += n_prob [i] * exp (x_ * static_cast <double> (n_score [i]));
      }
      return sum;
   }

   double n_meanPowerAssoc (double x_, long int power_ = 1L)
   {
      double sum = 0.0;
      for (size_t i = 0; i < n_dimension; i++) {
          sum += Integer::integerPower (static_cast <double> (n_score [i]), power_) * 
            n_prob [i] * exp (x_ * static_cast <double> (n_score [i]));
      }
      return sum;
   }

   double n_meanAssoc (double x_)
   {
      return n_meanPowerAssoc (x_);
   }

   void n_bracket (double *p_, double *q_)
   {
      const double FACTOR = 0.5;
      *p_ = -log (n_prob [n_dimension - 1]) / static_cast <double> (n_score [n_dimension - 1]);
      while (1.0 <= n_totalProbAssoc (*p_)) {
         *p_ *= FACTOR;
      }
      *q_ = *p_ / FACTOR;
   }


double LocalMaxStatUtil::mu (
size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *score_, // scores in increasing order
const double *prob_) // corresponding probabilities
{
   double mu = 0.0;
   for (size_t i = 0; i < dimension_; i++) {
      mu += static_cast <double> (score_ [i]) * prob_ [i];
   }
   return mu;
}

double LocalMaxStatUtil::lambda (
size_t dimension_, // #(distinct values)          
const long int *score_, // values 
const double *prob_) // probability of corresponding value  
{
   n_setParameters (dimension_, score_, prob_);

   double p = 0.0;
   double q = 0.0;

   n_bracket (&p, &q);

   return Root::bisection (1.0, n_totalProbAssoc, p, q, REL_TOL * fabs (p - q));
}

double LocalMaxStatUtil::muPowerAssoc (
size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *score_, // scores in increasing order
const double *prob_, // corresponding probabilities
double lambda_, // lambda
long int power_) // power
{
   n_setParameters (dimension_, score_, prob_);

   if (lambda_ == 0.0) lambda_ = lambda (dimension_, score_, prob_);

   return n_meanPowerAssoc (lambda_, power_);
}

double LocalMaxStatUtil::muAssoc (
size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *score_, // scores in increasing order
const double *prob_, // corresponding probabilities
double lambda_) // lambda
{
   return muPowerAssoc (dimension_, score_, prob_, lambda_);
}

double LocalMaxStatUtil::thetaMin (
size_t dimension_, // #(distinct values)          
const long int *score_, // values 
const double *prob_, // probability of corresponding value 
double lambda_) // lambda
// assumes logarithmic regime
{
   n_setParameters (dimension_, score_, prob_);

   if (lambda_ == 0.0) lambda_ = lambda (dimension_, score_, prob_);

   double p = 0.0;
   double q = 0.0;
   n_bracket (&p, &q);
   return Root::bisection (0.0, n_meanAssoc, 0.0, lambda_, REL_TOL * fabs (p - q));
}

double LocalMaxStatUtil::rMin (
size_t dimension_, // #(distinct values)          
const long int *score_, // values 
const double *prob_, // probability of corresponding value 
double lambda_, // lambda
double thetaMin_) // argument of rate
// assumes logarithmic regime
{
   n_setParameters (dimension_, score_, prob_);

   if (thetaMin_ == 0.0) thetaMin_ = thetaMin (dimension_, score_, prob_, lambda_);

   return n_totalProbAssoc (thetaMin_);
}

double LocalMaxStatUtil::r ( // r (theta)
size_t dimension_, // #(distinct values)          
const long int *score_, // scores in increasing order
const double *prob_, // probability of corresponding value 
double theta_) // argument of rate
// assumes logarithmic regime
{
   double sum = 0.0;
   for (size_t i = 0; i < dimension_; i++) {
      sum += prob_ [i] * exp (theta_ * static_cast <double> (score_ [i]));
   }
   return sum;
}

long int LocalMaxStatUtil::delta ( // theta [minus delta] for ungapped sequence comparison
size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *score_) // scores 
{
   size_t i = 0;

   long int delta = 0;
   for (i = 0; i < dimension_; i++) {
      delta = Integer::euclidAlgorithm <long int> (delta, score_ [i]);
   }
   return delta;
}

double LocalMaxStatUtil::thetaMinusDelta ( // theta [minus delta] for ungapped sequence comparison
double lambda_, // lambda, the exponential rate for the local maximum         
size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *score_) // scores 
{
   double del = static_cast <double> (delta (dimension_, score_));
   return (1.0 - exp (-lambda_ * del)) / del;
} 


   long int n_step (long int oldValue_, size_t state_)
   {
      assert (state_ < n_dimension);
      return n_morgue < oldValue_ ? oldValue_ + n_score [state_] : oldValue_;
   }

   long int n_bury (long int oldValue_, size_t state_)
   {
      assert (state_ < n_dimension);
      return n_entry < oldValue_ ? oldValue_ : n_morgue;
   }


void LocalMaxStatUtil::descendingLadderEpochRepeat (
size_t dimension_, // #(distinct values)          
const long int *score_, // values 
const double *prob_, // probability of corresponding value 
double *eSumAlpha_, // expectation (sum [alpha])
double *eOneMinusExpSumAlpha_, // expectation [1.0 - exp (sum [alpha])]
bool isStrict_, // ? is this a strict descending ladder epoch
double lambda_, // lambda for repeats : default is lambda0_ below
size_t endW_, // maximum w plus 1
double *pAlphaW_, // probability {alpha = w} : pAlphaW_ [0, wEnd)
double *eOneMinusExpSumAlphaW_, // expectation [1.0 - exp (sum [alpha]); alpha = w] : eOneMinusExpSumAlphaW_ [0, wEnd)
double lambda0_, // lambda for flattened distribution (avoid recomputation)
double mu0_, // mean of flattened distribution (avoid recomputation)
double muAssoc0_, // mean of associated flattened distribution (avoid recomputation)
double thetaMin0_, // thetaMin of flattened distribution (avoid recomputation)
double rMin0_, // rMin of flattened distribution (avoid recomputation)
double time_, // get time for the dynamic programming computation
bool *terminated_) // ? Was the dynamic programming computation terminated prematurely ?
// assumes logarithmic regime
{
    // Start dynamic programming probability calculation using notation in
    //
    // Mott R. and Tribe R. (1999)
    // J. Computational Biology 6(1):91-112
    //
    // Karlin S. and Taylor H.M.(1981)
    // A Second Course in Stochastic Processes, p. 480
    //
    // Note there is an error in Eq (6.19) there, which is corrected in Eq (6.20)
    //
    // This program uses departure into (-Inf, 0] not (-Inf, 0)

    // avoid recomputation
    double mu0 = 0.0 == mu0_ ? mu (dimension_, score_, prob_) : mu0_;
    assert (mu0 < 0.0);
    double lambda0 = 0.0 == lambda0_ ? lambda (dimension_, score_, prob_) : lambda0_;
    assert (0.0 < lambda0);
    if (lambda_ == 0.0) lambda_ = lambda0;
    assert (0.0 < lambda_);
    double muAssoc0 = 0.0 == muAssoc0_ ? muAssoc (dimension_, score_, prob_, lambda0) : muAssoc0_;
    assert (0.0 < muAssoc0);
    double thetaMin0 = 0.0 == thetaMin0_ ? thetaMin (dimension_, score_, prob_, lambda0) : thetaMin0_;
    assert (0.0 < thetaMin0);
    double rMin0 = 0.0 == rMin0_ ? rMin (dimension_, score_, prob_, lambda0, thetaMin0) : rMin0_;
    assert (0.0 < rMin0 && rMin0 < 1.0);

    const long int ITER_MIN = static_cast <long int> ((log (REL_TOL * (1.0 - rMin0)) / log (rMin0)));
    assert (0 < ITER_MIN);
    const long int ITER = static_cast <long int> (endW_) < ITER_MIN ? ITER_MIN : static_cast <long int> (endW_);
    assert (0 < ITER);
    const long int Y_MAX = static_cast <long int> (-log (REL_TOL) / lambda0);

    long int entry = isStrict_ ? -1 : 0;
    n_setParameters (dimension_, score_, prob_, entry);


    double time0 = 0.0;
    double time1 = 0.0;
	if (time_ > 0.0) Sls::sls_basic::get_current_time (time0);

    DynProgProbLim dynProgProb (n_step, dimension_, prob_, score_ [0] - 1, Y_MAX);

    if (pAlphaW_) pAlphaW_ [0] = 0.0;
    if (eOneMinusExpSumAlphaW_) eOneMinusExpSumAlphaW_ [0] = 0.0;

    dynProgProb.update (); // iterate random walk

    long int value = 0;

    if (eSumAlpha_) *eSumAlpha_ = 0.0;
    if (eOneMinusExpSumAlpha_) *eOneMinusExpSumAlpha_ = 0.0;

    for (size_t w = 1; w < static_cast <size_t> (ITER); w++) {

        if (w < endW_) { // sum pAlphaW_ [w] and eOneMinusExpSumAlphaW_ [w]

             if (pAlphaW_) pAlphaW_ [w] = 0.0;
             if (eOneMinusExpSumAlphaW_) eOneMinusExpSumAlphaW_ [w] = 0.0;

             for (value = score_ [0]; value <= entry; value++) {
                if (pAlphaW_) pAlphaW_ [w] += dynProgProb.getProb (value);
                if (eOneMinusExpSumAlphaW_) eOneMinusExpSumAlphaW_ [w] += 
                                               dynProgProb.getProb (value) * 
                                               (1.0 - exp (lambda_ * static_cast <double> (value)));
             }
        }

        for (value = score_ [0]; value <= entry; value++) {
         if (eSumAlpha_) *eSumAlpha_ += dynProgProb.getProb (value) * static_cast <double> (value);
         if (eOneMinusExpSumAlpha_) *eOneMinusExpSumAlpha_ += dynProgProb.getProb (value) * 
                                        (1.0 - exp (lambda_ * static_cast <double> (value)));
        }

        dynProgProb.setValueFct (n_bury); 
        dynProgProb.update (); // put probability into the morgue

        dynProgProb.setValueFct (n_step); 
        dynProgProb.update (); // iterate random walk

        if (time_ > 0.0)
        {
			Sls::sls_basic::get_current_time (time1);
            if (time1 - time0 > time_) 
            {
                *terminated_ = true;
                return;
            }
        }

    }

    for (value = score_ [0]; value <= entry; value++) {
      if (eSumAlpha_) *eSumAlpha_ += dynProgProb.getProb (value) * static_cast <double> (value);
      if (eOneMinusExpSumAlpha_) *eOneMinusExpSumAlpha_ += dynProgProb.getProb (value) * 
                                     (1.0 - exp (lambda_ * static_cast <double> (value)));
    }

    // check that not too much probability has been omitted
    double prob = 0.0;
    for (value = entry + 1; value < dynProgProb.getValueUpper (); value++) {
      prob += dynProgProb.getProb (value);
    }
    prob += dynProgProb.getProbLost ();

    const double FUDGE = 100.0;
    assert (prob <= FUDGE * static_cast <double> (dimension_) * REL_TOL);
}

void LocalMaxStatUtil::descendingLadderEpoch (
size_t dimension_, // #(distinct values)          
const long int *score_, // values 
const double *prob_, // probability of corresponding value 
double *eSumAlpha_, // expectation (sum [alpha])
double *eOneMinusExpSumAlpha_, // expectation [1.0 - exp (sum [alpha])]
bool isStrict_, // ? is this a strict descending ladder epoch
double lambda0_, // lambda for flattened distribution (avoid recomputation)
double mu0_, // mean of flattened distribution (avoid recomputation)
double muAssoc0_, // mean of associated flattened distribution (avoid recomputation)
double thetaMin0_, // thetaMin of flattened distribution (avoid recomputation)
double rMin0_, // rMin of flattened distribution (avoid recomputation)
double time_, // get time for the dynamic programming computation
bool *terminated_) // ? Was the dynamic programming computation terminated prematurely ?
{
   descendingLadderEpochRepeat (dimension_, score_, prob_, 
      eSumAlpha_, eOneMinusExpSumAlpha_, isStrict_, 0.0, 0, 0, 0, 
      lambda0_, mu0_, muAssoc0_, thetaMin0_, rMin0_, time_, terminated_); 
}

bool LocalMaxStatUtil::isProbDist (
size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
const double *prob_) // corresponding probabilities
{
   double sum = 0.0;
   for (size_t i = 0; i < dimension_; i++) {
      if (prob_ [i] < 0.0 || 1.0 < prob_ [i]) return false;
      sum += prob_ [i];
   }
   return Approx::relApprox (sum, 1.0, REL_TOL);
}

bool LocalMaxStatUtil::isScoreIncreasing (
size_t dimension_, // #(distinct values)          
const long int *score_) // scores in increasing order
{
   for (size_t i = 1; i < dimension_; i++) {
      if (score_ [i] <= score_ [i - 1]) return false;
   }
   return true;
}

bool LocalMaxStatUtil::isLogarithmic (
size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
const long int *score_, // scores in increasing order
const double *prob_) // corresponding probabilities
{
   assert (score_);
   assert (prob_);
   if (! isScoreIncreasing (dimension_, score_)) return false;
   if (! isProbDist (dimension_, prob_)) return false;
   if (0.0 <= mu (dimension_, score_, prob_)) return false;
   if (score_ [dimension_ - 1] <= 0.0) return false;
   return true;
}

