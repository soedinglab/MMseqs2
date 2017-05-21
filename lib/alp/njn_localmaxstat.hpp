#ifndef INCLUDED_NJN_LOCALMAXSTAT
#define INCLUDED_NJN_LOCALMAXSTAT

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

File name: njn_localmaxstat.hpp

Author: John Spouge

Contents: Random walk parameters

******************************************************************************/

#include <math.h>
#include <vector>
#include <assert.h>
#include <ostream>



namespace Njn {


    class LocalMaxStat { 

    // calculates the statistical parameters for the local maximum in a random walk
    //
    // The scores are uniqued and 
    //    with the correspondence to probabilities maintained, placed in ascending order.

      public:

      // The following subroutines control the time for the dynamic programming computation.
      // The default time_ = 0.0 permits the computation to run forever.
      static void setTime (double time_ = 0.0) {assert (time_ >= 0.0); s_time = time_;} // set time for the dynamic programming computation
      static double getTime () {return s_time;} // get time for the dynamic programming computation
      // For an object o, 
      //    if the computation is terminated before it finishes, 
      //    o.getTerminated () == true. 

      inline LocalMaxStat ( 
      size_t dimension_ = 0, // #(distinct values)          
      const long int *score_ = 0, // scores in increasing order
      const double *prob_ = 0) // probability of corresponding value  
      :  d_dimension (0), d_score_p (0), d_prob_p (0),
         d_lambda (0.0), d_k (0.0), d_c (0.0), d_thetaMin (0.0), d_rMin (0.0), 
         d_delta (0), d_thetaMinusDelta (0.0),
         d_mu (0.0), d_sigma (0.0), d_muAssoc (0.0), d_sigmaAssoc (0.0),
         d_meanWDLE (0.0), d_terminated (false)
      {
         copy (dimension_, score_, prob_);
      }

      inline LocalMaxStat (const LocalMaxStat &localMaxStat_) // random walk parameters
      :  d_dimension (0), d_score_p (0), d_prob_p (0),
         d_lambda (0.0), d_k (0.0), d_c (0.0), d_thetaMin (0.0), d_rMin (0.0), 
         d_delta (0), d_thetaMinusDelta (0.0),
         d_mu (0.0), d_sigma (0.0), d_muAssoc (0.0), d_sigmaAssoc (0.0),
         d_meanWDLE (0.0), d_terminated (false)
      {
         copy (localMaxStat_);
      }

      inline ~LocalMaxStat () {free2 ();}

      inline operator bool () // ? is the object ready for computation ?
      const {
         return d_dimension != 0;
      }

      inline LocalMaxStat &operator= (const LocalMaxStat &localMaxStat_) // random walk parameters
      {
         if (this != &localMaxStat_) copy (localMaxStat_);
         return *this;
      }

      void copy (
      size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
      const long int *score_, // scores in increasing order
      const double *prob_); // probabilities

      inline void copy (const LocalMaxStat &localMaxStat_)
      {
         copy (localMaxStat_.getDimension (), localMaxStat_.getScore (), localMaxStat_.getProb (),
               localMaxStat_.getLambda (), localMaxStat_.getK (), localMaxStat_.getC (), 
               localMaxStat_.getThetaMin (), localMaxStat_.getRMin (),
               localMaxStat_.getDelta (), localMaxStat_.getThetaMinusDelta (),
               localMaxStat_.getMu (), localMaxStat_.getSigma (), localMaxStat_.getMuAssoc (), localMaxStat_.getSigmaAssoc (),
               localMaxStat_.getMeanWDLE (), localMaxStat_.getTerminated ());  
      }

      void copy (
      size_t dimension_, // #(distinct values) of scores & probabilities (which are paired)         
      const long int *score_, // scores in increasing order
      const double *prob_, // probabilities
      double lambda_, // lambda for associated random walk
      double k_, // k for random walk : exponential prefactor
      double c_, // c for random walk : exponential prefactor (global alignment)
      double thetaMin_, // theta for minimum expectation (exp (theta * score))
      double rMin_, // minimum expectation (exp (theta * score))
      long int delta_, // span 
      double thetaMinusDelta_, // renewal span parameter
      double mu_, // step mean for random walk
      double sigma_, // step standard deviation for random walk
      double muAssoc_, // step mean for associated random walk (relative entropy)
      double sigmaAssoc_, // step standard deviation for associated random walk
      double meanDLE_, // expected renewal length
      bool terminated_ = false); // ? Was the dynamic programming computation terminated prematurely ?

      inline std::ostream &out (std::ostream &ostr_) const {return ostr_;} // output

      double getR (double theta_) const; // r (theta_) : dominant eigenvalue for theta_

      double getA () const {return getMuAssoc () == 0 ? HUGE_VAL : 1.0 / getMuAssoc ();} // expected [length / y] for achieving y

      double getAlpha () const {return getSigmaAssoc () * getSigmaAssoc () * getA () * getA () * getA ();} // var [length] / y for achieving y

      inline size_t getDimension () const {return d_dimension;} // #(distinct values) of scores & probabilities (which are paired)         
      inline const long int *getScore () const {return d_score_p;} // scores in increasing order
      inline const double *getProb () const {return d_prob_p;} // probabilities
      inline double getLambda () const {return d_lambda;} // lambda for associated random walk
      inline double getK () const {return d_k;} // k for random walk : exponential prefactor
      inline double getC () const {return d_c;} // c for random walk : exponential prefactor (global alignment)
      inline double getThetaMin () const {return d_thetaMin;} // theta for minimum expectation (exp (theta * score))
      inline double getRMin () const {return d_rMin;} // minimum expectation (exp (theta * score))
      inline long int getDelta () const {return d_delta;} // span
      inline double getThetaMinusDelta () const {return d_thetaMinusDelta;} // renewal span parameter
      inline double getMu () const {return d_mu;} // step mean for random walk
      inline double getSigma () const {return d_sigma;} // step standard deviation for random walk
      inline double getMuAssoc () const {return d_muAssoc;} // step mean for associated random walk (relative entropy)
      inline double getSigmaAssoc () const {return d_sigmaAssoc;} // step standard deviation for associated random walk
      inline double getMeanWDLE () const {return d_meanWDLE;} // expected renewal length for weak ladder epochs
      inline bool getTerminated () const {return d_terminated;} // ? Was the dynamic programming computation terminated prematurely ?

      private:

      // random walk distribution
      size_t d_dimension; // #(distinct values) of scores & probabilities (which are paired)         
      long int *d_score_p; // scores in increasing order
      double *d_prob_p; // probabilities
      
      // Karlin-Altschul parameters
      double d_lambda; // lambda for associated random walk
      double d_k; // k for random walk : exponential prefactor
      double d_c; // c for random walk : exponential prefactor (global alignment)
      double d_thetaMin; // theta for minimum expectation (exp (theta * score))
      double d_rMin; // minimum expectation (exp (theta * score))
      long int d_delta; // span 
      double d_thetaMinusDelta; // renewal span parameter
      double d_mu; // step mean for random walk
      double d_sigma; // step standard deviation for random walk
      double d_muAssoc; // step mean for associated random walk (relative entropy)
      double d_sigmaAssoc; // step standard deviation for associated random walk
      double d_meanWDLE; // expected renewal length for weak ladder epochs
      bool d_terminated; // ? Was the dynamic programming computation terminated prematurely ?

      void init (size_t dimension_);
      void free2 ();

      void clear ();

      void dynProgCalc ();
      // k for random walk : exponential prefactor 
      // expected renewal length for weak ladder epochs

      static double s_time;

    };

}

#endif 

