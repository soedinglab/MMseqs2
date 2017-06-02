#ifndef INCLUDED_NJN_LOCALMAXSTATMATRIX
#define INCLUDED_NJN_LOCALMAXSTATMATRIX

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

File name: njn_localmaxstatmatrix.hpp

Author: John Spouge

Contents: Random walk parameters

******************************************************************************/

 
#include "njn_localmaxstat.hpp"


namespace Njn {



     class LocalMaxStatMatrix : public LocalMaxStat { 

    // calculates the statistical parameters for the local maximum in a random walk
    //
    // The scores are uniqued and 
    //    with the correspondence to probabilities maintained, placed in ascending order.
    //
    // The default p2_ = 0 is equivalent to the symmetric probabilities p2_ = p_ (on the same alphabet).
    //
    // The default dimMatrix2_ = 0 is equivalent to the same alphabet dimMatrix2_ = dimMatrix_.

        public:

        inline LocalMaxStatMatrix ( 
        size_t dimMatrix_ = 0, // #(distinct values)          
        const long int *const *scoreMatrix_ = 0, // score matrix [0...dimMatrix_)[0...dimMatrix2_)
        const double *p_ = 0, // probability of "letters" p_ [0...dimMatrix_)
        const double *p2_ = 0, // probability of "letters" p2_ [0...dimMatrix2_), the second (j) set of letter-probabilities
        size_t dimMatrix2_ = 0, // #(distinct values) in the second alphabet         
		double time_=0)
        :  LocalMaxStat (), d_dimMatrix (0), d_scoreMatrix_p (0), d_p_p (0), d_p2_p (0), d_dimMatrix2 (0)
        {
			setTime(time_);
            copy (dimMatrix_, scoreMatrix_, p_, p2_, dimMatrix2_);
        }

        inline ~LocalMaxStatMatrix () {free2 ();}

        inline LocalMaxStatMatrix &operator= (const LocalMaxStatMatrix &localMaxStat_) // random walk parameters
        {
            if (this != &localMaxStat_) copy (localMaxStat_);
            return *this;
        }

        void copy (
        size_t dimMatrix_, // #(distinct values) of scores & probabilities (which are paired)         
        const long int *const *scoreMatrix_, // score matrix [0...dimMatrix_)[0...dimMatrix_)
        const double *p_, // probability of "letters" p_ [0...dimMatrix_)
        const double *p2_ = 0, // probability of "letters" p2_ [0...dimMatrix2_), the second (j) set of letter-probabilities
        size_t dimMatrix2_ = 0); // #(distinct letters) in the second alphabet         

        void copy (
        LocalMaxStat localMaxStat_, // base object 
        size_t dimMatrix_, // #(distinct values) of scores & probabilities (which are paired)         
        const long int *const *scoreMatrix_, // score matrix [0...dimMatrix_)[0...dimMatrix_)
        const double *p_, // probability of "letters" p_ [0...dimMatrix_)
        const double *p2_ = 0, // probability of "letters" p2_ [0...dimMatrix2_), the second (j) set of letter-probabilities
        size_t dimMatrix2_ = 0); // #(distinct letters) in the second alphabet         

        inline void copy (const LocalMaxStatMatrix &localMaxStatMatrix_)
        {
            copy (localMaxStatMatrix_, localMaxStatMatrix_.getDimMatrix (), localMaxStatMatrix_.getScoreMatrix (), localMaxStatMatrix_.getP (), localMaxStatMatrix_.getP2 (), localMaxStatMatrix_.getDimMatrix2 ());
        }

        using LocalMaxStat::operator bool; // ? is the object ready for computation ?
        using LocalMaxStat::out; // output
        using LocalMaxStat::getR; // r (theta_) : dominant eigenvalue for theta_
        using LocalMaxStat::getA; // lim expected [length] / y for achieving y
        using LocalMaxStat::getAlpha; // lim var [length] / y for achieving y
        using LocalMaxStat::getDimension; // #(distinct values) of scores & probabilities (which are paired)         
        using LocalMaxStat::getScore; // scores in increasing order
        using LocalMaxStat::getProb; // probabilities
        using LocalMaxStat::getLambda; // lambda for associated random walk
        using LocalMaxStat::getK; // k for random walk : exponential prefactor
        using LocalMaxStat::getC; // c for random walk : exponential prefactor (global alignment)
        using LocalMaxStat::getThetaMin; // theta for minimum expectation (exp (theta * score))
        using LocalMaxStat::getRMin; // minimum expectation (exp (theta * score))
        using LocalMaxStat::getDelta; // span
        using LocalMaxStat::getThetaMinusDelta; // renewal span parameter
        using LocalMaxStat::getMu; // step mean for random walk
        using LocalMaxStat::getSigma; // step standard deviation for random walk
        using LocalMaxStat::getMuAssoc; // step mean for associated random walk (relative entropy)
        using LocalMaxStat::getSigmaAssoc; // step standard deviation for associated random walk
        using LocalMaxStat::getMeanWDLE; // expected renewal length for weak ladder epochs

        inline size_t getDimMatrix () const {return d_dimMatrix;} // #(distinct values) of scores & probabilities (which are paired)
        inline const long int *const *getScoreMatrix () const {return d_scoreMatrix_p;} // score matrix [0...dimMatrix_)[0...dimMatrix_)
        inline const double *getP () const {return d_p_p;} // probability of "letters" d_p_p [0...dimMatrix_)
        inline const double *getP2 () const {return d_p2_p;}  // probability of "letters" p2_ [0...dimMatrix_), the second (j) set of letter-probabilities
        inline size_t getDimMatrix2 () const {return d_dimMatrix2;} // #(distinct letters) in the second alphabet

        private:

        size_t d_dimMatrix; // #(distinct values) of scores & probabilities (which are paired)         
        long int **d_scoreMatrix_p; // score matrix [0...dimMatrix_)[0...dimMatrix_)
        double *d_p_p; // probability of "letters" d_p_p [0...dimMatrix_)
        double *d_p2_p; // probability of "letters" p2_ [0...dimMatrix_), the second (j) set of letter-probabilities
        size_t d_dimMatrix2; // #(distinct letters) in the second alphabet

        void init (size_t dimMatrix_, size_t dimMatrix2_ = 0);
        void free2 ();
    };

	}


#endif //!INCLUDED_NJS_LOCALMAXSTAT 

