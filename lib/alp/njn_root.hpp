#ifndef INCLUDED_NJN_ROOT
#define INCLUDED_NJN_ROOT

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

File name: njn_root.hpp

Author: John Spouge

Contents: 

******************************************************************************/

#include <math.h>

#include "njn_approx.hpp"
#include "njn_function.hpp"
#include "njn_ioutil.hpp"


namespace Njn {
	namespace Root {


      const double FAILED = HUGE_VAL;

      // All routines find roots. 
      // They return FAILED if the root is not located within *itmax_ iterations.
      // If not a default 0 pointer, *itmax = iterations left.

      template <typename T>
      double newtonRaphson ( // finds root f_(x_) = y_ in [p_, q_] with derivative y_'
      double y_, // f_ (x_) = y_
      double (*f_) (double, const T &), // function : param_ contains parameters 
      double (*df_) (double, const T &), // derivative of function
      const T &param_, // parameters for function
      double p_, // end-point
      double x_, // initial guess : can be arbitrary
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance : set to 0.0 to ignore
      long int *itmax_); // maximum # of permitted iterations : default = 0
      // bisection locates x_ : f_ (x_) = y_ to within absolute error +-tol_, 
      //    then uses derivative information through a Newton-Raphson alternative
      //
      // asserts f_ (p_) <= y_ <= f_ (q_) or f_ (q_) <= y_ <= f_ (p_)

      template <typename T>
      double bisection ( // finds root f_ (x_) = y_ in [p_, q_]
      double y_, // f_ (x_) = y_
      double (*f_) (double, const T &), // function : param_ contains parameters 
      const T &param_, // parameters for function
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance : set to 0.0 to ignore
      long int *itmax_); // maximum # of permitted iterations : default = 0
      // bisection locates x_ : f_ (x_) = y to within absolute error +-tol_, 
      //
      // asserts f_ (p_) <= y_ <= f_ (q_) or f_ (q_) <= y_ <= f_ (p_)

      template <typename T>
      double hunt ( // finds root f_(x_) = y_ in (p_, q_) by looking until mesh = tol_
      double y_, // f_ (x_) = y_
      double (*f_) (double, const T &), // function : param_ contains parameters 
      const T &param_, // parameters for function
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance : set to 0.0 to ignore
      long int *itmax_); // maximum # of permitted iterations : default = 0
      // hunt locates x_ : f_ (x_) = y to within absolute error +-tol_, 

      template <typename T>
      double huntExtreme ( // finds root f_(x_) = y_ in (p_, q_) by looking until mesh = tol_
      double y_, // f_ (x_) = y_
      double (*f_) (double, const T &), // function : param_ contains parameters 
      const T &param_, // parameters for function
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance : set to 0.0 to ignore
      long int *itmax_, // maximum # of permitted iterations : default = 0
      bool isLargest_); // ? find the largest root ?
      // huntExtreme locates x_ : f_ (x_) = y to within absolute error +-tol_,
      //    finding either the largest or smallest root in the interval (p_, q_)
      // huntExtreme looks for a change in sign, so y_ should not be an extremum.

      //
      // Specializations of the template follow.
      //

      inline double newtonRaphson ( // finds root f_(x_) = y_ in [p_, q_] with derivative y_'
      double y_, // f_ (x_) = y_
      double (*f_) (double), // function
      double (*df_) (double), // derivative of function
      double p_, // end-point
      double x_, // initial guess : can be arbitrary
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_ = 0.0, // relative tolerance : set to 0.0 to ignore
      long int *itmax_ = 0); // maximum # of permitted iterations : default = 0
      // bisection routine locates x_ : f_ (x_) = y_ to within absolute error +-tol_, 
      //    then uses derivative information through a Newton-Raphson alternative
      //
      // asserts f_ (p_) <= y_ <= f_ (q_) or f_ (q_) <= y_ <= f_ (p_)

      inline double bisection ( // finds root f_(x_) = y_ in [p_, q_]
      double y_, // f_ (x_) = y_
      double (*f_) (double), // function 
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_ = 0.0, // relative tolerance : set to 0.0 to ignore
      long int *itmax_ = 0); // maximum # of permitted iterations : default = 0
      // bisection routine locates x_ : f_ (x_) = y to within absolute error +-tol_, 
      //
      // asserts f_ (p_) <= y_ <= f_ (q_) or f_ (q_) <= y_ <= f_ (p_)

      inline double hunt ( // finds root f_ (x_) = y_ in (p_, q_) by looking until mesh = tol_
      double y_, // f_ (x_) = y_
      double (*f_) (double), // function : param_ contains parameters 
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_ = 0.0, // relative tolerance : set to 0.0 to ignore
      long int *itmax_ = 0); // maximum # of permitted iterations : default = 0

      inline double huntExtreme ( // finds root f_(x_) = y_ in (p_, q_) by looking until mesh = tol_
      double y_, // f_ (x_) = y_
      double (*f_) (double), // function : param_ contains parameters 
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_ = 0.0, // relative tolerance : set to 0.0 to ignore
      long int *itmax_ = 0, // maximum # of permitted iterations : default = 0
      bool isLargest_ = false); // ? find the largest root ?
      // huntExtreme locates x_ : f_ (x_) = y to within absolute error +-tol_,
      //    finding either the largest or smallest root in the interval (p_, q_)
      // huntExtreme looks for a change in sign, so y_ should not be an extremum.

	}
}

//
// There are no more declarations beyond this point.
//

namespace Njn {
	namespace Root {

namespace {

         typedef double DoubleFct (double);

         DoubleFct *s_f = 0;
         DoubleFct *s_df = 0;
         const double ZERO = 0.0;

         double f (double x_, const double &/*sls deleted y_*/ ) {return (*s_f) (x_);}
         double df (double x_, const double &/*sls deleted y_*/ ) {return (*s_df) (x_);}

	}

      double newtonRaphson ( // finds root f_(x_) = y_ in [p_, q_] with derivative y_'
      double y_, // f_ (x_) = y_
      double (*f_) (double x_), // function
      double (*df_) (double x_), // derivative of function
      double p_, // end-point
      double x_, // initial guess : can be arbitrary
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance
      long int *itmax_) // maximum # of permitted iterations
      {
         s_f = f_;
         s_df = df_;
         return newtonRaphson (y_, f, df, ZERO, p_, x_, q_, tol_, rtol_, itmax_);
      }

      double bisection ( // finds root f_(x_) = y_ in [p_, q_]
      double y_, // f_ (x_) = y_
      double (*f_) (double x_), // function 
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance
      long int *itmax_) // maximum # of permitted iterations
      {
         s_f = f_;
         return bisection (y_, f, ZERO, p_, q_, tol_, rtol_, itmax_);
      }

      double hunt ( // finds root f_(x_) = y_ in (p_, q_)
      double y_, // f_ (x_) = y_
      double (*f_) (double x_), // function 
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance
      long int *itmax_) // maximum # of permitted iterations
      {
         s_f = f_;
         return hunt (y_, f, ZERO, p_, q_, tol_, rtol_, itmax_);
      }

      double huntExtreme ( // finds root f_(x_) = y_ in (p_, q_)
      double y_, // f_ (x_) = y_
      double (*f_) (double x_), // function 
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance
      long int *itmax_, // maximum # of permitted iterations
      bool isLargest_) // ? find the largest root ?
      {
         s_f = f_;
         return huntExtreme (y_, f, ZERO, p_, q_, tol_, rtol_, itmax_, isLargest_);
      }

      template <typename T>
      double newtonRaphson ( // finds root f_(x_) = y_ in [p_, q_] with derivative y_'
      double y_, // f_ (x_) = y_
      double (*f_) (double, const T &), // function : param_ contains parameters 
      double (*df_) (double, const T &), // derivative of function
      const T &param_, // parameters for function
      double p_, // end-point
      double x_, // initial guess : can be arbitrary
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance : set to 0.0 to ignore
      long int *itmax_) // maximum # of permitted iterations : default = 0
      {
         // #define F(x_)  ((*f_) (x_, param_) - y_)
         // #define DF(x_) ((*df_) (x_, param_))

         assert (p_ != HUGE_VAL && p_ != -HUGE_VAL);
         assert (q_ != HUGE_VAL && q_ != -HUGE_VAL);

         // checks for improper bracketing and end-point root.
         double fp = (*f_) (p_, param_) - y_;
         double fq = (*f_) (q_, param_) - y_;
         if (fp * fq > 0.0) 
            IoUtil::abort ("Root::newtonRaphson : root not bracketed");

         if (fp == 0.0) return p_;
         if (fq == 0.0) return q_;

         if (p_ == q_) IoUtil::abort ("Root::newtonRaphson : p_ == q_");

         double x = x_;
         // swaps end-points if necessary to make p_ < q_
         //if (q_ < p_) std::swap <double> (p_, q_);
		 if (q_ < p_) std::swap (p_, q_);/*sls deleted <double>*/

         // makes an initial guess within [p_, q_]
         if (x_ < p_ || q_ < x_) x = 0.5 * (p_ + q_);

         // swaps end-points if necessary to make F (p_) < 0.0 < F (q_)
         //if (fp > 0.0) std::swap <double> (p_, q_);
		 if (fp > 0.0) std::swap (p_, q_);/*sls deleted <double>*/
   
         // Set up the bisection & Newton-Raphson iteration.

         double dx; // present interval length
	      double dxold; // old interval length
	      double fx; // f_(x_)-y_
	      double dfx; // Df(x_)

	      dxold = dx = p_ - q_;

         long int iter = 100; // default iterations
         long int *itmax = itmax_ == 0 ? &iter: itmax_;
         
	      for ( ; 0 < *itmax; --*itmax) {

		      fx = (*f_) (x, param_) - y_; 
		      if (fx == 0.0) {           // Check for termination.
               return x;
		      } else if (fx < 0.0) {
			      p_ = x;
		      } else {
			      q_ = x;
		      }
		      dfx = (*df_) (x, param_) - y_;    
      
            // Is the root out of bounds, so bisection is faster than Newton-Raphson?		      
            if ((dfx * (x-p_) - fx) * (dfx * (x - q_) - fx) >= 0.0 || 
		           2.0 * fabs (fx) > fabs (dfx * dx)) { 
         
               // bisect
			      dx = dxold; 
			      dxold = 0.5 * (p_ - q_);
			      x = 0.5 * (p_ + q_);

			      if (fabs (dxold) <= tol_) return x;

		      } else {  
         
               // Newton-Raphson
			      dx = dxold; 
			      dxold = fx / dfx;
			      x -= dxold;

			      if (fabs (dxold) < tol_ || fabs (dxold) < rtol_ * fabs (x)) {
				      if (((*f_) ((x - Function::signum (dxold) * tol_), param_) - y_) * fx < 0.0) return x;
			      }
		      }
	      }
   
         return FAILED; // failure
      }

      template <typename T>
      double bisection ( // finds root f_ (x_) = y_ in [p_, q_]
      double y_, // f_ (x_) = y_
      double (*f_) (double, const T &), // function : param_ contains parameters 
      const T &param_, // parameters for function
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance : set to 0.0 to ignore
      long int *itmax_) // maximum # of permitted iterations : default = 0
      {
         assert (p_ != HUGE_VAL && p_ != -HUGE_VAL);
         assert (q_ != HUGE_VAL && q_ != -HUGE_VAL);

         // checks for improper bracketing and end-point root.
         double fp = (*f_) (p_, param_) - y_;
         double fq = (*f_) (q_, param_) - y_;
         if (fp * fq > 0.0) 
            IoUtil::abort ("Root::bisection : root not bracketed");

         if (fp == 0.0) return p_;
         if (fq == 0.0) return q_;

         if (p_ == q_) IoUtil::abort ("Root::bisection : p_ == q_");

         // swaps end-points if necessary to make F (p_) < 0.0 < F (q_)
         //if (fp > 0.0) std::swap <double> (p_, q_);
		 if (fp > 0.0) std::swap (p_, q_);/*sls deleted <double>*/

         double x = 0.0;
         double fx = 0.0;

         long int iter = 100; // default iterations
         long int *itmax = itmax_ == 0 ? &iter: itmax_;
         
         x = 0.5 * (p_ + q_);
	      for ( ; 0 < *itmax; --*itmax) {

		      fx = (*f_) (x, param_) - y_;    
            if (fx < 0.0) {
			      p_ = x;
            } else {
			      q_ = x;
            }
            x = 0.5 * (p_ + q_);

            if (Approx::absRelApprox <double> (p_, x, tol_, rtol_)) return x;
	      }

         return FAILED; // failure
      }

      template <typename T>
      double hunt ( // finds root f_(x_) = y_ in (p_, q_) by looking until mesh = tol_
      double y_, // f_ (x_) = y_
      double (*f_) (double, const T &), // function : param_ contains parameters 
      const T &param_, // parameters for function
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance : set to 0.0 to ignore
      long int *itmax_) // maximum # of permitted iterations : default = 0
      {
         assert (p_ != HUGE_VAL && p_ != -HUGE_VAL);
         assert (q_ != HUGE_VAL && q_ != -HUGE_VAL);

         if (p_ == q_) IoUtil::abort ("Root::hunt : p_ == q_");

         double x0 = 0.5 * (p_ + q_);
         double fx0 = (*f_) (x0, param_) - y_;
         if (fx0 == 0.0) return x0;

         // swaps end-points if necessary to make p_ < q_
         //if (q_ < p_) std::swap <double> (p_, q_);
		 if (q_ < p_) std::swap (p_, q_);/*sls deleted <double>*/


         size_t pts = 2;
         double del = 0.5 * (q_ - p_);
         double x = 0.0;
         double fx = 0.0;

         long int iter = 1000; // default iterations
         long int *itmax = itmax_ == 0 ? &iter: itmax_;
         
	      while (tol_ <= del) {

            x = p_ + 0.5 * del;
            for (size_t i = 0; i < pts && 0 < *itmax; i++, --*itmax) {

               fx = (*f_) (x, param_) - y_;
               if (fx * fx0 < 0.0) return bisection <T> (y_, f_, param_, x, x0, tol_, rtol_, itmax);
               x += del;
            }
            if (iter == 0) return FAILED;

            pts *= 2;
            del *= 0.5;
	      }

         return FAILED; // failure
      }

      template <typename T>
      double huntExtreme ( // finds root f_(x_) = y_ in (p_, q_) by looking until mesh = tol_
      double y_, // f_ (x_) = y_
      double (*f_) (double, const T &), // function : param_ contains parameters 
      const T &param_, // parameters for function
      double p_, // end-point
      double q_, // end-point
      double tol_, // absolute tolerance
      double rtol_, // relative tolerance : set to 0.0 to ignore
      long int *itmax_, // maximum # of permitted iterations : default = 0
      bool isLargest_) // ? find the largest root ?
      {
         long int iter = 1000; // default iterations
         long int *itmax = itmax_ == 0 ? &iter: itmax_;
         
         // swaps end-points if necessary to make p_ < q_
        // if (q_ < p_) std::swap <double> (p_, q_);
		 if (q_ < p_) std::swap (p_, q_);/*sls deleted <double>*/

         // check there is a root
         double x = hunt <T> (y_, f_, param_, p_, q_, tol_, rtol_, itmax);
         double x0 = x;

         // find the extreme root
         if (isLargest_) {
            while (0 < *itmax && x != FAILED) {
               x0 = x;
               x = hunt <T> (y_, f_, param_, x, q_, tol_, rtol_, itmax);
            }
         } else {
            while (0 < *itmax && x != FAILED) {
               x0 = x;
               x = hunt <T> (y_, f_, param_, p_, x, tol_, rtol_, itmax);
            }
         }

         return x0; 
      }


	}
}

#endif //! INCLUDED

