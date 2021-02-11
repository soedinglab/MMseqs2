#ifndef INCLUDED_SLS_BASIC
#define INCLUDED_SLS_BASIC

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

File name: sls_basic.hpp

Author: Sergey Sheetlin, Martin Frith

Contents: Some basic functions and types

******************************************************************************/

#ifndef _MSC_VER //UNIX program
#else
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef _MSC_VER //UNIX program
#include <sys/time.h>
#else
#include <sys/timeb.h>
#endif

#include <iomanip>
#include <cmath>  // ?
#include <math.h>
#include <string>

namespace Sls { 

	const double pi=3.1415926535897932384626433832795;
	const double const_val=1/sqrt(2.0*pi);
	const long int quick_tests_trials_number=100;

	struct error//struct to handle exceptions
	{
		std::string st;
		error(std::string st_,long int error_code_){st=st_;error_code=error_code_;};
		long int error_code;
	};

	//structure for user-defined Gumbel parameters without errors
	struct AlignmentEvaluerParameters
	{
		double d_lambda;
		double d_k;
		double d_a1;
		double d_b1;
		double d_a2;
		double d_b2;
		double d_alpha1;
		double d_beta1;
		double d_alpha2;
		double d_beta2;
		double d_sigma;
		double d_tau;
	};

	//structure for user-defined Gumbel parameters with errors
	struct AlignmentEvaluerParametersWithErrors
	{
		double d_lambda;
		double d_lambda_error;

		double d_k;
		double d_k_error;

		double d_a1;
		double d_a1_error;

		double d_b1;
		double d_b1_error;

		double d_a2;
		double d_a2_error;

		double d_b2;
		double d_b2_error;

		double d_alpha1;
		double d_alpha1_error;

		double d_beta1;
		double d_beta1_error;

		double d_alpha2;
		double d_alpha2_error;

		double d_beta2;
		double d_beta2_error;

		double d_sigma;
		double d_sigma_error;

		double d_tau;
		double d_tau_error;
	};


	class sls_basic{

	public:

		template<class T>
		static inline T Tmax(T i_, T j_)
		{
			if(i_>j_)
			{
				return i_;
			};
			return j_;
		}

		template<class T>
		static inline T Tmin(T i_, T j_)
		{
			if(i_<j_)
			{
				return i_;
			};
			return j_;
		}

		template<class T>
		static inline T Tmax(T x_,T y_,T z_)
		{
			return Tmax(Tmax(x_,y_),z_);
		}

		template<class T>
		static inline T Tmin(T x_,T y_,T z_)
		{
			return Tmin(Tmin(x_,y_),z_);
		}

		template<class T>
		static inline T Tmax(T x_,T y_,T z_,T w_)
		{
			return Tmax(Tmax(x_,y_),Tmax(z_,w_));
		}

		template<class T>
		static inline T Tmin(T x_,T y_,T z_,T w_)
		{
			return Tmin(Tmin(x_,y_),Tmin(z_,w_));
		}

		static inline void assert_mem(void *pointer_)
		{
			if(!pointer_)
			{
				throw error("Memory allocation error\n",41);
			};
		}

		static double round(//returns nearest integer to x_
		const double &x_);

		static void get_current_time(
		double &seconds_);

		static long int random_seed_from_time();

		static double one_minus_exp_function(
		double y_);

		static double normal_probability(double x_)
		{
			return 0.5*erfc(-sqrt(0.5)*x_);
		}

		static double normal_probability(
		double x_,
		double eps_);

		static double normal_probability(
		double a_,
		double b_,
		double h_,
		long int N_,
		double *p_,
		double x_,
		double eps_);

		static double ln_one_minus_val(
		double val_);



	};
}

#endif

