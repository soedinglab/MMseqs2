#ifndef INCLUDED_SLS_ALP_REGRESSION
#define INCLUDED_SLS_ALP_REGRESSION

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

File name: sls_alp_regression.hpp

Author: Sergey Sheetlin

Contents: Regression methods

******************************************************************************/
#include "sls_basic.hpp"

#include <complex>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <float.h>
#include <algorithm>

namespace Sls { 

	typedef double function_type(double x_,void* func_number_);


	class alp_reg{

	

	public:


		alp_reg(//constructor
			);


		~alp_reg();//destructor

		static void find_tetta_general(
		function_type *func_,
		void* func_pointer_,
		double a_,//[a,b] is the interval for search of equation solution
		double b_,
		long int n_partition_,
		double eps_,
		std::vector<double> &res_);

		static double find_single_tetta_general(
		function_type *func_,
		void* func_pointer_,
		double a_,//[a,b] is the interval for search of equation solution
		double b_,
		double eps_);

		static void correction_of_errors(
		double *errors_,
		long int number_of_elements_);


		static void robust_regression_sum_with_cut_LSM(
		long int min_length_,
		long int number_of_elements_,
		double *values_,
		double *errors_,
		bool cut_left_tail_,
		bool cut_right_tail_,
		double y_,
		double &beta0_,
		double &beta1_,
		double &beta0_error_,
		double &beta1_error_,
		long int &k1_opt_,
		long int &k2_opt_,
		bool &res_was_calculated_);

		static double function_for_robust_regression_sum_with_cut_LSM(
		double *values_,
		double *errors_,
		long int number_of_elements_,
		long int k_start_,
		double c_,
		double &beta0_,
		double &beta1_,
		double &beta0_error_,
		double &beta1_error_,
		bool &res_was_calculated_);

		static void robust_regression_sum_with_cut_LSM_beta1_is_defined(
		long int min_length_,
		long int number_of_elements_,
		double *values_,
		double *errors_,
		bool cut_left_tail_,
		bool cut_right_tail_,
		double y_,
		double &beta0_,
		double beta1_,
		double &beta0_error_,
		double beta1_error_,
		long int &k1_opt_,
		long int &k2_opt_,
		bool &res_was_calculated_);

		static double function_for_robust_regression_sum_with_cut_LSM_beta1_is_defined(
		double *values_,
		double *errors_,
		long int number_of_elements_,
		long int k_start_,
		double c_,
		double &beta0_,
		double beta1_,
		double &beta0_error_,
		double beta1_error_,
		bool &res_was_calculated_);


		inline static double sqrt_for_errors(
			double x_)
			{
				if(x_<=0)
				{
					return 0.0;
				}
				else
				{
					return sqrt(x_);
				};
			}


		static double median(
		long int dim_,
		double *array_);

		static double robust_sum(
		double *values,
		long int dim,
		long int N_points,
		bool *&remove_flag);


	};
}


#endif

