#ifndef INCLUDED_SLS_PVALUES
#define INCLUDED_SLS_PVALUES

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

File name: sls_pvalues.hpp

Author: Sergey Sheetlin

Contents: P-values calculation routines

******************************************************************************/

#include "sls_basic.hpp"

#include <vector>
#include <string>
#include <math.h>
#include <cstdlib>

#include <cmath>
#include <iostream>

namespace Sls {

	struct ALP_set_of_parameters
	{
		ALP_set_of_parameters()
		{
			d_params_flag=false;

			b_I=0;
			b_I_error=0;

			b_J=0;
			b_J_error=0;

			beta_I=0;
			beta_I_error=0;

			beta_J=0;
			beta_J_error=0;

			tau=0;
			tau_error=0;

		};

		double lambda;
		double lambda_error;

		double C;
		double C_error;


		double K;
		double K_error;

		double a_I;
		double a_I_error;

		double a_J;
		double a_J_error;

		double sigma;
		double sigma_error;

		double alpha_I;
		double alpha_I_error;

		double alpha_J;
		double alpha_J_error;

		double a;
		double a_error;

		double alpha;
		double alpha_error;

		double gapless_a;
		double gapless_a_error;

		double gapless_alpha;
		double gapless_alpha_error;

		long int G;
		long int G1;
		long int G2;

		std::vector<double > m_LambdaSbs;
		std::vector<double > m_KSbs;
		std::vector<double > m_CSbs;

		std::vector<double > m_SigmaSbs;

		std::vector<double > m_AlphaISbs;
		std::vector<double > m_AlphaJSbs;

		std::vector<double > m_AISbs;
		std::vector<double > m_AJSbs;

		double m_CalcTime;

		bool d_params_flag;//if true, then the parameters are defined and P-values can be calculated

		//intercepts
		double b_I;
		double b_I_error;

		double b_J;
		double b_J_error;

		double beta_I;
		double beta_I_error;

		double beta_J;
		double beta_J_error;

		double tau;
		double tau_error;

		std::vector<double > m_BISbs;
		std::vector<double > m_BJSbs;

		std::vector<double > m_BetaISbs;
		std::vector<double > m_BetaJSbs;

		std::vector<double > m_TauSbs;

		//tmp values
		double vi_y_thr;
		double vj_y_thr;
		double c_y_thr;
	};

	std::ostream &operator<<(std::ostream &s_,
	const ALP_set_of_parameters &gumbel_params_);

	std::istream &operator>>(std::istream &s_,
	ALP_set_of_parameters &gumbel_params_);

	class pvalues{

		public:


		pvalues();

		~pvalues();


		public:

		static void get_appr_tail_prob_with_cov(
		const ALP_set_of_parameters &par_,
		bool blast_,
		double y_,
		double m_,
		double n_,

		double &P_,
		double &P_error_,

		double &E_,
		double &E_error_,

		double &area_,

		double a_normal_,
		double b_normal_,
		double h_normal_,
		long int N_normal_,
		double *p_normal_,

		bool &area_is_1_flag_);

		static void compute_tmp_values(ALP_set_of_parameters &par_);

		static void get_appr_tail_prob_with_cov_without_errors(
		const ALP_set_of_parameters &par_,
		bool blast_,
		double y_,
		double m_,
		double n_,

		double &P_,

		double &E_,

		double &area_,

		double a_normal_,
		double b_normal_,
		double h_normal_,
		long int N_normal_,
		double *p_normal_,

		bool &area_is_1_flag_,
		bool compute_only_area_=false);

		static void get_P_error_using_splitting_method(
		const ALP_set_of_parameters &par_,
		bool blast_,
		double y_,
		double m_,
		double n_,

		double &P_,
		double &P_error_,

		double &E_,
		double &E_error_,

		double a_normal_,
		double b_normal_,
		double h_normal_,
		long int N_normal_,
		double *p_normal_,

		bool &area_is_1_flag_);


		public:

		static void compute_intercepts(
		ALP_set_of_parameters &par_);

		void calculate_P_values(
		long int Score1,
		long int Score2,
		double Seq1Len,
		double Seq2Len,
		const ALP_set_of_parameters &ParametersSet,
		std::vector<double> &P_values,
		std::vector<double> &P_values_errors,
		std::vector<double> &E_values,
		std::vector<double> &E_values_errors);

		void calculate_P_values(
		double Score,
		double Seq1Len,
		double Seq2Len,
		const ALP_set_of_parameters &ParametersSet,
		double &P_value,
		double &P_value_error,
		double &E_value,
		double &E_value_error,
		bool read_Sbs_par_flag=true);

		static inline double ran3()//generates the next random value
		{
			double rand_C=(double)((double)rand()/(double)RAND_MAX);
			return rand_C;	
		};

		static inline double standard_normal()//generates standard normal random value using the Box–Muller transform
		{
			double r1=0;
			while(r1==0)
			{
				r1=ran3();
			};
			double r2=0;
			while(r2==0)
			{
				r2=ran3();
			};

			double v1=-2*log(r1);
			if(v1<0)
			{
				v1=0;
			};
			return sqrt(v1)*cos(2*pi*r2);
		};

		static //returns "true" if the Gumbel parameters are properly defined and "false" otherwise
		bool assert_Gumbel_parameters(
		const ALP_set_of_parameters &par_);//a set of Gumbel parameters





		public:


		bool blast;
		double eps;
		double a_normal;
		double b_normal;
		long int N_normal;
		double h_normal;
		double *p_normal;


	};
}

#endif //! INCLUDED

