#ifndef INCLUDED_SLS_ALP_SIMULATION
#define INCLUDED_SLS_ALP_SIMULATION

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

File name: sls_alp_sim.hpp

Author: Sergey Sheetlin

Contents: Simulation of Gumbel parameters

******************************************************************************/


#include <complex>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <float.h>
#include <algorithm>

#include "sls_alp_data.hpp"
#include "sls_alp.hpp"
#include "sls_alp_regression.hpp"

namespace Sls {

	struct struct_for_lambda_calculation
	{
		void **d_alp_distr;
		void **d_alp_distr_errors;
		long int d_nalp;
		double d_f_error;

		double d_last_sum;
		double d_last_sum_error;

		bool d_calculate_alp_number;
		long int d_alp_number;
	};




	class alp_sim{

	

	public:


		alp_sim(//constructor
			alp_data *alp_data_
			);


		~alp_sim();//destructor

		void alp_sim_from_random_seed(//simulation from random seed
		alp_data *alp_data_);

		void quick_test(//runs quick tests to determine whether the scoring scheme is linear
		long int trials_number_,
		double max_time_);

		void get_minimal_simulation(
		long int ind1_,
		long int ind2_,
		long int &M_min_,
		long int &nalp_,
		long int &nalp_lambda_,
		bool C_calculation_,
		bool check_time_flag_);//simulation using [ind1_,ind2_] range of realizations with an estimation of parameters for accuracy, memory usage and calculation time


		void memory_release_for_get_minimal_simulation(
		long int nalp_,
		void **&alp_distr,
		void **&alp_distr_errors);

		bool the_criterion(//criteria of stopping of the simulating ALP
		//if the function returns true then calculates optimal M_min and ALP number
		//sets the flags M_min_flag_ and nalp_flag_ checking the corresponding condition
		long int upto_nalp_,
		long int &nalp_for_lambda_simulation_,
		long int ind1_,
		long int ind2_,
		void **&alp_distr,
		void **&alp_distr_errors,
		long int &M_min_,
		bool &M_min_flag_,
		bool &nalp_flag_,
		bool &inside_simulation_flag_,
		bool C_calculation_,
		double *lambda_=NULL,
		double *lambda_error_=NULL);

		void calculate_lambda(
		bool check_the_criteria_,
		long int nalp_,
		long int &nalp_thr_,
		bool &inside_simulation_flag_,
		void **alp_distr,
		void **alp_distr_errors,
		double &lambda_,
		double &lambda_error_,
		double &test_difference_,
		double &test_difference_error_);

		void calculate_C(
		long int starting_point,
		long int nalp_,
		void **alp_distr,
		void **alp_distr_errors,
		double lambda_,
		double lambda_error_,
		double &C_,
		double &C_error_,
		double &Sc_,
		double &Sc_error_);

		void memory_release_for_calculate_FSC(
		double *&exp_array,

		double *&delta_E,
		double *&delta_E_error,

		double *&delta_E_E,
		double *&delta_E_E_error,


		double *&delta_I,
		double *&delta_I_error,

		double *&delta_J,
		double *&delta_J_error,

		double *&delta_I_I,
		double *&delta_I_I_error,

		double *&delta_I_J,
		double *&delta_I_J_error,

		double *&delta_J_J,
		double *&delta_J_J_error,

		double *&cov_J_J,
		double *&cov_J_J_error,

		double *&cov_I_J,
		double *&cov_I_J_error,

		double *&cov_I_I,
		double *&cov_I_I_error,

		double *&cov_E_E,
		double *&cov_E_E_error);


		void calculate_FSC(
		long int nalp_,
		long int ind1_,
		long int ind2_,
		void **alp_distr,
		double lambda_,
		double Sc_,
		//double Sc_error_,

		double &a_I_,
		double &a_I_error_,
		double &a_J_,
		double &a_J_error_,
		double &sigma_,
		double &sigma_error_,
		double &alpha_I_,
		double &alpha_I_error_,
		double &alpha_J_,
		double &alpha_J_error_);

		void sigma_calculation(
		double delta_I_aver_,
		double delta_I_aver_error_,
		double delta_J_aver_,
		double delta_J_aver_error_,
		double delta_E_aver_,
		double delta_E_aver_error_,
		double cov_E_E_aver_,
		double cov_E_E_aver_error_,
		double cov_I_J_aver_,
		double cov_I_J_aver_error_,
		double &sigma_,
		double &sigma_error_);

		void get_and_allocate_alp_distribution(
		long int ind1_,
		long int ind2_,
		void **&alp_distr,
		void **&alp_distr_errors,
		long int nalp);

		bool check_K_criterion(
		long int nalp_,
		long int ind1_,
		long int ind2_,
		double lambda_,
		double eps_K_,
		long int &M_min_);

		void kill(
		bool check_time_,
		long int ind1_,
		long int ind2_,
		long int M_min_,
		double lambda_,
		double eps_K_,
		double &K_C_,
		double &K_C_error_,
		long int &level_,
		long int &diff_opt_);


		
		bool check_K_criterion_during_killing2(
		long int ind1_,
		long int ind2_,
		double lambda_,
		double eps_K_,
		long int current_level_,
		long int &recommended_level_,
		long int &diff_opt_,
		double &K_C_,
		double &K_C_error_);
		

		bool check_K_criterion_during_killing(
		long int ind1_,
		long int ind2_,
		double lambda_,
		double eps_K_,
		long int current_level_,
		long int &recommended_level_,
		long int &diff_opt_,
		double &K_C_,
		double &K_C_error_);




		inline static double lambda_exp(
		long int &i_,
		double *&exp_array_);

		void get_single_realization(
		bool check_time_,
		long int M_min_,
		long int nalp_,
		bool killing_flag_,
		long int level_,
		long int diff_opt_,
		alp *&obj_,
		bool &sucess_flag_,
		double &d_eps_);

		void calculate_main_parameters(
		long int final_realizations_number_lambda_,
		long int final_realizations_number_killing_,
		long int nalp,
		long int nalp_for_lambda_simulation,
		long int level,
		bool &inside_simulation_flag,
		double &lambda,
		double &lambda_error,
		double &test_difference,
		double &test_difference_error,
		double &C,
		double &C_error,
		double &C2,
		double &C2_error,
		double &C4,
		double &C4_error,
		double &K_C,
		double &K_C_error,
		double &a_I,
		double &a_I_error,
		double &a_J,
		double &a_J_error,
		double &sigma,
		double &sigma_error,
		double &alpha_I,
		double &alpha_I_error,
		double &alpha_J,
		double &alpha_J_error,
		double &K,
		double &K_error);

		void calculate_main_parameters2(
		long int final_realizations_number_lambda_,
		long int final_realizations_number_killing_,
		long int nalp,
		long int nalp_for_lambda_simulation,
		long int level,
		bool &inside_simulation_flag,
		double &lambda,
		double &lambda_error,
		double &test_difference,
		double &test_difference_error,
		double &C,
		double &C_error,
		double &C2,
		double &C2_error,
		double &C4,
		double &C4_error,
		double &K_C,
		double &K_C_error,
		double &a_I,
		double &a_I_error,
		double &a_J,
		double &a_J_error,
		double &sigma,
		double &sigma_error,
		double &alpha_I,
		double &alpha_I_error,
		double &alpha_J,
		double &alpha_J_error,
		double &K,
		double &K_error);

		void memory_release_for_calculate_main_parameters2m(
		long int nalp_for_lambda_simulation,
		long int *&d_mult_realizations,
		long int *&d_mult_K_realizations,

		double *&lambda_mult,
		double *&lambda_mult_error,

		double *&C_mult,
		double *&C_mult_error,

		double *&a_I_mult,
		double *&a_I_mult_error,

		double *&a_J_mult,
		double *&a_J_mult_error,

		double *&sigma_mult,
		double *&sigma_mult_error,

		double *&alpha_I_mult,
		double *&alpha_I_mult_error,

		double *&alpha_J_mult,
		double *&alpha_J_mult_error,

		double *&K_C_mult,
		double *&K_C_mult_error,

		double *&K_mult,
		double *&K_mult_error,

		double *&Sc_mult,
		double *&Sc_mult_error,


		void **&alp_distr,
		void **&alp_distr_errors,

		void ***&alp_mult_distr,
		void ***&alp_mult_distr_errors);

		void calculate_main_parameters2m(
		long int final_realizations_number_lambda_,
		long int final_realizations_number_killing_,
		long int nalp_for_lambda_simulation,
		long int level,
		bool &inside_simulation_flag,
		double &lambda,
		double &lambda_error,
		double &test_difference,
		double &test_difference_error,
		double &C,
		double &C_error,
		double &K_C,
		double &K_C_error,
		double &a_I,
		double &a_I_error,
		double &a_J,
		double &a_J_error,
		double &sigma,
		double &sigma_error,
		double &alpha_I,
		double &alpha_I_error,
		double &alpha_J,
		double &alpha_J_error,
		double &K,
		double &K_error,
		bool &flag_);

		void randomize_realizations(
		long int final_realizations_number_lambda_,
		long int final_realizations_number_killing_);

		void randomize_realizations_ind(
		long int ind1_,
		long int ind2_);

		void generate_random_permulation(
		long int *perm_,
		long int dim_);




		static void error_in_calculate_main_parameters2m(
		double C,
		double &C_error,
		double C_mult2,
		double C_mult2_error);




		void output_main_parameters(
		double time_,
		long int nalp,
		long int nalp_for_lambda_simulation,
		long int level,
		long int M_min_,
		bool &inside_simulation_flag,
		long int final_realizations_number_lambda_,
		long int final_realizations_number_killing_);

		void output_main_parameters2(
		double time_,
		long int nalp,
		long int nalp_for_lambda_simulation,
		long int level,
		long int M_min_,
		bool &inside_simulation_flag,
		long int final_realizations_number_lambda_,
		long int final_realizations_number_killing_);

		void output_main_parameters2m_new(
		long int nalp_for_lambda_simulation,
		long int level,
		bool &inside_simulation_flag,
		long int final_realizations_number_lambda_,
		long int final_realizations_number_killing_);

		void symmetric_parameters_for_symmetric_scheme();




		static double relative_error_in_percents(
		double val_,
		double val_error_);

		static double round_double(
		double val_,
		long int digits_);

		static long int get_number_of_subsimulations(
		long int number_of_realizations_);



		static double function_for_lambda_calculation(
		double lambda_,
		void * data_);

		static double get_root(
		const std::vector<double> &res_tmp_,
		double point_);
		



	public:


				
		alp_data *d_alp_data;//initial data
		array_positive<alp*> *d_alp_obj;//vector with the alp objects
		long int d_n_alp_obj;//number of alp objects
					
		array_positive<double> *d_lambda_tmp;
		array_positive<double> *d_lambda_tmp_errors;

		array_positive<double> *d_C_tmp;
		array_positive<double> *d_C_tmp_errors;

		//Subsimulations' parameters

		//number of subsimulations
		long int d_mult_number;


		//parameters estimations

		double m_Lambda;
		double m_LambdaError;
		double m_K;
		double m_KError;
		double m_C;
		double m_CError;

		double m_Sigma;
		double m_SigmaError;

		double m_GaplessAlpha;
		double m_GaplessAlphaError;

		double m_GaplessA;
		double m_GaplessAError;


		double m_AlphaI;
		double m_AlphaIError;
		double m_AlphaJ;
		double m_AlphaJError;

		double m_AI;
		double m_AIError;
		double m_AJ;
		double m_AJError;



		double m_CalcTime;

		long int m_G;
		long int m_G1;
		long int m_G2;

		std::vector<double> m_LambdaSbs;
		std::vector<double> m_KSbs;
		std::vector<double> m_CSbs;

		std::vector<double> m_SigmaSbs;

		std::vector<double> m_AlphaISbs;
		std::vector<double> m_AlphaJSbs;

		std::vector<double> m_AISbs;
		std::vector<double> m_AJSbs;



	};
}

#endif

