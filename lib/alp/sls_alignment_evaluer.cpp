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

File name: sls_alignment_evaluer.cpp

Author: Sergey Sheetlin

Contents: library functions of main routines

******************************************************************************/


#include "sls_alignment_evaluer.hpp"

#include "sls_alp.hpp"
#include "sls_alp_data.hpp"
#include "sls_alp_regression.hpp"
#include "sls_alp_sim.hpp"
#include "njn_localmaxstatmatrix.hpp"
#include "njn_localmaxstatutil.hpp"

using namespace std;


namespace Sls {

// Write the parameters:
std::ostream &operator<<(std::ostream &s_,
const AlignmentEvaluer &g_)
{

	if(!pvalues::assert_Gumbel_parameters(
	g_.d_params)||!g_.isGood())
	{
		throw error("Error - the Gumbel parameters are not defined properly in the function \"std::ostream &operator<<\"\n",1);
	};

	s_<<g_.d_params;
	return s_;
}

// Read the parameters:
std::istream &operator>>(std::istream &s_,
AlignmentEvaluer &g_)
{
	try
	{
		g_.d_params.d_params_flag=false;
		s_>>g_.d_params;
		g_.d_params.d_params_flag=true;

		//precompute intercepts
		pvalues::compute_intercepts(g_.d_params);

		if(!pvalues::assert_Gumbel_parameters(
		g_.d_params)||!g_.isGood())
		{
			g_.d_params.d_params_flag=false;
		};

		return s_;
	}
	catch (...)
	{ 
		g_.d_params.d_params_flag=false;
		throw;
	};
}

//check correctness of the input parameters for gapless alignment
void AlignmentEvaluer::assert_Gapless_input_parameters(
long alphabetSize_,//a number of letters in the alphabet
const double *letterFreqs1_,//background frequencies of letters in sequence #1
const double *letterFreqs2_,//background frequencies of letters in sequence #2
double *&letterFreqs1_normalized_,//normalized background frequencies of letters in sequence #1
double *&letterFreqs2_normalized_,//normalized background frequencies of letters in sequence #2
const string function_name_)//"assert_Gapless_input_parameters" is called from "function_name_" function
{
	if(!(alphabetSize_>0))
	{
		d_params.d_params_flag=false;
		throw error("Error - the parameter \"alphabetSize_\" in the function \""+function_name_+"\" must be positive\n",1);
	};

	long int i;

	double sum1=0;
	for(i=0;i<alphabetSize_;i++)
	{
		if(letterFreqs1_[i]<0)
		{
			d_params.d_params_flag=false;
			throw error("Error - the value \"letterFreqs1_["+alp_data::long_to_string(i)+"]\" in the function \""+function_name_+"\" must be non-negative\n",1);
		};
		sum1+=letterFreqs1_[i];
	};

	if(sum1<=0)
	{
		throw error("Error - sum of the frequencies \"letterFreqs1_\" is non-positive in the function \""+function_name_+"\"\n",1);
	};

	letterFreqs1_normalized_=new double[alphabetSize_];
	alp_data::assert_mem(letterFreqs1_normalized_);


	for(i=0;i<alphabetSize_;i++)
	{
		letterFreqs1_normalized_[i]=letterFreqs1_[i]/sum1;
	};

	double sum2=0;
	for(i=0;i<alphabetSize_;i++)
	{
		if(letterFreqs2_[i]<0)
		{
			d_params.d_params_flag=false;
			throw error("Error - the value \"letterFreqs2_["+alp_data::long_to_string(i)+"]\" in the function \""+function_name_+"\" must be non-negative\n",1);
		};
		sum2+=letterFreqs2_[i];
	};

	if(sum2<=0)
	{
		throw error("Error - sum of the frequencies \"letterFreqs2_\" is non-positive in the function \""+function_name_+"\"\n",1);
	};
	letterFreqs2_normalized_=new double[alphabetSize_];
	alp_data::assert_mem(letterFreqs1_normalized_);

	for(i=0;i<alphabetSize_;i++)
	{
		letterFreqs2_normalized_[i]=letterFreqs2_[i]/sum2;
	};

}


//Computes gapless Gumbel parameters:
void AlignmentEvaluer::initGapless(long alphabetSize_,//a number of letters in the alphabet
const long *const *substitutionScoreMatrix_,//scoring matrix
const double *letterFreqs1_,//background frequencies of letters in sequence #1
const double *letterFreqs2_,//background frequencies of letters in sequence #2
double max_time_)//maximum allowed calculation time in seconds
{

	try
	{
		

		double CurrentTime1;
		Sls::alp_data::get_current_time(CurrentTime1);

		//check correctness of the input parameters for gapless alignment
		string function_name="void AlignmentEvaluer::initGapless";
		double *letterFreqs1_normalized=NULL;//normalized background frequencies of letters in sequence #1
		double *letterFreqs2_normalized=NULL;//normalized background frequencies of letters in sequence #2

		assert_Gapless_input_parameters(
		alphabetSize_,//a number of letters in the alphabet
		letterFreqs1_,//background frequencies of letters in sequence #1
		letterFreqs2_,//background frequencies of letters in sequence #2
		letterFreqs1_normalized,//normalized background frequencies of letters in sequence #1
		letterFreqs2_normalized,//normalized background frequencies of letters in sequence #2
		function_name);//"assert_Gapless_input_parameters" is called from "function_name_" function

		
		if(max_time_<=0)
		{
			max_time_=60;
		};

		d_params.d_params_flag=false;

		Njn::LocalMaxStatMatrix local_max_stat_matrix(alphabetSize_,
							  substitutionScoreMatrix_,
							  letterFreqs1_normalized,
							  letterFreqs2_normalized,
							  alphabetSize_,
							  max_time_);

		if(local_max_stat_matrix.getTerminated()) 
		{
			throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
		};


		//calculation of a and sigma
		double calculation_error=1e-6;

		d_params.gapless_alpha = local_max_stat_matrix.getAlpha ();
		d_params.gapless_alpha=alp_data::Tmax(d_params.gapless_alpha,0.0);
		d_params.gapless_alpha_error = calculation_error;

		d_params.gapless_a = local_max_stat_matrix.getA ();
		d_params.gapless_a=alp_data::Tmax(d_params.gapless_a,0.0);
		d_params.gapless_a_error = calculation_error;

		//calculation of all required parameters for a gapless case
		d_params.G=0;
		d_params.G1=0;
		d_params.G2=0;

		d_params.lambda = local_max_stat_matrix.getLambda ();
		d_params.lambda_error = calculation_error;

		d_params.K = local_max_stat_matrix.getK ();
		d_params.K_error = calculation_error;
			
		d_params.C = local_max_stat_matrix.getC ();;
		d_params.C_error = calculation_error;

		d_params.sigma = d_params.gapless_alpha;
		d_params.sigma_error = calculation_error;

		d_params.alpha_I = d_params.gapless_alpha;
		d_params.alpha_I_error = calculation_error;

		d_params.alpha_J = d_params.gapless_alpha;
		d_params.alpha_J_error = calculation_error;

		d_params.a_I = d_params.gapless_a;
		d_params.a_I_error = calculation_error;

		d_params.a_J = d_params.gapless_a;
		d_params.a_J_error = calculation_error;


		std::vector<double > sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.lambda;
		sbs_arrays[1]=d_params.lambda + calculation_error;

		d_params.m_LambdaSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.K;
		sbs_arrays[1]=d_params.K+calculation_error;

		d_params.m_KSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.C;
		sbs_arrays[1]=d_params.C+calculation_error;

		d_params.m_CSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.sigma;
		sbs_arrays[1]=d_params.sigma + calculation_error;

		d_params.m_SigmaSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.alpha_I;
		sbs_arrays[1]=d_params.alpha_I + calculation_error;

		d_params.m_AlphaISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.alpha_J;
		sbs_arrays[1]=d_params.alpha_J + calculation_error;

		d_params.m_AlphaJSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.a_I;
		sbs_arrays[1]=d_params.a_I + calculation_error;

		d_params.m_AISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.a_J;
		sbs_arrays[1]=d_params.a_J + calculation_error;

		d_params.m_AJSbs=sbs_arrays;

		d_params.a 
			= (d_params.a_I + d_params.a_J) * 0.5;

		d_params.a_error = (d_params.a_I_error
						+ d_params.a_J_error)*0.5;
		
		d_params.alpha = (d_params.alpha_I
						+ d_params.alpha_J) * 0.5;

		d_params.alpha_error = (d_params.alpha_I_error
						+ d_params.alpha_J_error) * 0.5;

		d_params.d_params_flag=true;

		//precompute intercepts
		pvalues::compute_intercepts(d_params);

		if(!pvalues::assert_Gumbel_parameters(
		d_params)||!isGood())
		{
			d_params.d_params_flag=false;
			throw error("Error - computation of the Gumbel parameters is unsuccessful in the function \"void AlignmentEvaluer::initGapless\"\n",1);
		};

		delete[]letterFreqs1_normalized;
		delete[]letterFreqs2_normalized;

		double CurrentTime2;
		Sls::alp_data::get_current_time(CurrentTime2);
		d_params.m_CalcTime=CurrentTime2-CurrentTime1;

	}
	catch (...)
	{ 
		d_params.d_params_flag=false;
		throw;
	};

}

//Computes gapped Gumbel parameters:
//The NCBI convention is used for penalizing a gap:
//For example, a gap of length k is penalized as gapOpen1_+k*gapEpen1_ for sequence #1
//if max_time_<=0 then d_gapped_computation_parameters will be used;
//d_gapped_computation_parameters.d_parameters_flag must be true in this case
//every execution of initGapped with max_time_>0 rewrites d_gapped_computation_parameters by the actual values
void AlignmentEvaluer::initGapped(long alphabetSize_,//a number of letters in the alphabet
const long *const *substitutionScoreMatrix_,//scoring matrix
const double *letterFreqs1_,//background frequencies of letters in sequence #1
const double *letterFreqs2_,//background frequencies of letters in sequence #2
long gapOpen1_,//gap opening penalty for sequence #1
long gapEpen1_,//gap extension penalty for sequence #1
long gapOpen2_,//gap opening penalty for sequence #2
long gapEpen2_,//gap extension penalty for sequence #2
bool insertions_after_deletions_,//if true, then insertions after deletions are permitted
double eps_lambda_,//relative error for the parameter lambda
double eps_K_,//relative error for the parameter K
double max_time_,//maximum allowed calculation time in seconds
double max_mem_,//maximum allowed memory usage in Mb
long randomSeed_,//randomizaton seed
double temperature_)
{

	struct_for_randomization *randomization_parameters=NULL;
	try
	{

		double CurrentTime1;
		Sls::alp_data::get_current_time(CurrentTime1);

		//check correctness of the input parameters for gapless alignment
		string function_name="void AlignmentEvaluer::initGapped";
		double *letterFreqs1_normalized=NULL;//normalized background frequencies of letters in sequence #1
		double *letterFreqs2_normalized=NULL;//normalized background frequencies of letters in sequence #2
		assert_Gapless_input_parameters(
		alphabetSize_,//a number of letters in the alphabet
		letterFreqs1_,//background frequencies of letters in sequence #1
		letterFreqs2_,//background frequencies of letters in sequence #2
		letterFreqs1_normalized,//normalized background frequencies of letters in sequence #1
		letterFreqs2_normalized,//normalized background frequencies of letters in sequence #2
		function_name);//"assert_Gapless_input_parameters" is called from "function_name_" function

		if(!(gapEpen1_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"gapEpen1_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(gapEpen2_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"gapEpen2_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(eps_lambda_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"eps_lambda_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(eps_K_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"eps_K_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		if(!(max_mem_>0))
		{
			d_params.d_params_flag=false;
			throw error("Error - the parameter \"max_mem_\" in the function \""+function_name+"\" must be positive\n",1);
		};

		d_params.d_params_flag=false;

		//Gapless parameters calculation
		double GaplessTimePortion=0.5;
		double GaplessCalculationTime=max_time_;

		if(max_time_<=0)
		{
			GaplessCalculationTime=120;//the time is set to 120 seconds if not set as an input
		};

		//Gapless calculation may take only a portion of maximum allowed calculation time in the case of gapped calculation 
		GaplessCalculationTime*=GaplessTimePortion;


		
		Njn::LocalMaxStatMatrix local_max_stat_matrix(alphabetSize_,
							  substitutionScoreMatrix_,
							  letterFreqs1_normalized,
							  letterFreqs2_normalized,
							  alphabetSize_,
							  GaplessCalculationTime);

		if(local_max_stat_matrix.getTerminated()) 
		{
			throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
		};

		//calculation of a and sigma
		double calculation_error=1e-6;

		d_params.gapless_alpha = local_max_stat_matrix.getAlpha ();
		d_params.gapless_alpha=alp_data::Tmax(d_params.gapless_alpha,0.0);
		d_params.gapless_alpha_error = calculation_error;

		d_params.gapless_a = local_max_stat_matrix.getA ();
		d_params.gapless_a=alp_data::Tmax(d_params.gapless_a,0.0);
		d_params.gapless_a_error = calculation_error;


		double CurrentTimeGaplessPreliminary;
		Sls::alp_data::get_current_time(CurrentTimeGaplessPreliminary);
		double GaplessPreliminaryTime=CurrentTimeGaplessPreliminary-CurrentTime1;

		//the choice for the importance sampling
		//long int gapOpen=alp_data::Tmin(gapOpen1_,gapOpen2_);
		//long int gapEpen=alp_data::Tmin(gapEpen1_,gapEpen2_);

		long int gapEpen = alp_data::Tmin(gapEpen1_, gapEpen2_);
		long int gapOpen = alp_data::Tmin(gapOpen1_ + gapEpen1_, gapOpen2_ + gapEpen2_) - gapEpen;


		if(max_time_<=0)
		{
			if(d_gapped_computation_parameters.d_parameters_flag)
			{
				randomization_parameters=new struct_for_randomization;

				randomization_parameters->d_first_stage_preliminary_realizations_numbers_ALP=
					d_gapped_computation_parameters.d_first_stage_preliminary_realizations_numbers_ALP;

				randomization_parameters->d_preliminary_realizations_numbers_ALP=
					d_gapped_computation_parameters.d_preliminary_realizations_numbers_ALP;

				randomization_parameters->d_preliminary_realizations_numbers_killing=
					d_gapped_computation_parameters.d_preliminary_realizations_numbers_killing;

				randomization_parameters->d_random_seed=randomSeed_;

				randomization_parameters->d_total_realizations_number_with_ALP=
					d_gapped_computation_parameters.d_total_realizations_number_with_ALP;

				randomization_parameters->d_total_realizations_number_with_killing=
					d_gapped_computation_parameters.d_total_realizations_number_with_killing;

			}
			else
			{
				throw error("Error - d_gapped_computation_parameters must be defined before calling AlignmentEvaluer::initGapped with max_time_<=0\n",1);
			};
		};

		Sls::alp_data data_obj(//constructor
		randomSeed_,//randomization number
		randomization_parameters,//if not NULL, sets d_rand_flag to true and initializes d_rand_all

		gapOpen,//gap opening penalty
		gapOpen1_,//gap opening penalty for a gap in the sequence #1
		gapOpen2_,//gap opening penalty for a gap in the sequence #2

		gapEpen,//gap extension penalty
		gapEpen1_,//gap extension penalty for a gap in the sequence #1
		gapEpen2_,//gap extension penalty for a gap in the sequence #2

		alphabetSize_,
		substitutionScoreMatrix_,
		letterFreqs1_normalized,
		letterFreqs2_normalized,

		temperature_,
		max_time_,//maximum allowed calculation time in seconds
		max_mem_,//maximum allowed memory usage in MB
		eps_lambda_,//relative error for lambda calculation
		eps_K_,//relative error for K calculation
		insertions_after_deletions_,//if true, then insertions after deletions are allowed
		d_gapped_computation_parameters.d_max_time_for_quick_tests,//maximum allowed calculation time in seconds for quick tests
		d_gapped_computation_parameters.d_max_time_with_computation_parameters);//maximum allowed time in seconds for the whole computation



		data_obj.d_max_time=Sls::alp_data::Tmax((1.0-GaplessTimePortion)*data_obj.d_max_time,data_obj.d_max_time-GaplessPreliminaryTime);

		Sls::alp_sim sim_obj(&data_obj);

		if(max_time_>0)
		{
			d_gapped_computation_parameters.d_parameters_flag=true;
			
			d_gapped_computation_parameters.d_first_stage_preliminary_realizations_numbers_ALP=
				sim_obj.d_alp_data->d_rand_all->d_first_stage_preliminary_realizations_numbers_ALP;

			d_gapped_computation_parameters.d_preliminary_realizations_numbers_ALP=
				sim_obj.d_alp_data->d_rand_all->d_preliminary_realizations_numbers_ALP;

			d_gapped_computation_parameters.d_preliminary_realizations_numbers_killing=
				sim_obj.d_alp_data->d_rand_all->d_preliminary_realizations_numbers_killing;

			d_gapped_computation_parameters.d_total_realizations_number_with_ALP=
				sim_obj.d_alp_data->d_rand_all->d_total_realizations_number_with_ALP;

			d_gapped_computation_parameters.d_total_realizations_number_with_killing=
				sim_obj.d_alp_data->d_rand_all->d_total_realizations_number_with_killing;
		};

		sim_obj.m_GaplessAlpha = d_params.gapless_alpha;
		sim_obj.m_GaplessAlphaError = d_params.gapless_alpha_error;

		sim_obj.m_GaplessA = d_params.gapless_a;
		sim_obj.m_GaplessAError = d_params.gapless_a_error;

		
		sim_obj.m_G1=gapOpen1_+gapEpen1_;
		sim_obj.m_G2=gapOpen2_+gapEpen2_;
		sim_obj.m_G=alp_data::Tmin(sim_obj.m_G1,sim_obj.m_G2);

		//------------------------------------------------------------------

		d_params.G=sim_obj.m_G;
		d_params.G1=sim_obj.m_G1;
		d_params.G2=sim_obj.m_G2;

		d_params.lambda = sim_obj.m_Lambda;
		d_params.lambda_error = sim_obj.m_LambdaError;

		d_params.K = sim_obj.m_K;
		d_params.K_error = sim_obj.m_KError;
			
		d_params.C = sim_obj.m_C;
		d_params.C_error = sim_obj.m_CError;

		d_params.sigma = sim_obj.m_Sigma;
		d_params.sigma_error = sim_obj.m_SigmaError;

		d_params.alpha_I = sim_obj.m_AlphaI;
		d_params.alpha_I_error = sim_obj.m_AlphaIError;

		d_params.alpha_J = sim_obj.m_AlphaJ;
		d_params.alpha_J_error = sim_obj.m_AlphaJError;

		d_params.a_I = sim_obj.m_AI;
		d_params.a_I_error = sim_obj.m_AIError;

		d_params.a_J = sim_obj.m_AJ;
		d_params.a_J_error = sim_obj.m_AJError;


		d_params.m_LambdaSbs=sim_obj.m_LambdaSbs;

		d_params.m_KSbs=sim_obj.m_KSbs;

		d_params.m_CSbs=sim_obj.m_CSbs;

		d_params.m_SigmaSbs=sim_obj.m_SigmaSbs;

		d_params.m_AlphaISbs=sim_obj.m_AlphaISbs;

		d_params.m_AlphaJSbs=sim_obj.m_AlphaJSbs;

		d_params.m_AISbs=sim_obj.m_AISbs;

		d_params.m_AJSbs=sim_obj.m_AJSbs;


		d_params.a 
			= (d_params.a_I + d_params.a_J) * 0.5;

		d_params.a_error = (d_params.a_I_error
						+ d_params.a_J_error)*0.5;
		
		d_params.alpha = (d_params.alpha_I
						+ d_params.alpha_J) * 0.5;

		d_params.alpha_error = (d_params.alpha_I_error
						+ d_params.alpha_J_error) * 0.5;

		d_params.d_params_flag=true;

		//precompute intercepts
		pvalues::compute_intercepts(d_params);

		double CurrentTime2;
		Sls::alp_data::get_current_time(CurrentTime2);
		d_params.m_CalcTime=CurrentTime2-CurrentTime1;
		delete randomization_parameters;randomization_parameters=NULL;

		if(!pvalues::assert_Gumbel_parameters(
		d_params)||!isGood())
		{
			d_params.d_params_flag=false;
			throw error("Error - computation of the Gumbel parameters is unsuccessful in the function \"void AlignmentEvaluer::initGapped\"\n",1);
		};

		delete[]letterFreqs1_normalized;
		delete[]letterFreqs2_normalized;


	}
	catch (...)
	{ 
		delete randomization_parameters;randomization_parameters=NULL;
		d_params.d_params_flag=false;
		throw;
	};
}

//Initializes Gumbel parameters using precalculated values:
void AlignmentEvaluer::initParameters(
const AlignmentEvaluerParameters &parameters_)
{
	try
	{
		double CurrentTime1;
		Sls::alp_data::get_current_time(CurrentTime1);

		d_params.d_params_flag=false;

		double calculation_error=1e-6;

		d_params.lambda=parameters_.d_lambda;
		d_params.lambda_error=calculation_error;

		d_params.C=0;
		d_params.C_error=0;


		d_params.K=parameters_.d_k;
		d_params.K_error=calculation_error;

		d_params.a_I=parameters_.d_a2;
		d_params.a_I_error=calculation_error;

		d_params.a_J=parameters_.d_a1;
		d_params.a_J_error=calculation_error;

		d_params.sigma=parameters_.d_sigma;
		d_params.sigma_error=calculation_error;

		d_params.alpha_I=parameters_.d_alpha2;
		d_params.alpha_I_error=calculation_error;

		d_params.alpha_J=parameters_.d_alpha1;
		d_params.alpha_J_error=calculation_error;

		d_params.a=0.5*(parameters_.d_a1+parameters_.d_a2);
		d_params.a_error=calculation_error;

		d_params.alpha=0.5*(parameters_.d_alpha1+parameters_.d_alpha2);
		d_params.alpha_error=calculation_error;

		d_params.gapless_a=0;
		d_params.gapless_a_error=0;

		d_params.gapless_alpha=0;
		d_params.gapless_alpha_error=0;

		d_params.G=0;
		d_params.G1=0;
		d_params.G2=0;

		//intercepts
		d_params.b_I=parameters_.d_b2;
		d_params.b_I_error=calculation_error;

		d_params.b_J=parameters_.d_b1;
		d_params.b_J_error=calculation_error;

		d_params.beta_I=parameters_.d_beta2;
		d_params.beta_I_error=calculation_error;

		d_params.beta_J=parameters_.d_beta1;
		d_params.beta_J_error=calculation_error;

		d_params.tau=parameters_.d_tau;
		d_params.tau_error=calculation_error;


		//arrays initialization
		std::vector<double > sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.lambda;
		sbs_arrays[1]=d_params.lambda + calculation_error;

		d_params.m_LambdaSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.K;
		sbs_arrays[1]=d_params.K+calculation_error;

		d_params.m_KSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.C;
		sbs_arrays[1]=d_params.C+calculation_error;

		d_params.m_CSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.sigma;
		sbs_arrays[1]=d_params.sigma + calculation_error;

		d_params.m_SigmaSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.alpha_I;
		sbs_arrays[1]=d_params.alpha_I + calculation_error;

		d_params.m_AlphaISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.alpha_J;
		sbs_arrays[1]=d_params.alpha_J + calculation_error;

		d_params.m_AlphaJSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.a_I;
		sbs_arrays[1]=d_params.a_I + calculation_error;

		d_params.m_AISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.a_J;
		sbs_arrays[1]=d_params.a_J + calculation_error;

		d_params.m_AJSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.b_J;
		sbs_arrays[1]=d_params.b_J + calculation_error;

		d_params.m_BJSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.b_I;
		sbs_arrays[1]=d_params.b_I + calculation_error;

		d_params.m_BISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.beta_J;
		sbs_arrays[1]=d_params.beta_J + calculation_error;

		d_params.m_BetaJSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.beta_I;
		sbs_arrays[1]=d_params.beta_I + calculation_error;

		d_params.m_BetaISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=d_params.tau;
		sbs_arrays[1]=d_params.tau + calculation_error;

		d_params.m_TauSbs=sbs_arrays;


		d_params.d_params_flag=true;//if true, then the parameters are defined and P-values can be calculated

		pvalues::compute_tmp_values(d_params);

		double CurrentTime2;
		Sls::alp_data::get_current_time(CurrentTime2);
		d_params.m_CalcTime=CurrentTime2-CurrentTime1;

		if(!pvalues::assert_Gumbel_parameters(
		d_params)||!isGood())
		{
			d_params.d_params_flag=false;
			throw error("Error - computation of the Gumbel parameters is unsuccessful in the function \"void AlignmentEvaluer::initParameters\"\n",1);
		};
	}
	catch (...)
	{ 
		d_params.d_params_flag=false;
		throw;
	};
	
}

void AlignmentEvaluer::initParameters(
const AlignmentEvaluerParametersWithErrors &parameters_)
{
	try
	{
		double CurrentTime1;
		Sls::alp_data::get_current_time(CurrentTime1);

		d_params.d_params_flag=false;

		long int array_dim=20;

		long int seed_tmp=12345;
		srand(seed_tmp);


		d_params.lambda=parameters_.d_lambda;
		d_params.lambda_error=parameters_.d_lambda_error;

		d_params.C=0;
		d_params.C_error=0;


		d_params.K=parameters_.d_k;
		d_params.K_error=parameters_.d_k_error;

		d_params.a_I=parameters_.d_a2;
		d_params.a_I_error=parameters_.d_a2_error;

		d_params.a_J=parameters_.d_a1;
		d_params.a_J_error=parameters_.d_a1_error;

		d_params.sigma=parameters_.d_sigma;
		d_params.sigma_error=parameters_.d_sigma_error;

		d_params.alpha_I=parameters_.d_alpha2;
		d_params.alpha_I_error=parameters_.d_alpha2_error;

		d_params.alpha_J=parameters_.d_alpha1;
		d_params.alpha_J_error=parameters_.d_alpha1_error;

		d_params.a=0.5*(parameters_.d_a1+parameters_.d_a2);
		d_params.a_error=0.5*(parameters_.d_a1_error+parameters_.d_a2_error);

		d_params.alpha=0.5*(parameters_.d_alpha1+parameters_.d_alpha2);
		d_params.alpha_error=0.5*(parameters_.d_alpha1_error+parameters_.d_alpha2_error);

		d_params.gapless_a=0;
		d_params.gapless_a_error=0;

		d_params.gapless_alpha=0;
		d_params.gapless_alpha_error=0;

		d_params.G=0;
		d_params.G1=0;
		d_params.G2=0;

		d_params.m_CalcTime=0;


		//intercepts
		d_params.b_I=parameters_.d_b2;
		d_params.b_I_error=parameters_.d_b2_error;

		d_params.b_J=parameters_.d_b1;
		d_params.b_J_error=parameters_.d_b1_error;

		d_params.beta_I=parameters_.d_beta2;
		d_params.beta_I_error=parameters_.d_beta2_error;

		d_params.beta_J=parameters_.d_beta1;
		d_params.beta_J_error=parameters_.d_beta1_error;

		d_params.tau=parameters_.d_tau;
		d_params.tau_error=parameters_.d_tau_error;

		double sqrt_array_dim=sqrt((double)array_dim);

		d_params.m_LambdaSbs.clear();
		d_params.m_KSbs.clear();
		d_params.m_CSbs.clear();

		d_params.m_SigmaSbs.clear();

		d_params.m_AlphaISbs.clear();
		d_params.m_AlphaJSbs.clear();

		d_params.m_AISbs.clear();
		d_params.m_AJSbs.clear();

		d_params.m_BISbs.clear();
		d_params.m_BJSbs.clear();

		d_params.m_BetaISbs.clear();
		d_params.m_BetaJSbs.clear();

		d_params.m_TauSbs.clear();

		long int i;
		for(i=0;i<array_dim;i++)
		{
			d_params.m_LambdaSbs.push_back(d_params.lambda+d_params.lambda_error*pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_KSbs.push_back(d_params.K+d_params.K_error*pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_CSbs.push_back(0);

			d_params.m_SigmaSbs.push_back(d_params.sigma+d_params.sigma_error*pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_AlphaISbs.push_back(d_params.alpha_I+d_params.alpha_I_error*pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_AlphaJSbs.push_back(d_params.alpha_J+d_params.alpha_J_error*pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_AISbs.push_back(d_params.a_I+d_params.a_I_error*pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_AJSbs.push_back(d_params.a_J+d_params.a_J_error*pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_BISbs.push_back(d_params.b_I+d_params.b_I_error*pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_BJSbs.push_back(d_params.b_J+d_params.b_J_error*pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_BetaISbs.push_back(d_params.beta_I+d_params.beta_I_error*pvalues::standard_normal()*sqrt_array_dim);
			d_params.m_BetaJSbs.push_back(d_params.beta_J+d_params.beta_J_error*pvalues::standard_normal()*sqrt_array_dim);

			d_params.m_TauSbs.push_back(d_params.tau+d_params.tau_error*pvalues::standard_normal()*sqrt_array_dim);
		};

		d_params.d_params_flag=true;//if true, then the parameters are defined and P-values can be calculated

		pvalues::compute_tmp_values(d_params);

		double CurrentTime2;
		Sls::alp_data::get_current_time(CurrentTime2);
		d_params.m_CalcTime=CurrentTime2-CurrentTime1;

		if(!pvalues::assert_Gumbel_parameters(
		d_params)||!isGood())
		{
			d_params.d_params_flag=false;
			throw error("Error - computation of the Gumbel parameters is unsuccessful in the function \"void AlignmentEvaluer::initParameters\"\n",1);
		};
	}
	catch (...)
	{ 
		d_params.d_params_flag=false;
		throw;
	};
	
}

double AlignmentEvaluer::area(double score_,//pairwise alignment score
double seqlen1_,//length of sequence #1
double seqlen2_) const//length of sequence #2
{
	if(seqlen1_<=0||seqlen2_<=0)
	{
		throw error("Error - seqlen1_<=0 or seq2en1_<=0 in \"double AlignmentEvaluer::area\"\n",2);
	};

	if(!isGood())
	{
		throw error("Unexpected error - the Gumbel parameters are not defined properly in \"double AlignmentEvaluer::area\"\n",1);
	};

	static Sls::pvalues pvalues_obj;

	double P;
	double E;
	double area_res;
	bool area_is_1_flag=false;
	bool compute_only_area=true;

	pvalues_obj.get_appr_tail_prob_with_cov_without_errors(
	d_params,
	pvalues_obj.blast,
	score_,
	seqlen2_,
	seqlen1_,

	P,

	E,

	area_res,
	area_is_1_flag,
	compute_only_area);


	return area_res;

}

void AlignmentEvaluer::calc(double score_,//pairwise alignment score
double seqlen1_,//length of sequence #1
double seqlen2_,//length of sequence #2
double &pvalue_,//resulted P-value
double &pvalueErr_,//P-value error
double &evalue_,//resulted E-value
double &evalueErr_) const//E-value error
{
	if(seqlen1_<=0||seqlen2_<=0)
	{
		throw error("Error - seqlen1_<=0 or seqlen2_<=0 in \"double AlignmentEvaluer::calc\"\n",2);
	};

	if(!isGood())
	{
		throw error("Unexpected error - the Gumbel parameters are not defined properly in \"double AlignmentEvaluer::calc\"\n",1);
	};

	static Sls::pvalues pvalues_obj;

	pvalues_obj.calculate_P_values(
		score_, seqlen2_, seqlen1_,
		d_params, 
		pvalue_,
		pvalueErr_,
		evalue_,
		evalueErr_);
}

void AlignmentEvaluer::calc(double score_,//pairwise alignment score
double seqlen1_,//length of sequence #1
double seqlen2_,//length of sequence #2
double &pvalue_,//resulted P-value
double &evalue_) const//resulted E-value
{

	if(seqlen1_<=0||seqlen2_<=0)
	{
		throw error("Error - seqlen1_<=0 or seqlen2_<=0 in \"double AlignmentEvaluer::calc\"\n",2);
	};

	if(!isGood())
	{
		throw error("Unexpected error - d_params is not defined in \"double AlignmentEvaluer::calc\"\n",1);
	};

	static Sls::pvalues pvalues_obj;


	bool area_is_1_flag=false;

	double area;


	pvalues_obj.get_appr_tail_prob_with_cov_without_errors(
	d_params,
	pvalues_obj.blast,
	score_,
	seqlen2_,
	seqlen1_,

	pvalue_,

	evalue_,

	area,
	area_is_1_flag);

	

}

void AlignmentEvaluer::set_gapped_computation_parameters_simplified(
double max_time_,//maximum allowed time in seconds for the whole computation
long number_of_samples_,//number of realization for the main stage
long number_of_samples_for_preliminary_stages_)//number of realization for the preliminary stages
{

	if(number_of_samples_<=0||number_of_samples_for_preliminary_stages_<=0)
	{
		throw error("Error - number_of_samples_<=0 or number_of_samples_for_preliminary_stages_<=0 in \"void AlignmentEvaluer::set_gapped_computation_parameters_simplified\"\n",2);
	};

	d_gapped_computation_parameters.d_first_stage_preliminary_realizations_numbers_ALP.resize(1);
	d_gapped_computation_parameters.d_first_stage_preliminary_realizations_numbers_ALP[0]=number_of_samples_for_preliminary_stages_;
	d_gapped_computation_parameters.d_preliminary_realizations_numbers_ALP.resize(1);
	d_gapped_computation_parameters.d_preliminary_realizations_numbers_ALP[0]=number_of_samples_for_preliminary_stages_;
	d_gapped_computation_parameters.d_preliminary_realizations_numbers_killing.resize(1);
	d_gapped_computation_parameters.d_preliminary_realizations_numbers_killing[0]=number_of_samples_for_preliminary_stages_;
	d_gapped_computation_parameters.d_total_realizations_number_with_ALP=number_of_samples_;
	d_gapped_computation_parameters.d_total_realizations_number_with_killing=number_of_samples_;
	d_gapped_computation_parameters.d_parameters_flag=true;
	d_gapped_computation_parameters.d_max_time_with_computation_parameters=max_time_;
	if(max_time_>0)
	{
		d_gapped_computation_parameters.d_max_time_for_quick_tests=0.5*max_time_*(double)quick_tests_trials_number/(double)(quick_tests_trials_number+number_of_samples_+number_of_samples_for_preliminary_stages_);
	}
	else
	{
		d_gapped_computation_parameters.d_max_time_for_quick_tests=-1;
	};
}

void AlignmentEvaluer::gapped_computation_parameters_clear()
{
	d_gapped_computation_parameters.d_first_stage_preliminary_realizations_numbers_ALP.clear();
	d_gapped_computation_parameters.d_preliminary_realizations_numbers_ALP.clear();
	d_gapped_computation_parameters.d_preliminary_realizations_numbers_killing.clear();
	d_gapped_computation_parameters.d_parameters_flag=false;
	d_gapped_computation_parameters.d_max_time_for_quick_tests=-1;
	d_gapped_computation_parameters.d_max_time_with_computation_parameters=-1;
}


}

