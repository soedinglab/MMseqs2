#ifndef INCLUDED_SLS_ALIGNMENT_EVALUER
#define INCLUDED_SLS_ALIGNMENT_EVALUER

/* $Id: $
* ===========================================================================
*
*							PUBLIC DOMAIN NOTICE
*			   National Center for Biotechnology Information
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

File name: sls_alignment_evaluer.hpp

Authors: Martin Frith, Sergey Sheetlin

Contents: library functions of main routines

******************************************************************************/

#include "sls_pvalues.hpp"
#include <math.h>

namespace Sls {
	const double default_importance_sampling_temperature = 1.07;

	class AlignmentEvaluer {

	//Write the parameters:
	friend std::ostream &operator<<(std::ostream &s_,
			const AlignmentEvaluer &g_);

	//Read the parameters:
	friend std::istream &operator>>(std::istream &s_,
			AlignmentEvaluer &g_);

	public:

	//randomization parameters
	struct gapped_computation_parameters_struct
	{
		gapped_computation_parameters_struct()
		{
			d_parameters_flag=false;
			d_max_time_for_quick_tests=-1;
			d_max_time_with_computation_parameters=-1;
		};

		std::vector<long int> d_first_stage_preliminary_realizations_numbers_ALP;
		std::vector<long int> d_preliminary_realizations_numbers_ALP;
		std::vector<long int> d_preliminary_realizations_numbers_killing;
		long int d_total_realizations_number_with_ALP;
		long int d_total_realizations_number_with_killing;
		bool d_parameters_flag;//if false then the parameters are not defined
		double d_max_time_for_quick_tests;//maximum allowed calculation time in seconds for quick tests
		double d_max_time_with_computation_parameters;//maximum allowed time in seconds for the whole computation
	};


	//Computes gapless Gumbel parameters:
	void initGapless(long alphabetSize_,//a number of letters in the alphabet
			const long *const *substitutionScoreMatrix_,//scoring matrix
			const double *letterFreqs1_,//background frequencies of letters in sequence #1
			const double *letterFreqs2_,//background frequencies of letters in sequence #2
			double max_time_=60);//maximum allowed calculation time in seconds

	//Computes gapped Gumbel parameters:
	//The NCBI convention is used for penalizing a gap:
	//For example, a gap of length k is penalized as gapOpen1_+k*gapEpen1_ for sequence #1
	//if max_time_<=0 then d_gapped_computation_parameters will be used;
	//d_gapped_computation_parameters.d_parameters_flag must be true in this case
	//every execution of initGapped with max_time_>0 rewrites d_gapped_computation_parameters by the actual values
	void initGapped(long alphabetSize_,//a number of letters in the alphabet
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
			double max_time_,//maximum allowed calculation time in seconds; 
			double max_mem_,//maximum allowed memory usage in Mb
			long randomSeed_,//randomizaton seed
			double temperature_=default_importance_sampling_temperature);


	//Initializes Gumbel parameters using precalculated values:
	void initParameters(
	const AlignmentEvaluerParameters &parameters_);//parameters_ must be defined by the user

	void initParameters(
	const AlignmentEvaluerParametersWithErrors &parameters_);//parameters_ must be defined by the user

	//Computes P-values/E-values
	//Gumbel parameters must be defined via d_params
	void calc(double score_,//pairwise alignment score
			double seqlen1_,//length of sequence #1
			double seqlen2_,//length of sequence #2
			double &pvalue_,//resulted P-value
			double &pvalueErr_,//P-value error
			double &evalue_,//resulted E-value
			double &evalueErr_) const;//E-value error

	//The following routines compute P-values/E-values and related
	//quantities quickly, without calculating errors

	//Gumbel parameters must be defined via d_params
	void calc(double score_,//pairwise alignment score
		double seqlen1_,//length of sequence #1
		double seqlen2_,//length of sequence #2
		double &pvalue_,//resulted P-value
		double &evalue_) const;//resulted E-value

	double evalue(double score_,//pairwise alignment score
			double seqlen1_,//length of sequence #1
			double seqlen2_) const//length of sequence #2
	{
		return area(score_, seqlen1_, seqlen2_) * evaluePerArea(score_);
	}

	//Computes P-values from E-values
	static double pvalue(double evalue_)
	{
		return sls_basic::one_minus_exp_function(-evalue_);
	}

	//The "area" is approximately seqlen1_*seqlen2_, but it is
	//modified by a finite size correction
	double area(double score_,//pairwise alignment score
			double seqlen1_,//length of sequence #1
			double seqlen2_) const;//length of sequence #2

	double evaluePerArea(double score_) const
	{
	  return d_params.K*exp(-d_params.lambda*score_);
	}

	double bitScore(double score_, double logK) const
	{
		return (d_params.lambda*score_-logK)/log(2.0);
	}

	double bitScore(double score_) const
	{
		return (d_params.lambda*score_-log(d_params.K))/log(2.0);
	}

	//returns "true" if the set of parameters "d_params" is fully defined for P-value calculation
	bool isGood() const
	{
		return d_params.d_params_flag;
	}

	//provides access to the set of Gumbel parameters
	const ALP_set_of_parameters &parameters() const { return d_params; }

	//provides access to the set of randomization parameters
	const gapped_computation_parameters_struct &gapped_computation_parameters() const { return d_gapped_computation_parameters; }

	//simplified setting of the computation parameters
	void set_gapped_computation_parameters_simplified(
		double max_time_=-1.0,//maximum allowed time in seconds for the whole computation
		long number_of_samples_=1000,//number of realization for the main stage
		long number_of_samples_for_preliminary_stages_=200);//number of realization for the preliminary stages

	void set_gapped_computation_parameters(
		const gapped_computation_parameters_struct &gapped_computation_parameters_)
	{
		d_gapped_computation_parameters=gapped_computation_parameters_;
	};

	//clear the computation parameters
	void gapped_computation_parameters_clear();
	
	private:
	//check correctness of the input parameters for gapless alignment
	void assert_Gapless_input_parameters(
		long alphabetSize_,//a number of letters in the alphabet
		const double *letterFreqs1_,//background frequencies of letters in sequence #1
		const double *letterFreqs2_,//background frequencies of letters in sequence #2
		double *&letterFreqs1_normalized_,//normalized background frequencies of letters in sequence #1
		double *&letterFreqs2_normalized_,//normalized background frequencies of letters in sequence #2
		const std::string function_name_);//"assert_Gapless_input_parameters" is called from "function_name_" function

	private:

		ALP_set_of_parameters d_params;//set of Gumbel parameters

		//generalized randomization parameters for the gapped computation
		gapped_computation_parameters_struct d_gapped_computation_parameters;

	};
}

#endif //! INCLUDED

