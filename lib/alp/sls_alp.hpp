#ifndef INCLUDED_SLS_ALP
#define INCLUDED_SLS_ALP

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

File name: sls_alp.hpp

Author: Sergey Sheetlin

Contents: Ascending ladder points simulation

******************************************************************************/

#include <complex>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <float.h>


#include "sls_alp_data.hpp"
#include "sls_alp_regression.hpp"


namespace Sls {


	class state//struct to save a state of calculation
	{
	public:
		state();

	public:
		array<long int> *d_cells_counts;

		long int *d_HS_i_const_next;
		long int *d_HI_i_const_next;
		long int *d_HD_i_const_next;
		long int *d_H_i_const_next;

		long int *d_HS_j_const_next;
		long int *d_HI_j_const_next;
		long int *d_HD_j_const_next;
		long int *d_H_j_const_next;

		long int d_HS_ij_next;
		long int d_HI_ij_next;
		long int d_HD_ij_next;
		long int d_H_ij_next;

		long int d_H_matr_len;

		long int d_M;

		long int d_sentinel_i_next;
		long int d_sentinel_j_next;

	};


	class alp{

	

	public:


	alp(//constructor
		alp_data *alp_data_
		);


	~alp();//destructor

	long int random_AA1();//generates random AA for the sequence 1
	long int random_AA2();//generates random AA for the sequence 2

	bool one_step_of_importance_sampling_without_weight_calculation(
		long int d_dim1_,
		long int d_dim2_);

	void increment_sequences();

	void increment_W_matrix();
	void increment_H_matrix();

	void increment_W_weights();
	//the function calculates weigths for d_W_matr_len increased by 1
	//assumes that letters are defined for d_W_matr_len

	void increment_H_weights();
	//the function calculates alignment scores for d_H_matr_len increased by 1
	//assumes that letters are defined for d_H_matr_len

	void increment_H_weights_without_insertions_after_deletions();

	void increment_H_weights_with_insertions_after_deletions();


	void increment_H_weights_with_sentinels(
		long int diff_opt_);
	//the function calculates alignment scores for d_H_matr_len increased by 1
	//assumes that letters are defined for d_H_matr_len
	//uses sentinels

	void increment_H_weights_with_sentinels_without_insertions_after_deletions(
	long int diff_opt_);

	void increment_H_weights_with_sentinels_with_insertions_after_deletions(
	long int diff_opt_);



	void simulate_next_alp();//simulates next ALP

	void simulate_alp_upto_the_given_number(//simulates ALP upto the given number nalp_ including
	long int nalp_);

	void simulate_alp_upto_the_given_level(//simulates ALP upto the given level M_min_ including
	long int M_min_);





	template<typename T>
	inline void swap(
	T& a1_,
	T& a2_)
	{
		T tmp=a1_;
		a1_=a2_;
		a2_=tmp;
	}

	static double degree(//returns x_^n_
		double x_,
		double n_);


	double John2_weight_calculation(
		long int length_);//calculation of weigths for the importance sampling

	void save_state(
		state * &state_);

	void restore_state(
		state * &state_);

	void kill_upto_level(
		long int M_min_,
		long int M_level_,
		long int *M_upper_level_=NULL);

	void check_time_function();

	template<typename T>
	void release_and_calculate_memory(
	T *&pointer_,
	long int dim_)
	{
		if(pointer_==NULL)
		{
			return;
		};
		delete[]pointer_;pointer_=NULL;
		if(d_alp_data)
		{
			d_alp_data->d_memory_size_in_MB-=(double)(sizeof(T)*dim_)/mb_bytes;
		};
		
	}


	template<typename T>
	void release_and_calculate_memory(
	T *&pointer_)
	{
		if(pointer_==NULL)
		{
			return;
		};
		delete pointer_;pointer_=NULL;
		if(d_alp_data)
		{
			d_alp_data->d_memory_size_in_MB-=(double)(sizeof(T))/mb_bytes;
		};
		
	}

	void partially_release_memory();







	public:


				
	alp_data *d_alp_data;//initial data
	long int d_a_step;//increment for sequence length during memory allocation


	bool d_is_now;//true if the importance sampling is being used 

	



	//alignment data

	long int d_seqi_len;//current length of sequence 1
	long int d_seqj_len;//current length of sequence 2


	long int d_seq_a_len;//current length for memory allocation for the sequences
	long int d_H_matr_a_len;//current length for memory allocation for the matrices H
	long int d_W_matr_a_len;//current length for memory allocation for the matrices W

	long int *d_seqi;//AA from sequence 1 for the next step
	long int *d_seqj;//AA from sequence 2 for the next step

	long int d_H_matr_len;//length of the matrices H currently calculated
	long int d_W_matr_len;//length of the matrices W currently calculated

	//the importance sampling weights
	double *d_WS_i_const_pred;
	double *d_WI_i_const_pred;
	double *d_WD_i_const_pred;

	double *d_WS_i_const_next;
	double *d_WI_i_const_next;
	double *d_WD_i_const_next;

	double *d_WS_j_const_pred;
	double *d_WI_j_const_pred;
	double *d_WD_j_const_pred;

	double *d_WS_j_const_next;
	double *d_WI_j_const_next;
	double *d_WD_j_const_next;

	double d_WS_ij_pred;
	double d_WI_ij_pred;
	double d_WD_ij_pred;

	double d_WS_ij_next;
	double d_WI_ij_next;
	double d_WD_ij_next;


	//alignment matrix 
	long int *d_HS_i_const_pred;
	long int *d_HI_i_const_pred;
	long int *d_HD_i_const_pred;
	long int *d_H_i_const_pred;

	long int *d_HS_i_const_next;
	long int *d_HI_i_const_next;
	long int *d_HD_i_const_next;
	long int *d_H_i_const_next;

	long int *d_HS_j_const_pred;
	long int *d_HI_j_const_pred;
	long int *d_HD_j_const_pred;
	long int *d_H_j_const_pred;

	long int *d_HS_j_const_next;
	long int *d_HI_j_const_next;
	long int *d_HD_j_const_next;
	long int *d_H_j_const_next;

	long int d_HS_ij_pred;
	long int d_HI_ij_pred;
	long int d_HD_ij_pred;
	long int d_H_ij_pred;

	long int d_HS_ij_next;
	long int d_HI_ij_next;
	long int d_HD_ij_next;
	long int d_H_ij_next;

	bool d_success;

	//statistics
	long int *d_H_edge_max;
	long int d_M;

	long int d_nalp;
	long int d_nalp_killing;
	array_positive<long int> *d_alp;

	array_positive<long int> *d_H_I;
	array_positive<long int> *d_H_J;

	array_positive<long int> *d_alp_pos;
	array_positive<double> *d_alp_weights;

	array<long int> *d_cells_counts;

	array_positive<state*> *d_alp_states;

	long int d_sentinel_i_next;
	long int d_sentinel_j_next;

	long int d_sentinel_i_pred;
	long int d_sentinel_j_pred;

	long int d_diff_opt;
	bool d_sentinels_flag;

	bool d_check_time_flag;
	bool d_time_error_flag;
	bool d_time_limit_flag;

	bool d_single_realiztion_calculation_flag;



	//for the importance sampling
	char d_IS_state;




				


	};

}


#endif

