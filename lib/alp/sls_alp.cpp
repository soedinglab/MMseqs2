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

File name: sls_alp.cpp

Author: Sergey Sheetlin, Martin Frith

Contents: Ascending ladder points simulation

******************************************************************************/


#include "sls_alp.hpp"

using namespace Sls;

static long int small_long=(long int)((double)LONG_MIN/2.0);


alp::alp(//constructor
alp_data *alp_data_
)
{
	d_seqi=NULL;
	d_seqj=NULL;

	d_WS_i_const_pred=NULL;
	d_WI_i_const_pred=NULL;
	d_WD_i_const_pred=NULL;

	d_WS_i_const_next=NULL;
	d_WI_i_const_next=NULL;
	d_WD_i_const_next=NULL;

	d_WS_j_const_pred=NULL;
	d_WI_j_const_pred=NULL;
	d_WD_j_const_pred=NULL;

	d_WS_j_const_next=NULL;
	d_WI_j_const_next=NULL;
	d_WD_j_const_next=NULL;


	//alignment matrix 
	d_HS_i_const_pred=NULL;
	d_HI_i_const_pred=NULL;
	d_HD_i_const_pred=NULL;
	d_H_i_const_pred=NULL;

	d_HS_i_const_next=NULL;
	d_HI_i_const_next=NULL;
	d_HD_i_const_next=NULL;
	d_H_i_const_next=NULL;

	d_HS_j_const_pred=NULL;
	d_HI_j_const_pred=NULL;
	d_HD_j_const_pred=NULL;
	d_H_j_const_pred=NULL;

	d_HS_j_const_next=NULL;
	d_HI_j_const_next=NULL;
	d_HD_j_const_next=NULL;
	d_H_j_const_next=NULL;

	d_H_edge_max=NULL;
	d_H_I=NULL;
	d_H_J=NULL;

	d_alp=NULL;
	d_alp_pos=NULL;
	d_cells_counts=NULL;
	d_alp_weights=NULL;
	d_alp_states=NULL;

	d_success=true;

	d_check_time_flag=false;
	d_time_error_flag=false;
	d_time_limit_flag=false;
	d_single_realiztion_calculation_flag=false;

	d_alp_data=alp_data_;
	if(!d_alp_data)
	{
		throw error("Unexpected error\n",4);
	};

	d_a_step=30;

	try
	{

		d_is_now=true;
		d_seqi_len=0;
		d_seqj_len=0;
		d_seq_a_len=0;
		d_H_matr_a_len=0;
		d_W_matr_a_len=0;
		d_H_matr_len=-1;
		d_W_matr_len=-1;
		d_nalp=-1;

		d_H_edge_max=new long int[1];
		alp_data::assert_mem(d_H_edge_max);
		d_H_edge_max[0]=0;

		d_alp_data->d_memory_size_in_MB+=(double)(sizeof(long int))/mb_bytes;


		d_alp=new array_positive<long int>(d_alp_data);
		alp_data::assert_mem(d_alp);

		d_H_I=new array_positive<long int>(d_alp_data);
		alp_data::assert_mem(d_H_I);

		d_H_J=new array_positive<long int>(d_alp_data);
		alp_data::assert_mem(d_H_J);


		d_alp_pos=new array_positive<long int>(d_alp_data);
		alp_data::assert_mem(d_alp_pos);

		d_alp_data->d_memory_size_in_MB+=4*(double)(sizeof(array_positive<long int>))/mb_bytes;

		d_alp_states=new array_positive<state*>(d_alp_data);
		alp_data::assert_mem(d_alp_states);

		d_alp_data->d_memory_size_in_MB+=(double)(sizeof(array_positive<state*>))/mb_bytes;


		d_alp_weights=new array_positive<double>(d_alp_data);
		alp_data::assert_mem(d_alp_weights);

		d_alp_data->d_memory_size_in_MB+=(double)(sizeof(array_positive<double>))/mb_bytes;


		d_cells_counts=new array<long int>(d_alp_data);
		alp_data::assert_mem(d_cells_counts);

		d_alp_data->d_memory_size_in_MB+=(double)(sizeof(array<long int>))/mb_bytes;


		increment_W_weights();
		increment_H_weights_with_sentinels(0);
	}
	catch (...)
	{
		this->~alp();
		throw;
	};

}

alp::~alp()//destructor
{

	release_and_calculate_memory(d_seqi,d_seq_a_len);
	release_and_calculate_memory(d_seqj,d_seq_a_len);



	release_and_calculate_memory(d_WS_i_const_pred,d_W_matr_a_len);

	release_and_calculate_memory(d_WI_i_const_pred,d_W_matr_a_len);

	release_and_calculate_memory(d_WD_i_const_pred,d_W_matr_a_len);


	release_and_calculate_memory(d_WS_i_const_next,d_W_matr_a_len);

	release_and_calculate_memory(d_WI_i_const_next,d_W_matr_a_len);

	release_and_calculate_memory(d_WD_i_const_next,d_W_matr_a_len);



	release_and_calculate_memory(d_WS_j_const_pred,d_W_matr_a_len);

	release_and_calculate_memory(d_WI_j_const_pred,d_W_matr_a_len);

	release_and_calculate_memory(d_WD_j_const_pred,d_W_matr_a_len);


	release_and_calculate_memory(d_WS_j_const_next,d_W_matr_a_len);

	release_and_calculate_memory(d_WI_j_const_next,d_W_matr_a_len);

	release_and_calculate_memory(d_WD_j_const_next,d_W_matr_a_len);



	release_and_calculate_memory(d_HS_i_const_pred,d_H_matr_a_len);

	release_and_calculate_memory(d_HI_i_const_pred,d_H_matr_a_len);

	release_and_calculate_memory(d_HD_i_const_pred,d_H_matr_a_len);

	release_and_calculate_memory(d_H_i_const_pred,d_H_matr_a_len);


	release_and_calculate_memory(d_HS_i_const_next,d_H_matr_a_len);

	release_and_calculate_memory(d_HI_i_const_next,d_H_matr_a_len);

	release_and_calculate_memory(d_HD_i_const_next,d_H_matr_a_len);

	release_and_calculate_memory(d_H_i_const_next,d_H_matr_a_len);


	release_and_calculate_memory(d_HS_j_const_pred,d_H_matr_a_len);

	release_and_calculate_memory(d_HI_j_const_pred,d_H_matr_a_len);

	release_and_calculate_memory(d_HD_j_const_pred,d_H_matr_a_len);

	release_and_calculate_memory(d_H_j_const_pred,d_H_matr_a_len);


	release_and_calculate_memory(d_HS_j_const_next,d_H_matr_a_len);

	release_and_calculate_memory(d_HI_j_const_next,d_H_matr_a_len);

	release_and_calculate_memory(d_HD_j_const_next,d_H_matr_a_len);

	release_and_calculate_memory(d_H_j_const_next,d_H_matr_a_len);


	release_and_calculate_memory(d_H_edge_max,d_H_matr_a_len+1);



	
	release_and_calculate_memory(d_alp);

	release_and_calculate_memory(d_H_I);

	release_and_calculate_memory(d_H_J);

	release_and_calculate_memory(d_alp_pos);





	long int i;

	if(d_alp_states)
	{
		for(i=0;i<=d_nalp;i++)
		{
			if(i<=d_alp_states->d_dim)
			{
				if(d_alp_states->d_elem[i])
				{

					for(i=0;i<=d_nalp;i++)
					{


						release_and_calculate_memory(d_alp_states->d_elem[i]->d_HS_i_const_next,d_alp_states->d_elem[i]->d_H_matr_len);

						release_and_calculate_memory(d_alp_states->d_elem[i]->d_HI_i_const_next,d_alp_states->d_elem[i]->d_H_matr_len);

						release_and_calculate_memory(d_alp_states->d_elem[i]->d_HD_i_const_next,d_alp_states->d_elem[i]->d_H_matr_len);

						release_and_calculate_memory(d_alp_states->d_elem[i]->d_H_i_const_next,d_alp_states->d_elem[i]->d_H_matr_len);


						release_and_calculate_memory(d_alp_states->d_elem[i]->d_HS_j_const_next,d_alp_states->d_elem[i]->d_H_matr_len);

						release_and_calculate_memory(d_alp_states->d_elem[i]->d_HI_j_const_next,d_alp_states->d_elem[i]->d_H_matr_len);

						release_and_calculate_memory(d_alp_states->d_elem[i]->d_HD_j_const_next,d_alp_states->d_elem[i]->d_H_matr_len);

						release_and_calculate_memory(d_alp_states->d_elem[i]->d_H_j_const_next,d_alp_states->d_elem[i]->d_H_matr_len);





						release_and_calculate_memory(d_alp_states->d_elem[i]->d_cells_counts);

						release_and_calculate_memory(d_alp_states->d_elem[i]);

						
					};
				};
			};
		};
	};

	
	release_and_calculate_memory(d_alp_states);

	release_and_calculate_memory(d_alp_weights);

	
	release_and_calculate_memory(d_cells_counts);


}

void alp::partially_release_memory()
{

	
	release_and_calculate_memory(d_seqi,d_seq_a_len);
	release_and_calculate_memory(d_seqj,d_seq_a_len);


	release_and_calculate_memory(d_WS_i_const_pred,d_W_matr_a_len);
	release_and_calculate_memory(d_WI_i_const_pred,d_W_matr_a_len);
	release_and_calculate_memory(d_WD_i_const_pred,d_W_matr_a_len);

	release_and_calculate_memory(d_WS_i_const_next,d_W_matr_a_len);
	release_and_calculate_memory(d_WI_i_const_next,d_W_matr_a_len);
	release_and_calculate_memory(d_WD_i_const_next,d_W_matr_a_len);


	release_and_calculate_memory(d_WS_j_const_pred,d_W_matr_a_len);
	release_and_calculate_memory(d_WI_j_const_pred,d_W_matr_a_len);
	release_and_calculate_memory(d_WD_j_const_pred,d_W_matr_a_len);


	release_and_calculate_memory(d_WS_j_const_next,d_W_matr_a_len);
	release_and_calculate_memory(d_WI_j_const_next,d_W_matr_a_len);
	release_and_calculate_memory(d_WD_j_const_next,d_W_matr_a_len);



	release_and_calculate_memory(d_HS_i_const_pred,d_H_matr_a_len);
	release_and_calculate_memory(d_HI_i_const_pred,d_H_matr_a_len);
	release_and_calculate_memory(d_HD_i_const_pred,d_H_matr_a_len);
	release_and_calculate_memory(d_H_i_const_pred,d_H_matr_a_len);


	release_and_calculate_memory(d_HS_i_const_next,d_H_matr_a_len);
	release_and_calculate_memory(d_HI_i_const_next,d_H_matr_a_len);
	release_and_calculate_memory(d_HD_i_const_next,d_H_matr_a_len);
	release_and_calculate_memory(d_H_i_const_next,d_H_matr_a_len);



	release_and_calculate_memory(d_HS_j_const_pred,d_H_matr_a_len);
	release_and_calculate_memory(d_HI_j_const_pred,d_H_matr_a_len);
	release_and_calculate_memory(d_HD_j_const_pred,d_H_matr_a_len);
	release_and_calculate_memory(d_H_j_const_pred,d_H_matr_a_len);


	release_and_calculate_memory(d_HS_j_const_next,d_H_matr_a_len);
	release_and_calculate_memory(d_HI_j_const_next,d_H_matr_a_len);
	release_and_calculate_memory(d_HD_j_const_next,d_H_matr_a_len);
	release_and_calculate_memory(d_H_j_const_next,d_H_matr_a_len);



	release_and_calculate_memory(d_H_edge_max,d_H_matr_a_len+1);


	

	long int i;

	if(d_alp_states)
	{
		for(i=0;i<=d_nalp;i++)
		{
			if(i<=d_alp_states->d_dim)
			{
				if(d_alp_states->d_elem[i])
				{

					release_and_calculate_memory(d_alp_states->d_elem[i]->d_HS_i_const_next,d_alp_states->d_elem[i]->d_H_matr_len);
					release_and_calculate_memory(d_alp_states->d_elem[i]->d_HI_i_const_next,d_alp_states->d_elem[i]->d_H_matr_len);
					release_and_calculate_memory(d_alp_states->d_elem[i]->d_HD_i_const_next,d_alp_states->d_elem[i]->d_H_matr_len);
					release_and_calculate_memory(d_alp_states->d_elem[i]->d_H_i_const_next,d_alp_states->d_elem[i]->d_H_matr_len);

					release_and_calculate_memory(d_alp_states->d_elem[i]->d_HS_j_const_next,d_alp_states->d_elem[i]->d_H_matr_len);
					release_and_calculate_memory(d_alp_states->d_elem[i]->d_HI_j_const_next,d_alp_states->d_elem[i]->d_H_matr_len);
					release_and_calculate_memory(d_alp_states->d_elem[i]->d_HD_j_const_next,d_alp_states->d_elem[i]->d_H_matr_len);
					release_and_calculate_memory(d_alp_states->d_elem[i]->d_H_j_const_next,d_alp_states->d_elem[i]->d_H_matr_len);



					release_and_calculate_memory(d_alp_states->d_elem[i]->d_cells_counts);
				};
			};

			
		};
	};

}


long int alp::random_AA1()
{
	return d_alp_data->random_long(
		d_alp_data->ran2(),
		d_alp_data->d_number_of_AA,
		d_alp_data->d_RR1_sum,
		d_alp_data->d_RR1_sum_elements);
}

long int alp::random_AA2()
{
	return d_alp_data->random_long(
		d_alp_data->ran2(),
		d_alp_data->d_number_of_AA,
		d_alp_data->d_RR2_sum,
		d_alp_data->d_RR2_sum_elements);
}

bool alp::one_step_of_importance_sampling_without_weight_calculation(
	long int d_dim1_,
	long int d_dim2_)
{
	char &state_=d_IS_state;

	alp_data *&las_object_=d_alp_data;

	importance_sampling *&d_is_=d_alp_data->d_is;
	long int &length1_=d_seqi_len;
	long int &length2_=d_seqj_len;

	long int *&d_seqi_rglobal_=d_seqi;
	long int *&d_seqj_rglobal_=d_seqj;


	bool res=true;

	if(length1_==0&&length2_==0)
	{
		state_=alp_data::random_long(
			las_object_->ran2(),
			3,
			d_is_->d_for_S,
			d_is_->d_for_S_states);
	};


	if(state_=='D')
	{
		if(length1_==d_dim1_)
		{
			res=false;
			return res;
		};

		if(length1_>d_seq_a_len-1)
		{
			increment_sequences();
		};

		d_seqi_rglobal_[length1_]=random_AA1();
		length1_++;

		state_=alp_data::random_long(
			las_object_->ran2(),
			3,
			d_is_->d_for_D,
			d_is_->d_for_D_states);
		goto weight_calculation;
	};

	if(state_=='I')
	{
		if(length2_==d_dim2_)
		{
			res=false;
			return res;
		};

		if(length2_>d_seq_a_len-1)
		{
			increment_sequences();
		};


		d_seqj_rglobal_[length2_]=random_AA2();
		length2_++;

		state_=alp_data::random_long(
			las_object_->ran2(),
			2,
			d_is_->d_for_I,
			d_is_->d_for_I_states);
		goto weight_calculation;
	};

	if(state_=='S')
	{
		if(length1_==d_dim1_||length2_==d_dim2_)
		{
			res=false;
			return res;
		};
		q_elem pair=alp_data::random_long(
			las_object_->ran2(),
			d_is_->d_is_number_of_AA*d_is_->d_is_number_of_AA,
			d_is_->d_elements_values,
			d_is_->d_elements);

		if(length1_>d_seq_a_len-1||length2_>d_seq_a_len-1)
		{
			increment_sequences();
		};

		d_seqi_rglobal_[length1_]=pair.d_a;
		d_seqj_rglobal_[length2_]=pair.d_b;

		

		length1_++;
		length2_++;

		state_=alp_data::random_long(
			las_object_->ran2(),
			3,
			d_is_->d_for_S,
			d_is_->d_for_S_states);
		goto weight_calculation;
	};

weight_calculation:
//----deleted-------
return res;

}

void alp::increment_sequences()
{
	long int *d_seqi_new=NULL;
	long int *d_seqj_new=NULL;

	try
	{

		d_seq_a_len+=d_a_step;

		d_seqi_new=new long int[d_seq_a_len];
		alp_data::assert_mem(d_seqi_new);

		d_seqj_new=new long int[d_seq_a_len];
		alp_data::assert_mem(d_seqj_new);

		long int i;
		for(i=0;i<d_seqi_len;i++)
		{
			d_seqi_new[i]=d_seqi[i];
		};

		for(i=0;i<d_seqj_len;i++)
		{
			d_seqj_new[i]=d_seqj[i];
		};

		
		delete[]d_seqi;d_seqi=NULL;
		delete[]d_seqj;d_seqj=NULL;


		d_seqi=d_seqi_new;
		d_seqj=d_seqj_new;

		d_seqi_new=NULL;
		d_seqj_new=NULL;

		d_alp_data->d_memory_size_in_MB+=(double)(sizeof(long int)*d_a_step*2)/mb_bytes;

	}
	catch (...)
	{ 
		delete[]d_seqi_new;d_seqi_new=NULL;
		delete[]d_seqj_new;d_seqj_new=NULL;
		throw;
	};

}

void alp::increment_W_matrix()
{
	double *d_WS_i_const_pred_new=NULL;
	double *d_WI_i_const_pred_new=NULL;
	double *d_WD_i_const_pred_new=NULL;

	double *d_WS_i_const_next_new=NULL;
	double *d_WI_i_const_next_new=NULL;
	double *d_WD_i_const_next_new=NULL;

	double *d_WS_j_const_pred_new=NULL;
	double *d_WI_j_const_pred_new=NULL;
	double *d_WD_j_const_pred_new=NULL;

	double *d_WS_j_const_next_new=NULL;
	double *d_WI_j_const_next_new=NULL;
	double *d_WD_j_const_next_new=NULL;

	try
	{

		d_W_matr_a_len+=d_a_step;

		//the importance sampling weights
		d_WS_i_const_pred_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WS_i_const_pred_new);
		d_WI_i_const_pred_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WI_i_const_pred_new);
		d_WD_i_const_pred_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WD_i_const_pred_new);

		d_WS_i_const_next_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WS_i_const_next_new);
		d_WI_i_const_next_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WI_i_const_next_new);
		d_WD_i_const_next_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WD_i_const_next_new);

		d_WS_j_const_pred_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WS_j_const_pred_new);
		d_WI_j_const_pred_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WI_j_const_pred_new);
		d_WD_j_const_pred_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WD_j_const_pred_new);

		d_WS_j_const_next_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WS_j_const_next_new);
		d_WI_j_const_next_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WI_j_const_next_new);
		d_WD_j_const_next_new=new double[d_W_matr_a_len];
		alp_data::assert_mem(d_WD_j_const_next_new);



		long int i;
		for(i=0;i<d_W_matr_len;i++)
		{
			d_WS_i_const_next_new[i]=d_WS_i_const_next[i];
			d_WI_i_const_next_new[i]=d_WI_i_const_next[i];
			d_WD_i_const_next_new[i]=d_WD_i_const_next[i];

			d_WS_j_const_next_new[i]=d_WS_j_const_next[i];
			d_WI_j_const_next_new[i]=d_WI_j_const_next[i];
			d_WD_j_const_next_new[i]=d_WD_j_const_next[i];

		};

		

		for(i=0;i<d_W_matr_len-1;i++)
		{
			d_WS_i_const_pred_new[i]=d_WS_i_const_pred[i];
			d_WI_i_const_pred_new[i]=d_WI_i_const_pred[i];
			d_WD_i_const_pred_new[i]=d_WD_i_const_pred[i];

			d_WS_j_const_pred_new[i]=d_WS_j_const_pred[i];
			d_WI_j_const_pred_new[i]=d_WI_j_const_pred[i];
			d_WD_j_const_pred_new[i]=d_WD_j_const_pred[i];

		};


		delete[]d_WS_i_const_pred;d_WS_i_const_pred=NULL;
		delete[]d_WI_i_const_pred;d_WI_i_const_pred=NULL;
		delete[]d_WD_i_const_pred;d_WD_i_const_pred=NULL;

		delete[]d_WS_i_const_next;d_WS_i_const_next=NULL;
		delete[]d_WI_i_const_next;d_WI_i_const_next=NULL;
		delete[]d_WD_i_const_next;d_WD_i_const_next=NULL;

		delete[]d_WS_j_const_pred;d_WS_j_const_pred=NULL;
		delete[]d_WI_j_const_pred;d_WI_j_const_pred=NULL;
		delete[]d_WD_j_const_pred;d_WD_j_const_pred=NULL;

		delete[]d_WS_j_const_next;d_WS_j_const_next=NULL;
		delete[]d_WI_j_const_next;d_WI_j_const_next=NULL;
		delete[]d_WD_j_const_next;d_WD_j_const_next=NULL;

		d_alp_data->d_memory_size_in_MB+=(double)(sizeof(double)*d_a_step*12)/mb_bytes;


		d_WS_i_const_pred=d_WS_i_const_pred_new;d_WS_i_const_pred_new=NULL;
		d_WI_i_const_pred=d_WI_i_const_pred_new;d_WI_i_const_pred_new=NULL;
		d_WD_i_const_pred=d_WD_i_const_pred_new;d_WD_i_const_pred_new=NULL;

		d_WS_i_const_next=d_WS_i_const_next_new;d_WS_i_const_next_new=NULL;
		d_WI_i_const_next=d_WI_i_const_next_new;d_WI_i_const_next_new=NULL;
		d_WD_i_const_next=d_WD_i_const_next_new;d_WD_i_const_next_new=NULL;

		d_WS_j_const_pred=d_WS_j_const_pred_new;d_WS_j_const_pred_new=NULL;
		d_WI_j_const_pred=d_WI_j_const_pred_new;d_WI_j_const_pred_new=NULL;
		d_WD_j_const_pred=d_WD_j_const_pred_new;d_WD_j_const_pred_new=NULL;


		d_WS_j_const_next=d_WS_j_const_next_new;d_WS_j_const_next_new=NULL;
		d_WI_j_const_next=d_WI_j_const_next_new;d_WI_j_const_next_new=NULL;
		d_WD_j_const_next=d_WD_j_const_next_new;d_WD_j_const_next_new=NULL;

		

	}
	catch (...)
	{ 
		delete[]d_WS_i_const_pred_new;d_WS_i_const_pred_new=NULL;
		delete[]d_WI_i_const_pred_new;d_WI_i_const_pred_new=NULL;
		delete[]d_WD_i_const_pred_new;d_WD_i_const_pred_new=NULL;

		delete[]d_WS_i_const_next_new;d_WS_i_const_next_new=NULL;
		delete[]d_WI_i_const_next_new;d_WI_i_const_next_new=NULL;
		delete[]d_WD_i_const_next_new;d_WD_i_const_next_new=NULL;

		delete[]d_WS_j_const_pred_new;d_WS_j_const_pred_new=NULL;
		delete[]d_WI_j_const_pred_new;d_WI_j_const_pred_new=NULL;
		delete[]d_WD_j_const_pred_new;d_WD_j_const_pred_new=NULL;

		delete[]d_WS_j_const_next_new;d_WS_j_const_next_new=NULL;
		delete[]d_WI_j_const_next_new;d_WI_j_const_next_new=NULL;
		delete[]d_WD_j_const_next_new;d_WD_j_const_next_new=NULL;
		throw;
	};

}

void alp::increment_H_matrix()
{

	long int *d_HS_i_const_pred_new=NULL;
	long int *d_HI_i_const_pred_new=NULL;
	long int *d_HD_i_const_pred_new=NULL;
	long int *d_H_i_const_pred_new=NULL;

	long int *d_HS_i_const_next_new=NULL;
	long int *d_HI_i_const_next_new=NULL;
	long int *d_HD_i_const_next_new=NULL;
	long int *d_H_i_const_next_new=NULL;

	long int *d_HS_j_const_pred_new=NULL;
	long int *d_HI_j_const_pred_new=NULL;
	long int *d_HD_j_const_pred_new=NULL;
	long int *d_H_j_const_pred_new=NULL;

	long int *d_HS_j_const_next_new=NULL;
	long int *d_HI_j_const_next_new=NULL;
	long int *d_HD_j_const_next_new=NULL;
	long int *d_H_j_const_next_new=NULL;

	long int *d_H_edge_max_new=NULL;

	try
	{

		d_H_matr_a_len+=d_a_step;


		//alignment matrix 
		d_HS_i_const_pred_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HS_i_const_pred_new);
		d_HI_i_const_pred_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HI_i_const_pred_new);
		d_HD_i_const_pred_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HD_i_const_pred_new);
		d_H_i_const_pred_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_H_i_const_pred_new);

		d_HS_i_const_next_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HS_i_const_next_new);
		d_HI_i_const_next_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HI_i_const_next_new);
		d_HD_i_const_next_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HD_i_const_next_new);
		d_H_i_const_next_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_H_i_const_next_new);

		d_HS_j_const_pred_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HS_j_const_pred_new);
		d_HI_j_const_pred_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HI_j_const_pred_new);
		d_HD_j_const_pred_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HD_j_const_pred_new);
		d_H_j_const_pred_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_H_j_const_pred_new);

		d_HS_j_const_next_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HS_j_const_next_new);
		d_HI_j_const_next_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HI_j_const_next_new);
		d_HD_j_const_next_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_HD_j_const_next_new);
		d_H_j_const_next_new=new long int[d_H_matr_a_len];
		alp_data::assert_mem(d_H_j_const_next_new);


		d_H_edge_max_new=new long int[d_H_matr_a_len+1];
		alp_data::assert_mem(d_H_edge_max_new);


		long int i;
		for(i=0;i<d_H_matr_len;i++)
		{
			d_HS_i_const_next_new[i]=d_HS_i_const_next[i];
			d_HI_i_const_next_new[i]=d_HI_i_const_next[i];
			d_HD_i_const_next_new[i]=d_HD_i_const_next[i];
			d_H_i_const_next_new[i]=d_H_i_const_next[i];

			d_HS_j_const_next_new[i]=d_HS_j_const_next[i];
			d_HI_j_const_next_new[i]=d_HI_j_const_next[i];
			d_HD_j_const_next_new[i]=d_HD_j_const_next[i];
			d_H_j_const_next_new[i]=d_H_j_const_next[i];
		};

		for(i=0;i<d_H_matr_len-1;i++)
		{
			d_HS_i_const_pred_new[i]=d_HS_i_const_pred[i];
			d_HI_i_const_pred_new[i]=d_HI_i_const_pred[i];
			d_HD_i_const_pred_new[i]=d_HD_i_const_pred[i];
			d_H_i_const_pred_new[i]=d_H_i_const_pred[i];

			d_HS_j_const_pred_new[i]=d_HS_j_const_pred[i];
			d_HI_j_const_pred_new[i]=d_HI_j_const_pred[i];
			d_HD_j_const_pred_new[i]=d_HD_j_const_pred[i];
			d_H_j_const_pred_new[i]=d_H_j_const_pred[i];
		};


		for(i=0;i<=d_H_matr_len;i++)
		{
			d_H_edge_max_new[i]=d_H_edge_max[i];
		};

		

		delete[]d_HS_i_const_pred;d_HS_i_const_pred=NULL;
		delete[]d_HI_i_const_pred;d_HI_i_const_pred=NULL;
		delete[]d_HD_i_const_pred;d_HD_i_const_pred=NULL;
		delete[]d_H_i_const_pred;d_H_i_const_pred=NULL;

		delete[]d_HS_i_const_next;d_HS_i_const_next=NULL;
		delete[]d_HI_i_const_next;d_HI_i_const_next=NULL;
		delete[]d_HD_i_const_next;d_HD_i_const_next=NULL;
		delete[]d_H_i_const_next;d_H_i_const_next=NULL;

		delete[]d_HS_j_const_pred;d_HS_j_const_pred=NULL;
		delete[]d_HI_j_const_pred;d_HI_j_const_pred=NULL;
		delete[]d_HD_j_const_pred;d_HD_j_const_pred=NULL;
		delete[]d_H_j_const_pred;d_H_j_const_pred=NULL;

		delete[]d_HS_j_const_next;d_HS_j_const_next=NULL;
		delete[]d_HI_j_const_next;d_HI_j_const_next=NULL;
		delete[]d_HD_j_const_next;d_HD_j_const_next=NULL;
		delete[]d_H_j_const_next;d_H_j_const_next=NULL;

		delete[]d_H_edge_max;d_H_edge_max=NULL;

		d_alp_data->d_memory_size_in_MB+=(double)(sizeof(long int)*d_a_step*17)/mb_bytes;

		d_HS_i_const_pred=d_HS_i_const_pred_new;d_HS_i_const_pred_new=NULL;
		d_HI_i_const_pred=d_HI_i_const_pred_new;d_HI_i_const_pred_new=NULL;
		d_HD_i_const_pred=d_HD_i_const_pred_new;d_HD_i_const_pred_new=NULL;
		d_H_i_const_pred=d_H_i_const_pred_new;d_H_i_const_pred_new=NULL;

		d_HS_i_const_next=d_HS_i_const_next_new;d_HS_i_const_next_new=NULL;
		d_HI_i_const_next=d_HI_i_const_next_new;d_HI_i_const_next_new=NULL;
		d_HD_i_const_next=d_HD_i_const_next_new;d_HD_i_const_next_new=NULL;
		d_H_i_const_next=d_H_i_const_next_new;d_H_i_const_next_new=NULL;

		d_HS_j_const_pred=d_HS_j_const_pred_new;d_HS_j_const_pred_new=NULL;
		d_HI_j_const_pred=d_HI_j_const_pred_new;d_HI_j_const_pred_new=NULL;
		d_HD_j_const_pred=d_HD_j_const_pred_new;d_HD_j_const_pred_new=NULL;
		d_H_j_const_pred=d_H_j_const_pred_new;d_H_j_const_pred_new=NULL;

		d_HS_j_const_next=d_HS_j_const_next_new;d_HS_j_const_next_new=NULL;
		d_HI_j_const_next=d_HI_j_const_next_new;d_HI_j_const_next_new=NULL;
		d_HD_j_const_next=d_HD_j_const_next_new;d_HD_j_const_next_new=NULL;
		d_H_j_const_next=d_H_j_const_next_new;d_H_j_const_next_new=NULL;

		d_H_edge_max=d_H_edge_max_new;d_H_edge_max_new=NULL;

		
	}
	catch (...)
	{ 
		delete[]d_HS_i_const_pred_new;d_HS_i_const_pred_new=NULL;
		delete[]d_HI_i_const_pred_new;d_HI_i_const_pred_new=NULL;
		delete[]d_HD_i_const_pred_new;d_HD_i_const_pred_new=NULL;
		delete[]d_H_i_const_pred_new;d_H_i_const_pred_new=NULL;

		delete[]d_HS_i_const_next_new;d_HS_i_const_next_new=NULL;
		delete[]d_HI_i_const_next_new;d_HI_i_const_next_new=NULL;
		delete[]d_HD_i_const_next_new;d_HD_i_const_next_new=NULL;
		delete[]d_H_i_const_next_new;d_H_i_const_next_new=NULL;

		delete[]d_HS_j_const_pred_new;d_HS_j_const_pred_new=NULL;
		delete[]d_HI_j_const_pred_new;d_HI_j_const_pred_new=NULL;
		delete[]d_HD_j_const_pred_new;d_HD_j_const_pred_new=NULL;
		delete[]d_H_j_const_pred_new;d_H_j_const_pred_new=NULL;

		delete[]d_HS_j_const_next_new;d_HS_j_const_next_new=NULL;
		delete[]d_HI_j_const_next_new;d_HI_j_const_next_new=NULL;
		delete[]d_HD_j_const_next_new;d_HD_j_const_next_new=NULL;
		delete[]d_H_j_const_next_new;d_H_j_const_next_new=NULL;

		delete[]d_H_edge_max_new;d_H_edge_max_new=NULL;
		throw;
	};


}

void alp::increment_W_weights()
//the function calculates weigths for d_W_matr_len increased by 1
//assumes that letters are defined for d_W_matr_len
{
	if(d_W_matr_len==-1)
	{
		d_WS_ij_next=1.0;
		d_WI_ij_next=0.0;
		d_WD_ij_next=0.0;

		d_W_matr_len++;

		d_alp_weights->set_elem(0,1.0);

		return;
	};

	if(d_seqi_len<d_W_matr_len+1||d_seqj_len<d_W_matr_len+1)
	{
		throw error("Unexpected error in increment_W_weights\n",4);
	};

	if(d_W_matr_len+1>d_W_matr_a_len)
	{
		increment_W_matrix();
	};

	d_W_matr_len++;
	

	swap(d_WS_i_const_pred,d_WS_i_const_next);
	swap(d_WI_i_const_pred,d_WI_i_const_next);
	swap(d_WD_i_const_pred,d_WD_i_const_next);

	swap(d_WS_j_const_pred,d_WS_j_const_next);
	swap(d_WI_j_const_pred,d_WI_j_const_next);
	swap(d_WD_j_const_pred,d_WD_j_const_next);

	d_WS_ij_pred=d_WS_ij_next;
	d_WI_ij_pred=d_WI_ij_next;
	d_WD_ij_pred=d_WD_ij_next;

	long int d_W_matr_len_1=d_W_matr_len-1;
	long int d_W_matr_len_2=d_W_matr_len-2;

	//boundary conditions
	importance_sampling *&d_is_tmp=d_alp_data->d_is;

	d_WS_i_const_next[d_W_matr_len_1]=0;
	d_WS_j_const_next[d_W_matr_len_1]=0;

	d_WI_i_const_next[d_W_matr_len_1]=0;
	d_WD_j_const_next[d_W_matr_len_1]=0;

	double deg_tmp=degree(d_is_tmp->d_nu,d_W_matr_len_1);

	d_WD_i_const_next[d_W_matr_len_1]=d_is_tmp->d_mu_DS*deg_tmp;
	d_WI_j_const_next[d_W_matr_len_1]=d_is_tmp->d_mu_IS*deg_tmp;



	long int i;
	for(i=d_W_matr_len_2;i>=1;i--)
	{
		d_WS_i_const_next[i]=d_is_tmp->d_exp_s[d_seqi[d_W_matr_len_1]][d_seqj[d_W_matr_len_2-i]]*(d_is_tmp->d_eta*d_WS_i_const_pred[i]+d_is_tmp->d_mu_SI*d_WI_i_const_pred[i]+d_is_tmp->d_mu_SD*d_WD_i_const_pred[i]);
		d_WI_i_const_next[i]=d_is_tmp->d_mu_IS*d_WS_i_const_next[i+1]+d_is_tmp->d_nu*d_WI_i_const_next[i+1]+d_is_tmp->d_mu_ID*d_WD_i_const_next[i+1];
		d_WD_i_const_next[i]=d_is_tmp->d_mu_DS*d_WS_i_const_pred[i-1]+d_is_tmp->d_nu*d_WD_i_const_pred[i-1];

		d_WS_j_const_next[i]=d_is_tmp->d_exp_s[d_seqi[d_W_matr_len_2-i]][d_seqj[d_W_matr_len_1]]*(d_is_tmp->d_eta*d_WS_j_const_pred[i]+d_is_tmp->d_mu_SI*d_WI_j_const_pred[i]+d_is_tmp->d_mu_SD*d_WD_j_const_pred[i]);
		d_WI_j_const_next[i]=d_is_tmp->d_mu_IS*d_WS_j_const_pred[i-1]+d_is_tmp->d_nu*d_WI_j_const_pred[i-1]+d_is_tmp->d_mu_ID*d_WD_j_const_pred[i-1];
		d_WD_j_const_next[i]=d_is_tmp->d_mu_DS*d_WS_j_const_next[i+1]+d_is_tmp->d_nu*d_WD_j_const_next[i+1];
	};

	if(d_W_matr_len>1)
	{
		//copy of the previous lines with a modification for i-1
		i=0;
		d_WS_i_const_next[i]=d_is_tmp->d_exp_s[d_seqi[d_W_matr_len_1]][d_seqj[d_W_matr_len_2-i]]*(d_is_tmp->d_eta*d_WS_i_const_pred[i]+d_is_tmp->d_mu_SI*d_WI_i_const_pred[i]+d_is_tmp->d_mu_SD*d_WD_i_const_pred[i]);
		d_WI_i_const_next[i]=d_is_tmp->d_mu_IS*d_WS_i_const_next[i+1]+d_is_tmp->d_nu*d_WI_i_const_next[i+1]+d_is_tmp->d_mu_ID*d_WD_i_const_next[i+1];
		d_WD_i_const_next[i]=d_is_tmp->d_mu_DS*d_WS_ij_pred+d_is_tmp->d_nu*d_WD_ij_pred;

		d_WS_j_const_next[i]=d_is_tmp->d_exp_s[d_seqi[d_W_matr_len_2-i]][d_seqj[d_W_matr_len_1]]*(d_is_tmp->d_eta*d_WS_j_const_pred[i]+d_is_tmp->d_mu_SI*d_WI_j_const_pred[i]+d_is_tmp->d_mu_SD*d_WD_j_const_pred[i]);
		d_WI_j_const_next[i]=d_is_tmp->d_mu_IS*d_WS_ij_pred+d_is_tmp->d_nu*d_WI_ij_pred+d_is_tmp->d_mu_ID*d_WD_ij_pred;
		d_WD_j_const_next[i]=d_is_tmp->d_mu_DS*d_WS_j_const_next[i+1]+d_is_tmp->d_nu*d_WD_j_const_next[i+1];
	};


	d_WS_ij_next=d_is_tmp->d_exp_s[d_seqi[d_W_matr_len_1]][d_seqj[d_W_matr_len_1]]*(d_is_tmp->d_eta*d_WS_ij_pred+d_is_tmp->d_mu_SI*d_WI_ij_pred+d_is_tmp->d_mu_SD*d_WD_ij_pred);
	d_WI_ij_next=d_is_tmp->d_mu_IS*d_WS_i_const_next[0]+d_is_tmp->d_nu*d_WI_i_const_next[0]+d_is_tmp->d_mu_ID*d_WD_i_const_next[0];
	d_WD_ij_next=d_is_tmp->d_mu_DS*d_WS_j_const_next[0]+d_is_tmp->d_nu*d_WD_j_const_next[0];

	

}

double alp::degree(//returns x_^n_
double x_,
double n_)
{
	if(x_<0||n_<0)
	{
		throw error("Error - unexpected parameter in alp::degree\n",4);
	};

	if(x_==0)
	{
		if(n_==0)
		{
			return 1.0;
		}
		else
		{
			return 0.0;
		};
	};

	return exp(n_*log(x_));

}

void alp::increment_H_weights()
{
	if(d_alp_data->d_insertions_after_deletions)
	{
		increment_H_weights_with_insertions_after_deletions();
	}
	else
	{
		increment_H_weights_without_insertions_after_deletions();
	};
}

void alp::increment_H_weights_without_insertions_after_deletions()
//the function calculates alignment scores for d_H_matr_len increased by 1
//assumes that letters are defined for d_H_matr_len
{
	if(d_H_matr_len==-1)
	{
		d_HS_ij_next=0;
		d_HI_ij_next=0;
		d_HD_ij_next=0;
		d_H_ij_next=0;
		d_M=0;

		d_nalp=0;
		d_alp->set_elem(0,0);
		d_H_I->set_elem(0,0);
		d_H_J->set_elem(0,0);
		d_alp_pos->set_elem(0,0);

		d_cells_counts->increase_elem_by_1(0);

		d_H_matr_len++;

		d_alp_states->set_elem(d_nalp,NULL);
		save_state(d_alp_states->d_elem[d_nalp]);

		return;
	};

	if(d_seqi_len<d_H_matr_len+1||d_seqj_len<d_H_matr_len+1)
	{
		throw error("Unexpected error\n",4);
	};

	if(d_H_matr_len+1>d_H_matr_a_len)
	{
		increment_H_matrix();
	};

	d_H_matr_len++;
	

	swap(d_HS_i_const_pred,d_HS_i_const_next);
	swap(d_HI_i_const_pred,d_HI_i_const_next);
	swap(d_HD_i_const_pred,d_HD_i_const_next);
	swap(d_H_i_const_pred,d_H_i_const_next);

	swap(d_HS_j_const_pred,d_HS_j_const_next);
	swap(d_HI_j_const_pred,d_HI_j_const_next);
	swap(d_HD_j_const_pred,d_HD_j_const_next);
	swap(d_H_j_const_pred,d_H_j_const_next);

	d_HS_ij_pred=d_HS_ij_next;
	d_HI_ij_pred=d_HI_ij_next;
	d_HD_ij_pred=d_HD_ij_next;
	d_H_ij_pred=d_H_ij_next;

	long int d_H_matr_len_1=d_H_matr_len-1;
	long int d_H_matr_len_2=d_H_matr_len-2;

	//boundary conditions
	//long int gap_tmp=-d_alp_data->d_open-d_H_matr_len_1*d_alp_data->d_epen;
		long int gap_tmp1=-d_alp_data->d_open1-d_H_matr_len_1*d_alp_data->d_epen1;
		long int gap_tmp2=-d_alp_data->d_open2-d_H_matr_len_1*d_alp_data->d_epen2;


	d_HS_i_const_next[d_H_matr_len_1]=small_long;
	d_HS_j_const_next[d_H_matr_len_1]=small_long;

	d_HI_i_const_next[d_H_matr_len_1]=small_long;
	d_HD_j_const_next[d_H_matr_len_1]=small_long;


	d_HD_i_const_next[d_H_matr_len_1]=gap_tmp1;
	d_HI_j_const_next[d_H_matr_len_1]=gap_tmp2;

	d_H_i_const_next[d_H_matr_len_1]=gap_tmp1;
	d_H_j_const_next[d_H_matr_len_1]=gap_tmp2;

	long int i;
	for(i=d_H_matr_len_2;i>=1;i--)
	{
		d_HS_i_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_2-i]]+d_H_i_const_pred[i];
		d_HI_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i+1]-d_alp_data->d_open2,d_HI_i_const_next[i+1]-d_alp_data->d_epen2);
		d_HD_i_const_next[i]=alp_data::Tmax(d_HS_i_const_pred[i-1]-d_alp_data->d_open1,d_HD_i_const_pred[i-1]-d_alp_data->d_epen1);
		d_H_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i],d_HI_i_const_next[i],d_HD_i_const_next[i]);
	

		d_HS_j_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_2-i]][d_seqj[d_H_matr_len_1]]+d_H_j_const_pred[i];
		d_HI_j_const_next[i]=alp_data::Tmax(d_HS_j_const_pred[i-1]-d_alp_data->d_open2,d_HI_j_const_pred[i-1]-d_alp_data->d_epen2);
		d_HD_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i+1]-d_alp_data->d_open1,d_HD_j_const_next[i+1]-d_alp_data->d_epen1);
		d_H_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i],d_HI_j_const_next[i],d_HD_j_const_next[i]);
	};

	if(d_H_matr_len>1)
	{
		//copy of the previous lines with a modification for i-1
		i=0;
		d_HS_i_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_2-i]]+d_H_i_const_pred[i];
		d_HI_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i+1]-d_alp_data->d_open2,d_HI_i_const_next[i+1]-d_alp_data->d_epen2);
		d_HD_i_const_next[i]=alp_data::Tmax(d_HS_ij_pred-d_alp_data->d_open1,d_HD_ij_pred-d_alp_data->d_epen1);
		d_H_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i],d_HI_i_const_next[i],d_HD_i_const_next[i]);
	

		d_HS_j_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_2-i]][d_seqj[d_H_matr_len_1]]+d_H_j_const_pred[i];
		d_HI_j_const_next[i]=alp_data::Tmax(d_HS_ij_pred-d_alp_data->d_open2,d_HI_ij_pred-d_alp_data->d_epen2);
		d_HD_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i+1]-d_alp_data->d_open1,d_HD_j_const_next[i+1]-d_alp_data->d_epen1);
		d_H_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i],d_HI_j_const_next[i],d_HD_j_const_next[i]);
	};

	d_HS_ij_next=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_1]]+d_H_ij_pred;
	d_HI_ij_next=alp_data::Tmax(d_HS_i_const_next[0]-d_alp_data->d_open2,d_HI_i_const_next[0]-d_alp_data->d_epen2);
	d_HD_ij_next=alp_data::Tmax(d_HS_j_const_next[0]-d_alp_data->d_open1,d_HD_j_const_next[0]-d_alp_data->d_epen1);
	d_H_ij_next=alp_data::Tmax(d_HS_ij_next,d_HI_ij_next,d_HD_ij_next);

	d_cells_counts->increase_elem_by_1(d_H_ij_next);
	for(i=0;i<=d_H_matr_len_1;i++)
	{
		d_cells_counts->increase_elem_by_1(d_H_i_const_next[i]);
		d_cells_counts->increase_elem_by_1(d_H_j_const_next[i]);
	};


	long int tmp=d_H_ij_next;
	for(i=0;i<=d_H_matr_len_1;i++)
	{
		tmp=alp_data::Tmax(tmp,d_H_i_const_next[i]);
		tmp=alp_data::Tmax(tmp,d_H_j_const_next[i]);
	};

	d_H_edge_max[d_H_matr_len]=tmp;
	d_M=alp_data::Tmax(tmp,d_M);

	d_sentinel_i_next=d_H_matr_len_1;
	d_sentinel_j_next=d_H_matr_len_1;


	if(d_is_now)
	{
		long int i;
		if(tmp>d_alp->d_elem[d_nalp])
		{
			d_nalp++;
			d_alp->set_elem(d_nalp,tmp);
			d_alp_pos->set_elem(d_nalp,d_H_matr_len);

			d_alp_states->set_elem(d_nalp,NULL);
			save_state(d_alp_states->d_elem[d_nalp]);

		
			long int I=-1;
			long int J=-1;

			for(i=0;i<=d_H_matr_len_1;i++)
			{
				if(tmp==d_H_i_const_next[i])
				{
					I=i;
				};

				if(tmp==d_H_j_const_next[i])
				{
					J=i;
				};
			};

			d_H_I->set_elem(d_nalp,d_H_matr_len-I-1);
			d_H_J->set_elem(d_nalp,d_H_matr_len-J-1);


		};
	};


	check_time_function();


}


void alp::increment_H_weights_with_insertions_after_deletions()
//the function calculates alignment scores for d_H_matr_len increased by 1
//assumes that letters are defined for d_H_matr_len
{
	if(d_H_matr_len==-1)
	{
		d_HS_ij_next=0;
		d_HI_ij_next=0;
		d_HD_ij_next=0;
		d_H_ij_next=0;
		d_M=0;

		d_nalp=0;
		d_alp->set_elem(0,0);
		d_H_I->set_elem(0,0);
		d_H_J->set_elem(0,0);
		d_alp_pos->set_elem(0,0);

		d_cells_counts->increase_elem_by_1(0);

		d_H_matr_len++;

		d_alp_states->set_elem(d_nalp,NULL);
		save_state(d_alp_states->d_elem[d_nalp]);

		return;
	};

	if(d_seqi_len<d_H_matr_len+1||d_seqj_len<d_H_matr_len+1)
	{
		throw error("Unexpected error\n",4);
	};

	if(d_H_matr_len+1>d_H_matr_a_len)
	{
		increment_H_matrix();
	};

	d_H_matr_len++;
	

	swap(d_HS_i_const_pred,d_HS_i_const_next);
	swap(d_HI_i_const_pred,d_HI_i_const_next);
	swap(d_HD_i_const_pred,d_HD_i_const_next);
	swap(d_H_i_const_pred,d_H_i_const_next);

	swap(d_HS_j_const_pred,d_HS_j_const_next);
	swap(d_HI_j_const_pred,d_HI_j_const_next);
	swap(d_HD_j_const_pred,d_HD_j_const_next);
	swap(d_H_j_const_pred,d_H_j_const_next);

	d_HS_ij_pred=d_HS_ij_next;
	d_HI_ij_pred=d_HI_ij_next;
	d_HD_ij_pred=d_HD_ij_next;
	d_H_ij_pred=d_H_ij_next;

	long int d_H_matr_len_1=d_H_matr_len-1;
	long int d_H_matr_len_2=d_H_matr_len-2;

	//boundary conditions
	//long int gap_tmp=-d_alp_data->d_open-d_H_matr_len_1*d_alp_data->d_epen;
	long int gap_tmp1=-d_alp_data->d_open1-d_H_matr_len_1*d_alp_data->d_epen1;
	long int gap_tmp2=-d_alp_data->d_open2-d_H_matr_len_1*d_alp_data->d_epen2;


	d_HS_i_const_next[d_H_matr_len_1]=small_long;
	d_HS_j_const_next[d_H_matr_len_1]=small_long;

	d_HI_i_const_next[d_H_matr_len_1]=small_long;
	d_HD_j_const_next[d_H_matr_len_1]=small_long;


	d_HD_i_const_next[d_H_matr_len_1]=gap_tmp1;
	d_HI_j_const_next[d_H_matr_len_1]=gap_tmp2;

	d_H_i_const_next[d_H_matr_len_1]=gap_tmp1;
	d_H_j_const_next[d_H_matr_len_1]=gap_tmp2;

	long int i;
	for(i=d_H_matr_len_2;i>=1;i--)
	{
		d_HS_i_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_2-i]]+d_H_i_const_pred[i];
		d_HI_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i+1]-d_alp_data->d_open2,d_HI_i_const_next[i+1]-d_alp_data->d_epen2,d_HD_i_const_next[i+1]-d_alp_data->d_open2);
		d_HD_i_const_next[i]=alp_data::Tmax(d_HS_i_const_pred[i-1]-d_alp_data->d_open1,d_HD_i_const_pred[i-1]-d_alp_data->d_epen1);
		d_H_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i],d_HI_i_const_next[i],d_HD_i_const_next[i]);


		d_HS_j_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_2-i]][d_seqj[d_H_matr_len_1]]+d_H_j_const_pred[i];
		d_HI_j_const_next[i]=alp_data::Tmax(d_HS_j_const_pred[i-1]-d_alp_data->d_open2,d_HI_j_const_pred[i-1]-d_alp_data->d_epen2,d_HD_j_const_pred[i-1]-d_alp_data->d_open2);
		d_HD_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i+1]-d_alp_data->d_open1,d_HD_j_const_next[i+1]-d_alp_data->d_epen1);
		d_H_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i],d_HI_j_const_next[i],d_HD_j_const_next[i]);


	};

	if(d_H_matr_len>1)
	{
		//copy of the previous lines with a modification for i-1
		i=0;
		d_HS_i_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_2-i]]+d_H_i_const_pred[i];
		d_HI_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i+1]-d_alp_data->d_open2,d_HI_i_const_next[i+1]-d_alp_data->d_epen2,d_HD_i_const_next[i+1]-d_alp_data->d_open2);
		d_HD_i_const_next[i]=alp_data::Tmax(d_HS_ij_pred-d_alp_data->d_open1,d_HD_ij_pred-d_alp_data->d_epen1);
		d_H_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i],d_HI_i_const_next[i],d_HD_i_const_next[i]);
	


		d_HS_j_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_2-i]][d_seqj[d_H_matr_len_1]]+d_H_j_const_pred[i];
		d_HI_j_const_next[i]=alp_data::Tmax(d_HS_ij_pred-d_alp_data->d_open2,d_HI_ij_pred-d_alp_data->d_epen2,d_HD_ij_pred-d_alp_data->d_open2);
		d_HD_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i+1]-d_alp_data->d_open1,d_HD_j_const_next[i+1]-d_alp_data->d_epen1);
		d_H_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i],d_HI_j_const_next[i],d_HD_j_const_next[i]);


	};

	d_HS_ij_next=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_1]]+d_H_ij_pred;
	d_HI_ij_next=alp_data::Tmax(d_HS_i_const_next[0]-d_alp_data->d_open2,d_HI_i_const_next[0]-d_alp_data->d_epen2,d_HD_i_const_next[0]-d_alp_data->d_open2);
	d_HD_ij_next=alp_data::Tmax(d_HS_j_const_next[0]-d_alp_data->d_open1,d_HD_j_const_next[0]-d_alp_data->d_epen1);
	d_H_ij_next=alp_data::Tmax(d_HS_ij_next,d_HI_ij_next,d_HD_ij_next);



	d_cells_counts->increase_elem_by_1(d_H_ij_next);
	for(i=0;i<=d_H_matr_len_1;i++)
	{
		d_cells_counts->increase_elem_by_1(d_H_i_const_next[i]);
		d_cells_counts->increase_elem_by_1(d_H_j_const_next[i]);
	};


	long int tmp=d_H_ij_next;
	for(i=0;i<=d_H_matr_len_1;i++)
	{
		tmp=alp_data::Tmax(tmp,d_H_i_const_next[i]);
		tmp=alp_data::Tmax(tmp,d_H_j_const_next[i]);
	};

	d_H_edge_max[d_H_matr_len]=tmp;
	d_M=alp_data::Tmax(tmp,d_M);

	d_sentinel_i_next=d_H_matr_len_1;
	d_sentinel_j_next=d_H_matr_len_1;


	if(d_is_now)
	{
		long int i;
		if(tmp>d_alp->d_elem[d_nalp])
		{
			d_nalp++;
			d_alp->set_elem(d_nalp,tmp);
			d_alp_pos->set_elem(d_nalp,d_H_matr_len);

			d_alp_states->set_elem(d_nalp,NULL);
			save_state(d_alp_states->d_elem[d_nalp]);

		
			long int I=-1;
			long int J=-1;

			for(i=0;i<=d_H_matr_len_1;i++)
			{
				if(tmp==d_H_i_const_next[i])
				{
					I=i;
				};

				if(tmp==d_H_j_const_next[i])
				{
					J=i;
				};
			};

			d_H_I->set_elem(d_nalp,d_H_matr_len-I-1);
			d_H_J->set_elem(d_nalp,d_H_matr_len-J-1);


		};
	};


	check_time_function();


}

void alp::check_time_function()
{
	if(d_check_time_flag)
	{
		double time_after3;
		alp_data::get_current_time(time_after3);

		if((time_after3-d_alp_data->d_time_before1)>d_alp_data->d_max_time)
		{
			if(d_time_error_flag)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			}
			else
			{
				d_time_limit_flag=true;
				if(d_single_realiztion_calculation_flag)
				{
					throw error_for_single_realization();
				};
				return;
			};
		};

	};

	if(d_alp_data->d_max_time<=0&&d_alp_data->d_max_time_with_computation_parameters>0)
	{
		double time_after3;
		alp_data::get_current_time(time_after3);

		if((time_after3-d_alp_data->d_time_before1)>d_alp_data->d_max_time_with_computation_parameters)
		{
			throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
		};

	};
}

void alp::increment_H_weights_with_sentinels(
	long int diff_opt_)
//the function calculates alignment scores for d_H_matr_len increased by 1
//assumes that letters are defined for d_H_matr_len
//uses sentinels
{
	if(d_alp_data->d_insertions_after_deletions)
	{
		increment_H_weights_with_sentinels_with_insertions_after_deletions(diff_opt_);
	}
	else
	{
		increment_H_weights_with_sentinels_without_insertions_after_deletions(diff_opt_);
	};

}

void alp::increment_H_weights_with_sentinels_without_insertions_after_deletions(
	long int diff_opt_)
//the function calculates alignment scores for d_H_matr_len increased by 1
//assumes that letters are defined for d_H_matr_len
//uses sentinels
{
	if(d_H_matr_len==-1)
	{
		d_HS_ij_next=0;
		d_HI_ij_next=0;
		d_HD_ij_next=0;
		d_H_ij_next=0;
		d_M=0;

		d_nalp=0;
		d_alp->set_elem(0,0);
		d_H_I->set_elem(0,0);
		d_H_J->set_elem(0,0);
		d_alp_pos->set_elem(0,0);

		d_cells_counts->increase_elem_by_1(0);

		d_H_matr_len++;

		d_alp_states->set_elem(d_nalp,NULL);
		

		d_sentinel_i_next=0;
		d_sentinel_j_next=0;


		save_state(d_alp_states->d_elem[d_nalp]);

		

		return;
	};

	if(d_seqi_len<d_H_matr_len+1||d_seqj_len<d_H_matr_len+1)
	{
		throw error("Unexpected error\n",4);
	};


	if(d_H_matr_len+1>d_H_matr_a_len)
	{
		increment_H_matrix();
	};

	d_H_matr_len++;
	

	swap(d_HS_i_const_pred,d_HS_i_const_next);
	swap(d_HI_i_const_pred,d_HI_i_const_next);
	swap(d_HD_i_const_pred,d_HD_i_const_next);
	swap(d_H_i_const_pred,d_H_i_const_next);

	swap(d_HS_j_const_pred,d_HS_j_const_next);
	swap(d_HI_j_const_pred,d_HI_j_const_next);
	swap(d_HD_j_const_pred,d_HD_j_const_next);
	swap(d_H_j_const_pred,d_H_j_const_next);

	

	d_HS_ij_pred=d_HS_ij_next;
	d_HI_ij_pred=d_HI_ij_next;
	d_HD_ij_pred=d_HD_ij_next;
	d_H_ij_pred=d_H_ij_next;

	d_sentinel_i_pred=d_sentinel_i_next;
	d_sentinel_j_pred=d_sentinel_j_next;


	long int d_H_matr_len_1=d_H_matr_len-1;
	long int d_H_matr_len_2=d_H_matr_len-2;

	//boundary conditions
	long int sentinel_i_boundary=alp_data::Tmin(d_sentinel_i_pred+(long int)2,d_H_matr_len_1);
	long int sentinel_j_boundary=alp_data::Tmin(d_sentinel_j_pred+(long int)2,d_H_matr_len_1);



	d_HS_i_const_next[sentinel_i_boundary]=small_long;
	d_HS_j_const_next[sentinel_j_boundary]=small_long;

	d_HI_i_const_next[sentinel_i_boundary]=small_long;
	d_HD_j_const_next[sentinel_j_boundary]=small_long;


	d_HD_i_const_next[sentinel_i_boundary]=small_long;
	d_HI_j_const_next[sentinel_j_boundary]=small_long;

	d_H_i_const_next[sentinel_i_boundary]=small_long;
	d_H_j_const_next[sentinel_j_boundary]=small_long;

	long int i;
	for(i=sentinel_i_boundary-1;i>=1;i--)
	{
		d_HS_i_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_2-i]]+d_H_i_const_pred[i];
		d_HI_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i+1]-d_alp_data->d_open2,d_HI_i_const_next[i+1]-d_alp_data->d_epen2);
		d_HD_i_const_next[i]=alp_data::Tmax(d_HS_i_const_pred[i-1]-d_alp_data->d_open1,d_HD_i_const_pred[i-1]-d_alp_data->d_epen1);
		d_H_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i],d_HI_i_const_next[i],d_HD_i_const_next[i]);
	};

	for(i=sentinel_j_boundary-1;i>=1;i--)
	{
		d_HS_j_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_2-i]][d_seqj[d_H_matr_len_1]]+d_H_j_const_pred[i];
		d_HI_j_const_next[i]=alp_data::Tmax(d_HS_j_const_pred[i-1]-d_alp_data->d_open2,d_HI_j_const_pred[i-1]-d_alp_data->d_epen2);
		d_HD_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i+1]-d_alp_data->d_open1,d_HD_j_const_next[i+1]-d_alp_data->d_epen1);
		d_H_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i],d_HI_j_const_next[i],d_HD_j_const_next[i]);
	};


	if(d_H_matr_len>1)
	{
		//copy of the previous lines with a modification for i-1
		i=0;
		d_HS_i_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_2-i]]+d_H_i_const_pred[i];
		d_HI_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i+1]-d_alp_data->d_open2,d_HI_i_const_next[i+1]-d_alp_data->d_epen2);
		d_HD_i_const_next[i]=alp_data::Tmax(d_HS_ij_pred-d_alp_data->d_open1,d_HD_ij_pred-d_alp_data->d_epen1);
		d_H_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i],d_HI_i_const_next[i],d_HD_i_const_next[i]);
	

		d_HS_j_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_2-i]][d_seqj[d_H_matr_len_1]]+d_H_j_const_pred[i];
		d_HI_j_const_next[i]=alp_data::Tmax(d_HS_ij_pred-d_alp_data->d_open2,d_HI_ij_pred-d_alp_data->d_epen2);
		d_HD_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i+1]-d_alp_data->d_open1,d_HD_j_const_next[i+1]-d_alp_data->d_epen1);
		d_H_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i],d_HI_j_const_next[i],d_HD_j_const_next[i]);
	};

	d_HS_ij_next=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_1]]+d_H_ij_pred;
	d_HI_ij_next=alp_data::Tmax(d_HS_i_const_next[0]-d_alp_data->d_open2,d_HI_i_const_next[0]-d_alp_data->d_epen2);
	d_HD_ij_next=alp_data::Tmax(d_HS_j_const_next[0]-d_alp_data->d_open1,d_HD_j_const_next[0]-d_alp_data->d_epen1);
	d_H_ij_next=alp_data::Tmax(d_HS_ij_next,d_HI_ij_next,d_HD_ij_next);

	d_cells_counts->increase_elem_by_1(d_H_ij_next);
	for(i=0;i<=sentinel_i_boundary-1;i++)
	{
		d_cells_counts->increase_elem_by_1(d_H_i_const_next[i]);
	};

	for(i=0;i<=sentinel_j_boundary-1;i++)
	{
		d_cells_counts->increase_elem_by_1(d_H_j_const_next[i]);
	};

	long int tmp=d_H_ij_next;
	for(i=0;i<=sentinel_i_boundary-1;i++)
	{
		tmp=alp_data::Tmax(tmp,d_H_i_const_next[i]);
	};

	for(i=0;i<=sentinel_j_boundary-1;i++)
	{
		tmp=alp_data::Tmax(tmp,d_H_j_const_next[i]);
	};


	d_H_edge_max[d_H_matr_len]=tmp;
	d_M=alp_data::Tmax(tmp,d_M);


	{
		long int level=tmp-diff_opt_;
		long int i;
		d_sentinel_i_next=1;
		d_sentinel_j_next=1;
		for(i=sentinel_i_boundary-1;i>=1;i--)
		{
			if(d_H_i_const_next[i]>=level)
			{
				d_sentinel_i_next=i;
				break;
			};
		};

		for(i=sentinel_j_boundary-1;i>=1;i--)
		{
			if(d_H_j_const_next[i]>=level)
			{
				d_sentinel_j_next=i;
				break;
			};
		};
	};


	if(d_is_now)
	{
		long int i;
		if(tmp>d_alp->d_elem[d_nalp])
		{
			d_nalp++;
			d_alp->set_elem(d_nalp,tmp);
			d_alp_pos->set_elem(d_nalp,d_H_matr_len);

			d_alp_states->set_elem(d_nalp,NULL);
			save_state(d_alp_states->d_elem[d_nalp]);

		
			long int I=-1;
			long int J=-1;

			for(i=0;i<=sentinel_i_boundary-1;i++)
			{
				if(tmp==d_H_i_const_next[i])
				{
					I=i;
				};
			};

			for(i=0;i<=sentinel_j_boundary-1;i++)
			{
				if(tmp==d_H_j_const_next[i])
				{
					J=i;
				};
			};



			d_H_I->set_elem(d_nalp,d_H_matr_len-I-1);
			d_H_J->set_elem(d_nalp,d_H_matr_len-J-1);


		};
	};

	check_time_function();

}



void alp::increment_H_weights_with_sentinels_with_insertions_after_deletions(
	long int diff_opt_)
//the function calculates alignment scores for d_H_matr_len increased by 1
//assumes that letters are defined for d_H_matr_len
//uses sentinels
{
	if(d_H_matr_len==-1)
	{
		d_HS_ij_next=0;
		d_HI_ij_next=0;
		d_HD_ij_next=0;
		d_H_ij_next=0;
		d_M=0;

		d_nalp=0;
		d_alp->set_elem(0,0);
		d_H_I->set_elem(0,0);
		d_H_J->set_elem(0,0);
		d_alp_pos->set_elem(0,0);

		d_cells_counts->increase_elem_by_1(0);

		d_H_matr_len++;

		d_alp_states->set_elem(d_nalp,NULL);
		

		d_sentinel_i_next=0;
		d_sentinel_j_next=0;


		save_state(d_alp_states->d_elem[d_nalp]);

		

		return;
	};

	if(d_seqi_len<d_H_matr_len+1||d_seqj_len<d_H_matr_len+1)
	{
		throw error("Unexpected error\n",4);
	};


	if(d_H_matr_len+1>d_H_matr_a_len)
	{
		increment_H_matrix();
	};

	d_H_matr_len++;
	

	swap(d_HS_i_const_pred,d_HS_i_const_next);
	swap(d_HI_i_const_pred,d_HI_i_const_next);
	swap(d_HD_i_const_pred,d_HD_i_const_next);
	swap(d_H_i_const_pred,d_H_i_const_next);

	swap(d_HS_j_const_pred,d_HS_j_const_next);
	swap(d_HI_j_const_pred,d_HI_j_const_next);
	swap(d_HD_j_const_pred,d_HD_j_const_next);
	swap(d_H_j_const_pred,d_H_j_const_next);

	

	d_HS_ij_pred=d_HS_ij_next;
	d_HI_ij_pred=d_HI_ij_next;
	d_HD_ij_pred=d_HD_ij_next;
	d_H_ij_pred=d_H_ij_next;

	d_sentinel_i_pred=d_sentinel_i_next;
	d_sentinel_j_pred=d_sentinel_j_next;


	long int d_H_matr_len_1=d_H_matr_len-1;
	long int d_H_matr_len_2=d_H_matr_len-2;

	//boundary conditions
	long int sentinel_i_boundary=alp_data::Tmin(d_sentinel_i_pred+(long int)2,d_H_matr_len_1);
	long int sentinel_j_boundary=alp_data::Tmin(d_sentinel_j_pred+(long int)2,d_H_matr_len_1);



	d_HS_i_const_next[sentinel_i_boundary]=small_long;
	d_HS_j_const_next[sentinel_j_boundary]=small_long;

	d_HI_i_const_next[sentinel_i_boundary]=small_long;
	d_HD_j_const_next[sentinel_j_boundary]=small_long;


	d_HD_i_const_next[sentinel_i_boundary]=small_long;
	d_HI_j_const_next[sentinel_j_boundary]=small_long;

	d_H_i_const_next[sentinel_i_boundary]=small_long;
	d_H_j_const_next[sentinel_j_boundary]=small_long;

	long int i;
	for(i=sentinel_i_boundary-1;i>=1;i--)
	{
		d_HS_i_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_2-i]]+d_H_i_const_pred[i];
		d_HI_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i+1]-d_alp_data->d_open2,d_HI_i_const_next[i+1]-d_alp_data->d_epen2,d_HD_i_const_next[i+1]-d_alp_data->d_open2);
		d_HD_i_const_next[i]=alp_data::Tmax(d_HS_i_const_pred[i-1]-d_alp_data->d_open1,d_HD_i_const_pred[i-1]-d_alp_data->d_epen1);
		d_H_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i],d_HI_i_const_next[i],d_HD_i_const_next[i]);
	};

	for(i=sentinel_j_boundary-1;i>=1;i--)
	{
		d_HS_j_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_2-i]][d_seqj[d_H_matr_len_1]]+d_H_j_const_pred[i];
		d_HI_j_const_next[i]=alp_data::Tmax(d_HS_j_const_pred[i-1]-d_alp_data->d_open2,d_HI_j_const_pred[i-1]-d_alp_data->d_epen2,d_HD_j_const_pred[i-1]-d_alp_data->d_open2);
		d_HD_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i+1]-d_alp_data->d_open1,d_HD_j_const_next[i+1]-d_alp_data->d_epen1);
		d_H_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i],d_HI_j_const_next[i],d_HD_j_const_next[i]);
	};


	if(d_H_matr_len>1)
	{
		//copy of the previous lines with a modification for i-1
		i=0;
		d_HS_i_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_2-i]]+d_H_i_const_pred[i];
		d_HI_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i+1]-d_alp_data->d_open2,d_HI_i_const_next[i+1]-d_alp_data->d_epen2,d_HD_i_const_next[i+1]-d_alp_data->d_open2);
		d_HD_i_const_next[i]=alp_data::Tmax(d_HS_ij_pred-d_alp_data->d_open1,d_HD_ij_pred-d_alp_data->d_epen1);
		d_H_i_const_next[i]=alp_data::Tmax(d_HS_i_const_next[i],d_HI_i_const_next[i],d_HD_i_const_next[i]);
	

		d_HS_j_const_next[i]=d_alp_data->d_smatr[d_seqi[d_H_matr_len_2-i]][d_seqj[d_H_matr_len_1]]+d_H_j_const_pred[i];
		d_HI_j_const_next[i]=alp_data::Tmax(d_HS_ij_pred-d_alp_data->d_open2,d_HI_ij_pred-d_alp_data->d_epen2,d_HD_ij_pred-d_alp_data->d_open2);
		d_HD_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i+1]-d_alp_data->d_open1,d_HD_j_const_next[i+1]-d_alp_data->d_epen1);
		d_H_j_const_next[i]=alp_data::Tmax(d_HS_j_const_next[i],d_HI_j_const_next[i],d_HD_j_const_next[i]);
	};

	d_HS_ij_next=d_alp_data->d_smatr[d_seqi[d_H_matr_len_1]][d_seqj[d_H_matr_len_1]]+d_H_ij_pred;
	d_HI_ij_next=alp_data::Tmax(d_HS_i_const_next[0]-d_alp_data->d_open2,d_HI_i_const_next[0]-d_alp_data->d_epen2,d_HD_i_const_next[0]-d_alp_data->d_open2);
	d_HD_ij_next=alp_data::Tmax(d_HS_j_const_next[0]-d_alp_data->d_open1,d_HD_j_const_next[0]-d_alp_data->d_epen1);
	d_H_ij_next=alp_data::Tmax(d_HS_ij_next,d_HI_ij_next,d_HD_ij_next);

	d_cells_counts->increase_elem_by_1(d_H_ij_next);
	for(i=0;i<=sentinel_i_boundary-1;i++)
	{
		d_cells_counts->increase_elem_by_1(d_H_i_const_next[i]);
	};

	for(i=0;i<=sentinel_j_boundary-1;i++)
	{
		d_cells_counts->increase_elem_by_1(d_H_j_const_next[i]);
	};

	long int tmp=d_H_ij_next;
	for(i=0;i<=sentinel_i_boundary-1;i++)
	{
		tmp=alp_data::Tmax(tmp,d_H_i_const_next[i]);
	};

	for(i=0;i<=sentinel_j_boundary-1;i++)
	{
		tmp=alp_data::Tmax(tmp,d_H_j_const_next[i]);
	};


	d_H_edge_max[d_H_matr_len]=tmp;
	d_M=alp_data::Tmax(tmp,d_M);


	{
		long int level=tmp-diff_opt_;
		long int i;
		d_sentinel_i_next=1;
		d_sentinel_j_next=1;
		for(i=sentinel_i_boundary-1;i>=1;i--)
		{
			if(d_H_i_const_next[i]>=level)
			{
				d_sentinel_i_next=i;
				break;
			};
		};

		for(i=sentinel_j_boundary-1;i>=1;i--)
		{
			if(d_H_j_const_next[i]>=level)
			{
				d_sentinel_j_next=i;
				break;
			};
		};
	};


	if(d_is_now)
	{
		long int i;
		if(tmp>d_alp->d_elem[d_nalp])
		{
			d_nalp++;
			d_alp->set_elem(d_nalp,tmp);
			d_alp_pos->set_elem(d_nalp,d_H_matr_len);

			d_alp_states->set_elem(d_nalp,NULL);
			save_state(d_alp_states->d_elem[d_nalp]);

		
			long int I=-1;
			long int J=-1;

			for(i=0;i<=sentinel_i_boundary-1;i++)
			{
				if(tmp==d_H_i_const_next[i])
				{
					I=i;
				};
			};

			for(i=0;i<=sentinel_j_boundary-1;i++)
			{
				if(tmp==d_H_j_const_next[i])
				{
					J=i;
				};
			};



			d_H_I->set_elem(d_nalp,d_H_matr_len-I-1);
			d_H_J->set_elem(d_nalp,d_H_matr_len-J-1);


		};
	};

	check_time_function();

}



void alp::restore_state(
state * &state_)
{
	d_M=state_->d_M;
	d_H_matr_len=state_->d_H_matr_len;

	if(d_H_matr_len<0)
	{
		throw error("Unexpected error\n",4);
	};


	d_is_now=false;

	delete d_cells_counts;d_cells_counts=NULL;

	d_cells_counts=new array<long int>(d_alp_data);
	alp_data::assert_mem(d_cells_counts);


	d_cells_counts->set_elems(state_->d_cells_counts);


	d_HS_ij_next=state_->d_HS_ij_next;
	d_HI_ij_next=state_->d_HI_ij_next;
	d_HD_ij_next=state_->d_HD_ij_next;
	d_H_ij_next=state_->d_H_ij_next;

	long int i;
	for(i=0;i<d_H_matr_len;i++)
	{
		d_HS_i_const_next[i]=state_->d_HS_i_const_next[i];
		d_HI_i_const_next[i]=state_->d_HI_i_const_next[i];
		d_HD_i_const_next[i]=state_->d_HD_i_const_next[i];
		d_H_i_const_next[i]=state_->d_H_i_const_next[i];
		d_HS_j_const_next[i]=state_->d_HS_j_const_next[i];
		d_HI_j_const_next[i]=state_->d_HI_j_const_next[i];
		d_HD_j_const_next[i]=state_->d_HD_j_const_next[i];
		d_H_j_const_next[i]=state_->d_H_j_const_next[i];
	};

	d_sentinel_i_next=state_->d_sentinel_i_next;
	d_sentinel_j_next=state_->d_sentinel_j_next;


}

state::state()
{
	d_cells_counts=NULL;

	d_HS_i_const_next=NULL;
	d_HI_i_const_next=NULL;
	d_HD_i_const_next=NULL;
	d_H_i_const_next=NULL;

	d_HS_j_const_next=NULL;
	d_HI_j_const_next=NULL;
	d_HD_j_const_next=NULL;
	d_H_j_const_next=NULL;

}

void alp::save_state(
state * &state_)
{
	if(d_H_matr_len<0)
	{
		throw error("Unexpected error\n",4);
	};

	state_=new state;
	alp_data::assert_mem(state_);

	d_alp_data->d_memory_size_in_MB+=(double)(sizeof(state))/mb_bytes;


	state_->d_M=d_M;

	state_->d_cells_counts=new array<long int>(d_alp_data);
	alp_data::assert_mem(state_->d_cells_counts);

	d_alp_data->d_memory_size_in_MB+=(double)(sizeof(array<long int>))/mb_bytes;

	state_->d_cells_counts->set_elems(d_cells_counts);

	state_->d_H_matr_len=d_H_matr_len;

	state_->d_HS_ij_next=d_HS_ij_next;
	state_->d_HI_ij_next=d_HI_ij_next;
	state_->d_HD_ij_next=d_HD_ij_next;
	state_->d_H_ij_next=d_H_ij_next;

	if(d_H_matr_len==0)
	{
		state_->d_HS_i_const_next=NULL;;
		state_->d_HI_i_const_next=NULL;;
		state_->d_HD_i_const_next=NULL;;
		state_->d_H_i_const_next=NULL;;
		state_->d_HS_j_const_next=NULL;;
		state_->d_HI_j_const_next=NULL;;
		state_->d_HD_j_const_next=NULL;;
		state_->d_H_j_const_next=NULL;


	}
	else
	{
		state_->d_HS_i_const_next=new long int[d_H_matr_len];
		alp_data::assert_mem(state_->d_HS_i_const_next);

		state_->d_HI_i_const_next=new long int[d_H_matr_len];
		alp_data::assert_mem(state_->d_HI_i_const_next);

		state_->d_HD_i_const_next=new long int[d_H_matr_len];
		alp_data::assert_mem(state_->d_HD_i_const_next);

		state_->d_H_i_const_next=new long int[d_H_matr_len];
		alp_data::assert_mem(state_->d_H_i_const_next);

		state_->d_HS_j_const_next=new long int[d_H_matr_len];
		alp_data::assert_mem(state_->d_HS_j_const_next);

		state_->d_HI_j_const_next=new long int[d_H_matr_len];
		alp_data::assert_mem(state_->d_HI_j_const_next);

		state_->d_HD_j_const_next=new long int[d_H_matr_len];
		alp_data::assert_mem(state_->d_HD_j_const_next);

		state_->d_H_j_const_next=new long int[d_H_matr_len];
		alp_data::assert_mem(state_->d_H_j_const_next);

		d_alp_data->d_memory_size_in_MB+=8.0*(double)(d_H_matr_len*sizeof(long int))/mb_bytes;

		long int i;
		for(i=0;i<d_H_matr_len;i++)
		{
			state_->d_HS_i_const_next[i]=d_HS_i_const_next[i];
			state_->d_HI_i_const_next[i]=d_HI_i_const_next[i];
			state_->d_HD_i_const_next[i]=d_HD_i_const_next[i];
			state_->d_H_i_const_next[i]=d_H_i_const_next[i];
			state_->d_HS_j_const_next[i]=d_HS_j_const_next[i];
			state_->d_HI_j_const_next[i]=d_HI_j_const_next[i];
			state_->d_HD_j_const_next[i]=d_HD_j_const_next[i];
			state_->d_H_j_const_next[i]=d_H_j_const_next[i];
		};

	};

	state_->d_sentinel_i_next=d_sentinel_i_next;
	state_->d_sentinel_j_next=d_sentinel_j_next;

}

void alp::kill_upto_level(
long int M_min_,
long int M_level_,
long int *M_upper_level_)
{
	if(d_is_now)
	{
		while(d_alp->d_elem[d_nalp]<M_min_)
		{
			simulate_next_alp();
			if(!d_success)
			{
				return;
			};
		};
		d_is_now=false;

		long int i;
		d_nalp_killing=-1;
		for(i=0;i<=d_nalp;i++)
		{
			if(d_alp->d_elem[i]>=M_min_)
			{
				d_nalp_killing=i;
				break;
			};
		};

		if(d_nalp_killing==-1)
		{
			throw error("Unexpected error\n",4);
		};

		restore_state(d_alp_states->d_elem[d_nalp_killing]);

	};

	while(d_H_edge_max[d_H_matr_len]>=M_level_)
	{
		if(d_H_matr_len+1>=d_alp_data->d_dim1_tmp)
		{
			d_success=false;
			return;
		};

		if(M_upper_level_)
		{
			if(d_H_edge_max[d_H_matr_len]>(*M_upper_level_))
			{
				d_success=false;
				return;
			};
		};
		

		if(d_H_matr_len+1>d_seq_a_len)
		{
			increment_sequences();
		};


		d_seqi_len=d_seqj_len=d_H_matr_len+1;
		d_seqi[d_seqi_len-1]=random_AA1();
		d_seqj[d_seqj_len-1]=random_AA2();

		if(d_sentinels_flag)
		{
			increment_H_weights_with_sentinels(d_diff_opt);
		}
		else
		{
			increment_H_weights();
		};

		if(d_time_limit_flag)
		{
			d_success=false;
			return;
		};

	};

	d_success=true;
}

double alp::John2_weight_calculation(
long int length_)//calculation of weigths for the importance sampling
{
	if(length_==0)
	{
		return 1.0;
	};

	if(d_W_matr_len>length_)
	{
		throw error("Error - unexpected parameter in alp::John2_weight_calculation\n",4);
	};

	while(d_W_matr_len<length_)
	{
		increment_W_weights();
	};

	importance_sampling *&d_is_tmp=d_alp_data->d_is;

	long int d_W_matr_len_1=d_W_matr_len-1;

	double US=0;
	double UD=0;
	double UI=d_WI_j_const_next[d_W_matr_len_1]/(1-(d_is_tmp->d_nu));

	double VS=0;
	double VI=0;
	double VD=d_WD_i_const_next[d_W_matr_len_1]/(1-(d_is_tmp->d_nu));

	long int j;
	for(j=1;j<=length_-1;j++)
	{
		double US_next=d_alp_data->d_r_i_dot[d_seqi[j-1]]*((d_is_tmp->d_eta)*US+(d_is_tmp->d_mu_SI)*UI+(d_is_tmp->d_mu_SD)*UD)+d_WS_j_const_next[d_W_matr_len_1-j];
		double UD_next=((d_is_tmp->d_mu_DS)*US+(d_is_tmp->d_nu)*UD);
		double UI_next=((d_is_tmp->d_mu_IS)*US_next+(d_is_tmp->d_mu_ID)*UD_next+d_WI_j_const_next[d_W_matr_len_1-j])/(1-(d_is_tmp->d_nu));

		double VS_next=d_alp_data->d_r_dot_j[d_seqj[j-1]]*((d_is_tmp->d_eta)*VS+(d_is_tmp->d_mu_SI)*VI+(d_is_tmp->d_mu_SD)*VD)+d_WS_i_const_next[d_W_matr_len_1-j];
		double VI_next=((d_is_tmp->d_mu_IS)*VS+(d_is_tmp->d_mu_ID)*VD+(d_is_tmp->d_nu)*VI);
		double VD_next=((d_is_tmp->d_mu_DS)*VS_next+d_WD_i_const_next[d_W_matr_len_1-j])/(1-(d_is_tmp->d_nu));

		US=US_next;
		UD=UD_next;
		UI=UI_next;

		VS=VS_next;
		VD=VD_next;
		VI=VI_next;
	};

	//copy
	j=length_;
	double US_next=d_alp_data->d_r_i_dot[d_seqi[j-1]]*((d_is_tmp->d_eta)*US+(d_is_tmp->d_mu_SI)*UI+(d_is_tmp->d_mu_SD)*UD)+d_WS_ij_next;
	double UD_next=((d_is_tmp->d_mu_DS)*US+(d_is_tmp->d_nu)*UD);
	double UI_next=((d_is_tmp->d_mu_IS)*US_next+(d_is_tmp->d_mu_ID)*UD_next+d_WI_ij_next)/(1-(d_is_tmp->d_nu));

	double VS_next=d_alp_data->d_r_dot_j[d_seqj[j-1]]*((d_is_tmp->d_eta)*VS+(d_is_tmp->d_mu_SI)*VI+(d_is_tmp->d_mu_SD)*VD)+d_WS_ij_next;
	double VI_next=((d_is_tmp->d_mu_IS)*VS+(d_is_tmp->d_mu_ID)*VD+(d_is_tmp->d_nu)*VI);
	double VD_next=((d_is_tmp->d_mu_DS)*VS_next+d_WD_ij_next)/(1-(d_is_tmp->d_nu));

	US=US_next;
	UD=UD_next;
	UI=UI_next;

	VS=VS_next;
	VD=VD_next;
	VI=VI_next;


	double weight=-d_WS_ij_next+US+UD+VS+VI;




	if(weight==0)
	{
		throw error("Unexpected error\n",4);
	};
	weight=1.0/weight;

	return weight;

}

void alp::simulate_next_alp()//simulates next ALP
{
	if(!d_success)
	{
		return;
	};

	if(!d_is_now)
	{
		throw error("Unexpected error - ALP can be generated only in the importance sampling mode\n",4);
	};

	long int target_nalp=d_nalp+1;

	while(d_nalp<target_nalp)
	{
		long int k=alp_data::Tmin(d_seqi_len,d_seqj_len);
		//std::cout<<k<<"\t"<<this->d_H_edge_max[k]<<std::endl;

		while(alp_data::Tmin(d_seqi_len,d_seqj_len)!=k+1)
		{
			bool success=one_step_of_importance_sampling_without_weight_calculation(
			d_alp_data->d_dim1_tmp,
			d_alp_data->d_dim2_tmp);


			check_time_function();

			if(!success)
			{
				d_success=false;
				return;
			};
		};

		if(d_sentinels_flag)
		{
			increment_H_weights_with_sentinels(d_diff_opt);
		}
		else
		{
			increment_H_weights();
		};

		if(d_time_limit_flag)
		{
			d_success=false;
			return;
		};


		increment_W_weights();
	};

	double weight=John2_weight_calculation(alp_data::Tmin(d_seqi_len,d_seqj_len));
	if(weight<=0)
	{
		throw error("Unexpected error\n",4);
	};

	d_alp_weights->set_elem(d_nalp,weight);

}

void alp::simulate_alp_upto_the_given_number(//simulates ALP upto the given number nalp_ including
long int nalp_)
{
	d_sentinels_flag=false;
	while(d_nalp<nalp_)
	{
		simulate_next_alp();
		if(!d_success)
		{
			return;
		};
	};
}

void alp::simulate_alp_upto_the_given_level(//simulates ALP upto the given level M_min_ including
long int M_min_)
{
	d_sentinels_flag=false;
	while(d_alp->d_elem[d_nalp]<M_min_)
	{
		simulate_next_alp();
		if(!d_success)
		{
			return;
		};
	};
	d_nalp_killing=d_nalp;
}

