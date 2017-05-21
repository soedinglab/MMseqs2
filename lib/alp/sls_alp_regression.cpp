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

File name: sls_alp_regression.cpp

Author: Sergey Sheetlin

Contents: Regression methods

******************************************************************************/

#include "sls_alp_regression.hpp"

using namespace Sls;
using namespace std;


//-------------------------Solution of arbitrary equations----------------------------
void alp_reg::find_tetta_general(
function_type *func_,
void* func_pointer_,
double a_,//[a,b] is the interval for search of equation solution
double b_,
long int n_partition_,
double eps_,
std::vector<double> &res_
)
{
	res_.resize(0);
	vector<long int> intervals(0);

	if(n_partition_<=0)
	{
		throw error("Error in alp_reg::find_tetta_general\n",4);
		
	};
	
	long int i;
	double h=(b_-a_)/n_partition_;
	double x1,x2 = 0.0;
	for(i=0;i<n_partition_;i++)
	{
		if(i==0)
		{
		x1=(*func_)(//calculation of E(exp(tetta*X))
			a_+i*h,func_pointer_);
		if(fabs(x1)<eps_)
		{
			res_.push_back(a_+i*h);
		};

		}
		else
		{
			x1=x2;
		};


		x2=(*func_)(//calculation of E(exp(tetta*X))
			a_+(i+1)*h,
			func_pointer_);

		if(fabs(x2)<eps_)
		{
			res_.push_back(a_+(i+1)*h);
		};

		if((x1*x2<0)&&(fabs(x1)>=eps_&&fabs(x2)>=eps_))
		{
			intervals.push_back(i);
		};
	};

	for(i=0;i<(long int)intervals.size();i++)
	{
	double sol=find_single_tetta_general(
		func_,
		func_pointer_,
		a_+intervals[i]*h,//[a,b] is the interval for search of equation solution
		a_+(1+intervals[i])*h,
		eps_);
	res_.push_back(sol);
	};

	sort(res_.begin(),res_.end());

	return;
}

double alp_reg::find_single_tetta_general(
function_type *func_,
void* func_pointer_,
double a_,//[a,b] is the interval for search of equation solution
double b_,
double eps_
)
{
	if(b_<a_)
	{
		throw error("Error in alp_reg::find_single_tetta_general\n",4);
	};

	double x1=a_;
	double x2=b_;
	double precision=(x2-x1)/2;

	double y1=(*func_)(//calculation of E(exp(tetta*X))
			x1,
			func_pointer_);
	if(fabs(y1)<eps_)
	{
		return x1;
	};


	double y2=(*func_)(//calculation of E(exp(tetta*X))
			x2,
			func_pointer_);
	if(fabs(y2)<eps_)
	{
		return x2;
	};

	while(precision>eps_)
	{
		double x12=(x1+x2)/2;
		
		double y12=(*func_)(//calculation of E(exp(tetta*X))
		x12,
		func_pointer_);

		if(fabs(y12)<eps_)
		{
			return x12;
		};

		if(y12*y1<0)
		{
			x2=x12;
			y2=y12;
		}
		else
		{
			x1=x12;
			y1=y12;
		};
		precision=(x2-x1)/2;
	};

	return (x1+x2)/2;
}

//------------regression-----------------------------------

void alp_reg::correction_of_errors(
double *errors_,
long int number_of_elements_)
{

	if(number_of_elements_<=0)
	{
		throw error("Unexpected error\n",4);
	};
	

	double average_error=0;
	long int i;
	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]<0)
		{
			throw error("Error in alp_reg::correction_of_errors: input error in the regression model is less than 0\n",4);
		};

		average_error+=errors_[i];
	};

	average_error/=(double)number_of_elements_;

	double error_eps;
		
	if(average_error<=0)
	{
		error_eps=1e-50;
	}
	else
	{
		error_eps=average_error;
	};
		
	
	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]==0)
		{
			errors_[i]=error_eps;
		};
	};
}

void alp_reg::robust_regression_sum_with_cut_LSM(
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
bool &res_was_calculated_)
{

	if(number_of_elements_<2)
	{
		throw error("Unexpected error\n",4);
	};

	correction_of_errors(errors_,number_of_elements_);

	//minization of the function
	double c=y_*y_;


	long int k1_start,k1_end;
	long int k2_start,k2_end;

	if(cut_left_tail_&&cut_right_tail_)
	{
		k1_start=0;
		k1_end=number_of_elements_-1;

		k2_start=0;
		k2_end=number_of_elements_-1;

	}
	else
	{
		if(cut_left_tail_&&!cut_right_tail_)
		{
			k1_start=0;
			k1_end=number_of_elements_-1;

			k2_start=number_of_elements_-1;
			k2_end=number_of_elements_-1;

		}
		else
		{
			if(!cut_left_tail_&&cut_right_tail_)
			{
				k1_start=0;
				k1_end=0;

				k2_start=0;
				k2_end=number_of_elements_-1;
			}
			else
			{
				k1_start=0;
				k1_end=0;

				k2_start=number_of_elements_-1;
				k2_end=number_of_elements_-1;

			};
		};
	};
	
	long int k1_opt=0,k2_opt=0;

	double func_opt=DBL_MAX;
	double beta0_opt=0;
	double beta1_opt=0;
	double beta0_opt_error=0;
	double beta1_opt_error=0;

	long int k1,k2;

	
	res_was_calculated_=false;


	for(k1=k1_start;k1<=k1_end;k1++)
	{

		for(k2=sls_basic::Tmax(k1+(long int)1,sls_basic::Tmax(k1,k2_start)+min_length_);k2<=k2_end;k2++)
		{
			double beta0_opt_tmp,beta1_opt_tmp,beta0_opt_error_tmp,beta1_opt_error_tmp;
			bool res_was_calculated;

			double tmp=function_for_robust_regression_sum_with_cut_LSM(
				values_+k1,
				errors_+k1,
				k2-k1+1,
				k1,
				c,
				beta0_opt_tmp,
				beta1_opt_tmp,
				beta0_opt_error_tmp,
				beta1_opt_error_tmp,
				res_was_calculated);



			if(tmp<func_opt&&res_was_calculated)
			{
				func_opt=tmp;
				beta0_opt=beta0_opt_tmp;
				beta1_opt=beta1_opt_tmp;
				beta0_opt_error=beta0_opt_error_tmp;
				beta1_opt_error=beta1_opt_error_tmp;
				k1_opt=k1;
				k2_opt=k2;
				res_was_calculated_=true;
			};

		};
	};




	if(res_was_calculated_)
	{
		beta0_=beta0_opt;
		beta1_=beta1_opt;
		beta0_error_=beta0_opt_error;
		beta1_error_=beta1_opt_error;
		k1_opt_=k1_opt;
		k2_opt_=k2_opt;
	};

	
}

double alp_reg::function_for_robust_regression_sum_with_cut_LSM(
double *values_,
double *errors_,
long int number_of_elements_,
long int k_start_,
double c_,
double &beta0_,
double &beta1_,
double &beta0_error_,
double &beta1_error_,
bool &res_was_calculated_)
{
	long int i;
	double a11=0;
	double a12=0;
	double a21=0;
	double a22=0;
	double y1=0;
	double y2=0;

	double y1_error=0;
	double y2_error=0;

	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]!=0)
		{
			double tmp=1.0/(errors_[i]*errors_[i]);

			a11+=tmp;
			a12+=(double)(k_start_+i)*tmp;
			a22+=(double)((k_start_+i)*(k_start_+i))*tmp;
			y1+=values_[i]*tmp;
			y1_error+=tmp*tmp*errors_[i]*errors_[i];
			y2+=(double)(k_start_+i)*values_[i]*tmp;
			y2_error+=(double)(k_start_+i)*(double)(k_start_+i)*tmp*tmp*errors_[i]*errors_[i];
		};
	};

	a21=a12;
	y1_error=alp_reg::sqrt_for_errors(y1_error);
	y2_error=alp_reg::sqrt_for_errors(y2_error);

	double eps=1e-10*sls_basic::Tmax(fabs(a11*a22),fabs(a21*a12));

	double den=a11*a22-a21*a12;
	if(fabs(den)<=eps)
	{
		res_was_calculated_=false;
		return 0;
	}
	else
	{
		res_was_calculated_=true;
	};

	beta0_=(y1*a22-a12*y2)/den;
	beta1_=(a11*y2-a21*y1)/den;

	beta0_error_=sqrt(y1_error*y1_error*a22*a22+a12*a12*y2_error*y2_error)/den;
	beta1_error_=sqrt(a11*a11*y2_error*y2_error+a21*a21*y1_error*y1_error)/den;


	double res=0;
	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]!=0)
		{
			double tmp=(beta0_+beta1_*(i+k_start_)-values_[i])/errors_[i];
			res+=tmp*tmp-c_;
		};
	};

	return res;
}

void alp_reg::robust_regression_sum_with_cut_LSM_beta1_is_defined(
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
bool &res_was_calculated_)
{

	correction_of_errors(errors_,number_of_elements_);

	//minization of the function
	double c=y_*y_;


	long int k1_start,k1_end;
	long int k2_start,k2_end;

	if(cut_left_tail_&&cut_right_tail_)
	{
		k1_start=0;
		k1_end=number_of_elements_-1;

		k2_start=0;
		k2_end=number_of_elements_-1;

	}
	else
	{
		if(cut_left_tail_&&!cut_right_tail_)
		{
			k1_start=0;
			k1_end=number_of_elements_-1;

			k2_start=number_of_elements_-1;
			k2_end=number_of_elements_-1;

		}
		else
		{
			if(!cut_left_tail_&&cut_right_tail_)
			{
				k1_start=0;
				k1_end=0;

				k2_start=0;
				k2_end=number_of_elements_-1;
			}
			else
			{
				k1_start=0;
				k1_end=0;

				k2_start=number_of_elements_-1;
				k2_end=number_of_elements_-1;

			};
		};
	};
	
	long int k1_opt=0,k2_opt=0;

	double func_opt=DBL_MAX;
	double beta0_opt=0;
	double beta0_opt_error=0;

	long int k1,k2;

	res_was_calculated_=false;


	for(k1=k1_start;k1<=k1_end;k1++)
	{

		for(k2=sls_basic::Tmax(k1,k2_start)+min_length_;k2<=k2_end;k2++)
		{
			double beta0_opt_tmp,beta1_opt_tmp,beta0_opt_error_tmp,beta1_opt_error_tmp;
			bool res_was_calculated;

			beta1_opt_tmp=beta1_;
			beta1_opt_error_tmp=beta1_error_;

			double tmp=function_for_robust_regression_sum_with_cut_LSM_beta1_is_defined(
				values_+k1,
				errors_+k1,
				k2-k1+1,
				k1,
				c,
				beta0_opt_tmp,
				beta1_opt_tmp,
				beta0_opt_error_tmp,
				beta1_opt_error_tmp,
				res_was_calculated);

			if(tmp<func_opt&&res_was_calculated)
			{
				func_opt=tmp;
				beta0_opt=beta0_opt_tmp;
				beta0_opt_error=beta0_opt_error_tmp;
				k1_opt=k1;
				k2_opt=k2;
				res_was_calculated_=true;
			};

		};
	};

	if(res_was_calculated_)
	{
		beta0_=beta0_opt;
		beta0_error_=beta0_opt_error;
		k1_opt_=k1_opt;
		k2_opt_=k2_opt;
	};

	
}

double alp_reg::function_for_robust_regression_sum_with_cut_LSM_beta1_is_defined(
double *values_,
double *errors_,
long int number_of_elements_,
long int k_start_,
double c_,
double &beta0_,
double beta1_,
double &beta0_error_,
double beta1_error_,
bool &res_was_calculated_)
{
	long int i;
	double a11=0;
	double y1=0;

	double y1_error=0;


	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]!=0)
		{
			double tmp=1.0/(errors_[i]*errors_[i]);

			a11+=tmp;
			y1+=(values_[i]-(double)(k_start_+i)*beta1_)*tmp;
			double error_tmp=errors_[i]*errors_[i]+(double)(k_start_+i)*(double)(k_start_+i)*beta1_error_*beta1_error_;
			y1_error+=tmp*tmp*error_tmp;
		};
	};

	y1_error=sqrt(y1_error);

	double eps=1e-10*fabs(a11);

	double den=a11;
	if(fabs(den)<=eps)
	{
		res_was_calculated_=false;
		return 0;
	}
	else
	{
		res_was_calculated_=true;
	};

	beta0_=y1/den;

	beta0_error_=y1_error/den;


	double res=0;
	for(i=0;i<number_of_elements_;i++)
	{
		if(errors_[i]!=0)
		{
			double tmp=(beta0_+beta1_*(i+k_start_)-values_[i])/errors_[i];
			res+=tmp*tmp-c_;
		};
	};

	return res;

}

double alp_reg::median(
long int dim_,
double *array_)
{
	vector<double> array_vect(dim_);
	long int i;
	for(i=0;i<dim_;i++)
	{
		array_vect[i]=array_[i];
	};
	sort(array_vect.begin(),array_vect.end());
	if(dim_%2==0)
	{
		long int k=(long int)sls_basic::round((double)dim_/2.0);
		return 0.5*(array_vect[k-1]+array_vect[k]);
	}
	else
	{
		long int k=(long int)sls_basic::round((double)(dim_-1.0)/2.0);
		return array_vect[k];

	};
}

double alp_reg::robust_sum(
double *values,
long int dim,
long int N_points,
bool *&remove_flag)
{
	remove_flag=NULL;

	try
	{
		if(dim<=N_points)
		{
			throw error("Unexpected error\n",4);
		};

		long int i;
		remove_flag=new bool[dim];
		sls_basic::assert_mem(remove_flag);
		for(i=0;i<dim;i++)
		{
			remove_flag[i]=true;
		};


		double med_val=alp_reg::median(
		dim,
		values);

		vector<pair<double,long int> > array_vect(dim);

		for(i=0;i<dim;i++)
		{
			pair<double,long int> P;
			P.first=-fabs(values[i]-med_val);
			P.second=i;
			array_vect[i]=P;
		};

		sort(array_vect.begin(),array_vect.end());

		for(i=0;i<N_points;i++)
		{
			remove_flag[array_vect[i].second]=false;
		};

		double res=0.0;

		for(i=0;i<dim;i++)
		{
			if(remove_flag[i])
			{
				res+=values[i];
			};
		};

		res/=(double)(dim-N_points);
		return res;

	}
	catch (...)
	{ 
		delete[]remove_flag;remove_flag=NULL;
		throw;
	};

}

