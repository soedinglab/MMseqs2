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

File name: sls_pvalues.cpp

Author: Sergey Sheetlin

Contents: Calculation of P-values using precalculated Gumbel parameters

******************************************************************************/

#include "sls_pvalues.hpp"
#include "sls_alp_data.hpp"
#include <iomanip>      // std::setprecision

#include "sls_normal_distr_array.hpp"


using namespace Sls;
using namespace std;

const double nat_cut_off_in_max=2.0;//nat cut-off in max used in FSC


void pvalues::get_appr_tail_prob_with_cov(
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
bool &area_is_1_flag_)
{

	//to optimize performance
	blast_=false;

	double lambda_=par_.lambda;
	double lambda_error_=par_.lambda_error;
	double k_=par_.K;
	double k_error_=par_.K_error;

	double ai_hat_=par_.a_I;
	double ai_hat_error_=par_.a_I_error;
	double bi_hat_; 
	double bi_hat_error_; 
	double alphai_hat_=par_.alpha_I;
	double alphai_hat_error_=par_.alpha_I_error;
	double betai_hat_;
	double betai_hat_error_; 

	double aj_hat_=par_.a_J;
	double aj_hat_error_=par_.a_J_error;
	double bj_hat_; 
	double bj_hat_error_; 
	double alphaj_hat_=par_.alpha_J;
	double alphaj_hat_error_=par_.alpha_J_error;
	double betaj_hat_; 
	double betaj_hat_error_; 

	double sigma_hat_=par_.sigma;
	double sigma_hat_error_=par_.sigma_error;
	double tau_hat_;
 	double tau_hat_error_;

	{
		bi_hat_=par_.b_I;
		bi_hat_error_=par_.b_I_error;
		betai_hat_=par_.beta_I;
		betai_hat_error_=par_.beta_I_error;

		bj_hat_=par_.b_J;
		bj_hat_error_=par_.b_J_error;
		betaj_hat_=par_.beta_J;
		betaj_hat_error_=par_.beta_J_error;

		tau_hat_=par_.tau;
		tau_hat_error_=par_.tau_error;
	};


	if(blast_)
	{
		alphai_hat_=0;
		alphai_hat_error_=0;
		betai_hat_=0;
		betai_hat_error_=0;

		alphaj_hat_=0;
		alphaj_hat_error_=0;
		betaj_hat_=0;
		betaj_hat_error_=0;

		sigma_hat_=0;
		sigma_hat_error_=0;
		tau_hat_=0;
		tau_hat_error_=0;
	};

	double m_li_y_error=0;
	double m_li_y=0;

	double tmp=ai_hat_*y_+bi_hat_;

	m_li_y_error=alp_data::error_of_the_sum(fabs(y_)*ai_hat_error_,bi_hat_error_);
	m_li_y=m_-tmp;
	
	double vi_y_error=0;
	double vi_y=0;

	vi_y_error=alp_data::error_of_the_sum(fabs(y_)*alphai_hat_error_,betai_hat_error_);

	vi_y=alp_data::Tmax(par_.vi_y_thr,alphai_hat_*y_+betai_hat_);

	double sqrt_vi_y_error=alp_data::error_of_the_sqrt(vi_y,vi_y_error);

	double sqrt_vi_y=sqrt(vi_y);

	double m_F;
	double m_F_error;

	if(sqrt_vi_y==0.0||blast_)
	{
		m_F=1e100;
		m_F_error=0.0;
	}
	else
	{
		m_F_error=alp_data::error_of_the_ratio(m_li_y,m_li_y_error,sqrt_vi_y,sqrt_vi_y_error);
		m_F=m_li_y/sqrt_vi_y;
	};


	double P_m_F=sls_basic::normal_probability(m_F);
	double P_m_F_error=const_val*exp(-0.5*m_F*m_F)*m_F_error;

	double E_m_F=-const_val*exp(-0.5*m_F*m_F);
	double E_m_F_error=fabs(-E_m_F*m_F)*m_F_error;

	double m_li_y_P_m_F_error=alp_data::error_of_the_product(m_li_y,m_li_y_error,P_m_F,P_m_F_error);
	double m_li_y_P_m_F=m_li_y*P_m_F;

	double sqrt_vi_y_E_m_F_error=alp_data::error_of_the_product(sqrt_vi_y,sqrt_vi_y_error,E_m_F,E_m_F_error);
	double sqrt_vi_y_E_m_F=sqrt_vi_y*E_m_F;

	double p1_error=alp_data::error_of_the_sum(m_li_y_P_m_F_error,sqrt_vi_y_E_m_F_error);
	double p1=m_li_y_P_m_F-sqrt_vi_y_E_m_F;


	double n_lj_y_error=0;
	double n_lj_y=0;


	tmp=aj_hat_*y_+bj_hat_;

	n_lj_y_error=alp_data::error_of_the_sum(fabs(y_)*aj_hat_error_,bj_hat_error_);
	n_lj_y=n_-tmp;

	double vj_y_error=0;
	double vj_y=0;

	vj_y_error=alp_data::error_of_the_sum(fabs(y_)*alphaj_hat_error_,betaj_hat_error_);

	vj_y=alp_data::Tmax(par_.vj_y_thr,alphaj_hat_*y_+betaj_hat_);

	double sqrt_vj_y_error=alp_data::error_of_the_sqrt(vj_y,vj_y_error);

	double sqrt_vj_y=sqrt(vj_y);

	double n_F;
	double n_F_error;

	if(sqrt_vj_y==0.0||blast_)
	{
		n_F=1e100;
		n_F_error=0.0;
	}
	else
	{
		n_F_error=alp_data::error_of_the_ratio(n_lj_y,n_lj_y_error,sqrt_vj_y,sqrt_vj_y_error);

		n_F=n_lj_y/sqrt_vj_y;
	};

	double P_n_F=sls_basic::normal_probability(n_F);
	double P_n_F_error=const_val*exp(-0.5*n_F*n_F)*n_F_error;

	double E_n_F=-const_val*exp(-0.5*n_F*n_F);
	double E_n_F_error=fabs(-E_n_F*n_F)*n_F_error;

	double n_lj_y_P_n_F_error=alp_data::error_of_the_product(n_lj_y,n_lj_y_error,P_n_F,P_n_F_error);
	double n_lj_y_P_n_F=n_lj_y*P_n_F;

	double sqrt_vj_y_E_n_F_error=alp_data::error_of_the_product(sqrt_vj_y,sqrt_vj_y_error,E_n_F,E_n_F_error);
	double sqrt_vj_y_E_n_F=sqrt_vj_y*E_n_F;

	double p2_error=alp_data::error_of_the_sum(n_lj_y_P_n_F_error,sqrt_vj_y_E_n_F_error);
	double p2=n_lj_y_P_n_F-sqrt_vj_y_E_n_F;




	double c_y_error=0;
	double c_y=0;

	c_y_error=alp_data::error_of_the_sum(sigma_hat_error_*y_,tau_hat_error_);

	c_y=alp_data::Tmax(par_.c_y_thr,sigma_hat_*y_+tau_hat_);

	double P_m_F_P_n_F_error=alp_data::error_of_the_product(P_m_F,P_m_F_error,P_n_F,P_n_F_error);
	double P_m_F_P_n_F=P_m_F*P_n_F;

	double c_y_P_m_F_P_n_F_error=alp_data::error_of_the_product(c_y,c_y_error,P_m_F_P_n_F,P_m_F_P_n_F_error);
	double c_y_P_m_F_P_n_F=c_y*P_m_F_P_n_F;

	double p1_p2_error=alp_data::error_of_the_product(p1,p1_error,p2,p2_error);
	double p1_p2=p1*p2;


	double area_error=alp_data::error_of_the_sum(p1_p2_error,c_y_P_m_F_P_n_F_error);
	double area=p1_p2+c_y_P_m_F_P_n_F;




	if(!blast_)
	{
		//area=alp_data::Tmax(area,1.0);
	}
	else
	{
		if(area<=1.0)
		{
			area_is_1_flag_=true;
		};

		if(area_is_1_flag_)
		{
			area=1.0;
		};
	};


	double exp_lambda_y_error=fabs(lambda_error_*y_*exp(-lambda_*y_));
	double exp_lambda_y=exp(-lambda_*y_);

	double k_exp_lambda_y_error=alp_data::error_of_the_product(k_,k_error_,exp_lambda_y,exp_lambda_y_error);
	double k_exp_lambda_y=k_*exp_lambda_y;

	double area_k_exp_lambda_y_error=alp_data::error_of_the_product(area,area_error,k_exp_lambda_y,k_exp_lambda_y_error);
	double area_k_exp_lambda_y=-area*k_exp_lambda_y;

	E_=-area_k_exp_lambda_y;
	E_error_=area_k_exp_lambda_y_error;

	P_error_=exp(area_k_exp_lambda_y)*area_k_exp_lambda_y_error;

	P_=sls_basic::one_minus_exp_function(area_k_exp_lambda_y);
//	P_=1-exp(-k_*area*exp(-lambda_*y_));

	area_=area;


}

void pvalues::compute_intercepts(
ALP_set_of_parameters &par_)
{
	if(!par_.d_params_flag)
	{
		throw error("Unexpected error: pvalues::compute_intercepts is called for undefined parameters\n",1);
	};

	par_.b_I=2.0*par_.G*(par_.gapless_a-par_.a_I); 
	par_.b_I_error=2.0*par_.G*alp_data::error_of_the_sum(par_.gapless_a_error,par_.a_I_error); 
	par_.beta_I=2.0*par_.G*(par_.gapless_alpha-par_.alpha_I); 
	par_.beta_I_error=2.0*par_.G*alp_data::error_of_the_sum(par_.gapless_alpha_error,par_.alpha_I_error); 

	par_.b_J=2.0*par_.G*(par_.gapless_a-par_.a_J); 
	par_.b_I_error=2.0*par_.G*alp_data::error_of_the_sum(par_.gapless_a_error,par_.a_J_error); 
	par_.beta_J=2.0*par_.G*(par_.gapless_alpha-par_.alpha_J); 
	par_.beta_J_error=2.0*par_.G*alp_data::error_of_the_sum(par_.gapless_alpha_error,par_.alpha_J_error); 

	par_.tau=2.0*par_.G*(par_.gapless_alpha-par_.sigma);
	par_.tau_error=2.0*par_.G*alp_data::error_of_the_sum(par_.gapless_alpha_error,par_.sigma_error);

	long int vector_size=(long int)par_.m_LambdaSbs.size();

	par_.m_BISbs.resize(vector_size);
	par_.m_BJSbs.resize(vector_size);
	par_.m_BetaISbs.resize(vector_size);
	par_.m_BetaJSbs.resize(vector_size);
	par_.m_TauSbs.resize(vector_size);

	long int i;
	for(i=0;i<vector_size;i++)
	{
		par_.m_BISbs[i]=2.0*par_.G*(par_.gapless_a-par_.m_AISbs[i]); 
		par_.m_BetaISbs[i]=2.0*par_.G*(par_.gapless_alpha-par_.m_AlphaISbs[i]); 

		par_.m_BJSbs[i]=2.0*par_.G*(par_.gapless_a-par_.m_AJSbs[i]); 
		par_.m_BetaJSbs[i]=2.0*par_.G*(par_.gapless_alpha-par_.m_AlphaJSbs[i]); 

		par_.m_TauSbs[i]=2.0*par_.G*(par_.gapless_alpha-par_.m_SigmaSbs[i]);
	};

	compute_tmp_values(par_);

}

void pvalues::compute_tmp_values(ALP_set_of_parameters &par_)
{
	if(!par_.d_params_flag)
	{
		throw error("Unexpected call of pvalues::compute_tmp_values\n",1);
	};

	//tmp values
	if(par_.lambda>0)
	{
		par_.vi_y_thr=alp_data::Tmax(nat_cut_off_in_max*par_.alpha_I/par_.lambda,0.0);
		par_.vj_y_thr=alp_data::Tmax(nat_cut_off_in_max*par_.alpha_J/par_.lambda,0.0);
		par_.c_y_thr=alp_data::Tmax(nat_cut_off_in_max*par_.sigma/par_.lambda,0.0);
	}
	else
	{
		par_.vi_y_thr=0;
		par_.vj_y_thr=0;
		par_.c_y_thr=0;

		par_.d_params_flag=false;
	};
}

void pvalues::get_appr_tail_prob_with_cov_without_errors(
const ALP_set_of_parameters &par_,
bool blast_,
double y_,
double m_,
double n_,

double &P_,

double &E_,

double &area_,
bool &area_is_1_flag_,
bool compute_only_area_)
{

	//to optimize performance
	blast_=false;

	double lambda_=par_.lambda;
	double k_=par_.K;

	double ai_hat_=par_.a_I;
	double bi_hat_;
	double alphai_hat_=par_.alpha_I;
	double betai_hat_;

	double aj_hat_=par_.a_J;
	double bj_hat_;
	double alphaj_hat_=par_.alpha_J;
	double betaj_hat_;

	double sigma_hat_=par_.sigma;
	double tau_hat_;

	{
		bi_hat_=par_.b_I;
		betai_hat_=par_.beta_I;

		bj_hat_=par_.b_J;
		betaj_hat_=par_.beta_J;

		tau_hat_=par_.tau;
	};

	if(blast_)
	{
		alphai_hat_=0;
		betai_hat_=0;

		alphaj_hat_=0;
		betaj_hat_=0;

		sigma_hat_=0;
		tau_hat_=0;
	};

	double m_li_y=0;

	double tmp=ai_hat_*y_+bi_hat_;

	m_li_y=m_-tmp;
	
	double vi_y=0;

	vi_y=alp_data::Tmax(par_.vi_y_thr,alphai_hat_*y_+betai_hat_);

	double sqrt_vi_y=sqrt(vi_y);


	double m_F;

	if(sqrt_vi_y==0.0||blast_)
	{
		m_F=1e100;
	}
	else
	{
		m_F=m_li_y/sqrt_vi_y;
	};


	double P_m_F=sls_basic::normal_probability(m_F);

	double E_m_F=-const_val*exp(-0.5*m_F*m_F);

	double m_li_y_P_m_F=m_li_y*P_m_F;

	double sqrt_vi_y_E_m_F=sqrt_vi_y*E_m_F;

	double p1=m_li_y_P_m_F-sqrt_vi_y_E_m_F;


	double n_lj_y=0;

	tmp=aj_hat_*y_+bj_hat_;

	n_lj_y=n_-tmp;

	double vj_y=0;

	vj_y=alp_data::Tmax(par_.vj_y_thr,alphaj_hat_*y_+betaj_hat_);

	double sqrt_vj_y=sqrt(vj_y);

	double n_F;

	if(sqrt_vj_y==0.0||blast_)
	{
		n_F=1e100;
	}
	else
	{
		n_F=n_lj_y/sqrt_vj_y;
	};

	double P_n_F=sls_basic::normal_probability(n_F);

	double E_n_F=-const_val*exp(-0.5*n_F*n_F);

	double n_lj_y_P_n_F=n_lj_y*P_n_F;

	double sqrt_vj_y_E_n_F=sqrt_vj_y*E_n_F;

	double p2=n_lj_y_P_n_F-sqrt_vj_y_E_n_F;




	double c_y=0;

	c_y=alp_data::Tmax(par_.c_y_thr,sigma_hat_*y_+tau_hat_);

	double P_m_F_P_n_F=P_m_F*P_n_F;

	double c_y_P_m_F_P_n_F=c_y*P_m_F_P_n_F;

	double p1_p2=p1*p2;

	double area=p1_p2+c_y_P_m_F_P_n_F;




	if(!blast_)
	{
		//area=alp_data::Tmax(area,1.0);
	}
	else
	{
		if(area<=1.0)
		{
			area_is_1_flag_=true;
		};

		if(area_is_1_flag_)
		{
			area=1.0;
		};
	};

	area_=area;

	if(compute_only_area_)
	{
		return;
	};

	double exp_lambda_y=exp(-lambda_*y_);

	double k_exp_lambda_y=k_*exp_lambda_y;

	double area_k_exp_lambda_y=-area*k_exp_lambda_y;

	E_=-area_k_exp_lambda_y;

	P_=sls_basic::one_minus_exp_function(area_k_exp_lambda_y);
//	P_=1-exp(-k_*area*exp(-lambda_*y_));

}

void pvalues::get_P_error_using_splitting_method(
const ALP_set_of_parameters &par_,
bool blast_,
double y_,
double m_,
double n_,

double &P_,
double &P_error_,

double &E_,
double &E_error_,

bool &area_is_1_flag_)
{
	long int dim=par_.m_LambdaSbs.size();
	if(dim==0)
	{
		throw error("Unexpected error in get_P_error_using_splitting_method\n",1);
	};

	P_=0;
	P_error_=0;

	E_=0;
	E_error_=0;

	double exp_E_values_aver=0;
	double exp_E_values_error=0;


	vector<double> P_values(dim);
	vector<double> E_values(dim);
	vector<double> exp_E_values(dim);


	long int i;
	for(i=0;i<dim;i++)
	{
		ALP_set_of_parameters par_tmp;

		par_tmp.a_I=par_.m_AISbs[i];
		par_tmp.a_I_error=0;

		par_tmp.a_J=par_.m_AJSbs[i];
		par_tmp.a_J_error=0;

		

		par_tmp.gapless_a=par_.gapless_a;
		par_tmp.gapless_a_error=par_.gapless_a_error;

		par_tmp.a=0.5*(par_tmp.a_I+par_tmp.a_J);
		par_tmp.a_error=0;


		par_tmp.sigma=par_.m_SigmaSbs[i];
		par_tmp.sigma_error=0;

		par_tmp.gapless_alpha=par_.gapless_alpha;
		par_tmp.gapless_alpha_error=par_.gapless_alpha_error;


		par_tmp.C=par_.m_CSbs[i];
		par_tmp.C_error=0;

		par_tmp.K=par_.m_KSbs[i];
		par_tmp.K_error=0;


		par_tmp.lambda=par_.m_LambdaSbs[i];
		par_tmp.lambda_error=0;

		par_tmp.alpha_I=par_.m_AlphaISbs[i];
		par_tmp.alpha_I_error=0;

		par_tmp.alpha_J=par_.m_AlphaJSbs[i];
		par_tmp.alpha_J_error=0;

		par_tmp.alpha=0.5*(par_tmp.alpha_I+par_tmp.alpha_J);
		par_tmp.alpha_error=0;

		par_tmp.G=par_.G;
		par_tmp.G1=par_.G1;
		par_tmp.G2=par_.G2;

		//intercepts

		{
			par_tmp.b_I=par_.m_BISbs[i];
			par_tmp.b_I_error=0;

			par_tmp.b_J=par_.m_BJSbs[i];
			par_tmp.b_J_error=0;

			par_tmp.beta_I=par_.m_BetaISbs[i];
			par_tmp.beta_I_error=0;

			par_tmp.beta_J=par_.m_BetaJSbs[i];
			par_tmp.beta_J_error=0;

			par_tmp.tau=par_.m_TauSbs[i];
			par_tmp.tau_error=0;

			par_tmp.d_params_flag=true;

			compute_tmp_values(par_tmp);

		};


		double P_tmp,area_tmp,E_tmp;

		get_appr_tail_prob_with_cov_without_errors(
		par_tmp,
		blast_,
		y_,
		m_,
		n_,

		P_tmp,

		E_tmp,

		area_tmp,
		area_is_1_flag_);

		P_values[i]=P_tmp;

		P_+=P_tmp;

		E_values[i]=E_tmp;

		E_+=E_tmp;

		double exp_E_tmp=exp(-E_tmp);
		exp_E_values[i]=exp_E_tmp;
		exp_E_values_aver+=exp_E_tmp;


	};

	if(dim<=1)
	{
		return;
	};


	if(P_<=0)
	{
		return;
	};

	if(E_<=0)
	{
		return;
	};


	P_/=(double)dim;
	E_/=(double)dim;
	exp_E_values_aver/=(double)dim;

	for(i=0;i<dim;i++)
	{
		double tmp;
		
		if(P_>0)
		{
			tmp=P_values[i]/P_;
			P_error_+=tmp*tmp;
		};

		if(E_>0)
		{
			tmp=E_values[i]/E_;
			E_error_+=tmp*tmp;
		};

		if(exp_E_values_aver>0)
		{
			tmp=exp_E_values[i]/exp_E_values_aver;
			exp_E_values_error+=tmp*tmp;
		};

	};

	P_error_/=(double)dim;
	P_error_-=1;
	
	E_error_/=(double)dim;
	E_error_-=1;

	exp_E_values_error/=(double)dim;
	exp_E_values_error-=1;


	if(P_<1e-4)
	{
		P_error_=P_*alp_reg::sqrt_for_errors(P_error_/(double)dim);
	}
	else
	{
		P_error_=exp_E_values_aver*alp_reg::sqrt_for_errors(exp_E_values_error/(double)dim);
	};

	E_error_=E_*alp_reg::sqrt_for_errors(E_error_/(double)dim);

}


pvalues::pvalues()
{
	blast=false;
	eps=0.0001;
	a_normal=-10;
	b_normal=10;
	N_normal=NORMAL_DISTR_ARRAY_DIM;
	h_normal=(b_normal-a_normal)/(double)N_normal;
	p_normal=normal_distr_array_for_P_values_calculation;
}


pvalues::~pvalues()
{
	
}

void pvalues::calculate_P_values(
long int Score1,
long int Score2,
double Seq1Len,
double Seq2Len,
const ALP_set_of_parameters &ParametersSet,
vector<double> &P_values,
vector<double> &P_values_errors,
vector<double> &E_values,
vector<double> &E_values_errors)
{
	if(Score2<Score1)
	{
		throw error("Error - Score2<Score1\n",2);
	};

	if(Seq1Len<=0||Seq2Len<=0)
	{
		throw error("Error - Seq1Len<=0||Seq2Len<=0\n",2);
	};

	P_values.resize(Score2-Score1+1);
	P_values_errors.resize(Score2-Score1+1);

	E_values.resize(Score2-Score1+1);
	E_values_errors.resize(Score2-Score1+1);


	long int y;
	for(y=Score1;y<=Score2;y++)
	{
		calculate_P_values(
		y,
		Seq1Len,
		Seq2Len,
		ParametersSet,
		P_values[y-Score1],
		P_values_errors[y-Score1],
		E_values[y-Score1],
		E_values_errors[y-Score1]);
	};

}

void pvalues::calculate_P_values(
double Score,
double Seq1Len,
double Seq2Len,
const ALP_set_of_parameters &ParametersSet,
double &P_value,
double &P_value_error,
double &E_value,
double &E_value_error,
bool read_Sbs_par_flag)
{

	if(Seq1Len<=0||Seq2Len<=0)
	{
		throw error("Error - Seq1Len<=0||Seq2Len<=0\n",2);
	};

	double P;
	double P_error;
	double E;
	double E_error;
	double area;
	bool area_is_1_flag=false;


	if(read_Sbs_par_flag)
	{
		

		get_appr_tail_prob_with_cov_without_errors(
		ParametersSet,
		blast,
		Score,
		Seq1Len,
		Seq2Len,

		P,

		E,

		area,
		area_is_1_flag);


		
		double P_tmp,E_tmp;

		if(ParametersSet.m_LambdaSbs.size()>0)
		{
			get_P_error_using_splitting_method(
			ParametersSet,
			blast,
			Score,
			Seq1Len,
			Seq2Len,

			P_tmp,
			P_error,

			E_tmp,
			E_error,

			area_is_1_flag);


			if(P_tmp>0)
			{
				P_error=P_error/P_tmp*P;
			};

			P_value_error=P_error;

			if(E_tmp>0)
			{
				E_error=E_error/E_tmp*E;
			};

			E_value_error=E_error;

		}
		else
		{
			P_value_error=-DBL_MAX;
			E_value_error=-DBL_MAX;
		};

		
		
	}
	else
	{
		get_appr_tail_prob_with_cov(
		ParametersSet,
		blast,
		Score,
		Seq1Len,
		Seq2Len,

		P,
		P_error,

		E,
		E_error,

		area,
		area_is_1_flag);

		P_value_error=P_error;
		E_value_error=E_error;
	};

	P_value=P;
	E_value=E;
}


//input/output Gumbel parameters

namespace Sls {

std::ostream &operator<<(std::ostream &s_,
const ALP_set_of_parameters &gumbel_params_)
{


	s_<<"Lambda\tLambda error\tK\tK error\tC\tC error\ta\ta error\ta_1\ta_1 error\ta_2\ta_2 error\tsigma\tsigma error\talpha\talpha error\talpha_1\talpha_1 error\talpha_2\talpha_2 error\tGapless a\tGapless a error\tGapless alpha\tGapless alpha error\tG\tCalculation time\tArrays for error calculation\n";
	s_.precision(17);
	s_<<
		gumbel_params_.lambda<<"\t"<<gumbel_params_.lambda_error<<"\t"<<
		gumbel_params_.K<<"\t"<<gumbel_params_.K_error<<"\t"<<
		gumbel_params_.C<<"\t"<<gumbel_params_.C_error<<"\t"<<
		gumbel_params_.a<<"\t"<<gumbel_params_.a_error<<"\t"<<
		gumbel_params_.a_J<<"\t"<<gumbel_params_.a_J_error<<"\t"<<
		gumbel_params_.a_I<<"\t"<<gumbel_params_.a_I_error<<"\t"<<
		gumbel_params_.sigma<<"\t"<<gumbel_params_.sigma_error<<"\t"<<
		gumbel_params_.alpha<<"\t"<<gumbel_params_.alpha_error<<"\t"<<
		gumbel_params_.alpha_J<<"\t"<<gumbel_params_.alpha_J_error<<"\t"<<
		gumbel_params_.alpha_I<<"\t"<<gumbel_params_.alpha_I_error<<"\t"<<
		gumbel_params_.gapless_a<<"\t"<<gumbel_params_.gapless_a_error<<"\t"<<
		gumbel_params_.gapless_alpha<<"\t"<<gumbel_params_.gapless_alpha_error<<"\t"<<
		gumbel_params_.G<<"\t"<<
		gumbel_params_.m_CalcTime<<"\t";

	long int i;


	{
		const vector<double> &tmp=gumbel_params_.m_LambdaSbs;
		s_<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			s_<<tmp[i]<<"\t";
		};
	};

	{
		const vector<double> &tmp=gumbel_params_.m_KSbs;
		s_<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			s_<<tmp[i]<<"\t";
		};
	};

	{
		const vector<double> &tmp=gumbel_params_.m_CSbs;
		s_<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			s_<<tmp[i]<<"\t";
		};
	};

	{
		const vector<double> &tmp=gumbel_params_.m_AJSbs;
		s_<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			s_<<tmp[i]<<"\t";
		};
	};


	{
		const vector<double> &tmp=gumbel_params_.m_AISbs;
		s_<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			s_<<tmp[i]<<"\t";
		};
	};


	{
		const vector<double> &tmp=gumbel_params_.m_SigmaSbs;
		s_<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			s_<<tmp[i]<<"\t";
		};
	};

	{
		const vector<double> &tmp=gumbel_params_.m_AlphaJSbs;
		s_<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			s_<<tmp[i]<<"\t";
		};
	};


	{
		const vector<double> &tmp=gumbel_params_.m_AlphaISbs;
		s_<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			s_<<tmp[i]<<"\t";
		};
	};

	s_<<endl;

	return s_;

}

std::istream &operator>>(std::istream &s_,
ALP_set_of_parameters &gumbel_params_)
{

	gumbel_params_.d_params_flag=false;
	try
	{

		string st;
		getline(s_,st);
		s_>>
			gumbel_params_.lambda>>gumbel_params_.lambda_error>>
			gumbel_params_.K>>gumbel_params_.K_error>>
			gumbel_params_.C>>gumbel_params_.C_error>>
			gumbel_params_.a>>gumbel_params_.a_error>>
			gumbel_params_.a_J>>gumbel_params_.a_J_error>>
			gumbel_params_.a_I>>gumbel_params_.a_I_error>>
			gumbel_params_.sigma>>gumbel_params_.sigma_error>>
			gumbel_params_.alpha>>gumbel_params_.alpha_error>>
			gumbel_params_.alpha_J>>gumbel_params_.alpha_J_error>>
			gumbel_params_.alpha_I>>gumbel_params_.alpha_I_error>>
			gumbel_params_.gapless_a>>gumbel_params_.gapless_a_error>>
			gumbel_params_.gapless_alpha>>gumbel_params_.gapless_alpha_error>>
			gumbel_params_.G>>
			gumbel_params_.m_CalcTime;

		long int i;


		{
			vector<double> &tmp=gumbel_params_.m_LambdaSbs;
			long int tmp_size;
			s_>>tmp_size;
			if(tmp_size<=0)
			{
				throw error("Error in the input parameters\n",4);
			};
			tmp.resize(tmp_size);
			for(i=0;i<tmp_size;i++)
			{
				s_>>tmp[i];
			};
		};

		{
			vector<double> &tmp=gumbel_params_.m_KSbs;
			long int tmp_size;
			s_>>tmp_size;
			if(tmp_size<=0)
			{
				throw error("Error in the input parameters\n",4);
			};
			tmp.resize(tmp_size);
			for(i=0;i<tmp_size;i++)
			{
				s_>>tmp[i];
			};
		};


		{
			vector<double> &tmp=gumbel_params_.m_CSbs;
			long int tmp_size;
			s_>>tmp_size;
			if(tmp_size<=0)
			{
				throw error("Error in the input parameters\n",4);
			};
			tmp.resize(tmp_size);
			for(i=0;i<tmp_size;i++)
			{
				s_>>tmp[i];
			};
		};



		{
			vector<double> &tmp=gumbel_params_.m_AJSbs;
			long int tmp_size;
			s_>>tmp_size;
			if(tmp_size<=0)
			{
				throw error("Error in the input parameters\n",4);
			};
			tmp.resize(tmp_size);
			for(i=0;i<tmp_size;i++)
			{
				s_>>tmp[i];
			};
		};

		{
			vector<double> &tmp=gumbel_params_.m_AISbs;
			long int tmp_size;
			s_>>tmp_size;
			if(tmp_size<=0)
			{
				throw error("Error in the input parameters\n",4);
			};
			tmp.resize(tmp_size);
			for(i=0;i<tmp_size;i++)
			{
				s_>>tmp[i];
			};
		};



		{
			vector<double> &tmp=gumbel_params_.m_SigmaSbs;
			long int tmp_size;
			s_>>tmp_size;
			if(tmp_size<=0)
			{
				throw error("Error in the input parameters\n",4);
			};
			tmp.resize(tmp_size);
			for(i=0;i<tmp_size;i++)
			{
				s_>>tmp[i];
			};
		};


		{
			vector<double> &tmp=gumbel_params_.m_AlphaJSbs;
			long int tmp_size;
			s_>>tmp_size;
			if(tmp_size<=0)
			{
				throw error("Error in the input parameters\n",4);
			};
			tmp.resize(tmp_size);
			for(i=0;i<tmp_size;i++)
			{
				s_>>tmp[i];
			};
		};

		{
			vector<double> &tmp=gumbel_params_.m_AlphaISbs;
			long int tmp_size;
			s_>>tmp_size;
			if(tmp_size<=0)
			{
				throw error("Error in the input parameters\n",4);
			};
			tmp.resize(tmp_size);
			for(i=0;i<tmp_size;i++)
			{
				s_>>tmp[i];
			};
		};

		gumbel_params_.d_params_flag=true;
		return s_;
	}
	catch (...)
	{ 
		gumbel_params_.d_params_flag=false;
		throw;
	};

}

//returns "true" if the Gumbel parameters are properly defined and "false" otherwise
bool pvalues::assert_Gumbel_parameters(
const ALP_set_of_parameters &par_)//a set of Gumbel parameters
{
		if(!(par_.lambda>0)||
		par_.lambda_error<0||

		//the parameters C and K_C are not necessary for the P-value calculation
		//par_.C<0||
		//par_.C_error<0||

		!(par_.K>0)||
		par_.K_error<0||

		par_.a_I<0||
		par_.a_I_error<0||

		par_.a_J<0||
		par_.a_J_error<0||

		par_.sigma<0||
		par_.sigma_error<0||

		par_.alpha_I<0||
		par_.alpha_I_error<0||

		par_.alpha_J<0||
		par_.alpha_J_error<0||

		par_.gapless_a<0||
		par_.gapless_a_error<0||

		par_.gapless_alpha<0||
		par_.gapless_alpha_error<0||

		par_.G<0||
		par_.G1<0||
		par_.G2<0||

		//intercepts
		par_.b_I_error<0||

		par_.b_J_error<0||

		par_.beta_I_error<0||

		par_.beta_J_error<0||

		par_.tau_error<0

		)
		{
			return false;
		};



		size_t size_tmp=par_.m_LambdaSbs.size();
		if(
		par_.m_KSbs.size()!=size_tmp||
		//par_.m_CSbs.size()!=size_tmp||

		par_.m_SigmaSbs.size()!=size_tmp||

		par_.m_AlphaISbs.size()!=size_tmp||
		par_.m_AlphaJSbs.size()!=size_tmp||

		par_.m_AISbs.size()!=size_tmp||
		par_.m_AJSbs.size()!=size_tmp||

		par_.m_BISbs.size()!=size_tmp||
		par_.m_BJSbs.size()!=size_tmp||

		par_.m_BetaISbs.size()!=size_tmp||
		par_.m_BetaJSbs.size()!=size_tmp||

		par_.m_TauSbs.size()!=size_tmp)
		{
			return false;
		};


		return true;

}



}

