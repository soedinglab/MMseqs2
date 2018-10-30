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

File name: sls_alp_data.cpp

Author: Sergey Sheetlin, Martin Frith

Contents: Input data for the ascending ladder points simulation

******************************************************************************/

#include "sls_alp_data.hpp"

using namespace Sls;
using namespace std;

void alp_data::input_data_for_the_constructor(
string randout_,//if defined, then the program outputs complete randomization information into a file
string smatr_file_name_,//scoring matrix file name
string RR1_file_name_,//probabilities1 file name
string RR2_file_name_,//probabilities2 file name

struct_for_randomization &rand_all_,
bool &rand_flag_,
long int &rand_,

long int &alphabetSize_,
long int **&substitutionScoreMatrix_,
double *&letterFreqs1_,
double *&letterFreqs2_)
{
	ifstream frand;

	try
	{

		long int number_of_AA_RR1;
		long int number_of_AA_RR2;
		long int number_of_AA_smatr;

		read_smatr(
		smatr_file_name_,
		d_smatr,
		number_of_AA_smatr);

		d_number_of_AA_smatr=number_of_AA_smatr;
		
		read_RR(
		RR1_file_name_,
		d_RR1,
		d_RR1_sum,
		d_RR1_sum_elements,
		number_of_AA_RR1);


		read_RR(
		RR2_file_name_,
		d_RR2,
		d_RR2_sum,
		d_RR2_sum_elements,
		number_of_AA_RR2);


		if(number_of_AA_RR1==number_of_AA_smatr)
		{
			alphabetSize_=number_of_AA_smatr;
		}
		else
		{
			throw error("Number of letters is different in the files "+smatr_file_name_+" and "+RR1_file_name_+"\n",3);
		};

		if(number_of_AA_RR2!=number_of_AA_smatr)
		{
			throw error("Number of letters is different in the files "+smatr_file_name_+" and "+RR2_file_name_+"\n",3);
		};


		if(randout_!="")
		{
			rand_flag_=true;

			string rand_st=randout_;
			frand.open(rand_st.data(),ios::in);
			if(!frand)
			{
				rand_flag_=false;
			}
			else
			{

				long int i,size;
				frand>>rand_all_.d_random_seed;


				if(rand_all_.d_random_seed<0)
				{
					throw error("File "+rand_st+" is not correct\n",3);
				};

				rand_=rand_all_.d_random_seed;



				frand>>size;
				for(i=0;i<size;i++)
				{
					long int tmp;
					frand>>tmp;
					rand_all_.d_first_stage_preliminary_realizations_numbers_ALP.push_back(tmp);
					if(tmp<0)
					{
						throw error("File "+rand_st+" is not correct\n",3);
					};

				};

				frand>>size;
				for(i=0;i<size;i++)
				{
					long int tmp;
					frand>>tmp;
					rand_all_.d_preliminary_realizations_numbers_ALP.push_back(tmp);
					if(tmp<0)
					{
						throw error("File "+rand_st+" is not correct\n",3);
					};

				};

				frand>>size;
				for(i=0;i<size;i++)
				{
					long int tmp;
					frand>>tmp;
					rand_all_.d_preliminary_realizations_numbers_killing.push_back(tmp);
					if(tmp<0)
					{
						throw error("File "+rand_st+" is not correct\n",3);
					};

				};


				frand>>rand_all_.d_total_realizations_number_with_ALP;
				if(rand_all_.d_total_realizations_number_with_ALP<0)
				{
					throw error("File "+rand_st+" is not correct\n",3);
				};

				frand>>rand_all_.d_total_realizations_number_with_killing;
				if(rand_all_.d_total_realizations_number_with_killing<0)
				{
					throw error("File "+rand_st+" is not correct\n",3);
				};

				frand.close();
			};
		}
		else
		{
			rand_flag_=false;
		};


		if(!(alphabetSize_>0&&substitutionScoreMatrix_&&letterFreqs1_&&letterFreqs2_))
		{
			throw error("Incorrect input parameters\n",1);
		};

	}
	catch (...)
	{ 
		if(frand.is_open())
		{
			frand.close();
		};
		throw;
	};

}

void alp_data::init_main_class_members(
long int rand_,//randomization number
string randout_,//if defined, then the program outputs complete randomization information into a file

long int open_,//gap opening penalty
long int open1_,//gap opening penalty for a gap in the sequence #1
long int open2_,//gap opening penalty for a gap in the sequence #2

long int epen_,//gap extension penalty
long int epen1_,//gap extension penalty for a gap in the sequence #1
long int epen2_,//gap extension penalty for a gap in the sequence #2

double temperature_,
double max_time_,//maximum allowed calculation time in seconds
double max_mem_,//maximum allowed memory usage in MB
double eps_lambda_,//relative error for lambda calculation
double eps_K_,//relative error for K calculation
bool insertions_after_deletions_)//if true, then insertions after deletions are allowed
{
	try
	{
		d_randout=randout_;

		if(!d_rand_flag&&rand_<0)
		{
			rand_=sls_basic::random_seed_from_time();
			d_rand_flag=false;
			
		};


		d_random_seed=rand_;
		//cout<<"Random seed "<<d_random_seed<<endl;

		//srand(d_random_seed);

		Njn::Random::seed(d_random_seed);

		d_number_of_AA_smatr=d_number_of_AA;

		d_sentinels_flag=false;

		d_memory_size_in_MB=0;

		#ifndef _MSDOS_ //UNIX program

		#else
			_CrtMemCheckpoint( &d_s1 );
		#endif

		d_smatr_symmetric_flag=false;
		
		long int t;
		for(t=0;t<d_number_of_AA;t++)
		{
			if(d_RR1[t]!=d_RR2[t])
			{
				d_smatr_symmetric_flag=false;
				break;
			};
		};

		d_insertions_after_deletions=insertions_after_deletions_;



		d_open=open_+epen_;
		d_open1=open1_+epen1_;
		d_open2=open2_+epen2_;

		d_epen=epen_;
		d_epen1=epen1_;
		d_epen2=epen2_;

		d_max_time=max_time_;
		d_max_mem=max_mem_;
		d_eps_lambda=eps_lambda_;
		d_eps_K=eps_K_;
		d_minimum_realizations_number=40;



		


		d_is=new importance_sampling(
			this,
			d_open,

			d_epen,

			temperature_,
			d_number_of_AA,
			d_smatr,
			d_RR1,
			d_RR2);

		alp_data::assert_mem(d_is);

		d_memory_size_in_MB+=sizeof(d_is)/mb_bytes;

		d_r_i_dot=new double[d_number_of_AA];
		alp_data::assert_mem(d_r_i_dot);
		d_r_dot_j=new double[d_number_of_AA];
		alp_data::assert_mem(d_r_dot_j);
		long int k;
		for(k=0;k<d_number_of_AA;k++)
		{
			d_r_i_dot[k]=0;
			if(d_RR1[k]!=0)
			{
				long int i;
				for(i=0;i<d_number_of_AA;i++)
				{
					if(d_RR2[i]!=0)
					{
						d_r_i_dot[k]+=d_is->d_exp_s[k][i]*d_RR2[i];
					};
				};
			};
		};

		for(k=0;k<d_number_of_AA;k++)
		{
			d_r_dot_j[k]=0;
			if(d_RR2[k]!=0)
			{
				long int i;
				for(i=0;i<d_number_of_AA;i++)
				{
					if(d_RR1[i]!=0)
					{
						d_r_dot_j[k]+=d_is->d_exp_s[i][k]*d_RR1[i];
					};
				};
			};
		};


		d_memory_size_in_MB+=(double)(sizeof(double)*d_number_of_AA*2.0)/mb_bytes;

		double tmp_size1=LONG_MAX;

		double tmp_size=Tmin((double)(tmp_size1),
			(

			(double)mb_bytes*d_max_mem/(double)d_minimum_realizations_number
			)
			/(
			(double)(sizeof(double)*12)+(double)(sizeof(long int)*17)
			)
			);

		d_dim1_tmp=(long int)tmp_size;
		d_dim2_tmp=(long int)tmp_size;

	}
	catch (...)
	{ 
		release_memory();
		throw;
	};
}


alp_data::alp_data(//constructor
long int rand_,//randomization number
string randout_,//if defined, then the program outputs complete randomization information into a file

long int open_,//gap opening penalty
long int open1_,//gap opening penalty for a gap in the sequence #1
long int open2_,//gap opening penalty for a gap in the sequence #2

long int epen_,//gap extension penalty
long int epen1_,//gap extension penalty for a gap in the sequence #1
long int epen2_,//gap extension penalty for a gap in the sequence #2

string smatr_file_name_,//scoring matrix file name
string RR1_file_name_,//probabilities1 file name
string RR2_file_name_,//probabilities2 file name
double temperature_,
double max_time_,//maximum allowed calculation time in seconds
double max_mem_,//maximum allowed memory usage in MB
double eps_lambda_,//relative error for lambda calculation
double eps_K_,//relative error for K calculation
bool insertions_after_deletions_)//if true, then insertions after deletions are allowed
{


	d_smatr=NULL;
	d_RR1=NULL;
	d_RR1_sum=NULL;
	d_RR1_sum_elements=NULL;

	d_RR2=NULL;
	d_RR2_sum=NULL;
	d_RR2_sum_elements=NULL;

	d_is=NULL;
	d_r_i_dot=NULL;
	d_r_dot_j=NULL;

	d_rand_all=NULL;

	try
	{

		d_rand_flag=true;

		d_rand_all=new struct_for_randomization;
		alp_data::assert_mem(d_rand_all);
		d_memory_size_in_MB+=sizeof(struct_for_randomization)/mb_bytes;


		input_data_for_the_constructor(
		randout_,//if defined, then the program outputs complete randomization information into a file
		smatr_file_name_,//scoring matrix file name
		RR1_file_name_,//probabilities1 file name
		RR2_file_name_,//probabilities2 file name

		*d_rand_all,
		d_rand_flag,
		rand_,


		d_number_of_AA,
		d_smatr,
		d_RR1,
		d_RR2);

		init_main_class_members(
		rand_,//randomization number
		randout_,//if defined, then the program outputs complete randomization information into a file

		open_,//gap opening penalty
		open1_,//gap opening penalty for a gap in the sequence #1
		open2_,//gap opening penalty for a gap in the sequence #2

		epen_,//gap extension penalty
		epen1_,//gap extension penalty for a gap in the sequence #1
		epen2_,//gap extension penalty for a gap in the sequence #2

		temperature_,
		max_time_,//maximum allowed calculation time in seconds
		max_mem_,//maximum allowed memory usage in MB
		eps_lambda_,//relative error for lambda calculation
		eps_K_,//relative error for K calculation
		insertions_after_deletions_);//if true, then insertions after deletions are allowed

		d_max_time_with_computation_parameters=-1;

		if(max_time_>0)
		{
			d_max_time_for_quick_tests=0.25*max_time_;
		}
		else
		{
			d_max_time_for_quick_tests=1e99;
		};



	}
	catch (...)
	{
		release_memory();
		throw;
	};
}


alp_data::alp_data(//constructor
long int rand_,//randomization number
struct_for_randomization *randomization_parameters_,//if not NULL, sets d_rand_flag to true and initializes d_rand_all

long int open_,//gap opening penalty
long int open1_,//gap opening penalty for a gap in the sequence #1
long int open2_,//gap opening penalty for a gap in the sequence #2

long int epen_,//gap extension penalty
long int epen1_,//gap extension penalty for a gap in the sequence #1
long int epen2_,//gap extension penalty for a gap in the sequence #2

long alphabetSize_,
const long *const *substitutionScoreMatrix_,
const double *letterFreqs1_,
const double *letterFreqs2_,

double temperature_,
double max_time_,//maximum allowed calculation time in seconds
double max_mem_,//maximum allowed memory usage in MB
double eps_lambda_,//relative error for lambda calculation
double eps_K_,//relative error for K calculation
bool insertions_after_deletions_,//if true, then insertions after deletions are allowed
double max_time_for_quick_tests_,//maximum allowed calculation time in seconds for quick tests
double max_time_with_computation_parameters_)//maximum allowed time in seconds for the whole computation
{


	d_smatr=NULL;
	d_RR1=NULL;
	d_RR1_sum=NULL;
	d_RR1_sum_elements=NULL;

	d_RR2=NULL;
	d_RR2_sum=NULL;
	d_RR2_sum_elements=NULL;

	d_is=NULL;
	d_r_i_dot=NULL;
	d_r_dot_j=NULL;

	d_rand_all=NULL;

	try
	{

		d_rand_flag=false;

		d_number_of_AA=alphabetSize_;
		string randout="";


		get_memory_for_matrix(alphabetSize_,d_smatr);
		alp_data::assert_mem(d_smatr);
		d_RR1=new double[alphabetSize_];
		alp_data::assert_mem(d_RR1);
		d_RR2=new double[alphabetSize_];
		alp_data::assert_mem(d_RR2);

		long int i,j;
		for(i=0;i<alphabetSize_;i++)
		{
			d_RR1[i]=letterFreqs1_[i];
			d_RR2[i]=letterFreqs2_[i];
			for(j=0;j<alphabetSize_;j++)
			{
				d_smatr[i][j]=substitutionScoreMatrix_[i][j];
			};

		};

		d_rand_all=new struct_for_randomization;

		if(randomization_parameters_)
		{
			d_rand_flag=true;

			d_rand_all->d_first_stage_preliminary_realizations_numbers_ALP=
				randomization_parameters_->d_first_stage_preliminary_realizations_numbers_ALP;

			d_rand_all->d_preliminary_realizations_numbers_ALP=
				randomization_parameters_->d_preliminary_realizations_numbers_ALP;

			d_rand_all->d_preliminary_realizations_numbers_killing=
				randomization_parameters_->d_preliminary_realizations_numbers_killing;

			d_rand_all->d_random_seed=
				randomization_parameters_->d_random_seed;

			d_rand_all->d_total_realizations_number_with_ALP=
				randomization_parameters_->d_total_realizations_number_with_ALP;

			d_rand_all->d_total_realizations_number_with_killing=
				randomization_parameters_->d_total_realizations_number_with_killing;
		};

		alp_data::assert_mem(d_rand_all);
		d_memory_size_in_MB+=sizeof(struct_for_randomization)/mb_bytes;

		init_main_class_members(
		rand_,//randomization number
		randout,//if defined, then the program outputs complete randomization information into a file

		open_,//gap opening penalty
		open1_,//gap opening penalty for a gap in the sequence #1
		open2_,//gap opening penalty for a gap in the sequence #2

		epen_,//gap extension penalty
		epen1_,//gap extension penalty for a gap in the sequence #1
		epen2_,//gap extension penalty for a gap in the sequence #2

		temperature_,
		max_time_,//maximum allowed calculation time in seconds
		max_mem_,//maximum allowed memory usage in MB
		eps_lambda_,//relative error for lambda calculation
		eps_K_,//relative error for K calculation
		insertions_after_deletions_);//if true, then insertions after deletions are allowed

		if(max_time_for_quick_tests_>0)
		{
			d_max_time_for_quick_tests=max_time_for_quick_tests_;
		}
		else
		{
			if(max_time_>0)
			{
				d_max_time_for_quick_tests=0.25*max_time_;
			}
			else
			{
				d_max_time_for_quick_tests=1e99;
			};
		};

		if((max_time_with_computation_parameters_>0)&&(!(max_time_>0)))
		{
			d_max_time_with_computation_parameters=max_time_with_computation_parameters_;
		}
		else
		{
			d_max_time_with_computation_parameters=1e99;
		};

		calculate_RR_sum(
		d_RR1,
		alphabetSize_,
		d_RR1_sum,
		d_RR1_sum_elements);

		calculate_RR_sum(
		d_RR2,
		alphabetSize_,
		d_RR2_sum,
		d_RR2_sum_elements);


	}
	catch (...)
	{ 
		release_memory();
		throw;
	};
}

long int alp_data::random_long(
double value_,
long int dim_)
{
	if(value_<0||value_>1.0||dim_<=0)
	{
		throw error("Unexpected error\n",4);
	};

	if(dim_==1)
	{
		return 0;
	};

	long int tmp=(long int)floor(value_*(double)dim_);
	tmp=Tmin(tmp,dim_-1);
	return tmp;
}

void alp_data::release_memory()
{
	delete[]d_RR1;d_RR1=NULL;
	delete[]d_RR1_sum;d_RR1_sum=NULL;
	delete[]d_RR1_sum_elements;d_RR1_sum_elements=NULL;

	delete[]d_RR2;d_RR2=NULL;
	delete[]d_RR2_sum;d_RR2_sum=NULL;
	delete[]d_RR2_sum_elements;d_RR2_sum_elements=NULL;

	
	

	d_memory_size_in_MB-=2.0*(double)(2.0*sizeof(double)+sizeof(long int))*(double)d_number_of_AA/mb_bytes;



	if(d_smatr)
	{
		delete_memory_for_matrix(d_number_of_AA_smatr,d_smatr);
	};


	delete d_is;d_is=NULL;

	d_memory_size_in_MB-=sizeof(d_is)/mb_bytes;

	delete[]d_r_i_dot;d_r_i_dot=NULL;
	delete[]d_r_dot_j;d_r_dot_j=NULL;
	d_memory_size_in_MB-=(double)(sizeof(double)*d_number_of_AA*2.0)/mb_bytes;

	delete d_rand_all;d_rand_all=NULL;
	d_memory_size_in_MB-=sizeof(struct_for_randomization)/mb_bytes;

}

alp_data::~alp_data()//destructor
{
	release_memory();
}

void alp_data::check_out_file(
	string out_file_name_)
{
	ifstream f;
	char *str_ch=NULL;

	try
	{
		f.open(out_file_name_.data(),ios::in);
		if(!f)
		{
			return;
		};

		bool symmetric_case_flag;
		
		string str;
		getline(f,str);
		str_ch=new char[str.length()+1];
		alp_data::assert_mem(str_ch);

		long int k;
		for(k=0;k<(long int)str.length();k++)
		{
			str_ch[k]=str[k];
		};
		str_ch[str.length()]='\0';


		char str_for_test0[]="number of realizations with killing";
		char *test_flag0= strstr(str_ch,str_for_test0);

		if(!test_flag0)
		{
			throw error("The output file "+out_file_name_+" exists and does not have the correct format;\nplease delete the file and rerun the program\n",3);
		};

		char str_for_test[]="0.5*";

		char*test_flag= strstr(str_ch,str_for_test);
		if(test_flag)
		{
			symmetric_case_flag=true;
		}
		else
		{
			symmetric_case_flag=false;
		};


		

		if(symmetric_case_flag)
		{
			if(!d_smatr_symmetric_flag)
			{
				throw error("The output file "+out_file_name_+" exists and corresponds to symmetric case; \nthe current calculation uses non-symmetric parameters;\nplease define another output file name\n",3);
			};
		};

		if(!symmetric_case_flag)
		{
			if(d_smatr_symmetric_flag)
			{
				throw error("The output file "+out_file_name_+" exists and corresponds to non-symmetric case; \nthe current calculation uses symmetric parameters;\nplease define another output file name\n",3);
			};
		};

		f.close();
		delete[]str_ch;str_ch=NULL;
	}
	catch (...)
	{ 
		delete[]str_ch;str_ch=NULL;

		if(f.is_open())
		{
			f.close();
		};

		throw;
	};

}

double importance_sampling::lambda_equation(double x_,void* func_number_)
{
	data_for_lambda_equation *data=(data_for_lambda_equation*)func_number_;
	long int d_number_of_AA=data->d_number_of_AA;
	long int** d_smatr=data->d_smatr;
	double *d_RR1=data->d_RR1;
	double *d_RR2=data->d_RR2;

	double res=0;
	long int i,j;

	for(i=0;i<d_number_of_AA;i++)
	{
		for(j=0;j<d_number_of_AA;j++)
		{
			res+=d_RR1[i]*d_RR2[j]*exp(x_*d_smatr[i][j]);
		};
	};

	return res-1.0;
}

void alp_data::read_smatr(
string smatr_file_name_,
long int **&smatr_,
long int &number_of_AA_smatr_)
{
	ifstream f;

	try
	{

		long int i,j;
		f.open(smatr_file_name_.data(),ios::in);
		if(!f)
		{
			throw error("Error - file "+smatr_file_name_+" is not found\n",3);
		};

		f>>number_of_AA_smatr_;

		if(number_of_AA_smatr_<=0)
		{
			throw error("Error - number of letters in the scoring matrix file must be greater than 0\n",3);
		};

		get_memory_for_matrix(number_of_AA_smatr_,smatr_);


		for(i=0;i<number_of_AA_smatr_;i++)
		{
			for(j=0;j<number_of_AA_smatr_;j++)
			{
				f>>smatr_[i][j];
			};
		};

		f.close();

	}
	catch (...)
	{ 
		if(f.is_open())
		{
			f.close();
		};
		throw;
	};


}

void alp_data::read_RR(
string RR_file_name_,
double *&RR_,
double *&RR_sum_,
long int *&RR_sum_elements_,
long int &number_of_AA_RR_)
{

	read_RR(
	RR_file_name_,
	RR_,
	number_of_AA_RR_);

	calculate_RR_sum(
	RR_,
	number_of_AA_RR_,
	RR_sum_,
	RR_sum_elements_);
}

void alp_data::check_RR_sum(
double sum_tmp_,
long int number_of_AA_RR_,
string RR_file_name_)
{

	if(number_of_AA_RR_<=0)
	{
		throw error("Error - number of letters in the probabilities file must be greater than 0\n",3);
	};

	double diff_tmp=fabs(sum_tmp_-1.0);
	if(diff_tmp>0)
	{
		double lg_diff=-(log(diff_tmp)-log((double)number_of_AA_RR_))/log(10.0);
		double lg_eps=-log(DBL_EPSILON)/log(10.0)-1;
		if(lg_diff<lg_eps)
		{

			if(sum_tmp_<=0)
			{
				if(RR_file_name_!="")
				{
					throw error("Error: the sum of the probabilities from the file "+RR_file_name_+" is non-positive\n",3);
				}
				else
				{
					throw error("Error: the sum of the probabilities is non-positive\n",3);
				};

			};

			if(RR_file_name_!="")
			{
				static map<string, bool> flag_RR;

				if(!flag_RR[RR_file_name_])
				{
					cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
					cout<<"Warning: the sum of the probabilities from the file "<<RR_file_name_<<" is not equal to 1\n";
					cout<<"The probabilities will be normalized for the computation\n";
					cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n";

					flag_RR[RR_file_name_]=true;
				};
			}
			else
			{
				//no messages if called from the library functions
			};

		};

	};



}

void alp_data::calculate_RR_sum(
double *RR_,
long int number_of_AA_RR_,
double *&RR_sum_,
long int *&RR_sum_elements_)
{
	RR_sum_=NULL;
	RR_sum_elements_=NULL;

	try
	{

		long int i;
		if(number_of_AA_RR_<=0)
		{
			throw error("Error - number of letters in the probabilities file must be greater than 0\n",3);
		};
		
		RR_sum_=new double[number_of_AA_RR_];
		assert_mem(RR_sum_);

		RR_sum_elements_=new long int [number_of_AA_RR_];
		assert_mem(RR_sum_elements_);


		for(i=0;i<number_of_AA_RR_;i++)
		{
			if(RR_[i]<0)
			{
				throw error("Error - the frequencies must be non-negative\n",3);
			};

			if(i!=0)
			{
				RR_sum_[i]=RR_sum_[i-1]+RR_[i];
			}
			else
			{
				RR_sum_[i]=RR_[i];
			};
			RR_sum_elements_[i]=i;
		};

		double sum_tmp=RR_sum_[number_of_AA_RR_-1];

		check_RR_sum(
		sum_tmp,
		number_of_AA_RR_,
		"");

		if(sum_tmp>0)
		{
			long int i;
			for(i=0;i<number_of_AA_RR_;i++)
			{
				RR_[i]/=sum_tmp;
				RR_sum_[i]/=sum_tmp;
			};
		};

	}
	catch (...)
	{ 
		delete[]RR_sum_;RR_sum_=NULL;
		delete[]RR_sum_elements_;RR_sum_elements_=NULL;
		throw;
	};


}


void alp_data::read_RR(
string RR_file_name_,
double *&RR_,
long int &number_of_AA_RR_)
{

	ifstream f;
	RR_=NULL;

	try
	{

		long int i;
		f.open(RR_file_name_.data(),ios::in);
		if(!f)
		{
			throw error("Error - file "+RR_file_name_+" is not found\n",3);
		};

		f>>number_of_AA_RR_;

		if(number_of_AA_RR_<=0)
		{
			throw error("Error - number of letters in the probabilities file must be greater than 0\n",3);
		};
		
		RR_=new double[number_of_AA_RR_];
		assert_mem(RR_);


		double sum_tmp=0;
		for(i=0;i<number_of_AA_RR_;i++)
		{
			f>>RR_[i];

			if(RR_[i]<0)
			{
				throw error("Error - the frequencies defined in the file "+RR_file_name_+" must be non-negative\n",3);
			};

			sum_tmp+=RR_[i];

		};

		check_RR_sum(
		sum_tmp,
		number_of_AA_RR_,
		RR_file_name_);

		f.close();
	}

	catch (...)
	{ 
		delete[]RR_;RR_=NULL;
		if(f.is_open())
		{
			f.close();
		};
		throw;
	};

}



string alp_data::long_to_string(//convert interer ot string
long int number_)
{
	string res_="";
	string tmp_string;
	if(number_>0)
	{
		tmp_string="";
	}
	else
	{
		if(number_==0)
		{
			tmp_string="";
		}
		else
		{
			tmp_string="-";
		};
	};
	number_=abs(number_);

	for( ; ; )
	{
		long int reminder=number_%10;
		number_=(number_-reminder)/10;
		res_=digit_to_string(reminder)+res_;
		if (number_==0)
		{
			break;
		};
	};

	return tmp_string+res_;
}

char alp_data::digit_to_string(//convert interer ot string
long int digit_)
{
	switch(digit_)
	{
	case 0:return '0';
	case 1:return '1';
	case 2:return '2';
	case 3:return '3';
	case 4:return '4';
	case 5:return '5';
	case 6:return '6';
	case 7:return '7';
	case 8:return '8';
	case 9:return '9';
	default:return '?';
	};
}

importance_sampling::importance_sampling(
alp_data *alp_data_,
long int open_,

long int epen_,

double temperature_,
long int number_of_AA_,
long int **smatr_,
double *RR1_,
double *RR2_)
{
	d_elements=NULL;
	d_elements_values=NULL;

	d_exp_s=NULL;

	d_alp_data=alp_data_;
	if(!d_alp_data)
	{
		throw error("Unexpected error\n",4);
	};

	try
	{



		{

			//calculation of the importance sampling theta

			data_for_lambda_equation tmp_ptr;
			tmp_ptr.d_number_of_AA=number_of_AA_;
			tmp_ptr.d_RR1=RR1_;
			tmp_ptr.d_RR2=RR2_;
			tmp_ptr.d_smatr=smatr_;

			//calculate maximum of smatr_ elements
			long int smatr_max=smatr_[0][0];
			long int smatr_max_i=0;
			long int smatr_max_j=0;
			long int smatr_min=smatr_[0][0];

			long int smatr_pos_max=LONG_MIN;
			long int smatr_neg_min=LONG_MAX;

			double eps=0.00001;
			double threshold=DBL_MIN*10.0;

			double aver_score=0;
			long int i,j;
			for(i=0;i<number_of_AA_;i++)
			{
				for(j=0;j<number_of_AA_;j++)
				{
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					//if(RR1_[j]*RR2_[i]<=threshold)
					if(RR1_[i]*RR2_[j]<=threshold)
					{
						continue;
					};
										


					aver_score+=RR1_[i]*RR2_[j]*smatr_[i][j];

					if(smatr_max<smatr_[i][j])
					{
						smatr_max=smatr_[i][j];
						smatr_max_i=i;
						smatr_max_j=j;
					};
					smatr_min=alp_data::Tmin(smatr_min,smatr_[i][j]);
					

					if(smatr_[i][j]>0)
					{
						smatr_pos_max=alp_data::Tmax(smatr_pos_max,smatr_[i][j]);
					};

					if(smatr_[i][j]<0)
					{
						smatr_neg_min=alp_data::Tmin(smatr_neg_min,smatr_[i][j]);
					};

				};
			};

			if(aver_score>=-threshold)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};

			if(smatr_max<=0)
			{
				throw error("Error - at least one element of the scoring matrix must be positive\n",3);
			};

			

			double a=eps;

			while(importance_sampling::lambda_equation(a,(void*)(&tmp_ptr))>0)
			{
				a/=2.0;

				if(a<threshold*100.0)
				{
					throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
				};
			};

			if(a<threshold*100.0)
			{
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};

			eps=a/10.0;


			double tmp_pr=RR1_[smatr_max_i]*RR2_[smatr_max_j];
			double b=(log(1+10*eps)-log(tmp_pr))/(double)smatr_max;

			
			long int n_partition=2;
			vector<double> res_lambda;
			
			
			alp_reg::find_tetta_general(
			importance_sampling::lambda_equation,
			(void*)(&tmp_ptr),
			a,
			b,
			n_partition,
			eps,
			res_lambda);

			sort(res_lambda.begin(),res_lambda.end());

			if(res_lambda.size()==0)
			{
				//throw error("Error - the program is not able to find the ungapped lambda. The program does not work with the input scoring scheme. \n",3);
				throw error("Error - you have exceeded the calculation time or memory limit.\nThe error might indicate that the regime is linear or too close to linear to permit efficient computation.\nPossible solutions include changing the randomization seed, or increasing the allowed calculation time and the memory limit.\n",3);
			};

			d_lambda=res_lambda[res_lambda.size()-1];
			d_ungap_lambda=d_lambda;

			//cout<<"\nUngapped lambda is "<<d_ungap_lambda<<endl;

			d_lambda*=temperature_;
		};


		
		d_is_number_of_AA=number_of_AA_;

		d_elements=new q_elem[number_of_AA_*number_of_AA_];
		alp_data::assert_mem(d_elements);

		d_elements_values=new double[number_of_AA_*number_of_AA_];
		alp_data::assert_mem(d_elements_values);



		d_alp_data->get_memory_for_matrix(d_is_number_of_AA,d_exp_s);

		long int ind=0;
		double sum=0;
		long int a,b;
		for(a=0;a<number_of_AA_;a++)
		{
			for(b=0;b<number_of_AA_;b++)
			{
				d_exp_s[a][b]=exp(d_lambda*smatr_[a][b]);
				d_elements_values[ind]=RR1_[a]*RR2_[b]*d_exp_s[a][b];
				sum+=d_elements_values[ind];
				ind++;
			};
		};


		for(a=0;a<number_of_AA_;a++)
		{
			for(b=0;b<number_of_AA_;b++)
			{
				d_exp_s[a][b]/=sum;
			};

		};


		for(ind=0;ind<number_of_AA_*number_of_AA_;ind++)
		{
			d_elements_values[ind]/=sum;
		};

		
		for(ind=1;ind<number_of_AA_*number_of_AA_;ind++)
		{
			d_elements_values[ind]=d_elements_values[ind-1]+d_elements_values[ind];
		};

		
		ind=0;
		for(a=0;a<number_of_AA_;a++)
		{
			for(b=0;b<number_of_AA_;b++)
			{
				q_elem elem_tmp;

				elem_tmp.d_a=a;
				elem_tmp.d_b=b;

				d_elements[ind]=elem_tmp;
				d_elements_values[ind]=d_elements_values[ind];

				ind++;

			};
		};



		d_mu=exp(-fabs(d_lambda)*open_);
		d_nu=exp(-fabs(d_lambda)*epen_);

		double tmp=1+d_mu-d_nu;

		d_eta=(1-d_nu)*(1-d_nu)/(tmp*tmp);
		d_mu_SI=1-d_nu;
		d_mu_IS=d_mu*(1-d_nu)/(tmp*tmp);
		d_mu_DS=d_mu/tmp;
		d_mu_SD=(1-d_nu)*(1-d_nu)/tmp;
		d_mu_ID=d_mu*(1-d_nu)/tmp;


		d_for_D[0]=d_nu;				d_for_D_states[0]='D';
		d_for_D[1]=d_for_D[0]+d_mu_SD;	d_for_D_states[1]='S';
		d_for_D[2]=d_for_D[1]+d_mu_ID;	d_for_D_states[2]='I';

		d_for_I[0]=d_nu;				d_for_I_states[0]='I';
		d_for_I[1]=d_for_I[0]+d_mu_SI;	d_for_I_states[1]='S';

		d_for_S[0]=d_eta;				d_for_S_states[0]='S';
		d_for_S[1]=d_for_S[0]+d_mu_DS;	d_for_S_states[1]='D';
		d_for_S[2]=d_for_S[1]+d_mu_IS;	d_for_S_states[2]='I';

		d_alp_data->d_memory_size_in_MB+=sizeof(double)*number_of_AA_/mb_bytes;
		d_alp_data->d_memory_size_in_MB+=sizeof(q_elem)*number_of_AA_/mb_bytes;
	}
	catch (...)
	{
		this->~importance_sampling();
		throw;
	};

}

importance_sampling::~importance_sampling()
{
	delete []d_elements;d_elements=NULL;
	delete []d_elements_values;d_elements_values=NULL;

	if(d_alp_data)
	{
		d_alp_data->delete_memory_for_matrix(d_is_number_of_AA,d_exp_s);
		d_alp_data->d_memory_size_in_MB-=sizeof(double)*d_is_number_of_AA/mb_bytes;
		d_alp_data->d_memory_size_in_MB-=sizeof(q_elem)*d_is_number_of_AA/mb_bytes;
	};

}

double alp_data::error_of_the_sum(//v1_+v2_
double v1_error_,
double v2_error_)
{
	if(v1_error_>=1e100||v2_error_>=1e100)
	{
		return 1e100;
	};

	return sqrt(v1_error_*v1_error_+v2_error_*v2_error_);
}

double alp_data::error_of_the_product(//v1_*v2_
double v1_,
double v1_error_,
double v2_,
double v2_error_)
{
	if(v1_error_>=1e100||v2_error_>=1e100)
	{
		return 1e100;
	};

	double a1=(v1_+v1_error_)*(v2_+v2_error_);
	double a2=(v1_-v1_error_)*(v2_+v2_error_);
	double a3=(v1_+v1_error_)*(v2_-v2_error_);
	double a4=(v1_-v1_error_)*(v2_-v2_error_);

	double a=v1_*v2_;

	return Tmax(fabs(a1-a),fabs(a2-a),fabs(a3-a),fabs(a4-a));

}

double alp_data::error_of_the_sqrt(//sqrt(v1_)
double v1_,
double v1_error_)
{
	if(v1_error_>=1e100||v1_<0)
	{
		return 1e100;
	};

	double s=sqrt(v1_);
	double s1=sqrt(alp_data::Tmax(0.0,v1_-v1_error_));
	double s2=sqrt(alp_data::Tmax(0.0,v1_+v1_error_));

	return alp_data::Tmax(fabs(s-s1),fabs(s-s2));
}

double alp_data::error_of_the_ratio(//v1_/v2_
double v1_,
double v1_error_,
double v2_,
double v2_error_)
{
	if(v1_error_>=1e100||v2_error_>=1e100)
	{
		return 1e100;
	};

	if(v2_==0)
	{
		return 1e100;
	};

	if(v1_==0&&v1_error_==0)
	{
		return 0.0;
	};

	double a=v1_/v2_;

	if(((v2_+v2_error_)*v2_<=0))
	{
		double a3=(v1_+v1_error_)/(v2_-v2_error_);
		double a4=(v1_-v1_error_)/(v2_-v2_error_);
		return alp_data::Tmax(fabs(a-a3),fabs(a-a4));
	};

	if(((v2_-v2_error_)*v2_<=0))
	{
		double a1=(v1_+v1_error_)/(v2_+v2_error_);
		double a2=(v1_-v1_error_)/(v2_+v2_error_);
		return alp_data::Tmax(fabs(a-a1),fabs(a-a2));
	};


	double a1=(v1_+v1_error_)/(v2_+v2_error_);
	double a2=(v1_-v1_error_)/(v2_+v2_error_);
	double a3=(v1_+v1_error_)/(v2_-v2_error_);
	double a4=(v1_-v1_error_)/(v2_-v2_error_);

	return Tmax(fabs(a-a1),fabs(a-a2),fabs(a-a3),fabs(a-a4));
}

double alp_data::error_of_the_lg(//lg(v1_)
double v1_,
double v1_error_)
{
	if(v1_error_>=1e100||v1_<=0)
	{
		return 1e100;
	};

	return alp_data::Tmin(fabs(log(v1_)/log(10.0)),v1_error_/v1_/log(10.0));
}

bool alp_data::the_value_is_double(
string str_,
double &val_)
{
	if(str_=="")
	{
		return false;
	};

	bool res=false;
	errno=0;
	char *p;
	val_=strtod(str_.c_str(),&p);
	if(errno!=0)
	{
		res=false;
	}
	else
	{
		res=(*p==0);
	};
	return res;
}

bool alp_data::the_value_is_long(
string str_,
long int &val_)
{

	if(str_.length()==0)
	{
		return false;
	};
	if(!(str_[0]=='+'||str_[0]=='-'||isdigit(str_[0])))
	{
		return false;
	};

	long int start_digit=0;

	if(str_[0]=='+'||str_[0]=='-')
	{
		start_digit=1;
	};


	long int i;
	for(i=start_digit;i<(long int)str_.size();i++)
	{
		if(!isdigit(str_[i]))
		{
			return false;
		};
	};

	if(((long int)str_.size()-start_digit)<=0)
	{
		return false;
	};

	if(((long int)str_.size()-start_digit)>1)
	{
		/*
		if(str_[start_digit]=='0')
		{
			return false;
		};
		*/

		while(str_[start_digit]=='0')
		{
			string::iterator it=str_.begin()+start_digit;


			str_.erase(it);
			if((long int)str_.size()<=start_digit+1)
			{
				break;
			};
		};
	};

	if(((long int)str_.size()-start_digit>10)||((long int)str_.size()-start_digit)<=0)
	{
		return false;
	};


	if((long int)str_.size()-start_digit==10)
	{
		if(!(str_[start_digit]=='1'||str_[start_digit]=='2'))
		{
			return false;
		};

		if(str_[start_digit]=='2')
		{

			long int val2;
			string str2=str_.substr(start_digit+1,9);
			int flag=sscanf(str2.c_str(),"%ld",&val2);
			if(flag!=1)
			{
				return false;
			};

			bool positive=true;
			if(start_digit>0)
			{
				if(str_[0]=='-')
				{
					positive=false;
				};
			};

			if(positive)
			{
				if(val2>147483647)
				{
					return false;
				};
			}
			else
			{
				if(val2>147483648)
				{
					return false;
				};
			};

		};
	};

	int flag=sscanf(str_.c_str(),"%ld",&val_);
	if(flag!=1)
	{
		return false;
	};

	return true;
}

