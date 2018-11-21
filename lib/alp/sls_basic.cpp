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

File name: sls_basic.cpp

Author: Sergey Sheetlin, Martin Frith

Contents: Some basic functions and types

******************************************************************************/

// 2016: this voodoo is needed to compile on Cygwin, with g++ options
// such as -std=c++11 or -std=c++03, else it complains about gettimeofday
#define _DEFAULT_SOURCE 1

#include "sls_basic.hpp"
#include <cstdlib>  // std::abs
#include <ctime>

using namespace Sls;

double sls_basic::round(//returns nearest integer to x_
const double &x_)
{
	double x_floor=floor(x_);
	double x_ceil=ceil(x_);
	if(fabs(x_-x_floor)<0.5)
	{
		return x_floor;
	};
	return x_ceil;
}

void sls_basic::get_current_time(
double &seconds_)
{
#ifndef _MSC_VER //UNIX program
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv, &tz);
	seconds_=(double)(tv.tv_sec)+(double)(tv.tv_usec)/1000000.0;

#else

	struct _timeb timebuffer;
	_ftime( &timebuffer );
	seconds_=timebuffer.time+(double)(timebuffer.millitm)/1000.0;

#endif
}

long int sls_basic::random_seed_from_time()
{
	long int random_factor=(long int)std::time(NULL);
#ifndef _MSC_VER //UNIX program
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv, &tz);
	random_factor+=tv.tv_usec*10000000;
#else
	struct _timeb timebuffer;
	char *timeline;
	_ftime( &timebuffer );
	timeline = ctime( & ( timebuffer.time ) );
	random_factor+=timebuffer.millitm*10000000;
#endif
	return std::abs(random_factor);
}

double sls_basic::one_minus_exp_function(
double y_)
{
	if(fabs(y_)>1e-3)
	{
		return 1.0-exp(y_);
	}
	else
	{
		return -(y_*(120+y_*(60+y_*(20+y_*(5.0+y_))))/120.0);
	};
}

double sls_basic::normal_probability(
double x_,
double eps_)
{

	if(x_==0)
	{
		return 0.5;
	};


	eps_=Tmin(1.0,eps_);

	double x_max=10*eps_+sqrt(Tmax(0.0,-2*log(eps_)));


	if(x_>=x_max)
	{
		double x=x_/sqrt(2.0);
		return 1-0.5*exp(-x*x)/(x*sqrt(pi))*(1-1.0/(2*x*2*x));
	};

	if(x_<=-x_max)
	{
		double x=x_/sqrt(2.0);
		return 0.5*exp(-x*x)/(-x*sqrt(pi))*(1-1.0/(2*x*2*x));
	};


	double const_val=1/sqrt(2.0*pi);

	


	long int N=(long int)round(fabs(x_)/eps_)+1;
	double h=x_/(double)N;



	double res=0;
	long int i;
	for(i=0;i<=N;i++)
	{
		double y=h*i;
		double tmp=exp(-0.5*y*y);
		if(i==0||i==N)
		{
			res+=0.5*tmp;
		}
		else
		{
			res+=tmp;
		};
	};

	res*=h;

	return 0.5+const_val*(res);
}

double sls_basic::normal_probability(
double a_,
double b_,
double h_,
long int N_,
double *p_,
double x_,
double eps_)
{
	if(x_<a_||x_>b_)
	{
		return normal_probability(x_,eps_);
	};

	long int x_n=(long int)floor((x_-a_)/h_);
	x_n=Tmin(N_-1,x_n);
	return p_[x_n]+(p_[x_n+1]-p_[x_n])*(x_-(h_*x_n+a_))/h_;
}

double sls_basic::ln_one_minus_val(
double val_)
{
	if(val_>1e-8)
	{
		return log(1-val_);
	};

	return -val_-val_*val_/2.0-val_*val_*val_/3.0;
}

