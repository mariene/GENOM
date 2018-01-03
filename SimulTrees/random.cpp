/***
		file:     random.cpp
		function: all random functions
		author:   NumRec / amikezor
		date:     feb 04
		modif:    mar 09, add cauchy_dev
		          sept 09, add normal
***/



#include <iostream>
using namespace std;


#include "random.h"
#include <time.h>
#include <math.h>

#ifdef MACOSX
#include <float.h>
#elif defined(SGI)
#include <limits.h>
#elif defined(LINUX)
#include <values.h>
#endif

Random R;                          // a random object to be exported anywhere
long Random::idum_ran1=0;          // an internal random variable

char debug;

/**

	Private functions
	
**/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)

#define EPS FLT_MIN 
#define RNMX (1.0-EPS)              // RNMX should be the largest floating value less than 1 (numrec p 280)

/*
	Return a number between ]0.00,1.00] (( now certain because of a bug Apr 2007 ))
*/
double Random::ran1()
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (this->idum_ran1 <= 0 || !iy) {
		if (-this->idum_ran1 < 1) this->idum_ran1=1;
		else this->idum_ran1 = -this->idum_ran1;
		for (j=NTAB+7;j>=0;j--) {
			k=this->idum_ran1/IQ;
			this->idum_ran1=IA*(this->idum_ran1-k*IQ)-IR*k;
			if (this->idum_ran1 < 0) this->idum_ran1 += IM;
			if (j < NTAB) iv[j] = this->idum_ran1;
		}
		iy=iv[0];
	}
	k=this->idum_ran1/IQ;
	this->idum_ran1=IA*(this->idum_ran1-k*IQ)-IR*k;
	if (this->idum_ran1 < 0) this->idum_ran1 += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = this->idum_ran1;
	if ((temp=AM*iy) > RNMX) return (double)RNMX;
	else return (double)temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/**
	Seed the function ran1() with a negative long
	Should be used only ONCE before using ran1
*/
void Random::seed_ran1( long x )     
{

	if(x != 0)
		seed = x;
	else
		seed = (long)time(NULL);

	this->idum_ran1 = -seed; 
	this->ran1();                            /* seed the function */

}


/**

	Public functions
	
**/

// uniform_dev return a number [0,1[
double Random::uniform_dev( void ){

	double dum= ran1();
	
	if(dum == 0.0)cerr << "ran1() can also return 0\n";

	if(dum == 1.0)dum=0;
	
	return dum;
	
}

// exponen_dev is [0,inf[
// P( X>= x ) = exp(-x/mean)
//
double Random::exponential_dev( double mean )
{ 
	double U;
 
	U = this->uniform_dev();    // uniform_dev return a number [0,1[

 	U=(U==0.0)?1:U;
	                              // this is based on the inverse transformation 
	return -(mean*log(U));      // if X=F^-1(U), then X is F distributed see S.M. Ross (p.590)	
}	                              // numrec has a different approach interesting though, see Numrec(p.287)  
 
// cauchy_dev is [0,inf[
// P( X>= x ) = 1 / (1 + rate* x)
//
double Random::cauchy_dev( double rate )
{ 
	double U;
 
	U = this->uniform_dev();      // uniform_dev return a number [0,1[

	                                // this is based on the inverse transformation 
	return U/((1-U)*rate);      // if X=F^-1(U), then X is F distributed see S.M. Ross (p.590)	
}	                                // here we have x = U / [ (1-U)*rate ]
 

// birthdeath_dev is [0,inf[
// P( X>= x ) = r / (b exp^(rx) - d)
// with r = b-d
//
double Random::birthdeath_dev( double b, double  d )
{ 
	double U;
 
	U = this->uniform_dev();                              // uniform_dev return a number [0,1[
	                                
	return ( log( b - d*U ) - log( b*(1-U) ) ) / (b-d);   // this is based on the inverse transformation see S.M. Ross (p.590)	
}	                               
 

// Poisson_dev is [0,inf[
int Random::poisson_dev( double mean ){ 

	int em;
	double g,t; 
	 
	g=exp(-mean);
	
	if( g == 0.0 ){
		//cerr << "Watch out !! Random::poisson_dev(): cannot draw, so return mean\n";
		em = (int)round(  normal_dev( mean, sqrt(mean) ) );
		if(em<0){
			cerr << "Watch out !! Random::poisson_dev(): cannot draw -> use Normal that returned a neg. Set it to 0\n";
			em=0;
		}
		return em;
	}
	
	
	em = -1; 
	t=1.0; 
	do {
		em ++;
		t *= uniform_dev();   // cumulative draw until something happen (the unidev proba gets above the proba of no event)
		
	} while (t >= g); 


	return ( em );
} 



//inspired from numrec
//Normal dev is mean 0 and var 1
float Random::normal_CR_dev( )
{
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;


	if  (iset == 0) {
		do {
			v1=2.0*uniform_dev()-1.0;
			v2=2.0*uniform_dev()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

float Random::normal_dev( double mean, double stdev ){
	return normal_CR_dev( )*stdev + mean;
}


/*
	Pcik discrete numbers
*/

// uniform_int is [min , max]
int Random::uniform_int(int min, int max) 
{ 
	double dum;
 	int diff=max-min+1;
	
	dum = this->uniform_dev();
	 
	return (min + (int)floor(dum*(double)diff)); 
} 

void Random::uniform_2DiffInt(int min, int max, int &int1, int &int2){

	int1 = this->uniform_int( min,  max);
	do{
		int2 = this->uniform_int( min,  max);
	}while(int1 == int2);

}


int *Random::uniform_kDiffInt(int k, int min, int max ){

	int * from = new int[max-min+1];
	int * picked = new int[k];
	int r, npicked=0;;

	for(int i=0; i<max-min; i++)
		from[i] = min+i;

	while( npicked < k ){
		
		r = this->uniform_int( 0,  max-min-npicked);
		picked[npicked] = from [r];
		
		for(int i=r; i<max-min-npicked-1;i++)
			from[i] = from [i+1];
			
		from[max-min-npicked-1] = -1000000;
		
		npicked ++;
	}
	
	delete[] from;

	return picked;
}




/*
	Combinatorial
*/
long Random::n_choose_k( int n, int k){

	double res = 1;
	int min=-1;

	if( k < n-k )
		min = k;
	else
		min = n-k;
	
	for(int i=0; i< min; i++)
		res *= (n-i+0.0)/(i+1.0);

	return (long)res;
}
