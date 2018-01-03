/***
		file:     random.h
		function: header for the random functions
		author:   <amikezor>
		date:     sep 04
***/

#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <iostream>
using namespace std;


/***
	The Random class
	instantiate an object and you seed, otherwise you cannot.
***/
class Random{
	
	

	private:
	
		long seed;
	
		static long idum_ran1;                                    // the seed itself
		void seed_ran1( void ){seed_ran1( 0L );};                  // seed the function, only once
		void seed_ran1( long );                                   // seed the function, only once
		double ran1( void );                                      // ran1() returns a ]0,1]

	public:
		Random( void ){this->seed_ran1();};                       // constructor
		
		long get_seed(void ){return seed;};
		void set_seed( long x ){seed_ran1(x);};
		
		double uniform_dev( void );                               // run a uniform [0,1[ number
		double exponential_dev( double mean );                    // run an expornential with the mean m
		double cauchy_dev( double rate );                         // run an Cauchy 1/(1+rx), with a rate r
		double birthdeath_dev( double birth, double death );      // run a brith death r / (b exp^(rx) - d), with a rate r=b-d
		int poisson_dev( double mean );
		float normal_CR_dev(  );
		float normal_dev( double mean, double stdev );
		
		/*
			Pick discrete numbers
		*/
		int uniform_int(int min, int max);                        // uniform int in the range [min, max]
		void uniform_2DiffInt(int min, int max,
		                  int &int1, int &int2);                  // set int1 and int 2 as 2 different uniform int

		int *uniform_kDiffInt( int k, int min, int max);

		/*
			Combinatorial functions
		*/
		long n_choose_k( int n, int k);		
};

#endif
