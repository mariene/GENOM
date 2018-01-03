/***
		file:     distribution.h
		function: header of a distribution.cpp, a dedicated class for distribution handling
		author:   <amikezor>
		date:     sep 04
***/

#ifndef _DISTRIBUTION_H_
#define _DISTRIBUTION_H_


#include <stdio.h>

class distribution {


	private:
		bool full;            // If true, all values are stored. If false only the subset needed for CI is displayed

		float CI;             // Confidence Interval one wants to estimate, it conditions the array_size. Two tailed.
		int rep;              // total number of element in the distribution

		double *array;        // Values that need to be in memory
		double *w;            // weight vector
		int array_size;       // How many values.
		
		int c;                // current number of element in the distribution
		
		double sum;           // sum of all elements. Needed for mean
		double sum_sq;        // summ of all elements squared. needed for variance.
		
		double sum_w ;        // The total sum of the weight

		
	public:
		distribution(int rep, float initCI=0.95, bool full=false, bool weight=false);    // constructor
		~distribution(void);                                                             // destructor
		
		
		double get_mean( void ){return sum/(double)sum_w;};
		double get_variance( void ){ return sum_sq/(double)sum_w - get_mean()*get_mean(); };
		
		int get_array_size(void){return array_size;};
		
		void add_value( double x, double weight=-1 );
		
		double get_low_boundary( void );
		double get_high_boundary( void );
		
		void print_distribution( void );

		void clean_distribution( void );
		

};



#endif
