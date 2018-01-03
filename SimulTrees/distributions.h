/***
		file:     distribution.h
		function: header of a distribution.cpp, a dedicated class for distribution handling
		author:   <amikezor>
		date:     sep 04
***/

#ifndef _DISTRIBUTIONS_H_
#define _DISTRIBUTIONS_H_


#include <stdio.h>

class distributions {


	private:
		bool full;             // If true, all values are stored. If false only the subset needed for CI is displayed

		float CI;              // Confidence Interval one wants to estimate, it conditions the array_size. Two tailed.
		int   rep;             // total number of element in the distribution

		double **array;       // Values that need to be in memory
		int array_values;      // How many many different variables
		int array_size;        // size of the arrays.

		double *w;           // weight vector
		double sum_w ;        // The total sum of the weight
		
		int c;                 // current number of element in the distribution
		
		double *sum_X;         // sum of all elements of each col. Needed for mean
		double *sum_X2;        // summ of all elements squared. needed for variance.
		double *sum_XY;        // summ of all elements squared. needed for covariance.
		
		double **Covariance;    // The Covariance matrix
		
	public:
		distributions(int rep, float initCI=0.95, bool full=false, bool weight=false, int vals=2);    // constructor
		~distributions(void);                                                                         // destructor
		
		
		double get_value( int i, int j){ return this->array[i][j]; };
		double get_w( int i ){ return w[i]; };
		bool is_w( void ){ if(this->w==NULL){return false;}else{return true;} };
		
		void compute_Covariances( void );
		
		double get_mean( int i ){return sum_X[i]/(double)sum_w;};
		double get_variance( int i );
		double get_covariance( int i, int j );
		double **get_Covariances( void ){ return Covariance; };
		
		int get_array_size(void){return array_size;};
		int get_array_values(void){return array_values;};
		
		void add_values( double *X, double weight=-1 );
		
		double get_low_boundary( int i );
		double get_high_boundary( int i );
		
		void print_distribution( int i );

		void clean_distribution( void );
		

};



#endif
