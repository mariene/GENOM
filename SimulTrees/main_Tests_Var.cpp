/*
	file:     main_Tests_CI.S.n.cpp
	function: compute Confidence Interval for several test. Fixed S is used.
	author:   <amikezor>
	date:     May 06
	date:     Nov 06 - try without generating sequences (faster ?)
*/


#include <iostream>
using namespace std;

#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "population.h"
#include "sample.h"
#include "random.h"
#include "distribution.h"
#include "coalescent_tree_diversity.h"
#include "readfile.h"

#define _VERSION_ "v1.0 (May 06)"

#define ROUNDUP( A ) ( ((A)-(int)(A))>=0.5 )?( 1+(int)(A) ):( (int)(A) )


char opt_verbose;


static void printusage( void ){
		cerr<<"Syntax is Tests_Var w1_file w2_file [theta]\n\n";
		cerr<<"\t * w[12]_file     : file for the w1/w2 vector file\n";
		cerr<<"\t * [theta]        : set the population muta param. If not report alpha and beta\n";
}





int main(int argc, char *argv[]){
	

	readfile inputfiles;
	double theta=-1;
	double alpha=1, beta=-1;

	double *w_null;

	if(argc < 3)
		printusage( ),exit(1);	

	 	
	inputfiles.readfile_w( argv[1], 1);
	inputfiles.readfile_w( argv[2], 2);
	

	
	
	if(argc == 4)
		theta = atof( argv[3] );


	if( inputfiles.get_size_w1()  != inputfiles.get_size_w2() ){
		cerr << "w1 and w2 vectors should have the same size, bye\n";
		exit(1);
	}

	w_null = (double *)calloc( (size_t)  inputfiles.get_size_w1(), sizeof( double  ) );
	if(! w_null )fprintf(stderr, "main: cannot allocate w_null, bye\n"), exit(3);


	/*
		Set the population parameters
	*/
	
	population my_population;                                    // with Ne=1, mu=1 and haploid;
	 	
	
	/*
		Init the samples
	*/
	coalescent_tree_diversity *my_tree = new coalescent_tree_diversity( inputfiles.get_size_w2()+1, &my_population );
	


	/*
		Print out header
	*/

	cout << "/* Variance of the test "<< _VERSION_ <<"\n";
	if(theta != -1 ){
		cout << "/*    Population:\n";
		cout << "/*      * theta = "<< theta <<"\n";
	}
	cout << "/*    Sample:\n";
	cout << "/*      * n     = "<< inputfiles.get_size_w2()+1 <<" loci\n";
	cout << "/*    Vectors:\n";
	cout << "/*      * w1    = "<< argv[1]<<"\n";
	cout << "/*      * w2    = "<< argv[2]<<"\n";
	cout << "\n";


	/*
		Print out results
	*/
	if(theta != -1){
		cout <<"Variance(T): "<<my_tree->compute_VarT( inputfiles.get_w1(), inputfiles.get_w2(), theta)<<"\n";
	}
	else
	{
	
		my_tree->compute_VarTAlphaBeta( inputfiles.get_w1(), inputfiles.get_w2(), &alpha, &beta );
		printf("alpha: %.10f beta: %.10f\n", alpha, beta);
	}
	if(theta != -1){
	
		double c,v1,v2;
	
		c=my_tree->compute_CoVarw1w2( inputfiles.get_w1(), inputfiles.get_w2(), theta);
		v1=my_tree->compute_VarT( inputfiles.get_w1(), w_null, theta);
		v2=my_tree->compute_VarT( inputfiles.get_w2(), w_null, theta);
	
	
		cout <<"CovVariance(w1,w2): "<<c<<"\n";
		cout <<"Var(w1): "<<v1<<"\n";
		cout <<"Var(w2): "<<v2<<"\n";
		cout <<"r^2(w1,w2): "<<c*c/(v1*v2)<<"\n";
		cout <<"r(w1,w2): "<<c/sqrt(v1*v2)<<"\n";
		
	}
	
	
	
	free(w_null);
	
	delete my_tree;                                       // this would free also ancestral sample (in the destructor)

	return 0;
}
