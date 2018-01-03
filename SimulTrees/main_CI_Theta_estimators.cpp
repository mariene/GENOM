/*
	file:     main_Tests_estimators.cpp
	function: compute Confidence Interval for Theta estimators (for a given theta)
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

#define _VERSION_ "v1.0 (April 2010)"

#define ROUNDUP( A ) ( ((A)-(int)(A))>=0.5 )?( 1+(int)(A) ):( (int)(A) )

char opt_verbose;


static void printusage( void ){
		cerr<<"Syntax is 'Tests_estimators [-r #] [-x #] [-C #] n theta'\n";
		cerr<<"[options]\n";
		cerr<<"\t-C  : Confidence Interval, by default 0.95\n";
		cerr<<"\t-r  : number of replicates, default is 10^4\n";
		cerr<<"\t-x  : seed for random numbers, by default time() is used\n";
		cerr<<"[options]\n";
}





int main(int argc, char *argv[]){
	

	int n;
	double theta;
	
	double CI=0.95;
	long rep=10000;

	char ploidy=1;

	int optc;
	extern char *optarg;
	extern int optind;

	extern Random R;
	
	/*
		parse input
	*/
	while( (optc=getopt(argc, argv, "r:x:C:p:")) != -1 ){
	
		switch(optc){
		
			case 'r':
				rep = atol(optarg);
				break;

			case 'x':
				R.set_seed( atol(optarg) );
				break;

			case 'C':
				CI = atof(optarg);
				break;

		}	
	}

 	if(argc-optind != 2)
		printusage( ),exit(1);	


	n = atoi(argv[optind]);
	theta = atof(argv[optind+1]);

	/*
		Set the population parameters
	*/
	
	population my_population;                                    // with Ne=1, mu=1 and haploid;
	my_population.set_Theta( theta );
	my_population.set_ploidy( ploidy );
	 	
	
	/*
		Init the samples
	*/
	coalescent_tree_diversity *my_tree = new coalescent_tree_diversity( n, &my_population );
	


	/*
		Print out header
	*/

	cout << "/* CI of theta estimators "<< _VERSION_ <<"\n";

	cout << "/*    Population:\n";
	cout << "/*      * theta = "<< theta <<"\n";
	cout << "/*    Sample:\n";
	cout << "/*      * n     = "<< n <<" loci\n";
	cout << "/*    Distributions\n";
	cout << "/*      * CI    = "<< CI <<"\n";
	cout << "/*      * rep    = "<< rep <<"\n";
	cout << "\n";

	distribution *distrib_theta_pi;
	distribution *distrib_theta_S;


	distrib_theta_pi = new distribution(rep, CI );
	if( ! distrib_theta_pi  )fprintf(stderr, "main: cannot allocate distrib_pi, bye"), exit(3);

	distrib_theta_S = new distribution(rep, CI );
	if( ! distrib_theta_S  )fprintf(stderr, "main: cannot allocate distrib_pi, bye"), exit(3);


	for(int r=0; r<rep; r++){

		my_tree->standard_coalescent( );
		my_tree->mutate_whole_tree( );
		my_tree->compute_Xi_array();
		
		
	//	printf("an %f\n", my_tree->get_an() );
		
		distrib_theta_pi->add_value( my_tree->get_K() );
		distrib_theta_S->add_value( ((my_tree->get_S()+0.0)/my_tree->get_an()) );
		
		my_tree->clean_whole_tree( );
	}
	

	printf("Theta_pi\n");
	printf("E[pi]= %f ; Var[pi]= %f ; CI= [ %f , %f ]\n",
	       distrib_theta_pi->get_mean() , distrib_theta_pi->get_variance(),
	       distrib_theta_pi->get_low_boundary(), distrib_theta_pi->get_high_boundary() );


	printf("Theta_S\n");
	printf("E[S]= %f ; Var[S]= %f ; CI= [ %f , %f ]\n",
	       distrib_theta_S->get_mean() , distrib_theta_S->get_variance(),
	       distrib_theta_S->get_low_boundary(), distrib_theta_S->get_high_boundary() );


	delete distrib_theta_pi;
	delete distrib_theta_S;
	delete my_tree;                                       // this would free also ancestral sample (in the destructor)

	return 0;
}
