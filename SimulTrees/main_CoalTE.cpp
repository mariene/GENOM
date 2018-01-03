/*
	file:     main_CoalTE.cpp
	function: Perform Simulation n sequence where the MRCA had xx TE
	author:   <amikezor>
	date:     june 10
*/

#include <unistd.h>

#include <iostream>
using namespace std;

#include "population.h"
#include "sample.h"
#include "random.h"
#include "coalescent_tree_diversity.h"
#include "distribution.h"
#include "readfile.h"

#define _VERSION_ "June 17, 2010"

char opt_verbose;


static void printusage( char *argv0 ){

		cerr << "Syntax is "<<argv0<<" [option] (-m #TE_mrca | -t time_TE_inv) dupl_rate #seq\n\n";

		cerr << " [ Choice between ]\n";
		cerr << "\t* -m #TE_mrca         : Number of transposons in the MRCA\n";
		cerr << "\t* -t time_1st_TE      : set the time TE invasion began\n";

		cerr << " [ Required ]\n";
		cerr << "\t* dupl_rate           : duplication rate (of transposons)\n";
		cerr << "\t* #seq                : Number of sequences you have sampled in the clade\n";

		cerr << "\n";
		
		cerr << " [ Optionnal ]\n";
		cerr << "\t -T #                 : set Theta, the scaled mutation rate (2pNmu, N # chr and mu mut rate /gen /locus). default is 10.\n";
		cerr << "\t -p                   : set ploidy. Default is 1 (haploid).\n";
		cerr << "\t -r #                 : set number of replicates. default is 10^4.\n";
		cerr << "\t -x #                 : set seed for random generator (default is time())\n";
}


static void printheader( population *my_population, int rep, int n, int n_TE_mrca, float t_first_TE, float duplik_rate ){


	extern Random R;

	
	
	cout << "/* Perform simulations to generate a TE tree within a genome tree - NO RECOMBINATION - ("<<_VERSION_<<")\n";
	
	cout << "/*    Population:\n";
	cout << "/*      * Theta        = "<< my_population->get_Theta() <<" (2.Ne.p.mu)\n";
	cout << "/*      * Ploidy       = "<< (int)my_population->get_ploidy() <<" \n";

	cout << "/*    Genomes sample:\n";
	cout << "/*      * #seq         = "<< n <<"\n";

	cout << "/*    TEs:\n";
	if(n_TE_mrca != -1)
	cout << "/*      * #in MRCA     = "<< n_TE_mrca <<"\n";
	if(t_first_TE != -1)
	cout << "/*      * time_1st_TE  = "<< t_first_TE <<"\n";
	cout << "/*      * duplik_rate  = "<< duplik_rate <<"\n";

	
	cout << "/*    Iterations:\n";
	cout << "/*      * rep     = "<<rep<<" iterations\n";
	
	cout << "/*    seed random generator:\n";
	cout << "/*      * seed    = "<<R.get_seed()<<"\n";

	cout << "\n";
}



int main(int argc, char *argv[]){
	
	/*
		Parse options
	*/
	
	int optc;
	extern char *optarg;
	extern int optind;

	char ploidy=1;               // 1 for haploid, 2 for diploid
	char opt_verbose=0;          // become verbose when set to 1, option '-v'



	float dupl_rate;
	int n;
	
	int n_te_mrca=-1;
	float t_first_TE=-1;


	/*
		Default values
	*/
	
	float Theta=10;
	int rep=10000;

	extern Random R;

	/*
		parse input
	*/
	while( (optc=getopt(argc, argv, "T:p:r:vx:t:m:")) != -1 ){
	
		switch(optc){
		
			case 'T':
				Theta = atof(optarg);
				break;

			case 'r':
				rep = atol(optarg);
				break;

			case 'p':
				ploidy = (char)atoi(optarg);
				break;

			case 'v':
				opt_verbose=1;
				break;

			case 'x':
				R.set_seed( atol(optarg) );
				break;

			case 'm':
				n_te_mrca   = atol(optarg);
				if(n_te_mrca < 0){
					cout << "#TE in MRCA has be be >=0\n";
					exit(1);
				}
				break;

			case 't':
				t_first_TE = (float)atof(optarg);
				if(t_first_TE < 0){
					cout << "time of the first TE has be be >=0\n";
					exit(1);
				}
				break;


		}	
	}
	

 	if(argc-optind != 2)
		printusage( argv[0] ),exit(1);		 		
	
	dupl_rate   = atof(argv[optind]);
	n     	    = atoi(argv[optind+1]);

	if( (t_first_TE == -1 && n_te_mrca == -1) || ( t_first_TE != -1 && n_te_mrca != -1 ) ){
		cout << "Please set either time of the first TE (-t) or number of TE in the genome mrca\n";
		exit(1);
	}



	/*
		Set Pop, built sp tree
	*/

	population my_population;                                    // with Ne=1, mu=1 and haploid;
	my_population.set_Theta( Theta );
	my_population.set_ploidy( ploidy );

	/*
		Init the sequences samples
	*/
	coalescent_tree_diversity *my_tree = new coalescent_tree_diversity( n, &my_population );

	printheader( &my_population, rep, n, n_te_mrca, t_first_TE,  dupl_rate );

	
	for(int r=0; r<rep; r++){
	
		//cerr << "rep: "<<r<<"\n";
		
		locus * mrca_locus;

		my_tree->standard_coalescent( );
		mrca_locus = my_tree->get_oldest_sample()->get_loci()[0];
		
		if( n_te_mrca>0 )
			my_tree->coal_TE_TopDown_mrca( n_te_mrca, dupl_rate);
		else
			my_tree->coal_TE_TopDown_1stTE( t_first_TE, dupl_rate);
		
		my_tree->printout_tree( );
		
	//	for ( int i=0; i < my_tree->get_n() ; i++){
	//		my_tree -> get_all_samples()[i]->print_all_TE_connections();
	//	}
		
		my_tree->printout_tree_TE();

		my_tree->mutate_whole_tree( );

		
		
		
		my_tree->clean_whole_tree( );

	}


	delete my_tree;   

	return 0;
}
