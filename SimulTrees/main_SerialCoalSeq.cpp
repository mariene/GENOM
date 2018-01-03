/*
	file:     main_GenerCoalSeq.cpp
	function: Perform Simulation for s species with sequences within. Define by an entry a vector S.
	author:   <amikezor>
	date:     march 12
	modif:    
*/


/*
	C library
*/
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

/*
	C++ library
*/
#include <iostream>
#include <string>
using namespace std;

#include "population.h"
#include "sample.h"
#include "random.h"
#include "coalescent_tree_diversity.h"
#include "distribution.h"
#include "readfile.h"


char opt_verbose;


/*
	Print usage and help
*/
static void printusage( char *argv0 ){

	cerr << "Syntax is:\n";
	cerr <<  "\t"<< argv0<<" [ options ] n1 T12 n2 \n";

	cerr << "\t* arg n1 and n2 should be either a single value (1 supopulation) or\n"
	        "\t                 \"#:#:#..\" sampling_scheme ('2:10:4:1' translates to 4 subpopulations with respectively 2, 10, 4 and 1 sequences)\n";
	cerr << "\t* #:#:#...            : sampling_scheme ('2:10:4:1' translates to 4 subpopulations with respectively 2, 10, 4 and 1 sequences)\n";

	cerr << "\n";
	
	cerr << " [ options ]\n";
	
	cerr << "\t -N #:#:#             : set population size ('1:10:8' means that population are of size N, 10N, 8N). Default all population are 1\n";
	cerr << "\t                        !! it has to coincide with the number of species sampled. !!\n";
	cerr << "\t -M #                 : set the 'M'igration rate\n";
	cerr << "\t -T #                 : set Theta, the scaled mutation rate (2pNmu, N # chr and mu mut rate /gen /locus). Default is 10.\n";
	cerr << "\t -p #                 : set ploidy (default is 1)\n";
	cerr << "\t -x #                 : set seed for random generator (default is time())\n";
	cerr << "\t -q                   : 'q'uiet mode = only output sequences\n";
}


/*
	Print header with all selected parameters
*/
static void printheader( population *my_population, int rep, int n, int ns, int *n1 , int *n2, float T12, char *population_sizes ){


	extern Random R;

	cout << "/* Generate sequences with two serial samples\n";

	cout << "/*    Populations:\n";
	cout << "/*      * Theta    = "<< my_population->get_Theta() <<" (2.Ne.p.mu)\n";
	cout << "/*      * Ploidy   = "<< (int)my_population->get_ploidy() <<" \n";
	cout << "/*      * #subpop  = "<< ns <<"\n";
	
	if(population_sizes != NULL)
	cout << "/*      * Pop_size = "<< population_sizes <<"\n";

	cout << "/*    Samples:\n";
	
	cout << "/*      * sample1 = "<< n1[0];
	for(int i=1;i<ns;i++)
		cout << ":" << n1[i];
	cout << "\n";

	cout << "/*      * T12 = "<< T12 << "\n";
	
	cout << "/*      * sample2 = "<< n2[0];
	for(int i=1;i<ns;i++)
		cout << ":" << n2[i];
	cout << "\n";
	
	cout << "/*    Iterations:\n";
	cout << "/*      * rep      = "<<rep<<" iterations\n";
	
	cout << "/*    seed random generator:\n";
	cout << "/*      * seed     = "<<R.get_seed()<<"\n";

	cout << "\n";
}



/*
	sp_arg should be numbers separated by ':'
	return an array with those numbers
	and zeros up to 2*ns-1 --for ancestral species--
*/
int * parse_sparg(  char *sp_arg , int *n_sp ){

	
	char *p;
	int *ni;
	int i=0;
	
	*n_sp=1;
	p = sp_arg;
	
	while( *(p++) )
		if( *p == ':' )
			(*n_sp)++;

	ni = (int *)calloc( (size_t)(2*(*n_sp)-1), (size_t)sizeof(int) );
	if( !ni )fprintf(stderr, "parse_sparg: cannot allocate ni, bye\n"),exit(1);

	for (p = strtok(sp_arg, ":"); p; p= strtok(NULL, ":"))
		ni[i++] = atoi(p);
	
	p=sp_arg;
	for(i=0;i<(*n_sp)-1;i++){
		while(*(++p));
		*p=':';
	}
	
	
	return ni;
}

/*
	same as previous but for float *
*/
double * parse_sparg_f(  char *sp_arg , int *n_sp ){

	
	char *p;
	double *f;
	int i=0;
	
	*n_sp=1;
	p = sp_arg;
	
	while( *(p++) )
		if( *p == ':' )
			(*n_sp)++;

	f = (double *)calloc( (size_t)(*n_sp), (size_t)sizeof(double) );
	if( !f )fprintf(stderr, "parse_sparg: cannot allocate f, bye\n"),exit(1);

	for (p = strtok(sp_arg, ":"); p; p= strtok(NULL, ":")){
		f[i++] = atof(p);
	}
	
	p=sp_arg;
	for(i=0;i<(*n_sp)-1;i++){
		while(*(++p));
		*p=':';
	}
	
	
	return f;
}



/*
	Main funtion
*/
int main(int argc, char *argv[]){
	
	int tmp, i;
	

	char **sequences=NULL;

	/*
		Population(s) parameters
	*/

	double migration_rate=-1;      // migration rate 2Nm

	char ploidy=1;                 // 1 for haploid, 2 for diploid
	
	char *population_sizes=NULL;   // if set to NULL all population have the same size, otherwise, they can be of different size (e.g. "1:0.1:0.01").
	double *input_fN=NULL;         // population_sizes string turn into an array of doubles.
	double sum_fN=0;               // sum of all pop sizes
	
	double Theta=10;              // Population Mutation rate 2Nmu -> mutation rate on a branch is Theta/2



	/*
		Sample values
	*/
	
	int n=0;                        // total sample size (sum of all ni)

	char *sampling_scheme1=NULL;    // if set to NULL, species sampling is randomnly done, if not respect the sampling scheme
	char *sampling_scheme2=NULL;    // if set to NULL, species sampling is randomnly done, if not respect the sampling scheme
	
	int n1;                      // total number of loci in today's bpopulations
	int n2;                      // total number of loci in the population that must be added at T12
	int *sample1, *sample2;      // number of loci in today's subpopulations. Sum is n1. number of loci that were sampled previously at time T1
	int n_sp=0;                   // number of species

	double T12=0;

	/*
		runs and their results
	*/
	int rep=1;                // # of replicates


	/*
		Other values
	*/
	extern Random R;              // the random generator object
	extern char opt_verbose;             // become verbose when set to 1, option '-v'
	opt_verbose=2;

	/*
		Parse options
	*/

	int optc;
	extern char *optarg;
	extern int optind;

	while( (optc=getopt(argc, argv, "T:p:x:qr:DN:M:l:")) != -1 ){
	
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

			case 'q':
				opt_verbose=0;
				break;

			case 'x':
				R.set_seed( atol(optarg) );
				break;

			case 'N':
				population_sizes = optarg;
				break;

			case 'M':
				migration_rate = atof(optarg);
				break;

		}
	}
	
	

 	if(argc-optind != 3  )
		printusage( argv[0] ),exit(1);	
	 		
	/*
		Raw arguments
	*/
	sampling_scheme1 = argv[optind];
	T12 = atof( argv[optind+1 ]);
	sampling_scheme2 = argv[optind+2];
		
	
	/*
		Parse the sampling scheme to extract samples
	*/
	sample1 = parse_sparg( sampling_scheme1 , &n_sp  );
	sample2 = parse_sparg( sampling_scheme2 , &tmp );
	
	if( n_sp != tmp ){
		cerr << "both sampling scheme should have the same number of species\n";
		exit(5);
	}
	for( i=0, n1=0;  i<n_sp;  n1+=sample1[i++]);
	for( i=0, n2=0;  i<n_sp;  n2+=sample2[i++]);

	
	/*
		Built arrays for functions of coalescence
	*/
	
	int **all_samples_struct = new int *[2];
	all_samples_struct[0] = sample1;
	all_samples_struct[1] = sample2;

	double *serial_times = new double[2];
	serial_times[0]= 0.0;
	serial_times[1]= T12;

	int *all_samples = new int[2];
	all_samples[0] = n1;
	all_samples[1] = n2;

	n=n1+n2;

	
	/*
		If population sizes are not all equal, fill the array
	*/	
	if(population_sizes){
		
		input_fN = parse_sparg_f( population_sizes , &tmp );
		if(tmp != n_sp){
			cerr << "the species size are not compatible with the number of species, bye\n";
			exit(1);
		}
		
		for(int i=0;i<tmp;i++)
			sum_fN+=input_fN[i];
	}else{
	
		
	
	}

	


	/*
		Set Pop, built sp tree
	*/

	population my_population;                                    // with Ne=1, mu=1 and haploid;
	my_population.set_Theta( Theta );
	my_population.set_ploidy( ploidy );
	my_population.set_ploidy( ploidy );

	if(opt_verbose )
		printheader( &my_population, rep, n, n_sp, sample1 ,sample2, T12, population_sizes );



	/*
		Init the sequences samples
	*/
	coalescent_tree_diversity *my_tree = new coalescent_tree_diversity( n1+n2, &my_population );

	if(n_sp >1)
		my_population.set_structure( input_fN, n_sp, migration_rate );

	/*
		Run the simulations
	*/
	
	for(int r=0; r<rep; r++){
	
		cerr << "rep: "<<r<<"\n";
		

		/*
			Generate species tree
		*/
	
		
		if(migration_rate>0 && n_sp>1 && opt_verbose)
			cout << "/*\n\tMigration Event(s)\n*/\n";

		/*
			Sequence tree INSIDE the species tree
		*/
		if(n_sp>1)
			my_tree->structured_coalescent( 2, all_samples_struct, serial_times  );
		else
			my_tree->standard_coalescent( 2, all_samples, serial_times );



		/*
			Add mutations 
		*/
		my_tree->mutate_whole_tree_seq(  );


		if(opt_verbose)
			my_tree->printout_tree();
		

		sequences = my_tree->extract_sequences(  );


		if(opt_verbose)
			cout << "/*\n\tResulting Sequences\n*/\n";
			
		if( sequences == NULL )
		{
			cout <<"No polymorphic sites\n";
		}
		else
			for(int i=0; i<n1+n2; i++)
				if (i<n1)
					cout <<">current["<<i<<"]\n"<<sequences[i]<<"\n";  
				else
					cout <<">old["<<i-n1<<"]\n"<<sequences[i]<<"\n";


		if(n_sp>1)
			my_population.free_population_arrays( );

		/*
			FREE trees
		*/
		my_tree->clean_whole_tree( );
		
	}



	/*
		FREE
	*/
	free(sample1);
	free(sample2);
	
	delete all_samples;
	delete all_samples_struct;
	delete serial_times;
	
	free(input_fN);

	delete my_tree;   

	return 0;
}
