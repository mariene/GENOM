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


#define _VERSION_ "Jan 2013"

#define MAX( a, b )  (((a)>(b))?(a):(b))

char opt_verbose;

/*
	Print usage and help
*/
static void printusage( char *argv0 ){

	cerr << "Syntax is:\n";
	cerr <<  "\t"<< argv0<<" [ options ] sp_rate #:#:#   (given assignment)\n";
	cerr << "  OR\t"<<argv0<<" [ options ] sp_rate #species #sequences  (random assignment)\n\n";

	cerr << "\t* sp_rate             : speciation rate\n";
	cerr << "\t* #:#:#...            : sampling_scheme ('2:10:4:1' translates to 4 species with respectively 2, 10, 4 and 1 sequences)\n";
	cerr << "\t* #species            : Number of species in your clade\n";
	cerr << "\t* #sequences          : Number of sequences you have sampled in the clade\n";

	cerr << "\n";
	
	cerr << " [ options ]\n";
	cerr << "\t -S  M|C|Y|R|F        : Generate species tree using a 'M'oran model -default- (coal rate = sp_rate*i(i-1)/2).\n";
	cerr << "\t                                                      'Y'ule model (birth rate = sp_rate ; death rate=0).\n";
	cerr << "\t                                                      'C'ritical model (birth and death rates equal sp_rate).\n";
	cerr << "\t                                                      'R'adiation model (all species originate at a time drawn exp(sp_rate))\n";
	cerr << "\t                                                      'F'ixed model (rate is then read as a fixed time, and speciation occured at t, 2t, 3t, ...)\n";
	
	cerr << "\t -N #:#:#             : set population size ('1:10:8' means that population are of size N, 10N, 8N). Default all population are 1\n";
	cerr << "\t                        !! it has to coincide with the number of species sampled. !!\n";
	cerr << "\t -A m|M|s             : Ancestral populations sizes are either 'm'in, 'M'ax or the 's'um of the two daughter population sizes (default is min)\n";
	cerr << "\t -M #                 : set the 'M'igration rate\n";
	cerr << "\t -T #                 : set Theta, the scaled mutation rate (2pNmu, N # chr and mu mut rate /gen /locus). Default is 10.\n";
	cerr << "\t -p #                 : set ploidy (default is 1)\n";
	cerr << "\t -l #                 : sequence is of length l (default is infinite site model)\n";
	cerr << "\t -L #                 : sequence is the concatenation of L independent loci\n";
	cerr << "\t -x #                 : set seed for random generator (default is time())\n";
	cerr << "\t -q                   : 'q'uiet mode = only output sequences\n";
}


/*
	Print header with all selected parameters
*/
static void printheader( population *my_population, int rep, int ns, int n, char model, float sp_rate, int nLocus, char *sampling_scheme, char *population_sizes ){


	const char *charmodel=NULL;
	extern Random R;

	
	switch(model){
	
		case 'C':
			charmodel="Critic";
			break;
		case 'M':
			charmodel="Moran";
			break;
		case 'Y':
			charmodel="Yule";
			break;
		case 'R':
			charmodel="Radiation";
			break;
		default:
			charmodel="undef";
	}
	
	cout << "/* Perform simulations to generate homologous sequence with species ("<<_VERSION_<<")\n";
	
	cout << "/*    Species:\n";
	if(ns>1){
		cout << "/*      * Model    = "<< charmodel <<"\n";
		cout << "/*      * rate     = "<< sp_rate <<" (in Ne units)\n";
	}
	cout << "/*      * #sp      = "<< ns <<"\n";

	cout << "/*    Populations:\n";
	cout << "/*      * Theta    = "<< my_population->get_Theta() <<" (2.Ne.p.mu)\n";
	cout << "/*      * Ploidy   = "<< (int)my_population->get_ploidy() <<" \n";
	if(population_sizes != NULL)
	cout << "/*      * Pop_size = "<< population_sizes <<"\n";

	cout << "/*    Samples:\n";
	cout << "/*      * #seq     = "<< n <<"\n";
	if(sampling_scheme != NULL)
	cout << "/*      * sampling = "<< sampling_scheme <<"\n";
	cout << "/*      * locus = "<< nLocus <<"\n";
	
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
	IF input_fN is NULL, draw a sample uniformly
	OTHERWISE, assign the squences to the sample respecting the different population sizes
	
	The filled array is ni. It is first set to 0 and then increased one by one.
	
*/
void draw_sampling( double *input_fN, double sum_fN, int *ni, int n, int n_sp ){

	
	int i,j;                      // counters

	double  ran_doubl,            // [0,1[ random number
		cumul_fN;             // used to see what fraction of probabilities are in each population
	
	extern Random R;              // the random generator object

	for(i=0;i<n_sp;i++)                        /* draw it */
		ni[i]=0;
	
	if( input_fN==NULL){

		for(i=0;i<n;i++)
			ni[ R.uniform_int(0, n_sp-1) ]++;
			
	}
	else{
	
		for( i=0;i<n;i++){
	
			ran_doubl = R.uniform_dev();
			cumul_fN=0;
			
			for(j=0; j<n_sp; j++){
				
				cumul_fN += (input_fN[j] / sum_fN);

				if( ran_doubl < cumul_fN )
					break;
			}
			if( j == n_sp ){
				cerr << "draw_sampling: error in setting the sequences in species\n";
				exit(6);
			}
			else{
				ni[ j ] ++;
			}
		}
	}

}


/*
	Main funtion
*/
int main(int argc, char *argv[]){
	

	
	/*
		Population(s) parameters
	*/

	double migration_rate=0;      // migration rate 2Nm, FOR NOW always set to 0

	char ploidy=1;                 // 1 for haploid, 2 for diploid
	char *sampling_scheme=NULL;    // if set to NULL, species sampling is randomnly done, if not respect the sampling scheme
	
	int *ni;                      // number of loci in today's species. Sum is n.
	int n_sp=0;                   // number of species
	
	char *population_sizes=NULL;   // if set to NULL all population have the same size, otherwise, they can be of different size (e.g. "1:0.1:0.01").
	double *input_fN=NULL;         // population_sizes string turn into an array of doubles.
	double sum_fN=0;               // sum of all pop sizes
	
	char ancestral_type='m';       // if 'm' ancestral population size is the minimum of both desc, if 'M' it is
	                               // the maximum and if 's', it is the sum. Any other value is an error.

	double Theta=10;              // Population Mutation rate 2Nmu -> mutation rate on a branch is Theta/2

	/*
		Speciation model
	*/	
	double birth=-1,
	       death=-1;               // rate of birth and death 
	
	char SpeciationModel='M';          // which model 'C'is critic, 'Y' is Yule, 'R' is radiation and 'M' is moran model
	double sp_rate;                // the rate of speciation/branching in the SPECIES tree


	/*
		Sample values
	*/
	
	int n=0;                        // total sample size (sum of all ni)
	int l=0;
	
	int L=1;               // number of concatenated loci

	char **sequences=NULL;


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

	while( (optc=getopt(argc, argv, "T:p:x:qr:DN:A:M:l:S:L:")) != -1 ){
	
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

			case 'l':
				l = atoi(optarg);
				break;

			case 'A':
				ancestral_type = *optarg;
				if(ancestral_type != 's' && ancestral_type != 'M' && ancestral_type != 'm'){
					cerr << "Ancestral type is either 's'um, 'm'in or 'M'ax, bye\n";
					exit(1);
				}
				break;


			case 'S':
				SpeciationModel = *optarg;
				if(SpeciationModel != 'M' && SpeciationModel != 'C' && SpeciationModel != 'Y' && SpeciationModel != 'R'  && SpeciationModel != 'F'){
					cerr << "SpeciationModel is either 'M', moran, 'C', critical, 'Y', yule, 'R', radiation or 'F', fixed times, bye\n";
					exit(1);
				}
				break;

			case 'L':
				L = atoi(optarg);
				break;

		}	
	}
	
	

 	if(argc-optind != 3 && argc-optind != 2)
		printusage( argv[0] ),exit(1);	
	 		
	
	sp_rate 	= atof(argv[optind]);

	if( argc-optind == 3){
		
		n_sp  		= atoi(argv[optind+1]);
		n     		= atoi(argv[optind+2]);

	}
	else
	{
		sampling_scheme = argv[optind+1];
	}



	/*
		Parse the sampling scheme, if any
		otherwise get memory for the random assignation
	*/
	if (sampling_scheme == NULL){
		ni = (int *)calloc( (size_t) n_sp, (size_t)sizeof(int) );
		if(!ni )fprintf(stderr, "cannot allocate ni, bye\n"),exit(3);
	}
	else{
		int i;
		ni = parse_sparg(  sampling_scheme , &n_sp );
		for(n=0, i=0;i<n_sp;n+=ni[i++]);
	}

	
	/*
		If population sizes are not all equal, fill the array
	*/	
	if(population_sizes){

		int tmp;
		input_fN = parse_sparg_f( population_sizes , &tmp );
		if(tmp != n_sp){
			cerr << "the species size are not compatible with the number of species, bye\n";
			exit(1);
		}
		
		for(int i=0;i<tmp;i++)
			sum_fN+=input_fN[i];
	}

	


	/*
		Set Pop, built sp tree
	*/

	population my_population;                                    // with Ne=1, mu=1 and haploid;
	my_population.set_Theta( Theta / (L+0.0) );                  // theta is evenly distributed among all sequences
	my_population.set_ploidy( ploidy );

	/*
		Init the sequences samples
	*/
	coalescent_tree_diversity *my_tree = new coalescent_tree_diversity( n, &my_population );


	/*
		Give to user all parameters
	*/
	printheader( &my_population, rep, n_sp, n, SpeciationModel, sp_rate, L, sampling_scheme, population_sizes );

	/*
		Run the simulations
	*/
	
	for(int r=0; r<rep; r++){
	
		//cerr << "rep: "<<r<<"\n";
		
		if(!sampling_scheme)                                      /* if the sampling scheme is not set by the user */
			draw_sampling( input_fN, sum_fN, ni, n, n_sp );
			

		if(opt_verbose){
			cout << "/*\n\tHow many sequences in each species\n*/\n";
			for(int i=0;i<n_sp;i++)
				if(i==n_sp-1)
					cout <<ni[i]<<"\n";
				else
					cout <<ni[i]<<":";
		}
		
		/*
			Generate species tree
		*/
		if(n_sp>1)
			switch(SpeciationModel){
		
				case 'C':
					death=sp_rate;
				case 'Y':
					death=(death==-1)?0:death;
					birth=sp_rate;
					my_population.BirthDeath_species_tree(  n_sp, migration_rate, birth, death, input_fN, ancestral_type );
					break;
				case 'M':
					my_population.Moran_species_tree(  n_sp, migration_rate, sp_rate, input_fN, ancestral_type );
					break;
				case 'R':
					my_population.Radiation_species_tree(  n_sp, migration_rate, sp_rate, input_fN, ancestral_type );
					break;
				case 'F':
					my_population.Fixed_species_tree(  n_sp, migration_rate, sp_rate, input_fN, ancestral_type );
					break;
				default:
					cerr <<"Speciation Model undefined, bye\n";
					exit(1);
			}
	
		if(n_sp>1 && opt_verbose>1){
			cout << "/*\n\tSpecies tree\n*/\n";
			my_population.print_isolations();
		}
		
		for(int x=0; x<L; x++){
		
		
			if(migration_rate>0 && n_sp>1 && opt_verbose)
				cout << "/*\n\tMigration Event(s)\n*/\n";

			/*
				Sequence tree INSIDE the species tree
			*/
			if(n_sp>1)
				my_tree->isolation_coalescent(ni);
			else
				my_tree->standard_coalescent();

			if(opt_verbose){
				cout << "/*\n\tCoalescent tree (locus "<< x+1 <<")\n*/\n";
				my_tree->printout_tree();
			}
				
			/*
				Add mutations 
			*/
			if( x==L-1 )
				my_tree->mutate_whole_tree_seq( l - l/L*x );   /* add all the remaining sites */
			else
				my_tree->mutate_whole_tree_seq( l/L );
				
			
			sequences = my_tree->merge_sequences( sequences );
		
			/*
				FREE trees
			*/
			my_tree->clean_whole_tree( );
	

		}

		if(n_sp>1)
			my_population.free_population_arrays( );


		
		cout << "/*\n\tResulting Sequences\n*/\n";
		if( sequences == NULL )
		{
			cout <<"No polymorphic sites\n";
		}
		else
			for(int i=0; i<n; i++)
				cout <<">seq["<<i<<"]\n"<<sequences[i]<<"\n";     
		
	}



	/*
		FREE
	*/
	free(ni);
	free(input_fN);

	delete my_tree;   

	return 0;
}
