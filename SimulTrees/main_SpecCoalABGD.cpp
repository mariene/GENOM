/*
	file:     main_SpecCoal.cpp
	function: Perform Simulation for s species with sequences within. Define by an entry a vector S.
	author:   <amikezor>
	date:     march 09
	modif:    sept 09 - random partitions instead of a fix sampling scheme
	modif:    nov 09  - subs Kingman by Moran
	modif:    june 10  - use the ABGD with graph partition. This version is used in the Puillandre et al. ms. Add many comments and clean the code.
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

#include "abgd.c"


#define _VERSION_ "March 2010"
#define MAX_INTRAPOP_DIFF 50         // this is an excellent prior for theta=10 ONLY

#define MAX( a, b )  (((a)>(b))?(a):(b))

char opt_verbose;



/*
	Print usage and help
*/
static void printusage( char *argv0 ){

		cerr << "Syntax is "<<argv0<<" [opt] C|Y|K|R sp_rate #sp #seq\n\n";

		cerr << "\t  M|C|Y|R             : Generate species tree using a 'C'ritical model (birth and death rates equal sp_rate).\n";
		cerr << "\t                                                      'Y'ule model (birth rate = sp_rate ; death rate=0).\n";
		cerr << "\t                                                      'M'oran model (coal rate = sp_rate*i(i-1)/2).\n";
		cerr << "\t                                                      'R'adiation model (all species originate at a time drawn exp(sp_rate))\n";

		cerr << "\t* sp_rate             : speciation rate\n";
		cerr << "\t* #sp                 : Number of species in your clade\n";
		cerr << "\t* #seq                : Number of sequences you have sampled in the clade\n";

		cerr << "\n";
		
		cerr << " [ Options ]\n";
		cerr << "\t -s #:#:#             : set sampling scheme ('2:10:4' means 3 species with respectively 2, 10 and 4 sequences)\n";
		cerr << "\t                        !! #sp and #seq arguments are then ignored. !!\n";
		cerr << "\t -N #:#:#             : set population size ('1:10:8' means that population are of size N, 10N, 8N).\n";
		cerr << "\t                        !! #sp has to coincide with the number of sizes. !!\n";
		cerr << "\t -A m|M|s             : Ancestral population sizes are either 'm'in, 'M'ax or the 's'um of the two daughter population sizes (default is min)\n";
		cerr << "\t -M #                 : set the 'M'igration rate\n";
		cerr << "\t -T #                 : set Theta, the scaled mutation rate (2pNmu, N # chr and mu mut rate /gen /locus). Default is 10.\n";
		cerr << "\t -p #                 : set ploidy (default is 1)\n";
		cerr << "\t -r #                 : number of replicates (independant runs). default is 10^4\n";
		cerr << "\t -R #                 : set recurssive partition 1/0 (default is 1 (on). set to 0 to turn it off)\n";
		cerr << "\t -x #                 : set seed for random generator (default is time())\n";
		cerr << "\t -v                   : set verbose level to 1 (default is 0)\n";
		cerr << "\t -V                   : set verbose level to 2 (default is 0)\n";
}


/*
	Print header with all selected parameters
*/
static void printheader( population *my_population, int rep, int ns, int n, char model, float sp_rate, char *sampling_scheme, char *population_sizes, char opt_recursion ){


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
	
	cout << "/* Perform simulations to generate a distance matrix ("<<_VERSION_<<")\n";
	
	cout << "/*    Species:\n";
	cout << "/*      * Model    = "<< charmodel <<"\n";
	cout << "/*      * rate     =  "<< sp_rate <<" (in Ne units)\n";
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
	
	cout << "/*    Iterations:\n";
	cout << "/*      * rep      = "<<rep<<" iterations\n";
	
	cout << "/*    ABGD:\n";
	cout << "/*      * recur.   = "<< (int)opt_recursion <<"\n";

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
	This compares the groups built by abgd to the real ones
	it returns the number of correctly delineated species
	this is \in [0, n_sp]
*/
int Compare_Comp_To_Species( struct Composante final_partition, int *ni, int n_sp ){

	int i,j,k;                 // counters. k species, i & j sequences
	int n=0;                   // number of sequences

	int sp_min,                // first and last member of a species --they are arranged in ranges--
	    sp_max;
	    
	char identical_sp;         // set to 1 when two sequences belongs to the same species.
	
	int CorrectPairs;          // # of correct pairs 
	int CorrectSpecies=0;      // # correct species -- the returned value
	
	
	for(i=0;i<n_sp;i++)       // compute the total number of sequence
		n+=ni[i];

	sp_max=-1;
	for( k = 0 ; k< n_sp ;k++){               // loop on the species
		
		if(ni[k] == 0){
			CorrectSpecies++;         // if there are no member inside a species, it is considered successfully delineated.
			continue;
		}
		
		sp_min = sp_max+1;                // this species span from sp_min to sp_max
		sp_max = sp_min+ni[k]-1;
		
		//cerr << "k " <<k <<" spmin "<<sp_min <<" sp_max "<<sp_max<< "\n";

		CorrectPairs=0;
	
		for(i=0; i<= sp_max; i++){                       // loops on i,j: for all pairs involving at least one sequence rom the current species.
			for(j=MAX(i+1,sp_min); j<n; j++){ 
			
				if(i<sp_min && j>sp_max)
					break;				
				
				if(i<sp_min || j>sp_max)        // are they of the same species ?
					identical_sp = 0;
				else
					identical_sp = 1;

				if( identical_sp == (final_partition.node_compid[i]==final_partition.node_compid[j]) )   // count success if prediction equals truth
					CorrectPairs++;
				
			}
		
		}
				
		if( CorrectPairs == (ni[k]*(ni[k]-1))/2 + ni[k]*(n-ni[k]) )    // if all pairs involving one species member is correct, this species has been correctly delineated
			CorrectSpecies++;
		

	}
	
	return CorrectSpecies;
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
	
	char SpeciationModel;          // which model 'C'is critic, 'Y' is Yule, 'R' is radiation and 'M' is moran model
	double sp_rate;                // the rate of speciation/branching in the SPECIES tree


	/*
		Sample values
	*/
	
	int n;                        // total sample size (sum of all ni)
	int **Kmatrix;                // the [n n] square matrix of parwise differences between sequences


	/*
		runs and their results
	*/
	int rep=10000;                // # of replicates
	int ndef=0;                   // number of cases where only 1 species is found.
	int *HistoCorrectSp;          // Histogram of how many successfuly delinated species. \in [0, n_sp]


	/*
		ABGD parameters
	*/
	int stability=3;                     // how many successive window size should agree
	char opt_recursion=1;                // set it to 0 if you DO NOT want to perform recussive steps to define groups
	struct Composante final_partition;   // result of ABGD. A set of groups, also named composantes of a graph (where sequences are nodes)

	double Prior=MAX_INTRAPOP_DIFF;      // what is the PRIOR limit of divergence between intra and inter species divergence.
	                                     // MAX_INTRAPOP_DIFF is excellent as long as Theta=10 (see manuscript Puillandre et al.)


	/*
		Other values
	*/
	extern Random R;              // the random generator object
	extern char abgdDEBUG;            // debuggin purpose, set to 1 it is very very verbose.
	char opt_verbose=0;             // become verbose when set to 1, option '-v'
	abgdDEBUG=0;


	/*
		Parse options
	*/

	int optc;
	extern char *optarg;
	extern int optind;

	while( (optc=getopt(argc, argv, "R:T:p:vVx:r:Ds:N:A:M:")) != -1 ){
	
		switch(optc){
		
			case 'R':
				opt_recursion = (char)atoi(optarg);
				break;

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

			case 'V':
				opt_verbose=2;
				break;

			case 'x':
				R.set_seed( atol(optarg) );
				break;

			case 'D':
				abgdDEBUG=1;
				break;
				
			case 's':
				sampling_scheme = optarg;
				break;
				
			case 'N':
				population_sizes = optarg;
				break;

			case 'M':
				migration_rate = atof(optarg);
				break;

			case 'A':
				ancestral_type = *optarg;
				if(ancestral_type != 's' && ancestral_type != 'M' && ancestral_type != 'm'){
					cerr << "Ancestral type is either 's'um, 'm'in or 'M'ax, bye";
					exit(1);
				}
				break;

		}	
	}
	
	
	if(abgdDEBUG)opt_verbose=2;

 	if(argc-optind != 4)
		printusage( argv[0] ),exit(1);	
	 		
	
	SpeciationModel = argv[optind][0];
	sp_rate 	= atof(argv[optind+1]);
	n_sp  		= atoi(argv[optind+2]);
	n     		= atoi(argv[optind+3]);




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
		Get memory for the Results
	*/
	HistoCorrectSp = (int *)calloc( (size_t)(n_sp+1), (size_t)sizeof(int) );
	if( ! HistoCorrectSp )fprintf(stderr, "main: cannot allocate HistoCorrectSp, bye\n"), exit(3);


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


	/*
		Give to user all parameters
	*/
	printheader( &my_population, rep, n_sp, n, SpeciationModel, sp_rate, sampling_scheme, population_sizes, opt_recursion );

	/*
		Run
	*/
	
	for(int r=0; r<rep; r++){
	
		//cerr << "rep: "<<r<<"\n";
		
		if(!sampling_scheme)                                      /* if the sampling scheme is not set by the user */
			draw_sampling( input_fN, sum_fN, ni, n, n_sp );
			

		if(opt_verbose)
			for(int i=0;i<n_sp;i++)
				cerr << "ni["<<i<<"]= "<<ni[i]<<"\n";
		
		
		/*
			Generate species tree
		*/
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
			default:
				cerr <<"Speciation Model undefined, bye\n";
				exit(1);
		}
		

		/*
			Sequence tree INSIDE the species tree
		*/
		
		my_tree->isolation_coalescent(ni);



		if(opt_verbose == 2){
			my_population.print_isolations();
			my_tree->printout_tree();
		}

		
		
		/*
			Add mutation and generate [n, n] pairwise differences matrix 
		*/
		my_tree->mutate_whole_tree( );
		Kmatrix = my_tree->compute_Kmatrix();
	
		/*
			Run ABGD to find the partition
		*/
		final_partition = Find_ABGD_Partition( Kmatrix, n, stability, Prior , opt_recursion, opt_verbose );
		
		if(final_partition.nc == 1)              // count how many are a single species
			ndef++;
	
		if(opt_verbose == 2)
			cout << "final_partition.nc: "<< final_partition.nc<<"\n";


		/*
			Compare ABGD partition to true species and return the number of correctly predicted speciess
		*/
		int CorrectSpecies = Compare_Comp_To_Species( final_partition, ni, n_sp );		

		
		switch(opt_verbose)
		{
			case 2:
				for(int i=0; i<n; i++, printf("\n"))
					for(int j=0; j<=i; j++)
						cout << Kmatrix[i][j] << "\t";
			case 1:
				my_population.print_isolations( );
				my_tree->printout_treeWithMut(  );
				my_tree->printout_treeFromMut(  );
				cout << "CorrectSpecies " << CorrectSpecies << "\n"; 
		}

		
		HistoCorrectSp[CorrectSpecies]++;
			
		/*
			FREE Kmatrix
		*/		
		for(int i=1;i<n;i++)
			delete Kmatrix[i];
		delete Kmatrix;
		
		/*
			FREE trees
		*/
		my_tree->clean_whole_tree( );
		my_population.free_population_arrays( );
		free_composante(final_partition);
		

	}


	/*
		Print OUT results
	*/
	
	int HistCumul=0;
	printf("DistribCorrectSpecies: ");

	for( int i=0; i<=n_sp;i++ ){
		HistCumul += HistoCorrectSp[ i ];
		printf("%.3f\t", HistCumul/(float)rep);
	}
	printf("\n");



	/*
		FREE
	*/
	free(HistoCorrectSp);
	free(ni);
	free(input_fN);

	delete my_tree;   

	return 0;
}









/*
	The 2 following functions are deprecated --defined but NOT used
*/

float Compare_Kmin_Sp_Pairs( int *Kdistrib, float Kmin, int *ni, int n_sp ){

	int i,j, k;

	int n=0;                   // number of sequences
	
	float f_CorrectPairs=0;    // fractions of pairs that are correctly assigned (true identical species or true different species)
	
	
	for(i=0;i<n_sp;i++)
		n+=ni[i];
	
	for(i=0;i<n;i++)
		for(j=i+1;j<n;j++){                                                   // for all pairs of sequences
		
			int Kij = Kdistrib[  i*n - (i*(i+1))/2 + (j-i-1) ];           // 1. get dist in mutations
			 
			char identical_sp=0;                                          // 2. get status (1=same species ; 0=diff species)
			int nc=0;                                                     // # of current sequences since the begining (species by species)

			for(k=0; identical_sp==0 && k<n_sp; nc+=ni[k], k++)
				if( (i>=nc && i<nc+ni[k]) && (j>=nc && j<nc+ni[k]) )  // if numbers are within the same species
					identical_sp=1;
		
		
			if( identical_sp == (Kij<Kmin) )                                // check if predited status and true status equate
				f_CorrectPairs +=1.00;
		
		}
	
	return f_CorrectPairs/ (float)( (n*(n-1))/2 );
}

int Compare_Kmin_Sp_Species( int *Kdistrib, float Kmin, int *ni, int n_sp ){

	int i,j,k;                 // counters. k species, i & j sequences
	int n=0;                   // number of sequences

	int Kij;                   // distance between sequence i and j

	int sp_min, sp_max;
	char identical_sp;
	
	int CorrectPairs;          // # of correct pairs (ident sp with Kij<Kmin + diff sp with Kij >= Kmin)
	int CorrectSpecies=0;
	
	
	for(i=0;i<n_sp;i++)
		n+=ni[i];

	sp_max=-1;
	for( k = 0 ; k< n_sp ;k++){
		
		if(ni[k] == 0){
			CorrectSpecies++;
			continue;
		}
		
		sp_min = sp_max+1;
		sp_max = sp_min+ni[k]-1;
		
		//cerr << "k " <<k <<" spmin "<<sp_min <<" sp_max "<<sp_max<< "\n";

		CorrectPairs=0;
	
		for(i=0; i<= sp_max; i++){
			
			for(j=MAX(i+1,sp_min); j<n; j++){
			
				if(i<sp_min && j>sp_max)
					break;
				
				Kij = Kdistrib[  i*n - (i*(i+1))/2 + (j-i-1) ];           // 1. get dist in mutations
				
				
				if(i<sp_min || j>sp_max)
					identical_sp = 0;
				else
					identical_sp = 1;
				
				
				if( identical_sp == (Kij<Kmin || Kmin==-1) )
					CorrectPairs++;
				//else
				//	cerr <<"wrong pair "<<i<<" vs "<<j<<"\n";
				
				
			}
		
		}
		
		//cerr <<"CorrectPairs "<<CorrectPairs<<"\n";
		
		if( CorrectPairs == (ni[k]*(ni[k]-1))/2 + ni[k]*(n-ni[k]) )
			CorrectSpecies++;
		

	}
	
	return CorrectSpecies;
}



