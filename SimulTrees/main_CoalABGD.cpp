/*
	file:     main_CoalPGD.cpp
	function: Perform Simulation for 1 species with sequences within.
	author:   <amikezor>
	date:     sept 09
	modif:    mar 2010
*/


#include <iostream>
#include <string>
using namespace std;

#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "population.h"
#include "sample.h"
#include "random.h"
#include "coalescent_tree_diversity.h"
#include "distribution.h"
#include "readfile.h"

//#include "abgd.c"

#define _VERSION_ "Sept 09"

#define SIGN( a ) ( ( (a) > 0 )?1: ( ((a)==0)?0:-1)  )
static int Increase(const void *v1, const void *v2){  	return (int)SIGN( *((double *)v1) - *((double *)v2));  };
#undef SIGN

char DEBUG;

char opt_verbose;



struct Zscores {

	double zmax;
	double zmax_int;

};

static void printusage( char *argv0 ){

		cerr << "Syntax is "<<argv0<<" [opt] #seq\n\n";

		cerr << "\t* #seq                : Number of sequences you have sampled in the clade\n";

		cerr << "\n";
		
		cerr << " [ Options ]\n";
		cerr << "\t -T #                 : set Theta, the scaled mutation rate (2pNmu, N # chr and mu mut rate /gen /locus). Default is 10.\n";
		cerr << "\t -C 0.#               : set 'C'onfidence Interval (2-sided).\n";
		cerr << "\t -p #                 : set ploidy (default is 1)\n";
		cerr << "\t -r #                 : number of replicates (independant runs). default is 10^4\n";
		cerr << "\t -x #                 : set seed for random generator (default is time())\n";
		cerr << "\t -w #                 : winsiz min (default is 2)\n";
		cerr << "\t -W #                 : winsiz max (default is 100)\n";
		cerr << "\t -s #                 : set species number (independant replicates; default is 1)\n";
		cerr << "\t -v                   : set verbose level to 1 (default is 0)\n";
		cerr << "\t -V                   : set verbose level to 2 (default is 0)\n";
}



static void printheader( population *my_population, int rep, int n, float CI ){


	extern Random R;

	
	
	cout << "/* Perform simulations to generate a distance matrix ("<<_VERSION_<<")\n";
	
	cout << "/*    Population:\n";
	cout << "/*      * Theta   = "<< my_population->get_Theta() <<" (2.Ne.p.mu)\n";
	cout << "/*      * Ploidy  = "<< (int)my_population->get_ploidy() <<" \n";

	cout << "/*    Samples:\n";
	cout << "/*      * #seq    = "<< n <<"\n";

	cout << "/*    Confidence Interval:\n";
	cout << "/*      * CI one-tail = "<<CI<<"\n";

	
	cout << "/*    Iterations:\n";
	cout << "/*      * rep     = "<<rep<<" iterations\n";
	
	cout << "/*    seed random generator:\n";
	cout << "/*      * seed    = "<<R.get_seed()<<"\n";

	cout << "\n";
}


double Find_Val_of_MaxSlope(  int **Array, long n, int winsiz, int nsp ){
	
	int i=0,j=0,c;
	
	double *ValArray;        /* a double version of the values array */
	long N=nsp*n;            /* Total number of distance -- ValArray size */
	
	double  *Slope;           /* Slope <=> "derivative" of array */

	int  top=0;        
		
	int r,l;
	double Mean_dist;		
	float q;


	/*
		1. Create a sorted array of all distance: ValArray
	*/
	ValArray = ( double * ) malloc( N*sizeof(double)  );
	if(!ValArray)fprintf(stderr, "cannot allocate ValArry, bye\n"), exit(3);
	
	Slope = ( double * ) malloc( (N-winsiz+1)*sizeof(double)  );
	if(!Slope)fprintf(stderr, "cannot allocate ValArry, bye\n"), exit(3);


	for(i=0;i<nsp;i++)
		for(j=0;j<n;j++)
			ValArray[i*n+j] = (double)Array[i][j];
	
	qsort((void *) ValArray, (size_t) N, (size_t) sizeof(double), Increase );


	/*
		2. Compute the Slope Max
	*/
	
	top=0;
	for(i=0; i <= N-winsiz ; i++){
		Slope[i] = (ValArray[i+winsiz-1]-ValArray[i])/(double)(winsiz-1);
		top = (Slope[i]>Slope[top])?i:top;
	}

	r=l=top;
	Mean_dist=0;
	for(c=0, q=0.95; q>=0.20; q-=0.01){
	
		while( r<N-winsiz && Slope[r] > q*Slope[top] && Slope[r]<=Slope[top] ){
			r++;
		}
		
		while( l>=0 && Slope[l] > q*Slope[top] && Slope[l]<=Slope[top] ){
			l--;
		}
				
		Mean_dist += (ValArray[r+winsiz-1]+ValArray[l])/2;           /* Mean_dist is the sum of all mean values */
		c++;
			
	}

	if(DEBUG)printf("   final peak exploration == l: %d / top %d / r: %d / c: %d \n", l, top, r, c);

	Mean_dist /= c;       /* From sum to mean */


	if(DEBUG)printf("    Mean: dist %f \n", Mean_dist);

	
	free(ValArray);
	free(Slope);

	return Slope[top];
}



double get_MaxVal( int ws, int n, population *pop, float CI, long rep, int nsp){

	int **Kdistrib;                // the distribution of S between sequences
	
	coalescent_tree_diversity **ptrees;

	distribution *distrib_Val;


	double HighCI;
	double tmp;

	ptrees = (coalescent_tree_diversity **)malloc( (size_t) nsp*sizeof(coalescent_tree_diversity *) );
	if( ! ptrees )fprintf(stderr, "get_MaxVal: cannot allocate ptrees, bye"), exit(3);

	Kdistrib = (int **)malloc( (size_t) nsp*sizeof(int *) );
	if( ! Kdistrib )fprintf(stderr, "get_MaxVal: cannot allocate Kdistrib, bye"), exit(3);
	
	for(int i=0;i<nsp;i++)
		ptrees[i] = new coalescent_tree_diversity( n, pop );


	distrib_Val = new distribution(rep, 1.0-2.0*(1.0-CI) );                                                  // it is a one-tail CI that we are interested in
	if( ! distrib_Val  )fprintf(stderr, "get_MaxVal: cannot allocate distrib_Val, bye"), exit(3);

	for(int r=0; r<rep; r++){
	
		
		for(int i=0;i<nsp;i++){
			ptrees[i]->standard_coalescent();
			ptrees[i]->mutate_whole_tree( );                                                           // Compute distances (Sites/Mutations)
			Kdistrib[i] = ptrees[i]->compute_Kdistrib();
		}


		tmp = Find_Val_of_MaxSlope(  Kdistrib, (long)(n*(n-1))/2, ws, nsp );
		
		if(tmp == -1)
			r--;
		else
			distrib_Val->add_value( tmp );           // Find MaxSlope Increase

		
		for(int i=0;i<nsp;i++){
			delete Kdistrib[i];
			ptrees[i]->clean_whole_tree( );
		}
		
		
	}
	
	delete distrib_Val;
	
	HighCI  = distrib_Val->get_high_boundary();
	
	free(ptrees);
	free(Kdistrib);


	return HighCI;		
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
	
	/*
		For isolation event new version
	*/

	int n;                        //  sample size

	double Theta;                 // Population Mutation rate 2Nmu -> mutation rate on a branch is Theta/2
	
	int rep;                      // # of replicates


	extern Random R;              // the random generator object

	extern char DEBUG;


	float CI=0.950;              // confidence interval of the value
	

	long ws_min=2, ws_max=100;
	
	int nsp=1;

	/*
		Default values
	*/
	
	Theta=10;
	ploidy=1;
	rep=10000;

	opt_verbose=0;
	DEBUG=0;


	/*
		parse input
	*/
	while( (optc=getopt(argc, argv, "T:p:vVx:r:DC:W:w:s:")) != -1 ){
	
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

			case 'V':
				opt_verbose=2;
				break;

			case 'x':
				R.set_seed( atol(optarg) );
				break;

			case 'D':
				DEBUG=1;
				break;

			case 'C':
				CI=atof( optarg );
				break;
				
			case 'w':
				ws_min=atol(optarg);
				break;
				
			case 'W':
				ws_max=atol(optarg);
				break;

			case 's':
				nsp=atoi(optarg);
				break;

		}	
	}
	
 	if(argc-optind != 1)
		printusage( argv[0] ),exit(1);	
	 		
	
	n = atoi(argv[optind]);


	/*
		Set Pop, built sp tree
	*/

	population my_population;                                    // with Ne=1, mu=1 and haploid;
	my_population.set_Theta( Theta );
	my_population.set_ploidy( ploidy );
	


	/*
		Init the sequences samples
	*/
	


	printheader( &my_population, rep,  n, CI);	
	
	
	for(int w_size=ws_min;w_size<=ws_max && w_size<(n*(n-1)/2);w_size+=1){

		printf("win_siz %d --> SlopeMax %f\n", w_size, get_MaxVal( w_size, n, &my_population, CI, rep, nsp) );

	}

//	printf("Theta %f ; HighBound %f\n", Theta, get_MaxVal(  n, &my_population, CI, rep) );


	

	return 0;
}
