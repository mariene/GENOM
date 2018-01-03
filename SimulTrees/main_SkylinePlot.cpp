/*
	file:     main_SkylinePlot.cpp
	function: Perform simulations and output the Skyline Plot
	author:   <gachaz>
	date:     July 2014
*/


#include <iostream>
using namespace std;

#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "population.h"
#include "sample.h"
#include "random.h"
#include "coalescent_tree_diversity.h"
#include "distribution.h"


#define _VERSION_ "Oct 08"

#define KRONECKER( A, B ) (( (A) == (B) )?1:0)


char opt_verbose;



static void printusage( char *argv0 ){

		cerr << "Syntax is "<<argv0<<" [-h] [ options ] [ scenario ] n rep\n";
}

static void printhelp( char *argv0 ){

		printusage( argv0 );

		cerr << "\n[core values]\n";
		cerr << "\tn                   : total sample size.\n";
		cerr << "\trep                 : number of replicates (independant trees).\n";


		cerr << "\n[info]\n";
		cerr << "\t-v                  : be verbose\n";
		cerr << "\t-V                  : be extra verbose\n\n";
		cerr << "\t-D                  : debugging verbosity\n\n";

		cerr << "\n[option]\n";
		cerr << "\t-p 1|2              : ploidy - either 1 -haploid- or 2 -diploid- (default 1)\n";
		cerr << "\t-C #                : Confidence Interval (default 0.95)\n";
		cerr << "\t-x #                : set the seed for random numbers to # (otherwise time() is used)\n";
		
		cerr << "\n[scenario]\n!! At most, a single scenario is chosen !!\n";
		
		cerr << "\t-s \"Ts p R alpha\"       : perform a 's'weep using the parameters Ts time bwd to sweep, p freq at the sweep end\n"
		        "\t                          R recomb rate (p*2Nr) and alpha selection strength (p*Ns).\n"
		        "\t                          when p<1 and Ts=0, the sweep is ongoing. When Ts>0, p is set to 1\n";
		
		cerr << "\t-b \"Tb Tl f\"            : perform a 'b'ottleneck using the parameters Tb time bwd to bottleneck,\n"
		        "\t                          Tl time of the event and f bottleneck relative size.\n";

		cerr << "\t-e Gr                   : perform an 'e'xponential growth using the parameter Gr for growth rate.\n";

		cerr << "\t-l Gr                   : perform an 'l'inear growth using the parameter Gr for growth rate.\n";
		
		cerr << "\t-t MIN_trmca            : exclude any run where the TMRCA is smaller than MIN_trmca.\n";

		cerr << "\t-i \"Ti n0 f[0] f[1] M\"  : perform a 'i'solation with migration with 2 pop using the parameters Ti time bwd to isolation/speciation,\n"
		        "\t                          n0 the number of sample in species 0 and f[i] size of species i relative to the ancestral one\n"
			"\t                          and M the 2Nm migration rate (symetric for now).\n";

		cerr << "\t-n k                      : among the n sampled individual, report the spectrum within k (and n-k) individuals,\n"
		        "\t                          it mimics 'n'ested mutation, conditionned on the founder one to be at k in the sample\n";

		cerr << "\n";
}



static void printheader( population *my_population, int n , int *ni, int rep, int k, char chosen_scenario ){


	cout << "/* Perform simulations to get Freq and theta spectrum ("<<_VERSION_<<")\n";
	cout << "/*    Population:\n";
	cout << "/*      * ploidy = "<< 1*my_population->get_ploidy() <<" \n";
	cout << "/*    Sample:\n";
	cout << "/*      * n     = "<<n<<" loci\n";
	cout << "/*    Iterations:\n";
	cout << "/*      * rep   = "<<rep<<" iterations\n";
	
	switch(chosen_scenario){
	
		case 's':
			cout << "/*    [sweep]\n";
			cout << "/*      * Ts     = "<<my_population->get_sweep_TS()<<" (in "<<my_population->get_ploidy()<<"N generations)\n";
			cout << "/*      * alpha  = "<<my_population->get_sweep_alpha()<<" (select coef = "<<my_population->get_ploidy()<<"*Ns)\n";
			cout << "/*      * p      = "<<my_population->get_sweep_p()<<" (frequency of the selected allele at time 0)\n";
			cout << "/*      * R      = "<<my_population->get_sweep_R()<<" (recomb coef = "<<my_population->get_ploidy()<<"*2Nr)\n";
			break;
		
		case 'i':
			cout << "/*    [isolation]\n";
			cout << "/*      * Ti     = "<<my_population->get_isolation_Ti(0).time<<"\t(in "<< (int)my_population->get_ploidy() << "N generations)\n";
			cout << "/*      * fN[0]   = "<<my_population->get_isolation_fN(0)<<"\t(real size is "<< 1*my_population->get_ploidy() <<".f[0].N)\n";
			cout << "/*      * fN[1]   = "<<my_population->get_isolation_fN(1)<<"\t(real size is "<< 1*my_population->get_ploidy() <<".f[1].N)\n";
			cout << "/*      * M      = "<<my_population->get_isolation_M()<<" (in 2Nm units)\n";
			cout << "/*    Sample:\n";
			cout << "/*      * ni[0]  = "<<ni[0]<<" loci for species 0\n";
			cout << "/*      * ni[1]  = "<<ni[1]<<" loci for species 1\n";
			break;
		
		case 'b':
			cout << "/*     [bottleneck]\n";
			cout << "/*      * Tb     = "<<my_population->get_bottleneck_Tb()<<"\t(in "<< (int)my_population->get_ploidy()<<"N gen)\n";
			cout << "/*      * Tl/f   = "<<my_population->get_bottleneck_Tl()/my_population->get_bottleneck_f()<<"\t(strength of the bottleneck in coalescent time)\n";
			cout << "/*        * Tl   = "<<my_population->get_bottleneck_Tl()<<"\t(in "<< (int)my_population->get_ploidy()<<"N gen)\n";
			cout << "/*        * f    = "<<my_population->get_bottleneck_f()<<"\t(population at bottleneck is f."<<1*my_population->get_ploidy()<<".N)\n";
			break;
			
		case 'e':
			cout << "/*     [exp growth]\n";
			cout << "/*      * Gr     = "<<my_population->get_exponential_growth_Gr()<<"\n"
			        "                   (population multiplies forward "<< exp(my_population->get_exponential_growth_Gr())  <<" every "<< (int)my_population->get_ploidy()<<"Ne gen.)\n";
			break;
			
		case 'l':
			cout << "/*     [linear growth]\n";
			cout << "/*      * slope  = "<<my_population->get_linear_growth_Gr() <<"\n"
			        "                   (population adds forward "<< my_population->get_linear_growth_Gr()  <<" every "<< (int)my_population->get_ploidy()<<"Ne gen.)\n";
			break;
	
		case 'n':
			cout << "/*     [nested mutations]\n";
			cout << "/*      * k     = "<<k<<"\n";
			break;
	
		default:
			;
	}
	

	extern Random R;
	
	cout << "/*    seed random generator:\n";
	cout << "/*      * seed  = "<<R.get_seed()<<"\n";

	cout << "\n";
}


struct bin
{

	float sum;
	int n;
};


int main(int argc, char *argv[]){
	
	/*
		Parse options
	*/
	
	int optc;
	extern char *optarg;
	extern int optind;

	char ploidy=1;          // 1 for haploid, 2 for diploid
	
	char opt_verbose=0;     // print summary statistics of each rep
	
	char opt_scenario=0;    // it becomes 's' for sweep, 'S' for speciation, 'b' for bottleneck and 'e' for exponential growth, 'n' for nested mutations


	double TS, Rec, alpha, p;    // for a sweep
	double Tb, Tl, f;            // for a bottleneck
	double Gr=0;                 // for an exp/linear growth
	int k=0;                     // for nested mutations, size of the founder mutation

	struct isolation_event Ti[1];  // in struct, there are the time and which species are merged. Here only species 0 and 1 merge
	double TI;
	double fN[3];                  // relative size of each species. 0->(ns-1) today's species. >=ns: anceestral ones and 2ns-2 the FIRST species
	double M;                      // migration rate 2Nm
	int ni[3];                    // number of loci in today's species. Sum is n.
		
	
	/*
		Parse input
	*/
	

	extern char debug;
	extern Random R;
	
	float max_TMRCA=0;

	float binsize = 0.1;

	debug=0;

	while( (optc=getopt(argc, argv, "hp:vVs:b:i:e:x:uDn:t:l:B:")) != -1 ){
	
		switch(optc){
		
			case 'v':
				opt_verbose = 1;
				break;
			
			case 'D':
				debug = 1;
				break;
			
			case 'V':
				opt_verbose = 2;
				break;
			
			case 'p':
				ploidy = (char)atoi( optarg );
				if(ploidy != 1 && ploidy != 2)fprintf(stderr, "ploidy can only be set to 1 or 2\n"), exit(1);
				break;
				
			
			/*
				Optionnal scenario of non-regular coalescent
			*/
				
			case 's':
				if(opt_scenario)fprintf(stderr, "you cannot set several scenario at the same time\n"), exit(1);
				opt_scenario = 's';
				if(  sscanf(optarg, "%lf%lf%lf%lf", &TS, &p, &Rec, &alpha) != 4 )
					fprintf(stderr, "cannot parse sweep argument, bye\n"), exit(1);
				break;
				
			case 'b':
				if(opt_scenario)fprintf(stderr, "you cannot set several scenario at the same time\n"), exit(1);
				opt_scenario = 'b';
				if(  sscanf(optarg, "%lf%lf%lf", &Tb, &Tl, &f) != 3 )
					fprintf(stderr, "cannot parse bottleneck argument, bye\n"), exit(1);
				break;

			case 'i':
				if(opt_scenario)fprintf(stderr, "you cannot set several scenario at the same time\n"), exit(1);
				opt_scenario = 'i';
				if(  sscanf(optarg, "%lf%d%lf%lf%lf", &TI, ni, fN, fN+1, &M) != 5 )
					fprintf(stderr, "cannot parse speciation argument, bye\n"), exit(1);
				break;
				
			case 'e':
				if(opt_scenario)fprintf(stderr, "you cannot set several scenario at the same time\n"), exit(1);
				opt_scenario = 'e';
				Gr = atof(optarg);
				break;
				
			case 'l':
				if(opt_scenario)fprintf(stderr, "you cannot set several scenario at the same time\n"), exit(1);
				opt_scenario = 'l';
				Gr = atof(optarg);
				break;

			case 'n':
				if(opt_scenario)fprintf(stderr, "you cannot set several scenario at the same time\n"), exit(1);
				opt_scenario = 'n';
				k = atoi(optarg);
				break;

			case 'x':
				R.set_seed( atol(optarg) );
				break;

			case 'h':
				printhelp( argv[0] );
				exit(1);
				break;
						
			case 't':
				max_TMRCA = atof( optarg );
				break;
			
			case 'B':
				binsize = atof( optarg );
				break;
			
			
		}	
	}
	
 	if(argc-optind != 2)
		printusage( argv[0] ),exit(1);	
	 	
	int n         = atoi(argv[optind]);
	int rep       = atoi(argv[optind+1]);
		


	/*
		Set the population parameters
	*/
	
	population my_population;                                    // with Ne=1, mu=1 and haploid;
	my_population.set_ploidy( ploidy );
	 	
	
	int nbins = 2/binsize +1;
	struct bin * my_bins = new struct bin [ nbins ];
	
	for(int i=0;i<nbins;i++){
		my_bins[i].sum=0.0;
		my_bins[i].n=0;
	}
	
	switch(opt_scenario){
	
		case 's':
			my_population.set_sweep( TS, p, alpha, Rec );
			break;
		case 'b':
			my_population.set_bottleneck( Tb, Tl, f);
			break;
		case 'i':
			if(ni[0]>n)fprintf(stderr, "n[0] (number of sample in species 0) cannot exceed n (the total sample size), bye\n"), exit(1);
			
			Ti[0].time=TI;
			Ti[0].sp[0]=0;   // set the structure
			Ti[0].sp[1]=1;
			Ti[0].id_sp=2;
			ni[1]=n-ni[0];   // and the ni array
			fN[2]=1;          // the ancestral size (species '2' here) is always Ne*1

			my_population.set_isolation( Ti, fN, 2, M );
			
			break;
		case 'e':
			my_population.set_exponential_growth( Gr );
			break;
		case 'l':
			my_population.set_linear_growth( Gr );
			break;
		
	}
	
	/*
		Init the samples
	*/
	coalescent_tree_diversity *my_tree = new coalescent_tree_diversity( n, &my_population );

	

	/*
		Print out header
	*/

	printheader( &my_population, n , ni, rep, k,  opt_scenario );

	/*
		Run the process
		for each replicate,
			1. do a  tree
			2. compute Xi
			3. store it into the total Xi array
	*/
	
	
	long r;  // the replicates
	
	
	
	/*
		Run process
	*/
	
	for(r=0; r<rep; r++){
		
	
		switch(opt_scenario){
		
			case 's':
				my_tree->sweep_coalescent( );  
				break;
			case 'b':
				my_tree->bottleneck_coalescent();
				break;
			case 'i':
				my_tree->isolation_coalescent(ni);
				break;
			case 'n':
				my_tree->nested_coalescent( k );
				break;
			case 'e':
				if( Gr ){
					my_tree->growth_coalescent( 'e' ); 
					break;
				}
			case 'l':
				if( Gr ){
					my_tree->growth_coalescent( 'l' ); 
					break;
				}
			default:
				my_tree->standard_coalescent();
		
		}
		
		if(max_TMRCA){
		
			if( my_tree->get_Tmrca() > max_TMRCA ){
				
				my_tree->clean_whole_tree(  );
				r--;
				continue;
			}
		}
		

		if(opt_verbose>1){
			my_tree->printout_treeWithMut();                               // output all tree under the MRCA node
			my_tree->printout_sequences();                          // output all sequence under the MRCA node
		}

		
		
		for(int ncoal=1; ncoal < n; ncoal ++){
		
			int nlines = n - ncoal +1;  /* number of lineages that is in that step */
			
			float T_i = my_tree->get_Tcoal(ncoal)-my_tree->get_Tcoal(ncoal-1) ; 
			
			float T_med = (my_tree->get_Tcoal(ncoal)+my_tree->get_Tcoal(ncoal-1))/2;
			
		
			if(opt_verbose)
				printf("%d\t%f\t%f\t%f\n", nlines,(my_tree->get_Tcoal(ncoal)+my_tree->get_Tcoal(ncoal-1))/2, T_i, T_i*nlines*(nlines-1.0)/2.0 );
			
			
			int bin = (int) floor(T_med/binsize);
			if(bin>nbins-1)bin=nbins-1;
			
			my_bins[ bin ].n ++ ;
			my_bins[ bin ].sum += T_i*nlines*(nlines-1.0)/2.0;
			
		}
		
		
		

		my_tree->clean_whole_tree(  );                            // remove all mutations and species tag

		
	}


	for(int i=0; i < nbins; i ++)
		printf("%f %f\n", i*binsize, (my_bins[ i ].n>0)?my_bins[ i ].sum / my_bins[ i ].n :0);


	
	/*
		Empty memory
	*/


	delete my_tree;   

	return 0;
}

#undef KRONECKER
