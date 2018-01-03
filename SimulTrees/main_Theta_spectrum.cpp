/*
	file:     main_Theta_spectrum.cpp
	function: Perform simulations and output the Theta spectrum (Xi * i)
	author:   <gachaz>
	date:     sept-oct 2008
*/


#include <iostream>
using namespace std;

#include <stdlib.h>
#include <math.h>

#include "population.h"
#include "sample.h"
#include "random.h"
#include "coalescent_tree_diversity.h"
#include "distribution.h"


#define _VERSION_ "Oct 08"

#define KRONECKER( A, B ) (( (A) == (B) )?1:0)

static void printusage( char *argv0 ){

		cerr << "Syntax is "<<argv0<<" [-h] [ options ] [ scenario ] n Theta rep\n";
}

static void printhelp( char *argv0 ){

		printusage( argv0 );

		cerr << "\n[core values]\n";
		cerr << "\tn                   : total sample size.\n";
		cerr << "\tTheta               : Population mutation rate (2pNmu).\n";
		cerr << "\trep                 : number of replicates (independant trees).\n";


		cerr << "\n[info]\n";
		cerr << "\t-v                  : be verbose\n";
		cerr << "\t-V                  : be extra verbose\n\n";

		cerr << "\n[option]\n";
		cerr << "\t-p 1|2              : ploidy - either 1 -haploid- or 2 -diploid- (default 1)\n";
		cerr << "\t-C #                : Confidence Interval (default 0.95)\n";
		cerr << "\t-R #.#              : optionnal 'R'ate of sequencing/cloning errors (mutations added on the leaves)\n";
		cerr << "\t                      The total number of errors is drawn from Poisson with mean R*n (default 0)\n";
		cerr << "\t-x #                : set the seed for random numbers to # (otherwise time() is used)\n";
		cerr << "\t-S #                : set the number of mutations in the tree (instead of a Poisson random number using tree length and theta)\n";
		cerr << "\t-F                  : report a folded frequency spectrum\n";
		
		cerr << "\n[scenario]\n!! At most, a single scenario is chosen !!\n";
		
		cerr << "\t-s \"Ts p R alpha\"       : perform a sweep using the parameters Ts time bwd to sweep, p freq at the sweep end\n"
		        "\t                          R recomb rate (p*2Nr) and alpha selection strength (p*Ns).\n"
		        "\t                          when p<1 and Ts=0, the sweep is ongoing. When Ts>0, p is set to 1\n";
		
		cerr << "\t-b \"Tb Tl f\"            : perform a bottleneck using the parameters Tb time bwd to bottleneck,\n"
		        "\t                          Tl time of the event and f bottleneck relative size.\n";

		cerr << "\t-e Gr                   : perform an exponential growth using the parameter Gr for growth rate.\n";

		cerr << "\t-i \"Ti n0 f[0] f[1] M\"  : perform a isolation with migration with 2 pop using the parameters Ti time bwd to isolation/speciation,\n"
		        "\t                          n0 the number of sample in species 0 and f[i] size of species i relative to the ancestral one\n"
			"\t                          and M the 2Nm migration rate (symetric for now).\n";

		cerr << "\n";
}



static void printheader( population *my_population, int n , int *ni, int rep, float error_rate, char chosen_scenario ){


	cout << "/* Perform simulations to get Freq and theta spectrum ("<<_VERSION_<<")\n";
	cout << "/*    Population:\n";
	cout << "/*      * theta = "<<my_population->get_Theta()<<"\t("<<2*my_population->get_ploidy()<<".N.mu ; here E[T2]=Theta/2)\n";
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
			        "                   (population expands "<< exp(my_population->get_exponential_growth_Gr())  <<" every "<< (int)my_population->get_ploidy()<<"Ne gen.)\n";
	
		default:
			;
	}
	
	if( error_rate != 0 ){
		cout << "/*    Sequencing Errors:\n";
		cout << "/*      * rate  = "<< error_rate <<" (per sequence)\n";
		cout << "/*      * E[r]  = "<< error_rate * n <<" in the sample (Poisson dev)\n";
	}

	extern Random R;
	
	cout << "/*    seed random generator:\n";
	cout << "/*      * seed  = "<<R.get_seed()<<"\n";

	cout << "\n";
}



int main(int argc, char *argv[]){
	
	/*
		Parse options
	*/
	
	int optc;
	extern char *optarg;
	extern int optind;

	double error_rate=0.0;   // the error_rate per sequence 
	char ploidy=1;          // 1 for haploid, 2 for diploid
	
	char opt_verbose=0;     // print summary statistics of each rep
	
	char opt_scenario=0;    // it becomes 's' for sweep, 'S' for speciation, 'b' for bottleneck and 'e' for exponential growth


	double TS, Rec, alpha, p;    // for a sweep
	double Tb, Tl, f;       // for a bottleneck
	double Gr=0;            // for an exp growth

	struct isolation_event Ti[1];  // in struct, there are the time and which species are merged. Here only species 0 and 1 merge
	double TI;
	double fN[3];                  // relative size of each species. 0->(ns-1) today's species. >=ns: anceestral ones and 2ns-2 the FIRST species
	double M;                      // migration rate 2Nm
	int ni[3];                    // number of loci in today's species. Sum is n.
	

	int opt_Sites=-1;            // fixed number of sites instead of using theta
	int opt_Folded=0;            // by default report the full FS. if -F is set, report the folded FS
	
	
	/*
		Parse input
	*/
	

	extern char debug;
	extern Random R;


	debug=0;

	while( (optc=getopt(argc, argv, "hR:p:vVs:b:i:e:x:uS:F")) != -1 ){
	
		switch(optc){
		
			case 'v':
				opt_verbose = 1;
				break;
			
			case 'V':
				opt_verbose = 2;
				break;
			
			case 'R':
				error_rate = atof( optarg );
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

			case 'x':
				R.set_seed( atol(optarg) );
				break;

			case 'h':
				printhelp( argv[0] );
				exit(1);
				break;
			
			case 'S':
				opt_Sites = atoi( optarg );
				if(opt_Sites < 0)fprintf(stderr, "Number of mutations should be a positive value, bye\n"), exit(1);
				break;
			
			case 'F':
				opt_Folded = 1;
				break;
			
		}	
	}
	
 	if(argc-optind != 3)
		printusage( argv[0] ),exit(1);	
	 	
	int n         = atoi(argv[optind]);
	double Theta  = atof(argv[optind+1]);
	int rep       = atoi(argv[optind+2]);
		


	/*
		Set the population parameters
	*/
	
	population my_population;                                    // with Ne=1, mu=1 and haploid;
	my_population.set_Theta( Theta );
	my_population.set_ploidy( ploidy );
	 	
	
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
	
	}
	
	/*
		Init the samples
	*/
	coalescent_tree_diversity *my_tree = new coalescent_tree_diversity( n, &my_population );

	

	/*
		Print out header
	*/

	printheader( &my_population, n , ni, rep, error_rate, opt_scenario );

	/*
		Run the process
		for each replicate,
			1. do a  tree
			2. compute Xi
			3. store it into the total Xi array
	*/
	
	long r;  // the replicate counter
	
	long *Total_Xi = new long[n-1];
	long *Total_Xi_sq = new long[n-1];
	
	for(int x=1;x<n;x++)
		Total_Xi[x-1]=Total_Xi_sq[x-1]=0;

	
	distribution *distrib_TMRCA = new distribution(rep, 0.95);
	distribution *distrib_PL = NULL;

	if(opt_scenario == 's')
		distrib_PL = new distribution(rep, 0.95);
	
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
			case 'e':
				if( Gr ){
					my_tree->exponential_growth_coalescent(); 
					break;
				}
			default:
				my_tree->standard_coalescent();
		
		}
		if(opt_Sites == -1)
			my_tree->mutate_whole_tree( error_rate*(double)n );             // add some mutations
		else
			my_tree->mutate_whole_tree( opt_Sites, error_rate*(double)n );

		if(opt_verbose>1){
			my_tree->printout_tree();                               // output all tree under the MRCA node
			my_tree->printout_sequences();                          // output all sequence under the MRCA node
		}

		my_tree->compute_Xi_array( );                                   // K, S and Xi[i]

		for(int x=0;x<n-1;x++){
			Total_Xi[x] += my_tree->get_Xi_array()[x];
			Total_Xi_sq[x] += (my_tree->get_Xi_array()[x])*(my_tree->get_Xi_array()[x]);
		}
		
		if(opt_verbose>0){
		
			cout << "replicate #"<<r<<"\n";
		
			my_tree->print_summary_statistics();
			if(opt_scenario == 's')
				cout << "Surviving lineages= "<< my_tree->get_nsurvivors() <<"\n";
		}

		distrib_TMRCA->add_value( my_tree->get_all_samples()[n-1]->get_time()  );
		
		if(opt_scenario == 's')
			distrib_PL->add_value( my_tree->get_nsurvivors()  );


		my_tree->clean_whole_tree(  );                            // remove all mutations

		
	}

	printf("TMRCA\t");
	printf("Mean  = %f\tStdev = %f\t", distrib_TMRCA->get_mean(), sqrt(distrib_TMRCA->get_variance()) );
	printf("95%%CI = %f %f\n", distrib_TMRCA->get_low_boundary() , distrib_TMRCA->get_high_boundary() );
		
	if(opt_scenario == 's'){
		printf("Surv_n\t");
		printf("Mean  = %f\tStdev = %f\t", distrib_PL->get_mean(), sqrt(distrib_PL->get_variance()) );
		printf("95%%CI = %f %f\n", distrib_PL->get_low_boundary() , distrib_PL->get_high_boundary() );
	}

	
	if( ! opt_Folded ){
	
		
		
		printf("/*\n\tXI and THETA_i spectrum\n\n");
		printf("i\tE[Xi]     \tE[Th_i]      \tStdev[Xi]     \tStdev[Th_i]\n");
		for(int x=0;x<n-1;x++)
		      printf("%d\t%f\t%f\t%f\t%f\n",
		        	 x+1,
		         	Total_Xi[x]/(double)r, 
			 	(x+1.0)*(double)Total_Xi[x]/(double)r,
			 	sqrt( (Total_Xi_sq[x]/(double)r - Total_Xi[x]*Total_Xi[x]/((double)r*r) )* r / (r-1.0) ),
			 	sqrt( (x+1.0)*(x+1.0)* (Total_Xi_sq[x]/(double)r - Total_Xi[x]*Total_Xi[x]/((double)r*r) ) * r /  (r-1.0) )
				 );
		
		
		/*
		printf("/*\n\tUnfolded Frequency Spectrum\n\n");
		printf("i\tE[Xi]    \tStdev[Xi]\n");
		
		for(int x=0;x<n-1;x++)
		      printf("%d\t%f\t%f\n",
		        	 x+1,
		         	Total_Xi[x]/(double)r, 
			 	sqrt( (Total_Xi_sq[x]/(double)r - Total_Xi[x]*Total_Xi[x]/((double)r*r) )* r / (r-1.0) )
				 );
		*/		 
	}
	else{
	
		
		printf("/*\n\tETA and THETA_i spectrum\n\n");
		printf("i\tE[Eta_i]     \tE[Th_i]      \tStdev[Eta_i]     \tStdev[Th_i]\n");
		for(int x=0;x<=(n-2)/2;x++){
		
			int t_eta, t_eta_sq;
			float w;
			
			t_eta = (Total_Xi[x] + Total_Xi[n-2-x])/( 1+KRONECKER(x,n-2-x) ) ;
			t_eta_sq = (Total_Xi_sq[x] + Total_Xi_sq[n-2-x])/( 1+KRONECKER(x,n-2-x) ) ;
			w = ((x+1)*(n-1-x)*(1+KRONECKER(x+1,n-1-x)) / (n+0.0) );
					
			printf("%d\t%f\t%f\t%f\t%f\n",
		        	 x+1,
		         	t_eta / (r+0.0) ,
			 	w *  ( t_eta / (r+0.0) ),
			 	sqrt( (t_eta_sq/(double)r - t_eta*t_eta/((double)r*r) )* r / (r-1.0) ),
			 	sqrt( w*w*(t_eta_sq/(double)r - t_eta*t_eta/((double)r*r) ) * r /  (r-1.0) )
				 );
		
		/*	
		printf("/*\n\tFolded Frequency Spectrum\n\n");
		printf("i\tE[Eta]    \tStdev[Eta]\n");
		for(int x=0;x<=(n-2)/2;x++){
		
			int t_eta, t_eta_sq;
			
			t_eta = (Total_Xi[x] + Total_Xi[n-2-x])/( 1+KRONECKER(x,n-2-x) ) ;
			t_eta_sq = (Total_Xi_sq[x] + Total_Xi_sq[n-2-x])/( 1+KRONECKER(x,n-2-x) ) ;
					
			printf("%d\t%f\t%f\n",
		        	 x+1,
		         	t_eta / (r+0.0) ,
			 	sqrt( (t_eta_sq/(double)r - t_eta*t_eta/((double)r*r) )* r / (r-1.0) )
				 );
		*/
		}
	}
	
			

	delete Total_Xi;
	delete Total_Xi_sq;
	delete distrib_TMRCA;

	delete my_tree;   

	return 0;
}

#undef KRONECKER
