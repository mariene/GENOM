/*
	file:     main_generate_CI_coal.cpp
	function: Using various scheme, generate a CI for several values
	author:   <amikezor>
	date:     june 12
	modif:    
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
#include "readfile.h"

#include "get_theta_CI.c"

#define _VERSION_ "Feb 2012"

char opt_verbose;

long factorial( long x ){

	long f=1;
	
	for( int i=1; i<=x; i++)
		f *= i;
		
	return f;
}

float log_factorial( long x ){
	
	float log_f=0;
	
	for( int i=1; i<=x; i++)
		log_f += log(i);
		
	return log_f;
}

static void printusage( char *argv0 ){

		cerr << "Syntax is "<<argv0<<" [-h] [ options ] [ scenario ] n Theta S rep\n";
}

static void printhelp( char *argv0 ){

		printusage( argv0 );

		cerr << "\n[core values]\n";
		cerr << "\tn             : total sample size.\n";
		cerr << "\tTheta         : Population mutation rate (2pNmu).\n";
		cerr << "\tS             : sets the conditionnal segregating sites.\n";
		cerr << "\trep           : number of replicates (independant trees).\n";


		cerr << "\n[info]\n";
		cerr << "\t-v            : be verbose\n";

		cerr << "\n[option]\n";
		cerr << "\t-p 1|2        : ploidy - either 1 -haploid- or 2 -diploid- (default 1)\n";
		cerr << "\t-r            : use rejection algorithm (default use weight)\n";
		cerr << "\t-w #          : 'u': use uniform prior weight; 'l' log-uniform ; 'g'amma ; 't'heta value --given or theta_S\n";
		cerr << "\t-m #          : set mean (for gamma)\n";
		cerr << "\t-s #          : set stdev (for gamma)\n";

		cerr << "\n[scenario]\n!! At most, a single scenario is chosen !!\n";
		
		
		cerr << "\t-b \"Tb Tl f\"           : perform a bottleneck using the parameters Tb time bwd to bottleneck,\n"
		        "\t                          Tl time of the event and f bottleneck relative size.\n";
		cerr << "\t-S \"Ts p R alpha\"      : perform a sweep using the parameters Ts time bwd to sweep, p the frequency of the selected allele\n"
		        "\t                          R recomb rate (2pNr) and alpha selection strength (pNs)\n"
		        "\t                          if Ts>0, p is set to p=1 (sweep is over); when Ts=0, 0<=p<=1 (sweep is ongoing)\n";
		cerr << "\t-i \"Ti n0 f[0] f[1] M\" : perform a isolation with migration with 2 pop using the parameters Ti time bwd to isolation/speciation,\n"
		        "\t                          n0 the number of sample in species 0 and f[i] size of species i relative to the ancestral one\n"
			"\t                          and M the 2Nm migration rate (symetric for now).\n";
		cerr << "\t-e Gr                  : perform an exponential growth using the parameter Gr for growth rate.\n";

		cerr << "\n";
}



static void printheader( population *my_population, int n ,int rep, int S, char opt_rejection, char opt_weight ){


	cout << "/* Perform simulations to the P(\\pi | S) ("<<_VERSION_<<")\n";
	cout << "/*    Population:\n";
	cout << "/*      * theta = "<<my_population->get_Theta()<<"\t("<<2*my_population->get_ploidy()<<".N.mu ; here E[T2]=Theta/2)\n";
	cout << "/*      * ploidy = "<< 1*my_population->get_ploidy() <<" \n";
	cout << "/*    Sample:\n";
	cout << "/*      * n     = "<<n<<" loci\n";
	cout << "/*      * S     = "<<S<<" segregating sites\n";
	cout << "/*    Iterations:\n";
	cout << "/*      * rep   = "<<rep<<" iterations\n";
	
	cout << "/*    Algorithm:\n";
	cout << "/*      * type   = " << (opt_rejection?"rejection":(opt_weight?"weight":"theta")) << "\n";
	cout << "/*      * weight = "<<opt_weight<<"\n";
	
	extern Random R;
	
	cout << "/*    seed random generator:\n";
	cout << "/*      * seed  = "<<R.get_seed()<<"\n";

	cout << "\n";
}



static double Poisson( double rate, int x ){

	if(x<10)
		return exp( -rate )* pow( rate , (double) x ) / (double)factorial( (long)x );
	else
		return exp( -rate + x*log(rate) - log_factorial(x)  );


}

int main(int argc, char *argv[]){
	
	/*
		Parse options
	*/
	
	int optc;
	extern char *optarg;
	extern int optind;

	char ploidy=1;          // 1 for haploid, 2 for diploid
	
	char opt_verbose=0;     // print summary statistics of each rep
	char opt_rejection=0;   // if 1 use rejection instead of weight
	char opt_weight=0;      // if non-zero, specified a weighting scheme

	float CI=0.95;
	
	float g_mean=-1,
	      g_var=-1;     // when weight is gamma, it is mlean and variance

	int i;
	float an;


	/*
		Scenarios
	*/
	char opt_scenario=0;    // it becomes 's' for sweep, 'S' for speciation, 'b' for bottleneck and 'e' for exponential growth


	double TS, Rec, alpha, p;    // for a sweep
	double Tb, Tl, f;            // for a bottleneck
	double Gr=0;                 // for an exp growth

	/*
		For isolation event new version
	*/

	struct isolation_event Ti[1];  // in struct, there are the time and which species are merged. Here only species 0 and 1 merge
	double TI;
	double fN[3];                  // relative size of each species. 0->(ns-1) today's species. >=ns: anceestral ones and 2ns-2 the FIRST species
	double M;                      // migration rate 2Nm
	int ni[3];                    // number of loci in today's species. Sum is n.



	int nlocus=1;

	/*
		Parse input
	*/
	

	extern char debug;
	extern Random R;

	debug=0;

	while( (optc=getopt(argc, argv, "rhp:vVx:w:C:m:s:b:L:")) != -1 ){
	
		switch(optc){
		
			case 'r':
				opt_rejection=1;
				break;
		
			case 'v':
				opt_verbose = 1;
				break;
			
			case 'V':
				opt_verbose = 2;
				break;
			
			
			case 'p':
				ploidy = (char)atoi( optarg );
				if(ploidy != 1 && ploidy != 2)fprintf(stderr, "ploidy can only be set to 1 or 2\n"), exit(1);
				break;
				
			case 'x':
				R.set_seed( atol(optarg) );
				break;
			
			case 'h':
				printhelp( argv[0] );
				exit(1);
				
			case 'w':
				opt_weight = optarg[0];
				break;

			case 'C':
				CI= atof(optarg);
				if(CI<0 || CI>1)fprintf(stderr, "CI is in [0,1]\n"), exit(1);
				break;
				
			case 'm':
				g_mean = atof(optarg);
				break;

			case 's':
				g_var = pow(atof(optarg),2);
				break;

			case 'S':
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

			case 'L':
				nlocus = atoi(optarg);
				break;


		}	
	}
	
 	if(argc-optind != 4)
		printusage( argv[0] ),exit(1);	

	int n         = atoi(argv[optind]);
	double Theta  = atof(argv[optind+1]);
	int S         = atoi(argv[optind+2]);
	int rep       = atoi(argv[optind+3]);

	 
	an=1.0;
	for(i=2;i<n;i++)
		an += 1.0/i;
	 
	if(opt_weight == 'g' && (g_mean == -1 || g_var == -1) ){
		cerr << "when using a gamma, please provide a postive mean and variance (using options -m and -s)\n";
		exit(1);
	}
		 
	if(opt_weight != 0 && S <= 0){
		cerr << "when using weight scheme, please provide S>=1\n";
		exit(1);
	}
		 
		


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
		case 'e':
			my_population.set_exponential_growth( Gr );
			break;
	
		case 'i':
			if(ni[0]>n)fprintf(stderr, "n[0] (number of sample in species 0) cannot exceed n (the total sample size), bye\n"), exit(1);
			
			if(TI<=0)fprintf(stderr, "TI has to be >0, bye\n"), exit(1);
						
			Ti[0].time=TI;
			Ti[0].sp[0]=0;   // set the structure
			Ti[0].sp[1]=1;
			Ti[0].id_sp=2;
			
			ni[1]=n-ni[0];    // and the ni array (at the begining of the process)
			
			fN[2]=1;          // the ancestral size (species '2' here) is always Ne*1


			my_population.set_isolation( Ti, fN, 2, M );
			
			break;
	}
	
	/*
		Init the samples
	*/
	coalescent_tree_diversity *my_tree = new coalescent_tree_diversity( n, &my_population );


	       
	/*
		Print out header
	*/

	printheader( &my_population, n , rep, S, opt_rejection, opt_weight );

	/*
		Set tools for distributions
	*/	
	long r;
	bool full=false;

	if(opt_weight != 0)
		full=true;
		
	distribution * pi= new distribution(rep, CI, full, full);
	distribution * eta1= new distribution(rep, CI, full, full);
	distribution * D= new distribution(rep, CI, full, full);
	distribution * Du= new distribution(rep, CI, full, full);
	distribution * D2star= new distribution(rep, CI, full, full);
	distribution * Haplos= new distribution(rep, CI, full, full);
	distribution * TMRCA= new distribution(rep, CI, full, full);


	
	if(S>=0){
	
		double theta_min, theta_max;
	
		compute_theta_CI( S, n, &theta_min, &theta_max, 1-CI, 0.01, 0.0001);
		printf("theta in [%.3f, %.3f] (using P[S=%d | theta]) and CI=%f\n",  theta_min, theta_max, S,CI);
	
	}
	
	/*
		Run the process r times
		for each replicate,
			1. do a  tree
			2. compute pi if S is adequate
	*/
	for(r=0; r<rep; r++){
		
		int s;
		float weight=0;
		double eta1_c=0,
		        pi_c=0,
			 D_c=0,
			 Du_c=0,
			 TMRCA_c=0,
			 D2star_c=0,
			 Haplos_c=0,
			 weight_c=0;

		int X;
		
		for( X=0; X<nlocus; X++){
		
			/*
				Run the appropriate coalescent process
			*/
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
					my_tree->growth_coalescent( 'e' ); 
					break;
				default:
					my_tree->standard_coalescent();
		
			}
		

			if( opt_weight != 0 ){
		
				s=S;
		
				switch( opt_weight ){
					
					case 'u':
						weight = 1.0/my_tree->get_Ttot();
						break;
					case 'l':
						weight = 1.0;
						break;
					case 'g':
						weight = pow( my_tree->get_Ttot()/2, S)/pow( g_mean/g_var+my_tree->get_Ttot()/2, g_mean*g_mean/g_var+S);
						break;
					case 't':
						if(Theta <=0)
							weight = Poisson( (S/an)*my_tree->get_Ttot()/2.0, S);
						else
							weight = Poisson( Theta*my_tree->get_Ttot()/2.0, S);
						break;
				
				}
		
			}else{
		
				if(Theta <= 0)
					my_population.set_Theta( S/an ); /* set it to Theta_S */
			
				s = R.poisson_dev( my_tree->get_Ttot() * my_population.get_Theta() / 2.00  );

				if( opt_rejection && s != S){
					r--;
					continue;
				}
			
				weight = -1;
		
			}
		
			//printf("mutate tree with %d mutations\n", s);
		
			my_tree->mutate_whole_tree( s, 0 );
			my_tree->compute_Xi_array();
		
		
			eta1_c += (double) my_tree->get_Eta1();
			pi_c += my_tree->get_K();
			TMRCA_c += my_tree->get_all_samples()[n-1]->get_time();
	
			my_tree->compute_D(1, 1);
			D_c += my_tree->get_T();
			
			my_tree->compute_D(0, 1);
			Du_c += my_tree->get_T();
		
			my_tree->compute_Dstar_fl(1, 1);
			D2star_c += my_tree->get_T();
		
			my_tree->compute_Haplos();
			Haplos_c = (double)my_tree->get_Haplos();
	
		
			weight_c += weight;
		
			/*
				Verbose options
			*/
			if(opt_verbose>1){
				my_tree->printout_tree();                               // output all tree under the MRCA node
				//my_tree->printout_sequences();                        // output all sequence under the MRCA node
			}
			if(opt_verbose>0){
				my_tree->print_summary_statistics();
			}

		

			/*
				Update adequate distributions
			*/

			my_tree->clean_whole_tree(  );                            // remove all mutations
		}
		
		eta1->add_value( eta1_c/nlocus, weight_c/nlocus  );
		pi->add_value( pi_c/nlocus, weight_c/nlocus  );
		TMRCA->add_value( TMRCA_c/nlocus, weight_c/nlocus  );
		D->add_value( D_c/nlocus, weight_c/nlocus );
		Du->add_value( Du_c/nlocus, weight_c/nlocus );
		D2star->add_value( D2star_c/nlocus, weight_c/nlocus );
		Haplos->add_value( Haplos_c/nlocus, weight_c/nlocus  );
		

	}

	
	//TMRCA->print_distribution();
	
	/*
		Output
	*/
	printf("TMRCA\t");
	printf("Mean  = %f\tVar = %f\t", TMRCA->get_mean(), TMRCA->get_variance() );
	printf("%.1f%%CI = %f %f\n", CI*100, TMRCA->get_low_boundary() , TMRCA->get_high_boundary() );
		
	printf("Pi\t");
	printf("Mean  = %f\tVar = %f\t", pi->get_mean(), pi->get_variance() );
	printf("%.1f%%CI = %f %f\n", CI*100, pi->get_low_boundary() , pi->get_high_boundary() );
		
	printf("Eta1\t");
	printf("Mean  = %f\tVar = %f\t", eta1->get_mean(), eta1->get_variance() );
	printf("%.1f%%CI = %f %f\n", CI*100, eta1->get_low_boundary() , eta1->get_high_boundary() );
		
	printf("D\t");
	printf("Mean  = %f\tVar = %f\t", D->get_mean(), D->get_variance() );
	printf("%.1f%%CI = %f %f\n", CI*100, D->get_low_boundary() , D->get_high_boundary() );
		
	printf("Du\t");
	printf("Mean  = %f\tVar = %f\t", Du->get_mean(), Du->get_variance() );
	printf("%.1f%%CI = %f %f\n", CI*100, Du->get_low_boundary() , Du->get_high_boundary() );
		
	printf("D2*\t");
	printf("Mean  = %f\tVar = %f\t", D2star->get_mean(), D2star->get_variance() );
	printf("%.1f%%CI = %f %f\n", CI*100, D2star->get_low_boundary() , D2star->get_high_boundary() );
		
	printf("Haplos\t");
	printf("Mean  = %f\tVar = %f\t", Haplos->get_mean(), Haplos->get_variance() );
	printf("%.1f%%CI = %f %f\n", CI*100, Haplos->get_low_boundary() , Haplos->get_high_boundary() );
		

	/*
		Free Memory
	*/

	delete pi;
	delete eta1;
	delete D;
	delete D2star;
	delete Haplos;
	delete TMRCA;
	delete my_tree;   

	return 0;
}
