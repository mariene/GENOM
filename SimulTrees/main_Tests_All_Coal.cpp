/*
	file:     main_D_All_Coal.cpp
	function: Perform Similation and D, Ystar, Y tests under various coalescent scenarii
	author:   <amikezor>
	date:     march-may 07
	modif:    nov 08 - move to the generic tests
	                   feed the software, with a generic weight vector and its associated CI (for S given)
	modif:    feb 09 - add sampling in an ongoing sweep + debug when a CI file is given in argument + improve header
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

#define _VERSION_ "Feb 09"

char opt_verbose;


static void printusage( char *argv0 ){

		cerr << "Syntax is "<<argv0<<" [-h] [ options ] [ scenario ] n Theta mode rep\n";
}

static void printhelp( char *argv0 ){

		printusage( argv0 );

		cerr << "\n[core values]\n";
		cerr << "\tn                   : total sample size.\n";
		cerr << "\tTheta               : Population mutation rate (2pNmu).\n";
		cerr << "\tmode                : sets which test is used (1: Omega ; 2: D ; 3,4: Y,Y* ;\n";
		cerr << "\t                       5,6: Dfl,D*fl ; 7,8: F,F* ; 9: H; 10: Haplos; 11: Pi).\n";
		cerr << "\t                      mode can be several tests (e.g. '2:10' means D and Haploid test)\n";
		cerr << "\trep                 : number of replicates (independant trees).\n";


		cerr << "\n[info]\n";
		cerr << "\t-v                  : be verbose\n";
		cerr << "\t-V                  : be extra verbose\n\n";

		cerr << "\n[option]\n";
		cerr << "\t-p 1|2              : ploidy - either 1 -haploid- or 2 -diploid- (default 1)\n";
		cerr << "\t-f fileCI           : read CI file\n";
		cerr << "\t-w file_w1          : use w1 as a weight vector -- it should be the right size -- , w2 being uniform, it tuens into mode 1\n";
		cerr << "\t-C #                : Confidence Interval (default 0.95) -- ignored when a file is given (-f option)\n";
		cerr << "\t-R #.#              : optionnal 'R'ate of sequencing/cloning errors (mutations added on the leaves)\n";
		cerr << "\t                      The total number of errors is drawn from Poisson with mean R*n (default 0)\n";
		cerr << "\t-x #                : set the seed for random numbers to # (otherwise time() is used)\n";
		cerr << "\t-r #                : set # replicates when no CI file is given (default 2000)\n";
		cerr << "\t-u                  : compute unormalized tests (mode 1 to 9)\n";
		
		cerr << "\n[scenario]\n!! At most, a single scenario is chosen !!\n";
		
		cerr << "\t-s \"Ts p R alpha\"   : perform a sweep using the parameters Ts time bwd to sweep, p the frequency of the selected allele\n"
		        "\t                      R recomb rate (2pNr) and alpha selection strength (pNs)\n"
		        "\t                      if Ts>0, p is set to p=1 (sweep is over); when Ts=0, 0<=p<=1 (sweep is ongoing)\n";
		
		cerr << "\t-b \"Tb Tl f\"        : perform a bottleneck using the parameters Tb time bwd to bottleneck,\n"
		        "\t                      Tl time of the event and f bottleneck relative size.\n";

		
		cerr << "\t-i \"Ti n0 f[0] f[1] M\"  : perform a isolation with migration with 2 pop using the parameters Ti time bwd to isolation/speciation,\n"
		        "\t                            n0 the number of sample in species 0 and f[i] size of species i relative to the ancestral one\n"
			"\t                            and M the 2Nm migration rate (symetric for now).\n";

		cerr << "\t-e Gr               : perform an exponential growth using the parameter Gr for growth rate.\n";

		cerr << "\n";
}



static void printheader( population *my_population, int n , int * ni, int rep, float CI, float error_rate, char chosen_scenario, int mode, char inputCIfile ){


	cout << "/* Perform simulations to assess the CI and the power of neutrality test in various scenarios (Omega, D,...) ("<<_VERSION_<<")\n";
	cout << "/*    Population:\n";
	cout << "/*      * theta = "<<my_population->get_Theta()<<"\t("<<2*my_population->get_ploidy()<<".N.mu ; here E[T2]=Theta/2)\n";
	cout << "/*      * ploidy = "<< 1*my_population->get_ploidy() <<" \n";
	cout << "/*    Sample:\n";
	cout << "/*      * n     = "<<n<<" loci\n";
	cout << "/*    Iterations:\n";
	cout << "/*      * rep   = "<<rep<<" iterations\n";
	cout << "/*    Confidence Interval:\n";
	cout << "/*      * CI    = "<<CI<<"\n";
	if( inputCIfile == 0)
		cout << "/*      * from  = on the fly\n";
	else
		cout << "/*      * from  = an input file\n";

	switch(mode){
		case 1: 
			cout << "/*      * test = Omega ("<<mode<<") \n";
			break;
		case 2:
			cout << "/*      * test = D ("<<mode<<") \n";
			break;
		case 3:
			cout << "/*      * test = Y ("<<mode<<") \n";
			break;
		case 4: 
			cout << "/*      * test = Y* ("<<mode<<") \n";
			break;
		case 5:
			cout << "/*      * test = Dfl ("<<mode<<") \n";
			break;
		case 6:
			cout << "/*      * test = D*fl ("<<mode<<") \n";
			break;
		case 7:
			cout << "/*      * test = F ("<<mode<<") \n";
			break;
		case 8:
			cout << "/*      * test = F* ("<<mode<<") \n";
			break;
		case 9:
			cout << "/*      * test = H ("<<mode<<") \n";
			break;
		case 10:
			cout << "/*      * test = Haplos ("<<mode<<") \n";
			break;
			
		case 11:
			cout << "/*      * test = Pi ("<<mode<<") \n";
			break;
			
	}
	
	
	switch(chosen_scenario){
	
		case 's':
			cout << "/*    [sweep]\n";
			cout << "/*      * Ts     = "<<my_population->get_sweep_TS()<<" (in "<<my_population->get_ploidy()<<"N generations)\n";
			cout << "/*      * alpha  = "<<my_population->get_sweep_alpha()<<" (select coef = "<<1*my_population->get_ploidy()<<"Ns)\n";
			cout << "/*      * p      = "<<my_population->get_sweep_p()<<" (frequency of the selected allele at time 0)\n";
			cout << "/*      * R      = "<<my_population->get_sweep_R()<<" (recomb coef = "<<2*my_population->get_ploidy()<<"Nr)\n";
			break;
		
		case 'i':
			cout << "/*    [isolation]\n";
			cout << "/*      * Ti     = "<<my_population->get_isolation_Ti(0).time<<"\t(in "<<  1*my_population->get_ploidy() << "N generations)\n";
			cout << "/*      * fN[0]   = "<<my_population->get_isolation_fN(0)<<"\t(real size is "<< 1*my_population->get_ploidy() <<".f[0].N)\n";
			cout << "/*      * fN[1]   = "<<my_population->get_isolation_fN(1)<<"\t(real size is "<< 1*my_population->get_ploidy() <<".f[1].N)\n";
			cout << "/*      * M      = "<<my_population->get_isolation_M()<<" (in "<< 2*my_population->get_ploidy()<<"Nm units)\n";
			cout << "/*    Sample:\n";
			cout << "/*      * ni[0]  = "<<ni[0]<<" loci for species 0\n";
			cout << "/*      * ni[1]  = "<<ni[1]<<" loci for species 1\n";
			break;
		
		case 'b':
			cout << "/*     [bottleneck]\n";
			cout << "/*      * Tb     = "<<my_population->get_bottleneck_Tb()<<"\t(in "<< 1*my_population->get_ploidy()<<"N gen)\n";
			cout << "/*      * Tl/f   = "<<my_population->get_bottleneck_Tl()/my_population->get_bottleneck_f()<<"\t(strength of the bottleneck in coalescent time)\n";
			cout << "/*        * Tl   = "<<my_population->get_bottleneck_Tl()<<"\t(in "<< 1*my_population->get_ploidy()<<"N gen)\n";
			cout << "/*        * f    = "<<my_population->get_bottleneck_f()<<"\t(population at bottleneck is f."<<1*my_population->get_ploidy()<<"N)\n";
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
	
	float CI=0.950;          // the confidence interval
	double epsilon=1e-4;     // if the difference between T_CI and obs is more than epsilon, count it as over.
	
	
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

	char norm=1;                  // becomes 0 for unormalized tests
	
	readfile inputfiles;           // where the Ci file is stored

	
	/*
		Parse input
	*/
	

	extern char debug;
	extern Random R;

	int rep2 = 2000;


	debug=0;

	while( (optc=getopt(argc, argv, "hR:p:C:vVs:b:i:e:x:uf:w:r:")) != -1 ){
	
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
				
			case 'C':
				CI = atof(optarg);
				if(CI<0 || CI>1)fprintf(stderr, "Confidence Interval has to be 0<=CI<=1, bye"), exit(1);
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

			case 'u':
				norm=0;
				break;
				
			case 'f':
				inputfiles.readfile_CI_S(optarg);
				break;
				
			case 'w':
				inputfiles.readfile_w(optarg, 1);
				break;
			
			case 'r':
				rep2=atoi(optarg);
				break;
				
			case 'h':
				printhelp( argv[0] );
				exit(1);

		}	
	}
	
 	if(argc-optind != 4)
		printusage( argv[0] ),exit(1);	
	 	
	int n         = atoi(argv[optind]);
	double Theta  = atof(argv[optind+1]);
	int mode      = atoi(argv[optind+2]);
	int rep       = atoi(argv[optind+3]);
		

	if( inputfiles.get_size_w1() && mode != 1){
		cerr << "using input w1 vector, swicth to mode 1\n";
		mode = 1;
	}
	if( inputfiles.get_size_w1()  && inputfiles.get_size_w1() != n-1 ){
		cerr << "using input w1 vector, set n to "<< inputfiles.get_size_w1()+1 << "\n";
		n = inputfiles.get_size_w1()+1;
	}

	/*
		If mode 1 with no input w1 vector, use the freq spectrum
	*/
	if( inputfiles.get_size_w1() == 0 && mode == 1 ){
		cerr << "when using mode 1, please give an input w1 vector, bye\n";
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

	printheader( &my_population, n , ni, rep, CI, error_rate, opt_scenario,  mode, (inputfiles.get_size_CI_S())?1:0 );

	/*
		Set tools for distributions
	*/	
	long r;
	long power  = 0;

	distribution *distrib_Val = new distribution(rep, CI);
	
	distribution *distrib_TMRCA = new distribution(rep, CI);
	distribution *distrib_PL = NULL;

	if(opt_scenario == 's')
		distrib_PL = new distribution(rep, CI);
	
	
	

	/*
		The set tools for CI
	*/

	double Val_up=0,
	       Val_low=0;
	       
	
	coalescent_tree_diversity *std_tree = new coalescent_tree_diversity( n, &my_population );
	distribution *distrib_std = new distribution(rep2, CI);
	
	if(inputfiles.get_size_CI_S() == 0){
		cerr << "!! The CI is computed for each run --very slow-- !!\n";
		cerr << "!! Please consider feeding with a file of pre-computed CIs !!\n";
	}


	/*
		The Freq Spectrum weight -- mode 1
	*/

	double *w1=NULL;
	double *w2=NULL;
	if(mode == 1){
	
		w2 = new double [ n-1 ];

		for(int x=0;x<n-1;x++)      // whatever, the second vector is uniform
			w2[x]=1;


		if( inputfiles.get_size_w1() == 0){
			w1 = new double[ n-1 ];
			//cerr << "use mode 1 with empirik w1\n";
		}
		else{
			w1 = inputfiles.get_w1();
			//cerr << "use mode 1 with input w1\n";
			//for(int x=0;x<n-1;x++){
			//	printf("%d .. %f %f\n", x, w1[x], w2[x]);
			//}

		}
		
	}
	
	
	/*
		Run the process r times
		for each replicate,
			1. do a  tree
			2. compute T statistics ( D, Omega or any othe Stat value)
			3. compare it to its CI
	*/
	
	for(r=0; r<rep; r++){
		
	
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
				if( Gr ){
					my_tree->growth_coalescent( 'e' ); 
					break;
				}
			default:
				my_tree->standard_coalescent();
		
		}
		
		
		/*
			Drop mutations in the tree
		*/
		my_tree->mutate_whole_tree( 0, error_rate*(double)n );          // draw S with Poisson(Theta/2*L) and add extra singletons Poisson(n*mean_err)
		my_tree->compute_Xi_array( );                                   // K, S and Xi[i], the full freq spec


		/*
			Verbose options
		*/
		if(opt_verbose>1){
			my_tree->printout_tree();                               // output all tree under the MRCA node
			//my_tree->printout_sequences();                        // output all sequence under the MRCA node
		}
		if(opt_verbose>1){
			my_tree->print_summary_statistics();
			if(opt_scenario == 's')
				cout << "Surviving lineages= "<< my_tree->get_nsurvivors() <<"\n";
		}



		
		
		/*
			Compute	the appropriate statistics
		*/
		switch(mode){
			case 1:
				my_tree->compute_T( w1, w2, norm, 0);
				break;
			case 2:
				my_tree->compute_D( norm, 0 );
				break;
			case 3:
				my_tree->compute_Y( norm, 0 );
				break;
			case 4:
				my_tree->compute_Ystar(norm, 0 );
				break;
			case 5:
				my_tree->compute_D_fl(norm, 0);
				break;
			case 6:
				my_tree->compute_Dstar_fl(norm, 0);
				break;
			case 7:
				my_tree->compute_F(norm, 0);
				break;
			case 8:
				my_tree->compute_Fstar(norm, 0);
				break;
			case 9:
				my_tree->compute_H(norm, 0);
				break;
			case 10:
				my_tree->compute_Haplos();
				break;
		}
		

		/*
			Set the CI
			If no CI file is given, compute the CI using rep2 replicates
		*/
		if( inputfiles.get_size_CI_S() == 0  ){
						
//			cerr << "Round "<<r<<"\n"; 
//			printf("Nsites : %d\n", my_tree->get_S() );
			
						
			for(int r2=0; r2<rep2; r2++){
			
				std_tree->standard_coalescent();
				std_tree->mutate_whole_tree( my_tree->get_S(), 0 );
				std_tree->compute_Xi_array( );
								
				switch(mode){
					case 1:
						std_tree->compute_T( w1, w2, norm, 0);
						break;
					case 2:
						std_tree->compute_D( norm, 0 );
						break;
					case 3:
						std_tree->compute_Y( norm, 0 );
						break;
					case 4:
						std_tree->compute_Ystar(norm, 0 );
						break;
					case 5:
						std_tree->compute_D_fl(norm, 0);
						break;
					case 6:
						std_tree->compute_Dstar_fl(norm, 0);
						break;
					case 7:
						std_tree->compute_F(norm, 0);
						break;
					case 8:
						std_tree->compute_Fstar(norm, 0);
						break;
					case 9:
						std_tree->compute_H(norm, 0);
						break;
					case 10:
						std_tree->compute_Haplos();
						break;
				}
								
				if( mode <= 9 )
					distrib_std->add_value( std_tree->get_T() );
				else if( mode == 10 )
					distrib_std->add_value( std_tree->get_Haplos() );
				else if( mode == 11 )
					distrib_std->add_value( std_tree->get_K() );

				
				std_tree->clean_whole_tree(  );
				
			}
						
			Val_low = distrib_std->get_low_boundary();
			Val_up  = distrib_std->get_high_boundary();
			
			//printf("CI: S[%d] sin %f %f\n", my_tree->get_S(), Val_low, Val_up );
			
			distrib_std->clean_distribution();
			

		}else{
		
			if( my_tree->get_S() > inputfiles.get_size_CI_S() )
				fprintf(stderr, "warrning, cannot read a CI for S= %d; use CI of S=%d instead\n", my_tree->get_S(), inputfiles.get_size_CI_S());
				
			int S=( my_tree->get_S() > inputfiles.get_size_CI_S() )?inputfiles.get_size_CI_S():my_tree->get_S();
		
			Val_low = inputfiles.get_CI_S( S )[0];
			Val_up = inputfiles.get_CI_S( S )[1];
		
		}
	
		
		/*
			Update adequate distributions
		*/
		double myT=0;
		
		if(mode <= 9){
			myT=my_tree->get_T();
		}
		else if(mode == 10){
			myT = my_tree->get_Haplos();
		}
		else if(mode == 11)
			myT = my_tree->get_K();

		distrib_Val->add_value( myT );
		
		/*
			Compare the Statistics to its CI
		*/
		if( (myT - Val_up) > epsilon ||   (Val_low - myT) > epsilon  )
			power++;
		
		if(opt_verbose)
			printf("S: %d T: %f\n", my_tree->get_S(),my_tree->get_T() );
		


//		cout << "Stat is " << my_tree->get_T() << " with tmrca " << my_tree->get_all_samples()[n-1]->get_time() << " and S " <<  my_tree->get_S() << "\n";
		
		


		distrib_TMRCA->add_value( my_tree->get_all_samples()[n-1]->get_time()  );
		
		if(opt_scenario == 's')
			distrib_PL->add_value( my_tree->get_nsurvivors()  );




		my_tree->clean_whole_tree(  );                            // remove all mutations
		
	}

	/*
		Output
	*/
	printf("TMRCA\t");
	printf("Mean  = %f\tVar = %f\t", distrib_TMRCA->get_mean(), distrib_TMRCA->get_variance() );
	printf("95%%CI = %f %f\n", distrib_TMRCA->get_low_boundary() , distrib_TMRCA->get_high_boundary() );
		
	if(opt_scenario == 's'){
		printf("Surv_n\t");
		printf("Mean  = %f\tVar = %f\t", distrib_PL->get_mean(), distrib_PL->get_variance() );
		printf("95%%CI = %f %f\n", distrib_PL->get_low_boundary() , distrib_PL->get_high_boundary() );
	}
	printf("Stat\t");
	printf("Mean  = %f\tVar = %f\t", distrib_Val->get_mean(), distrib_Val->get_variance() );
	printf("95%%CI = %f %f\t", distrib_Val->get_low_boundary() , distrib_Val->get_high_boundary() );
	printf("Power = %.3f (%ld/%ld)\n",  power/(double)r, power, r );
		

	/*
		Free Memory
	*/

	delete distrib_Val;

	if(mode == 1){
		delete w2;
		if( inputfiles.get_size_w1() == 0)
			delete w1;
	}
	
	delete std_tree;
	delete distrib_std;

	delete my_tree;   

	return 0;
}
