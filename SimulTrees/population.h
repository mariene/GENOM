/***
		file:     population.h
		function: header of a population.cpp
		author:   <amikezor>
		date:     sep 04-05
		modif:    nov 08 - change to multispecies the isolation and allow symetric migration.
***/

#ifndef _POPULATION_H_
#define _POPULATION_H_

#include <iostream>
using namespace std;

#include <stdlib.h>

struct isolation_event {

	double time;   // the bckwd time to the event
	int sp[2];     // the species names that are merged. For now, we do not check what user sets.
	int id_sp;     // the id of the species in which sp[0] and sp[1] are merged

};


class population{


	private:
		double Ne;                 /* the effective population size */
		double mu;                 /* the mutation rate per generation */
		char ploidy;               /* either 1 (haploid) or 2 (diploid) */

		double Theta;              /* ploidy*2*Ne*mu -- scaled mutation/coalescence rate */


		/*
			The following are used for bottlenecks
			Notations follow Simonsen et al., Genetics, 1995
			with the exception that time is rescaled by mu
			
			|           |
 			|__       __|
			   |     |        Tl - time of the bottleneck
			 __|     |__      _   Reduction to size  f.N 
			|           |
			|           |     Tb - time to the bottleneck
			|           |
			
		*/

		double Tb;               // the time since the bottleneck ended  (in ploidy*N*mu=theta/2 generations)
		double Tl;               // time length of the bottleneck  (in ploidy*N*mu=theta/2 generations)
		double f;                // strength of the bottlenexk. The population was shrinked to f*2N individuals


		/*
			The following are used for isolation/speciation (with migration m)
			time is rescaled by mu

			     |          |         Ne - size of the true ancestral species, numbered 2*ns;
			     |    3     |  
			 ____|   _____   |       ----> isolation/speciation event Ts(0,4)
			|     |      |   |
			|   |     ___| 4 |__     ---> isolation/speciation event Ts(1,2) -> f[4]*Ne is the size of this anc species
			|   |    |    __    |     
			| 0 |   |  1 |  | 2 |     ns species at present time. individual migrate with M/2
		
		 Ne*    f[0]     f[1]   f[2]        - Ne*f[i] - size of species i

		*/

		long ns;                        // number of species - tagged from 0 to ns-1.
		struct isolation_event *Ti;     // time to the isolation event between two species. there will be at most (ns-1) such events.
		double *fN;                     // an array of 'relative' size when compared to the anc species that is Ne individual. len(f): ns+ns-2;
		double M;                       // coalescent migration rate 2pNm. symetric rate. P(migration bkwd) = M/2*(ns-1)/ns


		/*
			The following are used for exponential/linear growth // or decline (when Gr<0).

		        |\  
		       |  \            // size at time T was N(T) = N * exp( -Gr * T )  or N(T) = N( 1 - a t)
		     _|    \_ 
		  __|         \__
		 |               \      // size now is N

		*/

		double Gr;              // population growth rate  (when positive) or decline (when negative)
		double slope;           // in case of a linear growth (positive) or reduction (negative)


		/*
			The following are used for selective sweep with recombination


			  |   Ne   |      Ne - size of the pop
			  |\       |
			  |  \     |      sweep, with selec coef alpha=2pNs (s select coef)
			  |    \   |                  rec   coef  R = 2pNr  (r rec coef)
			  |     \  |
			  |       \| <---- TS - when sweep has started
			  |        |
			  |        |

		*/
		
		double TS;              // time back to the end of the sweep
		double R;               // recombination rate in p*Ne*r (r being per generation)
		double alpha;           // selection coef. in p*Ne*s (s being the usual 's')
		double p;               // frequency of the selected allele when sample is done
		bool sweep_over;        // becomes TRUE when sweep is over (when pS=0)
		
		

	public:
		population(double Ne=0, double mu=0, char ploidy=1);
		~population(void);
		
		void free_population_arrays( void );
		
		void set_Ne( double n ){Ne=n;};
		void set_mu( double m ){mu=m;};
		void set_ploidy( char p );

		
		double get_Ne( void ){return Ne;};
		double get_mu( void ){return mu;};
		char   get_ploidy( void ){return ploidy;};

		void set_Theta( double t ){Theta=t;};
		void compute_Theta( void );
		double get_Theta( void ){return Theta;};


		/*
			The following are used for bottlenecks
			Notations follow Simonsen et al., Genetics, 1995
		*/
		void set_bottleneck( double Tb, double l, double f ){this->Tb=Tb; this->Tl=l; this->f=f;};
		double get_bottleneck_Tb( void ){return this->Tb;};
		double get_bottleneck_Tl( void ){return this->Tl;};
		double get_bottleneck_f( void ){return this->f;};



		void set_structure(  double *fN, long ns, double M );     // Set a chosen scenario of species tree



		/*
			The following are used for isolation/speciation
			expanded into multi-species in Nov 2008
		*/
		
		void set_isolation( struct isolation_event *Ti, double *fN, long ns, double M );     // Set a chosen scenario of species tree
		
		void BirthDeath_species_tree(  long ns, double M, double b, double d, double *input_fN, char ancestral_type );    // built random isolation events (species tree) from birth & death process
		void Moran_species_tree(  long ns, double M, double r, double *input_fN, char ancestral_type );                  // built random isolation events (species tree) from Moran coalesc process
		void Radiation_species_tree(  long ns, double M, double r, double *input_fN, char ancestral_type );              // built random isolation events (species tree) from Radiation process
		void Fixed_species_tree(  long ns, double M, double t, double *input_fN, char ancestral_type );                  // built fixed isolation events (species tree) at t, 2t, 3t, etc...

		/*
			Assume you have a sample from all species in the clade
		*/
		void Kingman_species_tree(  long ns, double M, double r, double *input_fN, char ancestral_type );                 // built random isolation events (species tree) from Kingman coalesc process


		void print_isolations( void );
		
		
		struct isolation_event get_isolation_Ti( int i ){if(Ti == NULL){cerr <<"no Ti is defined\n"; exit(1);}else{return Ti[i];}};
		struct isolation_event *get_isolation_Ti( void ){return Ti;};
		double get_isolation_fN( int i ){if(fN == NULL){cerr <<"no fN is defined\n"; exit(1);}else{return this->fN[i];}};
		double *get_isolation_fN( void ){return this->fN;};
		long get_isolation_ns( void ){return this->ns;};
		double get_isolation_M( void ){return this->M;};
		
		int get_isolation_next_Ti( double time );                              // return the id of the next Ti > time




		/*
			The following are used for exponential/linear growth
		*/
		void set_exponential_growth( double Gr ){this->Gr=Gr;};
		double get_exponential_growth_Gr( void ){return this->Gr;};

		void set_linear_growth( double Gr ){this->slope=Gr;};
		double get_linear_growth_Gr( void ){return this->slope;};


		/*
			The following are used for distant sweep
			Notations follow Braverman et al., 1995, Genetics
		*/
		void set_sweep( double TS, double p, double alpha, double R ){this->TS=TS; this->p=(TS>0)?1:p; this->alpha=alpha; this->R=R; sweep_over=false;};
		double get_sweep_TS( void ){return this->TS;};
		double get_sweep_alpha( void ){return this->alpha;};
		double get_sweep_R( void ){return this->R;};
		double get_sweep_p( void ){return this->p;};
		bool is_sweep_over( void ){return this->sweep_over;};
		void set_sweep_over( bool state){this->sweep_over = state;};

};
#endif
