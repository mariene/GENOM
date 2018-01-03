/***
		file:     sample.h
		function: header for the sample functions
		author:   <amikezor>
		date:     sep 04
***/

#ifndef _SAMPLE_H_
#define _SAMPLE_H_

#include "locus.h"
#include "locus_TE.h"
#include "population.h"

/***
	The sample class
	each sample is a layer of the process (one for each time).
	There is a "new" sample each time 2 loci coalesces.
	So a sample has a FIXED number of loci 
***/

class sample{
	
	private:
		int id;                           // its id number
		double time;                      // time where this sample start existing. 0: present, >=0: backward in time;

		locus **loci;                     // a list of the current loci in that sample (at that time)
		int     nloci;                    // the number of them

		int     n_original_loci;          // the number of locus that were allocated here; usefull when freeing the memory
		
		int nspecial;                     // if one need to count some special lineages

		sample *ancestor;                 // its ancestral sample (otherwise NULL)  
		sample *descendant;               // its descendant sample (otherwise NULL) 

		population *pop;                  // the population this sample has been taken from
	
		static int Count_All_Samples;     // Keep track of all samples in use.

		void coalesce2loci( int n1, int n2, double coal_time);
		void coalesce2loci_species( int nl1, int nl2, double coal_time, int _chosen_type_);

		int get_n_alive_species( int *ni );
		int get_alive_species_x( int x, int *ni );
		int get_random_alive_species( int *ni );

		locus *get_species_locus_x( int species, int x );                        // select locus number x from a chosen species -- see bottom of this file
		locus *get_species_locus_x( int species, int x, double cuurent_time );   // select locus number x from a chosen species;
		
		locus *get_random_locus_species( int _chosen_species_ , double current_time); // pick randomnly a locus in species _chosen_species_
		locus *get_random_locus_species( int _chosen_species_ ); // pick randomnly a locus in species _chosen_species_
		

		int get_random_species_bysize( void );

	public:
		sample( population *p );                                // constructor;
		~sample( void );                                        // destructor;
	
		void set_time(double t){this->time=t;}                 // set time of sample initiation 
		double get_time(void){return this->time;}              // get this time
		
		int get_id(void){return this->id;}                     // get the id of the sample
 
		locus **get_loci(void){return this->loci;}             // get the current of them
		int get_nloci(void){return this->nloci;}               // how many loci in this sample
		void set_loci( int n, int n_memory );                     // set an array of n loci and "new" for the n_mem first ones
		int get_nloci_befT(double T);                            // how many loci in this sample with their time is < T
		double get_next_serial_sample( double time );


		int get_nloci_species( int species  ){ return get_nloci_species(species, this->time); };        // count the number of loci of a species ;
		int get_nloci_species( int species, double current_time  );                              // count the number of loci of a species ;

		locus **get_loci_species( int species  ){ return get_loci_species(species, this->time); };        // count the number of loci of a species ;
		locus **get_loci_species( int species, double current_time  );                   // count the number of loci of a species ;

		int get_nspecial(void){return  nspecial;};              //  set/get the nspecial; i.e. for now used only to count the sweep survivors 
		void set_nspecial(int x){ nspecial=x;};


		sample *get_ancestor( void );			        // return the ancestral sample   inline 
		sample *get_descendant( void ); 		        // return the descendant sample  inline 
		void connect_ancestor(sample *s);                       // connect this sample to an ancestor 
		
		sample *oldest_sample();                                // go back to the oldest sample and return it

		char **get_sequences();                                 // get all sequences (aligned) from this sample


		/*
			Coalescence scenarios
		*/
		void coalesce_sample( void  );                          // coalesce 2 loci into ancestor[0] then copy pointers to other locations;
		void bottleneck_reduce_time_sample( void );             // reduction in real time for a sample time + locus[0]
		void growth_reduce_time_sample( char type  );           // reduction in real time for a sample time + locus[0]. type is either 'e'xponential or 'l'inear
		void isolation_coalesce_sample( void  );                // coalesce 2 loci from a species into ancestor then copy pointers to other locations;
		void sweep_coalesce_sample( void );                     // coalesce either during sweep or in neutral
		void nested_coalesce_sample( void );                    // coalesce with one subtree having species tag set to 0
		sample *bs_coalesce_sample( void  );                    // coalesce k among n loci into ancestor[0] then copy pointers to other locations;
		void structured_coalesce_sample( void );                // a classical subpopulations setup

		void coalesce_sample_fixedtime( double time  );             // coalesce 2 loci into ancestor[0] then copy pointers to other locations;



		/*
			Pairwise differences
		*/
		int * compute_Kdistrib( void );                         // return an array with all 2 by 2 differences
		int **compute_Kmatrix( void );                          // return a nxn matrix with all differences
		int compute_K( int l1, int l2);
		
		/*
			TE --transposons managment functions
		*/
		void set_locus_TE( locus *l, int nTE, int nTE_anc );
		void coalesce_TE_locus( double time, int *alive_locus, int n_alive, locus_TE *array_TE, int size_array_TE, int te_anc );

		void coalesce_TE_mrca_BottomUp( float dupl_rate, int n_te_mrca, locus *mrca );
		void coalesce_1st_TE_locus_TopDown( float dupl_rate, float time_1st_TE, locus *l  );      
		
		void coalesce_TE_locus_TopDown( float dupl_rate, locus *l  );      
		void coalesce_TE_sample_TopDown( float dupl_rate );

		void print_all_TE_connections( void );

};

inline sample *sample::get_ancestor(void){ return this->ancestor; };
inline sample *sample::oldest_sample(void){ sample *p = this; while( p->get_ancestor() )p=p->get_ancestor(); return p; };
inline sample *sample::get_descendant( void ){ return this->descendant; };

inline locus *sample::get_species_locus_x( int species, int x ){return get_species_locus_x( species, x, this->time ); };
inline locus *sample::get_random_locus_species( int _chosen_species_  ){return get_random_locus_species(_chosen_species_, this->time ); };

#endif
