/***
		file:     coalescent_tree.h
		function: header of a coalescent_tree.cpp
		author:   <amikezor>
		date:     sep 04
		modif:    nov 08
***/

#ifndef _COALESCENT_TREE_H_
#define _COALESCENT_TREE_H_

#include "sample.h"
#include "locus.h"
#include "population.h"

class coalescent_tree{

	protected:
	
		population *pop;               // in what kind of population is this tree built
		sample **coalescent_samples;   // Each sample is a layer ; current one will be 0, the mrca will be (n-1)
		int n;                         // number of individuals in the sample
		
		int nsurvivors;                // number of individuals that survived an event (ie.e sweep/bottleneck, etc.)
		
		locus *TE_anc_locus;
	
	public:
	
		coalescent_tree( int n, population *p );
		~coalescent_tree();		

		int get_n(){return this->n;};
		population *get_pop(){return this->pop;};
		
		sample **get_all_samples(){return this->coalescent_samples;};
		sample *get_current_sample(){return this->coalescent_samples[0];};
		sample *get_oldest_sample(){return this->coalescent_samples[0]->oldest_sample();};

		locus *get_mrca_species( int sp );

		double get_Tmrca(void){return coalescent_samples[n-1]->get_loci()[0]->get_time();};
		double get_Tmrca_species( int sp ){return get_mrca_species( sp )->get_time();};
		
		locus * get_mrca_2loci( locus *l1, locus *l2 );
		locus * get_mrca_nloci( locus **list_locus, int size_list );
		
		double get_Ttot(void){return coalescent_samples[n-1]->get_loci()[0]->get_Ttime();};
		double get_Ttot_singletons(void);

		int get_nsurvivors(void){return nsurvivors;};
		void set_nsurvivors(int x){nsurvivors=x;};


		double get_Tcoal(int ncoal){return coalescent_samples[ncoal]->get_loci()[0]->get_time();};


		/*
			generate tree(s)
		*/
		void standard_coalescent( void ){ standard_coalescent(1, NULL, NULL); };  // the most classical standard coalescent
		void standard_coalescent( int nsamples, int *sample_sizes, double *serial_times );      // including serial samples aka time series

		void fixedtimes_coalescent( double *times );           // assuming ERM but with fixed times. Usefull for drawings
		
		
		void bottleneck_coalescent( void );
		void bottleneck_reduce_tree( void );
		void isolation_coalescent( int *ni );    // ni[0] contains the number of seq in species 0. Its size
		                                         // is given by mypop->get_ns. Its sum is n, the total sample size
							
		void sweep_coalescent( void );
		void growth_coalescent( char type );     // type is either 'l'inear or 'e'xponential
		void nested_coalescent( int k );         // k is the size of the founder mutation
		
		void bs_coalescent( void );              // Bolthausen-Sznitman

		void structured_coalescent( int *ni ){ structured_coalescent( 1, &ni, NULL ); };   // a classical structured coalescent
		void structured_coalescent( int nsamples, int **samples, double *serial_times );        // with serial samples


		/* add mutations */
		void mutate_whole_tree( void ){ mutate_whole_tree(0,0.0); };                            // default is draw S and set SeqErr to 0
		void mutate_whole_tree( double mean_err ){ mutate_whole_tree( (int)0 ,mean_err); };     // draw S from Poisson(Theta*L/2) and SeqErr from Poisson(mean_err)
		void mutate_whole_tree( int S, double mean_err );                                       // S is fix, draw SeqErr from Poisson(mean_err)
		void clean_whole_tree( void );                                                          // remove all mutations from the tree
		
		
		void mutate_whole_tree_TE( void ){ mutate_whole_tree(0,0.0); }; 
		void mutate_whole_tree_TE( double mean_err ){ mutate_whole_tree( (int)0 ,mean_err); };  // draw S from Poisson(Theta*L/2) and SeqErr from Poisson(mean_err)
		void mutate_whole_tree_TE( int S, double mean_err );                                    // S is fix, draw SeqErr from Poisson(mean_err)
		
		void mutate_whole_tree_seq( void ){mutate_whole_tree_seq(0, 0);};  
		void mutate_whole_tree_seq( int L ){mutate_whole_tree_seq(L, 0);};  
		void mutate_whole_tree_seq( int L, int errors );                                        // create mutated sequences (polymorphic sites)

		/* output */
		char **extract_sequences( void );
		char **merge_sequences( char **input_sequences );

		void printout_tree( void );
		void printout_tree_TE( void );
		void printout_sequences( void );
		void printout_all_sequences( void );
		void printout_nsites( void );
		
		
		/* T.E. */
		void coal_TE_TopDown_mrca( int n_TE_mrca, float dupl_rate);
		void coal_TE_TopDown_FromMrca( float dupl_rate );
		void coal_TE_TopDown_1stTE( float time_1st_TE, float dupl_rate);
		
};

#endif
