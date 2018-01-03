/***
		file:     locus.h
		function: header for locus functions
		author:   <amikezor>
		date:     sep 04
		modif:    mar 09 -- add the id_path
***/


#ifndef _LOCUS_H_
#define _LOCUS_H_

#include <stdlib.h>

class locus_TE;

/***
	The locus class
	each locus is a node in a tree --
	a locus is an object typically sharing a history with other loci
***/

class locus{
	
	protected:	
		int id;                             // id is unique and guarantee to be by 'number_of_loci'
		double time;                        // time to the present - 0=present, >=0=backward in time.
		double Ttime;                       // Total time to the bottom of the tree -- if this locus is a leaf, it is 0 --

		locus *ancestor;                    // the ancestral locus of this one
		locus *descendant1;                 // and its two descendant;
		locus *descendant2;

		static int Count_All_Loci;          // total number of loci in use
				
		int nleaves;                        // the tree level (how many leaves under this node ?)
		int newSites;                       // number of new polymorphic sites from this node

		int species;                        // a tag that specifies the species it belongs to (used e.g. for speciation_coalescent)
		int nleaves_sp;                     // the species tree level (how many leaves of a given species under this node ?)

		char *sequence;                     // the sequence of this locus (often only polymorphic site)
		int size_sequence;                  // the size of that seq
		int size_seq_mem;                   // the size of the memory allocation for that seq


		int *idpath;                        // path up to the root (including this one)
		int n_idpath;                       // number of nodes up to the root (including this one)
		
		locus *pick_node( double Tr );                           // pick a node according to a time Tr [between 0 and Ttime]


		locus_TE *array_TE;
		int size_array_TE;
		int n_TE;

	public:
		locus( double time=0);                                   // constructor
		~locus( void );                                          // destructor

		virtual int IsTE( void ){return 0;};
		virtual int get_locus(void ){ return -1;};

		int get_id( void ){return this->id;};                    // get the id;

		void set_time( double t ){this->time=t;};                // set the time where this locus exist;
		double get_time( void ){return this->time;}; 	         // get this time;

		void set_Ttime( double t ){Ttime=t;};                    // set the Total time to the bottom;
		double get_Ttime( void ){return this->Ttime;}; 	         // get this Ttime;
		void compute_Ttime( void );                              // compute T_time as a function of the descendants (if any)

		void set_species( int s ){this->species=s;};             // set the species;
		int get_species( void ){return this->species;}; 	 // get this species;


		int get_newSites(void){return this->newSites;};          // get the number of sites newly mutated at this node

		void compute_nleaves( void );                            // set the number of leaves under this node.
		int get_nleaves(void){return this->nleaves;};            // get the number of leaves under this node.

		int get_nleaves_sp( int sp );                            // compute and get the number of leaves under this node with species tag to sp.
		void compute_nleaves_sp( int sp );
		
		int get_n_idpath(void){return n_idpath;};
		int *get_idpath(void){return idpath;};
		
		void compute_idpath(void);
		
		void connect_ancestor( locus *anc , char descendant);    // connect with an ancestral locus.
		                                                         // descendant has to be 1 or 2 

		void reset_connections( void ){ancestor=descendant1=descendant2=NULL;}
		void reset_descendants( void ){descendant1=descendant2=NULL;}
		void reset_ancestor( void ){ancestor=NULL;}

		locus *get_ancestor(){ return this->ancestor; };         // extract the ancestral locus
		locus *get_descendant1(){ return this->descendant1; };   // extract the descendant locus
		locus *get_descendant2(){ return this->descendant2; };   // extract the other descendant

		bool is_coalesced();                                     // is this locus already coalesced ?
		
		void init_sequence( int len, char value );               // get a tring of size 'len' and init it to 'value'
		char *get_sequence( void ){return sequence;};            // get the sequence
		int get_size_sequence( void ){return size_sequence;};    // get the sequence

		void mutate_subtree( double theta, int S, double mean_err ); // add sites at all tree node under this node according to Ttime*theta/2
		void clean_subtree( void );                                  // set recursively all newSites to 0; 
		
		void mutate_subtree_seq( double theta, int L, int errors );         // create mutated sequences all tree under this node according to Ttime*theta/2
		void mutate_seq( int n, char value );                        // mutate the sequence at base n (as well as its descendants) into 'value'

		locus_TE *get_array_TE( void ){return array_TE;}            // set/get the sub-loci, i.e. transposon inside a haplotype
		void set_array_TE( locus_TE *l ){ array_TE = l; }
		int get_size_array_TE( void ){return size_array_TE;}        // set/get size of the array
		void set_size_array_TE( int a){size_array_TE = a;}		
		int get_n_TE( void ){return n_TE;}                // set/get the number of elements at the time of this locus (that is <= size of the array)
		void set_n_TE( int a ){n_TE = a;}

};


#endif
