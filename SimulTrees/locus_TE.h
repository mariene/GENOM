/***
		file:     locus.h
		function: header for locus functions
		author:   <amikezor>
		date:     sep 04
		modif:    mar 09 -- add the id_path
***/


#ifndef _LOCUS_TE_H_
#define _LOCUS_TE_H_

#include "locus.h"

/***
	The locus_transposon class
	each locus_transposon is a node in a transposon tree (nodes are either coalescent of genomes or duplications) --
***/

//class coalescent_tree_diversity : public coalescent_tree{

class locus_TE : public locus
{

	protected:
	
		static int Count_All_transposon;                   // total number of locus_transposon in use
		int locus;

	public:
	
		locus_TE( void );                               // constructor
		~locus_TE( void );                              // destructor
		
		void set_locus( int a ){this->locus=a;};
		int get_locus( void ){return this->locus;}; 

		int IsTE(void){return 1;};

};


#endif
