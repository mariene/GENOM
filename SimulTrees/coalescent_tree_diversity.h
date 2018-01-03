/***
		file:     coalescent_tree_diversity.h
		function: header of a coalescent_tree_diversity.cpp
		author:   <amikezor>
		date:     sep 04
***/

#ifndef _COALESCENT_TREE_DIVERSITY_H_
#define _COALESCENT_TREE_DIVERSITY_H_

#include "coalescent_tree.h"

class coalescent_tree_diversity : public coalescent_tree{

	protected:

		int Haplos;                        // total number of haplotypes
		
		int S;                             // total number of sites
		double K;                          // average number of pairwise difference

		double an;                         // the harmonic SUM(i=1 to n-1){1/i};

		long *Xi_array;                     // a full frequency spectrum

		double T;                          // a given Test
		double alpha;                      // alpha  -- only depends on w1, w2 and n
		double beta;                       // beta   -- can be computed once and stored
		
		double D;                          // a given Test

		double *HarmonicSums;               // Harmonic Sums, needed to improve computation of variances

		double beta_FU1995( int i );
		double sigma_ii_FU1995( int i );
		double sigma_ij_FU1995( int i, int j );

	public:
	
		coalescent_tree_diversity( int n, population *p );
		~coalescent_tree_diversity( void );

		double compute_VarT( double *w1, double *w2, double theta);                     // a function to compute the variance under neutrality of any test
		double compute_CoVarw1w2( double *w1, double *w2, double theta);

		void compute_VarTAlphaBeta( double *w1, double *w2, double *alpha, double *beta);      // a function to compute the alpha/beta of this variance

		void compute_Xi_array_species( int sp );                                              // restrict to a species. if sp==-1, compute for the whole tree
		void compute_Xi_array( void ){ compute_Xi_array_species(-1);};                         // K, S, an, Kno1 and Sno1

		void  compute_T( double *w1, double *w2, char norm);

		void compute_T( double *w1, double *w2, char norm, short compute_alphabeta );      // a generic function that can compute any test (Achaz, 2009)
		void compute_D( char norm, short compute_alphabeta );                              // Tajima 1989     --> K vs S
		void compute_F( char norm, short compute_alphabeta );                              // Fu and Li 1993  --> K vs Xi1
		void compute_Fstar( char norm, short compute_alphabeta );                          // ""              --> K vs Eta1
		void compute_D_fl( char norm, short compute_alphabeta );                           // "    "    "     --> S vs Xi1
		void compute_Dstar_fl( char norm, short compute_alphabeta );                       // "    "    "     --> S vs Eta1
		void compute_H( char norm, short compute_alphabeta );                              // Fay and Wu 2000 --> S vs Eta1
		void compute_Y( char norm, short compute_alphabeta );                              // Achaz      2008 --> K_noXi1 vs S_noXi1
		void compute_Ystar( char norm, short compute_alphabeta );                          // Achaz      2008 --> K_noEta1 vs S_noEta1


		void compute_D_classic( char norm );                                               // Tajima 1989     --> K vs S/an, orginal derivation

		void compute_Haplos( void );                                                        // Number of Haplotypes

		double get_K(void){return K;};
		int get_S(void){return S;};
		long * get_Xi_array(void){return Xi_array;};
		long get_Eta1(void){return Xi_array[0]+Xi_array[n-2];};
		double get_T(void){return T;};
		double get_D(void){return D;};
		double get_an(void){return an;};
		double get_Haplos(void){return Haplos;};
	

		int * compute_Kdistrib( void );
		int **compute_Kmatrix( void );

		void print_summary_statistics(void);           // print out S, Sno1, D, etc...
		void printout_treeWithMut( void  );
		void printout_treeFromMut( void  );
};


inline int * coalescent_tree_diversity::compute_Kdistrib( void ){return this->coalescent_samples[0]->compute_Kdistrib();};
inline int ** coalescent_tree_diversity::compute_Kmatrix( void ){return this->coalescent_samples[0]->compute_Kmatrix();};


#endif
