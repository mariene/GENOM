/***
		file:     coalescent_tree_diversity.cpp
		function: all functions associated with the computing diversity statictics without sequences
		author:   <amikezor>
		date:     nov 06-08
		          nov 08 - switch to generic omega tests
			  sept 09 - handle n=1 case (no tree)
			  dec 11 - bugs in compute_Xi_array()
***/


#include <iostream>
using namespace std;


#include <stdlib.h>
#include <stdio.h>
#include "sample.h"
#include "population.h"
#include "random.h"
#include "math.h"


#include "coalescent_tree_diversity.h"


/*
	Constructor
*/
coalescent_tree_diversity::coalescent_tree_diversity( int n, population *p ):coalescent_tree( n, p){

	Haplos = 0;
	S = 0;
	K = 0.0;
	T = 0;
	alpha = 0;
	beta = 0;
		
	Xi_array=(long *)calloc( (size_t)(n-1),sizeof(long) );
	if(Xi_array == 0){
		cerr << "Not enough memory to compute the frequency spectrum, bye\n";
		exit(1);
	}


	HarmonicSums = new double[n];
	if(HarmonicSums == NULL)
		fprintf(stderr, "coalescent_tree_diversity::coalescent_tree_diversity: cannot get memory for HarmonicSums, bye\n"), exit(3);

	HarmonicSums[0]=0;
	for(int i=1;i<n;i++)
		HarmonicSums[i] = HarmonicSums[i-1] + 1.0/i;

	an = HarmonicSums[n-1];
	
	
//	cerr << "//\nCreated a sample of id "<<this->id<<" and time "<<time <<"\n//\n";
}

coalescent_tree_diversity::~coalescent_tree_diversity( void){
	free(this->Xi_array);
	delete this->HarmonicSums;
}



/*
	When sp is set to -1, ignore the species tag. Otherwise
*/
void coalescent_tree_diversity::compute_Xi_array_species( int sp ){

	int i,j;
	int nleaves;   // on each node, how many descendants
	int nsites;    //  how many mutations in the branching leading to it


	int TreeSize;
	
	
	/*
		Simply go down and set the nleaves_sp to the number
		of leaves of species sp under this node
	*/
	this->get_mrca_species( -1 )->compute_nleaves_sp( sp );
	
	if(sp == -1)
		TreeSize = n;
	else
		TreeSize = get_mrca_species( sp )->get_nleaves_sp( sp );
		
	
	this->S=0;
	this->K=0;
			
	for(i=0;i<this->n-1;i++)
		this->Xi_array[i]=0;
		
	/*
		1. Do all internal nodes (sample 1 to n-1, thus excluding the ROOT node)
	*/
	for(j=1;j<this->n-1;j++){
		
		
		nleaves = coalescent_samples[j]->get_loci()[0]->get_nleaves_sp( sp );             /* this computes the number of leaves of same sp (or sp==-1) in the subtree */
		nsites  = coalescent_samples[j]->get_loci()[0]->get_newSites();                   /* number of sites on the branch that led to this node */

		if( sp != -1 && (nleaves == 0 || nsites == 0) )
			continue;
		
		this->S += nsites;
		this->K += nleaves*( TreeSize - nleaves) * nsites;

		this->Xi_array[  nleaves-1 ]+= nsites;



	}
	
//	printf("now singletons\n");

	/*
		Add the singletons (sample 0)
	*/
	for(i=0;i<this->n;i++){
		
		nsites =  coalescent_samples[0]->get_loci()[i]->get_newSites();
		
		if( sp != -1 && ( coalescent_samples[0]->get_loci()[i]->get_species() != sp || nsites == 0 ) )
			continue;
	
		this->S += nsites;
		this->K += nsites*(TreeSize-1);
		this->Xi_array[ 0 ]+= nsites;

	}

	
	this->K *= 2.0 / (double)( TreeSize*(TreeSize-1.0) ); 

}

/*
	This is used to compute variances and covariances
	of Xi in Fu 1995 Theor Pop Biol
*/
double coalescent_tree_diversity::beta_FU1995( int i  ){

	double ai=this->HarmonicSums[i-1];
	double beta=0;
	double nloci= (double)n;
	
	

	beta = 2.0 * nloci * ( an + (1.0/nloci) - ai )/(  (nloci-i+1.0 )*(nloci-i) ) - 2.0/(nloci - i);

	return beta;
}

double coalescent_tree_diversity::sigma_ii_FU1995( int i ){
	
	double nloci= (double)n;
	double sigma_ii=0;
	double ai=this->HarmonicSums[i-1];

	if( 2*i < n )
	{
		sigma_ii = beta_FU1995( i+1  );
	}
	else
	{
		if( 2*i == n  )
		{
		
			sigma_ii = 2.0*(an-ai)/(nloci - i) - 1.0/(i*i);

		}
		else
		{
			sigma_ii = beta_FU1995( i   ) - 1.0/(i*i);
		}

	}
	
	return sigma_ii;
}

double coalescent_tree_diversity::sigma_ij_FU1995( int i, int j ){
	
	double nloci= (double)n;

	double sigma_ij=0;

	if(i==j){
		cerr << "defined only for i>j, bye\n";
		printf("here i= %d, j=%d\n", i, j);
		exit(3);
	}
 	if(i<j){
		int tmp=i;
		i=j;
		j=tmp;
	}
 
	double ai=this->HarmonicSums[i-1],
	       aj=this->HarmonicSums[j-1];


	if( i+j < n )
	{
		sigma_ij = ( beta_FU1995( i+1  ) -  beta_FU1995( i  ) ) / 2.0;
	}
	else
	{
		if( i+j == n  )
		{
		
			sigma_ij  =  ((an - ai)/(nloci - i) + (an - aj)/(nloci - j)) 
			           - ( ( beta_FU1995( i   ) +  beta_FU1995( j+1  ) )  / 2.0 ) 
				   - (1.0/(i*j));

		}
		else
		{
			sigma_ij = (( beta_FU1995( j  ) -  beta_FU1995( j+1 ) )/2.0) - (1.0/(i*j));
		}

	}
	
	return sigma_ij;
}




/*
	Gives two weight vectors
	and if T should be normalized or not
	if compute_alphabeta is set to 1 set alpha and beta, otherwise use previous ones (if none compute)
*/

void coalescent_tree_diversity::compute_T( double *w1, double *w2, char norm, short compute_alphabeta ){


	double Sum_w1=0.0;
	double Sum_w2=0.0;
	int i,j;

	double alpha, beta;

	this->T = 0;
	
	if ( S == 0 )
		return;

	for(i=0;i<n-1;i++){
		 Sum_w1 += w1[i];
		 Sum_w2 += w2[i];
	}


	if(Sum_w1 == 0 || Sum_w2 == 0)
		return;
		
		
		
	for(i=0; i < this->n -1 ; i++)
		this->T += (double)(i+1) * Xi_array[i]* ( (w1[i]/Sum_w1) - (w2[i]/Sum_w2) );
	
	
//	printf("No Norm %f\n", this->T);
	
	if(norm == 0)
		return;

	/*
		Compute the alpha and beta of the variance
	*/

	//if( compute_alphabeta == 1 || this->alpha == 0 || this->beta == 0 ){
	
		alpha=0;
		beta=0;
		for(i=0; i < this->n - 1 ; i++){
	
			double omega = ( (w1[i]/Sum_w1) - (w2[i]/Sum_w2) );
		
			alpha += omega*omega*(i+1.0);
		}
	
		//cerr <<"alpha is "<<alpha<<"\n";
	
		beta=0;
	
	
		for(i=1; i < this->n ; i++){
			for( j=i; j< this->n ; j++){
		
				double omega_i = ( (w1[i-1]/Sum_w1) - (w2[i-1]/Sum_w2) );
				double omega_j = ( (w1[j-1]/Sum_w1) - (w2[j-1]/Sum_w2) );
		
		
				if( i == j)
					beta  +=  i*i* omega_i*omega_i * sigma_ii_FU1995( i  );
				else
					beta  +=  2*i*j* omega_i*omega_j * sigma_ij_FU1995( i, j  );
						
			}
		}
		//cerr <<"beta is "<<beta<<"\n";
	//}
	
 	double Estim_theta=0.0;         /* estimation of Theta */
 	double Estim_theta2=0.0;        /* estimation of Theta^2 */
	double a2=0;
	
	for(i=1;i<this->n;i++){
		a2 += 1.0/(double)(i*i);
	}


	Estim_theta =  (double) this->S / this->an;
	Estim_theta2 =  (double) this->S*( (double)this->S-1.0 )/(this->an*this->an +a2 );

	this->T /= sqrt( alpha*Estim_theta + beta*Estim_theta2   );

/*
	printf("Norm is %f\n", alpha*Estim_theta + beta*Estim_theta2 );
	printf("T norm is %f\n", this->T / sqrt( alpha*Estim_theta + beta*Estim_theta2   ) );
	printf("T norm is %f\n", this->T );
*/

}




/*
 *
 *
	All the following are
	std tests developped in the past
	
	D, F, F*, D_fl, D*_fl, H, Y, Y*
*
*/

// K - S
void coalescent_tree_diversity::compute_D( char norm, short compute_alphabeta ){

	double *w1 = new double[ this->n -1];
	double *w2 = new double[ this->n -1];


	for(int i=1; i<= this->n -1 ; i++){
	
		w1[i-1] = (double)this->n - i;
		w2[i-1] = 1.0/(double)i;
	
	}

	compute_T( w1, w2, norm, compute_alphabeta );

	delete w1;
	delete w2;
}



/*
	S is the number of polymorphic sites
	K is the mean differences between sequences
	D is the ratio between d and its std dev
	Computed as suggested by Tajima 83
	
	if norm == 0, then return unormalized value
	
*/
void coalescent_tree_diversity::compute_D_classic( char norm ){


	double d=0.0;            /* d is K - S/an (if it is neutral should be 0) */
	double Var_d=0.0;        /* Var_d is the variance of d */
	 
	double Nseq = (double)n;
		
	double b1=0.0,
	       c1=0.0;

	double a2=0.0,
	       b2=0.0,
	       c2=0.0;
	

 	double Estim_theta=0.0;         /* estimation of Theta */
 	double Estim_theta2=0.0;        /* estimation of Theta^2 */


	if( norm != 0 && norm != 1){
		cerr <<"coalescent_tree_diversity::compute_D takes either 0 or 1 as argument, bye\n";
		d=0;
		return;
	}
	
	/*
		no segregating sites
	*/
//      cerr << "S " << S << "\n";
	if ( S==0 ){
	
		d=0;
		return;
	}
//	else
//		cerr << "K "<<K<<" S "<<S<<" an "<<an << "\n";


	d = K - (double)S/an;
	

	/*
		If the statistics is not normalized, do not even compute the variance
	*/
	if( norm == 0 ){
		this->D = d ;
		return;
	}


	/*
		get a1 through e2
	*/
	for(int i=1;i< n;i++){
		a2 += 1.0/(double)(i*i);
	}
	
	b1 = (Nseq + 1.0) / (double)( 3.0 * (  Nseq - 1.0 )  );
	c1 = b1 - (1.0 / this->an);

	b2 = 2.0 * (double) (Nseq*Nseq + Nseq + 3.0) / ((double)(9.0 * Nseq * (Nseq-1.0) ) );
	c2 = b2 - (Nseq + 2.0)/(an * Nseq) + a2/( an * an);


	
//	printf("Classic D: alpha= %f ;  beta= %f\n", c1, c2);
		

	/*
		Get D
	*/

	
	Estim_theta =  (double)S / an;
	Estim_theta2 =  (double)S*( (double)S-1.0 )/(an*an +a2 );

	Var_d = c1* Estim_theta   + c2* Estim_theta2;        // estimate theta from S
	
	this->D = d / sqrt(Var_d);

/*
      cerr << "D " << D << "\t";
      cerr << "d " << d << "\t";
      cerr << "Estim_theta " << Estim_theta << "\t";
      cerr << "Estim_theta2 " << Estim_theta2 << "\t";
      cerr << "Var_S " << an*10+a2*100 << "\t";
      cerr << "Var_D " << Var_d << "\t";
      cerr << "sqrt(Var_d2) " << sqrt(Var_d) << "\n\n";
*/

}






// K - Xi1
void coalescent_tree_diversity::compute_F( char norm, short compute_alphabeta ){

	double *w1 = new double[ this->n -1];
	double *w2 = new double[ this->n -1];

	for(int i=0;i< this->n -1 ; i++){
	
		w1[i] = (double)(this->n-(i+1.0));
		w2[i] = 0;
	
	}
	w2[0]=1.0;

	compute_T( w1, w2, norm, compute_alphabeta );

	delete w1;
	delete w2;
}

// K - Eta1
void coalescent_tree_diversity::compute_Fstar( char norm, short compute_alphabeta ){

	double *w1 = new double[ this->n -1];
	double *w2 = new double[ this->n -1];

	for(int i=0;i< this->n -1 ; i++){
	
		w1[i] = (this->n-(i+1.0));
		w2[i] = 0.0;
	
	}
	w2[0]=(this->n-1.00);
	w2[this->n-2]=1.0;

	compute_T( w1, w2, norm, compute_alphabeta );

	delete w1;
	delete w2;
}

// S - Xi1
void coalescent_tree_diversity::compute_D_fl( char norm, short compute_alphabeta ){

	double *w1 = new double[ this->n -1];
	double *w2 = new double[ this->n -1];

	for(int i=0;i< this->n -1 ; i++){
	
		w1[i] = 1.0/(double)(i+1.0);
		w2[i] = 0;
	
	}
	w2[0]=1.0;

	compute_T( w1, w2, norm, compute_alphabeta );

	delete w1;
	delete w2;
}
// S - Eta1
void coalescent_tree_diversity::compute_Dstar_fl( char norm, short compute_alphabeta ){

	double *w1 = new double[ this->n -1];
	double *w2 = new double[ this->n -1];

	for(int i=0;i< this->n -1 ; i++){
	
		w1[i] = 1.0/(double)(i+1.0);
		w2[i] = 0;
	
	}
	w2[0]=(this->n-1.0);
	w2[this->n -2]=1.0;

	compute_T( w1, w2, norm, compute_alphabeta );

	delete w1;
	delete w2;
}
// S - Eta1
void coalescent_tree_diversity::compute_H( char norm, short compute_alphabeta ){

	double *w1 = new double[ this->n -1];
	double *w2 = new double[ this->n -1];

	for(int i=0;i< this->n -1 ; i++){
	
		w1[i] = (double)(this->n-(i+1.0));
		w2[i] = (double)(i+1.0);
	
	}

	compute_T( w1, w2, norm, compute_alphabeta );

	delete w1;
	delete w2;
}

// K_noXi1 - S_noXi1 
void coalescent_tree_diversity::compute_Y( char norm, short compute_alphabeta ){

	double *w1 = new double[ this->n -1];
	double *w2 = new double[ this->n -1];

	w1[0]=0;
	w2[0]=0;
	for(int i=1;i< this->n -1 ; i++){
	
		w1[i] = (this->n-(i+1.0));
		w2[i] = 1.0/(double)(i+1.0);
	
	}

	compute_T( w1, w2, norm, compute_alphabeta );

	delete w1;
	delete w2;
}

// K_noEta1 - S_noEta1 
void coalescent_tree_diversity::compute_Ystar( char norm, short compute_alphabeta ){

	double *w1 = new double[ this->n -1];
	double *w2 = new double[ this->n -1];

	w1[0]=w1[this->n - 2]=0;
	w2[0]=w2[this->n - 2]=0;
	for(int i=1;i< this->n -2 ; i++){
	
		w1[i] = (this->n-(i+1.0));
		w2[i] = 1.0/(double)(i+1.0);
	
	}

	compute_T( w1, w2, norm, compute_alphabeta );

	delete w1;
	delete w2;
}





void coalescent_tree_diversity::print_summary_statistics( void ){

	cout << "an: " << an << "\n";
	cout << "S: " << S << "\t K:" << K << "\n";

	for(int i=0;i<this->n-1;i++)
		cout << Xi_array[i] << "\t";

	cout << "\n\n";	
	
}



/*
	Compute the Variance of an Omega or an omega vector
*/
double coalescent_tree_diversity::compute_VarT( double *w1, double *w2, double theta){

	double Sum_w1=0.0;
	double Sum_w2=0.0;

	int i,j;

	double alpha=0,
	       beta=0;
	
	int tag=0;
	
	/*
		Compute the alpha and beta of the variance
	*/

	for(i=0;i<n-1;i++){
		 Sum_w1 += w1[i];
		 Sum_w2 += w2[i];
		 if(Sum_w2 != 0)tag=1;
	}

	if(tag == 0)   /* the second vector is null, so report the variance of the first vector (ie. a theta estimator */
		Sum_w2=1;

	alpha=0;
	beta=0;
	for(i=0; i < this->n - 1 ; i++){
	
		double omega = ( (w1[i]/Sum_w1) - (w2[i]/Sum_w2) );
		
		alpha += omega*omega*(i+1.0);
	}
	
//	cerr <<"alpha is "<<alpha<<"\n";
	
	beta=0;
	
	
	for(i=1; i < this->n ; i++){
		for( j=i; j< this->n ; j++){
		
			double omega_i = ( (w1[i-1]/Sum_w1) - (w2[i-1]/Sum_w2) );
			double omega_j = ( (w1[j-1]/Sum_w1) - (w2[j-1]/Sum_w2) );
		
		
			if( i == j)
				beta  +=  i*i* omega_i*omega_i * sigma_ii_FU1995( i  );
			else
				beta  +=  2*i*j* omega_i*omega_j * sigma_ij_FU1995( i, j  );
						
		}
	}
//	cerr <<"beta is "<<beta<<"\n";

	

	return alpha*theta + beta*theta*theta ;

}



/*
	Compute the Variance of an Omega or an omega vector
*/
double coalescent_tree_diversity::compute_CoVarw1w2( double *w1, double *w2, double theta){

	double Sum_w1=0.0;
	double Sum_w2=0.0;

	int i,j;

	double cov=0;
	
	
	/*
		Compute the alpha and beta of the variance
	*/

	for(i=0;i<n-1;i++){
		 Sum_w1 += w1[i];
		 Sum_w2 += w2[i];
	}

	if(Sum_w1 == 0 || Sum_w2 == 0){
		cerr << "cannot compute covariance, one vector sums to 0";
		return -1;
	}

	
	
	for(i=1; i < this->n ; i++){
		for( j=1; j< this->n ; j++){
		
			double omega_1 = (w1[i-1]/Sum_w1);
			double omega_2 = (w2[j-1]/Sum_w2);
		
		
			if( i == j )
				cov  +=  i*i* omega_1*omega_2 * sigma_ii_FU1995( i  );
			else
				cov  +=  i*j* omega_1*omega_2 * sigma_ij_FU1995( i, j  );
						
		}
	}
//	cerr <<"cov is "<<cov<<"\n";

	return cov * theta * theta ;

}




/*
	Compute the alpha and beta of an Omega or an omega vector
*/
void coalescent_tree_diversity::compute_VarTAlphaBeta( double *w1, double *w2, double *alpha, double *beta){

	double Sum_w1=0.0;
	double Sum_w2=0.0;

	int i,j;

	int tag=0;

	
	/*
		Compute the alpha and beta of the variance
	*/

	for(i=0;i<n-1;i++){
		 Sum_w1 += w1[i];
		 Sum_w2 += w2[i];
		 if(Sum_w2 != 0)tag=1;
	}

	if(tag == 0)   /* the second vector is null, so report alpha/beta of the first vector (ie. a theta estimator ) */
		Sum_w2=1;


	*alpha=0;
	for(i=0; i < this->n - 1 ; i++){
	
		double omega = ( (w1[i]/Sum_w1) - (w2[i]/Sum_w2) );
		
		*alpha += omega*omega*(i+1.0);
	}
	
	
	
//	cerr <<"alpha is "<<alpha<<"\n";
	
	*beta=0;
	
	for(i=1; i < this->n ; i++){
		for( j=i; j< this->n ; j++){
		
			double omega_i = ( (w1[i-1]/Sum_w1) - (w2[i-1]/Sum_w2) );
			double omega_j = ( (w1[j-1]/Sum_w1) - (w2[j-1]/Sum_w2) );
		
		
			if( i == j)
				*beta  +=  i*i* omega_i*omega_i * sigma_ii_FU1995( i  );
			else
				*beta  +=  2*i*j* omega_i*omega_j * sigma_ij_FU1995( i, j  );
						
		}
	}
//	cerr <<"beta is "<<beta<<"\n";


	return;
}





/*
	print out the tree using the number of Sites instead of the branch length
	adpated from my previous C code.
	to print the whole tree, the first pointer has to be at the top of the tree !!!
*/
static void treeout_FromMut( locus *ptree ){

	locus *ptmp;

	if( ptree == NULL )return;                                       // stop that recursive function

	if(ptree->get_descendant1() != NULL){                            // it is not an final node
		cout << "(";
		ptree = ptree->get_descendant1();
	} 
	else{                                                            // it is a final node

		printf("n[%d]:%d", ptree->get_id(),  ptree->get_newSites() );
		ptmp = ptree;
		
		while(ptmp->get_ancestor()->get_descendant2() == ptmp){   // are we a descendant2 ??
			ptmp=ptmp->get_ancestor();
	
			if( ! ptmp->is_coalesced() ){                          // we are at the the top of the tree
				cout << ");\n";
				ptree=NULL;
				return;
			}                                                     // in tree reconstruction, we want the branch length
			
			printf("):%d", ptmp->get_newSites()  );
		}
		
		ptree=ptmp->get_ancestor()->get_descendant2();          // switch to desc2
		cout << ",";
	}
	
	treeout_FromMut( ptree );
}
void coalescent_tree_diversity::printout_treeFromMut( void  ){
	if(this->n<2)
		return;
	treeout_FromMut( get_all_samples()[0]->oldest_sample()->get_loci()[0] );
}



/*
	print out the tree
	adpated from my previous C code.
	to print the whole tree, the first pointer has to be at the top of the tree !!!
*/
static void treeout_WithMut( locus *ptree ){

	locus *ptmp;

	if( ptree == NULL )return;                                       // stop that recursive function

	if(ptree->get_descendant1() != NULL){                            // it is not an final node
		cout << "(";
		ptree = ptree->get_descendant1();
	} 
	else{                                                            // it is a final node

		printf("n[%d]_m%d:%.10f", ptree->get_id(),  ptree->get_newSites(), ptree->get_ancestor()->get_time());
		ptmp = ptree;
		
		while(ptmp->get_ancestor()->get_descendant2() == ptmp){   // are we a descendant2 ??
			ptmp=ptmp->get_ancestor();
	
			if( ! ptmp->is_coalesced() ){                          // we are at the the top of the tree
				cout << ");\n";
				ptree=NULL;
				return;
			}                                                     // in tree reconstruction, we want the branch length
			
			printf(")m%d:%.10f", ptmp->get_newSites() , ptmp->get_ancestor()->get_time() - ptmp->get_time());  
		}
		
		ptree=ptmp->get_ancestor()->get_descendant2();          // switch to desc2
		cout << ",";
	}
	
	treeout_WithMut( ptree );
}

void coalescent_tree_diversity::printout_treeWithMut( void  ){
	if(this->n<2)
		return;
	treeout_WithMut( get_all_samples()[0]->oldest_sample()->get_loci()[0] );
}

void coalescent_tree_diversity::compute_Haplos( void  ){

	int **mat=NULL;
	int i,j;
	
	mat = compute_Kmatrix();
	
	this->Haplos = n;
	
	for(i=0;i<this->n;i++)
		for(j=0;j<i;j++)
			if( mat[i][j] == 0){
				this->Haplos--;
				break;
			}
	
	for(i=0 ; i<this->n ; i++ )
		delete mat[i];
	delete mat;
	
}
