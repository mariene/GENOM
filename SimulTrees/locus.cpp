/***
		file:     locus.cpp
		function: loci functions
		author:   <amikezor>
		date:     feb 04
		modif:    mar 09 -- add the id_path
		modif:    jan 2012 -- add JukesCantor
***/


#include <iostream>
using namespace std;

#include <string.h>
#include <stdio.h>


#include "locus_TE.h"
#include "locus.h"
#include "random.h"

int locus::Count_All_Loci=0;

/*
	Constructor
*/
locus::locus( double t ){
	
	id = Count_All_Loci++;
	time=t;                                         // initiate the locus at time t
	ancestor=descendant1=descendant2=NULL;          // all genealogy is still free
	Ttime=0.0;                                        // still unset

	sequence=NULL;                                  // sequence is set to unused;
	size_sequence=0;
	size_seq_mem=0;
	
	newSites=0;
	nleaves=1;
	nleaves_sp=0;
	species=0;

	idpath=NULL;       
	n_idpath=0;        

	array_TE=NULL;
	size_array_TE=0;
	n_TE=0;


	//cerr << "<-- create locus "<<this->id<<" with  time "<<time<<"\n";
}


/*
	Destructor
*/
locus::~locus(  ){
	
	//cerr <<"--> delete loci "<<id<<"\n";

	if(descendant1 != NULL)
		descendant1->ancestor=NULL;               // if it has any descendants set their anesctor as free
	
	if(descendant2 != NULL)
		descendant2->ancestor=NULL;
	
	if(sequence != NULL)                         // if there is a sequence delete it;
		delete sequence;
	
	if( idpath != NULL)
		delete idpath;
	
	if( array_TE != NULL )
		delete[] array_TE;
	
		
	Count_All_Loci--;
}



/*
	Connect the current locus to another, its ancestor
*/
void locus::connect_ancestor( locus *anc, char descendant ){


//	cerr <<"connect anc."<<anc->get_id()<<" to desc"<< descendant <<" "<<this->id<<"\n";

	if( descendant == 1 )
		anc->descendant1 = this;                         // set the descendant of the ancestor
	else if( descendant == 2 )
		anc->descendant2 = this;
	else{
		cerr <<"locus::connect_ancestor, a descendant is either 1 or 2\n";
		return;
	}


	this->ancestor = anc ;                              // set the ancestor

	if(anc->descendant1 && anc->descendant2 )
		anc->nleaves = anc->descendant1->nleaves + anc->descendant2->nleaves;

}



/*
	Ask a locus if it is coalesced with another one
*/

bool locus::is_coalesced( void ){
	if(this->ancestor)
		return true;
	else
		return false;
}




/*
	Mutate all tree beyond this node according to Poisson( Ttime * theta/2 ); 
	Ttime is in Ne generation, where theta mutations occurs (on average)
	mutation are stored into the newSites values
	if S is defined, it is fixed otherwise it is drawn as a poissondev
*/
void locus::mutate_subtree( double theta , int S, double mean_errors ){

	extern Random R;
	
	double Tr=0.0;                  // the random time at whihc a mutation occured
	double uni_r=0.0;               // a number between ]0-1]
	locus *node=NULL;               // the node picked for being mutated
	int i=0;                        // a counter
	int nmut=0;                     // how man mutations in this tree ?
	
	if( !S )
		nmut = R.poisson_dev( this->Ttime * theta / 2.00  );        // get the number of mutation event;
	else
		nmut=S;
		
	//cout << nmut<<" mutations on subtree with Ttot "<< this->Ttime << "\n";

	if(!nmut)return;                                            // if no mutation, just stay calm...

	/*
		First take care of real mutations
	*/
	for(i=0;i<nmut;i++){
	
		uni_r = R.uniform_dev();                 // get a uniform number between ]0,1[;

		Tr = Ttime * uni_r;                      // change a ]0,1] number into a ]0,Ttime[;

		node = this->pick_node( Tr );            // pick_node according to this branch length;
		                                         // by dichotomic search top down
		
		//cerr <<"pick node "<<node->get_id()<<" from Tr of "<<Tr<<"\n";
				
		node->newSites++;                      // mutate this node;
		
	}
	
	/*
		Then add PCR/sequencing/cloning mutations as external singletons
	*/
	
	if( mean_errors ){
		

		int rand_errors = R.poisson_dev( mean_errors );
		
		//cout << "rand_errors is "<<rand_errors<<"\n";
		
		for(i=nmut;i<nmut+rand_errors;i++){
	
			node=this;
		
			while(node->descendant1 != NULL){              // go to any leaf randomnly
				uni_r = R.uniform_dev();
				if(uni_r>0.5)
					node = node->descendant2;
				else
					node = node->descendant1;
			}

			node->newSites++; 	                      // mutate it
		}
	}

}




void locus::clean_subtree( void ){

	newSites=0;
	species=0;
	nleaves = 1;
	nleaves_sp = 0;
	
	ancestor=NULL;
	
	if(this->idpath != NULL){
		delete this->idpath;
		this->idpath=NULL;
		this->n_idpath=0;
	}
	if(this->sequence != NULL){
		delete this->sequence;
		this->sequence=NULL;
		this->size_seq_mem=0;
	}

	
	if( descendant1 )
		descendant1->clean_subtree();
	if( descendant2 )
		descendant2->clean_subtree();
		
		
	descendant1 = NULL;
	descendant2 = NULL;
}


static char JukesCantor(  char input ){

	extern Random R;
	char dna[4]={'A','C','G', 'T'};
	int i;
	
	switch(input){
		case 'A':
			i=0;
			break;
		case 'C':
			i=1;
			break;
		case 'G':
			i=2;
			break;
		case 'T':
			i=3;
			break;
		default:
			fprintf(stderr, "JukesCantor: unknown input base, bye"),exit(1);
	}	
	
	char n = dna[ (i+R.uniform_int(1,3))%4 ]; 
	
	return n;
	
}


/*
	Init a sequence of size 'len' and set it all to 'base'
	init all sequence in the desecdant
*/
void locus::init_sequence( int len, char base ){


	if(size_seq_mem < len)
	{
			
		if(sequence != NULL)
			delete[] sequence;
		
		if( (sequence=new char[len+1] ) == NULL)
			cerr << "locus::init_sequence, cannot allocate sequence, bye\n",exit(3);

		size_seq_mem = len;

	}
	
	memset( (void *)sequence, (int)base, (size_t)len);
	
	sequence[len]=0;
	size_sequence=len;
	

	if( descendant1 != NULL) descendant1->init_sequence( len, base );
	if( descendant2 != NULL) descendant2->init_sequence( len, base );
	
}



/*
	mutate the base m into 'value'
	and iterate to all descendant 
*/
void locus::mutate_seq( int n, char value ){

	if(sequence == NULL || size_sequence<n ){
		cerr << "locus::mutate_seq, cannot mutate the sequence. It is too short or it does not exist !\n";
		return;
	}

	while(sequence[n] == value){
		//cerr << "locus::mutate_seq, the mutation will not change the sequence !!\n";
		//printf("seq of node %d was %c and will be %c\n", id, sequence[n], value);
		value = JukesCantor(value);
	}
		
	sequence[n]=value;                                     // mutate
	
	if(descendant1 != NULL)descendant1->mutate_seq(n,value);   // spread the mutation
	if(descendant2 != NULL)descendant2->mutate_seq(n,value);

	return;
}


/*
	Mutate all tree beyond this node according to Poisson( Ttime * theta/2 ); 
	Ttime is in Ne generation, where theta mutations occurs (on average)
	Mutations are stored in sequences
*/
void locus::mutate_subtree_seq( double theta , int L, int errors ){

	extern Random R;
	
	double Tr=0.0;                  // the random time at whihc a mutation occured
	double uni_r=0.0;               // a number between ]0-1]
	locus *node=NULL;               // the node picked for being mutated
	int i=0;                        // a counter
	int nmut=0;                     // how man mutations in this tree ?
	
	char mutant;
		
	nmut = R.poisson_dev( this->Ttime * theta / 2.00  );        // get the number of mutation event;
	//cout << nmut<<" mutations on subtree with Ttot " << this->Ttime <<" \n";

	if(!nmut)return;                            // if no mutation, just stay calm...

	if(L == 0)
		this->init_sequence( nmut+errors, 'A' );    // init this locus and all it descendant to this string
	else{
		this->init_sequence( L, 'A' );
	}

	/*
		First take care of real mutations
	*/
	for(i=0;i<nmut;i++){
	
		uni_r = R.uniform_dev();                 // get a uniform number between ]0-1[;

		Tr = Ttime * uni_r;                      // turn thta 0-1 number into a ]0-Ttime[;

		node = this->pick_node( Tr );            // pick_node according to this branch length;
		
		//cerr <<"pick node "<<node->get_id()<<" from Tr of "<<Tr<<"\n";
		
		if(L == 0){
			
			mutant = JukesCantor( node->get_sequence()[i] );
			node->mutate_seq( i, mutant );                  // mutate this node;
		}
		else
		{
			int x=R.uniform_int(0, L-1);
			mutant = JukesCantor( node->get_sequence()[x] );
			
			//printf("at node %d site %d change %c to %c\n", x,  id,  node->get_sequence()[x], mutant);
			
			node->mutate_seq( x, mutant );                  // mutate this node;
		}
		
	}
	
	/*
		Then add PCR/sequencing/cloning mutations as singletons
	*/
	for(i=nmut;i<nmut+errors;i++){
	
		node=this;
		
		while(node->descendant1 != NULL){         // go to any tip
			uni_r = R.uniform_dev();
			if(uni_r>0.5)
				node = node->descendant2;
			else
				node = node->descendant1;
		}

		if(L == 0){
			mutant = JukesCantor( node->get_sequence()[i] );
			node->mutate_seq( i, mutant );                  // mutate this node;
		}
		else
		{
			int x=R.uniform_int(0, L-1);
			mutant = JukesCantor( node->get_sequence()[x] );
			node->mutate_seq( x, mutant );                  // mutate this node;
		}
	}
	

}


/*
	Compute the array which all nodes id drom this one to the root
	this is stored in an int array
	Very usefull when one wants to find a common ancestror (= find the common path)
*/
void locus::compute_idpath(void){

	locus *ploc;
	int i;

	ploc=this;
	this->n_idpath=1;
	
	while( (ploc = ploc->ancestor) != NULL)
		this->n_idpath++;
	
	this->idpath = new int[n_idpath];
	
	i=0;
	ploc=this;
	*(this->idpath)=this->id;

	while(  (ploc = ploc->ancestor)  && ++i ){
		*(this->idpath + i)=ploc->id;
	}
}





/**

	Private Function

**/

#define EPSILON 1e-4

locus *locus::pick_node( double Tr ){

	if( Tr > Ttime+EPSILON )
		fprintf(stderr,"locus::pick_node, pick Tr>Ttime (%.10f/%.10f), it is unexpected...\n",Tr,Ttime);
	if( Tr < 0 )
		cerr << "locus::pick_node, cannot pick Tr<0 ("<<Tr<<"), it is unexpected...\n";


	if(descendant1  == NULL || descendant2 == NULL)    // it should not happen but may due to
		return this;                               // the imprecision of floating numbers.

	double t1 = time - descendant1->get_time();        // time of the branch leading to desc 1
	double t2 = time - descendant2->get_time();        // time of the branch leading to desc 2
	double Tt1 = descendant1->get_Ttime();             // Total time under desc 1

	if( Tr < Tt1 )
		return descendant1->pick_node( Tr );
	else if( Tr < t1+Tt1 )
		return  descendant1;
	else if( Tr < t1+Tt1+t2 )
		return descendant2;
	else
		return descendant2->pick_node( Tr-(t1+Tt1+t2) );

}

#undef EPSILON






void locus::compute_Ttime( void ){

	this->Ttime=0;
	
	if( this->descendant1 != NULL )
		this->Ttime += this->get_time() - this->descendant1->get_time() + this->descendant1->get_Ttime();

	if( this->descendant2 != NULL )
		this->Ttime += this->get_time() - this->descendant2->get_time() + this->descendant2->get_Ttime();


}

void locus::compute_nleaves( void ){

	this->nleaves=0;
	
	if( this->descendant1 != NULL ){
		this->descendant1->compute_nleaves();
		this->nleaves +=  this->descendant1->nleaves;
	}

	if( this->descendant2 != NULL ){
		this->descendant2->compute_nleaves();
		this->nleaves +=  this->descendant2->get_nleaves();
	}
	
	if( this->descendant2 == NULL )
		this->nleaves = 1;

}

void locus::compute_nleaves_sp( int sp ){

	if( sp == -1 )
		this->compute_nleaves();

	this->nleaves_sp=0;
	
	if( this->descendant1 != NULL )
	{
		this->descendant1->compute_nleaves_sp( sp );
		this->nleaves_sp += this->descendant1->get_nleaves_sp( sp );
	}


	if( this->descendant2 != NULL )
	{
		this->descendant2->compute_nleaves_sp( sp );
		this->nleaves_sp +=  this->descendant2->get_nleaves_sp( sp );
	}


	if( this->descendant1 == NULL && this->get_species() ==  sp )
		this->nleaves_sp=1;
	
}

int locus::get_nleaves_sp( int sp ){


	if( sp == -1 )
		return this->get_nleaves();
	
	return  this->nleaves_sp;

}
