/***
		file:     coalescent_tree.cpp
		function: this would built several layers of sample objects that is the tree
		author:   <amikezor>
		date:     sept 06
		modif:    oct 08, nov 08
		          feb 2009 - change sweep so that sampling can be done when sweep is ongoing
			  sept 09 - handle n=1 case (no tree)
***/

#include "population.h"
#include "sample.h"
#include "coalescent_tree.h"
#include "locus.h"
#include "random.h"
#include "misc.h"

#include <iostream>
using namespace std;

#include <stdio.h>
#include <string.h>

coalescent_tree::coalescent_tree( int n, population *p ){

	pop=p;
	this->n=n;

	coalescent_samples = new sample * [n];                              // there will be n-1 coalescent events, so n layers (sample)
	
	coalescent_samples[0] = new sample( p );
	coalescent_samples[0]->set_loci( n, n );                             // an array with n loci and get memory for n
	
	for(int a=1;a<n;a++){
		
		coalescent_samples[a] = new sample( p );
		coalescent_samples[a]->set_loci(n-a, 1);                          // on each sample, one less loci, only one new locus
		coalescent_samples[a-1]->connect_ancestor( coalescent_samples[a] );
	}
	
}

coalescent_tree::~coalescent_tree(  ){

	delete coalescent_samples[0];                                       // this would free also ancestral samples (in the destructor of sample)
	delete[] coalescent_samples;                                        // free the pointer array


}


locus * coalescent_tree::get_mrca_2loci( locus *l1, locus *l2 ){


	int ncommon=0;               // how much they have in common. Shoud be >=1 (1: only the root)
	int i;                       // counter


	if(! l1->get_n_idpath() )
		l1->compute_idpath( );
	
	if( !l2->get_n_idpath() )
		l2->compute_idpath( );
	

	while( l1->get_idpath()[ l1->get_n_idpath()-1-ncommon ]  ==  l2->get_idpath()[ l2->get_n_idpath()-1-ncommon ] )
		ncommon++;
		
	locus * ploc=l1;
	
	for( i=0; i<l1->get_n_idpath()-ncommon ; i++ )
		ploc = ploc->get_ancestor();

	
	return ploc;
}


locus *coalescent_tree::get_mrca_nloci( locus **list_locus, int size_list ){


	int i;
	
	locus *plocus = list_locus[0];
	
	for(i=1;i<size_list;i++)
		plocus = get_mrca_2loci( plocus, list_locus[i] );
	
	
	return plocus;
}




/*
	check all coalescent events from the oldest until you find the species you are looking for.
*/

locus *coalescent_tree::get_mrca_species( int sp ){

	if( sp == -1 )
		return this->coalescent_samples[n-1]->get_loci()[0];
	
	locus *mrca;
	int nloc = this->get_current_sample()->get_nloci_species( sp  );
	
	locus **myloc = this->get_current_sample()->get_loci_species( sp  );

	mrca = this->get_mrca_nloci( myloc, nloc );

	delete myloc;
	
	return mrca;
	

}







/*
	Coalescent processes of all ancestral samples
*/

void coalescent_tree::standard_coalescent( int nsamples, int *sample_sizes, double *serial_times ){
	
	sample *psample = this->coalescent_samples[0];
	
	if( nsamples >1 )
	{
		int c=0;
		for(int s=0; s<nsamples; s++)
			for(int i=0;i<sample_sizes[s];i++)
				psample->get_loci()[ c++ ]->set_time(serial_times[s]);
	}


	while( psample->get_nloci() > 1 ){
	
		psample->coalesce_sample();
		psample = psample->get_ancestor();
	}
	
}


void coalescent_tree::fixedtimes_coalescent( double *times ){

	int tt=0;

	sample *psample = this->coalescent_samples[0];
	

	while( psample->get_nloci() > 1 ){
	
		psample->coalesce_sample_fixedtime( times[tt] );
		psample = psample->get_ancestor();
		tt++;
	}
	

}



/*
	Coalescent processes of all ancestral samples
*/

void coalescent_tree::bs_coalescent( void ){
	
	sample *psample = this->coalescent_samples[0];
	
	while( psample != NULL ){
		psample = psample->bs_coalesce_sample();
	}
	
	
}


/*
	Coalescent processes of all ancestral samples
*/

void coalescent_tree::bottleneck_coalescent( void ){
	

	/*
		Simulate a std coalescent
	*/
	this->standard_coalescent();
	
	/*
		Compact the tree to take the bottleneck into account
	*/
	this->bottleneck_reduce_tree();
	
}


void coalescent_tree::bottleneck_reduce_tree( void ){

	sample *psample = this->coalescent_samples[0];
	
	do{
		psample->bottleneck_reduce_time_sample();
		psample = psample->get_ancestor();

	}while( psample != NULL );
	


}

/*
	Coalescent processes of all ancestral samples
	type should be 'l'inear or 'e'xponential
*/

void coalescent_tree::growth_coalescent( char type ){
	

	/*
		Simulate a std coalescent
	*/
	this->standard_coalescent();
	
	/*
		Compact the tree to take the exponential/linear growth into account
	*/
	
	sample *psample = this->coalescent_samples[1];
	
	do{
		psample->growth_reduce_time_sample( type );
		psample = psample->get_ancestor();

	}while( psample != NULL );
	
}




/*
	Coalescent processes for isolation/speciation tree with ni[x] the sample x. there are pop->get_isolation_ns() species
*/

void coalescent_tree::isolation_coalescent( int *ni ){
	

	sample *psample = this->coalescent_samples[0];

	int sum_ni=0;
	for(int x=0;x<pop->get_isolation_ns(); x++)
		sum_ni+=ni[x];
	
		
	if( sum_ni != psample->get_nloci() ){
		cerr <<"coalescent_tree::speciation_coalescent the sum of all species loci ( "<< sum_ni <<" ) should be equal to "<<psample->get_nloci()<<", bye\n";
		exit(1);
	}
	
	sum_ni=0;
	for( int i = 0; i < pop->get_isolation_ns(); i++){                          // for the today sample, samples[0]
	
		//cout << ni[i]<<"\n";
	
		for( int j=0; j<ni[i]; j++){
			
			psample->get_loci()[sum_ni]->set_species(  i  );            // set the first ni[0] locis to species 0, the ni[1] next one to sp 1, ...
			psample->get_loci()[sum_ni]->reset_connections();
			
			sum_ni++;
		}
	}
//	cerr <<"sp all sets\n";
	
	while( psample->get_nloci() > 1 ){

//		cerr << "nloci is "<<psample->get_nloci()<<"\n";

		psample->isolation_coalesce_sample();
		psample = psample->get_ancestor();

	}

}





/*
	Coalescent processes for strutured-population tree with ni[x] the sample x.
	There are pop->get_isolation_ns() species
*/

void coalescent_tree::structured_coalescent( int nsamples,  int **samples, double *serial_times ){
	

	sample *psample = this->coalescent_samples[0];
	
	int *sample_sizes = new int [ nsamples ];
	int total_n = 0;
	
	if( sample_sizes == NULL){
		cerr << "coalescent_tree::structured_coalescent cannot allocate sample_sizes, bye\n";
		exit(3);
	}
	
	

	for(int y=0; y<nsamples; y++){
	
		sample_sizes[y] = 0;
		for( int x=0 ; x<pop->get_isolation_ns() ; x++ ){
			total_n += samples[y][x];
			sample_sizes[y] += samples[y][x];
		}
	}
		
	if( total_n != psample->get_nloci() ){
		cerr <<"coalescent_tree::structured_coalescent the sum of all species loci ( "<< total_n <<" ) should be equal to "<<psample->get_nloci()<<", bye\n";
		exit(1);
	}
	
	total_n=0;
	for (int s=0; s <  nsamples; s++){                                               // for all samples
	
		for( int i = 0; i < pop->get_isolation_ns(); i++){                          // for the today sample, samples[0]
	
			for( int j=0; j<samples[s][i]; j++){
			
				psample->get_loci()[total_n]->set_species(  i  );            // set the first ni[0] locis to species 0, the ni[1] next one to sp 1, ...
				psample->get_loci()[total_n]->reset_connections();
			
				total_n++;
			}
		}
	}
	
	if( nsamples >1 )
	{
		int c=0;
		for(int s=0; s<nsamples; s++)
			for(int i=0;i<sample_sizes[s];i++)
				psample->get_loci()[ c++ ]->set_time(serial_times[s]);
	}



	while( psample->get_nloci() > 1 ){

		//cerr << "nloci is "<<psample->get_nloci()<< "; current nloci are " <<  psample->get_nloci_befT( psample->get_time() ) << "\n";

		psample->structured_coalesce_sample();
		psample = psample->get_ancestor();

	}

	delete sample_sizes;

}




/*
	Coalescent processes for a neutral locus linked to a sweeping locus
*/

void coalescent_tree::sweep_coalescent(  ){


	sample *psample = this->coalescent_samples[0];
	
	extern Random R;             //  random assignment of 
	
	
	for( int i = 0; i<psample->get_nloci() ; i++){               // for the latest sample, samples[0]
	
		int random_sp;


		if( R.uniform_dev() <= pop->get_sweep_p() )
			random_sp=0;
		else
			random_sp=1;

		psample->get_loci()[i]->set_species( random_sp );    // set the species to 0 if linked to the selected locus, 1 otherwise
		psample->get_loci()[i]->reset_connections();         // set to null the connections
		
	}
	
	pop->set_sweep_over( false );
	
	/*
		Do the coalescent
	*/
	while( psample->get_nloci() > 1 ){
	
		psample->sweep_coalesce_sample(  );
		psample = psample->get_ancestor();
	}
	
		
	set_nsurvivors( get_oldest_sample()->get_nspecial() );
	
	//
	// if did not even reach the sweep, nsurvivors==0
	// if died in the sweep, it is set to 1
	// otherwise it is > 1
	
	if (nsurvivors < 0){
		cerr << "nsurvivors set to -1 !!\n";   
		exit(1);
	}

}


/*
	Coalescent processes of all ancestral samples
	with one subtree of size fixed to k
*/

void coalescent_tree::nested_coalescent( int k ){
	
	sample *psample = this->coalescent_samples[0];
	
	
	for( int i=0 ; i< psample->get_nloci() ; i++)
		if(i<k)
			psample->get_loci()[i]->set_species( 1 );
		else
			psample->get_loci()[i]->set_species( 2 );
	
	
	while( psample->get_nloci() > 1 ){
	
		psample->nested_coalesce_sample(  );
		
		psample = psample->get_ancestor();
	}
	
}




/*
	Mutation Funcstions
*/


void coalescent_tree::mutate_whole_tree_seq(  int L, int errors ){
	coalescent_samples[n-1]->get_loci()[0]->mutate_subtree_seq( pop->get_Theta(), L, errors ); 
}


void coalescent_tree::mutate_whole_tree( int S, double mean_err ){

	coalescent_samples[n-1]->get_loci()[0]->mutate_subtree( pop->get_Theta(), S, mean_err ); 
}

void coalescent_tree::clean_whole_tree( void ){
	coalescent_samples[n-1]->get_loci()[0]->clean_subtree( ); 
}



void coalescent_tree::mutate_whole_tree_TE( int S, double mean_err ){

	int oldest_te = get_all_samples()[0]->oldest_sample()->get_loci()[0]->get_size_array_TE() -1;

	get_all_samples()[0]->oldest_sample()->get_loci()[0]->get_array_TE()[   oldest_te ].mutate_subtree( pop->get_Theta(), S, mean_err ); 
}


/*
	Statistics on the tree
*/

double coalescent_tree::get_Ttot_singletons(void){

	double Ttot_s = 0.0;

	for(int x=0; x<n; x++){
		Ttot_s += coalescent_samples[0]->get_loci()[x]->get_ancestor()->get_time();
	}

	return Ttot_s;
}



void coalescent_tree::printout_nsites( void ){
	cout << get_all_samples()[0]->get_loci()[0]->get_size_sequence() << " polymorphic sites\n" ;
}






/*
	All functions for OUTPUT
*/




/*
	print out the tree
	adpated from my previous C code.
	to print the whole tree, the first pointer has to be at the top of the tree !!!
*/
static void treeout( locus *ptree ){

	locus *ptmp;

	if( ptree == NULL )return;                                   // stop that recursive function

	if(ptree->get_descendant1() != NULL){                        // it is not an final node
		cout << "(";
		ptree = ptree->get_descendant1();
	} 
	else{                                                        // it is a final node
		
		if( ptree->IsTE() == 1 )
			printf("n[%d]_l[%d]:%.10f", ptree->get_id(), ptree->get_locus(), ptree->get_ancestor()->get_time() - ptree->get_time());
		else
			printf("n[%d]:%.10f", ptree->get_id(), ptree->get_ancestor()->get_time() - ptree->get_time());
		
		
		ptmp = ptree;
		
		while(ptmp->get_ancestor()->get_descendant2() == ptmp){   // are we a descendant2 ??
			ptmp=ptmp->get_ancestor();
	
			if( ! ptmp->is_coalesced() ){                          // we are at the the top of the tree
				cout << ");\n";
				ptree=NULL;
				return;
			}                                                     // in tree reconstruction, we want the branch length
			
			printf("):%.10f",  ptmp->get_ancestor()->get_time() - ptmp->get_time());  
		}
		
		ptree=ptmp->get_ancestor()->get_descendant2();          // switch to desc2
		cout << ",";
	}
	
	treeout( ptree );
}



void coalescent_tree::printout_tree( void  ){

	if(this->n < 2)
		return;

	treeout( get_all_samples()[0]->oldest_sample()->get_loci()[0] );
}






/*
	print out the sequences from the current sample
*/
void coalescent_tree::printout_sequences( void ){

	if(  get_all_samples()[0]->get_loci()[0]->get_sequence() == NULL ){
		cout <<"No generated sequences (is S=0 ?)\n";
		return;
	}

	for(int i=0;i<n;i++){
		cout <<">seq["<<get_all_samples()[0]->get_loci()[i]->get_id()<<"]\n"<<get_all_samples()[0]->get_loci()[i]->get_sequence()<<"\n";   
	}
	
}


/*
	print out the sequences from the current sample
*/
char **coalescent_tree::extract_sequences( void ){


	int len;
	char **my_sequences=NULL;
	
	len = strlen( get_all_samples()[0]->get_loci()[0]->get_sequence() );
	
	if(len > 0){
		
		my_sequences = (char **)malloc( n*sizeof(char *) );
		for(int i=0;i<n;i++){
			my_sequences[i] = (char *)malloc( (len+1)*sizeof(char) );
			strcpy( my_sequences[i], get_all_samples()[0]->get_loci()[i]->get_sequence() );
		}
		
	}


	return my_sequences;
	
}

/*
	print out the sequences from the current sample
*/
char **coalescent_tree::merge_sequences( char **input_sequences ){

	if( input_sequences == NULL)
	{
		input_sequences = extract_sequences(  );
	}
	else
	{

		int l1 = strlen( input_sequences[0] );
		int l2 = (get_all_samples()[0]->get_loci()[0]->get_sequence())?strlen( get_all_samples()[0]->get_loci()[0]->get_sequence() ):0;
				
		for(int i=0;i<n;i++){
			input_sequences[i] = (char *)realloc( input_sequences[i], (l1+l2+2)*sizeof(char) );
			
			input_sequences[i][l1] = '-';
			if(l2)
				strcpy( input_sequences[i]+l1+1, get_all_samples()[0]->get_loci()[i]->get_sequence() );
		}

	
	
	}
	
	return input_sequences;
		

}

/*
	print out the sequences from this node to the leaves sequences
*/
void coalescent_tree::printout_all_sequences( void ){

	if(  get_all_samples()[0]->get_loci()[0]->get_sequence() == NULL ){
		cout <<"No polymorphic sites\n";
		return;
	}

	printout_sequences(  );              // sequences from all loci at level 0
	
	for(int i=1;i<n;i++){
		cout <<"seq["<<get_all_samples()[i]->get_loci()[0]->get_id()<<"]\t"<<get_all_samples()[i]->get_loci()[0]->get_sequence()<<"\n";   
	}
	
}


/*
 ***

	Transposable Elements
***
*/

/*
	This function computes TE tree using information from ancestral sample
	EXCEPT for sample mrca for which no ancestor exist
*/
void coalescent_tree::coal_TE_TopDown_FromMrca( float dupl_rate ){

	sample *psample = get_all_samples()[0]->oldest_sample();


	// built TE tree
	
	while ( psample->get_descendant() != NULL ){
	
		psample = psample->get_descendant();
		psample->coalesce_TE_sample_TopDown( dupl_rate );
	
	}

	// For the most recent sample, compute Ttot and nleaves for all loci.
	
	for( int i=0 ; i < psample->get_nloci() ; i++ ){
		
		for ( int j=0 ; j  < psample->get_loci()[i]->get_size_array_TE() ; j++ ){
			
			psample->get_loci()[i]->get_array_TE()[j].compute_nleaves();
			psample->get_loci()[i]->get_array_TE()[j].compute_Ttime();
		}
		
	}
	
	
	// For the other samples, compute Ttot and nleaves for loci[0].

	while(  psample->get_ancestor() != NULL  ){
		
		psample = psample->get_ancestor();

		for ( int j=0 ; j  < psample->get_loci()[0]->get_size_array_TE() ; j++ ){
		
			psample->get_loci()[0]->get_array_TE()[j].compute_nleaves();
			psample->get_loci()[0]->get_array_TE()[j].compute_Ttime();
		}
		
	}
}



/*
	This function sets a chosen number of TE in the MRCA locus (of the MRCA sample)
*/
void coalescent_tree::coal_TE_TopDown_mrca( int n_TE_mrca, float dupl_rate){

	sample *psample = get_all_samples()[0]->oldest_sample();

	// Top sample, generate the tree
	//


	this->TE_anc_locus = psample->get_loci()[0];

	psample->coalesce_TE_mrca_BottomUp( dupl_rate, n_TE_mrca, TE_anc_locus );
	
	this->coal_TE_TopDown_FromMrca( dupl_rate );

}



/*
	Using the Time of the first TE, this function chooses randomly a locus 
	from the appropriate sample to be in the lineage of the genome with the first TE.
	
*/
void coalescent_tree::coal_TE_TopDown_1stTE( float time_1st_TE, float dupl_rate){

	sample *psample = get_all_samples()[0]->oldest_sample();
	
	
	extern Random R;
	
	
	while( time_1st_TE < psample->get_time() ){
	
		cerr <<"first TE is younger than sample "<< psample->get_id() <<" (time "<< psample->get_time() <<" )\n";

		psample = psample->get_descendant();
		
	
	}

	this->TE_anc_locus = psample->get_loci()[ R.uniform_int(0, psample->get_nloci()-1) ];      // pick a random locus

	cerr <<"the selected locus is "<<TE_anc_locus->get_id()<<"\n";

	psample->coalesce_1st_TE_locus_TopDown( dupl_rate, time_1st_TE, TE_anc_locus );	      // set the first TE in its lineage

	this->coal_TE_TopDown_FromMrca( dupl_rate );

}
void coalescent_tree::printout_tree_TE( void  ){

	if(this->n < 2)
		return;

	int nTE_anc = TE_anc_locus->get_n_TE();
	
	if( (TE_anc_locus->get_array_TE() + 2*nTE_anc -2)->get_descendant1() == NULL ){
		cout <<"(te[0]_l["<<TE_anc_locus->get_id()<<"]:0.0)\n";
		return;
	}
	
//	cerr << "the root locus is " <<  get_all_samples()[0]->oldest_sample()->get_loci()[0]->get_array_TE()[2*nTE_mrca -2].get_id(  ) << "\n";

	treeout( (locus *) ( TE_anc_locus->get_array_TE() + 2*nTE_anc -2  ));
}
