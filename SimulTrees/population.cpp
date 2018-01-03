/***
		file:     popoulation.cpp
		function: all functions associated with a population model
		author:   <amikezor>
		date:     feb 04
***/

#include "population.h"
#include "sample.h"
#include "random.h"
#include <float.h>
#include <stdio.h>

#include <iostream>
using namespace std;

#define MIN2( a, b )  (((a)<(b))?(a):(b))
#define MAX2( a, b )  (((a)>(b))?(a):(b))

population::population(double Ne, double mu, char ploidy){

	this->Ne = Ne;
	this->mu = mu;
	this->set_ploidy(ploidy);
	this->compute_Theta( );

	this->Ti=NULL;
	this->fN=NULL;
	
}


population::~population( void ){

	if(this->Ti != NULL){
		delete this->Ti;
		this->Ti=NULL;
	}
		
	
	if(this->fN != NULL){
		delete this->fN;
		this->fN=NULL;
	}
	
}

void population::compute_Theta( void ){

	this->Theta = 2*this->ploidy*this->Ne*this->mu;
}


void population::set_ploidy( char p ){

	if(p==1 || p==2){
		this->ploidy=p;
	}
	else
		cerr << "population::set_ploidy: the ploidy can only be 1 or 2", exit(2);
}



void population::set_isolation( struct isolation_event *input_Ti, double *input_fN, long ns, double M ){

	if(ns<2){
		cerr << "population::set_speciation: this speciation/isolation makes only sense when ns >=2, bye\n";
		exit(4);
	}
	if(M<0){
		cerr  << "population::set_speciation: migration rate makes only sense when M >=0, bye\n";
		exit(4);
	}


	this->ns=ns;
	this->M=M;
	
	this->fN = new double [(int) (2*ns - 1) ];
	this->Ti = new struct isolation_event [(int) (ns - 1) ];
	
	for(int i=0 ; i< 2*ns - 1 ; i++)
		this->fN[i] = input_fN[i];
	
	for(int i=0; i<ns - 1;i++)
		this->Ti[i] = input_Ti[i];


}

void population::set_structure( double *input_fN, long ns, double M ){

	if(ns<2){
		cerr << "population::set_structure: this speciation/isolation makes only sense when ns >=2, bye\n";
		exit(4);
	}
	if(M<0){
		cerr  << "population::set_structure: migration rate makes only sense when M >=0, bye\n";
		exit(4);
	}


	this->ns=ns;
	this->M=M;
	
	this->fN = new double [(int) (2*ns - 1) ];
	
	for(int i=0 ; i< 2*ns - 1 ; i++)
		this->fN[i] = (input_fN == NULL)?1:input_fN[i];
	

}


/*
	Go through the whole Ti array
	and get the closest time
*/
int population::get_isolation_next_Ti( double time ){

	double min=-1;
	int i_next_event=-1;
	
	for( int i=0 ; i<this->ns-1 ; i++ ){
	
		if(min<0 && (this->Ti[i].time > time) ){
			min = this->Ti[i].time - time;
			i_next_event=i;
		}
		
		if(  this->Ti[i].time > time &&   this->Ti[i].time - time < min){
			min = this->Ti[i].time - time;
			i_next_event=i;
		}
	
	}
	
	//cerr <<"next event is "<<i_next_event<<" at time "<< this->Ti[i_next_event].time<<" between "<<this->Ti[i_next_event].sp[0] <<" & " <<this->Ti[i_next_event].sp[1] <<" into "<<this->Ti[i_next_event].id_sp <<"\n";

	return i_next_event;

}

/*
	ns - number of spcies
	M - igration rate (symetric)
	b - birth rate
	d - death rate
*/
void population::BirthDeath_species_tree(  long ns, double M, double b, double d, double *input_fN, char ancestral_type ){

	extern Random R;


	if(ns<2){
		cerr << "population::set_speciation: this speciation/isolation makes only sense when ns >=2, bye\n";
		exit(4);
	}
	if(M<0){
		cerr  << "population::set_speciation: migration rate makes only sense when M >=0, bye\n";
		exit(4);
	}
	if( d>b ){
		cerr  << "population::set_speciation: only sur-critic process are handled, when b>=d, bye\n";
		exit(4);
	}
	
//	cerr <<"rates are b:"<< b <<" d: "<<d<<"\n";

	this->ns=ns;
	this->M=M;
	
	/*
		For now all species have the same Ne
	*/

	this->fN = new double [(int) (2*ns - 1) ];
	
	for(int i=0 ; i< ns ; i++)
		this->fN[i] = (input_fN == NULL)?1:input_fN[i];   // set the leaf sizes



	this->Ti = new struct isolation_event [(int) (ns - 1) ];
	

	/*
		Built the species tree <=> isolation_events
	*/

	int *alive_species;
	
	alive_species = (int *)malloc( (size_t) ns*sizeof(int) );
	if( ! alive_species  )fprintf(stderr, "built_species_tree: cannot allocate alive_species, bye\n"), exit(1);

	for(int i=0; i<ns; i++)
		alive_species[i]=i;

	int exceed=0;
	for(int i=0; i<ns-1; i++){             // for all internal nodes
		
		/*
			Draw time
		*/
		if( d == 0 )
			this->Ti[i].time = R.exponential_dev( 1.0/b );              // Yule model - exponential distrib, death=0
		else
		{
			if( b == d )
				this->Ti[i].time = R.cauchy_dev( b );              // critical model birth rate = death rate
			else
				this->Ti[i].time = R.birthdeath_dev( b, d );       // Birth/Death general model, birth != death
		}
		

		if(this->Ti[i].time > 10000){
			this->Ti[i].time=10000+exceed;
			exceed++;
		}
	}
		
	
	
	int species=ns;
	int i_next_Ti=0;
	double base_time=0;
	
	while(species < 2*ns-1 ){
	
	
		// 1. find next isolation event (closest time)
		
		i_next_Ti = get_isolation_next_Ti( base_time );
		
		
		// 2. find connection (next larger time for Ti[i], i>i_next_Ti ; if none connect to alive_species[ns-1]

		int i = i_next_Ti + 1;
		
		while( i < ns-1 &&  this->Ti[ i_next_Ti ].time > this->Ti[i].time ){   // stop when reach a longer time (at most ns-1, the last alive_sp)
			i++;
		}

	
		// 3. set  alive_species[ connect ] to species
		
		Ti[ i_next_Ti ].sp[ 0 ] = alive_species[ i_next_Ti ];
		Ti[ i_next_Ti ].sp[ 1 ] = alive_species[ i ];
		Ti[ i_next_Ti ].id_sp = species;

		switch(ancestral_type){
		
			case 'm':
				fN[species] = ( fN[Ti[ i_next_Ti ].sp[ 0 ]] <  fN[Ti[ i_next_Ti ].sp[ 1 ]] )?fN[Ti[ i_next_Ti ].sp[ 0 ]]:fN[Ti[ i_next_Ti ].sp[ 1 ]] ;
				break;
			case 'M':
				fN[species] = ( fN[Ti[ i_next_Ti ].sp[ 0 ]] >  fN[Ti[ i_next_Ti ].sp[ 1 ]] )?fN[Ti[ i_next_Ti ].sp[ 0 ]]:fN[Ti[ i_next_Ti ].sp[ 1 ]] ;
				break;
			case 's':
				fN[species] = fN[Ti[ i_next_Ti ].sp[ 0 ]] + fN[ Ti[ i_next_Ti ].sp[ 1 ]] ;
				break;
		
		}
		
		alive_species[ i ] = species;
		
	

		// 4. species ++
		
		base_time = Ti[i_next_Ti].time;
		species++;

	}
	

	free( alive_species );

}



/*
	Built a species tree according to a Moran Process
	it is Kingman like although we assume we have the whole clade
	at the begining of the process

	ns - number of spcies
	M - igration rate (symetric)
	r - coalescent rate
*/
void population::Moran_species_tree(  long ns, double M, double r, double *input_fN, char ancestral_type ){

	extern Random R;


	if(ns<2){
		cerr << "population::Moran_species_tree: this speciation/isolation makes only sense when ns >=2, bye\n";
		exit(4);
	}
	if(M<0){
		cerr  << "population::Moran_species_tree: migration rate makes only sense when M >=0, bye\n";
		exit(4);
	}
	

	this->ns=ns;
	this->M=M;
	
	/*
		For now all species have the same Ne
	*/
	this->fN = new double [(int) (2*ns - 1) ];
	
	for(int i=0 ; i< ns ; i++)
		this->fN[i] = (input_fN == NULL)?1:input_fN[i];   // set the leaf sizes


	this->Ti = new struct isolation_event [(int) (ns - 1) ];
	

	/*
		Built the species tree <=> isolation_events
	*/

	int *alive_species;
	int n_alive=ns;
	
	double time_event = 0;
	int n1,
	    n2;
	    
	double coal_rate;
	
	alive_species = (int *)malloc( (size_t) ns*sizeof(int) );
	if( ! alive_species  )fprintf(stderr, "population::Moran_species_tree: cannot allocate alive_species, bye\n"), exit(1);

	for(int i=0; i<ns; i++)
		alive_species[i]=i;


	while(n_alive>1){
			
		double tt;
		
		coal_rate  = r * ( ( n_alive *(n_alive-1.00) ) / ns );   // scale rate
		tt =  R.exponential_dev( 1.0/coal_rate );
		
		time_event += tt;          // draw time

//		printf("-- nalive: %d tt %f time_event %f\n", n_alive, tt, time_event);

		R.uniform_2DiffInt(0, n_alive-1, n1, n2);                 // draw 2 species
	
		//cerr <<"n1: "<<n1 <<" ; n2: "<<n2<<"\n";
	
		if(n1>n2){
			int tmp=0; tmp=n1; n1=n2; n2=tmp;       /* swap n1,n2 */
		}
		

		Ti[ ns-n_alive ].sp[ 0 ] = alive_species[ n1 ];  // do the event
		Ti[ ns-n_alive ].sp[ 1 ] = alive_species[ n2 ];
		Ti[ ns-n_alive ].id_sp = 2*ns-n_alive;
		Ti[ ns-n_alive ].time = time_event;

		
		switch(ancestral_type){
		
			case 'm':
				fN[2*ns-n_alive] = ( fN[Ti[ ns-n_alive ].sp[ 0 ]] <  fN[Ti[ ns-n_alive ].sp[ 1 ]] )?fN[Ti[ ns-n_alive ].sp[ 0 ]]:fN[Ti[ ns-n_alive ].sp[ 1 ]] ;
				break;
			case 'M':
				fN[2*ns-n_alive] = ( fN[Ti[ ns-n_alive ].sp[ 0 ]] >  fN[Ti[ ns-n_alive ].sp[ 1 ]] )?fN[Ti[ ns-n_alive ].sp[ 0 ]]:fN[Ti[ ns-n_alive ].sp[ 1 ]] ;
				break;
			case 's':
				fN[2*ns-n_alive] = fN[Ti[ ns-n_alive ].sp[ 0 ]] + fN[ Ti[ ns-n_alive ].sp[ 1 ]] ;
				break;
		
		}
		
		
		
		
		
		alive_species[n1]=2*ns-n_alive;                   // update alive_species array

		while(n2<n_alive-1){
			alive_species[n2]=alive_species[n2+1];
			n2++;
		}
		
		n_alive--;

	}

	free( alive_species );

}



/*
	Built a species tree according to a Kingman Process

	ns - number of spcies
	M - igration rate (symetric)
	r - coalescent rate
*/
void population::Kingman_species_tree(  long ns, double M, double r, double *input_fN, char ancestral_type ){

	extern Random R;


	if(ns<2){
		cerr << "population::set_speciation: this speciation/isolation makes only sense when ns >=2, bye\n";
		exit(4);
	}
	if(M<0){
		cerr  << "population::set_speciation: migration rate makes only sense when M >=0, bye\n";
		exit(4);
	}
	

	this->ns=ns;
	this->M=M;
	
	/*
		For now all species have the same Ne
	*/
	this->fN = new double [(int) (2*ns - 1) ];
	
	for(int i=0 ; i< ns ; i++)
		this->fN[i] = (input_fN == NULL)?1:input_fN[i];   // set the leaf sizes


	this->Ti = new struct isolation_event [(int) (ns - 1) ];
	

	/*
		Built the species tree <=> isolation_events
	*/

	int *alive_species;
	int n_alive=ns;
	
	double time_event = 0;
	int n1,
	    n2;
	    
	double coal_rate;
	
	alive_species = (int *)malloc( (size_t) ns*sizeof(int) );
	if( ! alive_species  )fprintf(stderr, "built_species_tree: cannot allocate alive_species, bye\n"), exit(1);

	for(int i=0; i<ns; i++)
		alive_species[i]=i;


	while(n_alive>1){
			
		double tt;
		
		coal_rate  = r * ( ( n_alive *(n_alive-1.00) ) / 2.00 );   // scale rate
		tt =  R.exponential_dev( 1.0/coal_rate );
		
		time_event += tt;          // draw time

//		printf("-- nalive: %d tt %f time_event %f\n", n_alive, tt, time_event);

		R.uniform_2DiffInt(0, n_alive-1, n1, n2);                 // draw 2 species
	
		//cerr <<"n1: "<<n1 <<" ; n2: "<<n2<<"\n";
	
		if(n1>n2){
			int tmp=0; tmp=n1; n1=n2; n2=tmp;       /* swap n1,n2 */
		}
		

		Ti[ ns-n_alive ].sp[ 0 ] = alive_species[ n1 ];  // do the event
		Ti[ ns-n_alive ].sp[ 1 ] = alive_species[ n2 ];
		Ti[ ns-n_alive ].id_sp = 2*ns-n_alive;
		Ti[ ns-n_alive ].time = time_event;

		
		switch(ancestral_type){
		
			case 'm':
				fN[2*ns-n_alive] = ( fN[Ti[ ns-n_alive ].sp[ 0 ]] <  fN[Ti[ ns-n_alive ].sp[ 1 ]] )?fN[Ti[ ns-n_alive ].sp[ 0 ]]:fN[Ti[ ns-n_alive ].sp[ 1 ]] ;
				break;
			case 'M':
				fN[2*ns-n_alive] = ( fN[Ti[ ns-n_alive ].sp[ 0 ]] >  fN[Ti[ ns-n_alive ].sp[ 1 ]] )?fN[Ti[ ns-n_alive ].sp[ 0 ]]:fN[Ti[ ns-n_alive ].sp[ 1 ]] ;
				break;
			case 's':
				fN[2*ns-n_alive] = fN[Ti[ ns-n_alive ].sp[ 0 ]] + fN[ Ti[ ns-n_alive ].sp[ 1 ]] ;
				break;
		
		}
		
		
		
		
		
		alive_species[n1]=2*ns-n_alive;                   // update alive_species array

		while(n2<n_alive-1){
			alive_species[n2]=alive_species[n2+1];
			n2++;
		}
		
		n_alive--;

	}

	free( alive_species );

}



/*
	Built a species tree according to a Radiation Process (all species at the same time)

	ns - number of spcies
	M - igration rate (symetric)
	r - coalescent rate
*/
void population::Radiation_species_tree(  long ns, double M, double r, double *input_fN, char ancestral_type ){

	extern Random R;


	if(ns<2){
		cerr << "population::set_speciation: this speciation/isolation makes only sense when ns >=2, bye\n";
		exit(4);
	}
	if(M<0){
		cerr  << "population::set_speciation: migration rate makes only sense when M >=0, bye\n";
		exit(4);
	}
	

	this->ns=ns;
	this->M=M;
	
	/*
		For now all species have the same Ne
	*/
	this->fN = new double [(int) (2*ns - 1) ];
	for(int i=0 ; i< ns ; i++)
		this->fN[i] = (input_fN == NULL)?1:input_fN[i];   // set the leaf sizes


	this->Ti = new struct isolation_event [(int) (ns - 1) ];
	

	/*
		Built the species tree <=> isolation_events
	*/

	int sp_min=0;
	int sp_max=ns-1;
	int event=0;
	
	double time_event = 0;
	
	time_event = R.exponential_dev( 1.0/r );

	while( sp_max-sp_min >= 1  ){
	
		Ti[ event ].sp[ 0 ] = sp_min;  // do the event
		Ti[ event ].sp[ 1 ] = sp_min+1;
		Ti[ event ].id_sp = sp_max+1;
		Ti[ event ].time = time_event;

		// this takes care of ancestral size
		
		switch(ancestral_type){
		
			case 'm':
				fN[sp_max+1] = ( fN[Ti[ event ].sp[ 0 ]] <  fN[Ti[ event ].sp[ 1 ]] )?fN[Ti[ event ].sp[ 0 ]]:fN[Ti[ event ].sp[ 1 ]] ;
				break;
			case 'M':
				fN[sp_max+1] = ( fN[Ti[ event ].sp[ 0 ]] >  fN[Ti[ event ].sp[ 1 ]] )?fN[Ti[ event ].sp[ 0 ]]:fN[Ti[ event ].sp[ 1 ]] ;
				break;
			case 's':
				fN[sp_max+1] = fN[Ti[ event ].sp[ 0 ]] + fN[ Ti[ event ].sp[ 1 ]] ;
				break;
		
		}


		time_event += 1e-10;            // this seems to make a difference in numbers (implementation tricks), though it is still very close to time_event

		sp_min+=2;
		sp_max++;
		event++;

	}

}



/*
	Built a species tree according to a Radiation Process (all species at the same time)

	ns - number of spcies
	M - igration rate (symetric)
	t - speciation time (the first event was t ago, older are every t in Ne generations time
*/
void population::Fixed_species_tree(  long ns, double M, double t, double *input_fN, char ancestral_type ){


	if(ns<2){
		cerr << "population::Fixed_species_tree: this speciation/isolation makes only sense when ns >=2, bye\n";
		exit(4);
	}
	if(M<0){
		cerr  << "population::Fixed_species_tree: migration rate makes only sense when M >=0, bye\n";
		exit(4);
	}
	

	this->ns=ns;
	this->M=M;
	
	/*
		set species Ne
	*/
	this->fN = new double [(int) (2*ns - 1) ];
	for(int i=0 ; i< ns ; i++)
		this->fN[i] = (input_fN == NULL)?1:input_fN[i];   // set the leaf sizes


	this->Ti = new struct isolation_event [(int) (ns - 1) ];
	

	/*
		Built the species tree <=> isolation_events
	*/

	int sp_min=0;
	int sp_max=ns-1;
	int event=0;
	
	double time_event = 0;
	
	time_event = t;

	while( sp_max-sp_min >= 1  ){
	
		Ti[ event ].sp[ 0 ] = sp_min;  // do the event
		Ti[ event ].sp[ 1 ] = sp_min+1;
		Ti[ event ].id_sp = sp_max+1;
		Ti[ event ].time = time_event;

		// this takes care of ancestral size
		
		switch(ancestral_type){
		
			case 'm':
				fN[sp_max+1] = ( fN[Ti[ event ].sp[ 0 ]] <  fN[Ti[ event ].sp[ 1 ]] )?fN[Ti[ event ].sp[ 0 ]]:fN[Ti[ event ].sp[ 1 ]] ;
				break;
			case 'M':
				fN[sp_max+1] = ( fN[Ti[ event ].sp[ 0 ]] >  fN[Ti[ event ].sp[ 1 ]] )?fN[Ti[ event ].sp[ 0 ]]:fN[Ti[ event ].sp[ 1 ]] ;
				break;
			case 's':
				fN[sp_max+1] = fN[Ti[ event ].sp[ 0 ]] + fN[ Ti[ event ].sp[ 1 ]] ;
				break;
		
		}


		time_event += t;            // this seems to make a difference in numbers (implementation tricks), though it is still very close to time_event

		sp_min+=2;
		sp_max++;
		event++;

	}

}



void population::free_population_arrays( void ){

	if ( this->Ti != NULL){
		delete this->Ti ;
		 this->Ti=NULL;
	}
	
	if (  this->fN != NULL){
		 delete this->fN ;
		 this->fN=NULL;
	}
	

}



/*
	print out the tree
	adpated from my previous C code.
	to print the whole tree, the first pointer has to be at the top of the tree !!!
*/
static void treeout( locus *ptree, locus *pref ){

	locus *ptmp;

	if( ptree == NULL )return;                                   // stop that recursive function

	if(ptree->get_descendant1() != NULL){                        // it is not an final node
		cout << "(";
		ptree = ptree->get_descendant1();
	} 
	else{                                                        // it is a final node
		
		printf("sp_%ld:%.10f", ptree - pref, ptree->get_ancestor()->get_time());
		
		
		ptmp = ptree;
		
		while(ptmp->get_ancestor()->get_descendant2() == ptmp){   // are we a descendant2 ??
			ptmp=ptmp->get_ancestor();
	
			if( ! ptmp->is_coalesced() ){                          // we are at the the top of the tree
				cout << ");\n";
				ptree=NULL;
				return;
			}                                                     // in tree reconstruction, we want the branch length
			
			printf(")sp_%ld:%.10f",  ptmp - pref, ptmp->get_ancestor()->get_time() - ptmp->get_time());  
		}
		
		ptree=ptmp->get_ancestor()->get_descendant2();          // switch to desc2
		cout << ",";
	}
	
	treeout( ptree, pref );
}



void population::print_isolations( void ){

	
	locus *array = new locus[ 2*ns-1 ];


	for(int i=0;i< this->ns-1; i++){
	
		//printf("connect %d+%d -> %d at t= %f\n", Ti[i].sp[0], Ti[i].sp[1], Ti[i].id_sp, Ti[i].time );
	
		array[ Ti[i].sp[0] ].connect_ancestor( array+Ti[i].id_sp , 1 );
		array[ Ti[i].sp[1] ].connect_ancestor( array+Ti[i].id_sp , 2 );
		array[ Ti[i].id_sp ].set_time( Ti[i].time );
	
	}

	treeout( array+2*ns-2, array  );

	/*
	cout << "n_species:   "<< this->ns<<"\n";

	cout << "isol#\ttime\tsp1\tsp2\t->sp_anc\n";
	
	for(int i=0;i< this->ns-1; i++)
		printf("%d\t%.3f\t%d\t%d\t%d\n", i+1,this->Ti[i].time ,this->Ti[i].sp[0],this->Ti[i].sp[1],this->Ti[i].id_sp );
	
	*/
		
	delete[] array;

}


#undef MIN2
#undef MAX2
