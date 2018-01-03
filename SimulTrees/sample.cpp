/***
		file:     sample.cpp
		function: all functions associated with the samples
		author:   <amikezor>
		date:     feb 04,
		modif :   nov 06, add some extra check in set_loci()
		          oct 08 when the sweep is not even reached
		          nov 08, the isolation/structure
		          feb 09, sample in an ogoing sweep
***/

#include "sample.h"
#include "population.h"
#include "random.h"
#include "math.h"
#include "misc.h"

#include <iostream>
#include <vector>
using namespace std;

#include <stdio.h>

int sample::Count_All_Samples=0;            // number of sample objets in use


/*
	Constructor
*/
sample::sample( population *p ){

	this->id = Count_All_Samples++;          // increment sample objects
	this->pop = p;                           // this has been sampled of that population
		
	this->time=0;
	this->ancestor=NULL;                     // se to null/0 for now
	this->descendant=NULL;
	this->loci=NULL;
	this->nloci=0;

	this->nspecial=0;                         // a tag (for now, used for the # of lineage that has escaped the sweep in the oldest sample)
	
//	cerr << "//\nCreated a sample of id "<<this->id<<" and time "<<time <<"\n//\n";
}

/*
	Destructor
*/
sample::~sample( void ){
	
	//cerr << "// We will destroy a sample of id "<<this->id<<" and time "<<this->time <<"\n";
	//cerr << "// with "<<nloci<<" loci\n";

	if( this->ancestor != NULL )            // clean all ancestors
		delete this->ancestor;
	
	
	for(int i=0; i< this -> n_original_loci; i++)           // delete all first n_original_loci -- the ones that were allocated here
		delete loci[i];
			
			
	delete[] loci;                          // free the array of pointer
		
	if(descendant)                          // set the descendant ancestor as free
		descendant->ancestor=NULL;
		
	Count_All_Samples--;
}



/*
	Set number of loci of this sample
	n is the number of loci
	n_memory is the number of new loci object you want to be instanciated
	
	the first n_memory first ones will be initiated.
	
	n_memory shall be <= n
	
	the first one is 0
*/
void sample::set_loci( int n, int n_memory ){

	int a;

	if(n_memory>n)
		cerr<<"loci::set_loci, cannot allocate more loci (n_memory) than their numbers (n), bye\n", exit(3);
	

	this->loci = new locus * [n];            // get the array of pointer
	if(this->loci == NULL)
		cerr<<"loci::set_loci, allocate error, bye\n", exit(3);
	this->nloci=n;

	if(n_memory>0)
		for(a=0;a<n_memory;a++){
			this->loci[a] = new locus( 0 );      // get n_memory new loci objects
			if(this->loci[a] == NULL)
				cerr<<"loci::set_loci, loci[a] allocate error, bye\n", exit(3);
		}
	
	this -> n_original_loci = n_memory;
	
	
	for(a=n_memory;a<n;a++)
		this->loci[a] = NULL;                   // set the rest of pointer to free
}



/*
	Connect a sample to its ancestor
*/
void sample::connect_ancestor(sample *anc){
	this->ancestor = anc;
	anc->descendant = this;
}

/*
	How many locus within a chosen species
	and existing at the time of the sample
*/
int sample::get_nloci_species( int species, double current_time  ){

	int n=0;
	
	for(int i=0;i<nloci;i++)
		if (loci[i]->get_species() == species && loci[i]->get_time() <= current_time )
			n++;
	
	return n;
}



locus **sample::get_loci_species( int species, double current_time  ){

	int nloc = get_nloci_species( species  );
	int n=0;

	locus **mylocus = new locus *[ nloc ];
	
	for(int i=0;i<nloci;i++)
		if (loci[i]->get_species() == species && loci[i]->get_time() <=  current_time )
			mylocus[ n++ ] = loci[i];
	
	return mylocus;
}


/*
	within the chosen species, return the x° locus
	x is between 0 and n_alive-1
*/
locus *sample::get_species_locus_x( int _chosen_species_, int x, double current_time ){

	int tmp_n=0;
	
	if(  x >= get_nloci_species( _chosen_species_ , current_time ) ){
		cerr << "sample::get_species_locus_x cannot retrieve locus "<<x<<" of species "<<_chosen_species_<<", because this species has only "<<get_nloci_species( _chosen_species_  )<<" loci\n";
		return NULL;
	}

	for(int j=0; j<nloci ; j++ ){
				
		if( loci[j]->get_species() == _chosen_species_  && loci[j]->get_time() <= current_time ){
			
			if( tmp_n == x )
				return loci[j];
			
			tmp_n++;
		}
	}
					
	return NULL;
}


/*
	Draw randomnly a locus within a species
*/
locus *sample::get_random_locus_species( int _chosen_species_, double current_time ){

	int rand_loc_sp;
	extern Random R;                                               // an object to draw random variables declared in random.h
	
	rand_loc_sp = R.uniform_int(0, get_nloci_species( _chosen_species_, current_time  )-1 );
	return get_species_locus_x(  _chosen_species_, rand_loc_sp, current_time );

}

/*
	How many species are alive at that time
	in the int *ni vector
*/
int sample::get_n_alive_species( int *ni ){
	
	int i;
	int n_alive_species=0;

	for(i=0;i<2*this->pop->get_isolation_ns() -1; i++)
		if(ni[i]>0){
			n_alive_species++;
		}

	return n_alive_species;
}

/*
	From all alive species (in the int *ni vector)
	return the absolute number of the x-iest.
	(x in 0, n_alive-1)
*/
int sample::get_alive_species_x( int x, int *ni ){


	int tmp_x=0;
	
	if( x >= get_n_alive_species( ni ) || x< 0 ){
		cerr << "cannot retrieve alive species "<<x<<", because there are only "<< get_n_alive_species(ni) <<" species\n";
		return -1;
	}

	for(int i=0;  i<2*this->pop->get_isolation_ns() -1 ; i++){
	
		if( ni[i] > 0 ){
		
			 if(tmp_x == x)
			 	return i;
				
			 tmp_x++;
		}
			 
	}
	
	cerr << "sample::get_alive_species_x: Error ! It should not reach this point\n";
	return -1;

}
/*
	Return random a species alive at that time
*/
int sample::get_random_alive_species( int *ni ){

	extern Random R;                                               // an object to draw random variables declared in random.h
	int rand_sp;
	
	rand_sp = R.uniform_int(0, get_n_alive_species( ni )-1 );       // a random sp # between 1 and n_alive_sp
	
	return get_alive_species_x( rand_sp, ni );

}

/*
	Return random a species alive at that time
*/
int sample::get_random_species_bysize( void ){

	extern Random R;                                               // an object to draw random variables declared in random.h
	int rand_sp;
	double dum,c;
	double total_size=0;
	
		
	if( pop->get_isolation_fN() == NULL )
	{
	
		rand_sp = R.uniform_int(0, pop->get_isolation_ns()-1 );
	
	}
	else
	{
	
		dum = R.uniform_dev();
		
		total_size=0;
		for(int i=0;i<pop->get_isolation_ns();i++)
			total_size += pop->get_isolation_fN(i);
	
		c = pop->get_isolation_fN(0)/total_size;
		rand_sp=0;

		while( dum > c )
		{
			rand_sp++;
			c += pop->get_isolation_fN(rand_sp)/total_size;
		}
		
		if(c>1){
			cerr << "sample::get_random_species_bysize: this value should  never exceed 1.00 ; check code\n";
			exit(5);
		}
	
	
	}
	
	
	
	return rand_sp;

}


/*
	Change time for this sample and of the locus[0] - locus[0] being the loci that is new to this sample
	The time reduction is implemented as described in Simonsen et al., 1995, Genetics
*/
void sample::bottleneck_reduce_time_sample( void ){

	
	//cerr <<"for sample "<<get_id()<<" ; time is "<<time<<" ...";

	if( time > pop->get_bottleneck_Tb() && time <= (pop->get_bottleneck_Tb() + pop->get_bottleneck_Tl() / pop->get_bottleneck_f() ) )
		time =  pop->get_bottleneck_Tb() + (time - pop->get_bottleneck_Tb() )*pop->get_bottleneck_f();                             // time becomes Tb + (t-Tb)*f
	else{
		if( time > pop->get_bottleneck_Tb() + pop->get_bottleneck_Tl() / pop->get_bottleneck_f() )
			time =  time -  ( (1.00/pop->get_bottleneck_f() ) - 1.00 )*pop->get_bottleneck_Tl();                                // time becomes t - (1/f -1)*l
	}

	//cerr <<" and now, time is "<<time<<"\n";


	loci[0]->set_time( this->time );                                                                           // change also the time of the first locus (the one specific to this sample )


	double branch1 = ( loci[0]->get_descendant1() ) ? time - loci[0]->get_descendant1()->get_time() : 0;       // the branch length to desc1
	double branch2 = ( loci[0]->get_descendant2() ) ? time - loci[0]->get_descendant2()->get_time() : 0;       //   "           "      desc2

	double Ttime_desc1 =  ( loci[0]->get_descendant1() ) ? loci[0]->get_descendant1()->get_Ttime() : 0;	   // the branch length to desc1
	double Ttime_desc2 =  ( loci[0]->get_descendant2() ) ? loci[0]->get_descendant2()->get_Ttime() : 0;	   //	"	    "	   desc2


	loci[0]->set_Ttime( branch1+branch2+Ttime_desc1+Ttime_desc2 );                                             // reset Ttime for locus[0];

}



/*
	Change time for this sample and of the locus[0] - locus[0] being the loci that is new to this sample
	The time reduction is done according to exponential/linear growth time - scale (see any good book ;o)
	
	type is 'l'inear or 'e'xponential
	
*/
void sample::growth_reduce_time_sample( char type ){

	extern char debug;
	double Gr=0,
	       a=0;

	if(debug)cerr <<"for sample "<<get_id()<<" ; time is "<<time<<" ...";
	
	switch( type ){
	
		case 'e':
			Gr = pop->get_exponential_growth_Gr();
			time = log( 1.00 + Gr * time ) / Gr;
			break;
	
		case 'l':
			a = pop->get_linear_growth_Gr();
			time = ( 1.0 - exp( - a * time  ) ) / a;
			break;
		
		default:
			cerr <<"growth_reduce_time_sample: unknown type\n";
			exit(1);
	
	}


	if(debug)cerr <<" and now, time is "<<time<<"\n";


	loci[0]->set_time( this->time );                                                                           // change also the time of the first locus (the one specific to this sample )


	double branch1 = ( loci[0]->get_descendant1() ) ? time - loci[0]->get_descendant1()->get_time() : 0;       // the branch length to desc1
	double branch2 = ( loci[0]->get_descendant2() ) ? time - loci[0]->get_descendant2()->get_time() : 0;       //   "           "      desc2

	double Ttime_desc1 =  ( loci[0]->get_descendant1() ) ? loci[0]->get_descendant1()->get_Ttime() : 0;	   // the branch length to desc1
	double Ttime_desc2 =  ( loci[0]->get_descendant2() ) ? loci[0]->get_descendant2()->get_Ttime() : 0;	   //	"	    "	   desc2


	loci[0]->set_Ttime( branch1+branch2+Ttime_desc1+Ttime_desc2 );                                             // reset Ttime for locus[0];

}


void sample::coalesce2loci( int n1, int n2, double coal_time){

	locus *l1 = loci[n1];
	locus *l2 = loci[n2];
	
	extern char debug;


	/*
		Coalesce at the locus level
	*/
	locus *panc=ancestor->loci[0];                                      // The first locus of the ancestor loci is the ancestral one

	panc->reset_connections();
	panc->set_time(coal_time);                                         // this new node started at time_event

	l1->connect_ancestor( panc, 1 );                                   // and is connected to l1 and l2
	l2->connect_ancestor( panc, 2 );
	
	if(l1->get_species() == l2->get_species() )
		panc->set_species(l1->get_species());                      // inherits the linkage/species/type if both are compatible
	else{
		 
		panc->set_species(-1);                                      // otherwise set it to 0  --undefined--
	
		if(debug){
			cerr << "Warning: coalescing lineage with different types\n";
			cerr << "loci that will be coalesced "<<n1<<" " << n2 << "\n";
			cerr << "time is "<< coal_time << "\n";
			cerr << "nloci is "<< nloci << "\n";
			cerr << "species are "<<l1->get_species()<<" "<<l2->get_species()<<"\n";
		}
		
	}

	
		
	double branch1 = coal_time - l1->get_time();
	double branch2 = coal_time - l2->get_time();
	
	panc->set_Ttime(  branch1 + l1->get_Ttime()  + branch2 + l2->get_Ttime()); 


	/*
		Copy other loci
	*/
	ancestor->time = coal_time;
	int a=1;                                                            // a=0 is reserved for the ancestor
	for( int q=0;q<nloci;q++ ){
		if ( q != n2 && q!= n1 )
			this->ancestor->loci[a++] = loci[q];                          // copy pointer to the uncoalesces loci
	}

//	cerr<<"--> Coalesced nodes "<< n1 <<" & "<<n2 <<" at time " <<coal_time<<"\n";

}


void sample::coalesce2loci_species( int nl1, int nl2, double coal_time, int _chosen_type_){

	locus *l1 = get_species_locus_x( _chosen_type_, nl1, coal_time  );
	locus *l2 = get_species_locus_x( _chosen_type_, nl2, coal_time  );

	/*
		Coalesce at the locus level
	*/
	locus *panc=this->ancestor->loci[0];                          // The first locus of the ancestor sample is the ancestral one
	
	panc->reset_connections();

	panc->set_time(coal_time);                                   // this new node started at time_event
	
	panc->set_species(_chosen_type_);                            // inherits the linkage/species/type

	l1->connect_ancestor( panc, 1 );                             // and is connected to l1 and l2
	l2->connect_ancestor( panc, 2 );
	
	
	double branch1 = coal_time - l1->get_time();
	double branch2 = coal_time - l2->get_time();
	
	panc->set_Ttime(  branch1 + l1->get_Ttime()  + branch2 + l2->get_Ttime()); 


	/*
		Copy other loci from the sample
	*/
	ancestor->time = coal_time;
	int a=1;                                                            // a=0 is reserved for the ancestor
	for( int q=0;q<nloci;q++ ){
		if ( ! loci[q]->is_coalesced() )
			this->ancestor->loci[a++] = loci[q];                // copy pointer to the uncoalesces loci
	}

//	cerr<<"--> Coalesced nodes "<< l1->get_id() <<" & "<<l2->get_id() <<" at time " <<time_event<<"\n";

}


int sample::get_nloci_befT( double t ){

	int nloci_befT = this->nloci;
	
	for(int i=nloci-1;i>=0;i--)
		if( this->loci[i]->get_time() > t )
			nloci_befT--;
		else
			break;

	return nloci_befT;
}


double sample::get_next_serial_sample( double time ){

	double next_time=0;

	for(int i=0;i<nloci;i++ )
		if( this->loci[i]->get_time() > time ){
			next_time = this->loci[i]->get_time();
			break;
		}

	return next_time;

}


/*
	Coalesce 2 loci chosen at random
*/
void sample::coalesce_sample( void ){


	if(nloci == 1)return;                                          // if only 1 locus, it cannot coalesce

	double time_event = 0;
	double tmp_time=0;

	int nloci_befT = get_nloci_befT(  this->time );
	
	 
 	if(nloci_befT <= 1 && nloci>1 )
	{
		time_event = loci[1]->get_time( );             // the next with some time as loci[0] is the last with time smaller than this->time
		nloci_befT = get_nloci_befT( time_event );
	}
 	else
		time_event = this->time;
 
 
	extern Random R;                                               // an object to draw random variables declared in random.h
	extern char debug;

	/*
		Draw time
	*/
	
	
	double rate = (double) ( nloci_befT *(nloci_befT-1) ) / 2.00;
	tmp_time = R.exponential_dev( 1.0/rate );
	
	while( get_nloci_befT(  time_event+tmp_time ) > nloci_befT )
	{
		// extract the first time between time_event and time_event+tmp_time
		// set time_event at that time
		// update rates
		// and then redraw with the new rate
		
		time_event = get_next_serial_sample( time_event );
		nloci_befT = get_nloci_befT( time_event );
		rate = (double) ( nloci_befT *(nloci_befT-1) ) / 2.00;
		tmp_time = R.exponential_dev( 1.0/rate );
			
	}
	
	time_event += tmp_time;

	

	if(debug)cerr << "rate is " << rate << " and time_event is " << time_event << " ; nlocus_befT " << nloci_befT <<"\n";

	/*
		Draw locus
	*/
	int n1, n2;

	R.uniform_2DiffInt(0, nloci_befT-1, n1, n2);

	
	this->coalesce2loci(  n1,  n2, time_event);
	
				
}


/*
	Coalesce 2 loci chosen at random
*/
void sample::coalesce_sample_fixedtime( double time ){


	if(nloci == 1)return;                                          // if only 1 locus, it cannot coalesce


	/*
		Draw locus
	*/
	extern Random R;                                               // an object to draw random variables declared in random.h
	int n1, n2;

	R.uniform_2DiffInt(0, get_nloci()-1, n1, n2);
	this->coalesce2loci(  n1,  n2, time);
				
}





/*
	Coalesce k loci chosen at random
*/
sample * sample::bs_coalesce_sample( void ){

	if(nloci == 1)return NULL;                                          // if only 1 locus, it cannot coalesce
 
	extern Random R;                                               // an object to draw random variables declared in random.h
//	extern char debug;

	double epsilon = 1e-20;
	int n1, n2;
	double time_event;
	
	sample *psample = this;

	/*
		Set rate
	*/
	
	double *l_rate = new double[ this->nloci-1 ];
	double total_rate = 0;
	
	for(int i=2;i<=this->nloci;i++){
		l_rate [ i-2 ] =  (this->nloci+0.0) / ( i * (i-1) );
		total_rate += l_rate [ i-2 ];	
	}

	
	total_rate = this->nloci -1;
	
	/*
		Draw time
	*/
	time_event = this->time + R.exponential_dev( 1.0/total_rate );


	/*
		Draw # of coalesced loci
	*/

	double U = R.uniform_dev();
	
	int ncoal=2;
	double tmp = l_rate [ 0 ]/total_rate;

	while( tmp < U )
	{
		ncoal++;
		tmp += l_rate [ ncoal-2 ] / total_rate;
	}
	
	/*
		Draw the two first loci
	*/

	R.uniform_2DiffInt(0, nloci-1, n1, n2);
	psample->coalesce2loci(  n1,  n2, time_event);  // coalesce them
	
	/*
		If more than 2, coalescen the otehr ones
	*/
	
	ncoal -= 2;
	int x=1;
	
	while( ncoal > 0){
	
	
		time_event += epsilon;
		n1 = 0;
		n2 = R.uniform_int(1, nloci-1-x);
	
		psample = psample->get_ancestor();
		
		psample->coalesce2loci(  n1,  n2, time_event);
		
		x++;
		ncoal --;
	
	}

	
	return psample->get_ancestor();
			
}

/*
	When a coalesce must happen within a subtree
	from Wiuf and Donelly, 1999
*/
void sample::nested_coalesce_sample( void ){                   // coalesce with one subtree of size k

	
	if(nloci == 1)return;                          // if only 1 locus, it cannot coalesce
 
	extern Random R;                               // an object to draw random variables declared in random.h
//	extern char debug;

	int k=0;
	int i;
	int sp=-1;

	for( i=0 ; i<this->nloci ; i++)                     // count how many species of type 1 --this is k
		if( this->loci[i]->get_species() == 1 )
			k++;

	if(k<2){                                           // the subtree is of size 1 or 0, so ignore it
		coalesce_sample( );
		
		if(k==1 && nloci>2)
			this->ancestor->loci[0]->set_species( 2 );   // in this case the whole subtree had coalesced but not the outtree
		
		return;
	}


	/*
		Draw time
	*/
	
	double rate = (double) ( nloci *(nloci-1) ) / 2.00;
	double time_event = this->time + R.exponential_dev( 1.0/rate );


	/*
		Draw locus within one of the two sub-sample
	*/
	int n1, n2;

//	float P = k*(k-1.0)/( k*(k-1.0)+(nloci-k)*(nloci-k-1.0));  /* the probability is given by the number of pairs within k or within n-k */

	float P = (k+1.0)/nloci;  /* the probability is given by Wiuf and Donelly, 1999, TBP */
	
	//cerr<< "k/n: "<<k<<"/"<<nloci << " ; P: "<<P<<"\n";


	if( R.uniform_dev() <  P )
	{
		R.uniform_2DiffInt(0, k-1, n1, n2);
		sp=1;
	}
	else
	{
		R.uniform_2DiffInt(0, nloci-k-1, n1, n2);
		sp=2;
	}

	coalesce2loci_species( n1, n2, time_event, sp);

	return;
}



/*
	Neutral coalescence back to Ts
	then sweep - completely selection driven only even for small frequencies (solution comes fom infinite pop) -
	then back to neutral again
	Follow Braverman et al. (1995), Kim and Wolfgang 2002 (without intergenik recombination)
	febr 09 - it can start when p, freq of the selected allle is p<1 (ongoing sweep)
	
*/
void sample::sweep_coalesce_sample( void ){

	if(nloci == 1)return;                                          // if only 1 locus, it cannot coalesce

	if( this->time < this->pop->get_sweep_TS() ){
		pop->set_sweep_over(false);
		oldest_sample()->set_nspecial( -1 );                  // mark it as an impossible value
	}
		

	extern Random R;                                              // an object to draw random variables declared in random.h
	
	
	/*
		Draw time
	*/
	
	double rate = ( double ) ( nloci *(nloci-1) ) / 2.00;
	double time_event = this->time + R.exponential_dev( 1/rate );


	/*
		In the following case, we simply are in a neutral case
		(before/after the sweep)
		when sweep is over
		  If no more linked lineages are left, then proceed to neutral
		  else lets go inside the sweep
	*/


	
	if( (pop->is_sweep_over() == true && get_nloci_species(0)==0) || time_event < this->pop->get_sweep_TS() ){
	
		/*
			Draw locus and coalesce them, regular coalescent
		*/
		int n1, n2;

		R.uniform_2DiffInt(0, nloci-1, n1, n2);		   // n1 and n2 are given by adresss, although this is not cristal clear here :-)

		this->coalesce2loci(  n1, n2, time_event);
		
		if(nloci == 2 && time_event < this->pop->get_sweep_TS() ){
			oldest_sample()->set_nspecial( 0 );        // none even reach the sweep;
		}
		

		return;
	
	}

	
	/*
		From here after, we are IN the sweep
	*/
	
	double Dt = 1/(100*pop->get_sweep_alpha() );      // this would be the time slices
	double base_time=0;                               // In the sweep, time of each event is actually given by base_time
 

       /*
               Set base_time
       */
	if( time >= pop->get_sweep_TS() )                   // already in the sweep
		base_time = this->time;
	else
		if( time_event >= pop->get_sweep_TS() )
			base_time = pop->get_sweep_TS();    // first time in the sweep



	/*
		do sweep from base_time
	*/      

	int nlinked;                     // number of loci linked to the sweeping one, the type (encoded in a species_tag) is 0

	double x=0;		         // frequency of the sweeping alleles. It starts at 1-epsilon and finish at epsilon
	double epsilon = 1e-4;           // a small number. Sweep start at frequency 1-epsilon and end at epsilon.
  
	double P_coal_linked=0,	         // Proba that two of the linked loci coalesce
               P_rec_linked=0,	         // Proba that one linked loci becomes unlinked
               P_coal_unlinked=0,        // Proba that two of the unlinked loci coalesce
               P_rec_unlinked=0,	 // Proba that one unlinked loci becomes linked
               P_noevent=0,	         // 1 - All the above proba;
	       
               P_noevent_cumul=1.0;      // Proba that nothing happened since the the last event


	double dum;		         // a random variable uniformly drawn between 0 and 1

	bool _coalescent_event_=false;   // when it becomes true, we need to coalesce two lineages
	int _chosen_type_ = -1;          // this is the type of both loci that will coalesce - the alleles linked to the sweeping one are type '0', the other ones are type '1'


	double time_end = 0;            //  All this sweep freq and times simply comes from the resolution of the continuous model
	                                //  x(t) = e^mt/(e^mt + 1/x(0) -1)
					//  time is in N gen, therefore instead of m we have alpha
					//  extract t, with p=freq when stop ; x0 freq when starts
					//  t = ( ln(p) - ln(1-p) + ln(1-x0) -ln(x0) ) / m
					//  with x0=epsilon and p=1-epsilon, find Kim and Stephans results - there is a typo in Braverman et al. x(t)-

	if( pop->get_sweep_p()==1 )
		time_end = -( 2 / pop->get_sweep_alpha() )*log( epsilon ) + pop->get_sweep_TS();
	else{
		if( pop->get_sweep_TS() == 0 && pop->get_sweep_p()<1)
			time_end = ( 1.0/pop->get_sweep_alpha() )*(  log(pop->get_sweep_p()) - log(1.0-pop->get_sweep_p()) - log(epsilon) );
		else{
			cerr <<"sweep_coalesce_sample: before sweep we have p<1 and Ts>0: this should never happen !!, bye\n";
			exit(1);
		}

	}
	
	/*
		As long as there are no coalescent event, we are in the same 'sample' layer (of the tree structure)
	*/
	while( _coalescent_event_ == false ){


		do
			dum = R.uniform_dev();               // draw a unidev random
		while( dum == 1.00 );
		
		nlinked = get_nloci_species(0);              // update the number of linked loci
		
		
		P_noevent_cumul=1.0;
	
		/*
			For each time step (Dt), see if something happens
		*/
		
		while( P_noevent_cumul > dum ){

			_chosen_type_ = -1;                 // type is still unknown

			base_time += Dt;                    // time gets bigger
			
			if(base_time > time_end)            // but cannot go beyond the limit
				base_time = time_end;
			
			
			/*
				This also all comes from x=e^mt/(e^mt + 1/x0 -1) = x0/( x0+(1-x0)*e^-mt )
			*/
			x = epsilon / ( epsilon +  (1.00-epsilon) * exp(  pop->get_sweep_alpha() * ( base_time - time_end ) ) );

			if(x<0 || x>1){
				cerr << "!!! At time "<< base_time <<" ; freq is "<< x <<" which cannot be trusted, error\n";
				exit(1);
			}


			
			/*
				When frequency gets too small, it stops
			*/
			if(x <=  epsilon){
			
				if( oldest_sample()->nspecial == -1)
					oldest_sample()->set_nspecial( nloci-nlinked+1 );
				
				pop->set_sweep_over( true );                // it is over now

				if(nlinked>1){
					_chosen_type_ = 0;                  // do coalesce persisting linked lineages (just in case)
					_coalescent_event_=true;
					break;                              // jump to the last part of this function
					
				
				}else{
				                                             // if only 1 persist, change its tag to 1
					if(nlinked == 1){
						locus *my_locus = get_species_locus_x( 0, 0  );
						my_locus->set_species( 1 );
					}
	
					
					_chosen_type_ = 1;
					_coalescent_event_=true;                       // then do a coalescent event
					
					base_time += R.exponential_dev( 1/rate );      // after some regular neutral time

					break;

				}
						
				
			}
			
			
			
			/*
				Compute the probabilities of each kind of "interesting" event
			*/
			P_rec_linked	= nlinked * pop->get_sweep_R() * (1.00-x) * Dt;         // from link to unlink
			P_rec_unlinked  = (nloci-nlinked) * pop->get_sweep_R() * x * Dt;        // from unlink to link

			P_coal_linked   = nlinked*(nlinked-1) * Dt / (2.00 * x);
			P_coal_unlinked = (nloci-nlinked) * ((nloci-nlinked)-1) * Dt / (2.0 * (1.0-x) );
		
			P_noevent	= 1.00 - P_coal_linked -P_rec_linked-P_coal_unlinked-P_rec_unlinked;


			/*
				And the cumulative that something shall happen
			*/
			P_noevent_cumul*= P_noevent;
			
		
		}
		
		/*
			If we still do not know what kind of event it is, choose it
			If it is known, it means that the sweep is over and we just need to coalesce linked loci
			otherwise choose beteween REC and COAL
		*/
		if(  _coalescent_event_ == false ){
		
		
		       do
		               dum = R.uniform_dev();	// draw another uniform number 
		       while( dum == 1.00 );


			/*
				From the probabilities, set the kind of events (coal or rec) and which type of loci (0 or 1)
			*/

			if( dum <= (P_rec_linked+P_rec_unlinked)/(1.00-P_noevent) ){
	
				_chosen_type_ = ( dum <= P_rec_linked/(1.00-P_noevent) )?0:1;                   // select the type to rec -> that switch type
				
				int my_nlocus = R.uniform_int( 0, get_nloci_species(_chosen_type_)-1 );         // pick a random locus number within the chosen_type
				locus *my_locus = get_species_locus_x( _chosen_type_, my_nlocus  );            // set the pointer to the random chosen locus
				my_locus->set_species( ((_chosen_type_ == 0)?1:0) );                            // switch its type
				
														// and then restart the sweep.
	
			}else{
	
				_chosen_type_ = ( dum <= (P_rec_linked+P_rec_unlinked+P_coal_linked)/(1.00-P_noevent) )?0:1;   // select type to coalesce
				_coalescent_event_=true;
	
			}
			
			
		} // ENDOF if(  _coalescent_event_ == false ){
		

	}  // ENDOF while (coalescent == false)
		
	
	/*
		When we reach this point, a coalescent event HAS happened
	*/
		
	int nl1, nl2;
	
	R.uniform_2DiffInt( 0, get_nloci_species(_chosen_type_)-1, nl1, nl2);      // pick the two loci given the chosen type to coalesce
	this->coalesce2loci_species( nl1, nl2, base_time, _chosen_type_);          // and do coalesce them
	
	if( nloci == 2 && pop->is_sweep_over( )==false )
		oldest_sample()->set_nspecial( 1 );                                // when only 2 loci left and sweep is on, they coalesced to a single last one and there was NO "escaping" lineage
		
}


static void update_rates( int *ni, double **coal_rate, double **migr_rate, double *total_coal_rate, double *total_migr_rate, population *mypop ){

	(*total_coal_rate)=0;
	(*total_migr_rate)=0;

	if (*migr_rate != NULL)
		delete *migr_rate;
	
	if(*coal_rate != NULL)
		delete *coal_rate;

	(*coal_rate) = new double[ 2*(mypop->get_isolation_ns()) - 1 ];
	(*migr_rate) = new double[ 2*(mypop->get_isolation_ns()) - 1 ];

	if( *coal_rate== NULL || *migr_rate==NULL){
		cerr << "update_rates: memory too short for ni coal_rate or migr_rate, bye\n";
		exit(3);
	}
	

	
	for (int i=0; i< 2*mypop->get_isolation_ns() - 1 ; i++)
	{
	
		(*coal_rate)[i] = ni[i]*(ni[i]-1.0) / ( 2.0*mypop->get_isolation_fN(i));

		(*migr_rate)[i] = ni[i] * mypop->get_isolation_M() / 2.0;                       // a migrant can come from any species (including itself)

		(*total_coal_rate) += (*coal_rate)[i];
		(*total_migr_rate) += (*migr_rate)[i];
	}

}


/*
	With isolation/speciation and migration,
	
	either - coalesce 2 loci or migrate one loci
	there can be several species and several isolation
	nov 08
*/
void sample::isolation_coalesce_sample( void ){

	if(nloci == 1)return;                // if only 1 locus, it cannot coalesce
 
	extern Random R;                     // an object to draw random variables declared in random.h

	int *ni;                            // the array with # loci in each species it has the ns today species plus the ns-1 ancestral ones
	
	
	double *coal_rate=NULL;                   // this will contains all coalescent rates for each species
	double *migr_rate=NULL;                   // "       "              migration "        "        "
	
	double total_coal_rate=0;            // the total coal rate - for exponential draw
	double total_migr_rate=0;            // "         migration " - "

	
	int i_next_Ti=-1;                   // number of the next event. If all are done, it is set to -1;
	int i;                              // a counter
	
	
	int _chosen_species_=0;             // when something happen, in which species ?


	extern char opt_verbose;

	/*
		Set the # of the next isolation event
	*/
	i_next_Ti = this->pop->get_isolation_next_Ti( this->time );
	
//	cerr << "this_time is "<< this->time <<" i_next_Ti is "<<i_next_Ti <<"\n";
	
	if( i_next_Ti == -1 ){
		coalesce_sample();               // If time is before all isolation, species are merged into a single ancestral one
		return; 		         // Just do a regular coalescent
	}


	/*
		Get memory for the requiered arrays
	*/
	ni = new int[ 2*(this->pop->get_isolation_ns()) - 1 ];
		
	if(ni ==NULL ){
		cerr << "sample::isolation_coalesce_sample: memory too short for ni, bye\n";
		exit(3);
	}
		
	
	/*
		Set the ni array
	*/
	int sum_ni=0;
	for (i=0; i<2*(this->pop->get_isolation_ns()) - 1 ; i++){
		ni[i] =  get_nloci_species( i );
		sum_ni+=ni[i];
		
	}
	
	if(sum_ni != this->nloci ){
		cerr << "sample::isolation_coalesce_sample: Error, nloci should equal the sum of loci in each species\n";
		exit(4);
	}
		
	/*
		Compute all the rates
	*/
	
	update_rates( ni, &coal_rate, &migr_rate, &total_coal_rate, &total_migr_rate, this->pop );

//	cerr <<"** We have nloci: "<<sum_ni<<"\n";
//	cerr <<"** total_coal_rate: "<<total_coal_rate<<"\n";
	
	
	bool _coal_event_=false;              // become true when a coalescent event happen --> sample is then the previous one
	double time_event=this->time;         // when something happens
	
	int _destination_species_;            // for migration

	while( _coal_event_ == false ){
	
	
		/*
			 draw time of the event if any left
		*/
		if( (total_coal_rate+total_migr_rate) > 0 ){
		
		
			double t =R.exponential_dev( 1/(total_coal_rate+total_migr_rate) );
			time_event += t;
		
		}
		else{
			time_event = this->pop->get_isolation_Ti(i_next_Ti).time;
		}
		
		
//		cerr <<"time_event is "<<time_event<<" and isolation event is "<<  this->pop->get_isolation_Ti(i_next_Ti).time <<"\n";
//		if(this->nloci<90)
//			exit(5);
		
		if( i_next_Ti == -1 || time_event < this->pop->get_isolation_Ti(i_next_Ti).time ){     // we did not reach bckwd any isolation point or all are done
	
			double dum;
			
			dum = R.uniform_dev();          // draw a dum [0,1[ 
			

			if( dum < total_migr_rate/(total_coal_rate+total_migr_rate) ){       // a migration event
			
			//cerr << "a migration event\n";
			
				
				locus *rand_locus;                    // a random locus when the species us selected for migration
				double sum_P = 0;                     // a tmp value to pick which locus is migrating.
				
				/*
					choose the from species
				*/
				
				_chosen_species_=0;
				while(  dum >= sum_P+migr_rate[_chosen_species_]/(total_coal_rate+total_migr_rate)){
				
					sum_P+=migr_rate[_chosen_species_]/(total_coal_rate+total_migr_rate);
					
					if (sum_P >total_migr_rate/(total_coal_rate+total_migr_rate)){
						cerr << "sample::isolation_coalesce_sample: impossible Prob in migration. Look for a bug, bye\n";
						exit(5);
					}
					
					 _chosen_species_++;   	                               // migrate from this _chosen_species_
				}
				
				
				rand_locus = get_random_locus_species( _chosen_species_ );     // pick a locus in the selected species
				
				_destination_species_ = get_random_alive_species( ni );        // pick a destination species --based on non-null entries of ni

				rand_locus->set_species( _destination_species_ );              // migrate (change species tag)
				
				if(opt_verbose > 1 && _chosen_species_ != _destination_species_ )
					cerr << "Migration from species "<<_chosen_species_ <<" to species "<<_destination_species_ << " at time: "<< time_event <<"\n";
				

				/*
					Update ni and the rates
				*/
				
				ni[ _destination_species_ ] ++;
				ni[ _chosen_species_ ] --;
				
				update_rates( ni, &coal_rate, &migr_rate, &total_coal_rate, &total_migr_rate, this->pop );
				
				
				continue;                                                       // redraw another event (until a coal is drawn)

			}
			else{   // a coalescent event
				
				//cerr << "a coalescent event\n";

				double sum_P = total_migr_rate/(total_coal_rate+total_migr_rate);

				_chosen_species_=0;
				while(  dum >= sum_P+coal_rate[_chosen_species_]/(total_coal_rate+total_migr_rate) ){
								
					sum_P+=coal_rate[_chosen_species_]/(total_coal_rate+total_migr_rate);

					if (sum_P >1){
						cerr << "sample::isolation_coalesce_sample: impossible Prob in coalescence. Look for a bug, bye\n";
						exit(5);
					}
					
					_chosen_species_++;                                     //  coalesce within this _chosen_species_
				}
			
				_coal_event_ = true;                                           //  since it is a coal event, stop the process
				
				break;
			}
			
			
	
		}
		else{     // We reached a bcwd time point of isolation. Merge the species. and go back to draw a number
			
			int anc_species=-1;                        // the ancestral species into which the two isolated ones will be merged.
			locus *plocus;                             // a pointer to a locus
			
			/*
				Get destination ancestral_species
			*/
			anc_species = this->pop->get_isolation_Ti(i_next_Ti).id_sp; 
			
			/*
				Change all loci of both species to ancestral_species
			*/
			for(int species=0; species<2; species++){
			
				int current_sp = this->pop->get_isolation_Ti(i_next_Ti).sp[species];                   // a tmp current sp
			
				for(i=0; i<ni[ current_sp ] ; i++){
					plocus = get_species_locus_x(  current_sp, 0 );
					plocus->set_species(anc_species);
				}
			}


			/*
				Update ni and the rates
			*/
			
			ni[ anc_species ] = ni[ this->pop->get_isolation_Ti(i_next_Ti).sp[0] ]+ ni[ this->pop->get_isolation_Ti(i_next_Ti).sp[1] ];
			ni[ this->pop->get_isolation_Ti(i_next_Ti).sp[0] ]=0;
			ni[ this->pop->get_isolation_Ti(i_next_Ti).sp[1] ]=0;
			
			update_rates( ni, &coal_rate, &migr_rate, &total_coal_rate, &total_migr_rate, this->pop );
			
			/* 
				Set Time to next_Ti;
			*/
			
			time_event = this->pop->get_isolation_Ti(i_next_Ti).time;
			
			
			/*
				Update i_next_Ti;
			*/
			i_next_Ti = this->pop->get_isolation_next_Ti( time_event );

			if(i_next_Ti == -1){                                                      // all merging are done
				
				for(i=0;i<2*(this->pop->get_isolation_ns())-2 ; i++){              // therefore no species except the 2*(this->pop->get_isolation_ns())-1 should have loci
					
					//cerr << i << " " <<ni[i]<<"\n";
					
					if(ni[i]){
						cerr <<"sample::isolation_coalesce_sample. Error in species: "<< i<<" still have lineages (and it should not : BUG :o) !!\n";
						exit(5);
					}
				}
				
				total_migr_rate=0;                                                // no more migration are pertinent
			}
			
			/*
				Next round
			*/
			continue;
		
		}
	
	
	
	}
	
	/*
		At this point, we HAVE encountered a coalescent event
		within species = _chosen_species_.
	*/
	
	delete coal_rate;
	delete migr_rate;
	delete ni;
	
	int nl1, nl2;                         // two random loci

	R.uniform_2DiffInt( 0, get_nloci_species(_chosen_species_)-1, nl1, nl2);
	coalesce2loci_species( nl1, nl2, time_event, _chosen_species_);

		
//	cerr<<"--> Coalesced nodes "<< l1->get_id() <<" & "<<l2->get_id() <<" at time " <<time_event<<"\n";
		
}




/*
	With isolation/speciation and migration,
	
	either - coalesce 2 loci or migrate one loci
	there can be several species and several isolation
	nov 08
*/
void sample::structured_coalesce_sample( void ){

	if(nloci == 1)return;                // if only 1 locus, it cannot coalesce
 
	extern Random R;                     // an object to draw random variables declared in random.h

	int *ni;                            // the array with # loci in each species it has the ns today species plus the ns-1 ancestral ones
	
	
	double *coal_rate=NULL;              // this will contains all coalescent rates for each species
	double *migr_rate=NULL;              // "       "              migration "        "        "
	
	double total_coal_rate=0;            // the total coal rate - for exponential draw
	double total_migr_rate=0;            // "         migration " - "

	
	int i;                              // a counter
	
	
	int _chosen_species_=0;             // when something happen, in which species ?


	extern char opt_verbose;


	double time_event = this->time;         // when something happens

	int nloci_befT = get_nloci_befT(  this->time );



	if(nloci_befT <= 1 && nloci>1 )
	{
		time_event = loci[1]->get_time( );             // the next with some time as loci[0] is the last with time smaller than this->time
		nloci_befT = get_nloci_befT( time_event );
		
		cerr << "sample::structured_coalesce_sample only 1 locus was remaining --> update, now "<< nloci_befT << " \n";
		
	}


	/*
		Get memory for the requiered arrays
	*/
	ni = new int[ this->pop->get_isolation_ns() ];
	
	if(ni ==NULL ){
		cerr << "sample::structured_coalesce_sample: memory too short for ni , bye\n";
		exit(3);
	}
		
	
	/*
		Set the ni array - number loci per species
		that exist at that time (this last point is for serial coalescent)
	*/
	int sum_ni=0;
	
	for (i=0; i<this->pop->get_isolation_ns() ; i++){
		ni[i] =  this->get_nloci_species( i, time_event );
		sum_ni+=ni[i];
		
	}
	
	if(sum_ni != nloci_befT ){
		cerr << "sample::structured_coalesce_sample: Error, nloci_befT should equal the sum of loci in each subpopulations\n";
		exit(4);
	}

	//cerr <<"** lineages: "<<sum_ni<<"\n";
	//for(i=0; i<this->pop->get_isolation_ns() ; i++)
	//	cerr <<"**   ni["<<i<<"]: "<<ni[i]<<"\n";


	/*
		Compute all the rates
	*/
	
	update_rates( ni, &coal_rate, &migr_rate, &total_coal_rate, &total_migr_rate, this->pop );
	
	//cerr <<"** We have nloci: "<<sum_ni<<"\n";
	//cerr <<"** total_coal_rate: "<<total_coal_rate<<"\n";
	//cerr <<"** total_migr_rate: "<<total_migr_rate<<"\n";
	
	
	
	bool _coal_event_=false;              // become true when a coalescent event happen --> sample is then the previous one
	
	int _destination_species_;            // for migration

	while( _coal_event_ == false ){
	
	
		double t = R.exponential_dev( 1/(total_coal_rate+total_migr_rate) );		
	

		/*
			If there is an increase in the number of locus between
			time_event and time_event+t
		*/
		while( get_nloci_befT(  time_event+t ) > nloci_befT )
		{
			// extract the first time between time_event and time_event+tmp_time
			// set time_event at that time
			// update ni and rates
			// and then redraw with the new rate
		
			time_event = get_next_serial_sample( time_event );
			
			for (i=0; i<this->pop->get_isolation_ns() ; i++)     
				ni[i] =  this->get_nloci_species( i, time_event );
			
			update_rates( ni, &coal_rate, &migr_rate, &total_coal_rate, &total_migr_rate, this->pop );

			t = R.exponential_dev( 1.0/(total_coal_rate+total_migr_rate) );

			nloci_befT = get_nloci_befT( time_event );
			
		}
	
		time_event += t;

		
				
		//cerr <<"time_event is "<<time_event<<" and there are nlineages "<<  nloci_befT <<"\n";

		double dum;
		
		dum = R.uniform_dev();          // draw a dum [0,1[ 
		

		if( dum < total_migr_rate/(total_coal_rate+total_migr_rate) ){       // a migration event
		
			//cerr << "a migration event\n";
		
			
			locus *rand_locus;                    // a random locus when the species us selected for migration
			double sum_P = 0;                     // a tmp value to pick which locus is migrating.
			
			/*
				choose the from species
			*/
			
			_chosen_species_=0;
			
			//cerr << "Sum_P is "<<sum_P<<" vs dum "<< dum << " migr " << migr_rate[0] << " , " <<migr_rate[1] <<"\n";

			
			while(  dum >= sum_P+migr_rate[_chosen_species_]/(total_coal_rate+total_migr_rate)){
			
				sum_P+=migr_rate[_chosen_species_]/(total_coal_rate+total_migr_rate);
				
				
				if (sum_P >total_migr_rate/(total_coal_rate+total_migr_rate)){
					cerr << "sample::structured_coalesce_sample: impossible Prob in migration. Look for a bug, bye\n";
					exit(5);
				}
				
				 _chosen_species_++;   	                               // migrate from this _chosen_species_
			}

			//cerr << "  -> species "<<_chosen_species_ <<"\n";
						
			
			rand_locus = get_random_locus_species( _chosen_species_ , time_event);     // pick a locus in the selected species
			
			_destination_species_ = get_random_species_bysize();                     // pick a destination species -- based on their relative size

			rand_locus->set_species( _destination_species_ );              // migrate (change species tag)
			
			if(opt_verbose > 1 && _chosen_species_ != _destination_species_ ){
				cerr << "Migration from species "<<_chosen_species_ <<" to species "<<_destination_species_ << " at time: "<< time_event <<"\n";
			
			}
			

			/*
				Update ni and the rates
			*/
			
			ni[ _destination_species_ ] ++;
			ni[ _chosen_species_ ] --;
			
			update_rates( ni, &coal_rate, &migr_rate, &total_coal_rate, &total_migr_rate, this->pop );
						
			continue;                                                       // redraw another event (until a coal is drawn)

		}
		else{   // a coalescent event
			
			//cerr << "a coalescent event\n";

			double sum_P = total_migr_rate/(total_coal_rate+total_migr_rate);

			_chosen_species_=0;
			while(  dum >= sum_P+coal_rate[_chosen_species_]/(total_coal_rate+total_migr_rate) ){
							
				sum_P+=coal_rate[_chosen_species_]/(total_coal_rate+total_migr_rate);

				if (sum_P >1){
					cerr << "sample::structured_coalesce_sample: impossible Prob in coalescence. Look for a bug, bye\n";
					exit(5);
				}
				
				_chosen_species_++;                                     //  coalesce within this _chosen_species_
			}
			
			//cerr << "  -> species "<<_chosen_species_ <<"\n";

			if(opt_verbose > 1  ){
				cerr << "Coalescence within species "<<_chosen_species_ <<" at time: "<< time_event <<"\n";
			
			}
		
			_coal_event_ = true;                                           //  since it is a coal event, stop the process
			
			break;
		}
	
	}
	
	/*
		At this point, we HAVE encountered a coalescent event
		within species = _chosen_species_.
	*/
	
	delete coal_rate;
	delete migr_rate;
	delete ni;
	
	int nl1, nl2;                         // two random loci

	R.uniform_2DiffInt( 0, get_nloci_species(_chosen_species_ , time_event )-1, nl1, nl2);
	coalesce2loci_species( nl1, nl2, time_event, _chosen_species_);

	//cerr << "--> done\n";
	
//	cerr<<"--> Coalesced nodes "<< l1->get_id() <<" & "<<l2->get_id() <<" at time " <<time_event<<"\n";
		
}









/*
	How many mutations between locus l1 and locus l2
	Mutations in each nodes are mutation between the node
	and the anc. node.
	as a consequence, do not count the common node between l1 and l2
*/
int sample::compute_K( int l1, int l2 ){


	int ncommon=0;               // how much they have in common. Shoud be >=1 (1: only the root)
	int Mutations;               // How many mutations between both nodes.
	int i;                       // counter


	if(! this->loci[l1]->get_n_idpath() || !this->loci[l2]->get_n_idpath())
	{
		cerr << "sample::compute_K Compute idpath before using this function, bye\n";
		exit(4);
	}

	while( loci[l1]->get_idpath()[ loci[l1]->get_n_idpath()-1-ncommon ]  ==  loci[l2]->get_idpath()[ loci[l2]->get_n_idpath()-1-ncommon ] )
		ncommon++;
		
	locus * ploc=this->loci[l1];
	Mutations=0;
	
	for( i=0; i<loci[l1]->get_n_idpath()-ncommon ; i++ ){   // this one includes the common node
		Mutations += ploc->get_newSites();
		ploc = ploc->get_ancestor();
	}

	ploc=this->loci[l2];
	
	for( i=0; i<loci[l2]->get_n_idpath()-ncommon ; i++ ){  // this one exclude the common node -- do not count it twice
		Mutations += ploc->get_newSites();
		ploc = ploc->get_ancestor();
	}
	
	return Mutations;

}

/*
	Distribution of all distance
*/
int * sample::compute_Kdistrib( void ){

	int *K;
	int p=0;
	
	K = new int[((this->nloci-1)*this->nloci)/2];
	
	for(int i=0; i<this->nloci ; i++)
		this->loci[i]->compute_idpath();
	
	for(int i=0 ; i<this->nloci-1 ; i++)
		for(int j=i+1 ; j<this->nloci ; j++)
			K[p++] = this->compute_K( i, j);

	return K;

}

/*
	Distribution of all distance
*/
int ** sample::compute_Kmatrix( void ){

	int **K;
	
	K = new int * [this->nloci];
	
	for(int i=0;i<this->nloci; i++)
		K[i] = new int [this->nloci];
	
	for(int i=0; i<this->nloci ; i++)
		this->loci[i]->compute_idpath();

	for(int i=0; i<this->nloci ; i++)
		K[i][i]=0;
	
	for(int i=0 ; i<this->nloci-1 ; i++)
		for(int j=i+1 ; j<this->nloci ; j++)
			K[i][j]=K[j][i] = this->compute_K( i, j);

	return K;

}



/*
               Funtcions to manipulate TE, locus within locus

*/	     

void sample::set_locus_TE( locus *l, int nTE, int nTE_anc ){

	locus_TE *new_TE;
	int i;
	
	l->set_size_array_TE(2*nTE - nTE_anc);
	l->set_n_TE ( nTE );
	
	new_TE = new locus_TE[ l->get_size_array_TE() ];
	
	cerr << "create a TE array of size "<<l->get_size_array_TE()<<"\n";
	
	if( ! new_TE ){
		cerr << "sample::set_locus_TE: cannot allocate new_TE\n";
		exit(3);
	}
	l->set_array_TE ( new_TE );          // it contains the te in the mrca plus their ancestor
	

	for(i=0;i< l->get_size_array_TE() ;i++)
		l->get_array_TE()[i].set_locus(l->get_id());        // they are all connected to their locus

	for(i=0;i<l->get_n_TE();  i++)
		l->get_array_TE()[i].set_time(l->get_time());       // the first ones 
}



static void print_connection( locus_TE *x ){

	printf("Locus_TE of id %d has anc: %d ; desc1: %d ; desc2: %d\n",
	       x->get_id(),
	       x->get_ancestor()? x->get_ancestor()->get_id():-1,
	       x->get_descendant1()?x->get_descendant1()->get_id():-1,
	       x->get_descendant2()?x->get_descendant2()->get_id():-1
	       );

}
void sample::print_all_TE_connections( void ){

//	cerr << "++ sample id : " << get_id() << "\n";

	if(descendant == NULL){
	
		for(int i=0 ; i<get_nloci() ; i++ )
			for(int j=0 ; j < loci[i]->get_size_array_TE() ; j++ )
				print_connection( loci[i]->get_array_TE()+j );
				
	}
	else
		for(int j=0 ; j<loci[0]->get_size_array_TE() ; j++ )
			print_connection( loci[0]->get_array_TE()+j );        
	
}

void sample::coalesce_TE_locus( double time, int *alive_locus, int n_alive, locus_TE *array_TE, int size_array_TE, int n_te_anc ){

	
	/*
		Draw locus
	*/
	int n1, n2;
	extern Random R;                                               // an object to draw random variables declared in random.h


//	cerr << "draw numbers between 0, "<<n_alive-1<<"\n";
	
	R.uniform_2DiffInt(0, n_alive-1, n1, n2);
	
//	cerr <<"anc is the elt number: "<<size_array_TE <<"\n";
	
	locus_TE *panc  =  array_TE +  size_array_TE + n_te_anc - n_alive ;
	locus_TE *l1    =  array_TE +  alive_locus[ n1 ] ;
	locus_TE *l2    =  array_TE +  alive_locus[ n2 ] ;

//	print_connection( panc );
//	print_connection( l1 );
//	print_connection( l2 );
	
	
//	cerr << "sample::coalesce_TE_locus: coalesce TE "<<l1->get_id()<<" with "<<l2->get_id()<<" at time "<<time <<" to anc "<<panc->get_id()<<"\n";
			
	/*
		Coalesce at the locus level
	*/

	panc->reset_connections( );
	panc->set_time(time);                                              // this new node started at time_event

	l1->connect_ancestor( (locus *)panc, 1 );                                    // and is connected to l1 and l2
	l2->connect_ancestor(  (locus *)panc, 2 );
	
	/*
	print_connection( panc );
	print_connection( l1 );
	print_connection( l2 );
	*/
	
		
	/*
		update alive_locus
	*/
	if(n1>n2){
		int tmp=n2; n2=n1; n1=tmp;
	}
	
	alive_locus[ n1 ] = size_array_TE + n_te_anc - n_alive;
	alive_locus[ n2 ] = alive_locus[ n_alive-1 ];
	

//	cerr <<"nalive is now :"<< n_alive<<"\n";
		
//	for( i=0; i<n_alive-1; i++)
//		cerr << "alive_locus[ "<<i<<" ]: "<<alive_locus[ i ]<<"\n";

}


void sample::coalesce_TE_mrca_BottomUp( float dupl_rate, int n_te_mrca, locus *mrca ){                      // coalesce its transposons if mrca


	if(mrca->get_ancestor() != NULL){
		cerr <<"sample::coalesce_TE_mrca: this is not the ancestral sample, bye\n";
		exit(1);
	}
	
	
	set_locus_TE( mrca , n_te_mrca, 1 );
	
	int n_alive = mrca->get_n_TE();
	extern Random R;                                               // an object to draw random variables declared in random.h
	
	int i;
	double time_event = mrca->get_time();

	int * alive_array = new int[n_alive];
	
	
	
	if( alive_array ==  NULL ){
		cerr << "locus::coalesce_transposons_mrca: cannot allocate alive_array, bye\n";
		exit(3);
	}
	
	for( i=0;i<n_alive;i++)
		alive_array[ i ] = i;
	
//	cerr <<"BEGIN -- time is "<<time_event<<"\n";
//	cerr <<"BEGIN -- nalive is now :"<< n_alive<<"\n";
		
/*	for( i=0; i<n_alive; i++){
		cerr << "alive_array[ "<<i<<" ]: "<<alive_array[ i ]<<"\n";
	}
*/

	
	while( n_alive > 1){
	
		double rate = (double) ( n_alive-1 )*dupl_rate;

		time_event += R.exponential_dev( 1.0/rate );

//		cerr <<"sample::coalesce_TE_mrca_BottomUp: time is "<<time_event<<"\n";
		
		coalesce_TE_locus( time_event, alive_array, n_alive, mrca->get_array_TE(), mrca->get_size_array_TE(), 1 );
	
		n_alive --;
		
	}

//	cerr << "stop the coalescence\n";
	
	delete alive_array;
}

void sample::coalesce_1st_TE_locus_TopDown( float dupl_rate, float time_1st_TE, locus *l  ){    // built transposon tree from a locus l to THE ancestor (that conditioned it)

	
	double time_next_event = time_1st_TE;
	
	double rate;
	extern Random R;                                               // an object to draw random variables declared in random.h
	
	vector<double> duplication_times;                                     // I do not know if I need to check this worked properly ??
	
	/*
		do the TE expansion forward
	*/
	do{
	
		rate = (double) ( 1 + duplication_times.size() ) * dupl_rate;
		time_next_event -= R.exponential_dev( 1.0/rate );
		duplication_times.push_back( time_next_event );
	
	}while( time_next_event >  l->get_time() );

	duplication_times.pop_back( );

	cerr <<"sample::coalesce_1st_TE_locus_TopDown:  there was "<<duplication_times.size()<<" duplication from 1st TE to locus "<< l->get_id() <<"\n";
	
		
	int n_alive = 1 + (int) duplication_times.size() ;
	int * alive_array = new int[ n_alive ];
	if( alive_array ==  NULL ){
		cerr << "locus::coalesce_transposons_mrca: cannot allocate alive_array, bye\n";
		exit(3);
	}

	
	this->set_locus_TE( l, n_alive,  1 );

	for(int i=0;i<n_alive;i++)
		alive_array[ i ] = i;

	/*
		do the coalescent backward
	*/
	while( n_alive > 1 ){
	
		
		cerr << "now remains "<<n_alive<<" TE\n";
		
		coalesce_TE_locus( duplication_times.back(), alive_array, n_alive, l->get_array_TE(), l->get_size_array_TE(), 1 );
		duplication_times.pop_back( );
		n_alive --;
	
	}
	
	
	delete alive_array;
}




void sample::coalesce_TE_locus_TopDown( float dupl_rate, locus *l  ){         // built transposon tree from a locus l to its ancestor (that conditioned it)


	int nTE_anc = l->get_ancestor()->get_n_TE();
	
	if(nTE_anc == 0)
		return;
	
	
	double time_next_event = l->get_ancestor()->get_time();
	
	double rate;
	extern Random R;                                               // an object to draw random variables declared in random.h
	
	vector<double> duplication_times;                                     // I do not know if I need to check this worked properly ??
	
	/*
		do the TE expansion forward
	*/
	do{
	
		rate = (double) ( nTE_anc + duplication_times.size() )*dupl_rate;
		time_next_event -= R.exponential_dev( 1.0/rate );
		
//		cerr <<"sample::coalesce_TE_locus_TopDown:  time is "<<time_next_event<<"\n";

		duplication_times.push_back( time_next_event );
	
	}while( time_next_event >  l->get_time() );
	
	
	duplication_times.pop_back( );
	
	cerr <<"sample::coalesce_TE_locus_TopDown:  there was "<<duplication_times.size()<<" duplication from locus "<< l->get_ancestor()->get_id() << "to locus "<< l->get_id() <<"\n";
	
	
	this->set_locus_TE( l, (int)  nTE_anc + duplication_times.size(),  nTE_anc );
	
	int n_alive = nTE_anc + (int) duplication_times.size() ;
	int * alive_array = new int[ n_alive ];
	if( alive_array ==  NULL ){
		cerr << "locus::coalesce_transposons_mrca: cannot allocate alive_array, bye\n";
		exit(3);
	}


	for(int i=0;i<n_alive;i++)
		alive_array[ i ] = i;

	/*
	for(int i=0; i<n_alive; i++)
		cerr << "alive_array[ "<<i<<" ]: "<<alive_array[ i ]<<" locus_TE of id "<< l->get_array_TE()[ alive_array[i] ].get_id() <<"\n";
	cerr<<"\n\n";
	*/

	/*
		do the coalescent backward
	*/
	while( n_alive > nTE_anc ){
	
		
		coalesce_TE_locus( duplication_times.back(), alive_array, n_alive, l->get_array_TE(), l->get_size_array_TE(), nTE_anc );
		
		duplication_times.pop_back( );
		
		n_alive --;
	
		/*
		for(int i=0; i<n_alive; i++)
			cerr << "alive_array[ "<<i<<" ]: "<<alive_array[ i ]<<" locus_TE of id "<< l->get_array_TE()[ alive_array[i] ].get_id() <<"\n";
		cerr<<"\n\n";
		*/
	}
	
	/*
		Connect all remaining TE to the ancestral TE
	*/
		
		
	/*
	for(int i=0; i<n_alive; i++){
		cerr << "alive_array[ "<<i<<" ]: "<<alive_array[ i ]<<" locus_TE of id "<< l->get_array_TE()[ alive_array[i] ].get_id() <<"\n";
	}
	cerr<<"\n\n";
	*/

	for(int i=0; i<n_alive; i++){
			
		locus *desc=l->get_array_TE()+ alive_array[i] ;
	
		if( l->get_ancestor()->get_array_TE()[i].get_descendant1() == NULL  ){
			desc->connect_ancestor( (locus *) (l->get_ancestor()->get_array_TE()+i) , 1 );
		}
		else{
			desc->connect_ancestor( (locus *) (l->get_ancestor()->get_array_TE()+i) , 2 );
		}
		
/*
		print_connection( l->get_array_TE()+ alive_array[i]  );
		print_connection( l->get_ancestor()->get_array_TE()+i );		
*/		
	}
	
	
	//exit(1);
	
	delete alive_array;
}


void sample::coalesce_TE_sample_TopDown( float dupl_rate  ){


	if( this->ancestor == NULL){
	
		cerr << "this function cannot be used for the MRCA sample\n";
		exit(1);
	
	}
	
	if( this->descendant == NULL){                                  // this is the youngest sample --do it on all leaves--
	
		for(int i =0; i< get_nloci(); i++)
			coalesce_TE_locus_TopDown( dupl_rate, loci[i]  );
	
	}else
		coalesce_TE_locus_TopDown( dupl_rate, loci[0]  );


}






