/***
		file:     species.cpp
		function: all functions associated with the species, a tree in which loci can evolve
		author:   <amikezor>
		date:     mar 09,
		modif :   
***/



int sample::Count_All_Species=0;            // number of sample objets in use



/*
	Constructor
*/
species::species( double input_time=0, double input_N=0, double input_mu=0 ){

	this->id = Count_All_Species++;          // increment sample objects
	
	this->time=input_time;
	this->N=input_N;
	this->mu=input_mu;


	this->anc=NULL;
	this->desc1=NULL;
	this->desc2=NULL;

}


/*
	Set connections between species
	type is 0=anc0, 1=desc1 or 2=desc2
*/
void species::connect_sp( species * psp, char type ){


	switch( type ){
	
	
		case 0:
			this->anc=psp;
			break;
		case 1:
			this->desc1=psp;
			break;
		case 2:
			this->desc2=psp;
			break;
		default:
			cerr << "type is either 0-anc, 1-desc1 or 2-desc2, bye\n";
			exit(1);
		
	}
}
