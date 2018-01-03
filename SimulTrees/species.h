/***
		file:     species.h
		function: header of a species.cpp - handle species tree where lineage can exist
		author:   <amikezor>
		date:     mar 09
***/

#ifndef _SPECIES_H_
#define _SPECIES_H_


class species {

	private:
		int id;             // species id

		double time;        // the bckwd time when this species ended its existence
		double N;           // the species # of individuals
		double mu;          // the sp mutation rate

		species *anc;       // point to anc species.
		species *desc1;     // point to desc1 species
		species *desc2;     // point to desc2 species.
		
	public:
	
		species( double time=0, double N=0, double mu=0 );
		~species( void );
		
		/*
			get functions
		*/
		int get_id( void ){return this->id; };
		double get_time( void ){return this->time;};
		double get_N( void ){ return this->N;};
		double get_mu( void ){ return this->mu;};
		get_Theta( void ){return 2*this->N*this->mu;};
		
		species * get_desc1(void){return this->desc1;};
		species * get_desc2(void){return this->desc2;};
		species * get_anc(void){return this->anc;};
		
		
		/*
			Set functions
		*/
		void connect_sp( species * psp, char type );      // type is 0=anc, 1=desc1 or 2=desc2
		

};


#endif
