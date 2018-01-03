/***
		file:     locus_transposon.cpp
		function: loci functions
		author:   <amikezor>
		date:     feb 04
		modif:    mar 09 -- add the id_path
***/


#include <string.h>

#include <iostream>
using namespace std;

#include "locus_TE.h"
#include "random.h"

int locus_TE::Count_All_transposon=0;

/*
	Constructor
*/
locus_TE::locus_TE( void ){
	
	id = Count_All_transposon++;
	time=0;
	
	//cerr << "<-- create locus_TE "<<this->id<<" at time "<<this->time<<"\n";
}


/*
	Destructor
*/
locus_TE::~locus_TE( void  ){
		
	Count_All_transposon--;

	//cerr << "<-- delete locus_TE "<<this->id<<" at time "<<this->time<<"\n";

}
