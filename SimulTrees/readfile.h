/***
		file:     readfile.h
		function: header for readfile functions -- from file to memory
		author:   <gachaz>
		date:     nov 08
***/


#ifndef _READFILE_H_
#define _READFILE_H_

#include <stdlib.h>


/***
	The locus class
	each locus is a node in a tree --
	a locus is a string sharing an identical history
***/

class readfile{


	private:
	
		double (* CI_S)[2];  /* an array of pairs -- CI for a given S */
		int size_CI_S;      /* the size of this array */
	
		double *w1;        /* an array of weight for the w1 vector */
		int size_w1;      /* the size of this array */
	
		double *w2;        /* an array of weight for the w2 vector */
		int size_w2;      /* the size of this array */
	
	public:
	
		readfile(void);
		~readfile(void);
	
		void readfile_CI_S( char *file );
		void readfile_w( char *file, char w_num );
		
		int get_size_CI_S( void ){ return size_CI_S; };
		double * get_CI_S( int x ){ return CI_S[x]; };

		int get_size_w1( void ){ return size_w1; };
		double * get_w1( ){ return w1; };

		int get_size_w2( void ){ return size_w2; };
		double * get_w2( ){ return w2; };

};

#endif
