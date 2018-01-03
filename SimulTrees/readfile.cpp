/***
		file:     readfile.cpp
		function: functions that readfile and get them into memory
		author:   <gachaz>
		date:     nov 08
***/


#include <stdio.h>
#include "readfile.h"

readfile::readfile( void ){

	this->CI_S = NULL;
	this->size_CI_S = 0;
	
	this->w1 = NULL;
	this->size_w1 = 0;
}


readfile::~readfile( void ){

	if( CI_S != NULL)
		free(CI_S);
	
	if( w1 != NULL)
		free(w1);
	
}

void readfile::readfile_w( char *file, char w_num ){


	FILE *f_in;
	double w;
	int memory_size=30;
	
	double *w_ptr;
	int size_w;

	if(w_num != 1 && w_num != 2)
		fprintf(stderr, "readfile_w: w_num has to be 1 or 2\n"), exit(6);   

	w_ptr = (double *)malloc( (size_t)  memory_size*sizeof( double  ) );
	if(! w_ptr )fprintf(stderr, "readfile_w: cannot allocate w_ptr, bye\n"), exit(3);
     
	f_in = fopen(file, "r");
     
	if( !f_in )
		fprintf(stderr, "readfile_w: cannot open file %s, check it !\n", file), exit(2);


	/*
		then extract the three values per line S, Ci_low, CI_up
	*/
	size_w=0;
	while( fscanf(f_in, "%lf", &w) == 1){
     
		size_w++;

		if( size_w > memory_size){
	     
			memory_size*=2;
		     
			w_ptr = ( double *) realloc( (void *)w_ptr, (size_t)  memory_size*sizeof(double) );
			if(! w_ptr )fprintf(stderr, "readfile_w: cannot reallocate w_ptr, bye\n"), exit(3);
		     
		}
	     
		w_ptr[size_w-1]=w;
     
	}
     
	w_ptr = ( double *) realloc( (void *)w_ptr, (size_t)  size_w*sizeof(double) );
	if(! w_ptr )fprintf(stderr, "readfile_w: cannot reallocate w_ptr, bye\n"), exit(3);


	fclose(f_in);

	if(w_num==1){
		this->w1=w_ptr;
		this->size_w1=size_w;
	}
	else{
		this->w2=w_ptr;
		this->size_w2=size_w;
	}

}



void readfile::readfile_CI_S( char *file ){


     FILE *f_in;
     int S;
     double CI_low,
	   CI_up;
	   
     int memory_size=10;
     
     
     this->CI_S = (double (*)[2]) malloc( (size_t)  memory_size*sizeof( double [2] ) );
     
     if(! this->CI_S )fprintf(stderr, "readfile_CI_S: cannot allocate this->CI_S, bye\n"), exit(3);
   	    
     this->CI_S[0][0]=0;
     this->CI_S[0][1]=0;
     this->size_CI_S=1;
     
     
     f_in = fopen(file, "r");
     
     if( !f_in )
	     fprintf(stderr, "readfile_CI_S: cannot open file %s, check it !\n", file), exit(2);


     while( fgetc(f_in) != '\n')		 /* while the first letter is not '\n' */
	     while( fgetc(f_in) != '\n');	 /* remove the entire line */


     /*
	     then extract the three values per line S, Ci_low, CI_up
     */
     
     while( fscanf(f_in, "%d%lf%lf", &S, &CI_low, &CI_up) == 3){
     
	     this->size_CI_S++;

	     if( this->size_CI_S >= memory_size){
	     
		     memory_size*=2;
		     
		     this->CI_S = ( double (*)[2]) realloc( (void *)this->CI_S, (size_t)  memory_size*sizeof(double [2]) );
		     if(! this->CI_S )fprintf(stderr, "readfile_CI_S: cannot reallocate this->CI_S, bye\n"), exit(3);
		     
	     }
	     
	     this->CI_S[S][0] = CI_low;
	     this->CI_S[S][1] = CI_up;

     
     
     }
     
     this->CI_S = ( double (*)[2]) realloc( (void *)this->CI_S, (size_t)  size_CI_S*sizeof(double [2]) );
     if(! this->CI_S )fprintf(stderr, "readfile_CI_S: cannot reallocate this->CI_S, bye\n"), exit(3);


//     for(int i=0; i < this->size_CI_S; i++)
//     	printf("%d %f %f\n", i, this->CI_S[i][0],this->CI_S[i][1] );
     	

     fclose(f_in);
}
