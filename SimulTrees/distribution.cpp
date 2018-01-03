
#include "distribution.h"
#include <stdlib.h>
#include <math.h>

#include <iostream>
using namespace std;



/***
	sort two array by the order of the first one
	Qsort is implemented following the algorithm explained in Numrec C, pp 332-335
	Under a certain size, a part of the array is sorted by the naive method (straight insertion)
***/

#define NSTACK 1000                           /* the number of subarray that could be memorized */
#define SWITCH_METHOD 7                       /* under that size switch to a Insertion method */
#define SWAP(a,b) {tmp=(a); (a)=(b); (b)=tmp;}

static void qsort2(double *array1, double *array2, long n){

	
	long *Stack;       /* where remaining values of subarrays are temporarily stocked */
	long nStack=0;     /* How many of these values are in the Stack. Is never odd */
	
	long end=n-1,      /* the last value of the array to sort */
	     beg=0,        /* the first value ---  */
	     postbeg;      /* the second value ---  */

	double val_postbeg;  /* the value of the postbeg position - the one which is used for partion-exchange */

	long demi_len;     /* half of the distance between beg and end */
	
	long i,            /* counter from postbeg to right */
	     j;            /* counter from end to left */

	double    val1_i,       /* for insertion stock temporarily a value */ 
		  val2_i;
	
	double tmp;          /* used for the SWAP macro */
	
	
	if ( (Stack = (long *)malloc(NSTACK*sizeof(long))) == NULL )
	fprintf(stderr, "qsort2: not enough memory for Stack, bye\n"), exit(3); 
		
	while(1){ 

		if( end-beg+1 > SWITCH_METHOD){
	
			demi_len = (end-beg) >> 1 ; 
			postbeg = beg+1;
			
			SWAP( array1[beg+demi_len], array1[postbeg] );
			SWAP( array2[beg+demi_len], array2[postbeg] );
			
			if(array1[beg] > array1[postbeg]){            /* rearrange to have  beg <= postbeg <= end */ 
				SWAP( array1[beg], array1[postbeg] );
				SWAP( array2[beg], array2[postbeg] );
			}

			if(array1[beg] > array1[end]){
				SWAP( array1[beg], array1[end] );
				SWAP( array2[beg], array2[end] );
			}

			if(array1[postbeg] > array1[end]){
				SWAP( array1[postbeg], array1[end] );
				SWAP( array2[postbeg], array2[end] );
			}
			
			
			i = postbeg;
			j = end;
						
			val_postbeg =  array1[postbeg];
			
			while(1)                                   /* enter the partition exchange process */
				{
					do i++; while( array1[i] < val_postbeg );
					do j--; while( array1[j] > val_postbeg );
					
					if(j<i) break;
					
					SWAP( array1[i], array1[j] );
					SWAP( array2[i], array2[j] );
				}
						
			SWAP( array1[postbeg] , array1[j] );   /* place the postbeg value into j */
			SWAP( array2[postbeg] , array2[j] );
			
			if(nStack+2 > NSTACK)
				fprintf(stderr, "qsort2: not enough Stack... sorry bye\n"),exit(1);
			
			if(end-i+1 >= j-beg){
				Stack[nStack++] = i;                /* stock in Stack the largest and go with the smallest */
				Stack[nStack++] = end;
				
				end = j-1;
			}
			else{
				Stack[nStack++] = beg;
				Stack[nStack++] = j-1;
				
				beg = i;
			}
			
		}
		else{                                      /* Under a certain size, switch to the straight insertion method */
	
			
			for(i=beg+1; i <= end ; i++){
			
				val1_i = array1[i];
				val2_i = array2[i];
				
				for(j=i-1;j>=beg;j--){
				
					if(array1[j] < val1_i)break;
					
					array1[j+1] = array1[j];
					array2[j+1] = array2[j];
				}
				array1[j+1]=val1_i;
				array2[j+1]=val2_i;
			}

			if(nStack==0)break;          /* this si the end  - exit the process */
			
			end = Stack[--nStack];       /* go for the next round with the stacked parameters */
			beg = Stack[--nStack];
		}

	 }

	free(Stack);

}

/*	
	Sort many arrays as a function of a first array
	n : arrays size
	narrays: how many arrays
	X : what is the sorting array
*/
static void qsortn(double **array, long n, long narrays, int X){

	
	long *Stack;       /* where remaining values of subarrays are temporarily stocked */
	long nStack=0;     /* How many of these values are in the Stack. Is never odd */
	
	long end=n-1,      /* the last value of the array to sort */
	     beg=0,        /* the first value ---  */
	     postbeg;      /* the second value ---  */

	double val_postbeg;  /* the value of the postbeg position - the one which is used for partion-exchange */

	long demi_len;     /* half of the distance between beg and end */
	
	long i,            /* counter from postbeg to right */
	     j;            /* counter from end to left */

	double    *val_i;       /* for insertion stock temporarily a value */ 
	
	double tmp;          /* used for the SWAP macro */
	
	

	if ( (val_i = (double *)malloc(narrays*sizeof(double))) == NULL )
	fprintf(stderr, "qsort2: not enough memory for val_i, bye\n"), exit(3); 
	

	if ( (Stack = (long *)malloc(NSTACK*sizeof(long))) == NULL )
	fprintf(stderr, "qsort2: not enough memory for Stack, bye\n"), exit(3); 
		
	while(1){ 

		if( end-beg+1 > SWITCH_METHOD){
	
			demi_len = (end-beg) >> 1 ; 
			postbeg = beg+1;
			
			
			for(int a=0; a<narrays; a++)
				SWAP( array[a][beg+demi_len], array[a][postbeg] );
			
			if(array[X][beg] > array[X][postbeg]){            /* rearrange to have  beg <= postbeg <= end */ 

				for(int a=0; a<narrays; a++)
					SWAP( array[a][beg], array[a][postbeg] );
			}

			if(array[X][beg] > array[X][end]){
			
				for(int a=0; a<narrays; a++)
					SWAP( array[a][beg], array[a][end] );
			
			}

			if(array[X][postbeg] > array[X][end]){
				
				for(int a=0; a<narrays; a++)
					SWAP( array[a][postbeg], array[a][end] );

			}
			
			
			i = postbeg;
			j = end;
						
			val_postbeg =  array[X][postbeg];
			
			while(1)                                   /* enter the partition exchange process */
				{
					do i++; while( array[X][i] < val_postbeg );
					do j--; while( array[X][j] > val_postbeg );
					
					if(j<i) break;
					
					for(int a=0; a<narrays; a++)
						SWAP( array[a][i], array[a][j] );
					
				}
						
			for(int a=0; a<narrays; a++)
				SWAP( array[a][postbeg], array[a][j] );
			
			if(nStack+2 > NSTACK)
				fprintf(stderr, "qsort2: not enough Stack... sorry bye\n"),exit(1);
			
			if(end-i+1 >= j-beg){
				Stack[nStack++] = i;                /* stock in Stack the largest and go with the smallest */
				Stack[nStack++] = end;
				
				end = j-1;
			}
			else{
				Stack[nStack++] = beg;
				Stack[nStack++] = j-1;
				
				beg = i;
			}
			
		}
		else{                                      /* Under a certain size, switch to the straight insertion method */
	
			
			for(i=beg+1; i <= end ; i++){
			
				for(int a=0; a<narrays; a++)
					val_i[a] = array[a][i];
				
				for(j=i-1;j>=beg;j--){
				
					if(array[X][j] < val_i[X])break;
					
					for(int a=0; a<narrays; a++)
						array[a][j+1] = array[a][j];
				}
				for(int a=0; a<narrays; a++)
					array[a][j+1] = val_i[a];
						
			}

			if(nStack==0)break;          /* this si the end  - exit the process */
			
			end = Stack[--nStack];       /* go for the next round with the stacked parameters */
			beg = Stack[--nStack];
		}

	 }

	free(Stack);
	free(val_i);

}




#define SIGN( a ) ( ( (a) > 0 )?1: ( ((a)==0)?0:-1)  )
static int Increase(const void *v1, const void *v2){  	return (int)SIGN( *((double *)v1) - *((double *)v2));  };
#undef SIGN


distribution::distribution(int rep, float initCI, bool full, bool weight){


	if( weight == true )
		full = true;      /* do not know how to handle weights in incomplete distribution */

	this->CI=initCI;
	this->rep=rep;
	this->full=full;

	if(this->full == false){
		
		array_size = (int) (rep*(1.00-initCI));       /* only store a small amount of the distrib -- the extrema -- */
		if(array_size%2 == 1)
			array_size++;
	}
	else
		array_size = rep;
	
	array = (double *)calloc( array_size, sizeof(double) );
	if(!array)fprintf(stderr, "cannot allocate array of number \n"), exit(1);

	if(weight == true){
		this->w = (double *)calloc( array_size, sizeof(double) );
		if(! this->w )fprintf(stderr, "cannot allocate array of number \n"), exit(1);
	}
	else
		this->w = NULL;

	this->sum=0;
	this->sum_sq=0;
	this->sum_w=0;
	this->c=0;
	
}

distribution::~distribution(void){

	free(array);
	if(w)free(w);
}


/*
	If weight is unused, set it to -1
*/
void  distribution::add_value( double x, double weight ){



	if( isinf(x) )
		cerr << "Watch out ! There is an infinite number in the distribution\n";

	if(weight >= 0){
		sum    += weight*x;
		sum_w  += weight;
		sum_sq += weight*x*x;
	}
	else{
		sum    += x;
		sum_sq += x*x;
		sum_w  ++;
	}

	if( this->full == false ){
	
		if(c<array_size){
			array[c] = x;
		}
		else{
			if(c==array_size)
				qsort( (void *)array, (size_t) c, (size_t) sizeof(double), Increase );
			
			//for(int i=0;i<array_size;i++)
			//	printf("%f ",array[i]);
			//printf("\n");
		
			int a = array_size/2;			               // set it to the first large values (right part)

			if(  x > array[a]  ){
		
			
				while( a+1 < array_size &&  x > array[a+1] ){    // scan for the place to insert x on the right side
				
					array[a]=array[a+1];
					a++;
					
				}
				array[a]=x;
			}
			else if( a-- && x < array[a] ){                        // set it to the first small values (left part)
						
				while( a>0 && x < array[a-1] ){             // scan for the place to insert x on the left side
				
					array[a]=array[a-1];
					a--;
					
				}
				array[a]=x;
			
			}
		
		}
	
	}else{
	
		if(c >= array_size){
			cerr << "Cannot add further values to the distrib, sorry\n";
			return;
		}


		this->array[c] = x;
		if(weight >= 0 )
			this->w[c] = weight;

		if(c == array_size-1){
			if(this->w)
				qsort2(this->array, this->w, array_size);
			else
				qsortn(&(this->array), array_size, 1, 0);
		}
		

	}
		
	c++;
}


void  distribution::clean_distribution( void ) {

	this->sum=0;
	this->sum_sq=0;
	this->c=0;
	this->sum_w=0;

}


double distribution::get_low_boundary( void ){

	int index=-1;

	if( this->full == false )
		index = array_size/2-1;
	else
		if( this->w == NULL)
			index = (int) floor(  array_size * (1.0-this->CI)/2.0 ) -1;
	else
	{
		float w_t=w[0];
		int c=0;
		while( w_t/sum_w < (1.00-this->CI)/2.0 ){
			w_t += w[++c];
		}
		
		index=c-1;
	}
		
	if( index < 0){
		//cerr << "the CI is too small for your distribution, bye\n";
		//exit(1);
		index++;
	}
	
	return this->array[index];


}

double distribution::get_high_boundary( void ){
	
	int index=-1;

	if( this->full == false )
		index = array_size/2;
	else if( this->w == NULL)
		index = floor(  array_size * (1.0-((1.0-this->CI)/2.0)) ) -1;
	else
	{
		float w_t=w[array_size-1];
		int c=array_size-1;
		
		while( w_t/sum_w < (1.00-this->CI)/2.0 )
			w_t += w[--c];
		
		index=c+1;
	}
		
	if( index >= array_size){
		 index--;
		//cerr << "the CI is too small for your distribution, bye\n";
		//exit(1);
	}

	return this->array[index];


}

void distribution::print_distribution( void ){
	
	double XXX=0;
	int i;
	
	for(i=0;i<array_size; i++){
	
		if(this->w){
			XXX += w[i]/sum_w;
			printf("%d %.3f %.5f %.5f %.5f %.5f\n", i, array[i], w[i], w[i]/sum_w, XXX, 1-XXX+w[i]/sum_w );
		}
		else
			printf("%d %.3f\n", i, array[i]);
	}
	
}
