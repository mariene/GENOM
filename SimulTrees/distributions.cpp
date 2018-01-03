
#include "distributions.h"
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


distributions::distributions(int rep, float initCI, bool full, bool weight, int vals){


	if( weight == true )
		full = true;      /* do not know how to handle weights in incomplete distribution */

	this->CI=initCI;
	this->rep=rep;
	this->full=full;
	this->array_values=vals;


	/*
		array
	*/
	if(this->full == false){
		
		this->array_size = (int) (rep*(1.00-initCI));       /* only store a small amount of the distrib -- the extrema -- */
		if(this->array_size%2 == 1)
			this->array_size++;
	}
	else
		this->array_size = rep;
	
	this->array = (double **)malloc( this->array_values * sizeof(double *) );
	if(!this->array)fprintf(stderr, "distributions: cannot allocate array\n"), exit(1);
	
	for(int i=0; i<this->array_values; i++){
		array[i] = (double *)calloc( this-> array_size, sizeof(double) );
		if(!array[i])fprintf(stderr, "distributions: cannot allocate array[%d]\n", i), exit(1);
	}

	/*
		weight
	*/
	if(weight == true){
	
		this->w = (double *)calloc( this->array_size , sizeof(double ) );
		if(! this->w )fprintf(stderr, "distributions: cannot allocate w\n"), exit(1);
	}
	else
		this->w = NULL;

	/*
		and sums
	*/
	this->sum_X  = (double *)calloc( this->array_values , sizeof(double) );
	this->sum_X2 = (double *)calloc( this->array_values , sizeof(double) );
	this->sum_XY = (double *)calloc( (this->array_values*(this->array_values-1))/2 , sizeof(double) );
	if( ! this->sum_X || !this->sum_X2 || !this->sum_XY )fprintf(stderr, "distributions: cannot allocate sum_X or sum_X2 or sum_XY\n"), exit(1);
	
	/*
		and Covariance matrix
	*/
	this->Covariance = (double **)malloc( this->array_values * sizeof(double *) );
	if(! this->Covariance )fprintf(stderr, "distributions: cannot allocate Covariance\n"), exit(1);

	for(int i=0; i<this->array_values; i++){
		this->Covariance[i] = (double *)calloc( this-> array_values, sizeof(double) );
		if(!this->Covariance[i])fprintf(stderr, "distributions: cannot allocate Covariance[%d]\n", i), exit(1);
	}	
	
	
	
	this->c=0;
	this->sum_w = 0;
	
}

distributions::~distributions(void){

	for( int i=0; i<this->array_values; i++ )
		free( array[i] );
	
	free(array);
	free(sum_X);
	free(sum_X2);
	free(sum_XY);
	
	if(w) free(w);
}

double distributions::get_variance( int i ){

	double mean_w = sum_w/c;

	return (sum_w/(sum_w-mean_w))*(sum_X2[i]/(double)sum_w - get_mean(i)*get_mean(i));
};


double distributions::get_covariance( int i, int j ){
	
	double mean_w = sum_w/c;

	if( i == j )
		return  get_variance( i );
	
	else
	{
		if(i>j){
			int tmp=i;
			i=j;
			j=tmp;
		}
		
		int cell = (j-i-1);         // last row
		for( int a=1; a<=i; a++ ) cell += this->array_values -a; // al previous row 
	
		return (sum_w/(sum_w-mean_w))*( sum_XY[ cell ]/(double)sum_w - get_mean(i)*get_mean(j) );
	}
}


void distributions::compute_Covariances( void ){

	for( int i=0; i<this->array_values; i++)
		for(int j=0; j<=i; j++ )
			this->Covariance[i][j] = this->Covariance[j][i] = get_covariance( i, j );

}


/*
	If weight is unused, set it to -1
*/
void distributions::add_values( double *X, double weight ){

	double w;

	for( int a=0; a<this->array_values; a++ ) 
		if( isinf(X[a]) )
			cerr << "Watch out ! Adding an infinite number in the distribution "<<a<<"\n";

	if(weight >= 0)
		w = weight;
	else
		w = 1;
	
	sum_w  += w;

	for(int i=0; i<this->array_values; i++){	
		
		sum_X[i]  += w * X[i];
		sum_X2[i] +=  w * X[i] * X[i];
		
		for(int j=i+1;j<this->array_values; j++){
		
			int cell = (j-i-1);         // last row
			for( int a=1; a<=i; a++ ) cell += this->array_values -a; // al previous row 
								
			sum_XY[ cell ] += w * X[i] * X[j];

		}
	}


	if( this->full == false ){
	
		if(c<this->array_size){
		
			for(int i=0; i<this->array_values; i++)
				array[i][c] = X[i];
		}
		else{
			if(c==array_size){
			
				for(int i=0; i<this->array_values; i++)
					qsort( (void *)array[i], (size_t) c, (size_t) sizeof(double), Increase );
			}
			
			//for(int i=0;i<array_size;i++)
			//	printf("%f ",array[i]);
			//printf("\n");
		
			int a = array_size/2;			               // set it to the first large values (right part)

			for(int i=0; i<this->array_values; i++){
			
				if(  X[i] > array[i][a]  ){
		
					while( a+1 < array_size &&  X[i] > array[i][a+1] ){    // scan for the place to insert x on the right side
				
						array[i][a]=array[i][a+1];
						a++;
					
					}
					array[i][a]=X[i];
				}
				else if( a-- && X[i] < array[i][a] ){                        // set it to the first small values (left part)
						
					while( a>0 && X[i] < array[i][a-1] ){             // scan for the place to insert x on the left side
				
						array[i][a]=array[i][a-1];
						a--;
					
					}
					array[i][a]=X[i];
				
				}
			}
		}
	
	}else{
	
		if(c >= array_size){
			cerr << "Cannot add further values to the distrib, sorry\n";
			return;
		}


		for(int i=0; i<this->array_values; i++)
			this->array[i][c] = X[i];

		if(weight >= 0 )
			this->w[c] = weight;

//		if(c == array_size-1)
//			qsort2(this->array[i], this->w[i], array_size);
		

	}
		
	c++;
}


void  distributions::clean_distribution( void ) {


	for( int a=0; a < this->array_values; a++){
		this->sum_X[a]=0;
		this->sum_X2[a]=0;
	}

	for( int a=0; a < ((this->array_values-1)*this->array_values)/2; a++){
		this->sum_XY[a]=0;
	}
	
	
	this->c=0;
	this->sum_w=0;

}


double distributions::get_low_boundary( int X ){

	int index=-1;

	if( this->full == false ){
		index = array_size/2-1;
	}
	else{
	
		/*
			Sort the values of all arrays as well as the weight array by array[i]
		*/
		double **all_array = (double **)malloc( (this->array_values+1)*sizeof(double *) );
		if(!all_array)fprintf(stderr, "get_low_boundary: cannot allocate all_array, bye"), exit(3);
	
		
		if( this-> w ){
		
			for(int i=0; i<this->array_values; i++)
				all_array[i] = array[i];            // merge all arrays with w array into a single one (for qsortn)
			all_array[this->array_values] = w;
		
			qsortn( (double **)all_array, (long)this->array_size, (long) this->array_values+1, X );  // sort by array i
		
			free(all_array);

		}
		else
			qsortn( (double **)array, (long)this->array_size, (long) this->array_values, X );  // sort by array i
		
			/*for(int i=0; i<this->array_values; i++, printf("\n"))
				for(int j=0; j<this->array_size; j++)
					printf("%f\t",array[i][j]);            
			*/
	
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
		
	}
		
	if( index < 0){
		//cerr << "the CI is too small for your distribution, bye\n";
		//exit(1);
		index++;
	}
	
	
	
	return this->array[X][index];


}

double distributions::get_high_boundary( int X ){
	
	int index=-1;

	if( this->full == false ){
		index = array_size/2;
	}
	else{
		/*
			Sort the values of all arrays as well as the weight array by array[i]
		*/
		double **all_array = (double **)malloc( (this->array_values+1)*sizeof(double *) );
		if(!all_array)fprintf(stderr, "get_low_boundary: cannot allocate all_array, bye"), exit(3);
	
		if( this->w ){
		
			for(int i=0; i<this->array_values; i++)
				all_array[i] = array[i];            // merge all arrays with w array into a single one (for qsortn)
			all_array[this->array_values] = w;
		
			qsortn( (double **)all_array, (long)this->array_size, (long) this->array_values+1, X );  // sort by array i
		
			free(all_array);

		}
		else
			qsortn( (double **)array, (long)this->array_size, (long) this->array_values, X );  // sort by array i
	

		if( this->w == NULL)
			index = floor(  array_size * (1.0-((1.0-this->CI)/2.0)) ) -1;
		else
		{
			float w_t=w[array_size-1];
			int c=array_size-1;
		
			while( w_t/sum_w < (1.00-this->CI)/2.0 )
				w_t += w[--c];
		
			index=c+1;
		}
		
		free(all_array);
	}
		
	if( index >= array_size){
		 index--;
		//cerr << "the CI is too small for your distribution, bye\n";
		//exit(1);
	}

	

	return this->array[X][index];

}

void distributions::print_distribution( int X ){
	
	double XXX=0;
	int i;
	
	for(i=0;i<array_size; i++){
	
		if(this->w){
			XXX += w[i]/sum_w;
			printf("%d %.3f %.5f %.5f %.5f %.5f\n", i, array[X][i], w[i], w[i]/sum_w, XXX, 1-XXX+w[i]/sum_w );
		}
		else
			printf("%d %.3f\n", i, array[X][i]);
	}
	
}
