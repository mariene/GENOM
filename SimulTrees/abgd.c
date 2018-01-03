/******
        file     : abgg.c -- automatic barcod gap discovery
        function : rank values and find a gap in their density - devised to find the limit
	           between population and species in phylogeography (done with/for Nicolas Puillandre)
                                                 
        created  : April 2008
        modif    : Nov 09 --stable version-- (with a minimum of slope increase)
        modif    : April 10 --stable version-- (with (A) some minimum divergence and (B) at least 1.5 times of slope increase)
        modif    : June 10 --stable version-- 
	                   merged with Simultions from SimulTree
		  
        modif    : July 25 2010 --bug in find_gap
                   Aug 3rd - 3 bugs of over-reading an array (from sophieb using valgrind).
		  
        author   : gachaz
*****/


#define MIN_SLOPE_INCREASE 1.5

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

static char abgdDEBUG;
static char verbose;


#define SIGN( a ) ( ( (a) > 0 )?1: ( ((a)==0)?0:-1)  )
static int Increase(const void *v1, const void *v2){  	return (int)SIGN( *((double *)v1) - *((double *)v2));  };
#undef SIGN

#define ABS( x )  (((x)>0)?(x):(-x))




/*
	Ran1 from Numrec -- for bootstrap
*/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)

#define EPS FLT_MIN 
#define RNMX (1.0-EPS)              /* RNMX should be the largest floating value less than 1 (numrec p 280) */

/*
	Return a number between ]0.00,1.00] (( now certain because of a bug Apr 2007 ))
*/
long idum_ran1=0;
double ran1()
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (idum_ran1 <= 0 || !iy) {
		if (-idum_ran1 < 1) idum_ran1=1;
		else idum_ran1 = -idum_ran1;
		for (j=NTAB+7;j>=0;j--) {
			k=idum_ran1/IQ;
			idum_ran1=IA*(idum_ran1-k*IQ)-IR*k;
			if (idum_ran1 < 0) idum_ran1 += IM;
			if (j < NTAB) iv[j] = idum_ran1;
		}
		iy=iv[0];
	}
	k=idum_ran1/IQ;
	idum_ran1=IA*(idum_ran1-k*IQ)-IR*k;
	if (idum_ran1 < 0) idum_ran1 += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum_ran1;
	if ((temp=AM*iy) > RNMX) return (double)RNMX;
	else return (double)temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/* uniform_dev return a number [0,1[ */
double uniform_dev( void ){

	double dum= ran1();
	
	if(dum == 0.0)fprintf(stderr, "ran1() can also return 0\n");

	if(dum == 1.0)dum=0;
	
	return dum;
	
}


/*
	create a bootstrap array_out from an array_in of size n
	(sample without removing)
*/
void BootStrapArray(double * array_in, double * array_out, long n){

	long i=0;
	for(i=0;i<n;i++)
		array_out[i]=array_in[ ( int ) floor(uniform_dev()*n) ];

}






/********************

	Peak/Gap

*********************/


struct Peak {

	double Dist;
	double Rank;
	double theta_hat;
};


/*
	This function computes the derivative using
	winsiz and find its maximal value. Actually,
	it scans for a peak and find it summit using both
	left (low values) and right (high values) sides by
	averaging them.
*/
struct Peak FindFirstPeak( double *Array, long N, int winsiz, short output_slope, double *Pi, double MaxDist ){

	long i,             /* Indice of the slope */
	     top=0;         /* indices the summit of the preak */
	     
	     
	double *Slope;           /* Slope <=> derivative of array */
	double SlopeMax;         /* Current max of Slope */
	
	
	float Mean_i=0;          /* The a priori rank of the summit --> then it is set to j+0.5 */
	
	double Mean_dist=0.0;    /* The distance that corresponds to SlopeMax */
		
	struct Peak my_abgd;     /* The structure that contain both the estimated indice and distance */
		
	extern char abgdDEBUG;
	
			
	long wt;
	long ct;
	
	my_abgd.Dist = -1;
	my_abgd.Rank = -1;
	my_abgd.theta_hat = 0;
	
	
	/*
		Store Slope
	*/
	Slope = (double *)malloc(  (size_t) (N-winsiz+1) *sizeof(double) );
	if(!Slope )fprintf(stderr, "FindMaxDifferiental: cannot allocate SlopeS --%ld double--, bye\n", N-winsiz+1 ), exit(2);

	
	for(i=0; i <= N-winsiz ; i++)
		Slope[i] = (Array[i+winsiz-1]-Array[i])/(double)(winsiz-1);

	


	/*
		Compute Pi up to the maximum of interest
	*/

	

	for(i=1; i<N && Array[i] <= MaxDist ; i++);
	
	my_abgd.theta_hat = Pi[i-1];
	
/*	printf("theta is %f\n", Pi[i-1]);*/

	
	/*
		Print out stuf if needed
	*/

	
	if(output_slope)
		for(i=0; i <= N-winsiz ; i++)
			printf("slope %ld %.10f ; dist %ld %f\n", i, Slope[i], i,  Array[i] );

	if(abgdDEBUG){
		for(i=0;i<N;i++){
			if(i<N-winsiz)
				fprintf( stdout, "[%ld] ; Val= %.5f ; Slope= %f\n",i,Array[i], Slope[i] );
			else
				fprintf( stdout, "[%ld] ; Val= %.5f \n",i, Array[i]  );
		}
	}



	/*
		1. Find the first slope max
		
	*/

	i=0;                     /* start at the first Slope value */
	top=0;                   /* set the top of  peaks and the structure to default value */
	my_abgd.Dist = -1;
	my_abgd.Rank = -1;
	SlopeMax = 0;

	while( i<N-winsiz ){
		
			
		/*
			Find the first value with enough divergence --> get the peak
		*/
				
		while(   Slope[i+1] >= Slope[i] && i<N-winsiz )     /* Find a Local Maxima */
			i++;
		
		/*
			2. Explore Local Maxima on right (winsiz scale)
		*/
	
		top = i;
		while( i<N-winsiz  && ABS(i-top) <= winsiz/10  ){                 /* explore 10 from slope(top) on right */
			if(Slope[i] > Slope[top])
				top=i;
			i++;
		}


		if(abgdDEBUG)printf("  local maxima == top: %ld\n", top);


		/*
			2. Find Threshold Distance. Go from window_size up to 2, keeping track of the origin of the peak
		*/
		Mean_i=0;
		Mean_dist=0;
		
		ct = top;
			
		for( wt = winsiz-1; wt>=2; wt--){
			if(abgdDEBUG)
				printf("ct: %ld, wt: %ld ++  Diff[l]: %f -- Diff[r]: %f\n", ct, wt, Array[ct+wt-1]-Array[ct] , Array[ct+1+wt-1]-Array[ct+1]);

			if( Array[ct+wt-1]-Array[ct] <   Array[ct+1+wt-1]-Array[ct+1] && ct < N-wt-1)
				ct++;				
		}
		Mean_dist = (Array[ct]+Array[ct+1])/2.0;


		if(abgdDEBUG){
			printf(" final mean_dist= %f ; rank= [%ld, %ld]\n",Mean_dist, ct, ct+1 );
			printf("   final countdown == from top: %ld to ct: %ld ; MeanDist: %f\n", top, ct, (Array[ct]+Array[ct+1])/2.0);
			printf(  "(ct:%ld) Mean dist: %f vs %f in th  ++ Slope[top] %f vs SlopeMax %f\n",ct, Mean_dist  , 2.581 *  2*my_abgd.theta_hat, Slope[top], SlopeMax );

		}
		
		if( Mean_dist  > 2.581 *2* my_abgd.theta_hat &&           /* because theta estimates can be as half as true value */
		    Slope[top] > MIN_SLOPE_INCREASE*SlopeMax){                             /* really we need some slope jump ! */
		
		
/*			printf("++++ thetaF Pi[j]= %f, ws %d (prev: %f)\n", Pi[ct], winsiz, my_abgd.theta_hat); */
			

			my_abgd.Dist = Mean_dist;                       /* this the chosen candidate */
			my_abgd.Rank = ct+0.5;
			my_abgd.theta_hat = Pi[ct];


			/*
				if the Pi[ct] is smaller, set it as the new theta and run a next round
			*/
			if( Pi[ct] <= my_abgd.theta_hat && my_abgd.Rank == -1 ){ 
				
				SlopeMax = 0;
				i=0;
				top=0;
			}
			else{
				
				break;
			}

		}
		else{
			/* still not sure about this */
			
			SlopeMax= (Slope[top]>SlopeMax && Pi[ct] <= my_abgd.theta_hat )?Slope[top]:SlopeMax;   
			
			/*SlopeMax= (Slope[top]>SlopeMax  )?Slope[top]:SlopeMax; */ 
			
			i =  (top>ct)?top+1:ct+1;
			while( i<N-winsiz && Slope[i+1] <=  Slope[i] )           /* go downhill as far as we can */
				i++;

		}
				
	}

	free(Slope);

	return my_abgd;
}




struct Peak find_gap( double *Array, long N, long windsize_min, long windsize_max, short output_slope, double MaxDist  ){

	int c,i;
	int stable=0;
	
	int windsize_step = (windsize_min>10)?windsize_min/10:1;
	
	struct Peak my_abgd;

	double *Pi;              /* average from array[0] to array[i] */ 

	double stable_dist=-1;
	
	my_abgd.Dist = -1;
	my_abgd.Rank = -1;
	my_abgd.theta_hat = -1;

	Pi = (double *)malloc(  (size_t) N *sizeof(double) );
	if(!Pi )fprintf(stderr, "FindMaxDifferiental: cannot allocate Pi --%ld double--, bye\n", N ), exit(2);

	for(Pi[0]=Array[0], i=1; i<N ; i++)
		Pi[i] = (Array[i] + Pi[i-1]*i)/(i+1.0);


	for(c=windsize_min; c <= windsize_max && stable<3; c+=windsize_step){
		
					
			my_abgd = FindFirstPeak( Array, N, c, output_slope, Pi,  MaxDist);
			
			if(abgdDEBUG)
				fprintf(stderr,"abs( %f - %f )= %f vs %f\n", my_abgd.Dist,stable_dist, fabs(my_abgd.Dist-stable_dist) , 0.1*stable_dist );
			
			if( my_abgd.Dist != -1 && fabs( my_abgd.Dist-stable_dist) < 0.1*stable_dist ){
			
				stable++;
			}
			else{
				stable=1;
				stable_dist=my_abgd.Dist;
			}
			
			if(abgdDEBUG)
				printf("w: %d  d_peak: %f r_peak: %.1f\n", c, my_abgd.Dist, my_abgd.Rank);
			
		}

	if(my_abgd.Dist == -1){
		
		my_abgd.Dist=Array[N-1];
		my_abgd.Rank = N+0.5;
	
	}

	free(Pi);

	return my_abgd;
}



/********************

	Distance Matrix

*********************/


/*
	In mask 1 means yes take it, 0, don't use it
*/
double *matrix2list( int **Kmat, int n, char *mask, long *Nval ){

	int i,j;
	int nseq=0;
	
	double *Pairs;
	
	for(i=0;i<n;i++)
		nseq+=mask[i];          

	Pairs = (double *)malloc( ((nseq*(nseq-1))/ 2 )*sizeof(double) );
	if(!Pairs)fprintf(stderr, "matrix2list: cannot allocate Pairs, bye\n"), exit(4);


	*Nval=0;
	for(i=0; i<n;i++){
	
		if( mask[i] == 0 )
			continue;
		
		for(j=i+1;j<n; j++)
			if( mask[j] )
				Pairs[ (*Nval)++ ] = (double) Kmat[i][j];
	
	}

	return Pairs;

}



/********************

	  Groups (graph composantes) extraction

*********************/


struct Composante {

	int nc;           /* numnber of composantes */
	int nn;           /* number of nodes / sequences */
	int nm;           /* number of masked/excluded nodes (with -1 in node_compid) */

	int *node_compid;   /* the comp_id of each node -- this starts by generating this array */

	int *n_in_comp;     /* number of nodes in each composante */
	int **comp;         /* a list of node id in each composante */

};


/*
	recursive function that add a node to the current composante and 
	check for next ones in the composante
*/
void setcomp( int node, int compid, int * node_compid, int ** matrix, int n, double max_dist, char *mask){


	int j;
	
	node_compid[node] = compid;                                                    /* the node is added to the composante */

	for( j=0; j<n; j++  ){

		if( mask[j] == 0 )
			continue;
	
		if( matrix[node][j] < max_dist ){

			if( node_compid[j] == 0 )
				setcomp( j, compid, node_compid, matrix, n, max_dist, mask);     /* add an undiscovered connected node */
			else{
			
				if( node_compid[j] != compid ){
					printf("setcomp: really strange ??? This should not happen %d in comp %d && %d in comp %d\n", node, compid, j , node_compid[j] );
					exit(1);
				}
				else
					continue;
			}
		}
	}
}


/*
	take a distance matrix with a maximum distance
	and return a list with the composante id of each sequence
	-- all connections are propagated --
*/
struct Composante compute_node_compid(  int ** matrix, int n, double max_dist, char *mask ){


	int row=0;                      /* current row of the matrix while reading it */
	struct Composante my_comp;      /*  a structure that store most of composante features --see above-- */
	int i;

	
	my_comp.nm = 0;
	for(i=0; i<n; i++)
		if( mask[i] == 0 )
			my_comp.nm++;

	my_comp.nc        = 1;
	my_comp.nn        = n - my_comp.nm;
	my_comp.comp      = NULL;
	my_comp.n_in_comp = NULL;


	my_comp.node_compid = (int *)calloc( (size_t) n , (size_t)sizeof(int) );
	if( !my_comp.node_compid )fprintf(stderr, "compute_composante: cannot allocate comp, bye"),exit(4);


	for( row=0; row < n ; row++){
	
		if(my_comp.node_compid[row] > 0 || mask[row] == 0)                                          /* if stored already, just skip it */
			continue;
		
		setcomp( row, my_comp.nc, my_comp.node_compid, matrix, n, max_dist, mask );        /* otherwise create a new composante and propgate to all connections */
		
		my_comp.nc++;

	}
	

	/*
		Reset id numbers from [1,nc] to [0, nc[
		nb: masked entries will have the -1 composant value
	*/
	
	my_comp.nc--;

	for( i=0 ; i < n ; i++ )
		my_comp.node_compid[i]--;
	
	return my_comp;
}

/*
	From a list where each node has a composante, 
	generate the composante list
*/
struct Composante extract_composante(  int **matrix, int n, double max_dist, char *mask ){

	int i;
	struct Composante my_comp;
		
	
	/*
		First make an array where each node as a comp_id (from 0 to nc; masked ones are -1)
	*/
	my_comp = compute_node_compid(  matrix, n, max_dist, mask );


	/*
		From this array built composantes
	*/
	my_comp.n_in_comp = (int *)calloc( (size_t)my_comp.nc , (size_t)sizeof(int) );
	if(!my_comp.n_in_comp)fprintf(stderr, "extract_composante: cannot allocate my_comp.n_in_comp, bye\n"), exit(2);

	my_comp.comp = (int **)malloc( (size_t)my_comp.nc * (size_t)sizeof(int*) );
	if(!my_comp.comp)fprintf(stderr, "extract_composante: cannot allocate my_comp.comp, bye\n"), exit(2);
	
	
	/*
		Count how many nodes in each composante for
		for memory allocation optimization
	*/
	for( i=0;  i< n ;  i++ ){
	
		if( my_comp.node_compid[i] == -1 )continue;         /* a masked node */
		
		my_comp.n_in_comp[ my_comp.node_compid[i] ]++;

	}
	
	for( i=0;  i<my_comp.nc;  i++ ){		
		my_comp.comp[i] = (int *)malloc(  sizeof(int)*my_comp.n_in_comp[ i ] );
		if( !my_comp.comp[i] )fprintf(stderr, "compute_composante: cannot allocate my_comp.comp[%d], bye\n", i), exit(2);
	}
	
	
	/*
		reset the array to 0
	*/
	for( i=0;i<my_comp.nc;i++ )
		my_comp.n_in_comp[ i ] = 0;


	/*
		Store all nodes in its corresponding composante
		and re-count how many nodes in each composante
	*/
	for(i=0; i< n; i++ ){
	
		if( my_comp.node_compid[i] == -1)continue;         /* a masked entry */
	
		my_comp.comp[ my_comp.node_compid[i] ][ my_comp.n_in_comp[ my_comp.node_compid[i] ] ] = i;
		my_comp.n_in_comp[ my_comp.node_compid[i] ] ++;
	}
	

	return my_comp;
}
/*
	Use this function to split one composante (given by id) into several ones given by sub_comp)
*/
void update_composante(  struct Composante *main_comp, int id, struct Composante sub_comp ){

	int i;
	
	/*
		first a small check
	*/
	if( main_comp->n_in_comp[id] != sub_comp.nn )
		fprintf(stderr, "update_composante: make sure the size of the composante to split equals the number of nodes in the splitted one\n"), exit(1);
	
	
	/*
		Extend memory
	*/
	main_comp->n_in_comp = (int *)realloc( (void *) main_comp->n_in_comp, (size_t) (main_comp->nc+sub_comp.nc-1)*sizeof(int)  );
	main_comp->comp = (int **)realloc( (void *) main_comp->comp, (size_t) (main_comp->nc+sub_comp.nc-1)*sizeof(int *)  );
	if( !main_comp->n_in_comp || !main_comp->n_in_comp )
		fprintf(stderr, "update_composante: cannot reallocate main_comp->n_in_comp or main_comp->n_in_comp, bye\n"), exit(3);
	

	/*
		UPDATE COMP_ID 
		in sub_comp, -1:ignore, 0: leave it in the main_comp (ie. ignore), 1+, append it to the end of main_comp
	*/
	
	for(i =0; i < main_comp->nn+main_comp->nm ; i++){

		if( sub_comp.node_compid[ i ] > 0 )
			main_comp->node_compid[ i ] = main_comp->nc + sub_comp.node_compid[ i ];
	}
	
	
	
	/* UPDATE N_IN_COMP */	
	
	main_comp->n_in_comp[ id ] = sub_comp.n_in_comp[ 0 ];                      /* the first sub one replace the original comp */
	
	for(i =1; i < sub_comp.nc ; i++){
		main_comp->n_in_comp[ main_comp->nc + i-1 ] = sub_comp.n_in_comp[ i ];  /* the others are append to the main_com */
	}
	
	
	/* UPDATE COMP */	
	
	main_comp->comp[ id ] = sub_comp.comp[ 0 ];
	sub_comp.comp[ 0 ] = NULL;
	
	for(i =1; i < sub_comp.nc ; i++){
		
		main_comp->comp[ main_comp->nc+i-1 ] = sub_comp.comp[ i ];  /* the others are append to the main_com */
		sub_comp.comp[ i ] = NULL;
		
	}


	/* UPDATE NC */
	
	main_comp->nc = main_comp->nc + sub_comp.nc - 1;

}


void print_groups( struct Composante my_comp , int ** distmat, int n  ){


	int i,j;

	for(i=0; i<my_comp.nc; i++){
	
		printf("Group[ %d ] n: %d ; id:", i, my_comp.n_in_comp[i] );
		
		for(j=0; j< my_comp.n_in_comp[i]; j++)
			printf(" n%d ", my_comp.comp[i][j]);
		printf("\n");
	
	}
	
}

void free_composante(  struct Composante c  ){

	int i;

	free(c.node_compid);
	for(i=0;i<c.nc;i++)
		if( c.comp[i] != NULL  )
			free(c.comp[i]);
	free(c.n_in_comp);
	free(c.comp);
}

void reset_composante( struct Composante * c ){

	c->nc = 0;
	c->nn = 0;
	c->nm = 0;
	
	c->node_compid = NULL;
	c->n_in_comp   = NULL;
	c->comp        = NULL;

	
}



/******
	heuristics to get a min window size
******/
long min_ws( long nval ){

	if( nval > 10000 )
		return 1000;

	return (nval/10>1)?nval/10:1;
}




/*
	This is trhe CORE function
	it runs find_gap() --rank, derived and look for peak--,
	then it extracts the groups by extract_comp().
	If opt_recursion is set on (1), it recurssively apply find_gap and update_composante
	
*/
struct Composante Find_ABGD_Partition( int ** Matrix_dist, long n, int stability, double Prior , char opt_recursion, char opt_verbose ){


	double *ValArray=NULL;                /* a double array for the values */
	long NVal=n*(n-1)/2;                  /* array size */
		
	int windsize_min=0;                   /* the smallest wind_size */
	int windsize_max=0;                   /* the largest wind_size */

	long i,j;                             /* simple counting variable */
	
	struct Peak my_abgd;                  /* In this Structure, There is the Peak dist and the corresponding rank */

	int c;
	
	struct Composante comp;               /* group partition */

	char *mask;                           /* used to mask some row/col in the distance matrix -- you can extract only sub-part of the matrix */

	char output_groups=opt_verbose;       /* output groups is verbose is on */
	char output_slope=(abgdDEBUG);   /* only output slope if verbose level is very high */
	
	extern char verbose;                  /* a global variable used in the whole file */
  
	verbose = opt_verbose;                

	my_abgd.Rank = -1;                    /* these are the default value for a gap --> it means no significant gap */
	my_abgd.Dist = -1;
	my_abgd.theta_hat = 0;

	
	/*
		1.1 Built a mask
	*/
	mask=(char*)malloc( n*sizeof(char) );
	if(!mask)fprintf(stderr, "main: cannot allocate mask, bye\n");
	for(j=0; j<n; j++)
		mask[j]=1;

	/*
		1.2 from matrix, built a distrib using the mask
	*/
	ValArray = matrix2list( Matrix_dist, n, mask , &NVal);
	qsort((void *) ValArray, (size_t) NVal, (size_t) sizeof(double), Increase );

	

	/*
		2. Find the estimated peak of the derivative on windsize values
	*/

	if(windsize_min==0)windsize_min = min_ws( NVal );
	if(windsize_max==0 || windsize_max>NVal-1)windsize_max = NVal-1;
	

	my_abgd = find_gap( ValArray, NVal, windsize_min, windsize_max, output_slope, Prior  );

	if(verbose == 2)
		printf("First partition, Kmin=%f\n", my_abgd.Dist);
	
	
	
	/*
		3. Extract groups using the limit
	*/
	comp = extract_composante(  Matrix_dist, n, my_abgd.Dist, mask );

	

	/*
		Try to resplit each group using recursion startegy on already defined groups
	*/
	if(opt_recursion){

		int flag=1;                     /* if 0, do change in groups, if 1, need another round */
		int a,b;                        /* dummy counters */
		int nc;                         /* number of composantes from the first partition before sub-splitting */ 
		
		long nval=0;                    /* number of pairwise comparisons */
		double *vals;                   /* pairwise distances */
		
		struct Peak recursive_abgd;     /* structure for storing extra split distance */
		int round=1;                    /* how many recurssion round */
	
		while( flag ){
		
			flag=0;                 /* if no sub-split is done, do not start a new round */
			nc= comp.nc;            
			
			if(verbose)
				printf("** Round %d of sub-split **\n", round);
			
			for(a=0; a< nc; a++){
			
				struct Composante recursive_comp;
			
				reset_composante( &recursive_comp );                     /* needed for the free in case of no new group */
			
				bzero( (void *)mask, (size_t)n*sizeof(char) );   /* built mask used to only consider some cells of the matrix */ 
				for(b=0;b<comp.n_in_comp[a]; b++)
					mask[ comp.comp[a][b] ] = 1;
							
				vals = matrix2list( Matrix_dist, n, mask , &nval);                                /* built array of pairwise dist */
				qsort((void *) vals, (size_t) nval, (size_t) sizeof(double), Increase );	
				
				if( nval > 2 ){                                                           /* at least 3 sequences are needed */
	
					windsize_min = min_ws( nval );
					windsize_max= nval-1;
	
					recursive_abgd = find_gap( vals, nval, windsize_min, windsize_max, output_slope, Prior  );
					
					if(recursive_abgd.Rank != nval+0.5){
						
						/*printf("can recut group %d at dist %f\n", a, recursive_abgd.Dist);*/
						
						recursive_comp = extract_composante(  Matrix_dist, n, recursive_abgd.Dist, mask );
						
						if( recursive_comp.nc > 1 ){
							
							switch( verbose ){
							
								case 2:
									printf("Subsequent partition %s\n", (verbose)?"":"(details with -v)" );
									printf("theta_hat  : %g\n", recursive_abgd.theta_hat );
									printf("ABGD dist  : %f\n",  recursive_abgd.Dist);
									printf("ws         : [%d, %d]\n", windsize_min, windsize_max  );
									printf("Group id   : %d (%d nodes)\n",  a, recursive_comp.nn);
									printf("-> groups  : %d\n",  recursive_comp.nc);
														
									printf("Subgroups are:\n");
									print_groups( recursive_comp, Matrix_dist, n );
									printf("\n");
									break;
								
								case 1:
									printf("\tGroup[%d] --> %d groups\n",  a, recursive_comp.nc);
							}
	
							update_composante(  &comp, a, recursive_comp );
	
							/*print_groups( comp, distmat );*/
							
							flag=1;
	
						}
	
					}
				}
			
				free( vals );
				free_composante( recursive_comp );
			}
			round++;
		}	
	}
	
	if(verbose){
		printf("Final partition:\n");

		printf("Groups n   : %d\n",  comp.nc);
		i=j=comp.n_in_comp[0];
		for(c=1;c<comp.nc;c++){
			i=(comp.n_in_comp[c]<i)?comp.n_in_comp[c]:i;
			j=(comp.n_in_comp[c]>j)?comp.n_in_comp[c]:j;
		}
		printf(" comp min  : %ld\n",  i);
		printf(" comp max  : %ld\n",  j);
	}

	if(output_groups)
		print_groups( comp, Matrix_dist, n );


	free(ValArray);
	free(mask);

	return comp;

}

