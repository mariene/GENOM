/*
	test the function of interest
*/

#include <iostream>
using namespace std;


#include <stdlib.h>
#include <math.h>

#define ABS( a ) (( (a) > 0 )?(a):(-a))
#define SIGN( a ) (( (a) > 0 )?1:-1)


typedef struct {
          int s;               /* the number of obs sites */
	  int n;               /* the number of seq sampled */
	  double half_beta;     /* the 1/2 * beta-risk in estimating theta (should be small) */
} DfuncArgs;



/*
	All to compute N choose N
*/
static long double LOG_factoriel_step( int32_t n, int32_t steps ){

	if( steps > n )printf("step should be smaller than n, bye\n"), exit(2);
	if( steps == 0 )return 0;
	return logl(n) + LOG_factoriel_step(n-1, steps-1);

}
static long double LOG_factoriel( int32_t n ){

	if(n<=1)return 0;
	return logl(n) + LOG_factoriel(n-1);
}
static long double log_N_choose_n1(int n, int n1){

	long double f=0;
	int n2=n-n1;

	if(n1>n2)
		f = LOG_factoriel_step(n, n-n1) - LOG_factoriel(n2);
	else
		f = LOG_factoriel_step(n, n-n2) - LOG_factoriel(n1);
 
	return f;
}


/* -------------------------------------------- */
/* compute gradient for macopt                  */
/* -------------------------------------------- */
static void mydfunc(long double x, long double *g, void *dfuncargs)
{
	int r=0;
	long double sign=1.0;
	long double log_ratio=0;
	long double add=0;
	long double log_combi=0;
	
	DfuncArgs * myFuncArg = (DfuncArgs *)dfuncargs;
	
	
	*g = 1.0;

	for(r=1; r<= myFuncArg->n -1; r++){
	
		sign = (r%2 == 0) ? -1.0 : 1.0;
		
		
		log_ratio = (myFuncArg->s+1)*logl( x / (x+ (long double)r) );
		log_combi = log_N_choose_n1( myFuncArg->n - 1, r );
				
		add = sign * expl( log_combi + log_ratio);
		*g -= sign * expl( log_combi + log_ratio);

	}

	*g -= myFuncArg->half_beta;
		
	return;
}




/*
	This function find the x value that gives
	a an output y very close to zero
	dedicated to my theta problem
*/
static long double find_root(   
	void (* pFunc)(long double x, long double *y, void *funcarg),
	void *dfuncargs,
	long double x_init,
	double step,
	double precision   /* difference to zero tolerated   by user */
){

	
	long double x=x_init;   /* x is the input value of the function */
	long double y;     /* y is the output */
	
	long double y_init;

	long double x_inf, x_sup;
		
	pFunc( x, &y_init, dfuncargs);
	step = (y_init>0)? step : -step;


	/*
		First part
	*/
	do{
		x += step;
		if(x<0)x=0;
		
		pFunc( x, &y, dfuncargs);
				
	}while( x >= 0.0 && SIGN(y) == SIGN(y_init) );
	
	/*
		Second part - oscillate until precision needed is found
	*/
	
	if(x==0 && SIGN(y) == SIGN(y_init) )
		return x;
	
	if( y>0 ){
		x_sup = x;
		x_inf = x-step;
	}else{
		x_sup = x-step;
		x_inf = x;
	}
	
		
	while( ABS(y)>precision ){
	
		if( y>0 ){
			x_sup=x;
			x = (x+x_inf)/2.0;
			pFunc( x, &y, dfuncargs);
			
			
		}else{
			x_inf =x;
			x = (x+x_sup)/2.0;
			pFunc( x, &y, dfuncargs);
		}

	}
	
	return x;
	
}


/*
	This compute theta_lower and theta_upper for a given beta risk
	using the Tavaré (1984) F cumulative distribution
*/
void compute_theta_CI( int Sites, int n_seq, double * theta_min, double * theta_max, double beta_risk, double step, double precision){

	void (* pFunc)(long double , long double *, void *) = mydfunc;

	double an=0.0;
	int i;
	
	for(i=1; i<n_seq; i++)
		an+=1.0/(double)i;

	DfuncArgs myFuncarg;
		
		
	myFuncarg.s = Sites;
	myFuncarg.n = n_seq;

	/*
		Find theta minimum for this #sites, then find theta maximum
	*/
	myFuncarg.half_beta = beta_risk/2.0;
	*theta_max = (double) find_root( pFunc, (void *)&myFuncarg, (double)Sites/an, step, precision);
	

	if( ! Sites )
		*theta_min = 0.0;
	else{
		myFuncarg.s = Sites-1;
		myFuncarg.half_beta = 1.00 - beta_risk/2.0;
		*theta_min = (double)  find_root( pFunc, (void *)&myFuncarg, (double)Sites/an, step, precision);
	}
	return;
}

#undef ABS
#undef SIGN
