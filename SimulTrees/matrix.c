
/*
	Straight Import from Num recipies
*/



#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NRANSI
#define NR_END 1
#define FREE_ARG char*

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

static void eigsrt(double d[], double **v, int n)
{
	int k,j,i;
	double p;

	for (i=0;i<n-1;i++) {
		p=d[k=i];
		for (j=i+1;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}


static void free_vector(double *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


static double *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) fprintf(stderr, "allocation failure in vector(), bye"), exit(3);
	return v-nl+NR_END;
}


static void jacobi(double **a, int n, double d[], double **v, int *nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=vector(0,n-1);
	z=vector(0,n-1);

	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++){
			v[ip][iq]=0.0;
		}
		v[ip][ip]=1.0;
	}
	
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	
		
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_vector(z,0,n-1);
			free_vector(b,0,n-1);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				
				g=100.0*fabs(a[ip][iq]);
				
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
					
				else if (fabs(a[ip][iq]) > tresh) {
				
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	fprintf(stderr, "Too many iterations in routine jacobi"), exit(1);
}
#undef ROTATE
#undef NRANSI


struct eigensystem Eigen_Symetric( double ** sym_matrix, long n ){
	
	struct eigensystem E;
	int i;
	int nrot;
	
	E.values = (double *)malloc( n*sizeof(double) );
	if(!E.values)fprintf(stderr, "Eigen_Symetric: cannot allocate E.values, bye\n"), exit(3);

	E.vectors = (double **)malloc( n*sizeof(double *) );
	if(!E.vectors)fprintf(stderr, "Eigen_Symetric: cannot allocate E.vectors, bye\n"), exit(3);
	
	for(i=0; i<n; i++){
		E.vectors[i] = (double *)malloc( n*sizeof(double) );
		if(!E.vectors[i])fprintf(stderr, "Eigen_Symetric: cannot allocate E.vectors[%d], bye\n", i), exit(3);
	}
	E.n = n;
	
	jacobi(sym_matrix, n, E.values, E.vectors, &nrot);	
	eigsrt(E.values, E.vectors, n);
	
	return E;
	
}


void free_Eigen_System( struct eigensystem e){
	int i;
	free(e.values);
	for(i=0; i<e.n; i++)free(e.vectors[i]);
	free(e.vectors);
	
}

