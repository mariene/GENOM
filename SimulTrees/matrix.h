
#ifndef _MYMATRIX_H_
#define _MYMATRIX_H_

struct eigensystem
{
	double *values;
	double **vectors;	
	long n;
};

struct eigensystem Eigen_Symetric( double ** sym_matrix, long n );

void free_Eigen_System( struct eigensystem e);

#endif
