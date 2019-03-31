#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "auxs.h"

#define MAX_IT  10000
#define EPSILON  0.000001

int jacobi(double** A, double* b, double* x_,int size, int max_it, double epsilon);

int main(int argc, char *argv[])
{
  int rows, cols, i, size;
  int k;
  double **A;
  double *b;
  double *x_;

  char *filename_A, *filename_b, *filename_x_;

  switch(argc)
  {
  	case 1: // no arguments
  		filename_A = "test_A.txt";
  		filename_b = "test_b.txt";
  		printf("0 input arguments: I am using demi data \n");
		break;
	case 3: // 2 arguments
		filename_A = argv[1];
		filename_b = argv[2];
		printf("2 input arguments: I am using A matrix from '%s' file and b vector from '%s' file\n", filename_A, filename_b);
		break;
	default:
		printf("BAD INPUT! 1: file name with matrix A, 2: file name with vector b. There can be only 0 or 2 input arguments!\n");
		exit(0);
		break;
  }


// get size of the A matrix
  size = get_rowsNo(filename_A);
  if (size != get_rowsNo(filename_b)) //checking if A and b have the same size
  {
  	printf("DIMENSION ERROR!\n");
  	printf("Size of matrix A (%d) and length of vector b (%d) are different!\n", size, get_rowsNo(filename_b));
  	exit(1);
  }
  printf("Size of the problem is equal to %d\n",size);
  rows = size;
  cols = size;

/********MEMORY ALLOCATION*********/
// Memory allocation for matrix A
  A = malloc(rows * sizeof(*A));
	  if(A == NULL)
	  {
	    printf("Error! memory for A not allocated.\n");
	    exit(0);
	  }
  for (i=0; i<rows; i++)
  {
    A[i] = malloc(cols * sizeof(**A)); //x[i] = malloc(cols * c);
	    if(A[i] == NULL)
		{
		    printf("Error! memory for A not allocated.\n");
		    exit(0);
		}
	  }
// Memory allocation for  b
  b = malloc(rows * sizeof(*b));
	  if(b == NULL)
	  {
		    printf("Error! memory for b not allocated.\n");
		    exit(0);
	  }
// Memory allocation for  x
  x_ = malloc(rows * sizeof(*x_));
	  if(x_ == NULL)
	  {
		    printf("Error! memory for x_ not allocated.\n");
		    exit(0);
	  }
/********LOADING DATA FROM FILE*********/
// loading data from file
  printf("Loading A ...");
  load_matrix(A,  size, size, filename_A);
  printf("OK\n");
  printf("Loading b...");
  load_vector(b,  size, filename_b);
  printf("OK\n");


/********CALCULATE x_estimated WITH TIME MEASURE*********/
//  int max_it = 10000000;
//  double epsilon = 0.000001;
  printf("Start solving equations \n");
  clock_t start_time, end_time;
  double cpu_time_used;
  start_time = clock();

  k = jacobi(A, b, x_, size, MAX_IT, EPSILON);

  end_time = clock();
  cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

/********SAVE DATA *********/
  save_vector("x_result.txt", x_, size);
/********PRINT DATA *********/
//  printf("X:\n");
//  show_matrix(A, rows, cols);
//  printf("b:\n");
//  show_vector(b, rows);
//  printf("x_estimated_calc:\n");
//  show_vector(x_, rows);

//Summary
  if (k>=0)
  	printf("Solution achieved in %d iteration. Elapsed time %f seconds\n", k, cpu_time_used);
  else
	printf("Solution isn't converged in %d iteration. Elapsed time %f seconds\n", MAX_IT, cpu_time_used);
 /********MEMORY DEALLOCATION*********/
/* deallocate array A, b, x_ */
  for (i=0; i<rows; i++)
  {
    free(A[i]);
  }
  free(A);
  free(b);
  free(x_);

  return 0;
}

int jacobi(double** A, double* b, double* x_,int size, int max_it, double epsilon)
{
	int i, j, k, rows, cols;
	double *x_0;
	rows = size;
	cols = size;
	// Estimated x from k-1
	x_0 = malloc(rows * sizeof(*x_0));
	  	if(x_0 == NULL)
	  	{
		    printf("Error! memory for x_0 not allocated.\n");
		    exit(0);
	 	}
/* INITIALIZE VECTORS */

	for(i=0; i<rows; i++)
	{
		x_0[i] = 0;
		x_[i] = 0;
	}
/* SOLVE LINEAR EQUATIONS */
	double ax, serr, rserr;
	k = 0;

	do //do_while
	{
		k = k + 1;
		for(i=0; i<rows; i++)
		{
			ax = 0;
			for(j=0; j<cols; j++)
				if(i!=j)
					ax = ax + (A[i][j]*x_0[j]);

			x_[i] = (-ax + b[i])/A[i][i];
		}
		serr = 0;
		for(i=0; i<rows; i++)
		{
			serr += pow(x_[i]-x_0[i],2);
			x_0[i]= x_[i];
		}
		rserr = sqrt(serr);
	} while(rserr>epsilon & k < max_it);
	if (k == max_it)
		k = -1;
	return k;
	
/* deallocate temporary array x-0 */
	free(x_0);
}
