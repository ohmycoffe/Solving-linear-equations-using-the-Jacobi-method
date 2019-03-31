/******************************************************************************
* FILE: jacobi_mpi.c
* DESCRIPTION:  
*   The Jacobi method of solving linear equations - C Version with MPI implementation 
*   In this code, the master task distributes a matrix 
*   operations to nummtasks-1 worker tasks.
*   NOTE:  A array is allocated in contiguous block of memory
* AUTHOR: Krzysztof Krolikowski.
* LAST modified: 27/01/2018
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "auxs.h"
#include <mpi.h>
#include "auxs.h"

#define MAX_IT  10000
#define EPSILON  0.000001
#define MASTER_ID 0
#define FROM_MASTER_TAG 100          
#define FROM_WORKER_TAG 200          
#define COMPUTE_TAG 1
#define STOP_TAG 0
// #define MATSIZE  10

int jacobi(double** A, double* b, double* x_,int size, int max_it, double epsilon);
void assign_worker(int worker_id);
void calculate_x(double* x_, double** A, double* b, double* x_0, int size, int offset, int nrows_chunk);



int main(int argc, char *argv[])
{
/**************************** MPI INITIALIZATION ******************************/
/**************************** mpi variables ***********************************/		
	int
		taskid,
		nummtasks,
		rc;
		
	MPI_Init(&argc,&argv);
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	MPI_Comm_size(MPI_COMM_WORLD,&nummtasks);
	if (nummtasks < 2 )
	{
		printf("Need at least two MPI tasks. Quitting...\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		exit(1);
	}
/**************************** GENERAL VARIABLES *********************************/
	if(taskid == MASTER_ID)
	{

		int 
			nrows,
			ncols, 
			i,
			size,
			k;
		double 
			**A,
			*b,
			*x_;
		char 
			*filename_A, 
			*filename_b, 
			*filename_x_;
			
/**************************** MENU ********************************************/
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
/**************************** GET SIZE OF THE PROBLEM *************************/
	// get size of the A matrix
		size = get_rowsNo(filename_A);
		if (size != get_rowsNo(filename_b)) //checking if A and b have the same size
		{
			printf("DIMENSION ERROR!\n");
			printf("Size of matrix A (%d) and length of vector b (%d) are different!\n", size, get_rowsNo(filename_b));
			exit(1);
		}
		printf("Size of the problem is equal to %d\n",size);
		printf("Number of processors is equal to %d\n",nummtasks);
/**************************** AUXILIARY VARIABLES *****************************/
		nrows = size;
		ncols = size;
/**************************** MEMORY ALLOCATION *******************************/
	// Matrix A:
		A = malloc(nrows * sizeof(*A));
			if(A == NULL)
			{
				printf("Error! memory for A not allocated.\n");
				exit(0);
			}		
		A[0] = malloc(nrows * ncols * sizeof (**A));
			if(A[0] == NULL)
			{
				printf("Error! memory for A not allocated.\n");
				exit(0);
			}
		for (i=1; i<nrows; i++)
			A[i] = A[0] + i * (nrows);
	// Vector b:
		b = malloc(nrows * sizeof(*b));
		if(b == NULL)
			{
				printf("Error! memory for b not allocated.\n");
				exit(0);
			}
	// Vector x_:
		x_ = malloc(nrows * sizeof(*x_));
		if(x_ == NULL)
		{
			printf("Error! memory for x_ not allocated.\n");
			exit(0);
		}
/**************************** LOADING DATA FROM FILE **************************/
		printf("Loading A ...");
		load_matrix(A,  size, size, filename_A);
		printf("OK\n");
		printf("Loading b...");
		load_vector(b,  size, filename_b);
		printf("OK\n");
/**************************** RUN JACOBI **************************************/
		double start_time, end_time;
		double cpu_time_used;	
		printf("I am starting...\n");
		start_time = MPI_Wtime(); 
		k = jacobi(A, b, x_, size, MAX_IT, EPSILON);
		end_time = MPI_Wtime(); 
		cpu_time_used = end_time - start_time;     

/**************************** SAVING RESULTS TO FILE **************************/	
		save_vector("x_result.txt", x_, size);
/**************************** SUMMARY *****************************************/
	  if (k>=0)
		  printf("Solution achieved in %d iteration. Elapsed time %f seconds\n", k, cpu_time_used);
	  else
		  printf("Solution isn't converged in %d iteration. Elapsed time %f seconds\n", MAX_IT, cpu_time_used);	
/**************************** MEMORY DEALLOCATION *****************************/	
	  free(A[0]);
	  free(A);
	  free(b);
	  free(x_);
 /**************************** WORKERS TASKS ***********************************/	
	}
	else
	{
		assign_worker(taskid); 
	}
  
  
  
  MPI_Finalize();
  return 0;
}

int jacobi(double** A, double* b, double* x_,int size, int max_it, double epsilon)
{
	int
			i,
			j,
			k,
			nrows,
			ncols;
	double
			*x_0;
	int
			taskid,
			nummtasks,
			numworkers,
			worker_id,
			avg_nrow,
			remind_nrow,
			offset,
			nrows_chunk;
		
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	MPI_Comm_size(MPI_COMM_WORLD,&nummtasks);
	nrows = size;
	ncols = size;
	numworkers = nummtasks-1;
	avg_nrow = size/numworkers;
	remind_nrow = size%numworkers;
// x_0 ALLOCATION	
	x_0 = malloc(nrows * sizeof(*x_0));
	if(x_0 == NULL)
	{
		printf("Error! memory for x_0 not allocated.\n");
		exit(0);
	}
	for(i=0; i<nrows; i++)
	{
		x_0[i] = 0;
	}
/**************************** MASTER TASK ************************************/
	if (taskid == MASTER_ID)
	{
		double ax, serr, rserr;
		int nrows_chunk;
		offset = 0;
		for (worker_id=1; worker_id < nummtasks; worker_id++)
		{
			nrows_chunk = (worker_id <= remind_nrow) ? avg_nrow+1 : avg_nrow;
			MPI_Send(&ncols, 1, MPI_INT, worker_id, FROM_MASTER_TAG, MPI_COMM_WORLD);
			MPI_Send(&offset, 1, MPI_INT, worker_id, FROM_MASTER_TAG, MPI_COMM_WORLD);
			MPI_Send(&nrows_chunk, 1, MPI_INT, worker_id, FROM_MASTER_TAG, MPI_COMM_WORLD);
			MPI_Send(&b[offset], nrows_chunk, MPI_DOUBLE, worker_id, FROM_MASTER_TAG, MPI_COMM_WORLD);
			MPI_Send(&A[offset][0], nrows_chunk*ncols, MPI_DOUBLE, worker_id, FROM_MASTER_TAG, MPI_COMM_WORLD);
			
			offset = offset + nrows_chunk;
		}
	/* SOLVE LINEAR EQUATIONS */
		k = 0;
		do //do_while
		{
			k = k + 1;
			for (worker_id=1; worker_id < nummtasks; worker_id++)
			{	
				MPI_Send(&k, 1, MPI_INT, worker_id, COMPUTE_TAG, MPI_COMM_WORLD);
				MPI_Send(x_0, ncols, MPI_DOUBLE, worker_id, FROM_MASTER_TAG, MPI_COMM_WORLD);
			}
			for (worker_id=1; worker_id < nummtasks; worker_id++)
			{
				MPI_Recv(&offset, 1, MPI_INT, worker_id, FROM_WORKER_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&nrows_chunk, 1, MPI_INT, worker_id, FROM_WORKER_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&x_[offset], nrows_chunk, MPI_DOUBLE, worker_id, FROM_WORKER_TAG, MPI_COMM_WORLD, &status);		
								
			}
			serr = 0;
			for(i=0; i<nrows; i++)
			{
				serr += pow(x_[i]-x_0[i],2);
				x_0[i]= x_[i];
			}
			rserr = sqrt(serr);	
		} while(rserr>epsilon & k < max_it);
		for (worker_id=1; worker_id < nummtasks; worker_id++)
			MPI_Send(&k, 0, MPI_INT, worker_id, STOP_TAG, MPI_COMM_WORLD);
		if (k == max_it)
			k = -1;
		return k;
	}

/* deallocate temporary array x-0 */
	free(x_0);
}

void assign_worker(int worker_id)
{
		int i,j;
		int k, offset, nrows_chunk, ncols;
		double
			**A_tmp,
			*b_tmp,
			*x__tmp,
			*x_0_tmp;
		MPI_Status status;
		MPI_Recv(&ncols, 1, MPI_INT, MASTER_ID, FROM_MASTER_TAG, MPI_COMM_WORLD, &status);	
		MPI_Recv(&offset, 1, MPI_INT, MASTER_ID, FROM_MASTER_TAG, MPI_COMM_WORLD, &status);	
		MPI_Recv(&nrows_chunk, 1, MPI_INT, MASTER_ID, FROM_MASTER_TAG, MPI_COMM_WORLD, &status);
	// Matrix A:
		A_tmp = malloc(nrows_chunk * sizeof(*A_tmp));
			if(A_tmp == NULL)
			{
				printf("Error! memory for A not allocated.\n");
				exit(0);
			}
		A_tmp[0] = malloc(nrows_chunk * ncols * sizeof (**A_tmp));
			if(A_tmp[0] == NULL)
			{
				printf("Error! memory for A not allocated.\n");
				exit(0);
			}
		for (i=1; i<nrows_chunk; i++)
			A_tmp[i] = A_tmp[0] + i * (ncols);
	// Vector b:
		b_tmp = malloc(nrows_chunk * sizeof(*b_tmp));
		if(b_tmp == NULL)
			{
				printf("Error! memory for b not allocated.\n");
				exit(0);
			}
	// Vector x_:
		x__tmp = malloc(nrows_chunk * sizeof(*x__tmp));
		if(x__tmp == NULL)
		{
			printf("Error! memory for x_ not allocated.\n");
			exit(0);
		}
	// Vector x_0:
		x_0_tmp = malloc(ncols * sizeof(*x_0_tmp));		
		if(x_0_tmp == NULL)
		{
			printf("Error! memory for x_ not allocated.\n");
			exit(0);
		}
		MPI_Recv(&b_tmp[0], nrows_chunk, MPI_DOUBLE, MASTER_ID, FROM_MASTER_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&A_tmp[0][0], nrows_chunk*ncols, MPI_DOUBLE, MASTER_ID, FROM_MASTER_TAG, MPI_COMM_WORLD, &status);
		while (1)
		{
			MPI_Recv(&k, 1, MPI_INT, MASTER_ID, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if (status.MPI_TAG == STOP_TAG)
				break;
			MPI_Recv(x_0_tmp, ncols, MPI_DOUBLE, MASTER_ID, FROM_MASTER_TAG, MPI_COMM_WORLD, &status);
			calculate_x(x__tmp, A_tmp, b_tmp, x_0_tmp, ncols, offset, nrows_chunk);
			MPI_Send(&offset, 1, MPI_INT, MASTER_ID, FROM_WORKER_TAG, MPI_COMM_WORLD);
			MPI_Send(&nrows_chunk, 1, MPI_INT, MASTER_ID, FROM_WORKER_TAG, MPI_COMM_WORLD);
			MPI_Send(&x__tmp[0], nrows_chunk, MPI_DOUBLE, MASTER_ID, FROM_WORKER_TAG, MPI_COMM_WORLD);
		}
		free(A_tmp[0]);
		free(A_tmp);
		free(b_tmp);
		free(x__tmp);
		free(x_0_tmp);
}		

void calculate_x(double* x_, double** A, double* b, double* x_0, int size, int offset, int nrows_chunk)
{
	int i, j;
	double ax;	
	for(i=0; i< nrows_chunk; i++)
	{
		ax = 0;
		for(j=0; j<size; j++)
			if((i+offset)!=j)
				ax = ax + (A[i][j]*x_0[j]);
		x_[i] = (-ax + b[i])/A[i][i+offset];
	}
}
