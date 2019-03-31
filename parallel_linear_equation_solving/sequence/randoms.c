#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "auxs.h"
/* run this program using the console pauser or add your own getch, system("pause") or input loop */

int main(int argc, char *argv[])
{
	int rows, cols, i, size;
	double **A;
	double *b;
	double *x;
	double diagonal_domination_gain;
	srand(time(0));
	char *filename_A, *filename_b, *filename_x;

    switch(argc)
    {
  		case 1: // no arguments

			size = 10;                          		//1
			diagonal_domination_gain = 1.2;    			//2
			filename_A = "test_A.txt";          		//3
  			filename_b = "test_b.txt";          		//4
  			printf("Input arguments: 0\n I am going to generate demi data.\n Size of the problem will be %d\n Gain is equal %f\n", size, diagonal_domination_gain);
			break;
		case 2: // 1 argument
			size = atoi(argv[1]);               		//1
			diagonal_domination_gain = 1.2;			    //2
  			filename_A = "test_A.txt";          		//3
  			filename_b = "test_b.txt";          		//4
  			printf("Input arguments: 1\n Size of the problem is %d\n Gain is equal %f\n", size, diagonal_domination_gain);
  			break;
		case 3: // 2 argument
			size = atoi(argv[1]);               		//1
			diagonal_domination_gain = atof(argv[2]);   //2
  			filename_A = "test_A.txt";          		//3
  			filename_b = "test_b.txt";          		//4
  			printf("Input arguments: 2\n Size of the problem is %d\n Gain is equal %f\n", size, diagonal_domination_gain);
			break;
	  	case 5: // 4 arguments
			size = atoi(argv[1]);             			 //1
			diagonal_domination_gain = atof(argv[2]);    //2
			filename_A = argv[3];               		 //3
			filename_b = argv[4];               		 //4
			printf("Input arguments: 3\n I will generate A matrix called '%s', b vector called '%s'.\n Size of the problem is %d\n Gain is equal %f\n", filename_A, filename_b, size, diagonal_domination_gain);
			break;
		default:
			printf("BAD INPUT:\n 1: problem size,\n 2: gain of the diagonal domination, \n 3: file name where matrix A will be saved,\n 4: file name where vector b will be saved.\nThere can be only 0, 1 or 3 input arguments!");
			exit(0);
			break;
			
  }
	rows = size;
	cols = size;
/********MEMORY ALLOCATION*********/
/* Memory allocation for matrix A */
	A = malloc(rows * sizeof(*A));
	if(A == NULL)
		{
		    printf("Error! memory for A not allocated.");
		    exit(0);
		}
	for (i=0; i<rows; i++)
	{
	    A[i] = malloc(cols * sizeof(**A)); //x[i] = malloc(cols * c);
	    	if(A[i] == NULL)
		    {
		        printf("Error! memory for A not allocated.");
		        exit(0);
		    }
	}
	zeroing_matrix(A, rows, cols);
/* Memory allocation for  b */
  	b = malloc(rows * sizeof(*b));
    if(b == NULL)
    {
	    printf("Error! memory for b not allocated.");
	    exit(0);
	}
  	zeroing_vector(b, rows);
/* Memory allocation for  x */
  	x = malloc(rows * sizeof(*x));
	  	if(x == NULL)
	    {
		    printf("Error! memory for x not allocated.");
		    exit(0);
		}
  	zeroing_vector(x, rows);

	generate_vector(x, size, 1, 0);
 	generate_matrix(A, size, 1, 0);
 	if(diagonal_domination_gain != 0)
 		makeDiagDominant(A, rows, cols, diagonal_domination_gain);
 	save_matrix(filename_A, A, rows, cols);
 	printf("Matrix A has been saved successfully.\n");
 	calculate_b(A, x, b, size);
 	save_vector(filename_b, b, size);
 	printf("Vector b has been saved successfully.\n");
 	

/********PRINT DATA*********/
//  	printf("X:\n");
//  	show_matrix(A, rows, cols);
//  	printf("x:\n");
//  	show_vector(x, rows);
//  	printf("b_calc:\n");
//  	show_vector(b, rows);
//  show_vector(x_, rows);

/********MEMORY DEALLOCATION*********/
/* deallocate array A, b, x */
  for (i=0; i<rows; i++)
  {
    free(A[i]);
  }
  free(A);
  free(b);
  free(x);
  return 0;
}
