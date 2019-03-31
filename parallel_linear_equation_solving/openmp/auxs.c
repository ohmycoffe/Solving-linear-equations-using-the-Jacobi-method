#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "auxs.h"

void show_matrix(double** array, int rows, int cols)
{
  int i, j;
  for (i=0; i<rows; i++)
  {

    for (j=0; j<cols; j++)
    {
      printf("%.3f ",array[i][j]);
    }
	printf("\n");
  }
}

void show_vector(double* array, int rows)
{
  int i, j;
  for (i=0; i<rows; i++)
    printf("%f ",array[i]);
  printf("\n");
}

void save_matrix(char *filename,  double** matrix, int rows, int cols)
{
	FILE *fp;
    fp = fopen(filename, "w+");
    if (fp == NULL)
	{
   		printf("Error while opening file");
   		exit(1);
	}

	int i, j;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			fprintf(fp, "%f ", matrix[i][j]);
		}
		if(i < rows - 1) fprintf(fp, "\n");
	}
	fclose(fp);
}

void save_vector(char *filename,  double* vector, int size)
{
	FILE *fp;
    fp = fopen(filename, "w+");
    if (fp == NULL)
	{
   		printf("Error while opening file");
   		exit(1); //[TODO] Zastapic bledem
	}
 	int i;
	for(i=0;i<size;i++)
	{
		
		fprintf(fp, "%f ", vector[i]);
		if (i < size - 1)
			fprintf(fp, "\n");
	}

	fclose(fp);
}



void calculate_b(double** A, double* x, double* b, int size)
{
 	int rows = size;
 	int cols = size;
	int i, j;
    	for (i=0; i<cols; i++)
  		{
    		b[i]  = 0;
    	}
 	for (i=0; i<rows; i++)
  	{
    	for (j=0; j<cols; j++)
  		{
    		b[i]  = b[i] + A[i][j] * x[j];
    	}
  	}
}

void zeroing_matrix(double** array, int rows, int cols)
{
  int i, j;
  for (i=0; i<rows; i++)
  {
    for (j=0; j<cols; j++)
  	{
    	array[i][j] = 0;
    }
  }

}

void zeroing_vector(double* vector, int rows)
{
  int i;
  for (i=0; i<rows; i++)
      vector[i] = 0;
}


void generate_matrix(double **matrix,  int size, int range, float offset)
{
	//rand(); // without this, first generaterd element is always the same - what is going on?[TODO]
	int r,c;
	for(r=0;r<size;r++)
	{
		for(c=0;c<size;c++)
		{
			matrix[r][c] = ((float)rand()/(float)(RAND_MAX)+offset)*(float)range;
		}
	}
}
void generate_vector(double *vector,  int size, int range, float offset)
{
	rand(); // without this, first generaterd element is always the same - what is going on?[TODO]

	 int r;
	float rand_number;
	for(r=0;r<size;r++)
	{
		rand_number = ((float)rand()/(float)(RAND_MAX)+offset)*(float)range;
		vector[r] = rand_number;

	}
}

void makeDiagDominant(double** matrix, int rows, int cols, float gain_muliplier)
{
	int i,j;
	double sum_in_row;

	for (i=0; i<rows; i++)
  	{
		sum_in_row = 0;
    	for (j=0; j<cols; j++)
  		{
  			if(i!=j)
    			sum_in_row = sum_in_row + fabs(matrix[i][j])* gain_muliplier;
    	}
    	matrix[i][i] = sum_in_row;
  	}

}

 int get_rowsNo(char *filename)
{
	FILE *fp;
    int count_rows = 0;
    char chr;

    fp = fopen(filename, "r");
    if (fp == NULL)
	{
   		printf("Error while opening file (f:get_rowsNo)");
    	exit(1);
	}
	chr = getc(fp);
	while (chr != EOF)
	{
		if(count_rows == 0)
			count_rows = count_rows+1;
		if(chr == '\n')
			{
				count_rows = count_rows +1;
			}
	chr = getc(fp);
	}
	return count_rows;
	fclose(fp);
}

void load_matrix(double** array, int rows, int cols, char *filename)
{
// [TODO] Jak to zrobic sprtyniej?
  FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL)
	{
   		printf("Error while opening file");
   	exit(1); //[TODO] Zastapic bledem
	}
  int i, j, ret;


  for (i=0; i<rows; i++)
  {
    for (j=0; j<cols; j++)
  	{
    	array[i][j] = 0;
    }
  }


  for (i=0; i<rows; i++)
  {
    for (j=0; j<cols; j++)
    {
	  ret = fscanf(fp,"%lf",&array[i][j]);

    }
  }
  fclose(fp);
}

void load_vector(double* vector, int rows, char *filename)
{
// [TODO] Jak to zrobic sprtyniej?
  FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL)
	{
   		printf("Error while opening file");
   	exit(1); //[TODO] Zastapic bledem
	}
  int i, ret;
  for (i=0; i<rows; i++)
      vector[i] = 0;
  for (i=0; i<rows; i++)
	  ret = fscanf(fp,"%lf",&vector[i]);
  fclose(fp);
 }
























