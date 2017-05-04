/* Copyright 2016 by

  Laboratoire d'Informatique de Paris 6 - Équipe PEQUAN
  Sorbonne Universités
  UPMC Univ Paris 06
  UMR 7606, LIP6
  4, place Jussieu
  F-75252 Paris Cedex 05
  France

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France

  Contributors:
  Anastasia Volkova anastasia.volkova@lip6.fr


  This software is a mathematical library whose purpose is to provide
  functions to compute the Worst-Case Peak Gain measure of Linear
  Time-Invariant Digital Filters.

  This software is governed by the CeCILL-C license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-C
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-C license and that you accept its terms.

  This program is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


 */

/*
 ============================================================================
 Name        : exampleProgram.c
 Author      : Anastasia Volkova
 Version     :
 Copyright   : CeCILL-C license
 Description : Uses shared library to print greeting
               To run the resulting executable the LD_LIBRARY_PATH must be
               set to ${project_loc}/WCPG/.libs
               Alternatively, libtool creates a wrapper shell script in the
               build directory of this program which can be used to run it.
               Here the script will be called exampleProgram.
 ============================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "wcpg.h"
#include "mpfr_matrix.h"

int readMatrix(double *A,int m, int n, FILE *stream);
void printMatrix(FILE *stream, double *D, int m, int n);

int main(int argc, char *argv[] )
{

	/* Input files that we use for examples have following structure:
		n, n
		matrix_A, for wich column elements are separated with tabular and rows are each on new line
		n, p
		matrix_B
		q, n
		matrix_C
		q, p
		matrix_D
	*/



	/* File input interface */
	FILE *input;
	mpfr_t mpeps;
	mpfr_init(mpeps);


	if(argc < 3)
	{
		char s[100];
		char d[11];
		printf("Enter the path to an example file: ");
		scanf("%s", s);
		input = fopen(s, "r");
		printf("\nExample file: %s\n", argv[1]);
		printf("Enter the power of two to bound the WCPG 2^-X: X=");
		scanf("%s", d);

		mpfr_set_str(mpeps, d, 2, MPFR_RNDN);

	}
	else if(argc == 3)
	{
		input = fopen(argv[1], "r");
		printf("Example file: %s\n", argv[1]);
		mpfr_set_str(mpeps, argv[2], 2, MPFR_RNDN);

	}
	 if(input == NULL)
	 {
		 fprintf(stderr, "Error opening example file \n");
		 return 0;
	 }

	/* Reading matrices */
	int n;
	fscanf(input,"%d %d \n", &n, &n);
	double *A;
	A = (double*)wcpgSafeMalloc(n*n*sizeof(double));
	if(readMatrix(A,n, n, input))
	{
		fprintf(stderr, "Could not read matrix A from file. \n");
		return -1;
	}

	printf("Matrix A\n");
	printMatrix(stdout, A, n, n);

	double *B;
	int q;
	fscanf(input,"%d %d \n", &n, &q);
	B = (double*)wcpgSafeMalloc(n*q*sizeof(double));
	if(readMatrix(B,n, q, input))
	{
		fprintf(stderr, "Could not read matrix B from file. \n");
		return -1;
	}
	printf("Matrix B\n");
	printMatrix(stdout, B,n, q);

	double *C;
	int p;
	fscanf(input,"%d %d \n", &p, &n);
	C = (double*)wcpgSafeMalloc(p*n*sizeof(double));
	if(readMatrix(C,p, n, input))
	{
		fprintf(stderr, "Could not read matrix C from file. \n");
		return -1;
	}
	printf("Matrix C\n");
	printMatrix(stdout, C,p, n);

	double *D;
	fscanf(input,"%d %d \n", &p, &q);
	D = (double*)wcpgSafeMalloc(p*q*sizeof(double));
	if(readMatrix(D,p, q, input))
	{
		fprintf(stderr, "Could not read matrix D from file. \n");
		return -1;
	}
	printf("Matrix D\n");
	printMatrix(stdout, D,p, q);

	// Declaration and space allocation for the result wcpg approximation matrix
	//The default precision is 64, but the precision of the result variable is
	//modified by the algorithm in order to satisfy error bound eps.


	 /*-----------------------------------------------------------------------*/
	 //							Testing WCPG_ABCD
	 /*-----------------------------------------------------------------------*/
	printf("/*-----------------------------------------------------------------------*/\n \t\t\t\tTesting WCPG_ABCD \n/*-----------------------------------------------------------------------*/\n");
	double *W = (double*)wcpgSafeMalloc(p * q * sizeof(double*));
	 if (!WCPG_ABCD(W, A, B, C, D,n, p, q))
	 	printf("Could not compute WCPG \n");
	 else
	 {
	 	printf("\nWorst-case Peak Gain of WCPG_ABCD:\n");

	 	printMatrix(stdout, W, p, q);
	 }

	 /*-----------------------------------------------------------------------*/
	 //							Testing WCPG_ABCD_mprec
	 /*-----------------------------------------------------------------------*/
	printf("/*-----------------------------------------------------------------------*/\n \t\t\t\tTesting WCPG_ABCD_mprec \n/*-----------------------------------------------------------------------*/\n");
	mpfr_t *W_mprec;
	W_mprec = allocateMPFRMatrix(p,q, 64);
	 if (!WCPG_ABCD_mprec(W_mprec, A, B, C, D,n, p, q, mpeps))
	 	printf("Could not compute WCPG \n");
	 else
	 {
	 	printf("\nWorst-case Peak Gain of WCPG_ABCD_mprec:\n");

	 	printMatrix(stdout, W, p, q);
	 }
	 /*-----------------------------------------------------------------------*/
	 //							Testing WCPG_tf
	 /*-----------------------------------------------------------------------*/
	  printf("/*-----------------------------------------------------------------------*/\n \t\t\t\tTesting WCPG_tf \n/*-----------------------------------------------------------------------*/\n");

	  printf("This is a hardcoded example of a normalized transfer function.\n");


	  double num[3] = {1.0, 2.0, 3.0};
	  double denum[2] = {1.0/2.0, 1.0/3.0};
	  printf("Numerator:\n");
	  printMatrix(stdout, (double*)&num, 1, 3);
	  printf("Denumerator:\n");
	  printMatrix(stdout, (double*)&denum, 1, 2);

	  p = 1;
	  q = 1;
	  double *W2 = (double*)wcpgSafeMalloc(p * q * sizeof(double*));

	  if (!WCPG_tf(W2, num, denum, 3, 2, (int)0))
	  	printf("Could not compute WCPG \n");
	  else
	  {
	  	printf("\nWorst-case Peak Gain of WCPG_tf: %f \n", W2[0]);

		printMatrix(stdout, W2, 1, 1);
	  }



		 freeMPFRMatrix(W_mprec, p, q);
		 wcpgSafeFree(A);
		 wcpgSafeFree(B);
		 wcpgSafeFree(C);
		 wcpgSafeFree(D);
		 wcpgSafeFree(W);
		 wcpgSafeFree(W2);

  return 0;
}

int readMatrix(double *A,int m, int n, FILE *stream)
{

	int i, j;
	for(i = 0; i < m; ++i)
	{
	      for(j = 0; j < n; ++j)
	      {
		      if(!fscanf(stream, " %le ", &(A[i*n + j])))
		      {
		      		fprintf(stderr, "Problem reading matrix from file using scanf function. Process exit.\n");
		      		return 1;
		      }
		      if(isnan(A[i*n + j]))
		      {
		      		fprintf(stderr, "Problem reading matrix from file. One of the elements is a NaN. Process exit.\n");
		      		return 1;
		      }
	      }
	}
	return 0;

}

void printMatrix(FILE *stream, double *D, int m, int n)
{
	int i,j;
	fprintf(stream, "Matrix %d x %d \n", m, n);
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
		     fprintf(stream, "%e \t", (D[i*n + j]));
		}
		printf("\n");
	}
}

