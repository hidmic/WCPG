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
  Christoph Lauter

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

#include "mpfr_matrix.h"




/* Allocates a n * m matrix and initializes all entries to prec
   bits
*/
mpfr_t *allocateMPFRMatrix(uint64_t n,uint64_t m, mp_prec_t prec) {
  mpfr_t *A;
  uint64_t i, j;

  A = (mpfr_t *) wcpgSafeCalloc(n * m, sizeof(mpfr_t));
  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfr_init2(A[i * m + j], prec);
    }
  }

  return A;
}



/* Clears all entries of a n * m matrix A and frees the memory
   allocated to that matrix
*/
void freeMPFRMatrix(mpfr_t *A, uint64_t n, uint64_t m) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfr_clear(A[i * m + j]);
    }
  }

  wcpgSafeFree(A);
}

/* Allocates a n sized vector and initializes all entries to prec
   bits
*/
mpfr_t *allocateMPFRVector(uint64_t n, mp_prec_t prec) {
  mpfr_t *v;
  uint64_t i;

  v = (mpfr_t *) wcpgSafeCalloc(n, sizeof(mpfr_t));
  for (i=0;i<n;i++) {
    mpfr_init2(v[i], prec);
  }

  return v;
}

/* Clears all entries of a n sized vector v and frees the memory
   allocated to that vector
*/
void freeMPFRVector(mpfr_t *v, uint64_t n) {
  uint64_t i;

  for (i=0;i<n;i++) {
      mpfr_clear(v[i]);
  }

  wcpgSafeFree(v);
}


/* Sets a n * m matrix A to all zeros */
void setMatrixZero(mpfr_t *A, uint64_t n, uint64_t m) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfr_set_si(A[i * m + j], 0, GMP_RNDN); /* exact */
    }
  }
}

/* Sets a n * n matrix A to identity */
void setMatrixIdentity(mpfr_t *A, uint64_t n) {
  uint64_t i, j;

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      mpfr_set_si(A[i * n + j], 0, GMP_RNDN); /* exact */
    }
    mpfr_set_si(A[i * n + i], 1, GMP_RNDN); /* exact */
  }
}



/* Copies a matrix exactly, changing the precision of the elements of
   the output matrix where needed
*/
void copyMatrix(mpfr_t *B, mpfr_t *A, uint64_t n, uint64_t m) {
  uint64_t i,j;

  for (i=0;i<n;i++) {
    for (j=0;j<m;j++) {
      mpfr_set_prec(B[i * m + j], mpfr_get_prec(A[i * m + j]));
      mpfr_set(B[i * m + j], A[i * m + j], GMP_RNDN); /* exact */
    }
  }
}






void MPFRComplexMatrixPrint( mpfr_t *reA ,mpfr_t *imA, uint64_t m, uint64_t n)
{
	printf("MPFR matrix of size %llu x %llu \n", m, n);
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			printf("%f + i %f \t", mpfr_get_d(reA[i * n + j], MPFR_RNDN), mpfr_get_d(imA[i * n + j], MPFR_RNDN));
		}
		printf("\n");
	}
}

void MPFRComplexMatrixPrint2( FILE *file, mpfr_t *reA ,mpfr_t *imA, uint64_t m, uint64_t n)
{
	printf("MPFR matrix of size %llu x %llu \n", m, n);
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			fprintf(file, "%f + i %f \t", mpfr_get_d(reA[i * n + j], MPFR_RNDN), mpfr_get_d(imA[i * n + j], MPFR_RNDN));
		}
		fprintf(file, "\n");
	}
}



/*
Read a floating-point matrix A of size m * n from file stream, using rounding direction rnd.
Matrix A is assumed to be declared and initialized outside this function. Precision of matrix A is set outside this function.
Format of input: floating-point numbers must be in base 10 in form A@B or AeB, where A is mantissa and B is exponent.
*/
void readMPFRMatrix(mpfr_t *A, FILE *stream, uint64_t m, uint64_t n, mpfr_rnd_t rnd)
{

	int i, j;
	for(i = 0; i < m; ++i)
	{
	      for(j = 0; j < n; ++j)
	      {
		      mpfr_inp_str(A[i * n + j], stream, (int)10, rnd);
	      }
	}
}

/*
Write to file stream a complex m * n matrix rounded in the direction rnd with its real and imaginary parts in ReA and ImA respectively.
The function prints nmbr significant digits exactly, or if nmbr is 0, enough digits
so that matrix could be read back exactly.
Format of output: first line is two difits, representing size of matrix.
Then values are printed in form "ReAij + i*Imij", separated with tabulation.
The function prints matrix in base 10.
*/
void writeMPFRComplexMatrix(FILE *stream, mpfr_t *ReA, mpfr_t *ImA, uint64_t m, uint64_t n,size_t nmbr, mpfr_rnd_t rnd)
{
	fprintf(stream, "%d %d \n", (int)m, (int)n);
	int i, j;
	for(i = 0; i < m; ++i)
	{
	      for(j = 0; j < n; ++j)
	      {
		      mpfr_out_str(stream, (int)10, nmbr, ReA[i * n + j], rnd);
		      fprintf(stream, " + i* ");
		      mpfr_out_str(stream, (int)10, nmbr, ImA[i * n + j], rnd);
		      fprintf(stream, "\t");
	      }
	      fprintf(stream, "\n");
	}
}

void writeMPFRMatrix(FILE *stream, mpfr_t *A, uint64_t m, uint64_t n,size_t nmbr, mpfr_rnd_t rnd)
{
	fprintf(stream, "%d %d \n", (int)m, (int)n);
	int i, j;
	for(i = 0; i < m; ++i)
	{
	      for(j = 0; j < n; ++j)
	      {
		      mpfr_out_str(stream, (int)10, nmbr, A[i * n + j], rnd);
		      fprintf(stream, "\t");
	      }
	      fprintf(stream, "\n");
	}
}


/* Convert a real n * m matrix A represented in format clapack double to the MPFR matrix ReA.
Matrix ReA is assumed to be declared and pre-allocated outside the function.
THe function changes precision of ReA to 64.
*/
void doubleToMPFRMatrix(mpfr_t *ReA, double *A, int m, int n)
{
	int i, j;
	mp_prec_t prec = 64;

	for(i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			mpfr_set_prec(ReA[i * n + j], prec);
			mpfr_set_d(ReA[i * n + j], A[i * n + j], MPFR_RNDN);
		}
	}
}


void getMPFRMatrixPrecision(mp_prec_t *ReA_p, mp_prec_t *ImA_p, mpfr_t *ReA, mpfr_t *ImA, uint64_t m, uint64_t n)
{
	int i, j;
	mpfr_prec_t maxR, maxI;
	maxR = 0; maxI = 0;
	for(i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if(mpfr_zero_p(ReA[i * n + j])) ReA_p[i * n + j] = 0;
			// else ReA_p[i * n + j] = mpfr_get_prec(ReA[i * n + j]);
			if(mpfr_zero_p(ImA[i * n + j])) ImA_p[i * n + j] = 0;
			// else ImA_p[i * n + j] = mpfr_get_prec(ImA[i * n + j]);
			// printf("(%d, %d) \t",(int)ReA_p[i * n + j], (int)ImA_p[i * n + j]);

			if (mpfr_get_prec(ReA[i * n + j]) > maxR) maxR = mpfr_get_prec(ReA[i * n + j]);
			if (mpfr_get_prec(ImA[i * n + j]) > maxR) maxI = mpfr_get_prec(ImA[i * n + j]);
		}
		// printf("\n");

	}
	fprintf(stderr, "Max precisions: (%ld, %ld)\n", maxR, maxI);
}

mpfr_prec_t getMaxPrecision(mpfr_t *ReA, mpfr_t *ImA,  uint64_t m, uint64_t n)
{
	int i, j;
	mpfr_prec_t maxR, maxI;
	maxR = 0; maxI = 0;
	for(i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (mpfr_get_prec(ReA[i * n + j]) > maxR) maxR = mpfr_get_prec(ReA[i * n + j]);
			if (mpfr_get_prec(ImA[i * n + j]) > maxR) maxI = mpfr_get_prec(ImA[i * n + j]);
		}

	}
	int res = (maxR > maxI ? maxR : maxI);
	return res;
}


//The function computes the element-by-element abs of matrix A
void absMPFRMatrix(mpfr_t *Aabs,mpfr_t *A, uint64_t m, uint64_t n)
{
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			mpfr_abs(Aabs[i*n + j], A[i*n + j], MPFR_RNDN);
		}
	}

}


/* CHeck if any of the elements of a double m x n matrix is NaN.
Returns a non-zero value (true) if A has NaN value, and zero (false) otherwise. */
int matrixIsNan_mpfr(mpfr_t *A, uint64_t m, uint64_t n)
{
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
		    if(!mpfr_number_p(A[i*n + j])) return 1;
		}
	}
	return 0;

}


/* For a double m x n matrix A the function returns its maximum in absolute value
element, converted to MPFR. Output variable is assumed to be allocated outside the function and its
precision is not changes within the function. */
void getMaxInMPFR(mpfr_t max, double *A, uint64_t m, uint64_t n)
{
	double maxA = fabs(A[0]);
	double current = fabs(A[0]);
	int i,j;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			current = fabs(A[i * n + j]);
			if(current > maxA)
				maxA = current;
		}
	}
	mpfr_set_d(max, maxA, MPFR_RNDU);
}

