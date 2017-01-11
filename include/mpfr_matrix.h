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

#ifndef MPFR_MATRIX_H_
#define MPFR_MATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif


#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>

#include "aux_funcs.h"

/*
void *safeMalloc(size_t size);
void *safeCalloc(size_t nmemb, size_t size);
void *safeRealloc(void *ptr, size_t size);
void safeFree(void *ptr);
*/

/* Allocates a n * m matrix and initializes all entries to prec
   bits
*/
mpfr_t *allocateMPFRMatrix(uint64_t n,uint64_t m, mp_prec_t prec);

/* Clears all entries of a n * m matrix A and frees the memory
   allocated to that matrix
*/
void freeMPFRMatrix(mpfr_t *A, uint64_t n, uint64_t m);

/* Allocates a n sized vector and initializes all entries to prec
   bits
*/
mpfr_t *allocateMPFRVector(uint64_t n, mp_prec_t prec);

/* Clears all entries of a n sized vector v and frees the memory
   allocated to that vector
*/
void freeMPFRVector(mpfr_t *v, uint64_t n);

/* Sets a n * m matrix A to all zeros */
void setMatrixZero(mpfr_t *A, uint64_t n, uint64_t m);

/* Sets a n * n matrix A to identity */
void setMatrixIdentity(mpfr_t *A, uint64_t n);

/* Copies a matrix exactly, changing the precision of the elements of
   the output matrix where needed
*/
void copyMatrix(mpfr_t *B, mpfr_t *A, uint64_t n, uint64_t m);






/* Displays a complex m x n MPFR matrix to the standart output (mpfr floats are preliminarly
converted to doubles before display). */
void MPFRComplexMatrixPrint( mpfr_t *reA ,mpfr_t *imA, uint64_t m, uint64_t n);
void MPFRComplexMatrixPrint2( FILE *file, mpfr_t *reA ,mpfr_t *imA, uint64_t m, uint64_t n);



/* Returns maximum precision of all elements of a complex m x n MPFR matrix whose
complex and imaginary parts are in matrices ReA and ImA respectively. */
mpfr_prec_t getMaxPrecision(mpfr_t *ReA, mpfr_t *ImA,  uint64_t m, uint64_t n);

/*
Read a floating-point matrix A of size m * n from file stream, using rounding direction rnd.
Matrix A is assumed to be declared and initialized outside this function. Precision of matrix A is set outside this function.
Format of input: floating-point numbers must be in base 10 in form A@B or AeB, where A is mantissa and B is exponent.
*/
void readMPFRMatrix(mpfr_t *A, FILE *stream, uint64_t m, uint64_t n, mpfr_rnd_t rnd);

/* Convert a real n * m matrix A represented in format clapack double to the MPFR matrix ReA.
Matrix ReA is assumed to be declared and pre-allocated outside the function.
THe function changes precision of ReA to 64.
*/
void doubleToMPFRMatrix(mpfr_t *ReA, double *A, int m, int n);



/* For a double m x n matrix A the function returns its maximum in absolute value
element, converted to MPFR. Output variable is assumed to be allocated outside the function and its
precision is not changes within the function. */
void getMaxInMPFR(mpfr_t max, double *A, uint64_t m, uint64_t n);

/* CHeck if any of the elements of a double m x n matrix is NaN.
Returns a non-zero value (true) if A has NaN value, and zero (false) otherwise. */
int matrixIsNan_mpfr(mpfr_t *A, uint64_t m, uint64_t n);

/*
Write to file stream a complex m * n matrix rounded in the direction rnd with its real and imaginary parts in ReA and ImA respectively.
The function prints nmbr significant digits exactly, or if nmbr is 0, enough digits
so that matrix could be read back exactly.
Format of output: first line is two difits, representing size of matrix.
Then values are printed in form "ReAij + i*Imij", separated with tabulation.
The function prints matrix in base 10.
*/
void writeMPFRComplexMatrix(FILE *stream, mpfr_t *ReA, mpfr_t *ImA, uint64_t m, uint64_t n,size_t nmbr, mpfr_rnd_t rnd);


/*
Write to file stream a real m * n matrix rounded in the direction rnd.
The function prints nmbr significant digits exactly, or if nmbr is 0, enough digits
so that matrix could be read back exactly.
Format of output: first line is two difits, representing size of matrix.
Then values are separated with tabulation.
The function prints matrix in base 10.
*/
void writeMPFRMatrix(FILE *stream, mpfr_t *A, uint64_t m, uint64_t n,size_t nmbr, mpfr_rnd_t rnd);

/* Get matrices of element-by-element precision of a complex MPFR m x n matrix whose
real and imaginary parts are contained in matrices ReA and ImA respectively. */
void getMPFRMatrixPrecision(mp_prec_t *ReA_p, mp_prec_t *ImA_p, mpfr_t *ReA, mpfr_t *ImA, uint64_t m, uint64_t n);

/* Rerurn a matrix containing element-by-element absolute values of m x n MPFR matrix A. */
void absMPFRMatrix(mpfr_t *Aabs,mpfr_t *A, uint64_t m, uint64_t n);

#ifdef __cplusplus
}
#endif


#endif /* MPFR_MATRIX_H_ */
