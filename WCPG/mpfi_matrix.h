/*------------------------------------------------------------------------------------------------------*/
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
"lapack_linalg.h"
This is the sourse file, containing the code for
	- MPFI matrix allocation/freeing functions
	- matrix input/output functions
	- functions for basic matrix arithmetic (multiplication, substraction, etc)
	- some complex interval numbers arithmetic operations (absolute value, multiplication, etc)
	- some functions for interaction of MPFI and MPFR matrices
All the functions work with complex interval multiple precision matrices of type mpfi_t.
Complex matrices are represented with two matrices holding its real and imaginary parts.
All matrices are represented as one-dimension arrays, where for a n*m matrix A the element
A(i,j) is the element A[i * m + j] of the array.


*/
/*-----------------------------------------------------------------------------------------------------*/


#ifndef MPFI_MATRIX_ALG_H
#define MPFI_MATRIX_ALG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpfi.h>
#include <mpfi_io.h>
#include <mpfr.h>
#include <clapack_linalg.h>
#include "mpfr_matrix.h"
#include <math.h>


mpfi_t *allocateMPFIMatrix(uint64_t n,uint64_t m, mp_prec_t prec);
void freeMPFIMatrix(mpfi_t *A, uint64_t n, uint64_t m) ;

void MPFIComplexMatrixAdd(mpfi_t *reC, mpfi_t *imC, mpfi_t *reA ,mpfi_t *imA, mpfi_t *reB, mpfi_t *imB, uint64_t m, uint64_t n);
void MPFIComplexMatrixSub(mpfi_t *reC, mpfi_t *imC, mpfi_t *reA ,mpfi_t *imA, mpfi_t *reB, mpfi_t *imB, uint64_t m, uint64_t n);
void MPFIComplexMatrixMultiply(mpfi_t *reC,mpfi_t *imC, mpfi_t *reA ,mpfi_t *imA, mpfi_t *reB, mpfi_t *imB, uint64_t m, uint64_t n, uint64_t p, mpfi_t *scratch);
void MPFRComplexMatrixMultiplyMPFIComplexMatrix(mpfi_t *reC,mpfi_t *imC, mpfr_t *reA ,mpfr_t *imA, mpfi_t *reB, mpfi_t *imB, uint64_t m, uint64_t n, uint64_t p, mpfi_t *scratch);
void MPFIComplexMatrixMultiplyMPFRComplexMatrix(mpfi_t *reC, mpfi_t *imC, mpfi_t *reA ,mpfi_t *imA, mpfr_t *reB, mpfr_t *imB, uint64_t m, uint64_t n, uint64_t p, mpfi_t *scratch);


void MPFIIdentMatrix(mpfi_t *reA, uint64_t n);
void MPFIZeroMatrix(mpfi_t *reA, uint64_t m, uint64_t n);

void MPFIComplexMatrixSetD(mpfi_t *reA, mpfi_t *imA, complexdouble *B,  uint64_t m, uint64_t n);
void MPFIComplexMatrixSet(mpfi_t *reA, mpfi_t *imA, mpfi_t *reB, mpfi_t *imB, uint64_t m, uint64_t n);
void MPFIComplexMatrixNeg(mpfi_t *reA, mpfi_t *imA, uint64_t m, uint64_t n);

// void MPFIComplexMatrixAbs(mpfr_t *reA, mpfr_t *imA, uint64_t n, uint64_t m);

void MPFIMatrixCopy(mpfi_t *B, mpfi_t *A, uint64_t m, uint64_t n);
void MPFIComplexMatrixPrint( mpfi_t *reA ,mpfi_t *imA, uint64_t m, uint64_t n);

/* Concatanation functions */
void MPFIComplexMatrixVerticalConcat(mpfi_t *reC,mpfi_t *imC, mpfi_t *reA ,mpfi_t *imA, mpfi_t *reB, mpfi_t *imB, uint64_t n, uint64_t m, uint64_t q, mpfi_t *scratch);
void MPFIComplexMatrixHorizontalConcat(mpfi_t *reC,mpfi_t *imC, mpfi_t *reA ,mpfi_t *imA, mpfi_t *reB, mpfi_t *imB, uint64_t n, uint64_t m, uint64_t q, mpfi_t *scratch);


void mpfi_mul_complexdouble(mpfi_t reC, mpfi_t imC, mpfi_t reA, mpfi_t imA, complexdouble b);
void mpfi_mul_complex(mpfi_t reC, mpfi_t imC, mpfi_t reA, mpfi_t imA, mpfi_t reB, mpfi_t imB, mpfi_t scratch);
void mpfi_mul_fr_complex(mpfi_t reC, mpfi_t imC, mpfr_t reA, mpfr_t imA, mpfi_t reB, mpfi_t imB, mpfi_t scratch);

void mpfi_abs_complex(mpfi_t absA, mpfi_t reA, mpfi_t imA, mpfi_t scratch);

void MPFIComplexMatrixMidrad(mpfi_t *reA, mpfi_t *imA, mpfi_t *reC, mpfi_t *imC, uint64_t m, uint64_t n, mpfr_t eps, mpfr_t *scratch3);
void MPFIComplexMatrixMidradDouble(mpfi_t *reA, mpfi_t *imA, complexdouble *C, uint64_t m, uint64_t n, mpfr_t eps, mpfr_t *scratch3);
void MPFIComplexMatrixMidrad_fr(mpfi_t *reA, mpfi_t *imA, mpfr_t *reC, mpfr_t *imC, uint64_t m, uint64_t n, mpfr_t eps, mpfr_t *scratch3);

void ComplexScalarMultiplyMPFIMatrix(mpfi_t *reC, mpfi_t *imC, mpfi_t reK, mpfi_t imK, mpfi_t *reA, mpfi_t *imA, uint64_t m, uint64_t n, mpfi_t scratch);

void MPFIMatrixNormMax(mpfi_t max, mpfi_t *A, uint64_t m, uint64_t n, mpfi_t scratch);

void mpfi_maxabs(mpfi_t max, mpfi_t *reA, mpfi_t *imA, uint64_t m, uint64_t n, mpfi_t *scratch);

void MPFIConstructDiagonal(mpfi_t *recI, mpfi_t *imcI, mpfi_t rec, mpfi_t imc, uint64_t n);

void DoubleMatrixMultiplyMPFIMatrix(mpfi_t *reC,mpfi_t *imC, doublereal *A, mpfi_t *reB ,mpfi_t *imB, uint64_t m, uint64_t n, uint64_t p, mpfi_t *scratch);
void MPFIMatrixMultiplyDoubleMatrix(mpfi_t *reC,mpfi_t *imC, mpfi_t *reB ,mpfi_t *imB, doublereal *A, uint64_t m, uint64_t n, uint64_t p, mpfi_t *scratch);


#ifdef __cplusplus
}
#endif


#endif
