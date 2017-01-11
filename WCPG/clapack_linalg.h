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
// "lapack_linalg.h"
This is the header file, containing the headers for
	- matrix input/output functions
	- functions for basic matrix arithmetic (multiplication, substraction, etc)
	- Linear System Solver
	- Eigensolver
All the functions work with matrices of CLAPACK's types doublereal and complexdouble.
All matrices are represented as one-dimension arrays, where for a n*m matrix A the element
A(i,j) is the element A[i * m + j] of the array.
 */

#ifndef CLAPACK_LINALG_H
#define CLAPACK_LINALG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "clapack_functions_config.h"
#include "aux_funcs.h"
 

/* clapack matrix functions */
void 	clapack_matrix_neg			(complexdouble *op, int rows, int columns);
void 	clapack_matrix_diagonal		(complexdouble *rop, complexdouble v, int n);
void	clapack_matrix_ident		(complexdouble *Identity, int m, int n);
void 	clapack_matrix_sub			(complexdouble *rop, complexdouble *op1, complexdouble *op2, int rows, int cols);
void 	clapack_matrix_mul			(complexdouble *C, complexdouble *A, complexdouble *B, int n);
void 	clapack_rmatrix_as_zmatrix	(complexdouble *res, doublereal *a, int m, int n);
void 	clapack_matrix_transp_sqr	(complexdouble *T, complexdouble *A, int n);
void 	clapack_matrix_copy_d		(doublereal *cpy, doublereal *src, int m, int n);
void 	clapack_matrix_copy_z		(complexdouble *cpy, complexdouble *src, int m, int n);

complexdouble complexdouble_mul(complexdouble op1, complexdouble op2);
doublereal abs_complexdouble(complexdouble *z);

void 	complexdoubleCopy			(complexdouble *Acopy, complexdouble *A, int m, int n);


/*Matrix concat functions */
void 	clapack_matrix_ver_concat	(complexdouble *rop, complexdouble *op1, complexdouble *op2, int cols, int rows1, int rows2);
void 	clapack_matrix_hor_concat	(complexdouble *rop, complexdouble *op1, complexdouble *op2, int rows, int col1, int col2);

/* Eigendecomposition related functions */

void 	eigvect_to_matrix		(complexdouble *V, doublereal *v,int* flags, int n);
void 	wiwr_to_matrix			(complexdouble *p, doublereal *wr, doublereal *wi, int n);
void 	eigval_flags			(int *flags, doublereal *wi, doublereal *wr, int n);
int 	clapack_eigenSolver		(complexdouble *p, complexdouble *V,doublereal *verrbnd, doublereal *eerrbnd, doublereal *A, int n, double eps);

/* Linear system solution */
//void 	clapack_LSSolver					(complexdouble *X,doublereal *ferr,doublereal *berr, complexdouble *A, complexdouble *B, integer n, integer bcols);

/* Matrix inverse */
int 	clapack_complex_matrix_inverse					(complexdouble *U, complexdouble *A, integer n);

/*Input & Output functions */
void 	clapack_matrix_inp_str_z		(complexdouble *A, int m, int k, FILE *stream);
int 	clapack_matrix_inp_str_d		(doublereal *A,int m, int k, FILE *stream);
void 	clapack_matrix_print_d			(doublereal *D, int mD, int kD);
void 	clapack_matrix_print_z			(complexdouble *D, int m, int n);
/* CHeck if any of the elements of a double m x n matrix is NaN.
Returns a non-zero value (true) if A has NaN value, and zero (false) otherwise. */
int matrixIsNan_double(doublereal *A, uint64_t m, uint64_t n);

/* CHeck if any of the elements of a complexdouble m x n matrix is NaN.
Returns a non-zero value (true) if A has NaN value, and zero (false) otherwise. */
int matrixIsNan_complexdouble(complexdouble *A, uint64_t m, uint64_t n);



#ifdef __cplusplus
}
#endif


 #endif
