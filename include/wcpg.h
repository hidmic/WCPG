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
#ifndef WCPG_H
#define WCPG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <mpfr.h>

/** @brief For an LTI filter given in its State-Space representation {A,B,C,D},
where A is n*n, B is n*q, C is p*n and D is p*q real matrix the function 
returns integer value indicating if WCPG was successfully computed.
In p*q matrix W the Worst-Case peak gain is stored if algorithm successfully exited.
Input:
	A, B, C, D - pointers for double arrays representing filter in state-space realization
	n, p, q - order of filter, number of inputs and number of outputs respectively
	W (output) - if function succeeds, on the output will hold the p*q size WCPG matrix of the filter {A,B,C,D}
				space for W is assumed to be preallocated outside the function
Output:
	integer value equal to 1 if WCPG computation is successful and 0 otherwise.
 */
int WCPG_ABCD(double *W, double *A, double *B, double *C, double *D, uint64_t n, uint64_t p, uint64_t q);

/* For an LTI filter given in its State-Space representation {A,B,C,D},
where A is n*n, B is n*q, C is p*n and D is p*q double precision real matrix the function
returns integer value indicating if WCPG is successfully computed.
The function takes eps, a desired absolute error bound on the computed WCPG measure.
In p*q MPFR matrix W the Worst-Case peak gain is stored if algorithm successfully exited. 
Input:
	A, B, C, D - pointers for double arrays representing filter in state-space realization
	n, p, q - order of filter, number of inputs and number of outputs respectively
	W (output) - if function succeeds, on the output will hold the p*q MPFR WCPG matrix of the filter {A,B,C,D}
				matrix W is assumed to be preallocated outside the function
Output:
	integer value equal to 1 if WCPG computation is successful and 0 otherwise.*/
int WCPG_ABCD_mprec(mpfr_t *W, double *A, double *B, double *C, double *D, uint64_t n, uint64_t p, uint64_t q, mpfr_t mpeps);

/* For an LTI filter given in its State-Space representation {A,B,C,D},
where A is n*n, B is n*q, C is p*n and D is p*q MPFR real matrix the function
returns integer value indicating wether WCPG is successfully computed.
The function takes eps, a desired absolute error bound on the computed WCPG measure.
In p*q MPFR matrix W the Worst-Case peak gain is stored if algorithm successfully exited.
Input:
	A, B, C, D - pointers for double arrays representing filter in state-space realization
	n, p, q - order of filter, number of inputs and number of outputs respectively
	W (output) - if function succeeds, on the output will hold the p*q MPFR WCPG matrix of the filter {A,B,C,D}
				matrix W is assumed to be preallocated outside the function
Output:
	integer value equal to 1 if WCPG computation is successful and 0 otherwise.*/
int WCPG_mp(mpfr_t *Sk, mpfr_t *A, mpfr_t *B, mpfr_t *C, mpfr_t *D, uint64_t n, uint64_t p, uint64_t q, mpfr_t mpeps);


/** @brief Nth order LTI filter is represented by its transfer function numerator (array of size Nb) and
denumerator (array of size Na), where N := max(Na, Nb).
For such a filter, the function computes its WCPG in double precision
(i.e. such that the absolute error of WCPG computation is bounded by 2^53).
Input:
	num - pointer for the array holding numerator coefficients, sizeof(num) must be Nb
	denum - pointer for the array holding denumerator coefficients, sizeof(denum) must be Na
	W(output) - if function succeeds, on the output will hold the WCPG of filter represented with transfer function num/denum
				space for W is assumed to be preallocated outside the function
Output:
	integer value equal to 1 if WCPG computation is successful and 0 otherwise.
*/
int WCPG_tf(double *W, double *num, double *denum, uint64_t Nb, uint64_t Na);


/* Compute the lower bound on the WCPG matrix of a LTI filter
represented with state matrices A, B, C, D with a given 
absolute error bound 2^k.

W = W' + 2^k

where k < 0 is given in the argument.

Returns non-zero value (true) if computation is succesful and zero (false)
if could not compute lower bound on WCPG.
 */

int WCPGLowerBound(mpfr_t *W, mpfr_t* A, mpfr_t* B, mpfr_t* C, mpfr_t* D, uint64_t n, uint64_t p, uint64_t q, mpfr_exp_t k);

int WCPGLowerBound_double(mpfr_t *W, double* A, double* B, double* C, double* D, uint64_t n, uint64_t p, uint64_t q, mpfr_exp_t k);


/* For a real square MPFR matrix A of size n x n the function computes
 * the eigenvalues of the maitrx A and corresponding eigenvectors
 * in multiple precision, where all the internal computations are
 * performed with precision prec.
 *
 *  The function returns as well the approximate error bounds on
 *  the computed eigendecomposition: epsLambda for eigenvalues and
 *  epsV for eigenvectors.
 *
 *  All mpfr matrices are considered allocated and initilized outside the function!
 *
 * Function returns 1 in case of success and 0 therwise.
 *
 * */
int eigen(mpfr_t *re, mpfr_t *im, mpfr_t *Vre, mpfr_t *Vim, mpfr_t epsLambda, mpfr_t epsV, mpfr_t *A, uint64_t n, mpfr_prec_t prec);
void WCPG_welcome_message(int n);


/* Memory handling functions: setter and getter */
void wcpg_set_memory_funcs(void *(*)(size_t),
			   void *(*)(size_t, size_t),
			   void *(*)(void *, size_t),
			   void (*)(void *));

void wcpg_get_memory_funcs(void *(**)(size_t),
			   void *(**)(size_t, size_t),
			   void *(**)(void *, size_t),
			   void (**)(void *));

#ifdef __cplusplus
}
#endif

#endif  /* WCPG_H */
