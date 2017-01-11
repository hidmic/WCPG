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


#ifndef _WCPG_H_
#define _WCPG_H_


#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>

#include "glue_conversions.h"
#include "inclusion_verif.h"
#include "mpfr_linalg.h"
#include "clapack_functions_config.h"
#include "mpfr_eigendecomposition.h"
#include "aux_funcs.h"
#include "mpfi_matrix.h"


//Structure wcpg_result contains the information that will be saved after WCPG computation
typedef struct wcpg_result_struct
{
	int N;
	mpfr_t one_minus_rhoA;	
	mpfr_t maxSN;
	mpfr_t minSN;
	double time_overall;
	double time_Ncomp;
	double time_Summation;
	int inversion_Iter;
	mpfr_prec_t maxprec_PN;
	mpfr_prec_t maxprec_U;
	mpfr_prec_t maxprec_SN;


}wcpg_result;

void wcpg_result_init(wcpg_result *result);
void wcpg_result_clear(wcpg_result *result);
void wcpg_result_print(FILE *stream, wcpg_result *result, size_t ndigits);

/* THe function determines the minimal and maximum values of the
computed WCPG matrix and adds them to the corresonding fields
of the result structure. 
 */
void getDeltaOfWCPG(mpfr_t *reS, uint64_t p, uint64_t q, wcpg_result *result);

/* Support functions for N_lower_bound algorithm. Their placement to be duscussed later.
A priori, user does not need access to them (though possibly to lowerBoundN function).
However, circular dependency of headers makes result structure xcpg_result invisible, therefore
for the moment these functions are left here. */
void R_l(mpfi_t *reRl, mpfi_t *imRl, mpfi_t *rePhi, mpfi_t *imPhi, mpfi_t *rePsi, mpfi_t *imPsi, uint64_t n, uint64_t p, uint64_t q, int l);

int lowerBoundN(mpfi_t *rePhi, mpfi_t *imPhi, mpfi_t *rePsi, mpfi_t *imPsi, \
					 mpfi_t *reLambda, mpfi_t *imLambda,mpfr_t eps, uint64_t n, \
					  uint64_t p, uint64_t q, mpfr_prec_t prec, wcpg_result *context);
					  

/* Compute the lower bound on the WCPG matrix of a LTI filter
represented with state matrices A, B, C, D with a given 
absolute error bound 2^k:
	W = W' + 2^k
where k < 0 is given in the argument.

Returns non-zero value (true) if computation is succesful and zero (false)
if could not compute lower bound on WCPG. */
int WCPGLowerBound(mpfr_t *W, mpfr_t* A, mpfr_t* B, mpfr_t* C, mpfr_t* D, uint64_t n, uint64_t p, uint64_t q, mpfr_exp_t k);


/* For an LTI filter given in its State-Space representation {A,B,C,D}, 
where A is n*n, B is n*q, C is p*n and D is p*q real matrix,
and for an eps>0 the function returns a multi-precision n*q real matrix Sk of Worst-Case Peak Gains 
of the system such that the overall error of computation is bounded by eps. */
int WCPG(mpfr_t *Sk, mpfr_t *A, mpfr_t *B, mpfr_t *C, mpfr_t *D, mpfr_t mpeps, uint64_t n, uint64_t p, uint64_t q, wcpg_result *result);
int MPFItrustedEigendecompsition(mpfi_t *reLambdaint, mpfi_t *imLambdaint, mpfi_t *reVint, mpfi_t *imVint, mpfr_t *A, uint64_t n);

#ifdef __cplusplus
}
#endif

#endif  /* _WCPG_H */
