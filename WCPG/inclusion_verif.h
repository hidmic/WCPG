
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
 * inclusion_verif.h
 *
 *This is the header file, containing the prototypes for the eigensystem inclusion verification functions
 */

#ifndef INCLUSION_VERIF_H
#define INCLUSION_VERIF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include "clapack_linalg.h"
#include "glue_conversions.h"
#include "mpfi_matrix.h"
#include "mpfr_linalg.h"
#include "aux_funcs.h"

 typedef struct context_input_data_struct
 {
 	uint64_t n;
 	uint64_t k;
 	complexdouble *A;
 	complexdouble *v;
 	complexdouble lambda;
 	mpfi_t *reA, *imA;
 	mpfi_t *rev, *imv;
 	mpfi_t reLambda, imLambda;
 	mpfi_t *reEps, *imEps;
 } context_input_data;

 typedef struct context_algo_data_struct
 {
 	uint64_t n;
 	mpfr_prec_t prec;
 	mpfr_t *reR, *imR;
 	mpfi_t *reRint, *imRint;
 	mpfi_t *reC1, *imC1;
 	mpfi_t *reC, *imC;
 	mpfi_t *reZ, *imZ;

 	complexdouble *R;
 	//mpfi_t *reY, *imY;
 	//mpfi_t *reX, *imX;
 	//mpfi_t *reXX, *imXX;
 }context_algo_data;

void context_input_data_init(context_input_data *context_id, complexdouble *A, complexdouble *v, complexdouble lambda, uint64_t n,mpfr_prec_t prec);
void context_input_data_clear(context_input_data *context_id);
void context_algo_data_init(context_algo_data *context, uint64_t n, mpfr_prec_t prec);
void context_algo_data_clear(context_algo_data *context);

int checkEigensystemInclusion(mpfi_t *rev_corrected, mpfi_t *imv_corrected, mpfi_t relambda_corrected, mpfi_t imlambda_corrected, \
								complexdouble *A, complexdouble *v, complexdouble lambda, uint64_t n, doublereal eps_v, doublereal eps_lambda, mpfr_prec_t prec);


void compute_R(context_algo_data *algo, context_input_data *input);


/* Computes matrix [C1] = [A] - [lambda]*[I] */
void compute_C1(context_algo_data *algo, context_input_data *input);

/* Computes matrix [Z] = -R * ([C1] * v)
	where Z is of size n x 1 */
void compute_Z(context_algo_data *algo, context_input_data *input);

/* Computes matrix [C]:

	[Ctmp] = [C1]
	[Ctmp(:,k)] = [-v]
	[C] = [I] - R*[Ctmp]

	where C is of size n x n */
void compute_C(context_algo_data *algo, context_input_data *input);

#ifdef __cplusplus
}
#endif


 #endif
