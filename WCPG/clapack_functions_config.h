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

#ifndef SRC_LIB_CLAPACK_FUNCTIONS_CONFIG_H_
#define SRC_LIB_CLAPACK_FUNCTIONS_CONFIG_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include "config.h"
#include "aux_funcs.h"


#include <f2c.h>

#ifdef CLAPACK_HEADER
	#include <clapack.h>
#endif

#ifdef LAPACK_HEADER
	#include <lapack.h>
#endif

#ifdef LAPACKE_HEADER
	#include <lapacke.h>
#endif

#ifdef BLAS_HEADER
	#include <blas.h>
#endif

#ifdef CBLAS_HEADER
	#include <cblas.h>
#endif

#ifdef BLASWRAP_HEADER
	#include <blaswrap.h>
#endif
 
#include <complex.h>
typedef struct complexdouble
{
	double r;
	double i;
}complexdouble;

#define doublereal double
#define lapack_int int
#define integer long int 
#define lapacke_complex double _Complex




void my_zgetrf(int *m, int *n, complexdouble *A,
 	int *lda, int *ipiv, int *info);

void my_zgetri(int *n, complexdouble *A, int *lda,
		int *ipiv, complexdouble *work, int *lwork, int *info);
void my_dgeevx(int* n, double *A, int* lda, double *wr, double *wi, double *vl,  
	       int* ldvl, double *vr, int* ldvr, int *ilo, int *ihi,  
	       double *scale, double *abnrm, double *rconde, double *rcondv,  
	       double *work, int *lwork, long int *iwork, int *info);


#ifdef __cplusplus
}
#endif


#endif /* SRC_LIB_CLAPACK_FUNCTIONS_CONFIG_H_ */
