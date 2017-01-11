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


#ifndef GLUE_H
#define GLUE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include "clapack_functions_config.h"	


#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>

#include "mpfr_matrix.h"
#include "mpfi_matrix.h"
#include "clapack_linalg.h"



// extern mpfr_prec_t MaxPrec;


// void mpfr_init2_my(mpfr_t op, mpfr_prec_t prec);
// void mpfr_set_prec_my(mpfr_t op, mpfr_prec_t prec);


/* Convert a complex n * m matrix A represented in format clapack complexdouble to the format of two n * m mpfrt_t matrices 
ReA and ImA containing real and imaginary parts of A respectively, possibly changing their precision.
Matrices ReA and ImA are assumed to have been allocated outside this function. 
The function uses the MPFR variable scratch as scratch space, assuming they have been allocated, possibly changing its
precision and not clearing it.
*/
void complexdoubleToMPFRMatrix(mpfr_t *ReA, mpfr_t *ImA, complexdouble *A, int m, int n);

/* Convert a complex n * m matrix A represented in format clapack complexdouble to the format of two the interval n * m mpfi_t matrices 
ReA and ImA containing real and imaginary parts of A respectively, possibly changing their precision.
Matrices ReA and ImA are assumed to have been allocated outside this function. 
The function uses the MPFR variable scratch as scratch space, assuming they have been allocated, possibly changing its
precision and not clearing it.
*/
void complexdoubleToMPFIMatrix(mpfi_t *ReA, mpfi_t *ImA, complexdouble *A, int m, int n);


/* Convert a complex MPFR m x n matrix whose real and imaginary parts are in ReA and ImA to
a complexdouble matrix. The function converts element-by-element mpfr floats to double precision
using mpfr_get_d function */
void MPFRMatrixToDoublecomplex(mpfr_t *ReA, mpfr_t *ImA, complexdouble *A, int m, int n);

/* Convert a MPFR m x n matrix to a double matrix. The function converts element-by-element mpfr floats to double precision
using mpfr_get_d function */
void MPFRMatrixToDouble(double *Ad, mpfr_t *A, uint64_t m, uint64_t n);

/* Convert a complex MPFI m x n matrix whose real and imaginary parts are in ReA and ImA to
a complexdouble matrix. The function converts element-by-element mpfi floats to double precision */
void MPFIMatrixToDoublecomplex(mpfi_t *ReA, mpfi_t *ImA, complexdouble *A, int m, int n);

/* For a complex interval m x n matrix A the function returns its conversion to floating-point.
Output matrix B is assumed to be preallocated outside the function. Its precision may be changed.
 */

void MPFIMatrixtoMPFRMatrix(mpfr_t *reB, mpfr_t *imB, mpfi_t *reA, mpfi_t *imA, uint64_t m, uint64_t n);

 /* For a complex m x n MPFR matrix A the function returns its conversion to complex interval matrix. The result
 matrix B is assumed to be preallocated outside the function. Its precision may be changed. */

void MPFRMatrixToMPFIMatrix(mpfi_t *reB, mpfi_t *imB, mpfr_t *ReA, mpfr_t *ImA,uint64_t m, uint64_t n);

/* Convert a real n * m matrix A represented in format clapack doublereal to the MPFR matrix ReA.
Matrix ReA is assumed to be declared and pre-allocated outside the function.
THe function changes precision of ReA to 64.
*/
void doublerealToMPFRMatrix(mpfr_t *ReA, doublereal *A, int m, int n);

#ifdef __cplusplus
}
#endif


#endif 
