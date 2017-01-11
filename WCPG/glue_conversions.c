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


#include "glue_conversions.h"

// mpfr_prec_t MaxPrec = 0;

// void mpfr_init2_my(mpfr_t op, mpfr_prec_t prec)
// {
// 		mpfr_init2(op, prec);
// 		if (prec > MaxPrec)
// 			MaxPrec = prec;
		
// }

// void mpfr_set_prec_my(mpfr_t op, mpfr_prec_t prec)
// {
// 		mpfr_set_prec(op, prec);
// 		if (prec > MaxPrec)
// 			MaxPrec = prec;
		
// }



/* Convert a real n * m matrix A represented in format clapack doublereal to the MPFR matrix ReA.
Matrix ReA is assumed to be declared and pre-allocated outside the function.
THe function changes precision of ReA to 64.
*/
void doublerealToMPFRMatrix(mpfr_t *ReA, doublereal *A, int m, int n)
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

/* Convert a complex n * m matrix A represented in format clapack complexdouble to the format of two n * m mpfrt_t matrices 
ReA and ImA containing real and imaginary parts of A respectively, possibly changing their precision.
Matrices ReA and ImA are assumed to have been allocated outside this function. 
The function uses the MPFR variable scratch as scratch space, assuming they have been allocated, possibly changing its
precision and not clearing it.
*/
void complexdoubleToMPFRMatrix(mpfr_t *ReA, mpfr_t *ImA, complexdouble *A, int m, int n)
{
	int i, j;
	mp_prec_t prec = ceil(sizeof(complexdouble)/2);

	for(i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			prec = (prec >= mpfr_get_prec(ReA[i * n + j]) ? prec : mpfr_get_prec(ReA[i * n + j]));
			mpfr_set_prec(ReA[i * n + j], prec);
			mpfr_set_d(ReA[i * n + j], A[i * n + j].r, MPFR_RNDN);
			
			prec = (prec >= mpfr_get_prec(ImA[i * n + j]) ? prec : mpfr_get_prec(ReA[i * n + j]));
			mpfr_set_prec(ImA[i * n + j], prec);
			mpfr_set_d(ImA[i * n + j], A[i * n + j].i,MPFR_RNDN );
		}
	}

}


void MPFRMatrixToDoublecomplex(mpfr_t *ReA, mpfr_t *ImA, complexdouble *A, int m, int n)
{
	int i, j;

	for(i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			A[i * n + j].r = mpfr_get_d(ReA[i * n + j], MPFR_RNDN);
			A[i * n + j].i = mpfr_get_d(ImA[i * n + j], MPFR_RNDN);

		}
	}

}

/* Convert a MPFR m x n matrix to a double matrix. The function converts element-by-element mpfr floats to double precision
using mpfr_get_d function */
void MPFRMatrixToDouble(double *Ad, mpfr_t *A, uint64_t m, uint64_t n)
{
	int i, j;
	for(i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			Ad[i * n + j] = mpfr_get_d(A[i * n + j], MPFR_RNDN);
		}
	}
}

void MPFIMatrixToDoublecomplex(mpfi_t *ReA, mpfi_t *ImA, complexdouble *A, int m, int n)
{
	int i, j;

	for(i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			A[i * n + j].r = mpfi_get_d(ReA[i * n + j]);
			A[i * n + j].i = mpfi_get_d(ImA[i * n + j]);

		}
	}

}





/* Convert a complex n * m matrix A represented in format clapack complexdouble to the format of two the interval n * m mpfi_t matrices 
ReA and ImA containing real and imaginary parts of A respectively, possibly changing their precision.
Matrices ReA and ImA are assumed to have been allocated outside this function. 
The function uses the MPFR variable scratch as scratch space, assuming they have been allocated, possibly changing its
precision and not clearing it.
*/
void complexdoubleToMPFIMatrix(mpfi_t *ReA, mpfi_t *ImA, complexdouble *A, int m, int n)
{
	int i, j;
	// mp_prec_t prec = ceil(sizeof(complexdouble)/2);

	for(i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			mpfi_set_ui(ReA[i * n + j], 0);
			mpfi_set_ui(ImA[i * n + j], 0);
			mpfi_set_d(ReA[i * n + j], A[i * n + j].r);
			
			mpfi_set_d(ImA[i * n + j], A[i * n + j].i);
		}
	}

}

/* For a complex interval m x n matrix A the function returns its conversion to floating-point.
Output matrix B is assumed to be preallocated outside the function. Its precision may be changed.
 */
void MPFIMatrixtoMPFRMatrix(mpfr_t *reB, mpfr_t *imB, mpfi_t *reA, mpfi_t *imA, uint64_t m, uint64_t n)
{
	int i,j;
	mpfr_prec_t precA;
	mpfr_prec_t precB;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			precA = mpfi_get_prec(reA[i * n + j]);
			precB = mpfr_get_prec(reB[i * n + j]);
			if(precA > precB)
				mpfr_set_prec(reB[i * n + j], precA);
			mpfi_get_fr(reB[i * n + j], reA[i * n + j]);

			precA = mpfi_get_prec(imA[i * n + j]);
			precB = mpfr_get_prec(imB[i * n + j]);
			if(precA > precB)
				mpfr_set_prec(imB[i * n + j], precA);
			mpfi_get_fr(imB[i * n + j], imA[i * n + j]);

		}
	}
}


/* For a complex m x n MPFR matrix A the function returns its conversion to complex interval matrix. The result
matrix B is assumed to be preallocated outside the function. Its precision may be changed. */
void MPFRMatrixToMPFIMatrix(mpfi_t *reB, mpfi_t *imB, mpfr_t *reA, mpfr_t *imA,uint64_t m, uint64_t n)
{
	int i,j;
	mpfr_prec_t precA;
	mpfr_prec_t precB;
	for(i = 0; i < m; ++i)
	{
		for(j = 0; j < n; ++j)
		{
			precA = mpfr_get_prec(reA[i * n + j]);
			precB = mpfi_get_prec(reB[i * n + j]);
			if(precA > precB)
				mpfi_set_prec(reB[i * n + j], precA);
			mpfi_set_fr(reB[i * n + j], reA[i * n + j]);

			precA = mpfr_get_prec(imA[i * n + j]);
			precB = mpfi_get_prec(imB[i * n + j]);
			if(precA > precB)
				mpfi_set_prec(imB[i * n + j], precA);
			mpfi_set_fr(imB[i * n + j], imA[i * n + j]);

		}
	}

}









