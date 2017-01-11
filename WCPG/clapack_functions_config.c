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

#include "clapack_functions_config.h"



void my_zgetrf(int *m, int *n, complexdouble *A,
 	int *lda, int *ipiv, int *info)
{
		
	#ifdef LAPACKE_HEADER		//if we have non-standard clapack wrapper

		//we convert complexdouble *A into lapacke_complex type
		lapacke_complex *newA = (lapacke_complex*)wcpgSafeMalloc((*m) * (*n) * sizeof(lapacke_complex));
		int i,j;
		for(i = 0; i < *m; ++i)
		{
			for(j = 0; j < *n; ++j)
			{
				newA[i*(*n)+j] = A[i*(*n)+j].r + I * A[i*(*n)+j].i;
			}
		}

		//then use the lapacke function
	/*	lapack_int LAPACKE_zgetrf( int matrix_order, lapack_int m, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_int* ipiv ); */

		 *info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, *m, *n, newA, *lda, ipiv);

		//and convert back to complexdouble *A type
		for(i = 0; i < *m; ++i)
				{
					for(j = 0; j < *n; ++j)
					{
						A[i*(*n)+j].r = creal(newA[i*(*n)+j]);
						A[i*(*n)+j].i = cimag(newA[i*(*n)+j]);
					}
				}
		wcpgSafeFree(newA);
	#else
		zgetrf_((integer*)m, (integer*)n, (doublecomplex*)A, (integer*)lda, (integer*)ipiv, (integer*)info);
	#endif

}

void my_zgetri(int *n, complexdouble *A, int *lda,
		int *ipiv, complexdouble *work, int *lwork, int *info)
{
		
	#ifdef LAPACKE_HEADER 	//if we have non-standard clapack wrapper
		
		//we convert complexdouble *A into lapacke_complex type
		lapacke_complex *newA = (lapacke_complex*)wcpgSafeMalloc((*n) * (*n) * sizeof(lapacke_complex));
		
		int i,j;
		for(i = 0; i < *n; ++i)
		{
			for(j = 0; j < *n; ++j)
			{
				newA[i*(*n)+j] = A[i*(*n)+j].r + I * A[i*(*n)+j].i;
			}
		}

		//then use the lapacke function
		/* lapack_int LAPACKE_zgetri( int matrix_order, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           const lapack_int* ipiv ); */

		
		*info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, *n, newA, *lda, ipiv);

		//and convert back to complexdouble *A type
		for(i = 0; i < *n; ++i)
				{
					for(j = 0; j < *n; ++j)
					{
						A[i*(*n)+j].r = creal(newA[i*(*n)+j]);
						A[i*(*n)+j].i = cimag(newA[i*(*n)+j]);
					}
				}
		wcpgSafeFree(newA);
	
	#else
		zgetri_((integer*)n, (doublecomplex*)A, (integer*)lda, (integer*)ipiv, (doublecomplex*)work, (integer*)lwork, (integer*)info);
	#endif



}

void my_dgeevx(int *n, double *A, int *lda, double *wr, double *wi, double *vl,  \
	       int *ldvl, double *vr, int *ldvr, int *ilo, int *ihi,  \
	       double *scale, double *abnrm, double *rconde, double *rcondv,  \
	       double *work, int *lwork, long int *iwork, int *info)

{

	#ifdef LAPACKE_HEADER 	//if we have non-standard clapack wrapper

	/*lapack_int LAPACKE_dgeevx( int matrix_order, char balanc, char jobvl,
                           char jobvr, char sense, lapack_int n, double* a,
                           lapack_int lda, double* wr, double* wi, double* vl,
                           lapack_int ldvl, double* vr, lapack_int ldvr,
                           lapack_int* ilo, lapack_int* ihi, double* scale,
                           double* abnrm, double* rconde, double* rcondv );	*/

		*info = LAPACKE_dgeevx(LAPACK_COL_MAJOR, 'N', 'V', 'V','B',*n, A, *lda, wr, wi, vl, *ldvl, vr, *ldvr, ilo, ihi, scale, abnrm, rconde, rcondv);

	#else
		/* dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, real *a, integer *lda, real *wr, real *wi, real *
	vl, integer *ldvl, real *vr, integer *ldvr, integer *ilo, integer *
	ihi, real *scale, real *abnrm, real *rconde, real *rcondv, real *work, 
	 integer *lwork, integer *iwork, integer *info);*/
		char *jobvr = "V";
		char *jobvl = "V";
		char *balanc = "N";
		char *sense = "B";
		dgeevx_(balanc, jobvl, jobvr, sense,(integer*)n, (doublereal*)A, (integer*)lda, wr, wi, \
				 vl, (integer*)ldvl, vr, (integer*)ldvr, (integer*)ilo, (integer*)ihi,\
				 scale, abnrm, rconde, rcondv, work, (integer*)lwork, (integer*)iwork,\
				 (integer*)info); 
	
	#endif
}


