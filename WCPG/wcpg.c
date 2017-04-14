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

#include "_wcpg.h"
#include "wcpg.h"


/* Swaps two MPFR array pointers */
static inline void swapMPFRArrayPointers(mpfr_t **p, mpfr_t **q) {
  mpfr_t *temp;
  temp = *p;
  *p = *q;
  *q = temp;
}

void wcpg_result_init(wcpg_result *result)
{
	result->N = 0;
	result->time_overall = 0;
	result->time_Ncomp = 0;
	result->time_Summation = 0;
	mpfr_init2(result->one_minus_rhoA, 64);
	mpfr_init(result->maxSN);
	mpfr_init(result->minSN);
	result->inversion_Iter = 0;
	result->maxprec_PN = 0;
	result->maxprec_U = 0;
	result->maxprec_SN = 0;
}

void wcpg_result_clear(wcpg_result *result)
{
	mpfr_clear(result->maxSN);
	mpfr_clear(result->minSN);
	mpfr_clear(result->one_minus_rhoA);
}

void wcpg_result_print(FILE *stream, wcpg_result *result, size_t ndigits)
{
	fprintf(stream, "Truncation order N: %d\n", result->N);
	fprintf(stream, "1 - rho(A): \n");
	mpfr_out_str(stream, (int)10, ndigits, result->one_minus_rhoA, MPFR_RNDN );
	fprintf(stream, "\nInversion iterations: %d \n", result->inversion_Iter );
	fprintf(stream, "max(S_N): \n");
	mpfr_out_str(stream, (int)10, ndigits, result->maxSN, MPFR_RNDN );
	fprintf(stream, "\nmin(S_N): \n");
	mpfr_out_str(stream, (int)10, ndigits, result->minSN, MPFR_RNDN );
	fprintf(stream, "\nOverall time spent: %.2f \n", result->time_overall );
	fprintf(stream, "N computation time: %.2f \n", result->time_Ncomp );
	fprintf(stream, "Summation time: %.2f \n", result->time_Summation );
	fprintf(stream, "Maximum precision of U = inv(V): %ld \n", result->maxprec_U );
	fprintf(stream, "Maximum precision of P_N: %ld \n", result->maxprec_PN );
	fprintf(stream, "Maximum precision of S_N: %ld \n", result->maxprec_SN );

}

/* --------------Functions for N_lower_bound algorithm. To be replaced later. -----------------*/



/* For complex p x n matrix Phi, complex n x q matrix Psi the funciton constructs
l-th residue matrix R, wich is a complex p x q. */
void R_l(mpfi_t *reRl, mpfi_t *imRl, mpfi_t *rePhi, mpfi_t *imPhi, mpfi_t *rePsi, mpfi_t *imPsi, uint64_t n, uint64_t p, uint64_t q, int l)
{
	int i, j;
	// if (l >= n)
	// 	printf("Error! R_(%d)th residue matrix does not exist!\n", l);
	mpfi_t scratch;
	mpfi_init2(scratch, 64*2);
	for( i = 0; i < p; ++i)
	{
		for( j = 0; j < q; ++j)
		{
			mpfi_mul_complex(reRl[i*q + j], imRl[i*q + j], rePhi[i * n + l], imPhi[i * n + l], rePsi[l * q + j], imPsi[l * q + j], scratch);

		}
	}
	mpfi_clear(scratch);

}

/* This funciton coputes the lower bound on truncation order N for WCPG sum. The algorithm for 
this bound is described in details in report(article). Short exmplanation:
Suppose state-space [A,B,C,D] representation of LTI filter is given. 
Then,

N >= ceil{log(eps/norm(M)_min) / log(1 - rho(A))}

with 
	* eps - desired error of truncation;
	* norm(M)_min = min_i,j {M_ij}
	* M = sum_l=1:n {k_l * R_l} is p x q real matrix
	* k_l = abs(lambda_l) / ((1 - abs(lambda_l))*rho) is a scalar coefficient
	* R_l is l-th residue matrix 
 */

int lowerBoundN(mpfi_t *rePhi, mpfi_t *imPhi, mpfi_t *rePsi, mpfi_t *imPsi, \
					 mpfi_t *reLambda, mpfi_t *imLambda,mpfr_t eps, uint64_t n, uint64_t p, uint64_t q, mpfr_prec_t prec, wcpg_result *result)
{


	mpfi_t scratch1, scratch2, scratch3, scratchabs, *scratchVector2;
	mpfi_init2(scratch1, prec);
	mpfi_init2(scratch2, prec);
	mpfi_init2(scratch3, prec);
	mpfi_init2(scratchabs, prec);
	scratchVector2 = allocateMPFIMatrix(1,2,prec);

	mpfi_t rho;
	mpfi_init2(rho, prec);
	mpfi_maxabs(rho, reLambda, imLambda, 1, n, scratchVector2);
		if(!mpfi_bounded_p(rho))
	{
		fprintf(stderr, "ERROR: NaN/Inf occurred, could not compute WCPG. Exit with error --rho\n");
	}

	mpfi_t oneMinusRho;						
	mpfi_init2(oneMinusRho, prec);
	mpfi_ui_sub(oneMinusRho, (unsigned long int)1, rho);

	mpfr_t rhofr;
	mpfr_init(rhofr);
	mpfi_get_fr(rhofr, oneMinusRho);		// converting interval to floating point
	mpfr_set(result->one_minus_rhoA, rhofr, MPFR_RNDU);
	if(!mpfi_bounded_p(oneMinusRho))
	{
		fprintf(stderr, "ERROR: NaN/Inf occurred, could not compute WCPG. Exit with error\n");
	}



	mpfi_t denum;				
	mpfi_init2(denum, prec);
	mpfi_log2(denum, rho);		
	mpfi_abs(denum, denum);		//denum = abs(log2(rho));


	if(!mpfi_bounded_p(denum))
	{
		fprintf(stderr, "ERROR: NaN/Inf occurred, could not compute WCPG. Exit with error\n");
	}

	//We need to compute M = sum_1^n (K_i* R_i) with K = rho*abs(lambda) / (1 - abs(lambda))
	mpfi_t K;
	mpfi_init2(K, prec);

	mpfi_t *reRl, *imRl;
	reRl = allocateMPFIMatrix(p, q, prec);
	imRl = allocateMPFIMatrix(p, q, prec);

	mpfi_t Mij, min;
	mpfi_init2(min, prec);
	mpfi_init2(Mij, prec);
	mpfi_set_ui(min, (unsigned long int)INT_MAX);

	int l, i, j;
	for(l = 0; l < n; ++l)
	{
		mpfi_abs_complex(scratch1, reLambda[l], imLambda[l], scratchabs);    //scratch1 = abs(lambda_l)
		mpfi_ui_sub(scratch2, 1, scratch1); 						//scratch2 = 1 - abs(lambda_l);
		mpfi_mul(scratch3, scratch2, rho); 						// scratch3 = (1 - abs(lambda_l))*rho
		mpfi_div(K, scratch1, scratch3); 							// K = abs(lambda_l) / ((1 - abs(lambda_l))*rho)
		R_l(reRl, imRl, rePhi, imPhi, rePsi, imPsi, n, p, q, l);		//l-th residue matrix 
		
		for(i = 0; i < p; ++i)
		{
			for(j = 0; j < q; ++j)
			{
				
				mpfi_abs_complex(scratch1, reRl[i * q + j], imRl[i * q + j], scratchabs);	//scratch1 = abs(R_l(i,j))


				mpfi_mul(Mij, K, scratch1);		// M(i,j) = K * abs(R_l(i,j))

				if(mpfi_cmp(Mij, min) < 0 && mpfi_get_d(Mij))
				{
					mpfi_set(min, Mij);

				}

			}
		}
	}

	// so now we know the norm_min of matrix M: min = norm(M)_min;
	mpfi_fr_div(scratch1, eps, min);  //scratch1 = eps / norm(M)_min
	mpfi_log2(scratch1, scratch1);	// scratch1 = log2(eps / norm(M)_min)



	mpfi_div(scratch2, scratch1, denum); //scratch2 = scratch1 / denum = log2(eps / norm(M)_min) / log2(rho(A))
	mpfi_abs(scratch2, scratch2);



	mpfr_t Nsup;
	mpfr_init(Nsup);
	mpfi_get_right(Nsup, scratch2);


	int N = ceil(mpfr_get_d(Nsup, MPFR_RNDN));

	mpfi_clear(scratch1);
	mpfi_clear(scratch2);
	mpfi_clear(scratch3);
	mpfi_clear(scratchabs);
	mpfi_clear(rho);
	mpfi_clear(oneMinusRho);
	mpfi_clear(denum);
	mpfi_clear(K);
	mpfi_clear(Mij);
	mpfi_clear(min);
	mpfr_clear(Nsup);
	mpfr_clear(rhofr);
	freeMPFIMatrix(scratchVector2, 1, 2);
	freeMPFIMatrix(reRl, p, q);
	freeMPFIMatrix(imRl, p, q);

	return N;
}



/* ---------------- End of N_lower_bound algorithm functions ---------*/



/* THe function determines the minimal and maximum values of the
computed WCPG matrix and adds them to the corresonding fields
of the result structure. 
 */
void getDeltaOfWCPG(mpfr_t *reS, uint64_t p, uint64_t q, wcpg_result *result)
{
	int i,j;

	mpfr_prec_t e = getMaxPrecision(reS, reS, p, q);

	mpfr_t min;
	mpfr_init2(min, e);
	mpfr_set(min, reS[0], MPFR_RNDN);

	mpfr_t max;
	mpfr_init2(max, e);
	mpfr_set(max, reS[0], MPFR_RNDN);

	for( i = 0; i < p; ++i)
	{
		for(j = 0; j < q; ++j)
		{
			if(mpfr_cmp(reS[i * q + j],min) < 0 )
			{
				mpfr_set(min, reS[i * q + j], MPFR_RNDN);

			}
			if(mpfr_cmp(reS[i * q + j],max) > 0 )
				mpfr_set(max, reS[i * q + j], MPFR_RNDN);
		}
	}

	mpfr_set_prec(result->maxSN, e);
	mpfr_set_prec(result->minSN, e);

	mpfr_set(result->maxSN, max, MPFR_RNDN);
	mpfr_set(result->minSN, min, MPFR_RNDU);

	mpfr_clear(min);
	mpfr_clear(max);
}



/* For an LTI filter given in its State-Space representation {A,B,C,D},
where A is n*n, B is n*q, C is p*n and D is p*q real matrix the function
returns integer value indicating if WCPG is successfully computed.
The function takes eps, a desired absolute error bound on the computed WCPG measure.
In p*q MPFR matrix W the Worst-Case peak gain is stored if algorithm successfully exited. */
int WCPG_ABCD_mprec(mpfr_t *W, double *A, double *B, double *C, double *D, uint64_t n, uint64_t p, uint64_t q, mpfr_t mpeps)
{
	mpfr_prec_t prec = 106;
	int i;


	if(matrixIsNan_double(A, n, n))
		{
			fprintf(stderr, "Matrix A has a NaN value. Could not compute WCPG. Exit with error.\n");
			return 0;
		}
		if(matrixIsNan_double(B, n, q))
		{
			fprintf(stderr, "Matrix B has a NaN value. Could not compute WCPG. Exit with error.\n");
			return 0;
		}
		if(matrixIsNan_double(C, p, n))
		{
			fprintf(stderr, "Matrix C has a NaN value. Could not compute WCPG. Exit with error.\n");
			return 0;
		}
		if(matrixIsNan_double(D, p, q))
		{
			fprintf(stderr, "Matrix D has a NaN value. Could not compute WCPG. Exit with error.\n");
			return 0;
		}

		/* Convertion of the double matrices into multiprecision complex format */

		mpfr_t *reA;
		reA = allocateMPFRMatrix(n, n, prec);
		doublerealToMPFRMatrix(reA, A, n, n);

		mpfr_t *reB;
		reB = allocateMPFRMatrix(n, q, prec);
		doublerealToMPFRMatrix(reB, B, n, q);

		mpfr_t *reC;
		reC = allocateMPFRMatrix(p, n, prec);
		doublerealToMPFRMatrix(reC, C, p, n);

		mpfr_t *reD;
		reD = allocateMPFRMatrix(p, q, prec);
		doublerealToMPFRMatrix(reD, D, p, q);

			/*
			writeMPFRMatrix(stderr, reA, n, n, 20, MPFR_RNDN);
			writeMPFRMatrix(stderr, reB, n, q, 20, MPFR_RNDN);
			writeMPFRMatrix(stderr, reC, p, n, 20, MPFR_RNDN);
			writeMPFRMatrix(stderr, reD, p, q, 20, MPFR_RNDN);
		*/


		if(matrixIsNan_mpfr(reA, n, n))
		{
			fprintf(stderr, "MPFR Matrix A has a NaN value. Could not compute WCPG. Exit with error.\n");
			freeMPFRMatrix(reA, n, n);
			freeMPFRMatrix(reB, n, q);
			freeMPFRMatrix(reC, p, n);
			freeMPFRMatrix(reD, p, q);

			return 0;
		}
		if(matrixIsNan_mpfr(reB, n, q))
		{
			fprintf(stderr, "MPFR Matrix B has a NaN value. Could not compute WCPG. Exit with error.\n");
			freeMPFRMatrix(reA, n, n);
			freeMPFRMatrix(reB, n, q);
			freeMPFRMatrix(reC, p, n);
			freeMPFRMatrix(reD, p, q);

			return 0;
		}
		if(matrixIsNan_mpfr(reC, p, n))
		{
			fprintf(stderr, "MPFR Matrix C has a NaN value. Could not compute WCPG. Exit with error.\n");
			freeMPFRMatrix(reA, n, n);
			freeMPFRMatrix(reB, n, q);
			freeMPFRMatrix(reC, p, n);
			freeMPFRMatrix(reD, p, q);

			return 0;
		}
		if(matrixIsNan_mpfr(reD, p, q))
		{
			fprintf(stderr, "MPFR Matrix D has a NaN value. Could not compute WCPG. Exit with error.\n");
			freeMPFRMatrix(reA, n, n);
			freeMPFRMatrix(reB, n, q);
			freeMPFRMatrix(reC, p, n);
			freeMPFRMatrix(reD, p, q);

			return 0;
		}

		wcpg_result result;


		int res = WCPG(W, reA, reB, reC, reD, mpeps, n, p, q, &result);

		wcpg_result_clear(&result);
		freeMPFRMatrix(reA, n, n);
		freeMPFRMatrix(reB, n, q);
		freeMPFRMatrix(reC, p, n);
		freeMPFRMatrix(reD, p, q);
		return res;


}

int WCPG_mp(mpfr_t *Sk, mpfr_t *A, mpfr_t *B, mpfr_t *C, mpfr_t *D, uint64_t n, uint64_t p, uint64_t q, mpfr_t mpeps)
{
	wcpg_result result;

	int res = WCPG(Sk, A, B, C, D, mpeps, n, p, q, &result);

	wcpg_result_clear(&result);

	return res;

}

/* For an LTI filter given in its State-Space representation {A,B,C,D},
where A is n*n, B is n*q, C is p*n and D is p*q real matrix,
and for an eps>0 the function returns a multi-precision n*q real matrix Sk of Worst-Case Peak Gains
of the system such that the overall error of computation is bounded by eps. */
int WCPG(mpfr_t *Sk, mpfr_t *A, mpfr_t *B, mpfr_t *C, mpfr_t *D, mpfr_t mpeps, uint64_t n, uint64_t p, uint64_t q, wcpg_result *result)
{


	mpfr_prec_t prec = 106;
	int i;

	//wcpg_result *result;
	wcpg_result_init(result);
	mpfr_t *imA;
	imA = allocateMPFRMatrix(n, n, prec);
	setMatrixZero(imA, n, n);

	mpfr_t *imB;
	imB = allocateMPFRMatrix(n, q, prec);
	setMatrixZero(imB, n, q);

	mpfr_t *imC;
	imC = allocateMPFRMatrix(p, n, prec);
	setMatrixZero(imC, p, n);

	mpfr_t *imD;
	imD = allocateMPFRMatrix(p, q, prec);
	setMatrixZero(imD, p, q);

	clock_t begin, end, begin_Ncomp, end_Ncomp, begin_summation, end_summation;



	// mpfr_t *Wlow;
	// Wlow = allocateMPFRMatrix(p, q, prec);
	// if(!WCPGLowerBound(Wlow, reA, reB, reC, reD, n, p, q, -64))
	// {
	// 	printf("Inside WCPG function: Could not compute lower bound on the Worst Case Peak Gain matrix.\n");
	// 	freeMPFRMatrix(Wlow, p, q);
	// }
	// else
	// {
	// 	printf("Inside WCPG function: Lower bound on the WCPG matrix: \n");
	// 	writeMPFRMatrix(stderr, Wlow, p, q, 20, MPFR_RNDN);
	// 	freeMPFRMatrix(Wlow, p, q);
	// }



	/* starting the timer */

	begin = clock();				//we start timer here!
	begin_Ncomp = clock();


	/* Computing eps_1 = eps/6 */

	mpfr_t eps_1;
	mpfr_init2(eps_1, prec);				//eps_1 = eps / 6
	mpfr_set(eps_1, mpeps, MPFR_RNDD);
	mpfr_div_ui(eps_1, eps_1, 6, MPFR_RNDD);	


	/* Preparing scratch variables*/

	mpfi_t *scratch3intval;
	scratch3intval = allocateMPFIMatrix(1,3, prec);

	mpfr_t *scratch3;
	scratch3 = allocateMPFRVector(3, 64);

	mpfr_t *scratchT;
	scratchT = allocateMPFRVector(2 * n + 1, 64);

	mpfr_t scratch, scratch2;
	mpfr_init2(scratch, prec);
	mpfr_init2(scratch2, prec);


	/* -------------Step 1. Computing lower bound on N --------------*/

	/* If we want to separate this algorithm, then
	we must define local context structure for:
	1. Input: A,B,C,D, result, prec, n, p, q,mpeps
	2. Output: Nmax, reV, imV, reVint, imVint, Vinv

	*/

	
	/*Finding eigenvalues and eigenvectors*/

	
	// Variables for interval matrices of eigenvectors and eigenvalues
	mpfi_t *reVint, *imVint;
	reVint = allocateMPFIMatrix(n, n, prec);
	imVint = allocateMPFIMatrix(n, n, prec);

	mpfi_t *reLambdaint, *imLambdaint;
	reLambdaint = allocateMPFIMatrix(1, n, prec);
	imLambdaint = allocateMPFIMatrix(1, n, prec);


	if(!MPFItrustedEigendecompsition(reLambdaint, imLambdaint, reVint, imVint, A, n))
	{
		freeMPFRMatrix(imA, n, n);
		freeMPFRMatrix(imB, n, q);
		freeMPFRMatrix(imC, p, n);
		freeMPFRMatrix(imD, p, q);
		freeMPFIMatrix(reVint, n, n);
		freeMPFIMatrix(imVint, n, n);
		freeMPFIMatrix(reLambdaint, 1, n);
		freeMPFIMatrix(imLambdaint, 1, n);
		fprintf(stderr, "ERROR: could not compute the WCPG \n");
		return 0;
	}

	/* Preparing variables for inversion of V */
	mpfr_t *reV;
	mpfr_t *imV;
	reV = allocateMPFRMatrix(n, n, prec);
	imV = allocateMPFRMatrix(n, n, prec);
	MPFIMatrixtoMPFRMatrix(reV, imV, reVint, imVint, n, n);


	/* Compute inv(V) with clapack and check if we can use our inverse method */
	complexdouble *Vinv = (complexdouble*)wcpgSafeMalloc(n*n*sizeof(complexdouble));
	complexdouble *Vcpy = (complexdouble*)wcpgSafeMalloc(n*n*sizeof(complexdouble));		// we copy the matrix, because CLAPACK changes it in eigensolver

	MPFIMatrixToDoublecomplex(reVint, imVint, Vcpy, n, n);
	if(!clapack_complex_matrix_inverse(Vinv, Vcpy, n))
	{
		freeMPFRMatrix(imA, n, n);
		freeMPFRMatrix(imB, n, q);
		freeMPFRMatrix(imC, p, n);
		freeMPFRMatrix(imD, p, q);
		freeMPFIMatrix(reVint, n, n);
		freeMPFIMatrix(imVint, n, n);
		freeMPFIMatrix(reLambdaint, 1, n);
		freeMPFIMatrix(imLambdaint, 1, n);
		freeMPFRMatrix(reV, n, n);
		freeMPFRMatrix(imV, n, n);
		wcpgSafeFree(Vinv);
		wcpgSafeFree(Vcpy);

		fprintf(stderr, "ERROR: could not compute the WCPG \n");
		return 0;
	}


	if(matrixIsNan_complexdouble(Vinv, n, n))
	{
		freeMPFRMatrix(imA, n, n);
		freeMPFRMatrix(imB, n, q);
		freeMPFRMatrix(imC, p, n);
		freeMPFRMatrix(imD, p, q);
		freeMPFIMatrix(reVint, n, n);
		freeMPFIMatrix(imVint, n, n);
		freeMPFIMatrix(reLambdaint, 1, n);
		freeMPFIMatrix(imLambdaint, 1, n);
		freeMPFRMatrix(reV, n, n);
		freeMPFRMatrix(imV, n, n);
		wcpgSafeFree(Vinv);
		wcpgSafeFree(Vcpy);
		fprintf(stderr, "Complex Matrix inv(V) has a NaN value. Could not compute WCPG. Exit with error.\n");
		return 0;
	}

	/* Transform inv(V) into MPFR format matrix */
	mpfr_t *reVinv;
	mpfr_t *imVinv;
	reVinv = allocateMPFRMatrix(n,n, prec);
	imVinv = allocateMPFRMatrix(n,n, prec);
	complexdoubleToMPFRMatrix(reVinv, imVinv, Vinv, n,n);


	/* Check if we can inverse. 
	To update after frobeniusNormUpperBound basic brick implemented*/

	
	 //if(!inverseConditionsTest())
	 //{
	 //	printf("Could not compute WCPG, test for inversion possibility failed. Algorithm exits.\n");
	 //	return 0;
	 //}

	

	/* Calculating Psi = inv(V) * B */      

  	/* Here we inverse matrix V with just twice the precision, that would CLAPACK use, and after we inverse it with larger precision when compute matrix T.. */
  	/* Further we would like to do it only once, with max precision. */
	mpfr_t *reU1;				
	mpfr_t *imU1;
	imU1 = allocateMPFRMatrix(n,n, prec);
	reU1 = allocateMPFRMatrix(n,n, prec);
	setMatrixZero(imU1, n,n);
	setMatrixZero(reU1, n,n);

	invertComplexMatrix(reU1, imU1, reV, imV, reVinv, imVinv, n, -prec);


	mpfi_t *reU1_int, *imU1_int;
    reU1_int = allocateMPFIMatrix(n, n, prec);
    imU1_int = allocateMPFIMatrix(n, n, prec);
    MPFRMatrixToMPFIMatrix(reU1_int, imU1_int, reU1, imU1, n, n);	//convert computed inverse to mpfi

    freeMPFRMatrix(reU1, n, n);
    freeMPFRMatrix(imU1, n, n);

    /* Construct Psi = inv(V) * B */
	mpfi_t *rePsi, *imPsi;
    rePsi = allocateMPFIMatrix(n, q, prec);
    imPsi = allocateMPFIMatrix(n, q, prec);


    MPFIComplexMatrixMultiplyMPFRComplexMatrix(rePsi, imPsi, reU1_int, imU1_int, B, imB, n, n, q, scratch3intval);	//Psi = inv(V) * B

  	/* Computing Phi = C * V */
    mpfi_t *rePhi, *imPhi;
   	rePhi = allocateMPFIMatrix(p, n, prec);
   	imPhi = allocateMPFIMatrix(p, n, prec);

   	MPFRComplexMatrixMultiplyMPFIComplexMatrix(rePhi, imPhi, C, imC, reVint, imVint, p, n, n, scratch3intval);		//[Phi] = C*[V];


    /* Compute lower bound on N using formula (28) [see paper] */
    int Nmax = lowerBoundN(rePhi, imPhi, rePsi, imPsi, reLambdaint, imLambdaint,eps_1, n, p, q, prec, result);
    //fprintf(stderr, "Nmax = %d\n", Nmax);
  	if(Nmax <= 0)	//then we have a NaN
  	{

		freeMPFRMatrix(imA, n, n);
		freeMPFRMatrix(imB, n, q);
		freeMPFRMatrix(imC, p, n);
		freeMPFRMatrix(imD, p, q);
		freeMPFIMatrix(reVint, n, n);
		freeMPFIMatrix(imVint, n, n);
		freeMPFIMatrix(reLambdaint, 1, n);
		freeMPFIMatrix(imLambdaint, 1, n);
		freeMPFRMatrix(reV, n, n);
		freeMPFRMatrix(imV, n, n);
		wcpgSafeFree(Vinv);
		wcpgSafeFree(Vcpy);
		freeMPFIMatrix(imU1_int, n, n);
		freeMPFIMatrix(reU1_int, n, n);
		freeMPFIMatrix(rePsi, n, q);
		freeMPFIMatrix(imPsi, n, q);
		freeMPFIMatrix(rePhi, p, n);
		freeMPFIMatrix(imPhi, p, n);
  		fprintf(stderr, "ERROR: could not determine the lower bound on truncation order N. N = %d\n", Nmax);
  		return 0;
  	}

    result->N = Nmax;

    end_Ncomp = clock();



	/* Freeing complexdouble and doublereal matrices */

	/* TO BE COMPLETED */


	
	/* -------------------Summation stage------------------------ */



	begin_summation = clock();
	/* prepare Zero Matrices of different sizes that will be used for multiplyAndADDMPFRMatrices */
	mpfr_t *In;
	In = allocateMPFRMatrix(n,n, 64);
	setMatrixIdentity(In, n);

	mpfr_t *Zeronn;
	Zeronn = allocateMPFRMatrix(n,n, 64);
	setMatrixZero(Zeronn, n, n);

	mpfr_t *Zeropn;
	Zeropn = allocateMPFRMatrix(p,n, 64);
	setMatrixZero(Zeropn, p, n);

	mpfr_t *Zeronq;
	Zeronq = allocateMPFRMatrix(n,q, 64);
	setMatrixZero(Zeronq, n, q);

	mpfr_t *Zeropq;
	Zeropq = allocateMPFRMatrix(p,q, 64);
	setMatrixZero(Zeropq,p,q);

	/* Step 2. Constructing matrix T = inv(V)*A*V */

	/* Prepare variables for our inversion method */
	mpfr_t *reU;
	mpfr_t *imU;
	imU = allocateMPFRMatrix(n,n, prec);
	reU = allocateMPFRMatrix(n,n, prec);
	setMatrixZero(imU, n,n);
	setMatrixZero(reU, n,n);



	//Compute upper bounds on Frobenius norms of our matrices V, Vinv, B, C

	mpfr_t frobV;
	mpfr_init2(frobV, prec);
	frobeniusNormUpperBound(frobV, reV, imV, n, n, scratch);

	mpfr_t frobC;
	mpfr_init2(frobC, prec);
	frobeniusNormUpperBound(frobC, C, imC, p, n, scratch);

	mpfr_t frobB;
	mpfr_init2(frobB, prec);
	frobeniusNormUpperBound(frobB, B, imB, n, q, scratch);

	mpfr_t frobVinv;				//Frobenius norm of the inverse computed with CLAPACK
	mpfr_init2(frobVinv, prec);
	frobeniusNormUpperBound(frobVinv, reVinv, imVinv, n, n, scratch);


	/* Precisions for each step are computed via following formulas: */
	/*
		Step 2. eps_2 <= eps / { 6 * sqrt(n) * (Nmax + 1) * (Nmax + 2) * ||C||_F * ||V||_F * ||Vinv||_F * ||B||_F  }

		Step 3. eps_3 = min(eps_3C, eps_3B) = min {eps/(3 * sqrt(n) * 6 * (Nmax + 1) * ||C||_F * ||V||_F), 
			eps/(3 * sqrt(n) * 6 * (Nmax + 1) * ||B||_F * ||Vinv||_F)}

		Step 4. eps_4_k <= eps / {sqrt(n) * (Nmax + 1) * (Nmax -1) * 6 * ||C'||_F * ||B'||_F}

		Step 5. eps_5_k <= eps / {6 * (Nmax + 1)}

		Step 6. eps_6_k <= eps / {6 * (Nmax + 1)}

	 */

	mpfr_t sqrtn;
	mpfr_init2(sqrtn, prec);
	mpfr_sqrt_ui(sqrtn, n, MPFR_RNDU);


	mpfr_t eps_2;
	mpfr_init2(eps_2, prec);							//eps_2 = eps
	mpfr_set(eps_2, mpeps, MPFR_RNDD);
	mpfr_div(eps_2, eps_2, sqrtn, MPFR_RNDD);			//eps_2 = eps / sqrt(n)
	uint64_t beta = 6 * (Nmax + 1) * (Nmax + 2);	
	mpfr_div_ui(eps_2, eps_2, beta, MPFR_RNDD);			//eps_2 = eps / 6 * sqrt(n) * (Nmax + 1) * (Nmax + 2) 
	mpfr_div(eps_2, eps_2, frobVinv, MPFR_RNDD);
	mpfr_div(eps_2, eps_2, frobC, MPFR_RNDD);
	mpfr_div(eps_2, eps_2, frobB, MPFR_RNDD);
	mpfr_div(eps_2, eps_2, frobV, MPFR_RNDD);			//eps_2 = eps / { 6 * sqrt(n) * (Nmax + 1) * (Nmax + 2) * ||C||_F * ||V||_F * ||Vinv||_F * ||B||_F
	mpfr_exp_t eps_2_exp = mpfr_get_exp(eps_2) - 5;


	mpfr_t eps_3C;
	mpfr_init2(eps_3C, prec);				//eps_3C = eps
	mpfr_set(eps_3C, mpeps, MPFR_RNDD);
	mpfr_div(eps_3C, eps_3C, sqrtn, MPFR_RNDD);			//eps_3C = eps / sqrt(n)
	beta = 3 * 6 * (Nmax + 1);	
	mpfr_div_ui(eps_3C, eps_3C, beta, MPFR_RNDD);		//eps_3C = eps / 3* 6 * sqrt(n) * (Nmax + 1) 
	mpfr_div(eps_3C, eps_3C, frobC, MPFR_RNDD);
	mpfr_div(eps_3C, eps_3C, frobV, MPFR_RNDD);			//eps_3C = eps/(3 * sqrt(n) * 6 * (Nmax + 1) * ||C||_F * ||V||_F
	mpfr_exp_t eps_3C_exp = mpfr_get_exp(eps_3C) - 5;

	mpfr_t eps_3B;
	mpfr_init2(eps_3B, prec);				//eps_3C = eps
	mpfr_set(eps_3B, mpeps, MPFR_RNDD);
	mpfr_div(eps_3B, eps_3B, sqrtn, MPFR_RNDD);			//eps_3B = eps / sqrt(n)
	beta = 3 * 6 * (Nmax + 1);	
	mpfr_div_ui(eps_3B, eps_3B, beta, MPFR_RNDD);		//eps_3B = eps / 3* 6 * sqrt(n) * (Nmax + 1) 
	mpfr_div(eps_3B, eps_3B, frobB, MPFR_RNDD);
	mpfr_div(eps_3C, eps_3C, frobVinv, MPFR_RNDD);		//eps_3B = eps/(3 * sqrt(n) * 6 * (Nmax + 1) * ||B||_F * ||Vinv||_F
	mpfr_exp_t eps_3B_exp = mpfr_get_exp(eps_3B) - 5;


	mpfr_t eps_5;
	mpfr_init2(eps_5, prec);				//eps_5 = eps / {6 * (Nmax + 1)}
	mpfr_set(eps_5, mpeps, MPFR_RNDD);
	beta = 6 * (Nmax + 1);	
	mpfr_div_ui(eps_5, eps_5, beta, MPFR_RNDD);	
	mpfr_exp_t eps_5_exp = mpfr_get_exp(eps_5) - 5;

	mpfr_t eps_6;							//eps_6 = eps_5
	mpfr_init2(eps_6, prec);				
	mpfr_set(eps_6, eps_5, MPFR_RNDD);
	mpfr_exp_t eps_6_exp = mpfr_get_exp(eps_6) - 5;

	int invIterations = invertComplexMatrix(reU, imU, reV, imV, reVinv, imVinv, n, eps_2_exp);
	result->inversion_Iter = invIterations;

	/* Step 2 */

	// Construct T = inv(V)*A*V

	//Compute Ttmp = A*V
	mpfr_t *reTtmp, *imTtmp;
	reTtmp = allocateMPFRMatrix(n,n, 64);
	imTtmp = allocateMPFRMatrix(n,n, 64);
	setMatrixZero(Zeronn, n, n);
	complexMatrixMultiplyAndAdd(reTtmp, imTtmp, 1, A, imA, reV, imV, 1, Zeronn, Zeronn, n,n,n, eps_2_exp, scratchT);

	//Compute T = inv(V) * Ttmp
	mpfr_t *reT, *imT;
	reT = allocateMPFRMatrix(n,n, 64);
	imT = allocateMPFRMatrix(n,n, 64);
	complexMatrixMultiplyAndAdd(reT, imT, 1, reU, imU, reTtmp, imTtmp, 1, Zeronn, Zeronn, n,n,n, eps_2_exp, scratchT);

	/* Here we check if norm(T)_2 <= 1 */
	int normTcheck = checkTwoNorm(reT, imT, n);
	if(!normTcheck)
	{
		fprintf(stderr, "Error! norm(T)_2 > 1. Algorithms exit. \n");
		return 0;
	}

	/* Step 3 */

	//Compute C' = C * V
	mpfr_t *reCV, *imCV;
	reCV = allocateMPFRMatrix(p,n, 64);
	imCV = allocateMPFRMatrix(p,n, 64);
	complexMatrixMultiplyAndAdd(reCV, imCV, 1, C, imC, reV, imV, -1, Zeropn, Zeropn, p, n, n,eps_3C_exp, scratchT);

	//Compute B' = inv(V) * B
	mpfr_t *reBV, *imBV;
	reBV = allocateMPFRMatrix(n,q, 64);
	imBV = allocateMPFRMatrix(n,q, 64);
	complexMatrixMultiplyAndAdd(reBV, imBV, 1, reU, imU, B, imB, -1, Zeronq, Zeronq, n, n, q, eps_3B_exp, scratchT);

	//Compute Frobenius norms of C' and B'
	mpfr_t frobCV;
	mpfr_init2(frobCV, prec);
	frobeniusNormUpperBound(frobCV, reCV, imCV, p, n, scratch);

	mpfr_t frobBV;
	mpfr_init2(frobBV, prec);
	frobeniusNormUpperBound(frobBV, reBV, imBV, n, q, scratch);


	/* Prepare variables for summation */

	//Compute Sk = abs(D)
	absMPFRMatrix(Sk, D, p, q);

	// Set Pk = I

	mpfr_t *rePk, *imPk;
	rePk = allocateMPFRMatrix(n, n, 64);
	imPk = allocateMPFRMatrix(n, n, 64);
	setMatrixIdentity(rePk, n);
	setMatrixZero(imPk, n, n);

	// Set matrix P_k+1 to zero matrix

	mpfr_t *rePkp1, *imPkp1;
	rePkp1 = allocateMPFRMatrix(n, n, 64);
	imPkp1 = allocateMPFRMatrix(n, n, 64);
	setMatrixZero(rePkp1, n, n);
	setMatrixZero(imPkp1, n, n);

	//Set intermidiate matrix Ltmp to zero matrix

	mpfr_t *reLtmp, *imLtmp;
	reLtmp = allocateMPFRMatrix(n, q, 64);
	imLtmp = allocateMPFRMatrix(n, q, 64);
	setMatrixZero(reLtmp, n, q);
	setMatrixZero(imLtmp, n, q);


	/* Steps 4 - 5 - 6 */

	//Compute eps_4 for powerting matrix T
	mpfr_t eps_4;
	mpfr_init2(eps_4, prec);							//eps_4 = eps
	mpfr_set(eps_4, mpeps, MPFR_RNDD);
	mpfr_div(eps_4, eps_4, sqrtn, MPFR_RNDD);			//eps_4 = eps / sqrt(n)
	beta = 6 * (Nmax + 1) * (Nmax -1);	
	mpfr_div_ui(eps_4, eps_4, beta, MPFR_RNDD);			//eps_2 = eps / 6 * sqrt(n) * (Nmax + 1) * (Nmax -1) 
	mpfr_div(eps_4, eps_4, frobCV, MPFR_RNDD);
	mpfr_div(eps_4, eps_4, frobBV, MPFR_RNDD);			//eps_2 = eps / { 6 * sqrt(n) * (Nmax + 1) * (Nmax + 2) * ||C||_F * ||V||_F * ||Vinv||_F * ||B||_F
	mpfr_exp_t eps_4_exp = mpfr_get_exp(eps_4) - 5;

	//Compute sum for k = 0 before the main sycle: S_0 = abs(D) + abs(C'*B')

	// Set L = C' x B' 
	
	mpfr_t *reL, *imL;
	reL = allocateMPFRMatrix(p,q, 64);
	imL = allocateMPFRMatrix(p,q, 64);
	complexMatrixMultiplyAndAdd(reL, imL, 1, reCV, imCV, reBV, imBV, 1, Zeropq, Zeropq, p, n, q, eps_5_exp, scratchT);

	// Accumulate Sk = Sk + abs(Lk)
	accumulateAbsoluteValue(Sk, Sk, reL, imL, p, q, eps_6_exp, scratch3);
	int k;
	for(k = 1; k < Nmax + 1; ++k)
	{

		// P_k+1 = T * P_k
		complexMatrixMultiplyAndAdd(rePkp1, imPkp1, 1, reT, imT, rePk, imPk, 1, Zeronn, Zeronn, n, n, n, eps_4_exp, scratchT); // P_(k+1) = T x P_k

		// Ltmp = P_k+1 * ùùùùùùC'
		complexMatrixMultiplyAndAdd(reLtmp, imLtmp, 1, rePkp1, imPkp1, reBV, imBV, 1, Zeronq, Zeronq, n, n, q, eps_5_exp, scratchT);	
	
		// L_k+1 = B' * Ltmp
		complexMatrixMultiplyAndAdd(reL, imL, 1, reCV, imCV, reLtmp, imLtmp,1, Zeropq, Zeropq, p, n, q, eps_5_exp, scratchT);

		//S_k+1 = S_k + abs(L_k+1)
		accumulateAbsoluteValue(Sk, Sk, reL, imL, p, q, eps_6_exp, scratch3);

		// P_k = P_k+1
		swapMPFRArrayPointers(&rePk, &rePkp1);
    	swapMPFRArrayPointers(&imPk, &imPkp1);
	}

	end_summation = clock();
	end = clock();

	result->time_overall = (double)(end - begin) / CLOCKS_PER_SEC;
	result->time_Summation = (double)(end_summation - begin_summation) / CLOCKS_PER_SEC;
	result->time_Ncomp = (double)(end_Ncomp - begin_Ncomp) / CLOCKS_PER_SEC;
	
	mpfr_prec_t maxPN = getMaxPrecision(rePk, imPk, n, n);
	mpfr_prec_t maxSN = getMaxPrecision(Sk, Zeropq, p, q);
	mpfr_prec_t maxU = getMaxPrecision(reU, imU, n, n);

	result->maxprec_U = maxU;
	result->maxprec_PN = maxPN;
	result->maxprec_SN = maxSN;

	getDeltaOfWCPG(Sk, p, q, result);

	int j;
	for(i=0; i < p; ++i)
	{
		for(j = 0; j < q; ++j)
		{
			if(mpfr_nan_p(Sk[i * q + j]))
			{
				fprintf(stderr, "One of the WCPG matrix elements is NaN. Exit with error. \n");
				return 0;
			}
		}
	}

	//wcpg_result_clear(result);
	freeMPFRMatrix(Zeropn, p, n);
	freeMPFRMatrix(Zeronn, n, n);
	freeMPFRMatrix(Zeronq, n, q);
	freeMPFRMatrix(Zeropq, p, q);
	freeMPFRMatrix(reV, n, n);
	freeMPFRMatrix(imV, n, n);

	freeMPFRMatrix(imA, n, n);
	freeMPFRMatrix(imB, n, q);
	freeMPFRMatrix(imC, p, n);
	freeMPFRMatrix(imD, p, q);
	freeMPFRMatrix(reCV, p, n);
	freeMPFRMatrix(imCV, p, n);
	freeMPFRMatrix(reBV, n, q);
	freeMPFRMatrix(imBV, n, q);
	freeMPFRMatrix(reU, n, n);
	freeMPFRMatrix(imU, n, n);
	freeMPFRMatrix(reT, n, n);
	freeMPFRMatrix(imT, n, n);
	freeMPFRMatrix(reTtmp, n, n);
	freeMPFRMatrix(imTtmp, n, n);
	freeMPFRMatrix(rePk, n, n);
	freeMPFRMatrix(imPk, n, n);
	freeMPFRMatrix(rePkp1, n, n);
	freeMPFRMatrix(imPkp1, n, n);
	freeMPFRMatrix(reL, p, q);
	freeMPFRMatrix(imL, p, q);
	freeMPFRMatrix(reLtmp, n, q);
	freeMPFRMatrix(imLtmp, n, q);
	freeMPFRMatrix(In, n, n);

	freeMPFIMatrix(scratch3intval, 1, 3);
	freeMPFRMatrix(scratch3, 1, 3);
	freeMPFRMatrix(scratchT, 2 * n + 1, 1);
	freeMPFIMatrix(reVint, n, n);
	freeMPFIMatrix(imVint, n, n);
	freeMPFIMatrix(reLambdaint, 1, n);
	freeMPFIMatrix(imLambdaint, 1, n);

	freeMPFRMatrix(reVinv, n, n);
	freeMPFRMatrix(imVinv, n, n);

	freeMPFIMatrix(reU1_int, n, n);
	freeMPFIMatrix(imU1_int, n, n);

	freeMPFIMatrix(rePsi, n, q);
	freeMPFIMatrix(imPsi, n, q);
	freeMPFIMatrix(rePhi, p, n);
	freeMPFIMatrix(imPhi, p, n);



	wcpgSafeFree(Vinv);
	wcpgSafeFree(Vcpy);



	mpfr_clear(eps_1);
	mpfr_clear(frobCV);
	mpfr_clear(frobBV);
	mpfr_clear(frobV);
	mpfr_clear(frobC);
	mpfr_clear(frobB);
	mpfr_clear(frobVinv);
	mpfr_clear(sqrtn);
	mpfr_clear(eps_2);
	mpfr_clear(eps_3C);
	mpfr_clear(eps_3B);
	mpfr_clear(eps_5);
	mpfr_clear(eps_6);
	mpfr_clear(eps_4);
	mpfr_clear(scratch);
	mpfr_clear(scratch2);




	return 1;

}


/* Compute the lower bound on the WCPG matrix of a LTI filter
represented with state matrices A, B, C, D with a given 
absolute error bound eps.

W = |D| + |C*B| + |C*A*B|

We can compute only W' and we want 
||W - W'||_F <= eps

where eps>0 is given in the argument.

Since A,B,C,D  are not complex, then element-by-element absolute
value operation is performed exactly.

Then, computational errors in W' are due to two matrix multiplications and two summations.
Let us distribute error budget between two multiplications and summation in equal parts: eps/3 for each.

Step 1. We can compute (CB)' such that 
||C*B - (C*B)'||_F <= delta_CB,
where delta_CB>0 is the error matrix of the matrix multiplication.
Then, we give delta_CB=eps/3 in the argument to the multiplyAndAdd function.

Step 2. Suppose, ||A*B - (A*B)'||_F <= delta_AB.
Then, 
	||C*A*B - C*(A*B)'||_F <= ||C||_F*||A*B - (A*B)'||_F 
						   <= ||C||_F*delta_AB
Finally,
	||C*A*B - (C*(A*B)')'||_F <= ||C*A*B - C*(A*B)'||_F + ||C*(A*B)' - (C*(A*B)')'||_F
							  <= ||C||_F*delta_AB + delta_CAB,
	where delta_CAB is the error of multiplication of matrix C by (A*B)'.

Wa need ||C||_F*delta_AB + delta_CAB <= eps/3.
Then, for multiplication A*B we give
	delta_CAB <= eps/6.
	delta_AB=eps/(6*||C||_F) 
in the argument to the multiplyAndAdd function for A*B and C*(A*B) 
multiplications respectively.
NB: We use rigorous _lower_ bound on ||C||_F.

Step3. Since we need to perform two matrix summations, it is sufficient to
pass eps/6 for the sumAbs function for each summation in order to get 
overall summation error smaller than eps/3.

Summary:
	* compute |C*B|' with absolute error delta_CB = eps/3
	* compute rigorous lower bound on ||C||_F
	* compute (A*B)' with absolute error delta_AB = eps/(6*||C||_F) 
	* compute |C*(A*B)'|' with absolute error delta_CAB <= eps/6
	* compute sum S1':=|C*B|' + |C*(A*B)'|' with absolute error delta_S1 <= eps/6
	* compute sum S2':=|D| + S1' with absolute error delta_S2 <= eps/6

Since  delta_CAB=delta_S1=delta_S2=eps/6,
we can compute it just once as delta_eps6.

If eps = 2^k, then:
	* delta_CB <= 2^(k - 2)
	* delta_AB <= 2^(k - 3 - get_exp(||C||_F))
	* delta_eps6 <= 2^(k - 3)

Returns non-zero value (true) if computation is succesful and zero (false)
if could not compute lower bound on WCPG.

 */

int WCPGLowerBound(mpfr_t *W, mpfr_t* A, mpfr_t* B, mpfr_t* C, mpfr_t* D, uint64_t n, uint64_t p, uint64_t q, mpfr_exp_t k)
{


	/* Preparing Zero matrices of different sizes and scratch space*/
	if(matrixIsNan_mpfr(A, n, n))
	{
		fprintf(stderr, "Matrix A contains NaN. COuld not compute WCPGLowerBound. Exit with error\n");
		// writeMPFRMatrix(stderr, A, n, n, 10, MPFR_RNDN);
		return 0;
	}
	if(matrixIsNan_mpfr(B, n, q))
	{
		fprintf(stderr, "Matrix B contains NaN. COuld not compute WCPGLowerBound. Exit with error\n");
		// writeMPFRMatrix(stderr, B, n, q, 10, MPFR_RNDN);
		return 0;
	}
		if(matrixIsNan_mpfr(C, p, n))
	{
		fprintf(stderr, "Matrix C contains NaN. COuld not compute WCPGLowerBound. Exit with error\n");
		// writeMPFRMatrix(stderr, C, p, n, 10, MPFR_RNDN);
		return 0;
	}
	if(matrixIsNan_mpfr(D, p, q))
	{
		fprintf(stderr, "Matrix D contains NaN. COuld not compute WCPGLowerBound. Exit with error\n");
		// writeMPFRMatrix(stderr, D, p, q, 10, MPFR_RNDN);
		return 0;
	}

	mpfr_t scratch;
	mpfr_init(scratch);

	mpfr_t *scratchT;
	scratchT = allocateMPFRVector(2 * n + 1, 64);

	mpfr_t *scratch3;
	scratch3 = allocateMPFRVector(3, 64);

	mpfr_t *Zeronn;
	Zeronn = allocateMPFRMatrix(n,n, 64);
	setMatrixZero(Zeronn, n, n);

	mpfr_t *Zeropn;
	Zeropn = allocateMPFRMatrix(p,n, 64);
	setMatrixZero(Zeropn, p, n);

	mpfr_t *Zeronq;
	Zeronq = allocateMPFRMatrix(n,q, 64);
	setMatrixZero(Zeronq, n, q);

	mpfr_t *Zeropq;
	Zeropq = allocateMPFRMatrix(p,q, 64);
	setMatrixZero(Zeropq,p,q);

	/* Compute deltas for each FP matrix multiplication and summation */

	/* compute rigorous lower bound on ||C||_F */ 
	mpfr_t frobC;
	mpfr_init(frobC);
	frobeniusNormLowerBound(frobC, C, Zeropn, p, n, scratch);
	if(!mpfr_number_p(frobC))
	{
		fprintf(stderr, "||C||_F is NaN or Inf. Could not compute WCPGLowerBound. Exit with error\n");
		mpfr_clear(frobC);
		mpfr_clear(scratch);
		freeMPFRMatrix(scratchT, 1, 2 * n + 1);
		freeMPFRMatrix(scratch3, 1, 3);
		freeMPFRMatrix(Zeropn, p, n);
		freeMPFRMatrix(Zeronn, n, n);
		freeMPFRMatrix(Zeronq, n, q);
		freeMPFRMatrix(Zeropq, p, q);
		return 0;
	}
	mpfr_exp_t frbC_exp = mpfr_get_exp(frobC) - 1;
	
	/* delta_CB = eps/3 */	
	mpfr_exp_t delta_CB_exp = k - 2 - 3; 				//added 3 guard bits

	/* delta_AB = eps/(6*||C||_F)  */	
	mpfr_exp_t delta_AB_exp = k - 3 - frbC_exp - 3;		//added 3 guard bits

	/*Since  delta_CAB=delta_S1=delta_S2=eps/6,
	we can compute it just once as delta_eps6 */	
	mpfr_exp_t delta_eps6_exp = k - 3 - 3; 				//added 3 guard bits

	/* Compute |C*B|' with absolute error delta_CB */
	mpfr_t *CB, *imCB;
	CB = allocateMPFRMatrix(p,q, 64);
	imCB = allocateMPFRMatrix(p,q, 64);
	complexMatrixMultiplyAndAdd(CB, imCB, 1, C, Zeropn, B, Zeronq, 1, Zeropq, Zeropq, p, n, q, delta_CB_exp, scratchT);
	absoluteValue(CB, CB, p, q);

	if(matrixIsNan_mpfr(CB, p, q))
	{
		fprintf(stderr, "Matrix C*B contains NaN. COuld not compute WCPGLowerBound. Exit with error\n");
		mpfr_clear(frobC);
		mpfr_clear(scratch);
		freeMPFRMatrix(scratchT, 1, 2 * n + 1);
		freeMPFRMatrix(scratch3, 1, 3);
		freeMPFRMatrix(Zeropn, p, n);
		freeMPFRMatrix(Zeronn, n, n);
		freeMPFRMatrix(Zeronq, n, q);
		freeMPFRMatrix(Zeropq, p, q);
		freeMPFRMatrix(CB, p, q);
		freeMPFRMatrix(imCB, p, q);
		return 0;
	}

	/* Compute |A*B|' with absolute error delta_AB */
	mpfr_t *AB, *imAB;
	AB = allocateMPFRMatrix(n,q, 64);
	imAB = allocateMPFRMatrix(n,q, 64);
	complexMatrixMultiplyAndAdd(AB, imAB, 1, A, Zeronn, B, Zeronq, 1, Zeronq, Zeronq, n, n, q, delta_AB_exp, scratchT);

	if(matrixIsNan_mpfr(AB, n, q))
	{
		fprintf(stderr, "Matrix A*B contains NaN. COuld not compute WCPGLowerBound. Exit with error\n");
		mpfr_clear(frobC);
		mpfr_clear(scratch);
		freeMPFRMatrix(scratchT, 1, 2 * n + 1);
		freeMPFRMatrix(scratch3, 1, 3);
		freeMPFRMatrix(Zeropn, p, n);
		freeMPFRMatrix(Zeronn, n, n);
		freeMPFRMatrix(Zeronq, n, q);
		freeMPFRMatrix(Zeropq, p, q);
		freeMPFRMatrix(CB, p, q);
		freeMPFRMatrix(imCB, p, q);
		freeMPFRMatrix(AB, n, q);
		freeMPFRMatrix(imAB, n, q);
		return 0;
	}

	/* Compute |C*|A*B|'|' with absolute error delta_CAB = delta_eps6 */
	mpfr_t *CAB, *imCAB;
	CAB = allocateMPFRMatrix(p,q, 64);
	imCAB = allocateMPFRMatrix(p,q, 64);
	complexMatrixMultiplyAndAdd(CAB, imCAB, 1, C, Zeropn, AB, Zeronq, 1, Zeropq, Zeropq, p, n, q, delta_eps6_exp, scratchT);
	absoluteValue(CAB, CAB, p, q);

	if(matrixIsNan_mpfr(CAB, p, q))
	{
		fprintf(stderr, "Matrix C*A*B contains NaN. COuld not compute WCPGLowerBound. Exit with error\n");
		mpfr_clear(frobC);
		mpfr_clear(scratch);
		freeMPFRMatrix(scratchT, 1, 2 * n + 1);
		freeMPFRMatrix(scratch3, 1, 3);
		freeMPFRMatrix(Zeropn, p, n);
		freeMPFRMatrix(Zeronn, n, n);
		freeMPFRMatrix(Zeronq, n, q);
		freeMPFRMatrix(Zeropq, p, q);
		freeMPFRMatrix(CB, p, q);
		freeMPFRMatrix(imCB, p, q);
		freeMPFRMatrix(AB, n, q);
		freeMPFRMatrix(imAB, n, q);
		freeMPFRMatrix(CAB, p, q);
		freeMPFRMatrix(imCAB, p, q);
		return 0;
	}

	/* Compute sum  |C*B|' and |C*(A*B)'|' */ 

	absoluteValue(W, D, p, q);													//W = |D|
	accumulateAbsoluteValue(W, W, CB, Zeropq, p, q, delta_eps6_exp, scratch3);		//W = |D| + |C*B|'
	accumulateAbsoluteValue(W, W, CAB, Zeropq, p, q, delta_eps6_exp, scratch3);	//W = |D| + |C*B|' + |C*(A*B)'|'


	if(matrixIsNan_mpfr(W, p, q))
	{
		fprintf(stderr, "Result matrix W contains NaN. Could not compute WCPGLowerBound. Exit with error\n");
		mpfr_clear(frobC);
		mpfr_clear(scratch);
		freeMPFRMatrix(scratchT, 1, 2 * n + 1);
		freeMPFRMatrix(scratch3, 1, 3);
		freeMPFRMatrix(Zeropn, p, n);
		freeMPFRMatrix(Zeronn, n, n);
		freeMPFRMatrix(Zeronq, n, q);
		freeMPFRMatrix(Zeropq, p, q);
		freeMPFRMatrix(CB, p, q);
		freeMPFRMatrix(imCB, p, q);
		freeMPFRMatrix(AB, n, q);
		freeMPFRMatrix(imAB, n, q);
		freeMPFRMatrix(CAB, p, q);
		freeMPFRMatrix(imCAB, p, q);
		return 0;
	}


	mpfr_clear(frobC);
	mpfr_clear(scratch);
	freeMPFRMatrix(scratchT, 1, 2 * n + 1);
	freeMPFRMatrix(scratch3, 1, 3);
	freeMPFRMatrix(Zeropn, p, n);
	freeMPFRMatrix(Zeronn, n, n);
	freeMPFRMatrix(Zeronq, n, q);
	freeMPFRMatrix(Zeropq, p, q);
	freeMPFRMatrix(CB, p, q);
	freeMPFRMatrix(imCB, p, q);
	freeMPFRMatrix(AB, n, q);
	freeMPFRMatrix(imAB, n, q);
	freeMPFRMatrix(CAB, p, q);
	freeMPFRMatrix(imCAB, p, q);

	return 1;

}


int WCPGLowerBound_double(mpfr_t *W, double* A, double* B, double* C, double* D, uint64_t n, uint64_t p, uint64_t q, mpfr_exp_t k)
{
	if(matrixIsNan_double(A, n, n))
	{
		fprintf(stderr, "Matrix A has a NaN value. Could not compute WCPG. Exit with error.\n");
		return 0;
	}
	if(matrixIsNan_double(B, n, q))
	{
		fprintf(stderr, "Matrix B has a NaN value. Could not compute WCPG. Exit with error.\n");
		return 0;
	}
	if(matrixIsNan_double(C, p, n))
	{
		fprintf(stderr, "Matrix C has a NaN value. Could not compute WCPG. Exit with error.\n");
		return 0;
	}
	if(matrixIsNan_double(D, p, q))
	{
		fprintf(stderr, "Matrix D has a NaN value. Could not compute WCPG. Exit with error.\n");
		return 0;
	}

	int i;

	/* Convertion of the double matrices into multiprecision complex format */

	mp_prec_t prec = 64;
	mpfr_t *reA;
	reA = allocateMPFRMatrix(n, n, prec);
	doublerealToMPFRMatrix(reA, A, n, n);

	if(matrixIsNan_mpfr(reA, n, n))
	{
		fprintf(stderr, "MPFR Matrix A has a NaN value. Could not compute WCPG. Exit with error.\n");
		freeMPFRMatrix(reA, n, n);
		return 0;
	}

	mpfr_t *reB;
	reB = allocateMPFRMatrix(n, q, prec);
	doublerealToMPFRMatrix(reB, B, n, q);

	if(matrixIsNan_mpfr(reB, n, q))
	{
		fprintf(stderr, "MPFR Matrix B has a NaN value. Could not compute WCPG. Exit with error.\n");
		freeMPFRMatrix(reA, n, n);
		freeMPFRMatrix(reB, n, q);
		return 0;
	}

	mpfr_t *reC;
	reC = allocateMPFRMatrix(p, n, prec);
	doublerealToMPFRMatrix(reC, C, p, n);

	if(matrixIsNan_mpfr(reC, p, n))
	{
		fprintf(stderr, "MPFR Matrix C has a NaN value. Could not compute WCPG. Exit with error.\n");
		freeMPFRMatrix(reA, n, n);
		freeMPFRMatrix(reB, n, q);
		freeMPFRMatrix(reC, p, n);
		return 0;
	}

	mpfr_t *reD;
	reD = allocateMPFRMatrix(p, q, prec);
	doublerealToMPFRMatrix(reD, D, p, q);

	if(matrixIsNan_mpfr(reD, p, q))
	{
		fprintf(stderr, "MPFR Matrix D has a NaN value. Could not compute WCPG. Exit with error.\n");
		freeMPFRMatrix(reA, n, n);
		freeMPFRMatrix(reB, n, q);
		freeMPFRMatrix(reC, p, n);
		freeMPFRMatrix(reD, p, q);
		return 0;
	}

	// fprintf(stderr, "Filter order: %llu.\nNumber of inputs: %llu \nNumber of outputs: %llu \n",n, q, p);
	// printf("matrix A: \n");
	// writeMPFRMatrix(stderr, reA, n, n, 10, MPFR_RNDN);
	// printf("matrix B: \n");
	// writeMPFRMatrix(stderr, reB, n, q, 10, MPFR_RNDN);
	// printf("matrix C: \n");
	// writeMPFRMatrix(stderr, reC, p, n, 10, MPFR_RNDN);
	// printf("matrix D: \n");
	// writeMPFRMatrix(stderr, reD, p, q, 10, MPFR_RNDN);


	if(!WCPGLowerBound(W, reA, reB, reC, reD, n, p, q, k))
	{
		// printf("Could not compute lower bound on the Worst Case Peak Gain matrix.\n");
		freeMPFRMatrix(reA, n, n);
		freeMPFRMatrix(reB, n, q);
		freeMPFRMatrix(reC, p, n);
		freeMPFRMatrix(reD, p, q);
		return 0;
	}
	else
	{
		// printf("Lower bound on the WCPG matrix: \n");
		// writeMPFRMatrix(stderr, W, p, q, 10, MPFR_RNDN);
		freeMPFRMatrix(reA, n, n);
		freeMPFRMatrix(reB, n, q);
		freeMPFRMatrix(reC, p, n);
		freeMPFRMatrix(reD, p, q);
		return 1;
	}	

}



/**
 * @brief For an LTI filter given in its State-Space representation {A,B,C,D},
where A is n*n, B is n*q, C is p*n and D is p*q real matrix the function 
returns integer value indicating if WCPG is successfully computed.
In p*q matrix W the Worst-Case peak gain is stored if algorithm successfully exited. */
int WCPG_ABCD(double *W, double *A, double *B, double *C, double *D, uint64_t n, uint64_t p, uint64_t q)
{
	int flag = 0;



	mpfr_t *S;
	S = allocateMPFRMatrix(p, q, 64);


	mpfr_t mpeps;
	mpfr_init2(mpeps, 64);
	mpfr_set_str(mpeps, "1e-53", 2, MPFR_RNDN);


	if (!WCPG_ABCD_mprec(S, A, B, C, D, n, p, q, mpeps))
	{
		freeMPFRMatrix(S, p, q);
		mpfr_clear(mpeps);
		return 0;
	}
	else
	{

		int i, j;
		for(i = 0; i < p; ++i)
		{
			for(j = 0; j < q; ++j)
			{
				W[i * q + j] = mpfr_get_d(S[i * q + j], MPFR_RNDN);
			}
		}

		freeMPFRMatrix(S, p, q);
		mpfr_clear(mpeps);
		return 1;
	}

}

/* 
	This function computes the WCPG of a SISO filter represented with its transfer function.
	The filter is specified as following:

			sum_0^Nb b_i
		H=	-------------
			sum_0^Na a_i

	We consider only the case, when Nb <= Na.
	The order of the filter is equal then to fOrder = Na - 1
	
	In case of success, the function returns 1 and W contains the WCPG computed with double precision.
	Otherwise, the function returns 0 and W contains a NaN.
*/
int WCPG_tf(double *W, double *num, double *denum, uint64_t Nb, uint64_t Na)
{
	int i;

	if(Nb > Na)
	{
		fprintf(stderr, "\nCannot compute the WCPG_tf: numerator order is greater than denumerator order.\nReminer: we do not consider coeffa[0]=1 implicitly.\nTansfer function format is\n\tsum_0^Nb coeffb_i\n  H=\t-----------------\n\tsum_0^Na coeffa_i\n \n");
		//*W = -1;
		return 0;
	}

	uint64_t fOrder = Na - 1;

	double *a = wcpgSafeCalloc(Na, Na * sizeof(double));
	if(denum[0] != (double)1)
	{
		a[0] = 1;
		for(i=1; i < Na; ++i)
		{
			a[i] = denum[i]/denum[0];
		}
	}
	else
	{
		for(i=0; i < Na; ++i)
		{
			a[i] = denum[i];
		}
	}

	double *b = wcpgSafeCalloc(Na, Na * sizeof(double));
	for(i=0; i < Nb; ++i)
	{
			b[i] = num[i];
	}

	uint64_t p = 1;			//number of outputs
	uint64_t q = 1;			//number of inputs
	uint64_t n = fOrder;	//number of states

	double *A = (double*)wcpgSafeCalloc(n * n, n * n * sizeof(double));
	double *B = (double*)wcpgSafeCalloc(n * q, n * q * sizeof(double));
	double *C = (double*)wcpgSafeCalloc(p * n, p * n * sizeof(double));
	double *D = (double*)wcpgSafeCalloc(p * q, p * q * sizeof(double));

	// Computing matrix A

	// first line is [-a[1]  ...  -a[Na]]
	for(i = 0; i < n; ++i)
	{
		A[i] = -a[i+1];  
	}
	for(i = 1; i < n; ++i)
	{
		A[i * n + i - 1] = 1.0;	
	}

	// Matrix B
	B[0] = 1;

	//Matrix C
	for(i = 0; i < n; ++i)
	{
		C[i] = b[i+1] - b[0]*a[i+1];
	}

	//Matrix D
	D[0] = b[0];

/*
	printf("Matrix A: \n");
	clapack_matrix_print_d(A, n, n);
	printf("Matrix B: \n");
	clapack_matrix_print_d(B, n, 1);
	printf("Matrix C: \n");
	clapack_matrix_print_d(C, 1, n);
	printf("Matrix D: \n");
	clapack_matrix_print_d(D, 1, 1);
	printf("------> vector a: \n");
	clapack_matrix_print_d(a, 1, Na);
	printf("------> vector b: \n");
	clapack_matrix_print_d(b, 1, Na);
	*/ 
	
	if(!WCPG_ABCD(W, A, B, C, D, n, 1, 1))
	{
		wcpgSafeFree(a);
		wcpgSafeFree(b);
		wcpgSafeFree(A);
		wcpgSafeFree(B);
		wcpgSafeFree(C);
		wcpgSafeFree(D);
		return 0;
	}
	else
	{
		wcpgSafeFree(a);
		wcpgSafeFree(b);
		wcpgSafeFree(A);
		wcpgSafeFree(B);
		wcpgSafeFree(C);
		wcpgSafeFree(D);
		return 1;
	}

}



/* Nth order LTI filter is represented by its transfer function numerator (array of size Nb) and
denumerator (array of size Na), where N := max(Na, Nb).
For such a filter, the function computes its WCPG */
int WCPG_tf_initial(double *W, double *num, double *denum, uint64_t Nb, uint64_t Na)
{

	// FIrst, we need to convert tf representation to state-space
	//and then call the general WCPG function.


	uint64_t flag = 0;

	uint64_t n = (Nb-1) >= Na ? Nb - 1 : Na;		//n = max(Nb, Na)

	double *alpha = (double*)wcpgSafeCalloc(n, n * sizeof(double));
	double *beta = (double*)wcpgSafeCalloc(n, n * sizeof(double));
	double *b = (double*)wcpgSafeCalloc(n+1, (n+1) * sizeof(double));
	int i, j;
	for(i = 0; i < n; ++i)
	{
		if(i < (Na)) alpha[i] = denum[i];
		else alpha[i] = 0.0;	
	}

	for(i = 0; i < (n+1); ++i)
	{
		if(i < Nb + 1) b[i] = num[i];
		else b[i] = 0.0;
	}

	double b0 = num[0];

	for(i = 0; i < n; ++i)
	{
		beta[i] = b[i+1] - b0 * alpha[i];
	}
	uint64_t p = 1;
	uint64_t q = 1;


	double *A = (double*)wcpgSafeCalloc(n * n, n * n * sizeof(double));
	double *B = (double*)wcpgSafeCalloc(n * q, n * q * sizeof(double));
	double *C = (double*)wcpgSafeCalloc(p * n, p * n * sizeof(double));
	double *D = (double*)wcpgSafeCalloc(p * q, p * q * sizeof(double));

	//building matrix A;
	for(i = 0; i < n; ++i)
	{
		if(i == 0)
		{
			for(j = 0; j < n; ++j)
			{
				A[i * n + j] = -alpha[j];
			}
		}
		else
			A[i * n + i - 1] = 1.0;	

	}

	//building matrix B
	B[0] = 1.0;
	for(i = 1; i < n; ++i)
		B[i] = 0.0;

	//building matrix C
	for(i = 0; i < n; ++i)
		C[i] = beta[i];

	//building matrix D
	D[0] = b0;


	//printing result state-space
	/*
	 printf("Matrix A: \n");
	 clapack_matrix_print_d(A, n, n);
	 printf("Matrix B: \n");
	 clapack_matrix_print_d(B, n, 1);
	 printf("Matrix C: \n");
	 clapack_matrix_print_d(C, 1, n);
	 printf("Matrix D: \n");
	 clapack_matrix_print_d(D, 1, 1);
	 */

	//computing WCPG with double precision
	if(!WCPG_ABCD(W, A, B, C, D, n, 1, 1))
	{
		wcpgSafeFree(alpha);
		wcpgSafeFree(beta);
		wcpgSafeFree(b);
		wcpgSafeFree(A);
		wcpgSafeFree(B);
		wcpgSafeFree(C);
		wcpgSafeFree(D);
		return 0;
	}
	else
	{
		wcpgSafeFree(alpha);
		wcpgSafeFree(beta);
		wcpgSafeFree(b);
		wcpgSafeFree(A);
		wcpgSafeFree(B);
		wcpgSafeFree(C);
		wcpgSafeFree(D);
		return 1;
	}
	
}

int MPFItrustedEigendecompsition(mpfi_t *reLambdaint, mpfi_t *imLambdaint, mpfi_t *reVint, mpfi_t *imVint, mpfr_t *A, uint64_t n)
{
	mpfr_prec_t prec = 500;
	mp_exp_t prec_exp = -500;

	mpfr_t *re = allocateMPFRMatrix(1, n, prec);
	mpfr_t *im = allocateMPFRMatrix(1, n, prec);
	mpfr_t *reV = allocateMPFRMatrix(n, n, prec);
	mpfr_t *imV = allocateMPFRMatrix(n, n, prec);

	mpfr_t epsLambda, epsV;
	mpfr_init2(epsLambda, 64);
	mpfr_init2(epsV, 64);



	if(!eigen(re, im, reV, imV, epsLambda, epsV, A, n, prec))
	{
		fprintf(stderr, "ERROR: could not compute trusted intervals for eigendecomposition \n");
		return 0;
	}
	mpfr_t *scratch3 = allocateMPFRMatrix(1, 3, prec);
	mpfr_t scratch;
	mpfr_init(scratch);


	/* Here we need to check the inclusion, but this part is under reconstruction.*/


	/* Construct the interval matrices for eigenvalues and eigenvalues */
	MPFIComplexMatrixMidrad_fr(reVint, imVint, reV, imV, n, n, epsV, scratch3);
	MPFIComplexMatrixMidrad_fr(reLambdaint, imLambdaint, re, im, 1, n, epsLambda, scratch3);

/*	printf("Eigenvalues:\n");
	MPFIComplexMatrixPrint(reLambdaint, imLambdaint, n, 1);
*/
	/* CHeck if intervals abs(Lambda) contain 1 */
	int i;
	for(i = 0; i < n; ++i)
	{

		mpfr_set_prec(scratch3[0], prec);
		mpfr_set_prec(scratch3[1], prec);
		mpfr_set_prec(scratch3[2], prec);
		if(mpfi_cmp_ui(reLambdaint[i], (unsigned long int)0) >= 0)
		{
			mpfi_get_right(scratch3[0], reLambdaint[i]);
			mpfr_add(scratch3[0], scratch3[0], epsLambda, MPFR_RNDU);
		}
		else
		{
			mpfi_get_left(scratch3[0], reLambdaint[i]);
			mpfr_sub(scratch3[0], scratch3[0], epsLambda, MPFR_RNDU);
		}

		if(mpfi_cmp_ui(imLambdaint[i], (unsigned long int)0) >= 0)
		{
			mpfi_get_right(scratch3[1], imLambdaint[i]);
			mpfr_add(scratch3[1], scratch3[1], epsLambda, MPFR_RNDU);
		}
		else
		{
			mpfi_get_left(scratch3[1], imLambdaint[i]);
			mpfr_sub(scratch3[1], scratch3[1], epsLambda, MPFR_RNDU);
		}

		computeAbsoluteValue(scratch3[2], scratch3[0], scratch3[1], prec_exp, scratch);

		if(mpfr_cmp_ui(scratch3[2], (unsigned long int) 1) >= 0)	//if abs(lambda) > 1
		{
			fprintf(stderr, "ERROR: exists an eigenvalue with its absolute value larger than 1: abs(lambda) = \n");
			mpfr_out_str(stderr, 10, 10,scratch3[2], MPFR_RNDN);
			fprintf(stderr, "\n");
			return 0;
		}
	}

	mpfr_clear(scratch);
	freeMPFRMatrix(scratch3, 1, 3);
	freeMPFRMatrix(re, 1, n);
	freeMPFRMatrix(im, 1, n);
	freeMPFRMatrix(reV, n, n);
	freeMPFRMatrix(imV, n, n);

	return 1;
}


int eigen(mpfr_t *re, mpfr_t *im, mpfr_t *reV, mpfr_t *imV, mpfr_t epsLambda, mpfr_t epsV, mpfr_t *A, uint64_t n, mpfr_prec_t prec)
{
	mpfr_exp_t prec_result = -prec / 4;
	int flag;
    //flag = mpfr_eigen(n, A, re, im, (int)1, reV, imV, prec);

	flag = 0;

	if(flag == 2)
	{
		mpfr_set_ui_2exp(epsLambda, (unsigned long int) 1, prec_result, MPFR_RNDN);
		mpfr_set_ui_2exp(epsV, (unsigned long int) 1, prec_result, MPFR_RNDN);
		if(re == NULL || im == NULL || reV == NULL || imV == NULL)
		{
			fprintf(stderr, "ERROR in eigendecomposition computation: null pointers returned\n");
			return 0;
		}
		return 1;
	}
	else /* If flag is 0 or 1 it means that we need to compute at least eigenvectors with LAPACK */
	{

		complexdouble *V = (complexdouble*)wcpgSafeMalloc(n*n*sizeof(complexdouble));
		complexdouble *eigv = (complexdouble*)wcpgSafeMalloc(n*sizeof(complexdouble));

		doublereal *Acpy = (doublereal*)wcpgSafeMalloc(n*n*sizeof(doublereal));
		MPFRMatrixToDouble(Acpy, A, n, n);

		doublereal *eerrbnd = (doublereal*)wcpgSafeMalloc(n*sizeof(doublereal));
		doublereal *verrbnd = (doublereal*)wcpgSafeMalloc(n*sizeof(doublereal));
		doublereal eps_prec = DBL_EPSILON;



		if(!clapack_eigenSolver(eigv, V, verrbnd, eerrbnd, Acpy, n, eps_prec))	//eigensolver, which returns approx error bounds
		{
			fprintf(stderr, "ERROR: Could not compute eigensystem using LAPACK \n");
			wcpgSafeFree(V);
			wcpgSafeFree(eigv);
			wcpgSafeFree(Acpy);
			wcpgSafeFree(eerrbnd);
			wcpgSafeFree(verrbnd);
			return 0;
		}

		complexdoubleToMPFRMatrix(reV, imV, V, n, n);
		maxInVector(epsV, verrbnd, n);

		if(flag == 0)	/* If eigenvalues were not computed by the MPFR solver, we take the LAPACK-computed eigenvalues. */
		{
			complexdoubleToMPFRMatrix(re, im, eigv, 1, n);
			maxInVector(epsLambda, eerrbnd, n);
			//fprintf(stderr, "WARNING: could not compute the eigenvalues and eigenvectors with MPFR solver. Used LAPACK instead. \n");
		}
		else
		{
			//fprintf(stderr, "WARNING: could not compute the eigenvectors with MPFR solver. Used LAPACK instead. \n");
			mpfr_set_ui_2exp(epsLambda, (unsigned long int) 1, prec_result, MPFR_RNDN);
		}



		wcpgSafeFree(V);
		wcpgSafeFree(eigv);
		wcpgSafeFree(Acpy);
		wcpgSafeFree(eerrbnd);
		wcpgSafeFree(verrbnd);
		return 1;

	}

}


void WCPG_welcome_message(int n)
{
	printf("Hello, I'm the WCPG library! \n n = %d", n);
}
