#ifndef MPFR_EIGENDECOMPOSITION
#define MPFR_EIGENDECOMPOSITION

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <mpfr.h>
#include <stdint.h>

#include "aux_funcs.h"

mpfr_t* allocateMatrix(int n);
mpfr_t* allocateVector(int n);
void initMatrix(mpfr_t* A, int n, mpfr_prec_t prec);
void initVector(mpfr_t* v, int n, mpfr_prec_t prec);
void clearMatrix(mpfr_t* A, int n);
void clearVector(mpfr_t* v, int n);
mpfr_t* importMatrix(int* n, char* filename, mpfr_prec_t prec);
void exportMatrix(mpfr_t* A, int n, char* filename);
void writeMPFRMatrix(FILE *stream, mpfr_t *A, uint64_t m, uint64_t n,size_t nmbr, mpfr_rnd_t rnd);
void writeMPFRComplexMatrix(FILE *stream, mpfr_t *ReA, mpfr_t *ImA, uint64_t m, uint64_t n,size_t nmbr, mpfr_rnd_t rnd);
void writeMatrix(mpfr_t* A, int n, FILE *f);

void writeVectorComplex(mpfr_t* re, mpfr_t* im, int n, FILE *f);

void set0Matrix(mpfr_t* A, int n);
void set0Vector(mpfr_t* v, int n);
void setAleaMatrix(mpfr_t* A, int n, mpfr_t tmp);
void setAleaVector(mpfr_t* v, int n, mpfr_t tmp);
void identityMatrix(mpfr_t* A, int n);
void copyVector(mpfr_t* dest, mpfr_t* src, int n);
void eigcode_copyMatrix(mpfr_t* dest, mpfr_t* src, int n);
void trans(mpfr_t* A, int n, mpfr_t tmp);

void norm2Complex(mpfr_t renv, int n, mpfr_t* xre, mpfr_t* xim, mpfr_t tmp);
void normInfVectComplex(mpfr_t renv, mpfr_t* vre, mpfr_t* vim, int n, mpfr_t tmp);
void normInfMatrix(mpfr_t renv, mpfr_t* A, int n, mpfr_t tmp1, mpfr_t tmp2);
void norm2(mpfr_t renv, int n, mpfr_t* x, int incx);

void toHessenbergForm(int n, mpfr_t* A, mpfr_t* U, mpfr_prec_t prec);
void GivensCoefficients(int n, mpfr_t* A, int i, int j, int k, mpfr_t c, mpfr_t s, mpfr_t tmp);
void GivensRotationGH(int n, mpfr_t* A, int i, int j, mpfr_t c, mpfr_t s, mpfr_t tmp1, mpfr_t tmp2);
void GivensRotationDH(int n, mpfr_t* A, int i, int j, mpfr_t c, mpfr_t s, mpfr_t tmp1, mpfr_t tmp2);
void subDiagonalUpdate(int n, mpfr_t* A, mpfr_t tol, mpfr_t tmp1, mpfr_t tmp2);
void deflation(int n, mpfr_t* A, int* p, int* q);
void QRstepH(int n, mpfr_t* A, mpfr_t* rot, int p, int q, mpfr_t c, mpfr_t s, mpfr_t tmp1, mpfr_t tmp2, mpfr_t mu);
int iterQR(int n, mpfr_t* A, mpfr_prec_t prec);
int findEigenValues(int n, mpfr_t* A, mpfr_t* re, mpfr_t* im, mpfr_prec_t prec);

// Fonction de traitement de complexes et de calcul des vecteurs propres
void addComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2);
void subComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2);
void mulComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2, mpfr_t tmp);
void divComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp);
void luDecompositionH(mpfr_t* Re, mpfr_t* Im, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp1, mpfr_t tmp2, mpfr_t tmp3);
void upperTriangularSolve(mpfr_t* Re, mpfr_t* Im, mpfr_t* bre, mpfr_t* bim, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp1, mpfr_t tmp2, mpfr_t tmp3);
void lowerTriangularSolve(mpfr_t* Re, mpfr_t* Im, mpfr_t* bre, mpfr_t* bim, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp);
void systemSolveLU(mpfr_t* Re, mpfr_t* Im, mpfr_t* bre, mpfr_t* bim, int n, mpfr_t tmp1, mpfr_t tmp2, mpfr_t tmp3, mpfr_t tmp4, mpfr_t tmp5);
void mulMatrixComplexReal(mpfr_t* Cre, mpfr_t* Cim, mpfr_t* Are, mpfr_t* Aim, mpfr_t* Bre, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp1, mpfr_t tmp2);
int findEigenVectors(int n, mpfr_t* A, mpfr_t* eigVecre, mpfr_t* eigVecim, mpfr_t* re, mpfr_t* im, mpfr_prec_t prec);
////

int mpfr_eigen(int n, mpfr_t* A, mpfr_t* re, mpfr_t* im, int ask_for_Evector, mpfr_t* evre, mpfr_t* evim, mpfr_prec_t prec);
void exportVectorComplex(mpfr_t* re, mpfr_t* im, int n, char* filename);
void exportMatrixComplex(mpfr_t* re, mpfr_t* im, int n, char* filename);





#ifdef __cplusplus
}
#endif

#endif






















