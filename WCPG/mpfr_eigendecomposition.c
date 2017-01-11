#include "mpfr_eigendecomposition.h"

// Fonction d'allocation memoire d'une matrice carre nxn
mpfr_t* allocateMatrix(int n) {
	return wcpgSafeMalloc(n * n * sizeof(mpfr_t));
}

// Fonction d'allocation memoire d'un vecteur de taille n
mpfr_t* allocateVector(int n) {
	return wcpgSafeMalloc(n * sizeof(mpfr_t));
}

//// La bibliothèque mpfr nécessite l'initialisation des variables mpfr_t
//// qui contiennent les flottants multiprécisions avant de les utiliser.
// Initialisation des mpfr_t d'une matrice carre A de taille nxn
void initMatrix(mpfr_t* A, int n, mpfr_prec_t prec) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_init2(A[i*n + j], prec);
}

// Initialisation des mpfr_t d'un vecteur v de taille n
void initVector(mpfr_t* v, int n, mpfr_prec_t prec) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_init2(v[i], prec);
}
////

// Fonction de libération de mémoire d'une matrice carré A de taille nxn
void clearMatrix(mpfr_t* A, int n) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_clear(A[i*n + j]);
	wcpgSafeFree(A);
}



// Fonction de libération de mémoire d'un vecteur v de taille n
void clearVector(mpfr_t* v, int n) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_clear(v[i]);
	wcpgSafeFree(v);
}

// Importation d'une matrice carré à partir d'un fichier.
// Le format est le suivant :
// Le premier chiffre rencontré est la dimension n de la matrice (on aura A de taille nxn)
// Puis les coefficients de la matrice séparés par des blancs ligne par ligne
mpfr_t* importMatrix(int* n, char* filename, mpfr_prec_t prec) {
	FILE* f = fopen(filename, "r");
	int i, j;
	mpfr_t* renv;
	
	fscanf(f, " %d ", n);
	renv = allocateMatrix(*n);
	initMatrix(renv, *n, prec);
	
	for(i = 0 ; i < *n ; i++)
		for(j = 0 ; j < *n ;j++)
			mpfr_inp_str(renv[i*(*n) + j], f, 16, MPFR_RNDN);
	
	fclose(f);
	return renv;
}

// Exportation d'une matrice dans un fichier
// Ecrit la matrice A dans un fichier
// Le format de importMatrix est respecté, si bien qu'elle peut être relu par celle-ci
void exportMatrix(mpfr_t* A, int n, char* filename) {
	FILE* f = fopen(filename, "w");
	int i, j;
	
	fprintf(f, "%d\n", n);
	
	for(i = 0 ; i < n ; i++, fprintf(f, "\n"))
		for(j = 0 ; j < n ; j++, fprintf(f, "\t"))
			mpfr_out_str(f, 10, 0, A[i*n + j], MPFR_RNDN);
	
	fclose(f);
}

void writeMatrix(mpfr_t* A, int n, FILE *f) {
	int i, j;
	fprintf(f, "\n Matrix %d\n", n);

	for(i = 0 ; i < n ; i++, fprintf(f, "\n"))
		for(j = 0 ; j < n ; j++, fprintf(f, "\t"))
			mpfr_out_str(f, 10, 0, A[i*n + j], MPFR_RNDN);

	fclose(f);
}

// Donne à tout les éléments de la matrice A la valeur 0
void set0Matrix(mpfr_t* A, int n) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_set_d(A[i*n + j], 0, MPFR_RNDN);
}

// Donne à tout les éléments du vecteur v la valeur 0
void set0Vector(mpfr_t* v, int n) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_set_d(v[i], 0, MPFR_RNDN);
}

// Initialise la matrice A aleatoirement, chaque ligne etant unitaire
void setAleaMatrix(mpfr_t* A, int n, mpfr_t tmp) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_set_d(A[i*n + j], ((double) rand()) / RAND_MAX, MPFR_RNDN);
	for(i = 0 ; i < n ; i++) {
		norm2(tmp, n, A + i, 1);
		for(j = 0 ; j < n ; j++)
			mpfr_div(A[i*n + j], A[i*n + j], tmp, MPFR_RNDN); // A[i*n + j] /= tmp
	}
}

// Initialise le vecteur v aléatoirement. Le vecteur produit est unitaire
void setAleaVector(mpfr_t* v, int n, mpfr_t tmp) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_set_d(v[i], ((double) rand()) / RAND_MAX, MPFR_RNDN);
	norm2(tmp, n, v, 1);
	for(i = 0 ; i < n ; i++)
		mpfr_div(v[i], v[i], tmp, MPFR_RNDN); // v[i] /= tmp;
}

// A := In (matrice unité de taille n)
void identityMatrix(mpfr_t* A, int n) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_set_d(A[i*n + j], 0, MPFR_RNDN);
	for(i = 0 ; i < n ; i++)
		mpfr_set_d(A[i*n + i], 1, MPFR_RNDN);
}

void copyVector(mpfr_t* dest, mpfr_t* src, int n) {
	int i;
	for(i = 0 ; i < n ; i++)
		mpfr_set(dest[i], src[i], MPFR_RNDN);
}

void eigcode_copyMatrix(mpfr_t* dest, mpfr_t* src, int n) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++)
			mpfr_set(dest[i*n + j], src[i*n + j], MPFR_RNDN);
}

// A := A^^t
void trans(mpfr_t* A, int n, mpfr_t tmp) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < i ; j++) {
			mpfr_set(tmp, A[i*n + j], MPFR_RNDN);
			mpfr_set(A[i*n + j], A[j*n + i], MPFR_RNDN);
			mpfr_set(A[j*n + i], tmp, MPFR_RNDN);
		}
}

// Calcul de la norme 2 d'un vecteur x = xre + i * xim (vecteur complexe)
// Resultat placé dans renv
void norm2Complex(mpfr_t renv, int n, mpfr_t* xre, mpfr_t* xim, mpfr_t tmp) {
	int i;
	mpfr_set_d(renv, 0, MPFR_RNDN);
	for(i = 0 ; i < n ; i++) {
		mpfr_mul(tmp, xre[i], xre[i], MPFR_RNDN);
		mpfr_add(renv, renv, tmp, MPFR_RNDN);
		mpfr_mul(tmp, xim[i], xim[i], MPFR_RNDN);
		mpfr_add(renv, renv, tmp, MPFR_RNDN);
		// renv += xre[i] * xre[i] + xim[i] * xim[i];
	}
	mpfr_sqrt(renv, renv, MPFR_RNDN);
}

// Calcul de la norme 2 d'un vecteur v = vre + i * vim (vecteur complexe)
// Resultat placé dans renv
void normInfVectComplex(mpfr_t renv, mpfr_t* vre, mpfr_t* vim, int n, mpfr_t tmp) {
	int i;
	mpfr_set_d(renv, -1, MPFR_RNDN);
	for(i = 0 ; i < n ; i++) {
		mpfr_hypot(tmp, vre[i], vim[i], MPFR_RNDN);
		if(mpfr_cmp(renv, tmp) < 0)
			mpfr_set(renv, tmp, MPFR_RNDN);
	}
}

// Calcul de la norme 2 d'une matrice A (matrice réelle)
// Resultat placé dans renv
void normInfMatrix(mpfr_t renv, mpfr_t* A, int n, mpfr_t tmp1, mpfr_t tmp2) {
	int i, j;
	mpfr_set_d(renv, -1, MPFR_RNDN);
	for(i = 0 ; i < n ; i++) {
		mpfr_set_d(tmp1, 0, MPFR_RNDN);
		for(j = 0 ; j < n ; j++) {
			mpfr_abs(tmp2, A[i*n + j], MPFR_RNDN);
			mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
			// tmp += fabs(A[i*n + j]);
		}
		if(mpfr_cmp(renv, tmp1) < 0)
			mpfr_set(renv, tmp1, MPFR_RNDN);
	}
}

// Calcul de la norme 2 d'un vecteur x de taille n
// Le résultat est placé dans la variable renv
// Le vecteur x est vu de la sorte :
// x_0 = x[incx * 0]
// x_1 = x[incx * 1]
// ...
// x_k = x[incx * k]
void norm2(mpfr_t renv, int n, mpfr_t* x, int incx) {
	mpfr_set_d(renv, 0, MPFR_RNDN);
	int i, tmp_indx;
	for(i = 0, tmp_indx = 0 ; i < n ; i++, tmp_indx += incx)
		mpfr_hypot(renv, x[tmp_indx], renv, MPFR_RNDN);
}

// Réduction de la matrice A à une matrice d'Hessenberg par multiplications
// successives par des matrices d'householder
// Le résultat écrase le contenu de la matrice A
// Calcul la matrice de passage orthogonale U
// A := U^t * A * U
void toHessenbergForm(int n, mpfr_t* A, mpfr_t* U, mpfr_prec_t prec) {
	
	mpfr_t alpha, r, tmp_calc;
	int i, j, k, taille;
	mpfr_t* x = allocateVector(n); // Vecteur de la matrice d'Householder
	mpfr_t* tmp_A = allocateVector(n); // Vecteur de Travail
	mpfr_t* tmp_U = allocateVector(n); // Vecteur de Travail
	
	mpfr_init2(alpha, prec);
	mpfr_init2(r, prec);
	mpfr_init2(tmp_calc, prec);
	
	initVector(x, n, prec);
	initVector(tmp_A, n, prec);
	initVector(tmp_U, n, prec);
	
	// Réduction des différentes colonnes pour arriver à une forme d'Hessenberg
	for(i = 0 ; i < n - 2 ; i++) { 
		taille = n - i;
		
		// Le vecteur x va contenir le vecteur de la matrice d'Householder voulu
		for(j = i ; j < n ; j++) // Copie de la partie de colonne a supprimer
			mpfr_set(x[j - i], A[j*n + i], MPFR_RNDN);

		// Mise à zéro du premier élément du vecteur
		// car on ne touche pas au premier élément de la colonne à réduire
		mpfr_set_d(*x, 0, MPFR_RNDN);
		
		// Nouvel element sous-diagonal
		norm2(alpha, taille, x, 1);
		mpfr_copysign(alpha, alpha, x[1], MPFR_RNDN);
		mpfr_mul_d(alpha, alpha, -1, MPFR_RNDN);
		//// alpha = -copysign(norm2(taille, x, 1), x[1])
		
		mpfr_mul(r, alpha, alpha, MPFR_RNDN);
		mpfr_mul(tmp_calc, alpha, x[1], MPFR_RNDN);
		mpfr_sub(r, r, tmp_calc, MPFR_RNDN);
		mpfr_div_d(r, r, 2, MPFR_RNDN);
		mpfr_sqrt(r, r, MPFR_RNDN);
		//// r = sqrt((alpha * alpha - x[1] * alpha) / 2)
		
		mpfr_mul_d(tmp_calc, r, 2, MPFR_RNDN); // tmp_calc contient alors 2*r
		
		mpfr_sub(x[1], x[1], alpha, MPFR_RNDN);
		mpfr_div(x[1], x[1], tmp_calc, MPFR_RNDN);
		
		for(j = 2 ; j < taille ; j++) // Calcul du vecteur de la matrice d'Householder
			mpfr_div(x[j], x[j], tmp_calc, MPFR_RNDN);
		
		//// Calcul de HA
		for(j = 0 ; j < taille ; j++) { // Calcul de tmp_A := x^^t * A
			mpfr_set_d(tmp_A[j], 0, MPFR_RNDN);
			for(k = 0 ; k < taille ; k++) {
				mpfr_mul(tmp_calc, x[k], A[(i + k)*n + i + j], MPFR_RNDN);
				mpfr_add(tmp_A[j], tmp_A[j], tmp_calc, MPFR_RNDN);
			}
		}
		for(j = 0 ; j < taille ; j++) // Calcul de A - 2*x*tmp_A = A - 2*x*(x^^t * A)
			for(k = 0 ; k < taille ; k++) {
				mpfr_mul(tmp_calc, tmp_A[k], x[j], MPFR_RNDN);
				mpfr_mul_d(tmp_calc, tmp_calc, 2, MPFR_RNDN);
				mpfr_sub(A[(i + j)*n + i + k], A[(i + j)*n + i + k], tmp_calc, MPFR_RNDN);
			}
		////
		//// Calcul de (HA)H
		for(j = 0 ; j < n ; j++) { // Calcul de tmp_A := (HA) * x
			mpfr_set_d(tmp_A[j], 0, MPFR_RNDN);
			mpfr_set_d(tmp_U[j], 0, MPFR_RNDN);
			for(k = 0 ; k < taille ; k++) {
				mpfr_mul(tmp_calc, x[k], A[j*n + i + k], MPFR_RNDN);
				mpfr_add(tmp_A[j], tmp_A[j], tmp_calc, MPFR_RNDN);
				// tmp_A[j] += x[k] * A[j*n + i + k];
				mpfr_mul(tmp_calc, x[k], U[j*n + i + k], MPFR_RNDN);
				mpfr_add(tmp_U[j], tmp_U[j], tmp_calc, MPFR_RNDN);
				// tmp_U[j] += x[k] * U[j*n + i + k];
			}
		}
		for(j = 0 ; j < n ; j++) // Calcul de HAH = HA - 2 * tmp_A * x^^t = HA - 2 * ((HA) * x) * x^^t
			for(k = 0 ; k < taille ; k++) {
				mpfr_mul(tmp_calc, tmp_A[j], x[k], MPFR_RNDN);
				mpfr_mul_d(tmp_calc, tmp_calc, 2, MPFR_RNDN);
				mpfr_sub(A[j*n + i + k], A[j*n + i + k], tmp_calc, MPFR_RNDN);
				// A[j*n + i + k] -= 2 * x[k] * tmp_A[j];
				
				mpfr_mul(tmp_calc, tmp_U[j], x[k], MPFR_RNDN);
				mpfr_mul_d(tmp_calc, tmp_calc, 2, MPFR_RNDN);
				mpfr_sub(U[j*n + i + k], U[j*n + i + k], tmp_calc, MPFR_RNDN);
				// U[j*n + i + k] -= 2 * x[k] * tmp_U[j];
			}
		////
	}
	
	// Mise à zéro explicite des éléments sous la sous-diagonale
	for(i = 2 ; i < n ; i++)
		for(j = 0 ; j < i - 1 ; j++)
			mpfr_set_d(A[i*n + j], 0, MPFR_RNDN);
	
	// Libération des variables utilisées
	mpfr_clear(alpha);
	mpfr_clear(r);
	mpfr_clear(tmp_calc);
	
	clearVector(x, n);
	clearVector(tmp_A, n);
	clearVector(tmp_U, n);
}

// Donne les coefficients de la matrice de givens G qui effectue la transformation suivante :
// G * (A_ik) en (x)
//     (A_jk)    (0)
void GivensCoefficients(int n, mpfr_t* A, int i, int j, int k, mpfr_t c, mpfr_t s, mpfr_t tmp) {
	mpfr_hypot(tmp, A[i*n + k], A[j*n + k], MPFR_RNDN);
	mpfr_div(c, A[i*n + k], tmp, MPFR_RNDN);
	mpfr_div(s, A[j*n + k], tmp, MPFR_RNDN);
	mpfr_mul_d(s, s, -1, MPFR_RNDN);
}

// Applique une multiplication à gauche sur une matrice d'Hessenberg
// par une matrice de rotation de givens de coefficients c et s
// ! A doit être une matrice d'Hessenberg
void GivensRotationGH(int n, mpfr_t* A, int i, int j, mpfr_t c, mpfr_t s, mpfr_t tmp1, mpfr_t tmp2) {
	int ind;
	for(ind = i ; ind < n ; ind++) {
		mpfr_set(tmp1, A[i*n + ind], MPFR_RNDN);
		mpfr_mul(tmp2, A[j*n + ind], s, MPFR_RNDN);
		mpfr_mul(A[i*n + ind], A[i*n + ind], c, MPFR_RNDN);
		mpfr_sub(A[i*n + ind], A[i*n + ind], tmp2, MPFR_RNDN);
		//// A[i*n + ind] = A[i*n + ind] * c - A[j*n + ind] * s
		
		mpfr_mul(tmp2, tmp1, s, MPFR_RNDN);
		mpfr_mul(A[j*n + ind], A[j*n + ind], c, MPFR_RNDN);
		mpfr_add(A[j*n + ind], A[j*n + ind], tmp2, MPFR_RNDN);
		//// A[j*n + ind] = tmp1 * s + A[j*n + ind] * c
	}
}

// Applique une multiplication à droite sur une matrice d'Hessenberg
// par une matrice de rotation de givens de coefficients c et s
// ! A doit être une matrice d'Hessenberg
void GivensRotationDH(int n, mpfr_t* A, int i, int j, mpfr_t c, mpfr_t s, mpfr_t tmp1, mpfr_t tmp2) {
	int ind;
	for(ind = 0 ; ind < i + 2 ; ind++) {
		mpfr_set(tmp1, A[ind*n + i], MPFR_RNDN);
		mpfr_mul(tmp2, A[ind*n + j], s, MPFR_RNDN);
		mpfr_mul(A[ind*n + i], A[ind*n + i], c, MPFR_RNDN);
		mpfr_add(A[ind*n + i], A[ind*n + i], tmp2, MPFR_RNDN);
		//// A[ind*n + i] = A[ind*n + i] * c + A[ind*n + j] * s
		
		mpfr_mul(tmp2, tmp1, s, MPFR_RNDN);
		mpfr_mul(A[ind*n + j], A[ind*n + j], c, MPFR_RNDN);
		mpfr_sub(A[ind*n + j], A[ind*n + j], tmp2, MPFR_RNDN);
		//// A[ind*n + j] = -1 * tmp1 * s + A[ind*n + j] * c
	}
}

// Mise à jour de la sous-diagonal
// Si un élément de celle-ci est suffisament petit (|A_(i+1)i| <= tol * (|A_ii| + |A_(i + 1)(i + 1)|))
// alors il est explicitement mis à zéro
void subDiagonalUpdate(int n, mpfr_t* A, mpfr_t tol, mpfr_t tmp1, mpfr_t tmp2) {
	int i;
	for(i = 0 ; i < n - 1 ; i++) {
		mpfr_abs(tmp2, A[i*n + i], MPFR_RNDN);
		mpfr_abs(tmp1, A[(i + 1)*n + (i + 1)], MPFR_RNDN);
		mpfr_add(tmp2, tmp2, tmp1, MPFR_RNDN);
		mpfr_mul(tmp2, tmp2, tol, MPFR_RNDN);
		mpfr_abs(tmp1, A[(i+1)*n + i], MPFR_RNDN);
		if(mpfr_lessequal_p(tmp1, tmp2))
			mpfr_set_d(A[(i+1)*n + i], 0, MPFR_RNDN);
	}
}

// Calcul p et q, tel que A soit diagonal par bloc de dimension au plus 2x2 jusqu'à l'indice p et de l'indice q à n
void deflation(int n, mpfr_t* A, int* p, int* q) {
	int i;
	for(i = *p ; i < n - 1 ; i++, (*p)++)
		if(mpfr_cmp_ui(A[(i+1)*n + i], 0) != 0 && (i == n - 2 || mpfr_cmp_ui(A[(i+2)*n + i + 1], 0) != 0))
			break;
	for(i = *q ; i > 0 ; i--, (*q)--)
		if(mpfr_cmp_ui(A[i*n + (i - 1)], 0) != 0 && (i == 1 || mpfr_cmp_ui(A[(i - 1)*n + (i - 2)], 0) != 0))
			break;
}

// Une étape de l'algorithme QR avec simple shift d'une matrice d'Hessenberg A
void QRstepH(int n, mpfr_t* A, mpfr_t* rot, int p, int q, mpfr_t c, mpfr_t s, mpfr_t tmp1, mpfr_t tmp2, mpfr_t mu) {
	int i;
	
	mpfr_set(mu, A[q*n + q], MPFR_RNDN);
	
	// Calcul de A - mu*I
	for(i = p ; i <= q ; i++)
		mpfr_sub(A[i*n + i], A[i*n + i], mu, MPFR_RNDN);
	
	// Calcul de la décomposition QR
	// R est stocké dans A
	// Q est stocké sous la forme des coefficients des rotations successives appliquées a A dans rot
	for(i = p ; i <= q - 1 ; i++) {
		GivensCoefficients(n, A, i, i + 1, i, c, s, tmp1);
		GivensRotationGH(n, A, i, i + 1, c, s, tmp1, tmp2);
		mpfr_set(rot[2*i], c, MPFR_RNDN);
		mpfr_set(rot[2*i + 1], s, MPFR_RNDN);
	}
	// Calcul de RQ
	for(i = p ; i <= q - 1 ; i++) {
		mpfr_mul_d(rot[2*i + 1], rot[2*i + 1], -1, MPFR_RNDN);
		GivensRotationDH(n, A, i, i + 1, rot[2*i], rot[2*i + 1], tmp1, tmp2);
	}
	
	// A + mu*I
	for(i = p ; i <= q ; i++)
		mpfr_add(A[i*n + i], A[i*n + i], mu, MPFR_RNDN);
}

// Itération de l'algorithme QR sur une matrice d'Hessenberg A
// En sortie, la matrice A est équivalente à la matrice A d'origine passé en argument
// et elle est diagonal par bloc de dimension au plus 2x2
int iterQR(int n, mpfr_t* A, mpfr_prec_t prec) {
	mpfr_t c, s, tmp1, tmp2, tol, mu;
	int i = 0;
	int p = 0, q = n-1;
	int tmpp = p, tmpq = q;
	mpfr_t* rot = allocateVector(2*n);
	
	mpfr_init2(c, prec);
	mpfr_init2(s, prec);
	mpfr_init2(tmp1, prec);
	mpfr_init2(tmp2, prec);
	mpfr_init2(tol, prec);
	mpfr_init2(mu, prec);
	
	// Tolérence de mise à zéro
	// Le plus petit nombre x tel que 1 + x != 1 (déterminé en fonction de la précision)
	mpfr_set_d(tol, 2, MPFR_RNDN);
	mpfr_pow_si(tol, tol, -1 * (prec - 1), MPFR_RNDN);
	
	initVector(rot, 2*n, prec);
	
	subDiagonalUpdate(n, A, tol, tmp1, tmp2);
	while(p < q) {
		i++;
		QRstepH(n, A, rot, p, q, c, s, tmp1, tmp2, mu); // Etape QR
		subDiagonalUpdate(n, A, tol, tmp1, tmp2); // Mise à jour de la sous-diagonal
		deflation(n, A, &p, &q); // Mise à jour des indices p et q
		 //if(p != tmpp || q != tmpq)
		//	printf("p = %d et q = %d et i = %d\n", p, q, i);
		tmpp = p; tmpq = q;
		if(i > 1000000)
		{
			fprintf(stderr, "ERROR: tried over 1 000 000 iterations for QR decomposition without success. \n");
			return 0;
		}
	}
	
	mpfr_clear(c);
	mpfr_clear(s);
	mpfr_clear(tmp1);
	mpfr_clear(tmp2);
	mpfr_clear(tol);
	mpfr_clear(mu);
	
	clearVector(rot, 2*n);
	return 1;
}

//// FONCTION DE MANIPULATION DE COMPLEXES

//// Decomposition et resolution se font pour des matrices de complexes
//// Re est la matrice des parties reelles et Im des parties imaginaires

// !! Pour l'addition et la soustraction, les variables de retour peuvent être les mêmes que celles d'entrées
void addComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2) {
	mpfr_add(renvRe, re1, re2, MPFR_RNDN);
	mpfr_add(renvIm, im1, im2, MPFR_RNDN);
}

void subComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2) {
	mpfr_sub(renvRe, re1, re2, MPFR_RNDN);
	mpfr_sub(renvIm, im1, im2, MPFR_RNDN);
}

void mulComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2, mpfr_t tmp) {
	mpfr_mul(renvRe, re1, re2, MPFR_RNDN);
	mpfr_mul(tmp, im1, im2, MPFR_RNDN);
	mpfr_sub(renvRe, renvRe, tmp, MPFR_RNDN);
	// renvRe = re1 * re2  - im1 * im2;
	
	mpfr_mul(renvIm, re1, im2, MPFR_RNDN);
	mpfr_mul(tmp, re2, im1, MPFR_RNDN);
	mpfr_add(renvIm, renvIm, tmp, MPFR_RNDN);
	// renvIm = re1 * im2 + re2 * im1;
}

void divComplex(mpfr_t renvRe, mpfr_t renvIm, mpfr_t re1, mpfr_t im1, mpfr_t re2, mpfr_t im2, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp) {
	mpfr_mul(tmp, re2, re2, MPFR_RNDN);
	mpfr_mul(tmpre, im2, im2, MPFR_RNDN);
	mpfr_add(tmp, tmp, tmpre, MPFR_RNDN);
	// tmp = re2 * re2 + im2 * im2;
	
	mpfr_div(tmpre, re2, tmp, MPFR_RNDN);
	mpfr_div(tmpim, im2, tmp, MPFR_RNDN);
	mpfr_mul_d(tmpim, tmpim, -1, MPFR_RNDN);
	// tmpre = re2 / tmp;
	// tmpim = -1 * im2 / tmp;
	
	mpfr_mul(renvRe, re1, tmpre, MPFR_RNDN);
	mpfr_mul(tmp, im1, tmpim, MPFR_RNDN);
	mpfr_sub(renvRe, renvRe, tmp, MPFR_RNDN);
	// renvRe = re1 * tmpre  - im1 * tmpim;
	
	mpfr_mul(renvIm, re1, tmpim, MPFR_RNDN);
	mpfr_mul(tmp, tmpre, im1, MPFR_RNDN);
	mpfr_add(renvIm, renvIm, tmp, MPFR_RNDN);
	// renvIm = re1 * tmpim + tmpre * im1;
}

// Donne la decomposition LU de la Matrice complexe H = Re + i * Im
// !!! H doit etre une matrice d'Hessenberg
void luDecompositionH(mpfr_t* Re, mpfr_t* Im, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp1, mpfr_t tmp2, mpfr_t tmp3) {
	int i, k;
	for(i = 0 ; i < n - 1 ; i++) {
		divComplex(tmpre, tmpim, Re[(i+1)*n + i], Im[(i+1)*n + i], Re[i*n + i], Im[i*n + i], tmp1, tmp2, tmp3);
		mpfr_set(Re[(i+1)*n + i], tmpre, MPFR_RNDN);
		mpfr_set(Im[(i+1)*n + i], tmpim, MPFR_RNDN);
		for(k = i + 1 ; k < n ; k++) {
			mulComplex(tmpre, tmpim, Re[(i+1)*n + i], Im[(i+1)*n + i], Re[i*n + k], Im[i*n + k], tmp1);
			subComplex(Re[(i+1)*n + k], Im[(i+1)*n + k], Re[(i+1)*n + k], Im[(i+1)*n + k], tmpre, tmpim);
		}
	}
}

// Resout le système Ux = b avec U = Re + i * Im une matrice triangulaire superieure complexe
// La solution x est placée dans le vecteur complexe b = bre + i * bim
void upperTriangularSolve(mpfr_t* Re, mpfr_t* Im, mpfr_t* bre, mpfr_t* bim, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp1, mpfr_t tmp2, mpfr_t tmp3) {
	int i, j;
	for(i = n - 1 ; i >= 0 ; i--) {
		for(j = n - 1 ; j > i ; j--) {
			mulComplex(tmpre, tmpim, Re[i*n + j], Im[i*n + j], bre[j], bim[j], tmp1);
			subComplex(bre[i], bim[i], bre[i], bim[i], tmpre, tmpim);
		}
		divComplex(tmpre, tmpim, bre[i], bim[i], Re[i*n + i], Im[i*n + i], tmp1, tmp2, tmp3);
		mpfr_set(bre[i], tmpre, MPFR_RNDN);
		mpfr_set(bim[i], tmpim, MPFR_RNDN);
	}
}

// Resout le système Lx = b avec L = Re + i * Im une matrice triangulaire inférieure complexe avec des 1 sur la diagonale
// La solution x est placée dans le vecteur complexe b = bre + i * bim
void lowerTriangularSolve(mpfr_t* Re, mpfr_t* Im, mpfr_t* bre, mpfr_t* bim, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp) {
	int i, j;
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < i ; j++) {
			mulComplex(tmpre, tmpim, Re[i*n + j], Im[i*n + j], bre[j], bim[j], tmp);
			subComplex(bre[i], bim[i], bre[i], bim[i], tmpre, tmpim);
		}
}

// Resout Ax = b avec A = LU une matrice complexe à laquelle on a déjà appliqué la fonction luDecompositionH
// La solution x est placée dans le vecteur complexe b = bre + i * bim
void systemSolveLU(mpfr_t* Re, mpfr_t* Im, mpfr_t* bre, mpfr_t* bim, int n, mpfr_t tmp1, mpfr_t tmp2, mpfr_t tmp3, mpfr_t tmp4, mpfr_t tmp5) {
	lowerTriangularSolve(Re, Im, bre, bim, n, tmp1, tmp2, tmp3);
	upperTriangularSolve(Re, Im, bre, bim, n, tmp1, tmp2, tmp3, tmp4, tmp5);
}

// Calcul C = A*B
void mulMatrixComplexReal(mpfr_t* Cre, mpfr_t* Cim, mpfr_t* Are, mpfr_t* Aim, mpfr_t* Bre, int n, mpfr_t tmpre, mpfr_t tmpim, mpfr_t tmp1, mpfr_t tmp2) {
	int i, j, k;
	mpfr_set_d(tmp1, 0, MPFR_RNDN);
	for(i = 0 ; i < n ; i++)
		for(j = 0 ; j < n ; j++) {
			mpfr_set_d(Cre[i*n + j], 0, MPFR_RNDN);
			mpfr_set_d(Cim[i*n + j], 0, MPFR_RNDN);
			for(k = 0 ; k < n ; k++) {
				mulComplex(tmpre, tmpim, Are[i*n + k], Aim[i*n + k], Bre[k*n + j], tmp1, tmp2);
				addComplex(Cre[i*n + j], Cim[i*n + j], Cre[i*n + j], Cim[i*n + j], tmpre, tmpim);
			}
		}
}
////

// Trouve les valeurs propres de A, une matrice ayant subit l'algorithme QR
// re contient la partie réelle des valeurs propres et im la partie imaginaire 
int findEigenValues(int n, mpfr_t* A, mpfr_t* re, mpfr_t* im, mpfr_prec_t prec) {
	int i = 0;
	int success = 1;
	mpfr_t a, b, c, d, tmp, disc;
	
	mpfr_init2(a, prec);
	mpfr_init2(b, prec);
	mpfr_init2(c, prec);
	mpfr_init2(d, prec);
	mpfr_init2(tmp, prec);
	mpfr_init2(disc, prec);
	while(i < n) {
		if(i != n - 1 && mpfr_cmp_ui(A[(i+1)*n + i], 0) != 0) { // Si on se trouve sur un bloc diagonal 2x2, on calcul les racines du polynôme caractéristique de celui-ci
			mpfr_set(a, A[i*n + i], MPFR_RNDN);
			mpfr_set(b, A[i*n + i + 1], MPFR_RNDN);
			mpfr_set(c, A[(i+1)*n + i], MPFR_RNDN);
			mpfr_set(d, A[(i+1)*n + i + 1], MPFR_RNDN);
			
			mpfr_mul(disc, c, b, MPFR_RNDN);
			mpfr_mul(tmp, a, d, MPFR_RNDN);
			mpfr_sub(disc, tmp, disc, MPFR_RNDN);
			mpfr_mul_d(disc, disc, 4, MPFR_RNDN);
			mpfr_add(tmp, a, d, MPFR_RNDN);
			mpfr_pow_ui(tmp, tmp, 2, MPFR_RNDN);
			mpfr_sub(tmp, tmp, disc, MPFR_RNDN);

			
			if(mpfr_cmp_ui(tmp, 0) < 0) {
				mpfr_mul_d(tmp, tmp, -1, MPFR_RNDN);
				mpfr_sqrt(disc, tmp, MPFR_RNDN);
				
				mpfr_add(re[i], a, d, MPFR_RNDN);
				mpfr_div_d(re[i], re[i], 2, MPFR_RNDN);
				
				mpfr_set(re[i + 1], re[i], MPFR_RNDN);
				
				mpfr_div_d(im[i], disc, 2, MPFR_RNDN);
				mpfr_mul_d(im[i + 1], im[i], -1, MPFR_RNDN);

				i += 2;
			}
			else {
				mpfr_sqrt(disc, tmp, MPFR_RNDN);
				
				mpfr_add(re[i], a, d, MPFR_RNDN);
				mpfr_set(re[i + 1], re[i], MPFR_RNDN);
				
				mpfr_add(re[i], re[i], disc, MPFR_RNDN);
				mpfr_sub(re[i + 1], re[i + 1], disc, MPFR_RNDN);
				
				mpfr_div_d(re[i], re[i], 2, MPFR_RNDN);
				mpfr_div_d(re[i + 1], re[i + 1], 2, MPFR_RNDN);

				mpfr_set_d(im[i], 0, MPFR_RNDN);
				mpfr_set_d(im[i + 1], 0, MPFR_RNDN);

				i += 2;
			}
		}
		else { // Sinon, on lit la valeur propre sur le coefficient diagonal
			mpfr_set(re[i], A[i*n + i], MPFR_RNDN);
			mpfr_set_d(im[i], 0, MPFR_RNDN);
			i++;
		}
	}

	
	mpfr_clear(a);
	mpfr_clear(b);
	mpfr_clear(c);
	mpfr_clear(d);
	mpfr_clear(tmp);
	mpfr_clear(disc);

	return success;
}

// Calcul les vecteurs propre de la matrice A associés aux valeurs propres complexes val = re + i * im
int findEigenVectors(int n, mpfr_t* A, mpfr_t* eigVecre, mpfr_t* eigVecim, mpfr_t* re, mpfr_t* im, mpfr_prec_t prec) {

	int success = 0;
	int i, j, compt, essai;
	
	mpfr_t tol, tmp1, tmp2, tmp3, tmp4, tmp5;
	
	mpfr_t* Are = allocateMatrix(n);
	mpfr_t* Aim = allocateMatrix(n);
	
	mpfr_t* tmpre = allocateVector(n);
	mpfr_t* tmpim = allocateVector(n);
	
	mpfr_init2(tol, prec);
	mpfr_init2(tmp1, prec);
	mpfr_init2(tmp2, prec);
	mpfr_init2(tmp3, prec);
	mpfr_init2(tmp4, prec);
	mpfr_init2(tmp5, prec);
	
	initMatrix(Are, n, prec);
	initMatrix(Aim, n, prec);
	
	initVector(tmpre, n, prec);
	initVector(tmpim, n, prec);
	
	// Calcul de la tolérance sur le calcul des vecteurs propres
	mpfr_set_d(tol, 2, MPFR_RNDN);
	mpfr_pow_si(tol, tol, -1 * (prec - 3), MPFR_RNDN);
	normInfMatrix(tmp1, A, n, tmp2, tmp3);
	mpfr_mul(tol, tol, tmp1, MPFR_RNDN);
	// tol = 2^(3 - prec) * normInfMatrix(A, n);
	
	// Calcul de chaque vecteurs propres individuellement par une "inverse iteration"
	// On part d'un vecteur aléatoire
	// On effectue l'"inverse iteration" sur lui
	// Elle va écraser les composantes associées aux valeurs propres qui n'est pas celle du vecteur propre calculé
	// Cette méthode, si elle converge, converge très vite (1, voir 2 itérations suffisent)
	// Sinon, cela signifie que le vecteur aléatoire choisi en premier lieu ne possèdait pas une composante
	// Suffisament grande avec la valeur propre concernée et on réessaie avec un nouveau vecteur aléatoire
	for(i = 0 ; i < n ; i++) {
		essai = 0;
		
		// Copie de A pour faire les calculs sans modifier A
		eigcode_copyMatrix(Are, A, n);
		set0Matrix(Aim, n);
		for(j = 0 ; j < n ; j++) {
			mpfr_sub(Are[j*n + j], Are[j*n + j], re[i], MPFR_RNDN);
			// Are[j*n + j] -= re[i];
			mpfr_sub(Aim[j*n + j], Aim[j*n + j], im[i], MPFR_RNDN);
			// Aim[j*n + j] -= im[i];
		}

		// Decomposition pour la résolution de système linéaire
		luDecompositionH(Are, Aim, n, tmp1, tmp2, tmp3, tmp4, tmp5);
		do {
			compt = 0;
			essai++;
			
			// On choisi un vecteur aléatoire réel ou complexe selon que la valeur propre soit réelle ou complexe
			setAleaVector(eigVecre + i*n, n, tmp1);
			if(mpfr_cmp_d(im[i], 0) == 0)
				set0Vector(eigVecim + i*n, n);
			else
				setAleaVector(eigVecim + i*n, n, tmp1);
			
			// Iteration
			do {
				compt++;
				copyVector(tmpre, eigVecre + i*n, n);
				copyVector(tmpim, eigVecim + i*n, n);

				systemSolveLU(Are, Aim, eigVecre + i*n, eigVecim + i*n, n, tmp1, tmp2, tmp3, tmp4, tmp5);
				
				norm2Complex(tmp1, n, eigVecre + i*n, eigVecim + i*n, tmp2);

				for(j = 0 ; j < n ; j++) {
					mpfr_div(eigVecre[i*n + j], eigVecre[i*n + j], tmp1, MPFR_RNDN);
					// eigVecre[i*n + j] /= tmp;
					mpfr_div(eigVecim[i*n + j], eigVecim[i*n + j], tmp1, MPFR_RNDN);
					// eigVecim[i*n + j] /= tmp;
					mpfr_div(tmpre[j], tmpre[j], tmp1, MPFR_RNDN);
					// tmpre[j] /= tmp;
					mpfr_div(tmpim[j], tmpim[j], tmp1, MPFR_RNDN);
					// tmpim[j] /= tmp;
				}
				normInfVectComplex(tmp1, tmpre, tmpim, n, tmp2);
			} while(mpfr_cmp(tmp1, tol) > 0 && compt <= 3); // S'arrête si le vecteur propre a convergé ou si on déjà fait trop d'itérations
		}  while(compt > 3 && essai <= n*n); // S'arrête si on a pas fait trop d'itérations (donc que le vecteur propre a convergé) ou si on a fait trop d'essaie (pour éviter le programme de mouliner dans le vide)
		if(mpfr_cmp_d(im[i], 0) != 0) { // Si on avait une valeur propre complexe, on affecte à la valeur propre conjuguée, le vecteur propre conjugué
			copyVector(eigVecre + (i + 1)*n, eigVecre + i*n, n);
			for(j = 0 ; j < n ; j++)
				mpfr_mul_d(eigVecim[(i + 1)*n + j], eigVecim[i*n + j], -1, MPFR_RNDN);
			i++;
		}
		
		if(essai > n*n) {
			if(mpfr_cmp_d(im[i], 0) != 0)
				i += 2;
			else
				i++;
			//printf("Calcul du Vecteur propre de la valeur propre numéro %d n'a pas aboutie au bout de %d essais, sauté\n", i, essai);
			success = 0;
		}
		//printf("Vecteur propre de la valeur propre numéro %d calculé en %d essais\n", i, essai);
	}
	
	mpfr_clear(tol);
	mpfr_clear(tmp1);
	mpfr_clear(tmp2);
	mpfr_clear(tmp3);
	mpfr_clear(tmp4);
	mpfr_clear(tmp5);
	
	clearMatrix(Are, n);
	clearMatrix(Aim, n);
	
	clearVector(tmpre, n);
	clearVector(tmpim, n);
	return success;
}

// Calcul des valeurs propres d'une matrice quelconque et, si ask_for_Evector est différents de 0, des vecteurs propres
/*Function returns
 * 	0 if could not compute eigenvalues (and consequently, eigenvectors)
 * 	1 if computed eigenvalues but failed to compute eigenvectors
 * 	2 if computed eigenvalues and eigenvectors
 */
int mpfr_eigen(int n, mpfr_t* A, mpfr_t* re, mpfr_t* im, int ask_for_Evector, mpfr_t* evre, mpfr_t* evim, mpfr_prec_t prec) {
	srand(time(NULL));
	


	// Variables temporaires utilisées dans les fonctions
	mpfr_t tmp1, tmp2, tmp3, tmp4;
	mpfr_init2(tmp1, prec);
	mpfr_init2(tmp2, prec);
	mpfr_init2(tmp3, prec);
	mpfr_init2(tmp4, prec);
	//
	
	// Matrice de sauvegarde des transformations dans la réduction en matrice d'Hessenberg
	mpfr_t* U = allocateMatrix(n);
	initMatrix(U, n, prec);
	identityMatrix(U, n);

	mpfr_t* tmpevre = allocateMatrix(n);
	mpfr_t* tmpevim = allocateMatrix(n);
	mpfr_t* H = allocateMatrix(n);
	
	initMatrix(tmpevre, n, prec);
	initMatrix(tmpevim, n, prec);
	initMatrix(H, n, prec);
	


	toHessenbergForm(n, A, U, prec);

	trans(U, n, tmp1);
	
	eigcode_copyMatrix(H, A, n); // Copie de la reduction d'Hessenberg de A dans H avant de continuer les calculs

	if(!iterQR(n, A, prec)) // Réduction de A à une matrice quasi-triangulaire
	{
		fprintf(stderr, "ERROR: could not compute the QR decomposition of matrix A.\n");
		return 0;
	}

	if(!findEigenValues(n, A, re, im, prec)) // Lecture des valeurs propres
	{
		fprintf(stderr, "ERROR in MPFR eigensolver: could not compute eigenvalues. \n");
		mpfr_clear(tmp1);
		mpfr_clear(tmp2);
		mpfr_clear(tmp3);
		mpfr_clear(tmp4);

		clearMatrix(U, n);
		clearMatrix(tmpevre, n);
		clearMatrix(tmpevim, n);
		clearMatrix(H, n);
		return 0;
	}
	if(ask_for_Evector) { // Si les vecteurs propres sont demandés

		if(!findEigenVectors(n, H, tmpevre, tmpevim, re, im, prec)) // Calcul des vecteurs propres de H (forme d'Hessenberg de A)
		{
			fprintf(stderr, "WARNING in MPFR eigensolver: could not compute eigenvectors.\nWARNING: returning only the eigenvalues \n");
			mpfr_clear(tmp1);
			mpfr_clear(tmp2);
			mpfr_clear(tmp3);
			mpfr_clear(tmp4);

			clearMatrix(U, n);
			clearMatrix(tmpevre, n);
			clearMatrix(tmpevim, n);
			clearMatrix(H, n);
			return 1;
		}
		mulMatrixComplexReal(evre, evim, tmpevre, tmpevim, U, n, tmp1, tmp2, tmp3, tmp4); // Multiplication par U pour trouver les vecteurs propres de A
		
		// Les vecteurs sont organisés lignes par lignes
		// On va les organiser colonnes par colonnes
		trans(evre, n, tmp1);
		trans(evim, n, tmp1);
	}
	
	// Libération
	mpfr_clear(tmp1);
	mpfr_clear(tmp2);
	mpfr_clear(tmp3);
	mpfr_clear(tmp4);
	
	clearMatrix(U, n);
	clearMatrix(tmpevre, n);
	clearMatrix(tmpevim, n);
	clearMatrix(H, n);
	return 2;
}

// Fonction d'écriture de complex dans un fichier
void exportVectorComplex(mpfr_t* re, mpfr_t* im, int n, char* filename) {
	int i;
	FILE* f = fopen(filename, "w");
	for(i = 0 ; i < n ; i++) {
		mpfr_out_str(f, 10, 0, re[i], MPFR_RNDN);
		
		if(mpfr_cmp_ui(im[i], 0) != 0) {
			fprintf(f, " + i * ");
			mpfr_out_str(f, 10, 0, im[i], MPFR_RNDN);
		}
		
		fprintf(f, "\n");
	}
	fclose(f);
}

void writeVectorComplex(mpfr_t* re, mpfr_t* im, int n, FILE *f) {
	int i;
	for(i = 0 ; i < n ; i++) {
		mpfr_out_str(f, 10, 0, re[i], MPFR_RNDN);

		if(mpfr_cmp_ui(im[i], 0) != 0) {
			fprintf(f, " + i * ");
			mpfr_out_str(f, 10, 0, im[i], MPFR_RNDN);
		}

		fprintf(f, "\n");
	}
	fclose(f);
}

// Fonction d'écriture de complex dans un fichier
void exportMatrixComplex(mpfr_t* re, mpfr_t* im, int n, char* filename) {
	int i, j;
	FILE* f = fopen(filename, "w");
	for(i = 0 ; i < n ; i++, fprintf(f, "\n"))
		for(j = 0 ; j < n ; j++, fprintf(f, "\t")) {
			mpfr_out_str(f, 10, 0, re[i*n + j], MPFR_RNDN);
			if(mpfr_cmp_ui(im[i*n + j], 0) != 0) {
				fprintf(f, " + i * ");
				mpfr_out_str(f, 10, 0, im[i*n + j], MPFR_RNDN);
			}
		}
	fclose(f);
}














